# app/R/qc_module.R
# ---------------------------------------------------------------------------
# QC view: per-sample DT table + the pipeline's own interactive QC HTML
# plots embedded via iframe.
#
# We don't re-render QC plots from the underlying minfi QC object — the
# pipeline's `figures/qc/interactive_*.html` files are already plotly
# outputs self-contained with all their assets. Embedding them is faster,
# smaller, and matches what the shareable HTML report already shows.
# ---------------------------------------------------------------------------

# Plain-English explanations attached as native title attributes to the QC
# table's column headers. Hover-to-read; no Bootstrap init needed.
QC_COL_TOOLTIPS <- list(
  Sample_ID               = "Sentrix ID identifying this sample.",
  Mean_Detection_P        = "Average detection p-value across all probes for this sample. Lower is better; healthy samples typically sit near 0.001.",
  Failed_Probes_Count     = "Number of probes whose detection p-value exceeds the per-probe threshold (default 0.01).",
  Failed_Probes_Percent   = "Percent of probes that failed the detection-p check, within this sample. Compared against the failed-probe % threshold.",
  Median_Meth_Intensity   = "log2 median methylated-channel intensity (raw, pre-normalization).",
  Median_Unmeth_Intensity = "log2 median unmethylated-channel intensity (raw, pre-normalization).",
  Pass_QC                 = "TRUE = sample passed both the mean-detection-p and failed-probe checks.",
  Flag_Mean_DetP          = "TRUE when the sample's mean detection-p exceeds the cohort threshold (default 0.05).",
  Flag_Failed_Probes      = "TRUE when the failed-probe percent reaches the threshold (default 25%).",
  GCT_Score               = "Bisulfite-conversion control (GCT, Zhou et al. 2017). ~1.0 = complete conversion; higher = more incomplete. Computed for 450k/EPIC/EPICv2.",
  Flag_GCT                = "TRUE when the GCT score exceeds the conversion threshold (default 1.3) — incomplete bisulfite conversion. Contributes to Pass_QC (fails the sample).",
  Sesame_Sex              = "Predicted sex from sesame's curated X/Y probe model (MALE/FEMALE). Informational — compare against expected sex and the minfi prediction to catch sample swaps.",
  Karyotype               = "Coarse sex karyotype from the predicted sex plus X-inactivation heterozygosity (X_Het): XX or XY. 'XY (low Y - possible LOY)' = one X by methylation but depleted Y intensity (minfi yMed-xMed < -2) — likely somatic Loss-of-Y, common in tumours, and the usual cause of a sesame-MALE / minfi-FEMALE disagreement. A trailing '?' (XXY?, X0?) marks an UNCERTAIN aneuploidy guess to verify independently. Informational, never gates Pass_QC.",
  X_Het                   = "Fraction of X-linked probes with intermediate beta (0.3–0.7). Two active X chromosomes (X-inactivation) give ~0.45; a single X gives ~0.13. Drives the karyotype call.",
  Horvath_Age             = "Predicted epigenetic age (years) from the Horvath 353-CpG clock (Horvath 2013). Informational — a large gap from the known age can flag a mislabelled sample.",
  Leukocyte_Fraction      = "Estimated leukocyte (white-blood-cell) fraction from sesame's two-component model (0–1). Informational — gauges immune/normal-cell contamination in a tumour sample. EPICv2 is converted to EPIC space (Zhou Lab map) so it works across 450k/EPIC/EPICv2.",
  SNP_BestMatch           = "The other sample in this run with the most similar SNP genotype fingerprint — i.e. the closest genetic match. Pair this with SNP_Match_Pct.",
  SNP_Match_Pct           = "Genotype concordance (%) with SNP_BestMatch. ~95–100% = same individual (a true replicate / tumour-normal pair, OR a sample swap if they shouldn't match); ~35–50% = unrelated. See the Sample identity tab for the full pairwise heatmap.",
  SNP_Count               = "Number of usable rs SNP probes behind the identity match (typically ~59–65). Lower counts mean a less reliable match.",
  Note_Low_Intensity      = "Informational only — does NOT contribute to Pass_QC. Flags samples with low median intensities, often a scanner-gain issue that SWAN normalization can recover.",
  SWAN_Median_Meth        = "Median methylated intensity AFTER SWAN normalization. Computed only for low-intensity samples.",
  SWAN_Median_Unmeth      = "Median unmethylated intensity after SWAN normalization. Computed only for low-intensity samples.",
  SWAN_Recoverable        = "TRUE when SWAN normalization brings intensities above threshold — the low-intensity flag was a scanner-gain artifact, not a true failure."
)

qc_module_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::navset_tab(
    bslib::nav_panel(
      "Sample metrics",
      shiny::uiOutput(ns("summary_banner")),
      shinycssloaders::withSpinner(
        DT::DTOutput(ns("qc_table")),
        type = 7, color = COLORS$primary
      )
    ),
    bslib::nav_panel(
      "Sample identity",
      shiny::p(
        class = "text-muted",
        "Pairwise SNP genotype concordance (%). ~95–100% = same individual ",
        "(replicate / tumour-normal pair, or an unexpected swap); ~35–50% = ",
        "unrelated. Computed from the Infinium rs probes."
      ),
      shinycssloaders::withSpinner(
        plotly::plotlyOutput(ns("snp_heatmap"), height = "auto"),
        type = 7, color = COLORS$primary
      )
    ),
    bslib::nav_panel(
      "Dye bias",
      shiny::p(
        class = "text-muted",
        "Per-sample Red–Green QQ plot (sesame). Points on the diagonal mean ",
        "the two colour channels are balanced; a strong departure indicates ",
        "dye bias (the pipeline corrects it downstream via noob/dyeBiasNL)."
      ),
      shiny::uiOutput(ns("dye_bias_selector")),
      shiny::uiOutput(ns("dye_bias_frame"))
    ),
    bslib::nav_panel(
      "Density plot",
      shiny::uiOutput(ns("density_frame"))
    ),
    bslib::nav_panel(
      "MDS plot",
      shiny::uiOutput(ns("mds_frame"))
    )
  )
}

qc_module_server <- function(id, results) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$summary_banner <- shiny::renderUI({
      r <- results()
      if (is.null(r)) {
        return(shiny::div(class = "alert alert-secondary",
                          "No completed run yet. Run the pipeline first ",
                          "(Run tab), then this view populates automatically."))
      }
      if (is.null(r$qc_report)) {
        return(shiny::div(class = "alert alert-secondary",
                          "QC results not available yet — waiting for the ",
                          "QC step to finish."))
      }
      n_total <- nrow(r$qc_report)
      n_fail  <- length(r$qc_fail_ids)
      cls <- if (n_fail == 0) "alert alert-success" else "alert alert-warning"
      msg <- if (n_fail == 0) {
        sprintf("%d of %d sample(s) passed QC.", n_total, n_total)
      } else {
        sprintf("%d of %d sample(s) passed QC. %d failed: %s",
                n_total - n_fail, n_total, n_fail,
                paste(r$qc_fail_ids, collapse = ", "))
      }
      shiny::div(class = cls, role = "alert", msg)
    })

    output$qc_table <- DT::renderDT({
      r <- results()
      if (is.null(r) || is.null(r$qc_report)) return(NULL)
      df <- r$qc_report
      # Pass_QC can read as logical, character, or integer depending on the
      # upstream CSV. Normalize to character "TRUE"/"FALSE" so DT's JS-side
      # styleEqual match is deterministic.
      if ("Pass_QC" %in% colnames(df)) {
        df$Pass_QC <- ifelse(as.logical(df$Pass_QC), "TRUE", "FALSE")
      }
      numeric_cols <- which(vapply(df, is.numeric, logical(1)))
      dt <- DT::datatable(
        df,
        rownames = FALSE,
        selection = "none",
        class = "stripe hover compact",
        width = "100%",
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          autoWidth = FALSE,
          dom = "ltip",
          headerCallback = dt_header_tooltips(QC_COL_TOOLTIPS),
          columnDefs = if (length(numeric_cols)) list(list(
            className = "dt-right",
            targets = as.integer(numeric_cols - 1L)
          )) else list()
        )
      )
      if ("Pass_QC" %in% colnames(df)) {
        # Shade failed-sample rows pink. Under bslib's Bootstrap 5 theme DT
        # paints a per-cell background (--bs-table-bg, plus the stripe
        # box-shadow), which hides a row-level (<tr>) background. So style every
        # cell keyed off Pass_QC (target = "cell" across all columns) rather than
        # using target = "row", which renders invisible under BS5.
        dt <- DT::formatStyle(
          dt, columns = colnames(df), valueColumns = "Pass_QC", target = "cell",
          backgroundColor = DT::styleEqual(
            c("TRUE", "FALSE"), c(NA, "#fde2e4")
          )
        )
      }
      dt
    })

    output$snp_heatmap <- plotly::renderPlotly({
      r <- results()
      m <- if (!is.null(r)) r$snp_concordance else NULL
      if (is.null(m) || nrow(m) < 2) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") |>
          plotly::layout(title = list(
            text = "Sample-identity heatmap needs at least 2 samples with SNP data.",
            font = list(size = 13))))
      }
      ids <- rownames(m)
      # Plotly draws the y-axis bottom-up; reverse rows so the diagonal reads
      # top-left to bottom-right like the column order.
      plotly::plot_ly(
        x = ids, y = rev(ids), z = m[rev(seq_len(nrow(m))), , drop = FALSE],
        type = "heatmap", zmin = 0, zmax = 100,
        colors = grDevices::colorRampPalette(c("#f7fbff", COLORS$primary))(64),
        hovertemplate = "%{y}<br>vs %{x}<br>concordance: %{z}%<extra></extra>",
        colorbar = list(title = "% match")
      ) |>
        plotly::layout(
          xaxis = list(tickangle = -45, title = ""),
          yaxis = list(title = ""),
          margin = list(l = 140, b = 140)
        )
    })

    output$dye_bias_selector <- shiny::renderUI({
      r <- results()
      ids <- if (!is.null(r)) {
        d <- file.path(r$run_dir, "figures", "qc", "dye_bias")
        if (dir.exists(d)) sub("\\.png$", "", list.files(d, pattern = "\\.png$"))
        else character(0)
      } else character(0)
      if (!length(ids)) {
        return(shiny::tags$em(class = "text-muted",
                              "No dye-bias plots yet (run preprocessing)."))
      }
      shiny::selectInput(ns("dye_sample"), label = NULL,
                         choices = ids, selected = ids[1])
    })

    output$dye_bias_frame <- shiny::renderUI({
      r <- results()
      sid <- input$dye_sample
      if (is.null(r) || is.null(sid)) {
        return(shiny::div(class = "alert alert-secondary",
                          "Select a sample to view its dye-bias QQ plot."))
      }
      rel <- file.path("figures", "qc", "dye_bias", paste0(sid, ".png"))
      if (!file.exists(file.path(r$run_dir, rel))) {
        return(shiny::div(class = "alert alert-warning",
                          sprintf("Dye-bias plot not found for %s.", sid)))
      }
      shiny::img(
        src = paste0(r$run_url_base, "/", rel),
        style = "max-width: 700px; width: 100%; border: 1px solid #dee2e6; border-radius: 4px;"
      )
    })

    output$density_frame <- shiny::renderUI({
      qc_iframe(results(),
                rel_path = "figures/qc/interactive_density_plot.html",
                fallback = "Interactive density plot not found.")
    })

    output$mds_frame <- shiny::renderUI({
      qc_iframe(results(),
                rel_path = "figures/qc/interactive_mds_plot.html",
                fallback = "Interactive MDS plot not found.")
    })
  })
}

# Render an iframe pointing at a file under the run's URL base. If the
# file doesn't exist on disk, render a friendly placeholder instead.
qc_iframe <- function(results, rel_path, fallback) {
  if (is.null(results)) {
    return(shiny::div(class = "alert alert-secondary",
                      "No completed run yet."))
  }
  disk_path <- file.path(results$run_dir, rel_path)
  if (!file.exists(disk_path)) {
    return(shiny::div(class = "alert alert-warning", fallback))
  }
  url <- paste0(results$run_url_base, "/", rel_path)
  shiny::tags$iframe(
    src = url,
    style = "width: 100%; height: 75vh; border: 1px solid #dee2e6; border-radius: 4px;"
  )
}
