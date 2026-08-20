# app/R/deconv_module.R
# ---------------------------------------------------------------------------
# Deconvolution view: per-sample cell-type composition from deconvMe
# (deconv/cell_fractions.csv), plus a comparison of the summed immune fraction
# against sesame's Leukocyte_Fraction.
#
# These answer DIFFERENT questions: deconvMe = immune COMPOSITION (which cell
# types); sesame leukocyte = AMOUNT of immune (a purity scalar). On tumours the
# blood references over-estimate the immune magnitude, so the summed fraction is
# expected to differ from sesame's number — compare composition, not magnitude.
# ---------------------------------------------------------------------------

# Cell types treated as the non-immune residual (excluded from the immune sum).
DECONV_NONIMMUNE <- c("other", "unknown", "Other", "Unknown")

deconv_module_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::p(
      class = "text-muted",
      "Reference-based immune cell-type composition (deconvMe). Each bar is a ",
      "sample; segments are the estimated cell-type fractions. This is the ",
      "immune ", shiny::tags$em("composition"), " — sesame's Leukocyte_Fraction ",
      "(QC tab) is the immune ", shiny::tags$em("amount"), ". On tumours the ",
      "blood reference over-estimates the magnitude, so the summed fraction ",
      "below will not match sesame's number — that's expected, not a failure."
    ),
    shiny::uiOutput(ns("method_selector")),
    bslib::navset_tab(
      bslib::nav_panel(
        "Composition",
        shinycssloaders::withSpinner(
          plotly::plotlyOutput(ns("barplot"), height = "auto"),
          type = 7, color = COLORS$primary
        )
      ),
      bslib::nav_panel(
        "vs sesame leukocyte",
        shinycssloaders::withSpinner(
          DT::DTOutput(ns("compare_table")),
          type = 7, color = COLORS$primary
        )
      ),
      bslib::nav_panel(
        "Fractions table",
        shinycssloaders::withSpinner(
          DT::DTOutput(ns("fractions_table")),
          type = 7, color = COLORS$primary
        )
      )
    )
  )
}

deconv_module_server <- function(id, results) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    deconv <- shiny::reactive({
      r <- results()
      if (is.null(r) || is.null(r$deconv)) return(NULL)
      d <- r$deconv
      need <- c("method", "sample", "celltype", "value")
      if (!all(need %in% names(d))) return(NULL)
      d
    })

    output$method_selector <- shiny::renderUI({
      d <- deconv()
      if (is.null(d)) {
        return(shiny::div(
          class = "alert alert-secondary",
          "No deconvolution results. Run the deconvolution step ",
          "(--step deconvolution) on a 450k/EPIC run."
        ))
      }
      ms <- unique(d$method)
      # Show the aggregated result first if present.
      ms <- c(intersect(c("aggregated"), ms), setdiff(ms, "aggregated"))
      shiny::selectInput(ns("method"), "Method", choices = ms, selected = ms[1])
    })

    sel <- shiny::reactive({
      d <- deconv(); m <- input$method
      if (is.null(d) || is.null(m)) return(NULL)
      d[d$method == m, , drop = FALSE]
    })

    output$barplot <- plotly::renderPlotly({
      s <- sel()
      if (is.null(s) || !nrow(s)) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") |>
          plotly::layout(title = list(text = "No deconvolution data.",
                                      font = list(size = 13))))
      }
      n <- length(unique(s$sample))
      p <- ggplot2::ggplot(
        s, ggplot2::aes(x = sample, y = value, fill = celltype)) +
        ggplot2::geom_col(position = "stack") +
        ggplot2::labs(x = NULL, y = "Fraction", fill = "Cell type") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      plotly::ggplotly(p, height = max(360, 60 + 22 * n)) |>
        plotly::layout(margin = list(b = 140))
    })

    output$compare_table <- DT::renderDT({
      s <- sel(); r <- results()
      if (is.null(s)) return(NULL)
      imm <- stats::aggregate(
        value ~ sample,
        data = s[!s$celltype %in% DECONV_NONIMMUNE, , drop = FALSE], FUN = sum)
      names(imm) <- c("Sample_ID", "Deconv_Immune_Sum")
      imm$Deconv_Immune_Sum <- round(imm$Deconv_Immune_Sum, 3)
      # Join sesame Leukocyte_Fraction from the QC report, if loaded.
      if (!is.null(r) && !is.null(r$qc_report) &&
          all(c("Sample_ID", "Leukocyte_Fraction") %in% names(r$qc_report))) {
        leuk <- r$qc_report[, c("Sample_ID", "Leukocyte_Fraction")]
        imm <- merge(imm, leuk, by = "Sample_ID", all.x = TRUE)
        imm$Leukocyte_Fraction <- round(as.numeric(imm$Leukocyte_Fraction), 3)
      }
      DT::datatable(imm, rownames = FALSE, selection = "none",
                    class = "stripe hover compact", width = "100%",
                    options = list(pageLength = 25, dom = "ltip", scrollX = TRUE))
    })

    output$fractions_table <- DT::renderDT({
      s <- sel()
      if (is.null(s)) return(NULL)
      # Wide: samples x cell types for the selected method.
      w <- stats::reshape(
        s[, c("sample", "celltype", "value")],
        idvar = "sample", timevar = "celltype", direction = "wide")
      names(w) <- sub("^value\\.", "", names(w))
      num <- vapply(w, is.numeric, logical(1))
      w[num] <- lapply(w[num], function(x) round(x, 3))
      DT::datatable(w, rownames = FALSE, selection = "none",
                    class = "stripe hover compact", width = "100%",
                    options = list(pageLength = 25, dom = "ltip", scrollX = TRUE))
    })
  })
}
