# Visualization module for methylation array analysis pipeline

#' Generate comprehensive report
#'
#' @param qc_results QC results
#' @param dim_reduction Dimensionality reduction results
#' @param cnv_data CNV analysis results
#' @param reference_projection Reference-projection results (the rp_result
#'   list from run_reference_projection: dataset, projected, class_hints, ...)
#' @param sample_info Sample information data frame
#' @param output_dirs Named list of module output directories (from setup_directories())
#' @param output_dir Output directory for reports
#' @return List of generated report paths
generate_report <- function(qc_results = NULL,
                          dim_reduction = NULL,
                          cnv_data = NULL,
                          reference_projection = NULL,
                          sample_info = NULL,
                          output_dirs = NULL,
                          output_dir = ".") {

  report_files <- list()

  # ------------------------------------------------------------------
  # Check prerequisites for HTML rendering
  # ------------------------------------------------------------------
  has_rmarkdown <- requireNamespace("rmarkdown", quietly = TRUE)

  # R from CRAN ships no pandoc (only RStudio bundles one), so a plain
  # double-click install has none on PATH. setup.R provisions a managed copy
  # via the `pandoc` package in that case; wire it in here. find_pandoc()
  # checks RSTUDIO_PANDOC first, so point it at the managed binary's dir and
  # clear the cache so the new location takes effect.
  if (has_rmarkdown && nchar(rmarkdown::find_pandoc()$dir) == 0 &&
      requireNamespace("pandoc", quietly = TRUE) &&
      isTRUE(tryCatch(pandoc::pandoc_available(), error = function(e) FALSE))) {
    Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc::pandoc_bin()))
    rmarkdown::find_pandoc(cache = FALSE)
    message("Using managed pandoc: ", Sys.getenv("RSTUDIO_PANDOC"))
  }

  has_pandoc    <- has_rmarkdown && nchar(rmarkdown::find_pandoc()$dir) > 0
  has_DT        <- requireNamespace("DT", quietly = TRUE)
  has_knitr     <- requireNamespace("knitr", quietly = TRUE)

  render_html <- has_rmarkdown && has_pandoc && has_DT && has_knitr

  if (!has_rmarkdown) message("Package 'rmarkdown' not available — skipping HTML report.")
  if (has_rmarkdown && !has_pandoc) message("pandoc not found — skipping HTML report. ",
    "On HPC, try: module load pandoc  or set RSTUDIO_PANDOC.")
  if (!has_DT)    message("Package 'DT' not available — skipping HTML report.")
  if (!has_knitr) message("Package 'knitr' not available — skipping HTML report.")

  if (render_html) {
    message("Generating HTML report...")

    # Ensure output dir exists and resolve to an absolute path — rmarkdown::render
    # changes the working directory, so relative paths break mid-render.
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    abs_output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

    # Save all results into a single RData file loaded by the Rmd
    report_data <- list(
      qc_results           = qc_results,
      dim_reduction        = dim_reduction,
      cnv_data             = cnv_data,
      reference_projection = reference_projection,
      sample_info          = sample_info,
      output_dirs          = output_dirs,
      generation_time      = Sys.time()
    )
    save(report_data, file = file.path(abs_output_dir, "report_data.RData"))

    # Write the Rmd template (absolute path so pandoc's chdir doesn't break it)
    rmd_file  <- create_report_rmd(abs_output_dir)
    html_file <- file.path(abs_output_dir, "methylation_analysis_report.html")

    # Render — catch errors so a rendering failure never kills the pipeline
    render_ok <- tryCatch({
      rmarkdown::render(
        input       = rmd_file,
        output_file = basename(html_file),
        output_dir  = abs_output_dir,
        intermediates_dir = abs_output_dir,
        knit_root_dir     = abs_output_dir,
        quiet       = FALSE,          # show render progress/errors
        envir       = new.env()       # clean environment for reproducibility
      )
      TRUE
    }, error = function(e) {
      message("HTML rendering failed: ", conditionMessage(e))
      message("  Rmd file kept at: ", rmd_file)
      message("  Falling back to plain-text report.")
      FALSE
    })

    if (render_ok) {
      message("HTML report saved: ", html_file)
      report_files$html <- html_file
    }
  }

  # Always produce a plain-text report as a reliable fallback
  if (is.null(report_files$html)) {
    text_report <- generate_text_report(qc_results, dim_reduction, cnv_data,
                                        sample_info, reference_projection)
    text_file   <- file.path(output_dir, "methylation_analysis_report.txt")
    writeLines(text_report, text_file)
    message("Plain-text report saved: ", text_file)
    report_files$text <- text_file
  }

  return(report_files)
}

#' Create R Markdown report template
#'
#' @param output_dir Output directory
#' @return Path to created Rmd file
create_report_rmd <- function(output_dir) {
  # Create the YAML header with proper date formatting
  yaml_header <- paste0('---
title: "MeQTrack Analysis Report"
author: "Methylation Pipeline"
date: "', format(Sys.time(), "%d %B, %Y"), '"
output:
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
    highlight: tango
---')

  # Create the R Markdown content body
  rmd_body <- '
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
                      fig.width = 10, fig.height = 7)
library(ggplot2)
library(DT)
library(knitr)

# Load report data saved by generate_report()
load("report_data.RData")

# ------------------------------------------------------------------
# Helper: resolve a stored plot path to one that actually exists.
# Tries the stored path first, then looks for the file by basename
# inside the known output module directories, then falls back to
# relative paths from the reports directory.
# ------------------------------------------------------------------
resolve_plot_path <- function(stored_path, module_dir = NULL) {
  if (is.null(stored_path) || is.na(stored_path) || !nzchar(stored_path)) return(NULL)

  # 1. Stored path already valid
  if (file.exists(stored_path)) return(stored_path)

  fname <- basename(stored_path)

  # 2. Look for the file in the supplied module directory
  if (!is.null(module_dir) && nzchar(module_dir)) {
    candidate <- file.path(module_dir, fname)
    if (file.exists(candidate)) return(candidate)
    # Also try one level deeper (e.g. qc/plots/)
    candidate2 <- file.path(module_dir, "plots", fname)
    if (file.exists(candidate2)) return(candidate2)
  }

  # 3. Look relative to the reports directory (where we are rendering)
  for (rel in c(
      file.path("..", "figures", "qc",            fname),
      file.path("..", "figures", "dim_reduction", fname),
      file.path("..", "figures", "cnv",           fname),
      file.path("..", "qc", "plots",              fname),
      file.path("..", "qc",                       fname),
      file.path("..", "dimensionality_reduction", fname),
      file.path("..", "cnv",                      fname),
      file.path("..", "cnv", "plots",             fname),
      file.path("..", "reference_projection",     fname),
      fname
  )) {
    if (file.exists(rel)) return(rel)
  }

  # 4. Return NULL so callers can show "not available" gracefully
  return(NULL)
}

# Convenience wrapper: include a static plot (PDF or PNG) or show a message
show_plot <- function(path, caption = "") {
  path <- resolve_plot_path(path)
  if (!is.null(path)) {
    knitr::include_graphics(path)
  } else {
    cat("*Plot not available*\n\n")
  }
}

# Convenience: make a datatable or print a message
show_table <- function(df, caption = "", ...) {
  if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
    DT::datatable(df, caption = caption,
                  options = list(pageLength = 10, scrollX = TRUE), ...)
  } else {
    cat("*Table not available*\n\n")
  }
}

# Pull directories from report_data
odirs <- report_data$output_dirs   # may be NULL if pipeline is old
```

# Overview

This report summarises the results of methylation array analysis performed using MeQTrack.

## Sample Information

```{r sample_info}
show_table(report_data$sample_info, caption = "Sample Information")
```

# Quality Control {.tabset}

```{r qc_check}
qc <- report_data$qc_results
has_qc <- !is.null(qc)
```

```{r qc_summary, eval=has_qc}
n_pass <- length(qc$passed_samples)
n_fail <- length(qc$failed_samples)
cat(sprintf("**%d samples passed QC** | **%d samples failed QC**\n\n", n_pass, n_fail))
```

## QC Metrics Table

```{r qc_table, eval=has_qc}
show_table(qc$sample_qc, caption = "Per-sample QC Metrics")
```

```{r qc_unavailable, eval=!has_qc}
cat("No QC results available.")
```

## Mean Detection P-value

```{r qc_detp, eval=has_qc}
p <- resolve_plot_path(qc$plots$mean_detection_pvalue,
                       if (!is.null(odirs)) odirs$figures_qc else NULL)
if (!is.null(p)) knitr::include_graphics(p) else cat("*Plot not available*\n\n")
```

## Beta Density

```{r qc_density, eval=has_qc}
p <- resolve_plot_path(qc$plots$beta_density,
                       if (!is.null(odirs)) odirs$figures_qc else NULL)
if (!is.null(p)) knitr::include_graphics(p) else cat("*Plot not available*\n\n")
```

## Beta Bean Plot

```{r qc_bean, eval=has_qc}
p <- resolve_plot_path(qc$plots$beta_bean,
                       if (!is.null(odirs)) odirs$figures_qc else NULL)
if (!is.null(p)) knitr::include_graphics(p) else cat("*Plot not available*\n\n")
```

## MDS Plot

```{r qc_mds, eval=has_qc}
p <- resolve_plot_path(qc$plots$mds,
                       if (!is.null(odirs)) odirs$figures_qc else NULL)
if (!is.null(p)) knitr::include_graphics(p) else cat("*Plot not available*\n\n")
```

## Interactive Plots

```{r qc_interactive, eval=has_qc}
idens <- resolve_plot_path(qc$plots$interactive_density,
                           if (!is.null(odirs)) odirs$figures_qc else NULL)
imds  <- resolve_plot_path(qc$plots$interactive_mds,
                           if (!is.null(odirs)) odirs$figures_qc else NULL)
if (!is.null(idens)) cat(sprintf("[Interactive density plot](%s)\n\n", basename(idens)))
if (!is.null(imds))  cat(sprintf("[Interactive MDS plot](%s)\n\n",     basename(imds)))
if (is.null(idens) && is.null(imds)) cat("*Interactive plots not available*\n\n")
```

# Dimensionality Reduction {.tabset}

```{r dr_check}
dr     <- report_data$dim_reduction
has_dr <- !is.null(dr)
dr_dir <- if (!is.null(odirs)) odirs$figures_dim_reduction else NULL
```

## t-SNE

```{r tsne_results, eval=has_dr && !is.null(dr$tsne)}
tsne <- dr$tsne
# coords field (actual name from dim_reduction.R)
coords_df <- tsne$coords
if (!is.null(coords_df) && is.data.frame(coords_df)) {
  show_table(coords_df, caption = "t-SNE Coordinates")
}

p <- resolve_plot_path(file.path(if (!is.null(dr_dir)) dr_dir else ".", "tsne_plot.pdf"),
                       dr_dir)
if (!is.null(p)) knitr::include_graphics(p) else cat("*t-SNE plot not available*\n\n")

if (length(tsne$duplicates) > 0) {
  cat("\n**Duplicate samples removed before t-SNE:**",
      paste(tsne$duplicates, collapse = ", "), "\n\n")
}
```

```{r tsne_unavailable, eval=!has_dr || is.null(dr$tsne)}
cat("No t-SNE results available.")
```

## UMAP

```{r umap_results, eval=has_dr && !is.null(dr$umap)}
umap_res <- dr$umap
coords_df <- umap_res$coords
if (!is.null(coords_df) && is.data.frame(coords_df)) {
  show_table(coords_df, caption = "UMAP Coordinates")
}

p <- resolve_plot_path(file.path(if (!is.null(dr_dir)) dr_dir else ".", "umap_plot.pdf"),
                       dr_dir)
if (!is.null(p)) knitr::include_graphics(p) else cat("*UMAP plot not available*\n\n")
```

```{r umap_unavailable, eval=!has_dr || is.null(dr$umap)}
cat("No UMAP results available.")
```

## Hierarchical Clustering

```{r hclust_results, eval=has_dr && !is.null(dr$hclust)}
hcl <- dr$hclust
cat(sprintf("**Method:** %s | **Distance:** %s\n\n",
            hcl$method, hcl$distance))

p <- resolve_plot_path(file.path(if (!is.null(dr_dir)) dr_dir else ".", "hclust_dendrogram.pdf"),
                       dr_dir)
if (!is.null(p)) knitr::include_graphics(p) else cat("*Dendrogram not available*\n\n")
```

```{r hclust_unavailable, eval=!has_dr || is.null(dr$hclust)}
cat("No hierarchical clustering results available.")
```

# Reference Projection {.tabset}

```{r refproj_check}
rp     <- report_data$reference_projection
has_rp <- !is.null(rp)
rp_dir <- if (!is.null(odirs)) odirs$reference_projection else NULL
```

```{r refproj_summary, eval=has_rp}
rp_name <- if (!is.null(rp$label)) rp$label else rp$dataset
cat(sprintf("Query samples projected onto the **%s** reference (%s reference samples). Each sample is assigned its nearest reference tumour group by a k-NN vote in the embedding.\n\n",
            rp_name,
            if (!is.null(rp$ref_meta)) nrow(rp$ref_meta) else "?"))
```

## Class Hints

```{r refproj_table, eval=has_rp}
proj <- as.data.frame(rp$projected)
ch   <- rp$class_hints
if (!is.null(ch)) {
  ch  <- ch[match(proj$Sample, ch$Sample), , drop = FALSE]
  tbl <- data.frame(
    Sample          = proj$Sample,
    tSNE1           = round(proj$tSNE1, 2),
    tSNE2           = round(proj$tSNE2, 2),
    `Nearest class` = ch$nearest_class,
    Confidence      = sprintf("%.0f%%", 100 * ch$confidence),
    `Top classes`   = ch$top_classes,
    Ambiguous       = ifelse(ch$ambiguous %in% TRUE, "yes", ""),
    Distant         = ifelse(ch$distant_from_reference %in% TRUE, "yes", ""),
    check.names     = FALSE
  )
} else {
  tbl <- data.frame(Sample = proj$Sample,
                    tSNE1  = round(proj$tSNE1, 2),
                    tSNE2  = round(proj$tSNE2, 2),
                    check.names = FALSE)
}
show_table(tbl, caption = "Per-sample nearest reference tumour class")
```

```{r refproj_unavailable, eval=!has_rp}
cat("No reference projection results available.")
```

## Projection Plot

```{r refproj_plot, eval=has_rp}
p <- resolve_plot_path(rp$pdf, rp_dir)
if (!is.null(p)) knitr::include_graphics(p) else cat("*Projection plot not available*\n\n")
```

# Copy Number Variation Analysis {.tabset}

```{r cnv_check}
cnv_data <- report_data$cnv_data
has_cnv  <- !is.null(cnv_data)
cnv_dir  <- if (!is.null(odirs)) odirs$figures_cnv else NULL
```

## CNV Segments

```{r cnv_segments, eval=has_cnv}
cat(sprintf("CNV analysis method: **%s**\n\n",
            if (!is.null(cnv_data$method)) cnv_data$method else "unknown"))
show_table(cnv_data$segments, caption = "CNV Segments")
```

```{r cnv_unavailable, eval=!has_cnv}
cat("No CNV analysis results available.")
```

## Frequency Plot

```{r cnv_freq, eval=has_cnv}
p <- resolve_plot_path(cnv_data$frequency_plot, cnv_dir)
if (!is.null(p)) knitr::include_graphics(p) else cat("*CNV frequency plot not available*\n\n")
```

## Individual Sample Profiles

```{r cnv_samples, eval=has_cnv && !is.null(cnv_data$sample_results)}
plots_found <- 0L
for (res in cnv_data$sample_results) {
  if (!is.null(res$plot_file)) {
    p <- resolve_plot_path(res$plot_file, cnv_dir)
    if (!is.null(p)) {
      cat("### Sample:", res$sample_id, "\n\n")
      print(knitr::include_graphics(p))
      cat("\n\n")
      plots_found <- plots_found + 1L
    }
  }
}
if (plots_found == 0L) {
  cat(sprintf("*No individual CNV plots found (expected %d samples)*\n\n",
              length(cnv_data$sample_results)))
}
```

# Session Information

```{r session_info}
cat("Report generated on:", format(report_data$generation_time, "%Y-%m-%d %H:%M:%S"), "\n\n")
sessionInfo()
```
'

  # Combine header and body
  rmd_content <- paste0(yaml_header, rmd_body)

  # Write the Rmd file
  rmd_file <- file.path(output_dir, "methylation_analysis_report.Rmd")
  writeLines(rmd_content, rmd_file)

  return(rmd_file)
}

#' Generate simple text report
#'
#' @param qc_results QC results
#' @param dim_reduction Dimensionality reduction results
#' @param cnv_data CNV analysis results
#' @param sample_info Sample information data frame
#' @param reference_projection Reference-projection results (rp_result), optional
#' @return Character vector with report text
generate_text_report <- function(qc_results, dim_reduction, cnv_data,
                                 sample_info, reference_projection = NULL) {
  report_lines <- c(
    "=================================",
    "Methylation Array Analysis Report",
    "=================================",
    paste("Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Overview",
    "--------",
    paste("Number of samples:", if(!is.null(sample_info)) nrow(sample_info) else "Not available"),
    ""
  )

  # QC summary
  if (!is.null(qc_results)) {
    report_lines <- c(report_lines,
      "Quality Control Summary",
      "----------------------",
      paste("Total samples processed:", length(c(qc_results$passed_samples, qc_results$failed_samples))),
      paste("Samples passing QC:", length(qc_results$passed_samples)),
      paste("Samples failing QC:", length(qc_results$failed_samples)),
      ""
    )
  }

  # Dimensionality reduction summary
  if (!is.null(dim_reduction)) {
    report_lines <- c(report_lines,
      "Dimensionality Reduction Summary",
      "-------------------------------"
    )

    if (!is.null(dim_reduction$tsne)) {
      n_samp <- if (!is.null(dim_reduction$tsne$coords)) nrow(dim_reduction$tsne$coords) else "unknown"
      report_lines <- c(report_lines,
        paste("t-SNE: samples =", n_samp),
        ""
      )
    }

    if (!is.null(dim_reduction$umap)) {
      n_samp <- if (!is.null(dim_reduction$umap$coords)) nrow(dim_reduction$umap$coords) else "unknown"
      report_lines <- c(report_lines,
        paste("UMAP: samples =", n_samp),
        ""
      )
    }

    if (!is.null(dim_reduction$hclust)) {
      report_lines <- c(report_lines,
        paste("Hierarchical clustering: method =", dim_reduction$hclust$method,
              "| distance =", dim_reduction$hclust$distance),
        ""
      )
    }
  }

  # CNV analysis summary
  if (!is.null(cnv_data)) {
    report_lines <- c(report_lines,
      "Copy Number Variation Analysis Summary",
      "------------------------------------",
      paste("Method:", if (!is.null(cnv_data$method)) cnv_data$method else "unknown"),
      paste("Number of segments:", if (!is.null(cnv_data$segments)) nrow(cnv_data$segments) else "unknown"),
      ""
    )
  }

  # Reference projection summary
  if (!is.null(reference_projection) &&
      !is.null(reference_projection$class_hints)) {
    ch <- reference_projection$class_hints
    report_lines <- c(report_lines,
      "Reference Projection Summary",
      "---------------------------",
      paste("Reference dataset:",
            if (!is.null(reference_projection$label))
              reference_projection$label else reference_projection$dataset),
      paste("Samples projected:", nrow(ch)),
      ""
    )
    for (i in seq_len(nrow(ch))) {
      report_lines <- c(report_lines,
        sprintf("  %s -> %s (%.0f%% confidence)%s",
                ch$Sample[i], ch$nearest_class[i], 100 * ch$confidence[i],
                if (isTRUE(ch$ambiguous[i])) " [ambiguous]" else ""))
    }
    report_lines <- c(report_lines, "")
  }

  return(report_lines)
}

#' Create interactive plots for various analysis results
#'
#' @param results Analysis results
#' @param output_dir Output directory for plots
#' @return List of interactive plot paths
create_interactive_plots <- function(results, output_dir = ".") {
  plots <- list()

  # Check if required packages are available
  has_plotly <- requireNamespace("plotly", quietly = TRUE)
  has_htmlwidgets <- requireNamespace("htmlwidgets", quietly = TRUE)

  if (!has_plotly || !has_htmlwidgets) {
    message("Packages 'plotly' and 'htmlwidgets' required for interactive plots.")
    return(plots)
  }

  return(plots)
}
