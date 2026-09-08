# app/R/report_module.R
# ---------------------------------------------------------------------------
# Report view: in-app preview of the pipeline's self-contained HTML report,
# plus the "Open in browser" and "Show in Finder/Explorer" actions.
#
# The same actions exist on the Run tab (they come for free from the run
# controller). They're duplicated here because the Report tab is where a
# returning user would intuitively look for them once a run is in the past.
# ---------------------------------------------------------------------------

report_module_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "d-flex gap-2 mb-3",
      shiny::uiOutput(ns("actions"), inline = TRUE)
    ),
    shiny::uiOutput(ns("report_frame"))
  )
}

report_module_server <- function(id, results) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    report_rel_path <- function(r) {
      # Match discover_report in run_controller.R — check both the spec'd
      # layout and the one the pipeline actually produces.
      candidates <- c(
        "reports/methylation_analysis_report.html",
        "report/methylation_analysis_report.html"
      )
      for (c in candidates) {
        if (file.exists(file.path(r$run_dir, c))) return(c)
      }
      # Fall back to any .html under reports/ or report/.
      for (sub in c("reports", "report")) {
        d <- file.path(r$run_dir, sub)
        if (!dir.exists(d)) next
        fs <- list.files(d, pattern = "\\.html$", full.names = FALSE)
        if (length(fs)) return(file.path(sub, fs[1]))
      }
      NULL
    }

    output$report_frame <- shiny::renderUI({
      r <- results()
      if (is.null(r)) {
        return(shiny::div(
          class = "alert alert-secondary",
          "No completed run yet. Run the pipeline first; the report will ",
          "appear here automatically."
        ))
      }
      rel <- report_rel_path(r)
      if (is.null(rel)) {
        return(shiny::div(
          class = "alert alert-warning",
          "Run completed but the HTML report was not found under ",
          shiny::tags$code("reports/"), "."
        ))
      }
      shiny::tags$iframe(
        src = paste0(r$run_url_base, "/", rel),
        style = "width: 100%; height: 82vh; border: 1px solid #dee2e6; border-radius: 4px;"
      )
    })

    output$actions <- shiny::renderUI({
      r <- results()
      if (is.null(r)) return(NULL)
      rel <- report_rel_path(r)
      if (is.null(rel)) return(NULL)
      shiny::tagList(
        shiny::actionButton(
          ns("open_new_tab"), "Open report in new tab",
          icon = shiny::icon("arrow-up-right-from-square"),
          class = "btn-primary"
        ),
        shiny::actionButton(
          ns("reveal"), reveal_label(),
          icon = shiny::icon("folder-open"),
          class = "btn-outline-secondary"
        )
      )
    })

    shiny::observeEvent(input$open_new_tab, {
      r <- results()
      if (is.null(r)) return()
      rel <- report_rel_path(r)
      if (is.null(rel)) return()
      disk <- file.path(r$run_dir, rel)
      utils::browseURL(normalizePath(disk, winslash = "/"))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$reveal, {
      r <- results()
      if (is.null(r)) return()
      rel <- report_rel_path(r)
      if (is.null(rel)) return()
      disk <- file.path(r$run_dir, rel)
      reveal_in_file_manager(disk)
    }, ignoreInit = TRUE)
  })
}
