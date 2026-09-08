# app/R/samplesheet_module.R
# ---------------------------------------------------------------------------
# Shiny module for samplesheet picking, validation, array-type detection,
# and optional-metadata preview.
#
# Public exports:
#   samplesheet_ui(id)               UI pieces: sidebar inputs + main panel
#   samplesheet_server(id, workspace, project_root)
#     -> reactive list(valid, array_type, samplesheet_path, validated_df,
#                       optional_cols)
#
# The module does NOT talk to the pipeline. It is strictly a pre-flight
# input stage. The Run controller (Wave 3) consumes its outputs via the
# returned reactive.
# ---------------------------------------------------------------------------

# This file is sourced from app.R, which has already sourced validators.R
# and array_detect.R, so we don't re-source here.

#' UI pieces for the samplesheet module.
#'
#' Returns a list with three named slots that app.R places in the layout:
#'   $sidebar        controls (file picker, array-type dropdown)
#'   $main           validation table + summary
#'   $metadata       optional-metadata preview card
samplesheet_ui <- function(id) {
  ns <- shiny::NS(id)
  list(
    sidebar = shiny::tagList(
      shiny::h5("Samplesheet"),
      shiny::helpText(
        "Pick a CSV samplesheet from your workspace. ",
        "Required columns: Sentrix_ID, Sample_Name, Basename."
      ),
      shinyFiles::shinyFilesButton(
        ns("pick_samplesheet"), label = "Choose samplesheet...",
        title = "Select a samplesheet (.csv)", multiple = FALSE,
        icon = shiny::icon("folder-open")
      ),
      shiny::br(), shiny::br(),
      shiny::uiOutput(ns("selected_path")),
      shiny::hr(),
      shiny::h5("Array type"),
      shiny::selectInput(
        ns("array_type"),
        label = NULL,
        choices = c("auto", "450K", "EPIC", "EPICv2"),
        selected = "auto"
      ),
      shiny::uiOutput(ns("array_detect_msg"))
    ),
    main = shiny::tagList(
      shiny::uiOutput(ns("top_status")),
      shiny::uiOutput(ns("summary_panel")),
      DT::DTOutput(ns("rows_table"))
    ),
    metadata = shiny::tagList(
      shiny::h5("Optional metadata columns"),
      shiny::uiOutput(ns("metadata_preview"))
    )
  )
}

#' Server logic for the samplesheet module.
#'
#' @param id               module namespace id
#' @param workspace        reactive expression returning the workspace path
#' @param project_root_    function/path returning the project root dir
#'                         (used to build shinyFiles roots)
#' @return reactive list with the validated state the run controller will
#'         consume in Wave 3.
samplesheet_server <- function(id, workspace, project_root_) {
  shiny::moduleServer(id, function(input, output, session) {

    # Roots for shinyFiles: labeled shortcuts for the common locations,
    # plus every mounted volume so the user can browse anywhere on the
    # machine (external drives, network mounts, /, etc.).
    roots_fn <- function() {
      ws <- if (is.function(workspace)) workspace() else workspace
      pr <- if (is.function(project_root_)) project_root_() else project_root_
      c(
        Workspace = ws,
        `MeQTrack App` = pr,
        Home = path.expand("~"),
        shinyFiles::getVolumes()()
      )
    }

    shinyFiles::shinyFileChoose(
      input, "pick_samplesheet",
      roots = roots_fn,
      filetypes = c("csv"),
      defaultRoot = "Workspace"
    )

    # Reactive: selected samplesheet path (absolute).
    selected_path <- shiny::reactive({
      f <- input$pick_samplesheet
      if (is.null(f) || length(f) == 0L || is.integer(f)) return(NULL)
      parsed <- shinyFiles::parseFilePaths(roots = roots_fn(), f)
      if (nrow(parsed) == 0L) return(NULL)
      normalizePath(as.character(parsed$datapath[1]),
                    winslash = "/", mustWork = FALSE)
    })

    output$selected_path <- shiny::renderUI({
      p <- selected_path()
      if (is.null(p)) {
        shiny::tags$em("No samplesheet selected.")
      } else {
        shiny::tagList(
          shiny::tags$strong("Selected: "),
          shiny::tags$code(p)
        )
      }
    })

    # Validation: re-runs whenever selected_path changes.
    validated <- shiny::reactive({
      p <- selected_path()
      if (is.null(p)) return(NULL)
      parsed <- read_samplesheet(p)
      if (!parsed$ok) {
        return(list(top_error = parsed$error, df = NULL,
                    samplesheet_dir = NULL))
      }
      df <- validate_rows(parsed$df,
                          samplesheet_dir = parsed$samplesheet_dir,
                          pipeline_dir    = pipeline_dir())
      list(top_error = NULL, df = df, samplesheet_dir = parsed$samplesheet_dir)
    })

    # Top-level error banner (parse or missing-column failures).
    output$top_status <- shiny::renderUI({
      v <- validated()
      if (is.null(v) || is.null(v$top_error)) return(NULL)
      shiny::div(class = "alert alert-danger", role = "alert",
                 shiny::tags$strong("Samplesheet error: "), v$top_error)
    })

    # Summary above the table.
    output$summary_panel <- shiny::renderUI({
      v <- validated()
      if (is.null(v) || is.null(v$df)) return(NULL)
      s <- validation_summary(v$df)
      cls <- if (s$ok) "alert alert-success" else "alert alert-warning"
      msg <- if (s$ok) {
        sprintf("All %d row(s) OK. Ready to run.", s$n_total)
      } else {
        sprintf("%d of %d row(s) OK. %d need attention before running.",
                s$n_ok, s$n_total, s$n_bad)
      }
      shiny::div(class = cls, role = "alert", msg)
    })

    # Per-row validation table.
    output$rows_table <- DT::renderDT({
      v <- validated()
      if (is.null(v) || is.null(v$df)) return(NULL)
      display <- v$df[, c("Sentrix_ID", "Sample_Name", "Basename",
                          "ResolvedBasename", "Status", "StatusDetail")]
      DT::datatable(
        display,
        rownames = FALSE,
        selection = "none",
        options = list(pageLength = 25, dom = "tip", autoWidth = TRUE),
        class = "stripe hover compact"
      ) |>
        DT::formatStyle(
          "Status",
          target = "row",
          backgroundColor = DT::styleEqual(
            c(STATUS_OK,
              STATUS_MISSING_RED, STATUS_MISSING_GRN, STATUS_MISSING_BOTH,
              STATUS_DUPLICATE_ID, STATUS_MALFORMED),
            c("#ffffff",
              "#fff3cd", "#fff3cd", "#f8d7da",
              "#f8d7da", "#f8d7da")
          )
        )
    })

    # Array-type detection (runs after rows pass basic validation).
    detection <- shiny::reactive({
      v <- validated()
      if (is.null(v) || is.null(v$df)) return(NULL)
      detect_array_from_samplesheet(v$df)
    })

    output$array_detect_msg <- shiny::renderUI({
      d <- detection()
      if (is.null(d)) return(shiny::tags$em("Load a samplesheet first."))
      if (!is.na(d$array_type)) {
        shiny::tagList(
          shiny::tags$small(
            shiny::tags$strong("Detected: "),
            shiny::tags$code(d$array_type),
            " — ", d$reason
          )
        )
      } else {
        shiny::tags$small(
          shiny::tags$em("Could not auto-detect. "),
          d$reason,
          " Use the dropdown to set array type manually."
        )
      }
    })

    # Optional-metadata preview.
    output$metadata_preview <- shiny::renderUI({
      v <- validated()
      if (is.null(v) || is.null(v$df)) {
        return(shiny::tags$em("Load a samplesheet to see optional metadata."))
      }
      opt <- optional_metadata_columns(v$df)
      if (!length(opt)) {
        return(shiny::tags$em("No optional metadata columns detected."))
      }
      items <- lapply(opt, function(col) {
        vals <- v$df[[col]]
        uniq <- unique(vals[!is.na(vals)])
        preview <- if (length(uniq) <= 6) {
          paste(uniq, collapse = ", ")
        } else {
          sprintf("%s ... (%d unique values)",
                  paste(utils::head(uniq, 6), collapse = ", "),
                  length(uniq))
        }
        shiny::tags$li(
          shiny::tags$strong(col), ": ", preview
        )
      })
      shiny::tags$ul(items)
    })

    # Resolved array type: dropdown overrides detection unless "auto".
    resolved_array_type <- shiny::reactive({
      chosen <- input$array_type %||% "auto"
      if (chosen != "auto") return(chosen)
      d <- detection()
      if (is.null(d) || is.na(d$array_type)) "auto" else d$array_type
    })

    # Public reactive. This is what the run controller will read in Wave 3.
    shiny::reactive({
      v <- validated()
      s <- if (is.null(v) || is.null(v$df)) {
        list(ok = FALSE, n_total = 0L, n_ok = 0L, n_bad = 0L)
      } else {
        validation_summary(v$df)
      }
      list(
        valid             = isTRUE(s$ok),
        array_type        = resolved_array_type(),
        samplesheet_path  = selected_path(),
        samplesheet_dir   = v$samplesheet_dir,
        validated_df      = v$df,
        optional_cols     = if (is.null(v)) character(0)
                            else optional_metadata_columns(v$df)
      )
    })
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
