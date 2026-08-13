# app/R/past_runs_module.R
# ---------------------------------------------------------------------------
# Past runs library — Wave 6 Theme 6b.
#
# Lists every <workspace>/runs/<id>/ directory as a row in a DT table and
# lets the user "Open" one to attach it as the current session's run_dir.
# Once attached, the result modules (QC / Dim. reduction / CNV / Report)
# render the past run's artifacts via the existing results_loader pipeline,
# and the per-step Run buttons in run_controller work against it.
#
# Public exports:
#   past_runs_module_ui(id)
#   past_runs_module_server(id, workspace) -> reactive selected_run_dir | NULL
# ---------------------------------------------------------------------------

PAST_RUNS_COL_TOOLTIPS <- list(
  `Run ID`     = "Timestamped folder name under <workspace>/runs/. Format: YYYYMMDD-HHMMSS_<samplesheet>.",
  Started      = "When the run was launched (local time).",
  Samplesheet  = "Name of the samplesheet CSV the run was launched against.",
  Array        = "Methylation array type detected from the IDATs (450K / EPIC / EPICv2).",
  `Last step`  = "The pipeline step requested for this run. 'all' = full pipeline, otherwise the single per-step name.",
  Status       = "Final state of the run: completed (exit 0), failed (non-zero exit), running (still active), cancelled (user stopped it).",
  Exit         = "Process exit code. 0 means success; 137 typically means the OS killed the process for being out of memory.",
  `# samples`  = "Sample count derived from the QC report (preferred) or processed_data/sample_info.txt. NA if neither file exists yet.",
  Action       = "Click Open to attach this run as the active session — result tabs render its artifacts and per-step buttons start operating against it."
)

past_runs_module_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "d-flex justify-content-between align-items-center mb-2",
      shiny::tags$small(class = "text-muted",
        "Click 'Open' to load a past run into the result tabs. ",
        "Per-step Run buttons will then work against it."),
      shiny::actionButton(ns("refresh"), "Refresh",
                          icon = shiny::icon("rotate"),
                          class = "btn-sm btn-outline-secondary")
    ),
    shiny::uiOutput(ns("attached_banner")),
    shinycssloaders::withSpinner(
      DT::DTOutput(ns("runs_table")),
      type = 7, color = COLORS$primary
    )
  )
}

past_runs_module_server <- function(id, workspace) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive trigger so we can refresh the listing on demand.
    refresh_trigger <- shiny::reactiveVal(0L)
    shiny::observeEvent(input$refresh, {
      refresh_trigger(refresh_trigger() + 1L)
    })

    runs_df <- shiny::reactive({
      refresh_trigger()
      ws <- if (is.function(workspace)) workspace() else workspace
      list_runs(file.path(ws, "runs"))
    })

    output$runs_table <- DT::renderDT({
      df <- runs_df()
      if (is.null(df) || !nrow(df)) {
        return(DT::datatable(
          data.frame(`Past runs` = character(0), check.names = FALSE),
          options = list(language = list(emptyTable = "No past runs found.")),
          rownames = FALSE
        ))
      }

      # One Open button per row. The handler dispatches via input$open_row,
      # populated by the per-button click via Shiny.setInputValue.
      df$Action <- vapply(seq_len(nrow(df)), function(i) {
        sprintf(
          '<button class="btn btn-sm btn-primary open-run-btn" data-row="%d">Open</button>',
          i
        )
      }, character(1))

      display <- df[, c("run_id", "started_at", "samplesheet_name",
                        "array_type", "last_step", "status",
                        "exit_code", "n_samples", "Action"), drop = FALSE]
      colnames(display) <- c("Run ID", "Started", "Samplesheet",
                             "Array", "Last step", "Status",
                             "Exit", "# samples", "Action")

      action_col <- which(colnames(display) == "Action")

      DT::datatable(
        display,
        rownames = FALSE,
        escape = setdiff(seq_len(ncol(display)), action_col),
        selection = "none",
        class = "stripe hover compact",
        options = list(
          pageLength = 25,
          order = list(),
          autoWidth = FALSE,
          dom = "ltip",
          headerCallback = dt_header_tooltips(PAST_RUNS_COL_TOOLTIPS),
          columnDefs = list(
            list(className = "dt-right",
                 targets = which(colnames(display) %in%
                                   c("Exit", "# samples")) - 1L)
          ),
          drawCallback = DT::JS(sprintf(
            "function(settings) {
              var ns = '%s';
              $(this).find('.open-run-btn').off('click').on('click', function() {
                var row = parseInt($(this).data('row'), 10);
                Shiny.setInputValue(ns + 'open_row',
                  {row: row, _nonce: Math.random()},
                  {priority: 'event'});
              });
            }",
            ns("")
          ))
        )
      ) |>
        DT::formatStyle(
          "Status",
          color = DT::styleEqual(
            c("completed", "failed", "cancelled", "running", "(no manifest)"),
            c("#198754",  "#dc3545", "#fd7e14",   "#0d6efd", "#6c757d")
          ),
          fontWeight = "bold"
        )
    }, server = FALSE)

    # Click-to-open: returns the run_dir path for the selected row.
    selected_run_dir <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$open_row, {
      df <- runs_df()
      idx <- as.integer(input$open_row$row)
      if (is.na(idx) || idx < 1L || idx > nrow(df)) return()
      path <- df$run_dir[idx]
      if (!dir.exists(path)) {
        shiny::showNotification(
          sprintf("Run directory no longer exists: %s", path),
          type = "warning"
        )
        return()
      }
      selected_run_dir(path)
    })

    output$attached_banner <- shiny::renderUI({
      path <- selected_run_dir()
      if (is.null(path)) return(NULL)
      shiny::div(
        class = "alert alert-info py-2 mb-2",
        shiny::icon("link"),
        shiny::tags$strong(" Attached: "),
        shiny::tags$code(basename(path)),
        shiny::tags$small(class = "ms-2 text-muted",
          "Result tabs and per-step buttons now operate on this run.")
      )
    })

    selected_run_dir
  })
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Scan a workspace's runs/ directory and return one row per subdirectory,
# parsed from run_manifest.json when present. Sorted newest-first.
list_runs <- function(runs_dir) {
  if (is.null(runs_dir) || !dir.exists(runs_dir)) {
    return(empty_runs_df())
  }
  subdirs <- list.dirs(runs_dir, recursive = FALSE, full.names = TRUE)
  if (!length(subdirs)) return(empty_runs_df())

  rows <- lapply(subdirs, parse_run_dir)
  df <- do.call(rbind, rows)
  if (is.null(df) || !nrow(df)) return(empty_runs_df())

  ord <- order(df$sort_key, decreasing = TRUE)
  df[ord, , drop = FALSE]
}

empty_runs_df <- function() {
  data.frame(
    run_dir = character(0),
    run_id = character(0),
    started_at = character(0),
    samplesheet_name = character(0),
    array_type = character(0),
    last_step = character(0),
    status = character(0),
    exit_code = character(0),
    n_samples = integer(0),
    sort_key = character(0),
    stringsAsFactors = FALSE
  )
}

parse_run_dir <- function(path) {
  manifest_path <- file.path(path, "run_manifest.json")
  m <- if (file.exists(manifest_path)) {
    tryCatch(jsonlite::read_json(manifest_path, simplifyVector = TRUE),
             error = function(e) NULL)
  } else NULL

  run_id           <- if (!is.null(m$run_id))      m$run_id      else basename(path)
  started_raw      <- if (!is.null(m$started_at))  m$started_at  else
                        format(file.info(path)$mtime, "%Y-%m-%dT%H:%M:%S")
  samplesheet_name <- if (!is.null(m$samplesheet)) basename(m$samplesheet) else "—"
  array_type       <- if (!is.null(m$array_type))  m$array_type  else "—"
  last_step        <- if (!is.null(m$step))        m$step        else "—"
  status           <- if (!is.null(m$status))      m$status      else "(no manifest)"
  exit_code <- if (!is.null(m$exit_code) && length(m$exit_code) == 1L &&
                   !is.na(m$exit_code)) {
    as.character(m$exit_code)
  } else "—"

  data.frame(
    run_dir = path,
    run_id = run_id,
    started_at = pretty_time(started_raw),
    samplesheet_name = samplesheet_name,
    array_type = array_type,
    last_step = last_step,
    status = status,
    exit_code = exit_code,
    n_samples = count_samples(path),
    sort_key = started_raw,
    stringsAsFactors = FALSE
  )
}

count_samples <- function(run_dir) {
  qc_csv <- file.path(run_dir, "qc", "sample_qc_report.csv")
  if (file.exists(qc_csv)) {
    n <- tryCatch(nrow(utils::read.csv(qc_csv)), error = function(e) NA_integer_)
    if (!is.na(n)) return(n)
  }
  si <- file.path(run_dir, "processed_data", "sample_info.txt")
  if (file.exists(si)) {
    n <- tryCatch(
      nrow(utils::read.table(si, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)),
      error = function(e) NA_integer_
    )
    if (!is.na(n)) return(n)
  }
  NA_integer_
}

pretty_time <- function(ts) {
  if (is.null(ts) || is.na(ts) || !nzchar(ts)) return("—")
  parsed <- tryCatch(
    as.POSIXct(ts, format = "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
    error = function(e) NA
  )
  if (is.na(parsed)) {
    parsed <- tryCatch(
      as.POSIXct(ts, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
      error = function(e) NA
    )
  }
  if (is.na(parsed)) return(ts)
  format(parsed, "%Y-%m-%d %H:%M", tz = Sys.timezone())
}
