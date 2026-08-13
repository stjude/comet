# app/R/run_controller.R
# ---------------------------------------------------------------------------
# Shiny module that owns pipeline execution: Run/Cancel controls, stage
# progress tracking, live log tail, and post-run actions (open report,
# reveal in file manager).
#
# Wave 3 module. Consumes the samplesheet module's reactive state and talks
# to the pipeline through pipeline_bridge.R (already written).
#
# Design note (deviation from the original design spec):
#   the spec suggested ExtendedTask + promises + future. We use callr::r_bg
#   directly (a real OS subprocess) and poll its state via invalidateLater.
#   ExtendedTask is designed to move a synchronous R computation to a future
#   worker so it doesn't block the session; since callr::r_bg is already a
#   separate OS process the session is never blocked, and ExtendedTask adds
#   scaffolding without a payoff. If we ever need the promise-resolution
#   model (chaining then()/catch() on the run outcome), we can wrap the
#   handle in an ExtendedTask without changing the UI.
#
# Public exports:
#   run_controller_ui(id)                list of UI fragments
#                                        (controls / log / actions / header)
#   run_controller_server(id, ss_state, workspace, project_root_)
#     -> reactive list(state, stage, run_dir, exit_code)   (for header badge)
# ---------------------------------------------------------------------------

# User-facing pipeline stages. Probe filtering is folded into "qc" so users
# see one combined step instead of QC and filtering as separate rows; the
# pipeline still runs them as Step 2 + Step 3 internally, but they map onto
# the same UI row (see PIPELINE_STEP_TO_STAGE below).
RUN_STAGES <- c(
  preprocess           = "Preprocessing",
  qc                   = "QC and probe filtering",
  dim_reduction        = "Dimensionality reduction",
  reference_projection = "Reference projection",
  cnv                  = "Copy-number variation",
  visualization        = "Report generation"
)

# The pipeline emits "Step N:" log lines for the five user-facing steps.
# Probe filtering is a sub-step of QC and its log line deliberately omits
# the "Step N:" prefix so this regex does not match it.
PIPELINE_STEP_TO_STAGE <- c(
  "preprocess",     # Step 1
  "qc",             # Step 2 — quality control + probe filtering
  "dim_reduction",  # Step 3
  "cnv",            # Step 4
  "visualization"   # Step 5
)

# Regex that matches the pipeline's stage-start log line. The pipeline emits
# lines like: "[2026-04-22 14:03:11] Step 3: Dimensionality reduction analysis".
STAGE_LINE_RE <- "Step\\s+([1-5]):"

# Overall run states.
RUN_STATE_IDLE      <- "idle"
RUN_STATE_RUNNING   <- "running"
RUN_STATE_COMPLETED <- "completed"
RUN_STATE_FAILED    <- "failed"
RUN_STATE_CANCELLED <- "cancelled"

# Sentinel for the full-pipeline run.
RUN_STEP_ALL <- "all"

# Upstream artifact requirements for each per-step run. Each element is a
# list of "requirement groups" (paths relative to run_dir). A group is
# satisfied if ANY of its files exists; the step is runnable if ALL groups
# are satisfied. Preprocess and reference projection have no upstream
# prereqs — both read the IDATs directly, so their buttons are always
# enabled (when no run is in progress).
STEP_PREREQS <- list(
  preprocess           = list(),
  qc                   = list("processed_data/preprocessed_data.RData"),
  dim_reduction        = list(c("processed_data/preprocessed_data.RData",
                                "processed_data/filtered_beta_values.txt")),
  reference_projection = list(),
  cnv                  = list("processed_data/preprocessed_data.RData"),
  visualization        = list(c("qc/qc_results.RData",
                                "dimensionality_reduction/dim_reduction_results.RData",
                                "cnv/cnv_results.RData"))
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
run_controller_ui <- function(id) {
  ns <- shiny::NS(id)
  list(
    # Controls + stage progress — meant to sit beside the Settings card in
    # the Run tab (app.R composes the two-column layout).
    controls = shiny::tagList(
      shiny::div(
        class = "d-flex gap-2 align-items-center mb-3",
        shiny::actionButton(
          ns("run"), label = "Run analysis",
          icon = shiny::icon("play"),
          class = "btn-primary"
        ),
        shiny::actionButton(
          ns("cancel"), label = "Cancel",
          icon = shiny::icon("stop"),
          class = "btn-outline-danger"
        ),
        shiny::uiOutput(ns("state_badge"), inline = TRUE)
      ),
      shiny::uiOutput(ns("prereq_msg")),
      bslib::card(
        bslib::card_header("Stages"),
        bslib::card_body(
          shiny::uiOutput(ns("stages_panel"))
        )
      )
    ),
    # Live log — full-width card below the controls/settings row.
    log = bslib::card(
      bslib::card_header("Live log (last 20 lines)"),
      bslib::card_body(
        shiny::verbatimTextOutput(ns("log_tail"), placeholder = TRUE)
      )
    ),
    # Post-run actions — full-width card below the log.
    actions = bslib::card(
      bslib::card_header("Actions"),
      bslib::card_body(
        shiny::uiOutput(ns("post_run_actions"))
      )
    ),
    # Compact header slot — state badge in the navbar.
    header = shiny::uiOutput(ns("header_badge"), inline = TRUE)
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
run_controller_server <- function(id, ss_state, workspace, project_root_,
                                  attach_run = shiny::reactive(NULL),
                                  parameters = shiny::reactive(list())) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Core reactive state. We keep one authoritative reactiveValues bag so
    # the UI renderers only depend on this, not on the bridge handle shape.
    rv <- shiny::reactiveValues(
      state       = RUN_STATE_IDLE,
      handle      = NULL,                  # pipeline_bridge handle, or NULL
      run_dir     = NULL,                  # absolute path
      log_file    = NULL,
      step        = RUN_STEP_ALL,          # requested step on the active/last run
      stage       = NA_character_,         # currently running stage key
      stage_times = NULL,                  # named list: <stage> -> start POSIXct
      stage_ends  = NULL,                  # named list: <stage> -> end POSIXct (NULL while running)
      stage_state = NULL,                  # named chr:  <stage> -> pending|running|done|failed
      started_at  = NULL,
      ended_at    = NULL,
      exit_code   = NA_integer_,
      error_msg   = NULL,
      report_path = NULL
    )

    reset_stages <- function() {
      rv$stage_state <- stats::setNames(
        rep("pending", length(RUN_STAGES)),
        names(RUN_STAGES)
      )
      rv$stage_times <- stats::setNames(
        vector("list", length(RUN_STAGES)),
        names(RUN_STAGES)
      )
      rv$stage_ends <- stats::setNames(
        vector("list", length(RUN_STAGES)),
        names(RUN_STAGES)
      )
      rv$stage <- NA_character_
    }

    reset_stages()

    # --- Launch helper: shared by Run-all and per-step buttons -----------
    launch_step <- function(step_key) {
      message(sprintf("[run_controller] click step=%s", step_key))
      st <- ss_state()
      if (!isTRUE(st$valid)) {
        shiny::showNotification(
          "Samplesheet is not ready. Fix validation issues first.",
          type = "warning"
        )
        return(invisible())
      }
      if (identical(rv$state, RUN_STATE_RUNNING)) {
        shiny::showNotification("A run is already in progress.", type = "message")
        return(invisible())
      }

      ws <- if (is.function(workspace)) workspace() else workspace
      pr <- if (is.function(project_root_)) project_root_() else project_root_

      # Per-step (non-"all") runs reuse the current run_dir so they can pick
      # up upstream artifacts. Run-all always creates a fresh run_dir.
      reuse_dir <- !identical(step_key, RUN_STEP_ALL) &&
                   !is.null(rv$run_dir) &&
                   dir.exists(rv$run_dir)

      # preprocess and reference_projection are self-contained — they read
      # the IDATs directly — so they can launch without an existing run_dir.
      # The other per-step runs need upstream artifacts and a run_dir to
      # find them; block the launch and tell the user why.
      self_contained <- step_key %in% c("preprocess", "reference_projection")
      if (!identical(step_key, RUN_STEP_ALL) && !self_contained && !reuse_dir) {
        shiny::showNotification(
          sprintf("Cannot run %s alone — start with a full run (or just preprocess) first.",
                  RUN_STAGES[[step_key]]),
          type = "warning", duration = 8
        )
        return(invisible())
      }
      if (reuse_dir) {
        pr_check <- check_step_prereqs(step_key, rv$run_dir)
        if (!isTRUE(pr_check$ok)) {
          shiny::showNotification(
            sprintf("Cannot run %s — missing upstream artifact(s): %s",
                    RUN_STAGES[[step_key]],
                    paste(basename(unlist(pr_check$missing)), collapse = ", ")),
            type = "warning", duration = 8
          )
          return(invisible())
        }
      }

      if (!reuse_dir) {
        stem <- tools::file_path_sans_ext(basename(st$samplesheet_path))
        run_id <- sprintf("%s_%s", format(Sys.time(), "%Y%m%d-%H%M%S"), stem)
        run_dir <- file.path(ws, "runs", run_id)
        reset_stages()
        rv$run_dir <- run_dir
      } else {
        # Reuse: clear prior state for the target stage only; preserve the
        # done/failed/pending state of every other stage from prior runs.
        rv$stage_state[[step_key]] <- "pending"
        rv$stage_ends[[step_key]]  <- NULL
        rv$stage_times[[step_key]] <- NULL
      }

      # Pre-flight disk-space check on the workspace volume. We don't
      # know exactly how much each run will consume but a rough estimate
      # of ~150 MB/sample + 200 MB buffer catches the obvious "tens of MB
      # free" case before the pipeline crashes mid-stage with ENOSPC.
      free <- free_bytes(rv$run_dir)
      if (!is.na(free)) {
        n_samples <- if (!is.null(st$validated_df)) {
          nrow(st$validated_df)
        } else 1L
        estimated <- 150 * 1024^2 * n_samples + 200 * 1024^2
        if (free < 200 * 1024^2) {
          shiny::showNotification(
            sprintf(
              paste("Less than 200 MB free on the workspace volume",
                    "(%.0f MB). Free space and try again."),
              free / 1024^2
            ),
            type = "error", duration = NULL
          )
          return(invisible())
        }
        if (free < estimated) {
          shiny::showNotification(
            sprintf(
              paste("Workspace volume has %.1f GB free; this run may need",
                    "around %.1f GB. Pipeline will crash mid-run if it",
                    "fills up."),
              free / 1024^3, estimated / 1024^3
            ),
            type = "warning", duration = 12
          )
        }
      }

      # Seed the first stage as "running" right now so its clock starts
      # ticking at click time. The pipeline takes a few seconds to bootstrap
      # before emitting its first "Step N:" log line; without this seed,
      # the user would see no clock movement on the row they clicked, while
      # any prior row with stale stage_times appears to tick instead.
      first_stage <- if (identical(step_key, RUN_STEP_ALL)) {
        "preprocess"
      } else {
        step_key
      }
      now <- Sys.time()
      rv$stage <- first_stage
      rv$stage_state[[first_stage]] <- "running"
      rv$stage_times[[first_stage]] <- now
      rv$stage_ends[[first_stage]]  <- NULL

      rv$step        <- step_key
      rv$started_at  <- now
      rv$ended_at    <- NULL
      rv$exit_code   <- NA_integer_
      rv$error_msg   <- NULL
      rv$report_path <- NULL

      handle <- tryCatch(
        bridge_launch(
          samplesheet     = st$samplesheet_path,
          output_dir      = rv$run_dir,
          data_dir        = file.path(pr, "pipeline", "data"),
          array_type      = st$array_type %||% "auto",
          threads         = 4L,
          step            = step_key,
          pipeline_script = file.path(pr, "pipeline", "methylation_pipeline.R"),
          parameters      = if (is.function(parameters)) parameters() else parameters
        ),
        error = function(e) {
          shiny::showNotification(
            sprintf("Failed to launch pipeline: %s", conditionMessage(e)),
            type = "error", duration = NULL
          )
          NULL
        }
      )
      if (is.null(handle)) return(invisible())

      rv$handle   <- handle
      rv$log_file <- handle$log_file
      rv$state    <- RUN_STATE_RUNNING

      message(sprintf("[run_controller] launched run_id=%s step=%s dir=%s",
                      handle$run_id, step_key, rv$run_dir))
      invisible()
    }

    # --- Run button (full pipeline) --------------------------------------
    shiny::observeEvent(input$run, launch_step(RUN_STEP_ALL),
                        ignoreInit = TRUE)

    # --- Per-step Run buttons --------------------------------------------
    # Observers spelled out one-per-stage. An earlier version registered
    # these via `for (.k in names(RUN_STAGES)) local({...})` which is the
    # textbook Shiny closure-in-loop pattern, but in practice the wrong
    # handler was firing — clicking Step N would launch step N-1 — so we
    # take the safe, explicit route here.
    shiny::observeEvent(input$run_preprocess,
                        launch_step("preprocess"),    ignoreInit = TRUE)
    shiny::observeEvent(input$run_qc,
                        launch_step("qc"),            ignoreInit = TRUE)
    shiny::observeEvent(input$run_dim_reduction,
                        launch_step("dim_reduction"), ignoreInit = TRUE)
    shiny::observeEvent(input$run_reference_projection,
                        launch_step("reference_projection"), ignoreInit = TRUE)
    shiny::observeEvent(input$run_cnv,
                        launch_step("cnv"),           ignoreInit = TRUE)
    shiny::observeEvent(input$run_visualization,
                        launch_step("visualization"), ignoreInit = TRUE)

    # --- Cancel button ---------------------------------------------------
    shiny::observeEvent(input$cancel, {
      if (!identical(rv$state, RUN_STATE_RUNNING)) return()
      if (is.null(rv$handle)) return()

      bridge_kill(rv$handle)
      rv$state     <- RUN_STATE_CANCELLED
      rv$ended_at  <- Sys.time()

      # Mark the in-progress stage as failed (it was interrupted) and
      # freeze its elapsed time.
      if (!is.na(rv$stage) && !is.null(rv$stage_state)) {
        rv$stage_state[[rv$stage]] <- "failed"
        rv$stage_ends[[rv$stage]]  <- rv$ended_at
      }

      # Only quarantine on a full-run cancel. For per-step cancels, the prior
      # stage outputs in this dir are still valid; only the current step's
      # half-finished files are suspect, and the user can re-run that step.
      if (identical(rv$step, RUN_STEP_ALL)) {
        quarantine_partial_outputs(rv$run_dir)
      }
      message(sprintf("[run_controller] cancelled run dir=%s step=%s",
                      rv$run_dir, rv$step))
    }, ignoreInit = TRUE)

    # --- Attach a past run (Theme 6b) ------------------------------------
    # When the past-runs module emits a path, adopt it as our current
    # run_dir. Stage state is synthesized from artifacts on disk: a stage
    # is "done" if its expected output file exists, "pending" otherwise.
    # We don't have historical durations, so stage_times/stage_ends stay
    # NULL and the elapsed-time column shows "—".
    shiny::observeEvent(attach_run(), {
      path <- attach_run()
      if (is.null(path)) return()
      if (identical(rv$state, RUN_STATE_RUNNING)) {
        shiny::showNotification(
          paste("Cannot attach a past run while a pipeline is in progress.",
                "Cancel the current run first."),
          type = "warning", duration = 8
        )
        return()
      }
      if (!dir.exists(path)) {
        shiny::showNotification(
          sprintf("Run directory not found: %s", path),
          type = "warning"
        )
        return()
      }
      rv$run_dir     <- path
      rv$stage_state <- synthesize_stage_state_from_disk(path)
      rv$stage_times <- stats::setNames(vector("list", length(RUN_STAGES)),
                                        names(RUN_STAGES))
      rv$stage_ends  <- stats::setNames(vector("list", length(RUN_STAGES)),
                                        names(RUN_STAGES))
      rv$stage       <- NA_character_
      rv$state       <- RUN_STATE_COMPLETED
      rv$step        <- RUN_STEP_ALL
      rv$exit_code   <- 0L
      rv$error_msg   <- NULL
      rv$report_path <- discover_report(path)
      rv$started_at  <- NULL
      rv$ended_at    <- NULL
      rv$handle      <- NULL
      rv$log_file    <- file.path(path, "logs", "pipeline.log")
      message(sprintf("[run_controller] attached past run dir=%s", path))
    }, ignoreInit = TRUE)

    # --- Poll loop: refreshes every second while a run is active ---------
    shiny::observe({
      if (!identical(rv$state, RUN_STATE_RUNNING)) return()
      if (is.null(rv$handle)) return()

      # Re-trigger this observer every second so long as we're running.
      shiny::invalidateLater(1000, session)

      # 1. Parse log tail for stage transitions.
      lines <- bridge_log_tail(rv$handle, n = 500L)
      update_stage_progress(rv, lines)

      # 2. Detect process exit and finalize state.
      if (!bridge_is_running(rv$handle)) {
        code <- bridge_exit_code(rv$handle)
        rv$exit_code <- if (is.na(code)) NA_integer_ else as.integer(code)
        rv$ended_at  <- Sys.time()

        if (isTRUE(rv$exit_code == 0L)) {
          rv$state <- RUN_STATE_COMPLETED
          # Close out the last running stage first, then mark final state.
          if (!is.na(rv$stage) && !is.null(rv$stage_state)) {
            rv$stage_state[[rv$stage]] <- "done"
            rv$stage_ends[[rv$stage]]  <- rv$ended_at
          }
          if (identical(rv$step, RUN_STEP_ALL)) {
            # Full run: every stage executed.
            rv$stage_state[] <- "done"
          } else if (!is.null(rv$stage_state)) {
            # Per-step run: ensure the requested step is marked done even
            # if we missed its "Step N:" log line. Other stages keep their
            # prior state (done from earlier runs, or pending).
            rv$stage_state[[rv$step]] <- "done"
            if (is.null(rv$stage_ends[[rv$step]])) {
              rv$stage_ends[[rv$step]] <- rv$ended_at
            }
          }
          rv$stage <- NA_character_
          rv$report_path <- discover_report(rv$run_dir)
        } else {
          rv$state <- RUN_STATE_FAILED
          if (!is.na(rv$stage) && !is.null(rv$stage_state)) {
            rv$stage_state[[rv$stage]] <- "failed"
            rv$stage_ends[[rv$stage]]  <- rv$ended_at
          }
          rv$error_msg <- extract_last_error(lines, rv$exit_code)
        }
      }
    })

    # --- UI outputs ------------------------------------------------------
    output$prereq_msg <- shiny::renderUI({
      st <- ss_state()
      if (isTRUE(st$valid)) return(NULL)
      shiny::div(
        class = "alert alert-secondary", role = "alert",
        shiny::tags$strong("Waiting on samplesheet. "),
        "Load a valid samplesheet on the Samplesheet tab to enable Run."
      )
    })

    output$state_badge <- shiny::renderUI({
      render_state_badge(rv$state, rv$stage, rv$started_at, rv$ended_at)
    })

    output$header_badge <- shiny::renderUI({
      # Header version is identical for now; app.R decides whether to show
      # this (run state) or the samplesheet state.
      render_state_badge(rv$state, rv$stage, rv$started_at, rv$ended_at)
    })

    output$stages_panel <- shiny::renderUI({
      # Re-render every second while running so the elapsed-time of the
      # currently-running stage ticks. Finished stages have a frozen end
      # time, so their elapsed value is stable across re-renders.
      if (identical(rv$state, RUN_STATE_RUNNING)) {
        shiny::invalidateLater(1000, session)
      }
      render_stages(rv$stage_state, rv$stage_times, rv$stage_ends,
                    current = rv$stage, run_state = rv$state,
                    run_dir = rv$run_dir, ns = ns)
    })

    output$log_tail <- shiny::renderText({
      if (identical(rv$state, RUN_STATE_RUNNING)) {
        shiny::invalidateLater(1000, session)
      }
      if (is.null(rv$handle)) return("(no run started yet)")
      lines <- bridge_log_tail(rv$handle, n = 20L)
      if (!length(lines)) return("(waiting for first log line…)")
      paste(lines, collapse = "\n")
    })

    output$post_run_actions <- shiny::renderUI({
      render_post_run_actions(ns, rv$state, rv$report_path, rv$log_file,
                              rv$error_msg, rv$run_dir)
    })

    # --- Action handlers -------------------------------------------------
    shiny::observeEvent(input$open_report, {
      if (is.null(rv$report_path) || !file.exists(rv$report_path)) {
        shiny::showNotification("Report file not found.", type = "warning")
        return()
      }
      utils::browseURL(normalizePath(rv$report_path, winslash = "/"))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$reveal_report, {
      if (is.null(rv$report_path) || !file.exists(rv$report_path)) {
        shiny::showNotification("Report file not found.", type = "warning")
        return()
      }
      reveal_in_file_manager(rv$report_path)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$open_log, {
      if (is.null(rv$log_file) || !file.exists(rv$log_file)) {
        shiny::showNotification("Log file not found.", type = "warning")
        return()
      }
      utils::browseURL(normalizePath(rv$log_file, winslash = "/"))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$reveal_run_dir, {
      if (is.null(rv$run_dir) || !dir.exists(rv$run_dir)) {
        shiny::showNotification("Run directory not found.", type = "warning")
        return()
      }
      reveal_in_file_manager(rv$run_dir)
    }, ignoreInit = TRUE)

    # --- Public reactive (for header badge in app.R) --------------------
    shiny::reactive({
      list(
        state     = rv$state,
        stage     = rv$stage,
        run_dir   = rv$run_dir,
        exit_code = rv$exit_code
      )
    })
  })
}

# ---------------------------------------------------------------------------
# Helpers (pure, module-internal)
# ---------------------------------------------------------------------------

# Parse the log tail and advance rv$stage / rv$stage_state / rv$stage_times.
# We only look for "Step N:" markers emitted by methylation_pipeline.R, and
# remap N onto a user-facing stage key via PIPELINE_STEP_TO_STAGE (so the
# pipeline's six internal steps collapse to five UI rows — Step 2 (QC) and
# Step 3 (filtering) both belong to the same "qc" stage).
update_stage_progress <- function(rv, lines) {
  if (!length(lines)) return(invisible())
  matches <- regmatches(lines, regexec(STAGE_LINE_RE, lines))
  step_nums <- vapply(matches, function(m) {
    if (length(m) < 2L) NA_integer_ else as.integer(m[2])
  }, integer(1))
  step_nums <- step_nums[!is.na(step_nums)]
  if (!length(step_nums)) return(invisible())

  latest <- max(step_nums)
  if (latest < 1L || latest > length(PIPELINE_STEP_TO_STAGE)) return(invisible())

  new_stage_key <- PIPELINE_STEP_TO_STAGE[latest]
  # Already on this stage — could be a) we seeded it at click time, or
  # b) we already transitioned to it on a previous poll. Either way, keep
  # the existing start time so the visible clock isn't reset.
  if (identical(rv$stage, new_stage_key)) return(invisible())

  # Guard against backward transitions. launch_step seeds rv$stage at click
  # time; if a poll then catches a stale log line that maps to an earlier UI
  # stage (e.g. a leftover "Step 2:" from a prior per-step QC run before the
  # bridge truncates the log), don't regress — that would reset the clicked
  # stage's clock and re-seed an already-completed stage as running.
  new_idx <- match(new_stage_key, names(RUN_STAGES))
  cur_idx <- if (is.na(rv$stage)) NA_integer_ else match(rv$stage, names(RUN_STAGES))
  if (!is.na(cur_idx) && !is.na(new_idx) && new_idx < cur_idx) {
    return(invisible())
  }

  now <- Sys.time()
  # Close out the previous stage (if any) as done and freeze its elapsed time.
  if (!is.na(rv$stage)) {
    rv$stage_state[[rv$stage]] <- "done"
    rv$stage_ends[[rv$stage]]  <- now
  }
  # Mark every earlier UI stage as done (in case we missed its log line).
  if (!is.na(new_idx) && new_idx > 1L) {
    for (i in seq_len(new_idx - 1L)) {
      k <- names(RUN_STAGES)[i]
      if (identical(rv$stage_state[[k]], "pending")) {
        rv$stage_state[[k]] <- "done"
      }
    }
  }
  rv$stage <- new_stage_key
  rv$stage_state[[new_stage_key]] <- "running"
  rv$stage_times[[new_stage_key]] <- now
  rv$stage_ends[[new_stage_key]]  <- NULL
}

# Render the 6-row stage list. Each row carries a small per-step "Run"
# button; the button is disabled while a run is active or while the step's
# upstream artifacts are missing.
render_stages <- function(stage_state, stage_times, stage_ends, current,
                          run_state, run_dir, ns) {
  # When there's no run yet, show an empty stage list so the per-step
  # buttons are still discoverable (preprocess will be enabled).
  if (is.null(stage_state)) {
    stage_state <- stats::setNames(rep("pending", length(RUN_STAGES)),
                                   names(RUN_STAGES))
    stage_times <- stats::setNames(vector("list", length(RUN_STAGES)),
                                   names(RUN_STAGES))
    stage_ends  <- stage_times
  }
  now <- Sys.time()
  is_running <- identical(run_state, RUN_STATE_RUNNING)

  rows <- lapply(names(RUN_STAGES), function(key) {
    state <- stage_state[[key]] %||% "pending"
    start <- stage_times[[key]]
    end   <- stage_ends[[key]]
    endpoint <- if (!is.null(end)) end else now
    elapsed <- if (is.null(start)) NA else as.numeric(difftime(endpoint, start, units = "secs"))

    badge <- status_pill(state)
    elapsed_txt <- if (is.na(elapsed)) "—" else format_elapsed(elapsed)

    pr_check <- check_step_prereqs(key, run_dir)
    btn_disabled <- is_running ||
      (!identical(key, "preprocess") && !isTRUE(pr_check$ok))
    btn_title <- if (is_running) {
      "A run is already in progress"
    } else if (!isTRUE(pr_check$ok) && !identical(key, "preprocess")) {
      sprintf("Requires upstream output: %s",
              paste(basename(unlist(pr_check$missing)), collapse = ", "))
    } else {
      sprintf("Run only the %s step", RUN_STAGES[[key]])
    }

    btn <- shiny::actionButton(
      inputId = ns(paste0("run_", key)),
      label   = "Run",
      icon    = shiny::icon("play"),
      class   = "btn-sm btn-outline-secondary",
      title   = btn_title
    )
    if (btn_disabled) btn$attribs$disabled <- "disabled"

    shiny::div(
      class = "d-flex align-items-center gap-3 py-1 border-bottom",
      shiny::div(style = "width: 110px;", badge),
      shiny::div(style = "flex: 1;", RUN_STAGES[[key]]),
      shiny::div(class = "text-muted small", style = "width: 90px; text-align: right;",
                 elapsed_txt),
      shiny::div(style = "width: 100px; text-align: right;", btn)
    )
  })
  shiny::tagList(rows)
}

# Synthesize a stage_state vector from artifacts present in a run_dir.
# Used when attaching a past run (Theme 6b) — the original in-memory
# durations are gone, but presence of each stage's output file is enough
# to drive the badges and the per-step prereq checks.
synthesize_stage_state_from_disk <- function(run_dir) {
  ss <- stats::setNames(rep("pending", length(RUN_STAGES)),
                        names(RUN_STAGES))
  if (is.null(run_dir) || !dir.exists(run_dir)) return(ss)
  has <- list(
    preprocess           = file.exists(file.path(run_dir, "processed_data",
                                                 "preprocessed_data.RData")),
    qc                   = file.exists(file.path(run_dir, "qc",
                                                 "qc_results.RData")),
    dim_reduction        = file.exists(file.path(run_dir, "dimensionality_reduction",
                                                 "dim_reduction_results.RData")),
    reference_projection = file.exists(file.path(run_dir, "reference_projection",
                                                 "reference_projection_results.RData")),
    cnv                  = file.exists(file.path(run_dir, "cnv",
                                                 "cnv_results.RData")),
    visualization        = length(list.files(file.path(run_dir, "reports"),
                                              pattern = "\\.html$")) > 0L
  )
  for (k in names(ss)) {
    if (isTRUE(has[[k]])) ss[[k]] <- "done"
  }
  ss
}

# Check whether a per-step run can launch against the given run_dir. Returns
# list(ok = logical, missing = list-of-groups). A "missing" group is a set
# of paths where none exist (i.e. that requirement is unsatisfied).
check_step_prereqs <- function(step_key, run_dir) {
  reqs <- STEP_PREREQS[[step_key]]
  if (is.null(reqs) || length(reqs) == 0L) {
    return(list(ok = TRUE, missing = list()))
  }
  if (is.null(run_dir) || !dir.exists(run_dir)) {
    return(list(ok = FALSE, missing = reqs))
  }
  missing <- list()
  for (group in reqs) {
    paths <- file.path(run_dir, group)
    if (!any(file.exists(paths))) {
      missing[[length(missing) + 1L]] <- group
    }
  }
  list(ok = length(missing) == 0L, missing = missing)
}

render_state_badge <- function(state, stage, started_at, ended_at) {
  switch(
    state,
    idle      = status_pill("idle", "Idle"),
    running   = {
      label <- if (!is.na(stage)) sprintf("Running: %s", RUN_STAGES[[stage]])
               else "Running"
      status_pill("running", label, icon = "spinner")
    },
    completed = status_pill("completed", "Completed", icon = "circle-check"),
    failed    = status_pill("failed",    "Failed",    icon = "triangle-exclamation"),
    cancelled = status_pill("cancelled", "Cancelled", icon = "ban"),
    status_pill("neutral", state)
  )
}

render_post_run_actions <- function(ns, state, report_path, log_file,
                                    error_msg, run_dir) {
  if (identical(state, RUN_STATE_IDLE)) {
    return(shiny::tags$em("Actions appear here after a run finishes."))
  }
  if (identical(state, RUN_STATE_RUNNING)) {
    return(shiny::tags$em("Available once the run finishes."))
  }

  # Completed: report actions + reveal run dir.
  if (identical(state, RUN_STATE_COMPLETED)) {
    has_report <- !is.null(report_path) && file.exists(report_path)
    return(shiny::tagList(
      if (has_report) shiny::tagList(
        shiny::actionButton(ns("open_report"),
                            "Open report in browser",
                            icon = shiny::icon("arrow-up-right-from-square"),
                            class = "btn-success me-2"),
        shiny::actionButton(ns("reveal_report"),
                            reveal_label(),
                            icon = shiny::icon("folder-open"),
                            class = "btn-outline-secondary me-2")
      ) else shiny::div(
        class = "alert alert-warning",
        "Run completed but no HTML report was found under ",
        shiny::tags$code(file.path(basename(run_dir), "report")), "."
      ),
      shiny::actionButton(ns("reveal_run_dir"),
                          "Show run folder",
                          icon = shiny::icon("folder"),
                          class = "btn-outline-secondary me-2"),
      shiny::actionButton(ns("open_log"),
                          "Open log",
                          icon = shiny::icon("file-lines"),
                          class = "btn-outline-secondary")
    ))
  }

  # Failed or cancelled: friendly hint (when available) + raw error line +
  # log + run dir buttons. error_msg is now a list(raw, hint); guard for
  # the legacy bare-string shape just in case.
  raw  <- if (is.list(error_msg)) error_msg$raw  else error_msg
  hint <- if (is.list(error_msg)) error_msg$hint else NULL

  shiny::tagList(
    if (!is.null(hint) && nzchar(hint)) {
      shiny::div(
        class = "alert alert-warning",
        shiny::icon("circle-info"), " ",
        shiny::tags$strong("Likely cause: "), hint
      )
    },
    if (!is.null(raw) && nzchar(raw)) {
      shiny::div(
        class = "alert alert-danger",
        shiny::tags$strong("Error: "), raw
      )
    },
    shiny::actionButton(ns("open_log"),
                        "Open log",
                        icon = shiny::icon("file-lines"),
                        class = "btn-danger me-2"),
    shiny::actionButton(ns("reveal_run_dir"),
                        "Show run folder",
                        icon = shiny::icon("folder"),
                        class = "btn-outline-secondary")
  )
}

reveal_label <- function() {
  # Finder on macOS, Explorer on Windows, file manager on Linux.
  sys <- Sys.info()[["sysname"]]
  if (sys == "Darwin") "Show in Finder"
  else if (sys == "Windows") "Show in Explorer"
  else "Show in file manager"
}

# Cross-platform "reveal" — select the file in the OS file manager when we
# can, otherwise open the containing directory.
reveal_in_file_manager <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  sys  <- Sys.info()[["sysname"]]
  if (sys == "Darwin") {
    system2("open", args = c("-R", shQuote(path)))
  } else if (sys == "Windows") {
    # explorer /select,<path> highlights the file.
    shell(sprintf('explorer /select,"%s"', gsub("/", "\\\\", path, fixed = TRUE)),
          wait = FALSE)
  } else {
    # xdg-open works on the directory, not individual-file selection.
    dir <- if (dir.exists(path)) path else dirname(path)
    system2("xdg-open", args = shQuote(dir))
  }
  invisible()
}

# Find the generated report HTML. The pipeline writes to <run_dir>/reports/
# (plural) as methylation_analysis_report.html; the original design spec
# called out <run_dir>/report/meqtrack_*.html. Check both layouts so either works.
discover_report <- function(run_dir) {
  if (is.null(run_dir)) return(NULL)
  candidates <- c(
    file.path(run_dir, "reports"),
    file.path(run_dir, "report")
  )
  candidates <- candidates[dir.exists(candidates)]
  if (!length(candidates)) return(NULL)
  files <- unlist(lapply(candidates, function(d) {
    list.files(d, pattern = "\\.html$", full.names = TRUE)
  }))
  if (!length(files)) return(NULL)
  # Prefer canonical pipeline filenames; fall back to any .html.
  preferred <- grep("methylation_analysis_report|^meqtrack_", basename(files))
  if (length(preferred)) files <- files[preferred]
  files[order(file.info(files)$mtime, decreasing = TRUE)][1]
}

# Pull the most useful error line from the log tail and pair it with a
# human-readable hint when the line matches a known failure pattern.
# Returns a list(raw, hint) where `hint` may be NULL. The post-run
# actions panel renders both.
extract_last_error <- function(lines, exit_code = NA_integer_) {
  raw <- NULL
  if (length(lines)) {
    hits <- grep("^Error\\b|Error in ", lines, value = TRUE)
    raw <- if (length(hits)) utils::tail(hits, 1L) else utils::tail(lines, 1L)
  }
  hint <- diagnose_failure(lines, raw, exit_code)
  list(raw = raw, hint = hint)
}

# Map common pipeline failure patterns onto a one-line hint. The order
# matters — exit-code signals are checked first, then log content.
diagnose_failure <- function(lines, raw, exit_code = NA_integer_) {
  text <- paste(c(lines, raw), collapse = "\n")

  # 137 = 128 + SIGKILL on Unix; the kernel's OOM killer is the most
  # common source. macOS jetsam behaves similarly. 9 (raw SIGKILL) too.
  if (isTRUE(exit_code %in% c(137L, 9L))) {
    return(paste(
      "Process was killed by the operating system, almost certainly out",
      "of memory. Try running with fewer samples, or on a machine with",
      "more RAM."
    ))
  }

  if (grepl("No space left on device|ENOSPC|Disk quota exceeded",
            text, ignore.case = TRUE)) {
    return(paste(
      "The pipeline ran out of disk space mid-run. Free space on the",
      "workspace volume and retry — partial outputs from this run can",
      "be deleted from the run directory."
    ))
  }

  if (grepl("cannot allocate vector of size", text, ignore.case = TRUE)) {
    return(paste(
      "R ran out of memory allocating an array. This usually means the",
      "sample × probe matrix is too large for available RAM. Try fewer",
      "samples per run, or close other applications."
    ))
  }

  if (grepl("BiocParallel errors|error reading from connection",
            text, ignore.case = TRUE)) {
    return(paste(
      "A parallel worker crashed — usually a fork-related issue with",
      "sesame on macOS. Re-running often clears it; if it persists, set",
      "MEQTRACK_SERIAL=1 in the environment to force serial execution."
    ))
  }

  if (grepl("could not find function", text, ignore.case = TRUE)) {
    return(paste(
      "A package the pipeline expects isn't installed. Re-run setup.R",
      "from the project root to refresh the renv library."
    ))
  }

  if (grepl("there is no package called", text, ignore.case = TRUE)) {
    return(paste(
      "A required R package is missing. Re-run setup.R from the project",
      "root."
    ))
  }

  if (grepl("cannot open .* No such file or directory|cannot find file",
            text, ignore.case = TRUE)) {
    return(paste(
      "The pipeline tried to read a file that isn't on disk. The most",
      "common cause is a samplesheet pointing at an IDAT path that has",
      "moved — re-validate the samplesheet before re-running."
    ))
  }

  NULL
}

# Move partial stage subdirs into <run_dir>/_cancelled/ so a cancelled run
# cannot be mistaken for a successful one. Keeps logs/ and the manifest at
# the top level so the user can still inspect what happened.
quarantine_partial_outputs <- function(run_dir) {
  if (is.null(run_dir) || !dir.exists(run_dir)) return(invisible())
  quarantine <- file.path(run_dir, "_cancelled")
  dir.create(quarantine, recursive = TRUE, showWarnings = FALSE)
  keep <- c("logs", "_cancelled", "run_manifest.json")
  entries <- list.files(run_dir, full.names = FALSE, include.dirs = TRUE,
                        no.. = TRUE)
  for (e in setdiff(entries, keep)) {
    src <- file.path(run_dir, e)
    dst <- file.path(quarantine, e)
    tryCatch(
      file.rename(src, dst),
      warning = function(w) {
        # file.rename can fail across filesystems; fall back to copy+delete.
        file.copy(src, quarantine, recursive = TRUE)
        unlink(src, recursive = TRUE, force = TRUE)
      }
    )
  }
  invisible()
}

format_elapsed <- function(secs) {
  if (is.na(secs)) return("—")
  if (secs < 60) return(sprintf("%ds", round(secs)))
  if (secs < 3600) return(sprintf("%dm %02ds", secs %/% 60, round(secs %% 60)))
  sprintf("%dh %02dm", secs %/% 3600, (secs %% 3600) %/% 60)
}
