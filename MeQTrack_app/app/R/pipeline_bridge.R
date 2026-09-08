# app/R/pipeline_bridge.R
# ---------------------------------------------------------------------------
# Thin wrapper around the existing R pipeline, launched as a background
# process via `callr::r_bg`.
#
# This module is the single boundary between the Shiny UI (to be built in
# Wave 2/3) and the existing pipeline/methylation_pipeline.R driver.
#
# Design goals:
#   * Fire-and-forget launch: returns a handle immediately.
#   * Non-blocking: UI stays responsive while the pipeline runs.
#   * Isolation: pipeline crash cannot take down the Shiny session.
#   * Live log tailing: the child's stdout + stderr are redirected to a
#     single log file that the UI can read on a poll timer.
#   * Stage tracking: we don't reinvent stages; the pipeline already names
#     them in its log output, and the UI can parse the tail for transitions.
#
# Public API (exported by sourcing this file):
#   bridge_launch(samplesheet, output_dir, data_dir, array_type, threads, step)
#     -> list(proc = <callr process>, log_file = <path>, run_id = <id>,
#             output_dir = <path>)
#   bridge_is_running(handle) -> logical
#   bridge_exit_code(handle) -> integer | NA (NA while running)
#   bridge_log_tail(handle, n = 100) -> character vector
#   bridge_kill(handle) -> invisibly the handle (marked cancelled on disk)
# ---------------------------------------------------------------------------

# Required (declared here for clarity; callers are responsible for install).
# - callr      (CRAN)
# - jsonlite   (CRAN)

#' Launch a pipeline run in a background R process.
#'
#' @param samplesheet  Absolute path to the samplesheet CSV.
#' @param output_dir   Absolute path to the run output directory. Created if missing.
#' @param data_dir     Absolute path to pipeline/data (for keep-probes, yamapData).
#' @param array_type   One of "450K", "EPIC", "EPICv2", or "auto".
#' @param threads      Integer number of threads to pass to --threads.
#' @param step         One of "preprocess", "qc", "filtering", "dim_reduction",
#'                     "cnv", "visualization", or "all".
#' @param pipeline_script  Path to methylation_pipeline.R. Defaults to
#'                     <project_root>/pipeline/methylation_pipeline.R.
#' @return A handle list: proc, log_file, run_id, output_dir.
bridge_launch <- function(samplesheet,
                          output_dir,
                          data_dir,
                          array_type = "auto",
                          threads = 4L,
                          step = "all",
                          pipeline_script = NULL,
                          parameters = list()) {

  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package 'callr' is required. Run setup.R.")
  }

  # Resolve the pipeline script relative to the current working directory
  # if not given. Callers (test_bridge.R, Shiny) are expected to cd into the
  # project root, which matches the documented project layout.
  if (is.null(pipeline_script)) {
    pipeline_script <- file.path(
      normalizePath(".", winslash = "/"),
      "pipeline", "methylation_pipeline.R"
    )
  }

  stopifnot(
    file.exists(samplesheet),
    file.exists(pipeline_script),
    dir.exists(data_dir)
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  logs_dir <- file.path(output_dir, "logs")
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(logs_dir, "pipeline.log")

  # Truncate the log file synchronously before returning. system2 in the
  # worker eventually truncates via shell redirect (>), but there is a gap
  # between bridge_launch returning and the worker actually getting there;
  # during that gap a UI poll reading the file would see stale "Step N:"
  # lines from a prior per-step run targeting this same run_dir, which
  # would mislead the stage-progress parser into transitioning backwards.
  cat("", file = log_file)

  # If the user has set any tunable parameters, write a small per-run
  # config file the pipeline will deep-merge over its defaults. Empty or
  # all-default parameters → no config file, no --config arg.
  config_path <- .write_run_config(output_dir, parameters)

  # Write a tiny run manifest so the UI can reason about the run later.
  # The parameters block is what the past-runs library reads back to
  # populate the Settings inputs on attach.
  run_id <- basename(output_dir)
  manifest <- list(
    run_id = run_id,
    started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    samplesheet = samplesheet,
    output_dir = output_dir,
    data_dir = data_dir,
    array_type = array_type,
    threads = threads,
    step = step,
    parameters = parameters,
    pipeline_script = pipeline_script,
    status = "running"
  )
  .write_manifest(output_dir, manifest)

  # Build the argv. We invoke the pipeline via Rscript so the child is a
  # fresh R process (no Shiny reactive state, no package load conflicts).
  argv <- c(
    pipeline_script,
    "--input", samplesheet,
    "--output", output_dir,
    "--data_dir", data_dir,
    "--array_type", array_type,
    "--threads", as.character(threads),
    "--step", step
  )
  if (!is.null(config_path)) {
    argv <- c(argv, "--config", config_path)
  }

  # callr::r_bg runs an R function in a background process. Inside that
  # child we shell out to Rscript to drive the existing CLI pipeline and
  # redirect its stdout + stderr to a single tailable log file.
  handle_proc <- callr::r_bg(
    func = function(script, args, log_file) {
      result <- tryCatch(
        {
          system2(
            file.path(R.home("bin"), "Rscript"),
            args = c(shQuote(script), args),
            stdout = log_file,   # capture child's stdout to log file
            stderr = log_file    # merge stderr into the same file
          )
        },
        error = function(e) {
          cat("[pipeline_bridge] Child exception: ", conditionMessage(e),
              "\n", file = log_file, append = TRUE)
          1L
        }
      )
      as.integer(result)
    },
    args = list(script = pipeline_script, args = argv[-1], log_file = log_file),
    stdout = file.path(logs_dir, "bridge.out"),
    stderr = file.path(logs_dir, "bridge.err"),
    supervise = TRUE
  )

  handle <- list(
    proc = handle_proc,
    log_file = log_file,
    run_id = run_id,
    output_dir = output_dir,
    manifest_path = file.path(output_dir, "run_manifest.json")
  )
  class(handle) <- c("pipeline_bridge_handle", "list")
  handle
}

#' @rdname bridge_launch
bridge_is_running <- function(handle) {
  if (!inherits(handle, "pipeline_bridge_handle")) {
    stop("Not a pipeline_bridge_handle")
  }
  handle$proc$is_alive()
}

#' @rdname bridge_launch
bridge_exit_code <- function(handle) {
  if (bridge_is_running(handle)) return(NA_integer_)
  res <- tryCatch(handle$proc$get_exit_status(), error = function(e) NA_integer_)
  .update_manifest_on_exit(handle, res)
  as.integer(res)
}

#' @rdname bridge_launch
bridge_log_tail <- function(handle, n = 100L) {
  if (!file.exists(handle$log_file)) return(character(0))
  lines <- tryCatch(readLines(handle$log_file, warn = FALSE),
                    error = function(e) character(0))
  tail(lines, n)
}

#' @rdname bridge_launch
bridge_kill <- function(handle) {
  if (bridge_is_running(handle)) {
    handle$proc$kill()
  }
  .update_manifest_on_exit(handle, NA_integer_, override_status = "cancelled")
  invisible(handle)
}

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Translate the UI's flat-keyed parameters list into the nested config
# structure the pipeline expects, then write it as run_config.R into the
# run_dir. Returns the file path (or NULL if there's nothing to write).
.write_run_config <- function(output_dir, parameters) {
  if (is.null(parameters) || !length(parameters)) return(NULL)
  cfg <- list()
  if (!is.null(parameters$qc.detection_p_threshold)) {
    cfg$qc$detection_p_threshold <- parameters$qc.detection_p_threshold
  }
  if (!is.null(parameters$qc.failed_probe_percent_threshold)) {
    cfg$qc$failed_probe_percent_threshold <-
      parameters$qc.failed_probe_percent_threshold
  }
  if (!is.null(parameters$qc.max_gct_score)) {
    cfg$qc$max_gct_score <- parameters$qc.max_gct_score
  }
  if (!is.null(parameters$dim.variable_probes)) {
    cfg$dim_reduction$variable_probes <- parameters$dim.variable_probes
  }
  if (!is.null(parameters$dim.tsne_perplexity)) {
    cfg$dim_reduction$tsne$perplexity <- parameters$dim.tsne_perplexity
  }
  if (!is.null(parameters$dim.umap_n_neighbors)) {
    cfg$dim_reduction$umap$n_neighbors <- parameters$dim.umap_n_neighbors
  }
  if (!is.null(parameters$cnv.gain_threshold)) {
    cfg$cnv$gain_threshold <- parameters$cnv.gain_threshold
  }
  if (!is.null(parameters$cnv.loss_threshold)) {
    cfg$cnv$loss_threshold <- parameters$cnv.loss_threshold
  }
  # Back-compat: a legacy single `cnv.threshold` still flows through
  # untouched. The pipeline applies it symmetrically when neither
  # gain_threshold nor loss_threshold is set.
  if (!is.null(parameters$cnv.threshold)) {
    cfg$cnv$threshold <- parameters$cnv.threshold
  }
  if (!is.null(parameters$refproj.dataset)) {
    cfg$reference_projection$dataset <- parameters$refproj.dataset
  }
  if (!is.null(parameters$refproj.knn_k)) {
    cfg$reference_projection$knn_k <- parameters$refproj.knn_k
  }
  if (!is.null(parameters$refproj.perplexity)) {
    cfg$reference_projection$perplexity <- parameters$refproj.perplexity
  }
  if (!length(cfg)) return(NULL)
  path <- file.path(output_dir, "run_config.R")
  # Serialize via dput so we don't hand-roll R syntax. Wrapping in a named
  # assignment is what load_config() expects (it sources the file and
  # reads the `config` variable).
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  cat("# Generated by MeQTrack settings_module — overrides only.\n", file = con)
  cat("config <- ", file = con)
  dput(cfg, file = con)
  path
}

.write_manifest <- function(output_dir, manifest) {
  jsonlite::write_json(
    manifest,
    file.path(output_dir, "run_manifest.json"),
    auto_unbox = TRUE, pretty = TRUE
  )
}

.update_manifest_on_exit <- function(handle, exit_code, override_status = NULL) {
  if (!file.exists(handle$manifest_path)) return(invisible(NULL))
  m <- tryCatch(
    jsonlite::read_json(handle$manifest_path, simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(m)) return(invisible(NULL))
  # Only write once; don't repeatedly flip status on poll.
  if (!is.null(m$ended_at)) return(invisible(NULL))
  m$ended_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  m$exit_code <- exit_code
  m$status <- if (!is.null(override_status)) {
    override_status
  } else if (isTRUE(exit_code == 0L)) {
    "completed"
  } else {
    "failed"
  }
  .write_manifest(handle$output_dir, m)
}
