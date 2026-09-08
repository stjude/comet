# scripts/test_bridge.R — Wave 1 smoke test for pipeline_bridge.R.
#
# Exercises the bridge API (launch, poll, log-tail, exit code) against the
# bundled example samplesheet. No Shiny; just proves the bridge works from a
# plain R script.
#
# Usage (from project root):
#   Rscript scripts/test_bridge.R

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

source(file.path("app", "R", "pipeline_bridge.R"))

project_root <- normalizePath(".", winslash = "/")
samplesheet <- file.path(
  project_root, "pipeline", "data", "example", "samplesheet_epic.csv"
)
data_dir <- file.path(project_root, "pipeline", "data")

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
output_dir <- file.path(project_root, "runs", paste0("bridge_test_", timestamp))

message("Launching pipeline via bridge...")
handle <- bridge_launch(
  samplesheet = samplesheet,
  output_dir  = output_dir,
  data_dir    = data_dir,
  array_type  = "EPIC",
  threads     = 4L,
  step        = "all"
)

message(sprintf("Run id: %s", handle$run_id))
message(sprintf("Log file: %s", handle$log_file))

# Poll until done, printing a tail every 10 seconds.
while (bridge_is_running(handle)) {
  Sys.sleep(10)
  tail_lines <- bridge_log_tail(handle, n = 5L)
  if (length(tail_lines) > 0) {
    cat(sprintf("--- tail (%s) ---\n", format(Sys.time(), "%H:%M:%S")))
    cat(paste(tail_lines, collapse = "\n"), "\n")
  }
}

code <- bridge_exit_code(handle)
message(sprintf("Pipeline exited with code: %s", code))

if (!isTRUE(code == 0L)) {
  message("Bridge test FAILED. Recent log lines:")
  cat(paste(bridge_log_tail(handle, n = 40L), collapse = "\n"), "\n")
  quit(status = 1L)
}

message("Bridge test OK. Output at: ", output_dir)
