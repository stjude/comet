# scripts/run_example.R — Wave 1 smoke test.
#
# Runs the pipeline end-to-end on the bundled example samplesheet:
#   pipeline/data/example/samplesheet_epic.csv
#
# Usage (from the project root):
#   Rscript scripts/run_example.R
#
# Expected outcome:
#   runs/example_<timestamp>/... is created with subfolders:
#     preprocess/ qc/ filtering/ dim_reduction/ cnv/ report/ logs/
#   and runs/example_<timestamp>/report/ contains an .html report.

# Activate renv so we get the project's library.
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
} else {
  stop("renv is not activated. Run `Rscript setup.R` first.")
}

# Paths (all anchored at project root).
project_root <- normalizePath(".", winslash = "/")
pipeline_script <- file.path(project_root, "pipeline", "methylation_pipeline.R")
samplesheet <- file.path(
  project_root, "pipeline", "data", "example", "samplesheet_epic.csv"
)
data_dir <- file.path(project_root, "pipeline", "data")

if (!file.exists(pipeline_script)) {
  stop("Pipeline script not found at ", pipeline_script)
}
if (!file.exists(samplesheet)) {
  stop("Example samplesheet not found at ", samplesheet)
}

# Run output goes into runs/example_<timestamp>/
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
output_dir <- file.path(project_root, "runs", paste0("example_", timestamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf("Running example pipeline into %s ...", output_dir))

# We invoke the CLI driver rather than sourcing it, because methylation_pipeline.R
# does setwd(script_dir) at startup and uses optparse to parse its own args.
result <- system2(
  "Rscript",
  args = c(
    shQuote(pipeline_script),
    "--input", shQuote(samplesheet),
    "--output", shQuote(output_dir),
    "--data_dir", shQuote(data_dir),
    "--array_type", "EPIC",
    "--threads", "4",
    "--step", "all"
  ),
  stdout = "",   # stream to our own stdout
  stderr = ""
)

if (result != 0) {
  stop(sprintf(
    "Pipeline exited with status %d. See stderr above for details.",
    result
  ))
}

# Verify the expected output tree exists.
expected_subdirs <- c(
  "preprocess", "qc", "filtering", "dim_reduction", "cnv", "report", "logs"
)
missing <- expected_subdirs[!vapply(
  file.path(output_dir, expected_subdirs),
  dir.exists,
  logical(1)
)]

if (length(missing) > 0) {
  message(sprintf(
    "Warning: expected subfolders missing from %s: %s",
    output_dir, paste(missing, collapse = ", ")
  ))
}

# Locate the report.
report_files <- list.files(
  file.path(output_dir, "report"),
  pattern = "\\.html$",
  full.names = TRUE
)

if (length(report_files) == 0) {
  message("No HTML report found in report/. This may be OK if rmarkdown/pandoc are missing; ",
          "check logs/ for details.")
} else {
  message(sprintf("Report: %s", report_files[1]))
}

message("\nDone. Open the report above in a browser to inspect the result.")
