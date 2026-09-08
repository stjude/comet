# app/R/results_loader.R
# ---------------------------------------------------------------------------
# Load all artifacts from a completed pipeline run and expose them as a
# single reactive bundle consumed by the Wave 4 result modules.
#
# Contract (return value of results_loader_server()): a reactive that yields
# NULL when there is no completed run, or a list with:
#   run_dir         absolute path to the run directory
#   qc_report       data.frame — per-sample QC table (sample_qc_report.csv)
#   conversion_qc   data.frame — GCT bisulfite-conversion QC (conversion_qc.csv) or NULL
#   qc_fail_ids     character  — Sample_IDs that failed QC
#   sample_info     data.frame — sample metadata (processed_data/sample_info.txt)
#   tsne            list(coords, sample_info, duplicates)   or NULL
#   umap            list(coords, sample_info, umap_object)  or NULL
#   hclust          list(hclust, distance, method)          or NULL
#   cnv             list(segments, sample_results, ...)     or NULL
#   reference_projection  list(dataset, label, projected, ref_meta, ...) or NULL
#   metadata_cols   character  — samplesheet columns eligible as color
#                                variables (Sample_Group, Diagnosis, etc.)
#   run_url_base    character  — web path (served via addResourcePath) from
#                                which the run's files are reachable, e.g.
#                                "runs/<id>" — modules concat with "figures/..."
# ---------------------------------------------------------------------------

# File layout produced by the pipeline (relative to run_dir):
RESULTS_PATHS <- list(
  qc_report       = file.path("qc", "sample_qc_report.csv"),
  conversion_qc   = file.path("qc", "conversion_qc.csv"),
  snp_concordance = file.path("qc", "snp_concordance.csv"),
  deconv          = file.path("deconv", "cell_fractions.csv"),
  qc_rdata        = file.path("qc", "qc_results.RData"),
  sample_info     = file.path("processed_data", "sample_info.txt"),
  tsne_rdata      = file.path("dimensionality_reduction", "tsne_results.RData"),
  umap_rdata      = file.path("dimensionality_reduction", "umap_results.RData"),
  hclust_rdata    = file.path("dimensionality_reduction", "hclust_results.RData"),
  cnv_rdata       = file.path("cnv", "cnv_results.RData"),
  refproj_rdata   = file.path("reference_projection", "reference_projection_results.RData")
)

# Columns we never offer as a "coloring variable" on t-SNE/UMAP — they are
# either identifiers or pipeline-internal.
METADATA_EXCLUDED_COLS <- c(
  "Sample_ID", "Sentrix_ID", "Sentrix_Position",
  "Basename", "Array", "Slide", "Pool_ID"
)

#' Load all artifacts present under a run directory and return them as a
#' bundle. Pieces that haven't been produced yet are returned as NULL so
#' downstream modules can render placeholders for the missing slots —
#' this is what enables Theme 6f's incremental result tabs.
#'
#' Returns NULL only when the run directory itself is missing or empty;
#' a bundle with every piece NULL is still preferable to a NULL bundle
#' because it lets modules surface "still running, nothing yet" placeholders
#' instead of "no completed run."
load_results_bundle <- function(run_dir, run_url_base = NULL) {
  if (is.null(run_dir) || !dir.exists(run_dir)) return(NULL)

  qc_report_path <- file.path(run_dir, RESULTS_PATHS$qc_report)
  qc_report <- if (file.exists(qc_report_path)) {
    tryCatch(
      utils::read.csv(qc_report_path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  } else NULL

  qc_fail_ids <- if (!is.null(qc_report) && "Pass_QC" %in% colnames(qc_report)) {
    as.character(qc_report$Sample_ID[!as.logical(qc_report$Pass_QC)])
  } else character(0)

  conversion_qc_path <- file.path(run_dir, RESULTS_PATHS$conversion_qc)
  conversion_qc <- if (file.exists(conversion_qc_path)) {
    tryCatch(
      utils::read.csv(conversion_qc_path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  } else NULL

  # Pairwise SNP genotype-concordance matrix (sample-identity heatmap). First
  # column holds the Sample_IDs (row names); keep names verbatim.
  snp_concordance_path <- file.path(run_dir, RESULTS_PATHS$snp_concordance)
  snp_concordance <- if (file.exists(snp_concordance_path)) {
    tryCatch(
      as.matrix(utils::read.csv(snp_concordance_path, row.names = 1,
                                check.names = FALSE)),
      error = function(e) NULL
    )
  } else NULL

  # Cell-type deconvolution (deconvMe) — tidy long: method, sample, celltype, value.
  deconv_path <- file.path(run_dir, RESULTS_PATHS$deconv)
  deconv <- if (file.exists(deconv_path)) {
    tryCatch(
      utils::read.csv(deconv_path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
  } else NULL

  sample_info_path <- file.path(run_dir, RESULTS_PATHS$sample_info)
  sample_info <- if (file.exists(sample_info_path)) {
    tryCatch(
      utils::read.table(sample_info_path,
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                        check.names = FALSE),
      error = function(e) NULL
    )
  } else NULL

  tsne   <- .load_rdata(file.path(run_dir, RESULTS_PATHS$tsne_rdata),   "tsne_results")
  umap   <- .load_rdata(file.path(run_dir, RESULTS_PATHS$umap_rdata),   "umap_results")
  hclust <- .load_rdata(file.path(run_dir, RESULTS_PATHS$hclust_rdata), "hclust_results")
  cnv    <- .load_rdata(file.path(run_dir, RESULTS_PATHS$cnv_rdata),    "cnv_results")
  reference_projection <- .load_rdata(
    file.path(run_dir, RESULTS_PATHS$refproj_rdata), "rp_result")

  # If literally nothing is on disk yet (very early in a run), don't return
  # a hollow bundle — let consumers keep showing the "no run yet" empty state.
  if (is.null(qc_report) && is.null(conversion_qc) && is.null(sample_info) &&
      is.null(tsne) && is.null(umap) && is.null(hclust) && is.null(cnv) &&
      is.null(reference_projection)) {
    return(NULL)
  }

  metadata_cols <- detect_metadata_cols(sample_info)

  list(
    run_dir       = run_dir,
    qc_report     = qc_report,
    conversion_qc = conversion_qc,
    snp_concordance = snp_concordance,
    deconv        = deconv,
    qc_fail_ids   = qc_fail_ids,
    sample_info   = sample_info,
    tsne          = tsne,
    umap          = umap,
    hclust        = hclust,
    cnv           = cnv,
    reference_projection = reference_projection,
    metadata_cols = metadata_cols,
    run_url_base  = run_url_base
  )
}

#' Shiny module that watches the run controller's state reactive and emits
#' a results bundle whenever artifacts are available. Polls every 2s while
#' a run is in progress so the QC/Dim-reduction/CNV tabs populate as each
#' pipeline stage finishes (Theme 6f), instead of only at end-of-run.
#'
#' @param id           module id
#' @param run_state    reactive — output of run_controller_server()
#' @param runs_url_prefix  the addResourcePath prefix for runs/. App.R
#'                     registers "runs" pointing at <workspace>/runs/.
results_loader_server <- function(id, run_state, runs_url_prefix = "runs") {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::reactive({
      rs <- run_state()
      if (is.null(rs) || is.null(rs$run_dir)) return(NULL)
      # Show partial results during RUNNING and final results after
      # COMPLETED. For other states (idle / failed / cancelled) only emit
      # if there's still a run_dir we can read from — failed-mid-run still
      # has whatever earlier stages produced.
      if (identical(rs$state, RUN_STATE_IDLE)) return(NULL)
      if (identical(rs$state, RUN_STATE_RUNNING)) {
        # Re-trigger this reactive every 2 seconds so newly-produced
        # artifacts are picked up while the pipeline is still running.
        shiny::invalidateLater(2000, session)
      }
      url_base <- paste0(runs_url_prefix, "/", basename(rs$run_dir))
      load_results_bundle(rs$run_dir, run_url_base = url_base)
    })
  })
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Pick columns from sample_info that could reasonably color a scatter plot:
# anything except known identifier columns, and with at least one non-NA
# value, and with at most 50 distinct values (avoid coloring by cg probe
# names if they accidentally leak in).
detect_metadata_cols <- function(sample_info) {
  if (is.null(sample_info) || !nrow(sample_info)) return(character(0))
  cols <- setdiff(colnames(sample_info), METADATA_EXCLUDED_COLS)
  keep <- vapply(cols, function(c) {
    v <- sample_info[[c]]
    if (all(is.na(v)) || all(!nzchar(as.character(v)))) return(FALSE)
    length(unique(v)) <= 50L
  }, logical(1))
  cols[keep]
}

.load_rdata <- function(path, obj_name) {
  if (!file.exists(path)) return(NULL)
  tryCatch({
    e <- new.env()
    load(path, envir = e)
    if (!exists(obj_name, envir = e, inherits = FALSE)) return(NULL)
    get(obj_name, envir = e)
  }, error = function(e) NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
