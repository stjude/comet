# ---------------------------------------------------------------------------
# Cross-platform cohort merging
#
# Public API:
#   merge_platform_betas(sources, policy = "intersection")
#   merge_platform_sample_info(sources)
#   merge_platform_segments(sources, output_file)
#   merge_platforms(sources, output_dir, policy)
#   correct_platform_batch(beta, platform, covariates = NULL)
#   platform_batch_diagnostics(beta, platform, control = NULL)
#
# A "source" is one completed single-platform run directory plus the platform
# label it was run as. Preprocessing must be per-platform: minfi's
# read.metharray.exp(force = TRUE) silently truncates to the smallest common
# probe set when IDAT sizes differ, so a mixed cohort read in one shot yields a
# corrupted probe space. We therefore never merge upstream of beta values.
#
# Note on probe IDs: the preprocessing step already strips EPICv2 design
# suffixes (cg00000029_TC21 -> cg00000029) and collapses replicate probes to
# per-CpG means, so beta_values.txt carries plain CpG IDs on every platform.
# .normalize_probe_ids() below is retained as a defensive no-op for the case
# where a caller passes a matrix that has not been through preprocessing.
# ---------------------------------------------------------------------------

#' Strip EPICv2 design suffixes and collapse replicate CpGs to per-CpG means.
#'
#' Mirrors reference_projection.R's .normalize_probe_ids(). A no-op when
#' rownames carry no suffix (450K / EPIC v1, and post-preprocessing EPICv2).
#'
#' @param beta Numeric matrix, probes x samples.
#' @return Matrix with plain-CpG rownames, replicates averaged.
.normalize_probe_ids <- function(beta) {
  beta <- as.matrix(beta)
  pid  <- rownames(beta)
  base <- sub("_.*$", "", pid)
  if (identical(base, pid)) return(beta)

  if (anyDuplicated(base) == 0L) {
    rownames(beta) <- base
    return(beta)
  }
  sums   <- rowsum(beta, base, na.rm = TRUE)
  counts <- rowsum(1 * !is.na(beta), base)
  counts[counts == 0] <- NA_real_          # all-NA group -> NA, not 0/0
  sums / counts
}

#' Read a beta matrix written by write_beta_values().
#'
#' @param path Path to beta_values.txt (or .txt.gz), first column ProbeID.
#' @return Numeric matrix, probes x samples.
.read_beta_file <- function(path) {
  if (!file.exists(path)) stop("Beta file not found: ", path)
  df <- utils::read.delim(path, row.names = 1, check.names = FALSE,
                          stringsAsFactors = FALSE)
  as.matrix(df)
}

#' Normalise the `sources` argument into a validated data frame.
#'
#' Accepts either a named character vector (names = platform labels, values =
#' run directories) or a data frame with `platform` and `dir` columns.
.as_sources <- function(sources) {
  if (is.data.frame(sources)) {
    stopifnot(all(c("platform", "dir") %in% names(sources)))
    out <- sources[, c("platform", "dir")]
  } else {
    if (is.null(names(sources)) || any(!nzchar(names(sources)))) {
      stop("`sources` must be a named vector: names = platform, values = run dir.")
    }
    out <- data.frame(platform = names(sources),
                      dir      = unname(as.character(sources)),
                      stringsAsFactors = FALSE)
  }
  if (nrow(out) < 2) {
    stop("Merging needs at least 2 platform sources; got ", nrow(out), ".")
  }
  if (anyDuplicated(out$platform)) {
    stop("Duplicate platform labels in `sources`: ",
         paste(out$platform[duplicated(out$platform)], collapse = ", "))
  }
  missing <- !dir.exists(out$dir)
  if (any(missing)) {
    stop("Run directory not found for platform(s): ",
         paste(sprintf("%s (%s)", out$platform[missing], out$dir[missing]),
               collapse = "; "))
  }
  out
}

#' Merge per-platform beta matrices into one cohort matrix.
#'
#' Policy "intersection" keeps only CpGs present on every platform, so no value
#' is ever imputed and every sample carries real measurements. This is the safe
#' default for joint dimensionality reduction: dim_reduction.R drops any probe
#' with an NA (see run_tsne/run_umap), and it does so *after*
#' select_variable_probes() has already ranked probes — so a union-style merge
#' would let imputed or missing values steer variable-probe selection and then
#' silently collapse back toward the intersection anyway.
#'
#' @param sources Named vector or data frame of platform -> run directory.
#' @param policy  Currently only "intersection".
#' @param which   "raw" (beta_values.txt) or "filtered"
#'                (filtered_beta_values.txt).
#' @return List: beta (matrix), probes_per_platform, n_common, platform_of_sample.
merge_platform_betas <- function(sources, policy = "intersection",
                                 which = "raw") {
  src <- .as_sources(sources)
  policy <- match.arg(policy, c("intersection"))
  fname <- switch(match.arg(which, c("raw", "filtered")),
                  raw      = "beta_values.txt",
                  filtered = "filtered_beta_values.txt")

  mats <- list()
  n_probes <- integer(nrow(src))
  for (i in seq_len(nrow(src))) {
    p    <- src$platform[i]
    path <- file.path(src$dir[i], "processed_data", fname)
    if (!file.exists(path) && file.exists(paste0(path, ".gz"))) {
      path <- paste0(path, ".gz")
    }
    message(sprintf("merge_platform_betas: reading %s (%s)...", p, basename(path)))
    m <- .normalize_probe_ids(.read_beta_file(path))
    if (anyDuplicated(colnames(m))) {
      stop("Duplicate sample columns within platform ", p, ".")
    }
    mats[[p]]   <- m
    n_probes[i] <- nrow(m)
    message(sprintf("  %s: %d probes x %d samples", p, nrow(m), ncol(m)))
  }
  names(n_probes) <- src$platform

  common <- Reduce(intersect, lapply(mats, rownames))
  if (length(common) == 0) {
    stop("No probes shared across platforms — cannot merge.")
  }
  common <- sort(common)

  # Guard against the same sample appearing under two platforms.
  all_cols <- unlist(lapply(mats, colnames), use.names = FALSE)
  if (anyDuplicated(all_cols)) {
    dup <- unique(all_cols[duplicated(all_cols)])
    stop("Sample(s) present in more than one platform source: ",
         paste(dup, collapse = ", "))
  }

  beta <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
  rownames(beta) <- common

  platform_of_sample <- rep(names(mats), vapply(mats, ncol, integer(1)))
  names(platform_of_sample) <- all_cols

  # Judge the intersection against the SMALLEST platform, not the largest: the
  # smallest manifest is the hard ceiling on what any intersection can retain,
  # so measuring against the largest would flag healthy merges as broken. A
  # 450K+EPICv2 merge, for instance, keeps ~81% of 450K but only ~40% of
  # EPICv2 -- entirely expected, since 450K carries ~485K probes to EPICv2's
  # ~924K. Reference values for the manifests (EPICv2 suffix-stripped):
  #   450K+EPIC 452K (93% of 450K), 450K+EPICv2 394K (81%),
  #   EPIC+EPICv2 721K (83%), all three 370K (76%).
  pct <- 100 * length(common) / min(n_probes)
  smallest <- names(n_probes)[which.min(n_probes)]
  message(sprintf(
    "merge_platform_betas: %d common probes (%.1f%% of the smallest platform, %s at %d); %d samples total.",
    length(common), pct, smallest, min(n_probes), ncol(beta)))
  if (pct < 60) {
    warning(sprintf(
      "merge_platform_betas: intersection retains only %.1f%% of the smallest platform (%s) — check that platforms were preprocessed consistently.",
      pct, smallest))
  }

  # NA audit: intersection removes platform-specific probes but not probes that
  # simply failed detection, so report what dim_reduction will later drop.
  n_na_rows <- sum(apply(beta, 1, anyNA))
  if (n_na_rows > 0) {
    message(sprintf(
      "  note: %d/%d common probes carry >=1 NA and will be dropped by dim reduction.",
      n_na_rows, length(common)))
  }

  list(beta               = beta,
       probes_per_platform = n_probes,
       n_common           = length(common),
       n_na_rows          = n_na_rows,
       platform_of_sample = platform_of_sample)
}

#' Concatenate per-platform sample_info, adding a Platform column.
#'
#' Optionally joins extra columns (e.g. Group) from each platform's
#' samplesheet, which sample_info.txt does not carry.
#'
#' @param sources      Named vector or data frame of platform -> run directory.
#' @param samplesheets Optional named vector platform -> samplesheet CSV.
#' @param join_cols    Columns to lift from the samplesheet when available.
#' @return Data frame with Sample_ID, Platform, and any joined columns.
merge_platform_sample_info <- function(sources, samplesheets = NULL,
                                       join_cols = c("Group")) {
  src <- .as_sources(sources)
  out <- list()

  for (i in seq_len(nrow(src))) {
    p    <- src$platform[i]
    path <- file.path(src$dir[i], "processed_data", "sample_info.txt")
    if (!file.exists(path)) stop("sample_info.txt not found for platform ", p)
    si <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    si$Platform <- p

    if (!is.null(samplesheets) && p %in% names(samplesheets) &&
        file.exists(samplesheets[[p]])) {
      ss <- utils::read.csv(samplesheets[[p]], stringsAsFactors = FALSE,
                            check.names = FALSE)
      key <- intersect(c("Sample_Name", "Sample_ID", "Sentrix_ID"), names(ss))
      have <- intersect(join_cols, names(ss))
      if (length(key) && length(have) && "Sample_Name" %in% names(si)) {
        idx <- match(si$Sample_Name, ss[[key[1]]])
        for (cn in have) si[[cn]] <- ss[[cn]][idx]
      }
    }
    out[[p]] <- si
  }

  # rbind tolerating differing column sets across platforms.
  all_cols <- unique(unlist(lapply(out, names), use.names = FALSE))
  out <- lapply(out, function(d) {
    for (cn in setdiff(all_cols, names(d))) d[[cn]] <- NA
    d[, all_cols, drop = FALSE]
  })
  merged <- do.call(rbind, out)
  rownames(merged) <- NULL

  if ("Sample_ID" %in% names(merged) && anyDuplicated(merged$Sample_ID)) {
    warning("Duplicate Sample_ID values across platforms in merged sample_info.")
  }
  message(sprintf("merge_platform_sample_info: %d samples across %d platforms.",
                  nrow(merged), nrow(src)))
  merged
}

#' Concatenate per-platform CNV segment files.
#'
#' Segments are genomic intervals, so they are platform-independent once called
#' — this is a plain row concatenation, no probe harmonisation involved. Each
#' platform must therefore be segmented by its own conumee2 annotation first
#' (cnv_analysis.R builds one annotation per run).
#'
#' @param sources     Named vector or data frame of platform -> run directory.
#' @param output_file Optional path to write the merged .seg to.
#' @return Data frame of all segments, with a Platform column.
merge_platform_segments <- function(sources, output_file = NULL) {
  src <- .as_sources(sources)
  segs <- list()

  for (i in seq_len(nrow(src))) {
    p    <- src$platform[i]
    path <- file.path(src$dir[i], "cnv", "segments", "combined_cnv_segments.seg")
    if (!file.exists(path)) {
      warning("No combined_cnv_segments.seg for platform ", p, " — skipping.")
      next
    }
    df <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(df) == 0) {
      warning("Empty segment file for platform ", p, " — skipping.")
      next
    }
    df$Platform <- p
    segs[[p]] <- df
    message(sprintf("merge_platform_segments: %s -> %d segments.", p, nrow(df)))
  }

  if (length(segs) == 0) {
    warning("No segment files found for any platform; nothing merged.")
    return(NULL)
  }
  if (length(segs) < nrow(src)) {
    warning(sprintf("Merged segments cover %d of %d platforms.",
                    length(segs), nrow(src)))
  }

  common_cols <- Reduce(intersect, lapply(segs, names))
  if (length(common_cols) <= 1) {
    stop("Segment files share no common columns — check .seg formats.")
  }
  dropped <- setdiff(unique(unlist(lapply(segs, names), use.names = FALSE)),
                     common_cols)
  if (length(dropped)) {
    message("  note: dropping non-shared segment column(s): ",
            paste(dropped, collapse = ", "))
  }
  merged <- do.call(rbind, lapply(segs, function(d) d[, common_cols, drop = FALSE]))
  rownames(merged) <- NULL

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(merged, output_file, sep = "\t", quote = FALSE,
                       row.names = FALSE)
    message(sprintf("merge_platform_segments: wrote %d segments -> %s",
                    nrow(merged), output_file))
  }
  merged
}

#' Merge beta, sample_info and segments from per-platform runs into one cohort.
#'
#' Writes a merged/ tree mirroring a normal run's processed_data + cnv layout,
#' so the existing dim-reduction and reference-projection steps can consume it
#' unchanged. Deliberately does NOT merge RGChannelSets: minfi::combineArrays
#' down-converts to the lowest common probe space and does not support EPICv2
#' alongside EPIC/450k, and none of dim reduction, reference projection, or CNV
#' plotting consumes an RGSet post-segmentation.
#'
#' @param sources      Named vector or data frame of platform -> run directory.
#' @param output_dir   Destination for the merged cohort.
#' @param policy       Probe reconciliation policy; "intersection".
#' @param which        "raw" or "filtered" beta matrices.
#' @param samplesheets Optional platform -> samplesheet CSV, to lift Group etc.
#' @param compress     Gzip the merged beta matrix (default TRUE).
#' @return List with beta, sample_info, segments and the paths written.
merge_platforms <- function(sources, output_dir, policy = "intersection",
                            which = "raw", samplesheets = NULL,
                            compress = TRUE) {
  src <- .as_sources(sources)
  message("=== Merging ", nrow(src), " platforms: ",
          paste(src$platform, collapse = ", "), " ===")

  proc_dir <- file.path(output_dir, "processed_data")
  seg_dir  <- file.path(output_dir, "cnv", "segments")
  dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(seg_dir,  recursive = TRUE, showWarnings = FALSE)

  bm <- merge_platform_betas(src, policy = policy, which = which)
  si <- merge_platform_sample_info(src, samplesheets = samplesheets)

  # Align sample_info to the beta column order; carry Platform from the merge
  # so the two can never disagree.
  if ("Sample_ID" %in% names(si)) {
    ord <- match(colnames(bm$beta), si$Sample_ID)
    if (anyNA(ord)) {
      warning(sprintf("%d beta column(s) absent from sample_info: %s",
                      sum(is.na(ord)),
                      paste(utils::head(colnames(bm$beta)[is.na(ord)], 5),
                            collapse = ", ")))
    }
    si <- si[ord[!is.na(ord)], , drop = FALSE]
    rownames(si) <- NULL
  }
  si$Platform <- unname(bm$platform_of_sample[
    match(si$Sample_ID, names(bm$platform_of_sample))])

  beta_path <- file.path(proc_dir,
                         if (compress) "beta_values.txt.gz" else "beta_values.txt")
  con <- if (compress) gzfile(beta_path, "w") else file(beta_path, "w")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  utils::write.table(
    data.frame(ProbeID = rownames(bm$beta), bm$beta, check.names = FALSE),
    con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  on.exit(NULL)
  message("merge_platforms: wrote ", beta_path)

  si_path <- file.path(proc_dir, "sample_info.txt")
  utils::write.table(si, si_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("merge_platforms: wrote ", si_path)

  segs <- merge_platform_segments(
    src, output_file = file.path(seg_dir, "combined_cnv_segments.seg"))

  summary_df <- data.frame(
    platform = src$platform,
    n_samples = as.integer(table(bm$platform_of_sample)[src$platform]),
    n_probes  = as.integer(bm$probes_per_platform[src$platform]),
    stringsAsFactors = FALSE)
  summary_df$n_common_probes <- bm$n_common
  utils::write.csv(summary_df, file.path(output_dir, "merge_summary.csv"),
                   row.names = FALSE)

  message("=== Merge complete: ", bm$n_common, " probes x ",
          ncol(bm$beta), " samples ===")
  invisible(list(beta = bm$beta, sample_info = si, segments = segs,
                 summary = summary_df, beta_path = beta_path))
}

# ---------------------------------------------------------------------------
# Platform batch correction
# ---------------------------------------------------------------------------

#' Beta <-> M-value conversion.
#'
#' Beta is bounded [0,1] and strongly heteroscedastic (variance collapses at
#' both ends), which violates ComBat's assumption of roughly homoscedastic
#' normal errors. M-values are the standard remedy. Values are clamped away
#' from 0/1 so the logit stays finite.
#'
#' @param b   Beta matrix.
#' @param eps Clamp distance from 0 and 1.
#' @return M-value matrix.
beta_to_m <- function(b, eps = 1e-4) {
  b <- pmin(pmax(as.matrix(b), eps), 1 - eps)
  log2(b / (1 - b))
}

#' @rdname beta_to_m
m_to_beta <- function(m) 2^m / (1 + 2^m)

#' Remove array-platform batch effects from a merged beta matrix.
#'
#' Uses ComBat (empirical-Bayes location/scale adjustment) on M-values. ComBat
#' is preferred over limma::removeBatchEffect here because per-probe variance
#' estimates are poor at the small per-platform sample sizes typical of a mixed
#' cohort; ComBat shrinks those estimates toward a common prior, which
#' removeBatchEffect does not do.
#'
#' Pass biological covariates you want protected via \code{covariates} — ComBat
#' will not absorb them into the batch term. Omit it when the cohort is a single
#' biological group. A covariate collinear with platform is not identifiable and
#' will make ComBat fail or mis-adjust; the design must be at least partly
#' balanced across platforms.
#'
#' Probes containing any NA, and probes with zero within-platform variance, are
#' dropped — ComBat handles neither.
#'
#' @param beta       Merged beta matrix, probes x samples.
#' @param platform   Platform label per sample (length = ncol(beta)).
#' @param covariates Optional data frame of covariates to preserve.
#' @return List: beta (corrected), n_probes_in, n_probes_out, dropped_na,
#'         dropped_invariant.
correct_platform_batch <- function(beta, platform, covariates = NULL) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for ComBat batch correction.")
  }
  beta <- as.matrix(beta)
  platform <- as.character(platform)
  if (length(platform) != ncol(beta)) {
    stop("`platform` length (", length(platform),
         ") must equal ncol(beta) (", ncol(beta), ").")
  }
  tb <- table(platform)
  if (length(tb) < 2) {
    stop("Batch correction needs >= 2 platforms; got ", length(tb), ".")
  }
  if (any(tb < 2)) {
    stop("Every platform needs >= 2 samples for ComBat; got: ",
         paste(sprintf("%s=%d", names(tb), tb), collapse = ", "))
  }

  n_in    <- nrow(beta)
  na_rows <- apply(beta, 1, anyNA)
  beta    <- beta[!na_rows, , drop = FALSE]

  # Zero within-platform SD breaks ComBat's scale term.
  min_sd <- apply(beta, 1, function(x) min(tapply(x, platform, stats::sd)))
  inv    <- is.na(min_sd) | min_sd <= 1e-8
  beta   <- beta[!inv, , drop = FALSE]

  if (nrow(beta) == 0) stop("No usable probes left for batch correction.")
  message(sprintf(
    "correct_platform_batch: %d probes in; dropped %d with NA, %d within-platform-invariant; %d remain.",
    n_in, sum(na_rows), sum(inv), nrow(beta)))

  mod <- if (is.null(covariates)) {
    stats::model.matrix(~ 1, data = data.frame(row.names = colnames(beta)))
  } else {
    cv <- as.data.frame(covariates)
    if (nrow(cv) != ncol(beta)) stop("`covariates` must have one row per sample.")
    stats::model.matrix(~ ., data = cv)
  }

  corrected <- sva::ComBat(dat = beta_to_m(beta), batch = platform, mod = mod,
                           par.prior = TRUE, prior.plots = FALSE)

  list(beta = m_to_beta(corrected), n_probes_in = n_in,
       n_probes_out = nrow(beta), dropped_na = sum(na_rows),
       dropped_invariant = sum(inv))
}

#' Quantify how strongly platform structures a merged cohort.
#'
#' Two complementary read-outs, computed on the top-variance probes:
#'   * ANOVA of PC1/PC2 against platform — a *small* p means platform drives
#'     the dominant axes of variation (bad).
#'   * mean within- minus between-platform sample correlation — a *positive*
#'     gap means samples resemble their platform-mates more than anything else
#'     (bad). A successful correction drives the gap toward zero or negative.
#'
#' \code{control} is an optional known biological grouping (e.g. predicted sex)
#' used as a negative control: correction should leave it detectable. Without
#' such a control, "platform separation disappeared" cannot be distinguished
#' from "all between-sample structure was flattened".
#'
#' @param beta     Beta matrix, probes x samples.
#' @param platform Platform label per sample.
#' @param control  Optional biological label per sample to check is preserved.
#' @param n_probes Number of top-variance probes to use.
#' @return List of diagnostics (also messaged).
platform_batch_diagnostics <- function(beta, platform, control = NULL,
                                       n_probes = 10000) {
  beta <- as.matrix(beta)
  beta <- beta[!apply(beta, 1, anyNA), , drop = FALSE]
  if (nrow(beta) < 2 || ncol(beta) < 3) {
    stop("Too few probes or samples for diagnostics.")
  }
  v   <- apply(beta, 1, stats::sd)
  top <- beta[order(v, decreasing = TRUE)[seq_len(min(n_probes, nrow(beta)))], ,
              drop = FALSE]

  pca <- stats::prcomp(t(top), center = TRUE, scale. = FALSE)
  pv  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  p_of <- function(g, i) tryCatch(
    summary(stats::aov(pca$x[, i] ~ factor(g)))[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_)
  p_plat <- c(p_of(platform, 1), p_of(platform, 2))
  p_ctrl <- if (is.null(control)) c(NA_real_, NA_real_)
            else c(p_of(control, 1), p_of(control, 2))

  cm   <- stats::cor(top, use = "pairwise.complete.obs")
  same <- outer(platform, platform, "==")
  diag(same) <- NA
  within  <- mean(cm[which(same)],  na.rm = TRUE)
  between <- mean(cm[which(!same)], na.rm = TRUE)

  message(sprintf(
    "platform_batch_diagnostics: PC1 %.1f%% PC2 %.1f%% | platform p: PC1=%.3g PC2=%.3g | corr within=%.4f between=%.4f (gap %+.4f)",
    pv[1], pv[2], p_plat[1], p_plat[2], within, between, within - between))
  if (!is.null(control)) {
    message(sprintf("  control signal p: PC1=%.3g PC2=%.3g (small = preserved)",
                    p_ctrl[1], p_ctrl[2]))
  }

  list(pc_var = pv[1:2], p_platform = p_plat, p_control = p_ctrl,
       corr_within = within, corr_between = between,
       corr_gap = within - between, pca = pca)
}
