# reference_projection.R
# ---------------------------------------------------------------------------
# Reference-projection module (MeQTrack v2.0.0).
#
# Projects user methylation samples onto a pre-trained reference t-SNE
# embedding, so a query sample is positioned against thousands of labelled
# reference methylomes instead of being clustered only against itself.
#
# Projection engine: snifter::project() (an R wrapper around openTSNE).
# snifter::project(x, new, old) recomputes the affinities of `old` and
# checks them against the embedding's stored affinities — so `old` MUST be
# the *complete, unmodified* reference beta matrix the embedding was built
# on. Probe harmonisation therefore happens on the QUERY side only: the
# query is reshaped to exactly the reference probe set (same probes, same
# order), imputing any it lacks with the reference per-probe mean.
#
# The reference embedding was built from SWAN-normalised betas, so the
# query is preprocessed with preprocessSWAN() to match.
#
# Functions:
#   load_reference(dataset, reference_dir)
#   preprocess_query_swan(samplesheet)
#   harmonize_query(query_beta, ref_beta)
#   project_onto_reference(reference, query_beta_harmonized, perplexity)
#   plot_reference_projection(reference, projected, output_dir, file)
#   run_reference_projection(samplesheet, dataset, reference_dir, output_dir, perplexity)
# ---------------------------------------------------------------------------

# Registry of available reference datasets. Each entry names the files under
# reference_dir and the metadata columns carrying the Sentrix ID, the
# tumour-class label, and the per-class plotting colour. v2.0.0 shipped the
# COMET primary-diagnostic set; v2.1 adds the Capper et al. CNS-tumour
# reference (GSE90496) and the Koelsche et al. sarcoma reference (GSE140686).
#
# The reference beta matrix ships pre-built as a compact .rds (`beta_rds`) —
# the canonical format for every dataset (smaller and far faster to load
# than a CSV). `beta_csv` is an optional on-ramp: if only a CSV is present,
# load_reference() parses it once and caches the .rds. Both shipped datasets
# set beta_csv = NA, so the .rds is the single source of truth.
#
# A query can be projected onto any of these — the app's Settings panel
# exposes the choice (settings_module.R reads this registry directly so the
# two never drift), and load_reference(dataset) resolves a key to its files.
.REFERENCE_DATASETS <- list(
  COMET_1915 = list(
    label       = "Pour et al. COMET pediatric solid-tumor (1915 samples, GSE305405)",
    embedding   = "tSNE_embedding_GSE305405_top10K.RData",
    beta_rds    = "beta_GSE305405_1915samples_top10K.rds",
    beta_csv    = NA_character_,
    metadata    = "COMET_Labkey_August_12_2025.csv",
    sentrix_col = "X850k Tumor File Name",
    class_col   = "Tumor Group For Clustering",
    color_col   = "Col"
  ),

  # Capper et al. (2018) CNS-tumour methylation reference, GSE90496 — 2801
  # samples across 91 methylation classes. The metadata is per-sample (one
  # row per GSM sample): GSE90496_MC_MCF_color_labels_key.csv carries the
  # methylation class (`meth.class`) and its colour (`col`). The sibling
  # GSE90496_color_key_allRef.csv is only a class -> colour lookup, not
  # per-sample, so load_reference() (which matches the embedding's samples
  # to metadata rows) cannot consume it directly — and it need not, since
  # the per-sample `col` values already match that key exactly.
  Capper_GSE90496 = list(
    label       = "Capper et al. CNS-tumor reference (2801 samples, GSE90496)",
    embedding   = "tSNE_embedding_GSE90496_top10K.RData",
    beta_rds    = "beta_GSE90496_top10K.rds",
    beta_csv    = NA_character_,
    metadata    = "GSE90496_MC_MCF_color_labels_key.csv",
    sentrix_col = "Sample",
    class_col   = "meth.class",
    color_col   = "col"
  ),

  # Koelsche et al. (2021) sarcoma methylation reference, GSE140686 — 1077
  # tumours across 65 sarcoma methylation classes. The metadata is per-sample
  # (one row per reference tumour): `IDAT` carries the Sentrix ID matching
  # the embedding, `Meth_Class` the short class code (e.g. RMS-ALV), and
  # `col` its hex plotting colour. The sibling `Methylation_Class_Name`
  # column holds the same classes spelled out in full.
  Sarcoma_GSE140686 = list(
    label       = "Koelsche et al. sarcoma reference (1077 samples, GSE140686)",
    embedding   = "tSNE_embedding_GSE140686_top10K.RData",
    beta_rds    = "beta_GSE140686_1077Sarcoma_top10K.rds",
    beta_csv    = NA_character_,
    metadata    = "GSE140686_sarcoma_methylation_labels.csv",
    sentrix_col = "IDAT",
    class_col   = "Meth_Class",
    color_col   = "col"
  )
)

# ---------------------------------------------------------------------------
# Load a reference dataset
# ---------------------------------------------------------------------------

#' Load a reference embedding, its training beta matrix, and its metadata.
#'
#' @param dataset       Key into \code{.REFERENCE_DATASETS} (default "COMET_1915").
#' @param reference_dir Directory holding the reference files.
#' @return A list with: dataset, label, embedding (snifter object),
#'   ref_beta (probes x samples, aligned to the embedding), and ref_meta
#'   (data frame: Sentrix, tSNE1, tSNE2, tumor_group, color).
load_reference <- function(dataset = "COMET_1915", reference_dir = "reference") {
  spec <- .REFERENCE_DATASETS[[dataset]]
  if (is.null(spec)) {
    stop("Unknown reference dataset '", dataset, "'. Known: ",
         paste(names(.REFERENCE_DATASETS), collapse = ", "))
  }

  # --- embedding (a single snifter object inside the .RData) ---
  emb_path <- file.path(reference_dir, spec$embedding)
  if (!file.exists(emb_path)) stop("Reference embedding not found: ", emb_path)
  ee <- new.env()
  load(emb_path, envir = ee)
  embedding <- get(ls(ee)[1], envir = ee)
  if (is.null(rownames(embedding))) {
    stop("Reference embedding '", dataset, "' has no rownames — cannot map ",
         "rows to samples.")
  }

  # --- reference beta matrix ---
  # The .rds is the canonical format and is what every shipped dataset
  # provides. beta_csv is an optional fallback (NA when absent): when a CSV
  # is present but the .rds is not, it is parsed once and cached as the .rds
  # for every subsequent load.
  rds     <- file.path(reference_dir, spec$beta_rds)
  has_csv <- !is.null(spec$beta_csv) && !is.na(spec$beta_csv)
  csv     <- if (has_csv) file.path(reference_dir, spec$beta_csv) else NA_character_
  if (file.exists(rds)) {
    ref_beta <- readRDS(rds)
  } else if (has_csv && file.exists(csv)) {
    message("load_reference: parsing reference beta CSV (slow) and caching as .rds ...")
    suppressPackageStartupMessages(library(data.table))
    d <- as.data.frame(data.table::fread(csv))
    rownames(d) <- d[[1]]
    ref_beta <- as.matrix(d[, -1, drop = FALSE])
    saveRDS(ref_beta, rds)
  } else {
    stop("Reference beta matrix not found (looked for ", rds,
         if (has_csv) paste(" and", csv) else "", ")")
  }
  ref_beta <- as.matrix(ref_beta)

  # Align beta columns to the embedding's sample order.
  if (!identical(colnames(ref_beta), rownames(embedding))) {
    if (!setequal(colnames(ref_beta), rownames(embedding))) {
      stop("Reference beta samples do not match the embedding samples.")
    }
    ref_beta <- ref_beta[, rownames(embedding), drop = FALSE]
  }

  # --- metadata: pull Sentrix / tumour-group / colour, aligned to embedding ---
  md_path <- file.path(reference_dir, spec$metadata)
  if (!file.exists(md_path)) stop("Reference metadata not found: ", md_path)
  suppressPackageStartupMessages(library(data.table))
  md <- as.data.frame(data.table::fread(md_path))
  for (col in c(spec$sentrix_col, spec$class_col, spec$color_col)) {
    if (!col %in% colnames(md)) {
      stop("Metadata is missing expected column '", col, "'.")
    }
  }
  idx <- match(rownames(embedding), as.character(md[[spec$sentrix_col]]))
  if (anyNA(idx)) {
    warning(sprintf("load_reference: %d embedding samples have no metadata row.",
                    sum(is.na(idx))))
  }
  ref_meta <- data.frame(
    Sentrix     = rownames(embedding),
    tSNE1       = as.numeric(embedding[, 1]),
    tSNE2       = as.numeric(embedding[, 2]),
    tumor_group = as.character(md[[spec$class_col]])[idx],
    color       = as.character(md[[spec$color_col]])[idx],
    stringsAsFactors = FALSE
  )
  ref_meta$tumor_group[is.na(ref_meta$tumor_group) | ref_meta$tumor_group == ""] <- "Unknown"

  message(sprintf("load_reference: '%s' — %d reference samples, %d probes, %d tumour groups.",
                  dataset, nrow(embedding), nrow(ref_beta),
                  length(unique(ref_meta$tumor_group))))

  list(dataset = dataset, label = spec$label, embedding = embedding,
       ref_beta = ref_beta, ref_meta = ref_meta)
}

# ---------------------------------------------------------------------------
# Preprocess a query (IDATs -> SWAN-normalised betas)
# ---------------------------------------------------------------------------

#' Read query IDATs and SWAN-normalise, to match the reference's preprocessing.
#'
#' @param samplesheet A data frame with a \code{Basename} column, a path to
#'   such a CSV, or a directory of IDAT files.
#' @return Beta matrix, probes x samples.
preprocess_query_swan <- function(samplesheet) {
  suppressPackageStartupMessages(library(minfi))

  if (is.character(samplesheet) && length(samplesheet) == 1L && dir.exists(samplesheet)) {
    rgset <- read.metharray.exp(base = samplesheet, recursive = TRUE, force = TRUE)
  } else {
    if (is.character(samplesheet)) {
      suppressPackageStartupMessages(library(data.table))
      ss <- as.data.frame(data.table::fread(samplesheet))
    } else {
      ss <- as.data.frame(samplesheet)
    }
    if (!"Basename" %in% colnames(ss)) {
      stop("Query samplesheet needs a 'Basename' column.")
    }
    rgset <- read.metharray.exp(targets = ss, recursive = TRUE, force = TRUE)
  }

  message(sprintf("preprocess_query_swan: %d query sample(s); running SWAN ...",
                  ncol(rgset)))
  mset <- preprocessSWAN(rgset)
  getBeta(mset)
}

# ---------------------------------------------------------------------------
# Harmonise a query beta matrix to the reference probe set
# ---------------------------------------------------------------------------

#' Normalise query probe IDs to plain CpG IDs.
#'
#' EPICv2 manifest probe IDs carry a design suffix (e.g. cg00000029_TC21),
#' while the reference embeddings were built on plain 450K/EPIC IDs
#' (cg00000029). Without this step an EPICv2 query intersects almost none
#' of the reference probes and collapses to the reference mean. We strip
#' the suffix and, where EPICv2 ships several replicate probes for one CpG,
#' average them per CpG. A no-op for 450K / EPIC v1 (their IDs have no
#' suffix), so it is safe to call unconditionally.
#'
#' @param beta Query beta matrix, probes x samples.
#' @return Beta matrix with plain-CpG rownames; replicate CpGs collapsed.
.normalize_probe_ids <- function(beta) {
  beta <- as.matrix(beta)
  pid  <- rownames(beta)
  base <- sub("_.*$", "", pid)
  if (identical(base, pid)) return(beta)            # no suffixes — 450K / EPIC v1

  n_suffixed <- sum(base != pid)
  if (anyDuplicated(base) == 0L) {
    rownames(beta) <- base
    message(sprintf(
      "harmonize_query: EPICv2-style probe IDs — stripped suffix from %d probe(s).",
      n_suffixed))
    return(beta)
  }
  # EPICv2 ships replicate probes for some CpGs — average them per CpG.
  sums   <- rowsum(beta, base, na.rm = TRUE)
  counts <- rowsum(1 * !is.na(beta), base)
  counts[counts == 0] <- NA_real_                   # all-NA group -> NA, not 0/0
  out <- sums / counts
  message(sprintf(
    "harmonize_query: EPICv2-style probe IDs — stripped suffix from %d probe(s); collapsed %d replicate row(s) into per-CpG means (%d -> %d probes).",
    n_suffixed, nrow(beta) - nrow(out), nrow(beta), nrow(out)))
  out
}

#' Reshape a query beta matrix to exactly the reference probe set.
#'
#' snifter::project() requires \code{ncol(new) == ncol(old)} with matching
#' features, and \code{old} cannot be modified — so the query must carry
#' every reference probe, in the reference's order. Probes the query lacks,
#' and any residual NA values, are imputed with the reference per-probe mean
#' (a neutral choice that does not pull the sample toward any class).
#' Query probe IDs are first normalised by \code{.normalize_probe_ids()} so
#' EPICv2 arrays match the 450K/EPIC-style reference probe names.
#'
#' @param query_beta Query beta matrix, probes x samples.
#' @param ref_beta   Reference beta matrix, probes x samples.
#' @return Query beta matrix, reference-probes x query-samples.
harmonize_query <- function(query_beta, ref_beta) {
  query_beta <- .normalize_probe_ids(as.matrix(query_beta))
  ref_probes <- rownames(ref_beta)

  common  <- intersect(ref_probes, rownames(query_beta))
  pct     <- 100 * length(common) / length(ref_probes)
  n_imput <- length(ref_probes) - length(common)
  message(sprintf("harmonize_query: %d/%d reference probes found in query (%.1f%%); %d probe(s) imputed.",
                  length(common), length(ref_probes), pct, n_imput))
  if (pct < 80) {
    warning(sprintf("harmonize_query: only %.1f%% of reference probes present — ",
                    pct), "projection may be unreliable.")
  }

  ref_mean <- rowMeans(ref_beta, na.rm = TRUE)

  # Start every probe at the reference mean, then overwrite with query values.
  out <- matrix(ref_mean, nrow = length(ref_probes), ncol = ncol(query_beta),
                dimnames = list(ref_probes, colnames(query_beta)))
  out[common, ] <- query_beta[common, , drop = FALSE]

  # Replace residual NA (failed probes) with the reference mean.
  na_idx <- which(is.na(out), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    out[na_idx] <- ref_mean[na_idx[, 1]]
    message(sprintf("harmonize_query: imputed %d residual NA value(s).", nrow(na_idx)))
  }
  out
}

# ---------------------------------------------------------------------------
# Project the query onto the reference embedding
# ---------------------------------------------------------------------------

#' Project a harmonised query onto a reference t-SNE embedding.
#'
#' @param reference              A list from \code{load_reference()}.
#' @param query_beta_harmonized  Output of \code{harmonize_query()}.
#' @param perplexity             Perplexity for the projection's own kNN
#'                               (snifter::project default = 5).
#' @return Data frame: Sample, tSNE1, tSNE2.
project_onto_reference <- function(reference, query_beta_harmonized,
                                   perplexity = 5) {
  suppressPackageStartupMessages(library(snifter))

  # snifter expects samples x features for both `new` and `old`.
  old <- t(as.matrix(reference$ref_beta))     # reference samples x probes
  new <- t(as.matrix(query_beta_harmonized))  # query samples x probes

  stopifnot(identical(colnames(old), colnames(new)))   # same probes, same order
  stopifnot(nrow(old) == nrow(reference$embedding))    # old aligns to embedding

  message(sprintf("project_onto_reference: projecting %d query sample(s) onto '%s' ...",
                  nrow(new), reference$dataset))
  coords <- snifter::project(reference$embedding, new = new, old = old,
                             perplexity = perplexity)
  coords <- as.matrix(coords)[, 1:2, drop = FALSE]

  data.frame(
    Sample = rownames(new),
    tSNE1  = as.numeric(coords[, 1]),
    tSNE2  = as.numeric(coords[, 2]),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Plot the projection
# ---------------------------------------------------------------------------

# Replace any colour string R can't parse with a generated fallback.
.safe_color_map <- function(col_map) {
  ok <- vapply(col_map, function(cc) {
    !is.na(cc) && tryCatch({ grDevices::col2rgb(cc); TRUE },
                           error = function(e) FALSE)
  }, logical(1))
  if (any(!ok)) col_map[!ok] <- grDevices::rainbow(sum(!ok))
  col_map
}

#' Plot the reference embedding with the projected query samples overlaid.
#'
#' @param reference  A list from \code{load_reference()}.
#' @param projected  Data frame from \code{project_onto_reference()}.
#' @param output_dir Directory for the PDF.
#' @param file       Optional explicit output path.
#' @return The output PDF path (invisibly).
plot_reference_projection <- function(reference, projected, output_dir = ".",
                                      file = NULL) {
  suppressPackageStartupMessages(library(ggplot2))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (is.null(file)) {
    file <- file.path(output_dir,
                      sprintf("reference_projection_%s.pdf", reference$dataset))
  }

  rm_ <- reference$ref_meta
  rm_$tumor_group <- factor(rm_$tumor_group)
  col_map <- tapply(rm_$color, rm_$tumor_group, function(x) x[1])
  col_map <- .safe_color_map(col_map)

  # A long tumour-group legend wraps into several columns and eats into the
  # plot panel, squishing the t-SNE map. COMET (16 groups) fits one column;
  # GSE140686 (65) and GSE90496 (91) need more. Cap each legend column at
  # ~32 entries, then widen the page by one column's worth for each legend
  # column beyond the first so the map keeps a stable, roughly-square size.
  # Multi-column datasets also get a taller page (10 in vs 8) so the wider
  # canvas stays proportional rather than letterboxed.
  n_groups    <- nlevels(rm_$tumor_group)
  legend_cols <- max(1L, ceiling(n_groups / 32))
  plot_width  <- 10 + 3 * (legend_cols - 1L)
  plot_height <- if (legend_cols > 1L) 10 else 8

  p <- ggplot() +
    geom_point(data = rm_,
               aes(x = tSNE1, y = tSNE2, fill = tumor_group),
               shape = 21, size = 3.2, stroke = 0.1,
               colour = "grey50", alpha = 0.8) +
    scale_fill_manual(values = col_map, name = "Reference\ntumour group") +
    guides(fill = guide_legend(ncol = legend_cols)) +
    geom_point(data = projected,
               aes(x = tSNE1, y = tSNE2),
               shape = 23, size = 2.5, stroke = 0.6,
               fill = "#111111", colour = "white") +
    # No per-sample text labels — they crowd the plot. Sample names are
    # available on hover in the interactive (HTML/Shiny) projection view.
    labs(title = sprintf("Query projected onto %s", reference$label),
         subtitle = sprintf("%d query sample(s); black diamonds = query, coloured = reference",
                            nrow(projected)),
         x = "t-SNE 1", y = "t-SNE 2") +
    theme_bw(base_size = 12) +
    theme(legend.text = element_text(size = rel(0.7)),
          panel.grid.minor = element_blank())

  ggsave(file, p, width = plot_width, height = plot_height)
  message("plot_reference_projection: wrote ", file)
  invisible(file)
}

# ---------------------------------------------------------------------------
# Nearest-class diagnostic
# ---------------------------------------------------------------------------

#' Assign each projected query sample a nearest reference tumour class by
#' k-NN vote in the embedding space.
#'
#' For each query point the k nearest reference points (Euclidean distance
#' in the 2-D t-SNE space) vote with their tumour group. The winning group
#' is the class hint and the vote fraction is the confidence. A sample is
#' flagged `ambiguous` when no group wins a majority, or the top two are
#' within 15 percentage points — the between-classes case to surface
#' rather than forcing a single label. `distant_from_reference` flags a
#' sample that lands far from any reference cluster (its methylome is
#' unlike anything in the reference), where the class hint is unreliable.
#'
#' @param projected  Data frame from project_onto_reference()
#'                   (Sample, tSNE1, tSNE2).
#' @param ref_meta   Reference metadata (tSNE1, tSNE2, tumor_group).
#' @param k          Number of nearest reference neighbours (default 25).
#' @param top_n      Classes to list in the `top_classes` summary string.
#' @return Data frame, one row per query sample.
nearest_reference_class <- function(projected, ref_meta, k = 25, top_n = 3) {
  ref_xy  <- as.matrix(ref_meta[, c("tSNE1", "tSNE2")])
  ref_grp <- as.character(ref_meta$tumor_group)
  k       <- max(1L, min(as.integer(k), nrow(ref_xy) - 1L))

  # Typical reference neighbourhood spread (mean distance to the k nearest
  # reference points), from a subsample for speed — the scale against which
  # a query is judged "distant from the reference".
  set.seed(1L)
  sidx <- sample(nrow(ref_xy), min(300L, nrow(ref_xy)))
  ref_spread <- stats::median(vapply(sidx, function(j) {
    d <- sqrt((ref_xy[, 1] - ref_xy[j, 1])^2 + (ref_xy[, 2] - ref_xy[j, 2])^2)
    mean(sort(d)[2:(k + 1L)])          # [2:] excludes the point itself
  }, numeric(1)))

  rows <- lapply(seq_len(nrow(projected)), function(i) {
    qx  <- projected$tSNE1[i]; qy <- projected$tSNE2[i]
    d   <- sqrt((ref_xy[, 1] - qx)^2 + (ref_xy[, 2] - qy)^2)
    ord <- order(d)[seq_len(k)]
    votes <- sort(table(ref_grp[ord]), decreasing = TRUE)
    cls   <- names(votes)
    frac  <- as.numeric(votes) / k

    runner_conf <- if (length(cls) >= 2L) frac[2] else 0
    n_top       <- min(top_n, length(cls))
    data.frame(
      Sample                 = projected$Sample[i],
      nearest_class          = cls[1],
      confidence             = round(frac[1], 3),
      runner_up_class        = if (length(cls) >= 2L) cls[2] else NA_character_,
      runner_up_confidence   = round(runner_conf, 3),
      ambiguous              = frac[1] < 0.5 || (frac[1] - runner_conf) < 0.15,
      distant_from_reference = mean(d[ord]) > 2.5 * ref_spread,
      mean_knn_distance      = round(mean(d[ord]), 3),
      top_classes            = paste(sprintf("%s (%.0f%%)",
                                             cls[seq_len(n_top)],
                                             100 * frac[seq_len(n_top)]),
                                     collapse = ", "),
      stringsAsFactors       = FALSE
    )
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

#' Run the full reference projection: query IDATs -> projected coordinates.
#'
#' @param samplesheet   Query samplesheet (data frame / CSV path / IDAT dir).
#' @param dataset       Reference dataset key (default "COMET_1915").
#' @param reference_dir Directory holding the reference files.
#' @param output_dir    Directory for the coordinate CSV and plot PDF.
#' @param perplexity    Projection perplexity (default 5).
#' @param knn_k         Neighbours for the nearest-class diagnostic (default 25).
#' @return A list: dataset, label, projected, class_hints, ref_meta, csv,
#'   pdf, hints_csv.
run_reference_projection <- function(samplesheet,
                                     dataset       = "COMET_1915",
                                     reference_dir = "reference",
                                     output_dir    = ".",
                                     perplexity    = 5,
                                     knn_k         = 25) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  reference  <- load_reference(dataset, reference_dir)
  query_beta <- preprocess_query_swan(samplesheet)
  query_h    <- harmonize_query(query_beta, reference$ref_beta)
  projected  <- project_onto_reference(reference, query_h, perplexity)

  # Nearest-class diagnostic — k-NN vote in the embedding space.
  class_hints <- nearest_reference_class(projected, reference$ref_meta, k = knn_k)
  for (i in seq_len(nrow(class_hints))) {
    ch <- class_hints[i, ]
    message(sprintf("  %s -> %s (%.0f%% confidence)%s%s",
                    ch$Sample, ch$nearest_class, 100 * ch$confidence,
                    if (isTRUE(ch$ambiguous)) " [ambiguous]" else "",
                    if (isTRUE(ch$distant_from_reference))
                      " [distant from reference]" else ""))
  }

  csv <- file.path(output_dir, sprintf("reference_projection_%s.csv", dataset))
  write.csv(projected, csv, row.names = FALSE)
  hints_csv <- file.path(output_dir,
                         sprintf("reference_projection_class_hints_%s.csv", dataset))
  write.csv(class_hints, hints_csv, row.names = FALSE)
  pdf <- plot_reference_projection(reference, projected, output_dir)

  message(sprintf("run_reference_projection: done — %d sample(s) projected onto '%s'.",
                  nrow(projected), dataset))
  # ref_meta (reference embedding coords + tumour groups + colours) is small
  # (~1900 rows) and is saved alongside the projected query coords + class
  # hints so the Shiny UI can render everything without the 137 MB beta.
  invisible(list(
    dataset     = dataset,
    label       = reference$label,
    projected   = projected,
    class_hints = class_hints,
    ref_meta    = reference$ref_meta,
    csv         = csv,
    pdf         = pdf,
    hints_csv   = hints_csv
  ))
}
