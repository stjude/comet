#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Build a serialized CNV control reference from a directory of control IDATs.
#
# Run once per control set; the pipeline then loads the .rds and never needs
# the IDATs (or the share they live on) again. This mirrors how yamapData
# ships a prepared EPIC reference, but for platforms where we have our own
# controls.
#
#   Rscript scripts/build_cnv_controls.R \
#     --idat_dir /path/to/EPICv2_CNV_Controls \
#     --platform EPICv2 \
#     --output Anno/EPICv2/EPICv2_CNV_controls.rds \
#     [--exclude 207071180012_R01C01] [--min_cor 0.85] [--keep_all]
#
# Quality gate: a CNV reference is the baseline every query sample is fit
# against, so a single noisy control degrades every result. Each control is
# correlated against the median log2 profile of the others and anything below
# --min_cor is dropped, with the reason recorded in the saved metadata.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--idat_dir", type = "character",
              help = "Directory of control IDATs (searched recursively)"),
  make_option("--platform", type = "character", default = "EPICv2",
              help = "Array platform label [default %default]"),
  make_option("--output", type = "character", default = NULL,
              help = "Destination .rds path"),
  make_option("--exclude", type = "character", default = NULL,
              help = "Comma-separated sample IDs to drop before building"),
  make_option("--min_cor", type = "double", default = 0.85,
              help = "Drop controls correlating below this with the others [default %default]"),
  make_option("--keep_all", action = "store_true", default = FALSE,
              help = "Skip the correlation gate and keep every control"),
  make_option("--normalization", type = "character", default = "noob",
              help = "noob | raw [default %default]"),
  make_option("--skip_flatness", action = "store_true", default = FALSE,
              help = "Skip the leave-one-out flatness gate"),
  make_option("--seg_thresh", type = "double", default = 0.20,
              help = "Flatness: |seg.mean| above this counts as altered [default %default]"),
  make_option("--frac_thresh", type = "double", default = 0.02,
              help = "Flatness: altered genome fraction above this drops the control [default %default]"),
  make_option("--flatness_pdf", type = "character", default = NULL,
              help = "Optional PDF of per-control leave-one-out genome plots")
)))

if (is.null(opt$idat_dir) || is.null(opt$output)) {
  stop("--idat_dir and --output are both required.")
}
if (!dir.exists(opt$idat_dir)) {
  stop("IDAT directory not found (is the volume mounted?): ", opt$idat_dir)
}

suppressPackageStartupMessages({
  library(minfi)
  library(conumee2)
})

# Shared leave-one-out flatness scan (also used by check_cnv_controls_flatness.R).
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])),
  "cnv_flatness.R"))

# ---- collect IDAT stems ---------------------------------------------------
# Case-insensitive: Illumina exports vary (_Grn.idat, _Grn.IDAT, .idat.gz).
grn <- list.files(opt$idat_dir, pattern = "_Grn\\.idat(\\.gz)?$",
                  full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
# macOS AppleDouble sidecars (._foo_Grn.idat) look like IDATs but are 4 KB
# resource forks; minfi will happily try to read them and fail.
grn <- grn[!grepl("(^|/)\\._", grn)]
stems <- unique(sub("_Grn\\.idat(\\.gz)?$", "", grn, ignore.case = TRUE))
if (!length(stems)) stop("No IDATs found under ", opt$idat_dir)

if (!is.null(opt$exclude) && nzchar(opt$exclude)) {
  drop <- trimws(strsplit(opt$exclude, ",", fixed = TRUE)[[1]])
  before <- length(stems)
  stems <- stems[!basename(stems) %in% drop]
  message(sprintf("Excluded %d sample(s) by request: %s",
                  before - length(stems), paste(drop, collapse = ", ")))
}
message(sprintf("Found %d control sample(s).", length(stems)))
if (length(stems) < 2) stop("Need at least 2 controls to build a reference.")

# ---- read and normalise ---------------------------------------------------
rgset <- read.metharray.exp(targets = data.frame(Basename = stems,
                                                 stringsAsFactors = FALSE),
                            force = TRUE)
message(sprintf("RGSet: %d probes x %d samples; annotation %s",
                nrow(rgset), ncol(rgset),
                paste(rgset@annotation, collapse = "/")))

mset <- if (identical(opt$normalization, "noob")) preprocessNoob(rgset) else
  preprocessRaw(rgset)
controls <- CNV.load(mset)

# EPICv2 arrays come out of CNV.load with per-replicate rownames
# (cg25324105_BC11, cg25324105_TC11, ...). process_conumee_sample() collapses
# the *query* to plain CpG IDs before CNV.fit, so the reference has to be
# collapsed identically or the two probe spaces will not intersect at all.
raw_ids <- rownames(controls@intensity)
if (any(grepl("_(BC|TC)\\d+$", utils::head(raw_ids, 200)))) {
  message("Collapsing EPICv2 replicate probes in the reference (mean)...")
  base_ids <- sub("_.*$", "", raw_ids)
  int_mat  <- as.matrix(controls@intensity)
  not_na   <- !is.na(int_mat)
  int_mat[!not_na] <- 0
  sums   <- rowsum(int_mat, base_ids)
  counts <- rowsum(matrix(as.numeric(not_na), nrow(int_mat), ncol(int_mat)),
                   base_ids)
  agg <- sums / counts
  agg[counts == 0] <- NA
  controls@intensity <- as.data.frame(agg)
  message(sprintf("  %d -> %d probes after collapse.",
                  length(raw_ids), nrow(controls@intensity)))
}

# ---- quality gate ---------------------------------------------------------
x  <- as.matrix(controls@intensity)
lx <- log2(x + 1)
cor_with_others <- vapply(seq_len(ncol(lx)), function(i) {
  stats::cor(lx[, i], apply(lx[, -i, drop = FALSE], 1, stats::median))
}, numeric(1))
names(cor_with_others) <- colnames(lx)
med_intensity <- apply(x, 2, stats::median, na.rm = TRUE)

message("\nPer-control QC:")
for (s in colnames(lx)) {
  message(sprintf("  %-22s median intensity %8.1f   cor-with-others %.4f",
                  s, med_intensity[s], cor_with_others[s]))
}

dropped <- character(0)
if (!opt$keep_all) {
  bad <- names(cor_with_others)[cor_with_others < opt$min_cor]
  if (length(bad)) {
    if (ncol(x) - length(bad) < 2) {
      stop("Correlation gate would leave fewer than 2 controls; ",
           "inspect the set or pass --keep_all.")
    }
    dropped <- bad
    message(sprintf("\nDropping %d control(s) below --min_cor %.2f: %s",
                    length(bad), opt$min_cor, paste(bad, collapse = ", ")))
    keep <- setdiff(colnames(x), bad)
    controls@intensity <- controls@intensity[, keep, drop = FALSE]
    x <- as.matrix(controls@intensity)
  }
}

# ---- flatness gate --------------------------------------------------------
# The correlation gate above catches technical outliers but is blind to
# copy-number events: a CNV shared across controls (or one large enough to
# survive) still passes, gets baked into the baseline, and then shows up
# inverted in every query. So fit each surviving control against the others
# (leave-one-out), segment, and drop any whose genome is not flat.
flatness <- NULL
nonflat  <- character(0)
if (!opt$skip_flatness) {
  if (ncol(controls@intensity) < 3) {
    message("\nFewer than 3 controls after the correlation gate; skipping the ",
            "flatness scan (leave-one-out needs >= 3).")
  } else {
    message("\nLeave-one-out flatness scan (fits every control; may take a while)...")
    data("exclude_regions"); data("detail_regions")
    anno <- conumee2::CNV.create_anno(
      exclude_regions = exclude_regions, detail_regions = detail_regions,
      bin_minprobes = 15, bin_minsize = 50000, array_type = opt$platform
    )
    common <- intersect(names(anno@probes), rownames(controls@intensity))
    if (!length(common)) stop("No probes shared between annotation and controls.")
    anno@probes <- anno@probes[names(anno@probes) %in% common]
    # Score on a probe-aligned COPY; the saved reference keeps its full probe
    # set (the pipeline intersects it against its own annotation at run time).
    ctrl_aligned <- controls
    ctrl_aligned@intensity <- controls@intensity[common, , drop = FALSE]

    plot_fn <- NULL
    if (!is.null(opt$flatness_pdf)) {
      dir.create(dirname(opt$flatness_pdf), recursive = TRUE, showWarnings = FALSE)
      pdf(opt$flatness_pdf, width = 10, height = 5)
      plot_fn <- function(seg, title) conumee2::CNV.genomeplot(
        seg, main = title,
        cols = c("#018571", "#80cdc1", "#f5f5f5", "#dfc27d", "#a6611a"))
    }
    flatness <- cnv_flatness_scan(ctrl_aligned, anno,
                                  seg_thresh  = opt$seg_thresh,
                                  frac_thresh = opt$frac_thresh,
                                  plot_fn     = plot_fn)
    if (!is.null(opt$flatness_pdf)) {
      dev.off(); message("Per-control genome plots: ", opt$flatness_pdf)
    }

    message("\nFlatness result (worst first):")
    for (j in seq_len(nrow(flatness))) {
      f <- flatness[j, ]
      message(sprintf("  %-22s max|seg.mean| %.3f  altered %5.1f%%  -> %s",
                      f$sample, f$max_abs_segmean, 100 * f$frac_genome_alt,
                      f$verdict))
    }
    nonflat <- flatness$sample[flatness$verdict == "NOT-FLAT"]
    if (length(nonflat)) {
      if (ncol(x) - length(nonflat) < 2)
        stop("Flatness gate would leave fewer than 2 controls; ",
             "inspect the set or pass --skip_flatness.")
      message(sprintf("\nDropping %d non-flat control(s): %s",
                      length(nonflat), paste(nonflat, collapse = ", ")))
      keep <- setdiff(colnames(controls@intensity), nonflat)
      controls@intensity <- controls@intensity[, keep, drop = FALSE]
      x <- as.matrix(controls@intensity)
    } else {
      message("\nAll controls are flat at these thresholds. ✔")
    }
  }
}

overall_median <- stats::median(x, na.rm = TRUE)
message(sprintf("\nFinal reference: %d probes x %d samples, median intensity %.1f",
                nrow(x), ncol(x), overall_median))
if (overall_median < 5000) {
  warning("Median intensity is below conumee2's 5000 threshold; ",
          "CNV calls from this reference may be unreliable.")
}

# ---- save -----------------------------------------------------------------
meta <- list(
  platform      = opt$platform,
  annotation    = paste(rgset@annotation, collapse = "/"),
  n_samples     = ncol(x),
  n_probes      = nrow(x),
  sample_ids    = colnames(x),
  excluded      = c(if (!is.null(opt$exclude)) trimws(strsplit(opt$exclude, ",")[[1]]),
                    dropped, nonflat),
  median_intensity = round(overall_median, 1),
  cor_with_others  = round(cor_with_others, 4),
  flatness         = flatness,          # per-control leave-one-out scan (NULL if skipped)
  flatness_dropped = nonflat,
  flatness_thresholds = c(seg_thresh = opt$seg_thresh, frac_thresh = opt$frac_thresh),
  normalization = opt$normalization,
  source_dir    = normalizePath(opt$idat_dir, mustWork = FALSE),
  built         = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(controls = controls, meta = meta), opt$output, compress = "xz")
message(sprintf("\nWrote %s (%.1f MB)", opt$output,
                file.info(opt$output)$size / 1e6))
