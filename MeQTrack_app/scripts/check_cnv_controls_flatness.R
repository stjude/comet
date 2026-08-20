#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Verify that a set of CNV "control" IDATs are TRUE controls, i.e. each has a
# flat genome (no real copy-number gains/losses).
#
# WHY a separate check: build_cnv_controls.R's only gate is a
# correlation-with-the-other-controls test. That catches technical outliers
# but is BLIND to copy-number events — a CNV shared across controls (or one
# large enough to survive the correlation) still passes and gets baked into
# the reference baseline, which then shows up as a spurious *inverse* event in
# every query sample. Flatness has to be tested directly.
#
# HOW: leave-one-out. Each control is fit (conumee2) against the OTHER
# controls, segmented, and scored on how far its segment log2-ratios stray
# from zero. This mirrors the pipeline's own recipe (Noob -> CNV.load ->
# EPICv2 replicate collapse -> CNV.create_anno -> CNV.fit -> CNV.segment) so
# the verdict reflects what the pipeline would actually do.
#
#   Rscript scripts/check_cnv_controls_flatness.R \
#     --idat_dir /Volumes/qtran/EPICv2_CNV_Controls \
#     --platform EPICv2 \
#     --out_pdf  Anno/EPICv2/EPICv2_controls_flatness.pdf \
#     [--seg_thresh 0.20] [--frac_thresh 0.02] [--exclude 2070..._R01C01]
#
# A control is flagged NOT-FLAT if it has any segment with |seg.mean| above
# --seg_thresh (default 0.20, the pipeline's own loss threshold) covering a
# meaningful stretch, or if the altered fraction of the genome exceeds
# --frac_thresh (default 2%).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--idat_dir",  type = "character", help = "Directory of control IDATs (recursive)"),
  make_option("--platform",  type = "character", default = "EPICv2",
              help = "Array platform for CNV.create_anno [default %default]"),
  make_option("--out_pdf",   type = "character", default = NULL,
              help = "Optional PDF of per-control genome plots"),
  make_option("--exclude",   type = "character", default = NULL,
              help = "Comma-separated sample IDs to drop before checking"),
  make_option("--seg_thresh",  type = "double", default = 0.20,
              help = "|seg.mean| above this counts as an alteration [default %default]"),
  make_option("--frac_thresh", type = "double", default = 0.02,
              help = "Altered genome fraction above this = NOT flat [default %default]"),
  make_option("--normalization", type = "character", default = "noob",
              help = "noob | raw [default %default]")
)))

if (is.null(opt$idat_dir)) stop("--idat_dir is required.")
if (!dir.exists(opt$idat_dir))
  stop("IDAT directory not found (is the volume mounted?): ", opt$idat_dir)

suppressPackageStartupMessages({
  library(minfi)
  library(conumee2)
})

# Shared leave-one-out flatness scan (also used by build_cnv_controls.R).
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])),
  "cnv_flatness.R"))

# ---- collect IDAT stems (same filtering as build_cnv_controls.R) ----------
# Case-insensitive: Illumina exports vary (_Grn.idat, _Grn.IDAT, .idat.gz).
grn <- list.files(opt$idat_dir, pattern = "_Grn\\.idat(\\.gz)?$",
                  full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
grn <- grn[!grepl("(^|/)\\._", grn)]                 # skip macOS AppleDouble
stems <- unique(sub("_Grn\\.idat(\\.gz)?$", "", grn, ignore.case = TRUE))
if (length(stems) < 3)
  stop("Need at least 3 controls for a leave-one-out flatness check; found ", length(stems))

if (!is.null(opt$exclude) && nzchar(opt$exclude)) {
  drop  <- trimws(strsplit(opt$exclude, ",", fixed = TRUE)[[1]])
  stems <- stems[!basename(stems) %in% drop]
  message(sprintf("Excluded %d sample(s) by request.", length(drop)))
}
message(sprintf("Found %d control sample(s).", length(stems)))

# ---- read + normalise + load intensities ----------------------------------
rgset <- read.metharray.exp(targets = data.frame(Basename = stems,
                                                 stringsAsFactors = FALSE),
                            force = TRUE)
message(sprintf("RGSet: %d probes x %d samples; annotation %s",
                nrow(rgset), ncol(rgset), paste(rgset@annotation, collapse = "/")))
mset <- if (identical(opt$normalization, "noob")) preprocessNoob(rgset) else preprocessRaw(rgset)
controls <- CNV.load(mset)

# ---- collapse EPICv2 replicate probes (identical to build + pipeline) ------
raw_ids <- rownames(controls@intensity)
if (any(grepl("_(BC|TC)\\d+$", utils::head(raw_ids, 200)))) {
  message("Collapsing EPICv2 replicate probes (mean)...")
  base_ids <- sub("_.*$", "", raw_ids)
  int_mat  <- as.matrix(controls@intensity)
  not_na   <- !is.na(int_mat)
  int_mat[!not_na] <- 0
  sums   <- rowsum(int_mat, base_ids)
  counts <- rowsum(matrix(as.numeric(not_na), nrow(int_mat), ncol(int_mat)), base_ids)
  agg <- sums / counts; agg[counts == 0] <- NA
  controls@intensity <- as.data.frame(agg)
  message(sprintf("  %d -> %d probes after collapse.", length(raw_ids), nrow(controls@intensity)))
}

# ---- annotation (same params as cnv_analysis.R matched-controls path) ------
data("exclude_regions"); data("detail_regions")
anno <- conumee2::CNV.create_anno(
  exclude_regions = exclude_regions, detail_regions = detail_regions,
  bin_minprobes = 15, bin_minsize = 50000, array_type = opt$platform
)
common <- intersect(names(anno@probes), rownames(controls@intensity))
if (!length(common)) stop("No probes shared between annotation and controls.")
anno@probes        <- anno@probes[names(anno@probes) %in% common]
controls@intensity <- controls@intensity[common, , drop = FALSE]
message(sprintf("Annotation + controls aligned on %d probes.", length(common)))

# ---- leave-one-out fit + flatness scoring (shared helper) -----------------
plot_fn <- NULL
if (!is.null(opt$out_pdf)) {
  dir.create(dirname(opt$out_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(opt$out_pdf, width = 10, height = 5)
  plot_fn <- function(seg, title) CNV.genomeplot(
    seg, main = title,
    cols = c("#018571", "#80cdc1", "#f5f5f5", "#dfc27d", "#a6611a"))
}

message("\nLeave-one-out flatness scan:")
res <- cnv_flatness_scan(controls, anno,
                         seg_thresh  = opt$seg_thresh,
                         frac_thresh = opt$frac_thresh,
                         plot_fn     = plot_fn)
if (!is.null(opt$out_pdf)) dev.off()

cat("\n================ CNV control flatness ================\n")
print(res, row.names = FALSE)
cat("\nThresholds: |seg.mean| >", opt$seg_thresh,
    "AND altered-fraction >", opt$frac_thresh, "=> NOT-FLAT\n")
bad <- res$sample[res$verdict == "NOT-FLAT"]
if (length(bad)) {
  cat(sprintf("\n%d control(s) are NOT flat and should be excluded:\n  %s\n",
              length(bad), paste(bad, collapse = ", ")))
  cat("Rebuild with:  --exclude ", paste(bad, collapse = ","), "\n", sep = "")
} else {
  cat("\nAll controls look flat at these thresholds. ✔\n")
}
if (!is.null(opt$out_pdf)) cat("\nPer-control genome plots:", opt$out_pdf, "\n")
