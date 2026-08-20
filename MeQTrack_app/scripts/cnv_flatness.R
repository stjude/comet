# scripts/cnv_flatness.R
# ---------------------------------------------------------------------------
# Shared leave-one-out CNV flatness scan.
#
# Used by:
#   - build_cnv_controls.R          (gate: drop non-flat controls at build)
#   - check_cnv_controls_flatness.R (report: diagnose a control set)
# One implementation, so the build drops EXACTLY what the check would flag.
#
# Why this exists: a CNV reference must be built from samples with FLAT
# genomes. Any real copy-number event shared by the controls is absorbed into
# the baseline and then appears, inverted, in every query sample. The
# correlation-with-others gate in build_cnv_controls.R cannot see this — a
# shared CNV correlates just fine. Only fitting each control against the
# others and segmenting reveals it.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(conumee2))

# Score every control by leave-one-out: fit control i against the other
# controls, segment, and measure how far its segments stray from log2 0.
#
#   controls    CNV.data; @intensity (probes x samples), already replicate-
#               collapsed and aligned to `anno`'s probes.
#   anno        CNV.anno (CNV.create_anno) aligned to the same probes.
#   seg_thresh  |seg.mean| above this = an altered segment.
#   frac_thresh altered fraction of the segmented genome above this = NOT-FLAT.
#   plot_fn     optional function(seg_obj, title) for a per-control genomeplot
#               (device must already be open).
#
# Returns a data.frame ordered worst-first:
#   sample, n_seg, n_seg_altered, max_abs_segmean, sd_segmean,
#   frac_genome_alt, verdict ("flat" | "NOT-FLAT")
cnv_flatness_scan <- function(controls, anno, seg_thresh = 0.20,
                              frac_thresh = 0.02, plot_fn = NULL) {
  samples <- colnames(controls@intensity)
  if (length(samples) < 3)
    stop("Need >= 3 controls for a leave-one-out flatness scan; have ",
         length(samples), ".")

  score_one <- function(i) {
    q <- controls; q@intensity <- controls@intensity[, i,  drop = FALSE]
    r <- controls; r@intensity <- controls@intensity[, -i, drop = FALSE]
    seg <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(q, r, anno))))
    # conumee2 stores one DNAcopy::segments.summary() data frame per sample in
    # @seg$summary; a single-query fit has exactly one element. Columns:
    # ID, chrom, loc.start, loc.end, num.mark, seg.mean, seg.sd, seg.median...
    sdf <- seg@seg$summary[[1]]
    if (is.null(sdf) || !nrow(sdf))
      return(data.frame(sample = samples[i], n_seg = 0L, n_seg_altered = 0L,
                        max_abs_segmean = NA_real_, sd_segmean = NA_real_,
                        frac_genome_alt = NA_real_, stringsAsFactors = FALSE))
    len     <- pmax(sdf$loc.end - sdf$loc.start, 0)
    altered <- abs(sdf$seg.mean) > seg_thresh
    covered <- sum(len)
    if (!is.null(plot_fn)) plot_fn(seg, sprintf("%s (LOO)", samples[i]))
    data.frame(
      sample          = samples[i],
      n_seg           = nrow(sdf),
      n_seg_altered   = sum(altered),
      max_abs_segmean = round(max(abs(sdf$seg.mean)), 3),
      sd_segmean      = round(stats::sd(sdf$seg.mean), 3),
      frac_genome_alt = round(if (covered > 0) sum(len[altered]) / covered else 0, 4),
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, lapply(seq_along(samples), function(i) {
    message(sprintf("  [%d/%d] %s", i, length(samples), samples[i]))
    tryCatch(score_one(i), error = function(e) {
      message("    flatness scan FAILED: ", conditionMessage(e)); NULL
    })
  }))
  if (is.null(res)) stop("Flatness scan produced no results (all fits failed).")

  res$verdict <- ifelse(
    !is.na(res$max_abs_segmean) &
      res$max_abs_segmean > seg_thresh & res$frac_genome_alt > frac_thresh,
    "NOT-FLAT", "flat")
  res[order(-res$frac_genome_alt, -res$max_abs_segmean), ]
}
