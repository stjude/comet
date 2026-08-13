# app/R/array_detect.R
# ---------------------------------------------------------------------------
# Detect Illumina methylation array type from an IDAT file.
#
# We use file size as a coarse heuristic. Real IDATs for each array have
# known, stable sizes:
#
#   450K   ~ 8 MB     ( 622K probes, ~13 bytes per probe)
#   EPIC   ~11 MB     ( 866K probes)
#   EPICv2 ~15 MB     (~935K probes, plus v2 extensions)
#
# A size-based heuristic is fast (no file parsing, no Bioconductor dep),
# good enough for the UI preview, and is never load-bearing at run time
# because the pipeline itself will read the actual manifest before any
# analysis. The user can always override via the dropdown.
# ---------------------------------------------------------------------------

# Size thresholds in bytes, expressed in decimal MB to match what
# macOS Finder / Windows Explorer report (1 MB = 1,000,000 bytes).
# Observed real _Grn.idat sizes from the bundled examples:
#   EPIC v1 (GSM3735546_201465940014_R01C01): ~13.68 MB
#   EPICv2  (209251850130_R03C01):            ~14.37 MB
# These differ by ~0.7 MB; boundary sits at 14.0 MB (midpoint). If
# real-world IDATs land in the dead zone around 14.0 MB, switch to a
# manifest-header parse (illuminaio::readIDAT) instead of size.
ARRAY_SIZE_THRESHOLDS <- list(
  `450K`   = c(min =  5e6, max = 10e6),
  `EPIC`   = c(min = 10e6, max = 14e6),
  `EPICv2` = c(min = 14e6, max = 20e6)
)

#' Guess array type from a Basename stem (expects _Grn.idat to exist).
#'
#' @return list(array_type = "450K"|"EPIC"|"EPICv2"|NA, size = integer|NA,
#'              reason = character)
detect_array_from_idat <- function(basename_stem) {
  grn <- paste0(basename_stem, "_Grn.idat")
  if (!file.exists(grn)) {
    return(list(array_type = NA_character_, size = NA_integer_,
                reason = "IDAT not found"))
  }
  sz <- file.info(grn)$size
  for (arr in names(ARRAY_SIZE_THRESHOLDS)) {
    th <- ARRAY_SIZE_THRESHOLDS[[arr]]
    if (sz >= th["min"] && sz < th["max"]) {
      return(list(array_type = arr, size = sz,
                  reason = sprintf("IDAT size %.1f MB", sz / 1e6)))
    }
  }
  list(array_type = NA_character_, size = sz,
       reason = sprintf("IDAT size %.1f MB outside known ranges",
                        sz / 1e6))
}

#' Infer the array type for a validated samplesheet.
#'
#' Consensus across samples: if every detectable row agrees, return that
#' array type. If they disagree, return `NA_character_` and note the
#' disagreement so the UI can prompt the user to override.
#'
#' @param validated_df  output of validators.R::validate_rows()
#' @return list(array_type, reason, per_sample = data.frame)
detect_array_from_samplesheet <- function(validated_df) {
  if (is.null(validated_df) || nrow(validated_df) == 0L) {
    return(list(array_type = NA_character_,
                reason = "No rows in samplesheet.",
                per_sample = NULL))
  }
  per_sample <- do.call(rbind, lapply(seq_len(nrow(validated_df)), function(i) {
    stem <- validated_df$ResolvedBasename[i]
    if (is.na(stem) || !nzchar(stem)) {
      return(data.frame(Sample_Name = validated_df$Sample_Name[i],
                        array_type = NA_character_, size = NA_integer_,
                        reason = "unresolved basename",
                        stringsAsFactors = FALSE))
    }
    d <- detect_array_from_idat(stem)
    data.frame(Sample_Name = validated_df$Sample_Name[i],
               array_type  = d$array_type %||% NA_character_,
               size        = d$size %||% NA_integer_,
               reason      = d$reason,
               stringsAsFactors = FALSE)
  }))
  detectable <- per_sample$array_type[!is.na(per_sample$array_type)]
  if (!length(detectable)) {
    return(list(array_type = NA_character_,
                reason = "Could not detect from any sample.",
                per_sample = per_sample))
  }
  tab <- table(detectable)
  if (length(tab) == 1L) {
    return(list(array_type = names(tab)[1],
                reason = sprintf("Detected from %d/%d samples.",
                                 length(detectable), nrow(per_sample)),
                per_sample = per_sample))
  }
  # Disagreement across samples — unusual, let user choose.
  list(array_type = NA_character_,
       reason = sprintf("Samples disagree: %s. Please override.",
                        paste(sprintf("%s=%d", names(tab), tab),
                              collapse = ", ")),
       per_sample = per_sample)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
