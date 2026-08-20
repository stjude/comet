# app/R/array_detect.R
# ---------------------------------------------------------------------------
# Array-type detection for the Shiny app.
#
# ARRAY_SIZE_THRESHOLDS, canonical_array_type() and detect_array_from_idat()
# now live in pipeline/pipeline_modules/array_detect.R so that the app and the
# CLI pipeline share one set of size thresholds. This file keeps only the
# app-specific consensus wrapper below.
#
# The size heuristic is fast (no file parsing, no Bioconductor dep), good
# enough for the UI preview, and never load-bearing at run time because the
# pipeline reads the actual manifest before any analysis. The user can always
# override via the dropdown.
# ---------------------------------------------------------------------------

# Resolve the shared module relative to the app, then fall back to the
# locations the pipeline's own find_anno_file() checks, so this works whether
# the app is launched from the repo root or from app/.
local({
  candidates <- c(
    file.path("pipeline", "pipeline_modules", "array_detect.R"),
    file.path("..", "pipeline", "pipeline_modules", "array_detect.R"),
    file.path(getwd(), "pipeline", "pipeline_modules", "array_detect.R")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("array_detect.R: cannot locate pipeline_modules/array_detect.R; ",
         "expected one of: ", paste(candidates, collapse = ", "))
  }
  source(hit[1], local = FALSE)
})

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
