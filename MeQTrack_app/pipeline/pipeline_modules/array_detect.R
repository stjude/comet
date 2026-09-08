# ---------------------------------------------------------------------------
# Detect Illumina methylation array type per sample.
#
# Shared by the Shiny app (app/R/array_detect.R sources this) and the CLI
# pipeline, so the size thresholds live in exactly one place.
#
# Two detection routes, in priority order:
#   1. An Array_Type / ArrayType / Platform column in the samplesheet. Explicit
#      beats inferred, and it has no dead zone.
#   2. IDAT file size. Fast (no parsing, no Bioconductor dependency) but the
#      EPIC v1 / EPICv2 boundary is narrow -- see the note on thresholds below.
#
# Unlike the app's consensus helper, detect_sample_platforms() deliberately
# returns the full per-sample table and does NOT collapse disagreement to NA:
# a mixed-platform cohort is a supported input, not an error.
#
# Public API:
#   detect_array_from_idat(basename_stem)
#   detect_sample_platforms(sample_sheet, prefer_column = TRUE)
#   split_sheet_by_platform(sample_sheet, ...)
# ---------------------------------------------------------------------------

# Size thresholds in bytes, in decimal MB to match what Finder / Explorer
# report (1 MB = 1,000,000 bytes). Observed real _Grn.idat sizes:
#   450K   ~  8 MB  (622K probes)
#   EPIC   ~ 11-13.7 MB  (866K probes)
#   EPICv2 ~ 14.4-15 MB  (~935K probes plus v2 extensions)
# The EPIC/EPICv2 gap is only ~0.7 MB, so the 14.0 MB boundary is a midpoint
# guess. Samples landing within DEAD_ZONE_MB of it are reported as ambiguous
# rather than silently assigned -- prefer an explicit samplesheet column there.
ARRAY_SIZE_THRESHOLDS <- list(
  `450K`   = c(min =  5e6, max = 10e6),
  `EPIC`   = c(min = 10e6, max = 14e6),
  `EPICv2` = c(min = 14e6, max = 20e6)
)

# Half-width of the ambiguous band around the EPIC/EPICv2 boundary, in bytes.
ARRAY_SIZE_DEAD_ZONE <- 0.35e6

#' Canonicalise a platform label.
#'
#' Accepts the spellings that appear in real samplesheets (450k, HM450,
#' EPICv1, EPIC_v2, MSA, ...) and maps them onto the pipeline's labels.
#'
#' @param x Character vector of raw labels.
#' @return Character vector of "450K" / "EPIC" / "EPICv2", or NA if unknown.
canonical_array_type <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- gsub("[^A-Z0-9]", "", y)                 # EPIC_v2 / EPIC-V2 -> EPICV2
  out <- rep(NA_character_, length(y))
  out[y %in% c("450K", "HM450", "450", "HUMANMETHYLATION450", "ILLUMINA450K")] <- "450K"
  out[y %in% c("EPIC", "EPICV1", "EPIC1", "850K", "HUMANMETHYLATIONEPIC")] <- "EPIC"
  out[y %in% c("EPICV2", "EPIC2", "935K", "HUMANMETHYLATIONEPICV2")] <- "EPICv2"
  out
}

#' Guess array type from a Basename stem (expects _Grn.idat to exist).
#'
#' @param basename_stem IDAT path without the _Grn.idat / _Red.idat suffix.
#' @return list(array_type, size, reason, ambiguous)
detect_array_from_idat <- function(basename_stem) {
  grn <- paste0(basename_stem, "_Grn.idat")
  if (!file.exists(grn)) {
    # Some archives ship gzipped IDATs.
    if (file.exists(paste0(grn, ".gz"))) {
      grn <- paste0(grn, ".gz")
    } else {
      return(list(array_type = NA_character_, size = NA_integer_,
                  reason = "IDAT not found", ambiguous = FALSE))
    }
  }
  sz <- file.info(grn)$size

  for (arr in names(ARRAY_SIZE_THRESHOLDS)) {
    th <- ARRAY_SIZE_THRESHOLDS[[arr]]
    if (sz >= th["min"] && sz < th["max"]) {
      # Flag the narrow EPIC/EPICv2 boundary band as ambiguous.
      amb <- abs(sz - ARRAY_SIZE_THRESHOLDS$EPICv2["min"]) < ARRAY_SIZE_DEAD_ZONE
      return(list(array_type = arr, size = sz,
                  reason = sprintf("IDAT size %.2f MB%s", sz / 1e6,
                                   if (amb) " (near EPIC/EPICv2 boundary)" else ""),
                  ambiguous = unname(amb)))
    }
  }
  list(array_type = NA_character_, size = sz,
       reason = sprintf("IDAT size %.2f MB outside known ranges", sz / 1e6),
       ambiguous = FALSE)
}

#' Determine the platform of every sample in a samplesheet.
#'
#' Prefers an explicit platform column when one is present, falling back to the
#' IDAT size heuristic per sample. Returns one row per sample so that a mixed
#' cohort can be split; it never collapses disagreement into an error.
#'
#' @param sample_sheet  Data frame with a Basename column (and optionally
#'                      Array_Type / ArrayType / Platform, Sample_Name).
#' @param prefer_column Use the samplesheet column when available.
#' @return Data frame: Sample_Name, Basename, platform, source, size, reason,
#'         ambiguous.
detect_sample_platforms <- function(sample_sheet, prefer_column = TRUE) {
  if (is.null(sample_sheet) || nrow(sample_sheet) == 0) {
    stop("detect_sample_platforms: empty sample sheet.")
  }
  base_col <- intersect(c("Basename", "ResolvedBasename"), names(sample_sheet))
  if (!length(base_col)) {
    stop("detect_sample_platforms: sample sheet needs a Basename column.")
  }
  base_col <- base_col[1]

  name_col <- intersect(c("Sample_Name", "Sample_ID", "Sentrix_ID"),
                        names(sample_sheet))
  nm <- if (length(name_col)) as.character(sample_sheet[[name_col[1]]]) else
    basename(as.character(sample_sheet[[base_col]]))

  col <- intersect(c("Array_Type", "ArrayType", "array_type", "Platform"),
                   names(sample_sheet))
  declared <- if (prefer_column && length(col)) {
    canonical_array_type(sample_sheet[[col[1]]])
  } else {
    rep(NA_character_, nrow(sample_sheet))
  }
  if (prefer_column && length(col) && all(is.na(declared))) {
    warning("detect_sample_platforms: column '", col[1],
            "' present but no value recognised; falling back to IDAT size.")
  }

  out <- do.call(rbind, lapply(seq_len(nrow(sample_sheet)), function(i) {
    stem <- as.character(sample_sheet[[base_col]][i])
    if (!is.na(declared[i])) {
      return(data.frame(Sample_Name = nm[i], Basename = stem,
                        platform = declared[i], source = "samplesheet",
                        size = NA_real_,
                        reason = sprintf("declared as '%s'",
                                         sample_sheet[[col[1]]][i]),
                        ambiguous = FALSE, stringsAsFactors = FALSE))
    }
    d <- detect_array_from_idat(stem)
    data.frame(Sample_Name = nm[i], Basename = stem,
               platform = if (is.null(d$array_type)) NA_character_ else d$array_type,
               source = "idat_size",
               size = if (is.null(d$size)) NA_real_ else as.numeric(d$size),
               reason = d$reason, ambiguous = isTRUE(d$ambiguous),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL

  tab <- table(out$platform[!is.na(out$platform)])
  message(sprintf("detect_sample_platforms: %d sample(s) -> %s%s",
                  nrow(out),
                  if (length(tab)) paste(sprintf("%s=%d", names(tab), tab),
                                         collapse = ", ") else "none detected",
                  if (anyNA(out$platform))
                    sprintf(", %d undetermined", sum(is.na(out$platform))) else ""))
  if (any(out$ambiguous)) {
    warning(sprintf("%d sample(s) sized near the EPIC/EPICv2 boundary; ",
                    sum(out$ambiguous)),
            "add an Array_Type column to the samplesheet to be certain.")
  }
  out
}

#' Split a samplesheet into one sub-sheet per platform.
#'
#' @param sample_sheet   Data frame as accepted by detect_sample_platforms().
#' @param prefer_column  Passed through.
#' @param drop_unknown   Drop samples whose platform could not be determined
#'                       (otherwise they are grouped under "unknown").
#' @return Named list of data frames, one per platform, each carrying an
#'         attribute "platform".
split_sheet_by_platform <- function(sample_sheet, prefer_column = TRUE,
                                    drop_unknown = FALSE) {
  det <- detect_sample_platforms(sample_sheet, prefer_column = prefer_column)
  pl  <- det$platform
  if (anyNA(pl)) {
    if (drop_unknown) {
      message(sprintf("split_sheet_by_platform: dropping %d undetermined sample(s).",
                      sum(is.na(pl))))
      keep <- !is.na(pl)
      sample_sheet <- sample_sheet[keep, , drop = FALSE]
      pl <- pl[keep]
    } else {
      pl[is.na(pl)] <- "unknown"
    }
  }
  if (!nrow(sample_sheet)) stop("split_sheet_by_platform: no samples left.")

  out <- lapply(split(seq_along(pl), pl), function(idx) {
    s <- sample_sheet[idx, , drop = FALSE]
    rownames(s) <- NULL
    s
  })
  for (p in names(out)) attr(out[[p]], "platform") <- p
  message(sprintf("split_sheet_by_platform: %d platform(s): %s",
                  length(out),
                  paste(sprintf("%s (n=%d)", names(out),
                                vapply(out, nrow, integer(1))), collapse = ", ")))
  out
}
