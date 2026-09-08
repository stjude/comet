# app/R/validators.R
# ---------------------------------------------------------------------------
# Pure, side-effect-free samplesheet validation.
#
# These functions know nothing about Shiny. They take file paths (or a
# data.frame) and return a validated structure the Shiny layer renders.
#
# Keeping the logic here, isolated from reactives, means we can (a) unit
# test it without spinning up Shiny, and (b) reuse it in the CLI pipeline
# bridge if we ever want pre-flight checks before a run.
# ---------------------------------------------------------------------------

# Required columns per pipeline/README.md and pipeline/pipeline_modules/
# preprocess.R:34.
REQUIRED_COLUMNS <- c("Sentrix_ID", "Sample_Name", "Basename")

# Row-level status codes surfaced in the UI validation table.
# Cross-platform "free bytes on the filesystem holding `path`". Returns
# NA_real_ when the platform tool isn't available or output can't be
# parsed. Used by the run controller for a pre-flight disk-space check
# before kicking off a pipeline run.
free_bytes <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NA_real_)
  if (!file.exists(path)) {
    parent <- dirname(path)
    if (!file.exists(parent)) return(NA_real_)
    path <- parent
  }
  if (.Platform$OS.type == "windows") {
    drive <- toupper(substr(normalizePath(path, mustWork = FALSE), 1, 2))
    out <- tryCatch(
      system2("fsutil", c("volume", "diskfree", drive),
              stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    # fsutil prints a line with a "free bytes" number — pick the first
    # large integer in the line containing "Avail" or "free".
    line <- grep("free|avail", out, ignore.case = TRUE, value = TRUE)
    if (!length(line)) return(NA_real_)
    m <- regmatches(line[1], regexpr("\\d{6,}", line[1]))
    if (!length(m)) return(NA_real_)
    suppressWarnings(as.numeric(m[1]))
  } else {
    out <- tryCatch(
      system2("df", c("-k", shQuote(path)),
              stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    if (length(out) < 2) return(NA_real_)
    fields <- strsplit(trimws(out[length(out)]), "\\s+")[[1]]
    if (length(fields) < 4L) return(NA_real_)
    avail_kb <- suppressWarnings(as.numeric(fields[4]))
    if (is.na(avail_kb)) return(NA_real_)
    avail_kb * 1024
  }
}

STATUS_OK               <- "OK"
STATUS_MISSING_RED      <- "Missing _Red.idat"
STATUS_MISSING_GRN      <- "Missing _Grn.idat"
STATUS_MISSING_BOTH     <- "Missing IDAT pair"
STATUS_DUPLICATE_ID     <- "Duplicate Sentrix_ID"
STATUS_MALFORMED        <- "Malformed row"
STATUS_CORRUPT_IDAT     <- "Corrupt IDAT"

# Verify a file looks like an IDAT by reading the first 4 bytes and
# checking for the "IDAT" magic header. Cheap (<1 ms per file) and
# catches truncated downloads, files with the wrong extension, or
# placeholders that minfi would otherwise crash on mid-preprocess.
idat_magic_ok <- function(path) {
  if (!file.exists(path)) return(FALSE)
  con <- tryCatch(file(path, "rb"), error = function(e) NULL)
  if (is.null(con)) return(FALSE)
  on.exit(close(con), add = TRUE)
  bytes <- tryCatch(readBin(con, "raw", n = 4L),
                    error = function(e) raw(0))
  if (length(bytes) < 4L) return(FALSE)
  identical(rawToChar(bytes), "IDAT")
}

#' Read a samplesheet CSV and return a data.frame + top-level error status.
#'
#' Top-level errors (missing file, unparseable CSV, missing required columns)
#' short-circuit row-level validation. The caller's UI should render the
#' top-level error prominently and hide the per-row table.
#'
#' @param path   Absolute path to a .csv samplesheet.
#' @return list(ok = logical, error = character|NULL, df = data.frame|NULL,
#'              samplesheet_dir = character|NULL)
read_samplesheet <- function(path) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, error = paste("Samplesheet not found:", path),
                df = NULL, samplesheet_dir = NULL))
  }
  # fileEncoding = "UTF-8-BOM" strips an optional UTF-8 BOM (common on
  # Excel/Numbers-exported CSVs); harmless for plain ASCII files.
  # Don't swallow warnings into NULL — read.csv routinely warns on
  # incomplete-final-line and similar non-fatal quirks; treat the parse
  # as successful as long as it returns a data frame with rows.
  df <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE,
                    fileEncoding = "UTF-8-BOM"),
    error = function(e) NULL
  )
  if (is.null(df)) {
    df <- tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE),
      error = function(e) NULL
    )
  }
  if (is.null(df) || nrow(df) == 0L) {
    return(list(ok = FALSE,
                error = "Samplesheet is empty or could not be parsed as CSV.",
                df = NULL, samplesheet_dir = NULL))
  }
  # Defensive BOM strip on first column name in case the encoding hint
  # didn't catch it (rare older R versions).
  colnames(df)[1] <- sub("^\\xef\\xbb\\xbf|^﻿", "", colnames(df)[1])
  missing_cols <- setdiff(REQUIRED_COLUMNS, colnames(df))
  if (length(missing_cols)) {
    return(list(
      ok = FALSE,
      error = sprintf(
        "Samplesheet is missing required column(s): %s. Required columns are: %s.",
        paste(missing_cols, collapse = ", "),
        paste(REQUIRED_COLUMNS, collapse = ", ")
      ),
      df = NULL, samplesheet_dir = NULL
    ))
  }
  list(ok = TRUE, error = NULL, df = df,
       samplesheet_dir = normalizePath(dirname(path), winslash = "/",
                                       mustWork = TRUE))
}

#' Resolve a Basename to the absolute stem path where IDATs should live.
#'
#' If Basename is absolute, return it as-is. If relative, try resolving
#' relative to the samplesheet's own directory first (most intuitive for
#' user-authored samplesheets) and fall back to the pipeline folder
#' (matches the CLI pipeline's setwd(script_dir) behaviour and supports
#' the bundled example samplesheet unchanged).
#'
#' Always returns a normalized path. Existence is NOT checked here —
#' validate_rows() does the file-existence check per row.
#'
#' @return list(resolved = character, tried = character)
resolve_basename <- function(basename, samplesheet_dir, pipeline_dir) {
  basename <- as.character(basename)
  if (!nzchar(basename)) {
    return(list(resolved = NA_character_, tried = character(0)))
  }
  # Absolute path — try as-is first; if the IDATs aren't there (common
  # when a samplesheet was authored on a different machine), fall back
  # to looking for the same filename stem in the samplesheet dir or the
  # pipeline dir.
  is_abs <- startsWith(basename, "/") || grepl("^[A-Za-z]:[/\\\\]", basename)
  if (is_abs) {
    abs_red <- paste0(basename, "_Red.idat")
    abs_grn <- paste0(basename, "_Grn.idat")
    if (file.exists(abs_red) && file.exists(abs_grn)) {
      return(list(
        resolved = normalizePath(basename, winslash = "/", mustWork = FALSE),
        tried = basename
      ))
    }
    # Fall through to the same dir-fallback used for relative paths,
    # using just the stem (final path component).
    stem <- basename(basename)
    fallback <- c(
      file.path(samplesheet_dir, stem),
      file.path(pipeline_dir,    stem)
    )
    for (cand in fallback) {
      if (file.exists(paste0(cand, "_Red.idat")) &&
          file.exists(paste0(cand, "_Grn.idat"))) {
        return(list(
          resolved = normalizePath(cand, winslash = "/", mustWork = FALSE),
          tried = c(basename, fallback)
        ))
      }
    }
    return(list(
      resolved = normalizePath(basename, winslash = "/", mustWork = FALSE),
      tried = c(basename, fallback)
    ))
  }
  # Relative — try samplesheet dir, then pipeline dir.
  candidates <- c(
    file.path(samplesheet_dir, basename),
    file.path(pipeline_dir,    basename)
  )
  for (cand in candidates) {
    red <- paste0(cand, "_Red.idat")
    grn <- paste0(cand, "_Grn.idat")
    if (file.exists(red) && file.exists(grn)) {
      return(list(
        resolved = normalizePath(cand, winslash = "/", mustWork = FALSE),
        tried = candidates
      ))
    }
  }
  list(
    resolved = normalizePath(candidates[1], winslash = "/", mustWork = FALSE),
    tried = candidates
  )
}

#' Per-row validation. Takes the data.frame from read_samplesheet() and
#' returns the same frame with three added columns: ResolvedBasename,
#' Status, StatusDetail.
#'
#' @param df             parsed samplesheet
#' @param samplesheet_dir  directory containing the samplesheet
#' @param pipeline_dir     pipeline/ folder for fallback resolution
#' @return data.frame with added validation columns
validate_rows <- function(df, samplesheet_dir, pipeline_dir) {
  n <- nrow(df)
  resolved    <- character(n)
  status      <- character(n)
  status_det  <- character(n)

  dup_ids <- df$Sentrix_ID[duplicated(df$Sentrix_ID)]

  for (i in seq_len(n)) {
    bn     <- df$Basename[i]
    sid    <- df$Sentrix_ID[i]
    sname  <- df$Sample_Name[i]

    # Malformed check: any required field empty or NA.
    if (is.na(bn) || !nzchar(bn) ||
        is.na(sid) || !nzchar(sid) ||
        is.na(sname) || !nzchar(sname)) {
      resolved[i]   <- NA_character_
      status[i]     <- STATUS_MALFORMED
      status_det[i] <- "Missing Sentrix_ID, Sample_Name, or Basename."
      next
    }

    r <- resolve_basename(bn, samplesheet_dir, pipeline_dir)
    resolved[i] <- r$resolved
    red <- paste0(r$resolved, "_Red.idat")
    grn <- paste0(r$resolved, "_Grn.idat")
    has_red <- file.exists(red)
    has_grn <- file.exists(grn)

    # Duplicate Sentrix_ID takes priority so the user sees which row is
    # problematic even if the IDATs also happen to be missing.
    if (sid %in% dup_ids) {
      status[i]     <- STATUS_DUPLICATE_ID
      status_det[i] <- sprintf("Sentrix_ID '%s' appears on more than one row.", sid)
    } else if (!has_red && !has_grn) {
      status[i]     <- STATUS_MISSING_BOTH
      status_det[i] <- sprintf("Neither %s_Red.idat nor %s_Grn.idat found.",
                               basename(r$resolved), basename(r$resolved))
    } else if (!has_red) {
      status[i]     <- STATUS_MISSING_RED
      status_det[i] <- sprintf("Expected %s, not found.", red)
    } else if (!has_grn) {
      status[i]     <- STATUS_MISSING_GRN
      status_det[i] <- sprintf("Expected %s, not found.", grn)
    } else {
      # Both files present — verify they look like real IDATs. minfi will
      # crash mid-preprocess on a truncated download or a placeholder file
      # with the .idat extension, with an error far less obvious than this.
      bad_red <- !idat_magic_ok(red)
      bad_grn <- !idat_magic_ok(grn)
      if (bad_red && bad_grn) {
        status[i]     <- STATUS_CORRUPT_IDAT
        status_det[i] <- "Both IDATs fail the magic-header check (corrupt or truncated)."
      } else if (bad_red) {
        status[i]     <- STATUS_CORRUPT_IDAT
        status_det[i] <- sprintf("%s fails the IDAT magic-header check (corrupt or truncated).",
                                 basename(red))
      } else if (bad_grn) {
        status[i]     <- STATUS_CORRUPT_IDAT
        status_det[i] <- sprintf("%s fails the IDAT magic-header check (corrupt or truncated).",
                                 basename(grn))
      } else {
        status[i]     <- STATUS_OK
        status_det[i] <- ""
      }
    }
  }

  df$ResolvedBasename <- resolved
  df$Status           <- status
  df$StatusDetail     <- status_det
  df
}

#' Summary of a validation result: overall ok/not and row counts by status.
validation_summary <- function(validated_df) {
  if (is.null(validated_df) || nrow(validated_df) == 0L) {
    return(list(ok = FALSE, n_total = 0L, n_ok = 0L, n_bad = 0L,
                counts = integer(0)))
  }
  counts <- table(validated_df$Status)
  n_ok  <- as.integer(counts[STATUS_OK] %||% 0L)
  n_bad <- nrow(validated_df) - n_ok
  list(ok = (n_bad == 0L),
       n_total = nrow(validated_df),
       n_ok = n_ok,
       n_bad = n_bad,
       counts = counts)
}

#' Optional metadata columns: everything in the samplesheet that isn't
#' required and isn't one of our synthesized validation columns.
optional_metadata_columns <- function(validated_df) {
  if (is.null(validated_df)) return(character(0))
  synthesized <- c("ResolvedBasename", "Status", "StatusDetail")
  setdiff(colnames(validated_df), c(REQUIRED_COLUMNS, synthesized))
}

# local null-coalesce to avoid sourcing pipeline_bridge.R just for this
`%||%` <- function(a, b) if (is.null(a)) b else a
