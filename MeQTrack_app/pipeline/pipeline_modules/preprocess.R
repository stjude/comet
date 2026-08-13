# Preprocessing module for methylation array analysis pipeline

#' Read sample sheet
#'
#' @param sample_sheet_path Path to sample sheet CSV
#' @return Data frame containing sample information
read_sample_sheet <- function(sample_sheet_path) {
  if (!file.exists(sample_sheet_path)) {
    stop("Sample sheet not found: ", sample_sheet_path)
  }
  
  # Function to read data files
  read_data <- function(file) {
      message("Reading data from: ", file)
      data <- fread(file, header = TRUE)
      data <- as.data.frame(data)
      
      # Clean column names if they start with X
      if (grepl("X", colnames(data)[1])) {
        colnames(data) <- gsub("X", "", colnames(data))
      }
      
      # Set row names and remove the first column if it's an index
      if (colnames(data[1]) != "Sentrix_ID") {
        rownames(data) <- data[, 1]
        data <- data[, -1]
      }
      
      return(data)
    }
  sample_sheet <- read_data(sample_sheet_path)
  
  # Check required columns
  required_cols <- c("Sentrix_ID", "Sample_Name", "Basename")
  missing_cols <- required_cols[!required_cols %in% colnames(sample_sheet)]
  
  if (length(missing_cols) > 0) {
    stop("Sample sheet is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  # Check for and remove duplicate Basename entries (would crash read.metharray)
  dup_mask <- duplicated(sample_sheet$Basename)
  if (any(dup_mask)) {
    dup_basenames <- unique(sample_sheet$Basename[dup_mask])
    warning(
      length(dup_basenames), " duplicate Basename(s) found and removed (keeping first occurrence):\n",
      paste0("  ", dup_basenames, collapse = "\n")
    )
    sample_sheet <- sample_sheet[!dup_mask, ]
    message("Sample sheet after deduplication: ", nrow(sample_sheet), " samples remain.")
  }

  return(sample_sheet)
}

#' Determine array type from data
#'
#' @param rgset Raw data RGChannelSet
#' @return String indicating array type (450k, EPIC, EPICv2)
determine_array_type <- function(rgset) {
  # Classify by RGChannelSet row count. These are raw red/green intensity
  # rows (not CpG loci), so they're much larger than probe counts:
  #   450K:    622,399 rows
  #   EPIC:  1,051,943 rows
  #   EPICv2:1,105,209 rows
  # The thresholds below bracket these three known values with ~50k headroom.
  n_probes <- nrow(rgset)

  if (n_probes < 700000) {
    return("450k")
  } else if (n_probes < 1080000) {
    return("EPIC")
  } else {
    return("EPICv2")
  }
}

#' Load appropriate annotation package based on array type
#'
#' @param array_type Array type (450k, EPIC, EPICv2)
#' @return Annotation object
load_array_annotation <- function(array_type) {
  # Handle auto-detection
  if (array_type == "auto") {
    warning("Array type 'auto' is not valid for loading annotations directly.
            Please determine the specific array type first or specify it explicitly.")
    array_type <- "EPIC"  # Default to EPIC as fallback
  }
  # Canonicalize 450K -> 450k (minfi manifest convention). Case mismatch
  # between the UI ("450K") and this pipeline ("450k") was a repeat
  # source of "Unsupported array type" errors.
  if (toupper(array_type) == "450K") array_type <- "450k"

  # Check for valid array types
  if (!array_type %in% c("450k", "EPIC", "EPICv2")) {
    stop("Unsupported array type: ", array_type, 
         ". Must be one of: '450k', 'EPIC', 'EPICv2'")
  }

  # Create a safe getAnnotation function
  safe_get_annotation <- function(pkg_name) {
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      tryCatch({
        # Get the annotation object from the package
        pkg_env <- asNamespace(pkg_name)
        anno_obj <- data(pkg_name)
        return(anno_obj)
      }, error = function(e) {
        # If there's an error, provide details and return NULL
        warning("Error getting annotation from package '", pkg_name, "': ", e$message)
        return(NULL)
      })
    } else {
      return(NULL)
    }
  }
  
  # Try to load the appropriate annotation package
  anno_obj <- NULL
  
  if (array_type == "450k") {
    message("Loading 450k annotation...")
    pkg_name <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    anno_obj <- safe_get_annotation(pkg_name)
    
    if (is.null(anno_obj)) {
      # Try to fall back to EPIC if 450k package is not available
      warning("Could not load 450k annotation. Trying EPIC annotation as fallback...")
      array_type <- "EPIC"
    } else {
      return(anno_obj)
    }
  }
  
  if (array_type == "EPIC") {
    message("Loading EPIC annotation...")
    # Try the newer version first
    pkg_name <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
    anno_obj <- safe_get_annotation(pkg_name)
    
    if (is.null(anno_obj)) {
      # Try older version
      pkg_name <- "IlluminaHumanMethylationEPICanno.ilm10b3.hg19"
      anno_obj <- safe_get_annotation(pkg_name)
      
      if (!is.null(anno_obj)) {
        warning("Using older EPIC annotation (ilm10b3) as ilm10b4 package is not installed")
        return(anno_obj)
      }
    } else {
      return(anno_obj)
    }
  }
  
  if (array_type == "EPICv2" || is.null(anno_obj)) {
    if (array_type == "EPICv2") {
      message("Loading EPICv2 annotation...")
      message("Loading sesame package ...")
      
    } else {
      message("Trying EPICv2 annotation as last resort...")
    }
    
    pkg_name <- "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
    anno_obj <- safe_get_annotation(pkg_name)
    
    if (!is.null(anno_obj)) {
      return(anno_obj)
    }
  }
}
#' Read methylation array data
#'
#' @param sample_sheet Sample sheet data frame
#' @param array_type Array type
#' @return RGChannelSet object
read_methylation_data <- function(sample_sheet, array_type = "auto") {
  # Read IDAT files
  message(" array data...")
  rgset <- read.metharray.exp(base = NULL, 
                             targets = sample_sheet, 
                             recursive = TRUE,
                             force = TRUE)
  
  # Set annotation
  if (array_type == "auto") {
    array_type <- determine_array_type(rgset)
    message("Auto-detected array type: ", array_type)
  }
  # Canonicalize 450K -> 450k so callers using either capitalization work.
  if (toupper(array_type) == "450K") array_type <- "450k"

  if (array_type == "450k") {
    rgset@annotation <- c(array = "IlluminaHumanMethylation450k", 
                         annotation = "ilmn12.hg19")
  } else if (array_type == "EPIC") {
    rgset@annotation <- c(array = "IlluminaHumanMethylationEPIC", 
                         annotation = "ilm10b4.hg19")
  } else if (array_type == "EPICv2") {
    rgset@annotation <- c(array = "IlluminaHumanMethylationEPICv2", 
                         annotation = "20a1.hg38")
  } else {
    stop("Unsupported array type: ", array_type)
  }
  
  return(rgset)
}

#' Compute per-sample GCT bisulfite-conversion control score
#'
#' GCT (Zhou et al. 2017) quantifies residual *incomplete* bisulfite
#' conversion from the Infinium-I C/T-extension probes. A score near 1.0
#' indicates complete conversion; higher values indicate more incomplete
#' conversion. This is the specific QC metric sesame provides that minfi does
#' not, and the reason preprocessing reads IDATs through sesame.
#'
#' sesame::bisConversionControl() auto-fetches the required extension probes
#' only for EPIC and HM450 (450k). EPICv2/MSA need extR/extA supplied manually
#' (Phase 2), so those array types return NA with an explanatory note rather
#' than erroring. The metric is informational only — it never gates Pass_QC.
#'
#' @param basenames IDAT basename prefixes (sample_sheet$Basename)
#' @param sample_ids Per-sample identifiers, aligned to basenames order
#' @param array_type Array type ("450k", "EPIC", "EPICv2", ...)
#' @param bpparam BiocParallel backend (reused from the caller)
#' @return data.frame with Sample_ID, GCT_Score, Array_Type, Note (one row/sample)
compute_gct_scores <- function(basenames, sample_ids, array_type, bpparam) {
  at <- toupper(array_type)

  # bisConversionControl auto-fetches its C/T-extension probes for EPIC/HM450
  # only. For EPICv2 (no sesameData EPICv2.probeInfo) we supply extR/extA from a
  # vendored probe list derived from the Zhou Lab EPICv2 manifest (type-I probes
  # split by extension base: nextBase "R" -> ext-C, "A" -> ext-T). This was
  # validated to reproduce sesame's native GCT exactly on EPIC. Other platforms
  # (MSA, HM27) stay NA until their ext probes are sourced.
  gct_func <- sesame::bisConversionControl
  if (at == "EPICV2") {
    ext <- load_epicv2_ext_probes()
    if (is.null(ext)) {
      message("EPICv2 GCT ext-probe list unavailable — emitting NA.")
      return(data.frame(
        Sample_ID  = sample_ids,
        GCT_Score  = NA_real_,
        Array_Type = array_type,
        Note       = "GCT skipped: EPICv2 ext-probe list not found under Anno/EPICv2/",
        stringsAsFactors = FALSE
      ))
    }
    gct_func <- function(sdf) {
      sesame::bisConversionControl(sdf, extR = ext$extC, extA = ext$extT)
    }
  } else if (!at %in% c("EPIC", "450K")) {
    message("GCT bisulfite-conversion control not yet supported for array type '",
            array_type, "' — emitting NA.")
    return(data.frame(
      Sample_ID  = sample_ids,
      GCT_Score  = NA_real_,
      Array_Type = array_type,
      Note       = "GCT not yet supported for this array type",
      stringsAsFactors = FALSE
    ))
  }

  message("Computing GCT bisulfite-conversion control scores...")
  # prep = "C" (inferInfiniumIChannel) is the minimal prep bisConversionControl
  # needs: it reads InfIR(sdf), which requires Infinium-I channel assignment.
  # Heavier prep (noob/dyebias) would distort the raw extension-probe signal.
  scores <- tryCatch(
    sesame::openSesame(basenames, prep = "C",
                       func = gct_func,
                       BPPARAM = bpparam),
    error = function(e) {
      warning("GCT computation failed (", conditionMessage(e),
              "); recording NA for all samples.")
      NULL
    }
  )

  if (is.null(scores)) {
    return(data.frame(
      Sample_ID  = sample_ids,
      GCT_Score  = NA_real_,
      Array_Type = array_type,
      Note       = "GCT computation failed; see pipeline log",
      stringsAsFactors = FALSE
    ))
  }

  # openSesame returns a named numeric vector (names = IDAT prefixes/basenames),
  # in the same order as `basenames`. Align by position to sample_ids.
  scores <- as.numeric(scores)
  if (length(scores) != length(sample_ids)) {
    warning("GCT score count (", length(scores), ") != sample count (",
            length(sample_ids), "); results may be misaligned.")
    length(scores) <- length(sample_ids)  # pad/truncate with NA to keep shape
  }

  data.frame(
    Sample_ID  = sample_ids,
    GCT_Score  = round(scores, 4),
    Array_Type = array_type,
    Note       = ifelse(is.na(scores), "GCT not computed for this sample", ""),
    stringsAsFactors = FALSE
  )
}

#' Dye-bias Red/Green QQ plots (sesame) — one PNG per sample.
#'
#' sesame::sesameQC_plotRedGrnQQ draws a Red-vs-Green quantile-quantile plot
#' that exposes dye bias: a strong departure from the diagonal means the two
#' colour channels are imbalanced. Computed on the channel-inferred SigDF
#' (prep = "C") — i.e. BEFORE the pipeline's dye-bias correction — so the plot
#' shows the bias that downstream noob/dyeBiasNL then corrects.
#'
#' Writes one rasterized PNG per sample to figures/qc/dye_bias/<Sample_ID>.png
#' (~40 KB each) rather than a single vector PDF: the QQ plots every address, so
#' a vector page is ~1 MB and a 200-sample PDF would be hundreds of MB. The app
#' shows them one at a time behind a sample selector. Each sample is wrapped so
#' one unreadable IDAT can't abort the rest.
#'
#' @param basenames IDAT basenames (sample_sheet$Basename)
#' @param sample_ids Per-sample IDs (drive the file names + titles)
#' @param figures_qc_dir The run's figures/qc directory (dirs$figures_qc); PNGs
#'   land in a dye_bias/ subfolder of it
plot_dye_bias_qq <- function(basenames, sample_ids, figures_qc_dir) {
  fig_dir <- file.path(figures_qc_dir, "dye_bias")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  message("Plotting dye-bias R/G QQ (", length(basenames), " samples)...")
  safe <- gsub("[^A-Za-z0-9._-]", "_", sample_ids)
  n_ok <- 0L
  for (i in seq_along(basenames)) {
    png_path <- file.path(fig_dir, paste0(safe[i], ".png"))
    ok <- tryCatch({
      sdf <- sesame::inferInfiniumIChannel(sesame::readIDATpair(basenames[i]))
      grDevices::png(png_path, width = 700, height = 700, res = 110)
      tryCatch(
        sesame::sesameQC_plotRedGrnQQ(sdf, main = paste0("R-G QQ: ", sample_ids[i])),
        finally = grDevices::dev.off()
      )
      TRUE
    }, error = function(e) {
      warning("Dye-bias QQ failed for ", sample_ids[i], ": ", conditionMessage(e))
      if (file.exists(png_path)) unlink(png_path)
      FALSE
    })
    if (isTRUE(ok)) n_ok <- n_ok + 1L
  }
  message("Dye-bias QQ written for ", n_ok, "/", length(basenames),
          " samples: figures/qc/dye_bias/")
  invisible(fig_dir)
}

#' Preprocess methylation data
#'
#' @param sample_sheet Sample sheet data frame
#' @param array_type Array type
#' @param normalization Normalization method
#' @param threads Number of CPU threads to use
#' @param output_dir Output directory
#' @return List containing preprocessed data
preprocess_methylation <- function(sample_sheet, array_type = "auto",
                                   normalization = "swan", 
                                   threads = 6,
                                   output_dir = ".") {
  # Build a BiocParallel backend appropriate for the current environment.
  #
  # MulticoreParam (fork-based) causes two known failures on LSF/bsub:
  #   1. "wrong args for environment subassignment" – forked workers inherit
  #      large S4 objects and crash during environment subassignment.
  #   2. "reached CPU time limit" – LSF tracks cumulative CPU time across all
  #      forked child processes and kills the job when the queue limit is hit.
  #
  # Solution: force SerialParam when running inside an LSF job (detected via
  # the LSB_JOBID environment variable), otherwise use MulticoreParam.
  in_lsf <- nchar(Sys.getenv("LSB_JOBID")) > 0
  bpparam <- tryCatch({
    if (in_lsf) {
      message("LSF job detected (LSB_JOBID=", Sys.getenv("LSB_JOBID"),
              "): using SerialParam to avoid fork-based CPU time limit issues.")
      BiocParallel::SerialParam()
    } else if (threads > 1 && .Platform$OS.type == "unix") {
      message("Parallel backend: MulticoreParam (", threads, " workers)")
      BiocParallel::MulticoreParam(threads)
    } else {
      message("Parallel backend: SerialParam (threads = ", threads, ")")
      BiocParallel::SerialParam()
    }
  }, error = function(e) {
    message("Backend creation failed (", conditionMessage(e), "); falling back to SerialParam.")
    BiocParallel::SerialParam()
  })
  BiocParallel::register(bpparam)

  # Read data
  rgset <- read_methylation_data(sample_sheet, array_type)
  print(rgset)
  message("Check rgSet: ", class(rgset)[1], " with ", ncol(rgset), " samples and ", nrow(rgset), " probes")
  save(rgset, file = file.path(output_dir, "rgset.RData"))

  # Get sample information
  sample_info <- data.frame(
    Sample_ID = sample_sheet$Sentrix_ID,
    Sample_Name = sample_sheet$Sample_Name
  )

  # Check normalization method and use default if invalid
  if (is.null(normalization) || !normalization %in% c("raw", "illumina", "functional", "quantile", "swan", "sesame")) {
    normalization <- "swan"
    warning("Invalid normalization method specified. Using 'SWAN' normalization instead.")
  }

  ####If EPICv2, then use sesame to prepare idats and extract beta values
  ## prep = "QCDB" applies these technical corrections (left-to-right):
  ## Q = qualityMask, C = infer Infinium-I channel, D = dyeBiasNL, B = Noob.
  if (array_type == "EPICv2"){
    beta_v2 <- sesame::openSesame(sample_sheet$Basename, prep = "QCDB", platform= "EPICv2", func = getBetas, BPPARAM=bpparam,
                                  collapseToPfx = TRUE, collapseMethod = "mean")
    # Reorder beta table to be same order as sample_table
  }else {
    beta_v2 <- sesame::openSesame(sample_sheet$Basename, prep = "QCDB", func = getBetas, BPPARAM=bpparam)
  }
  beta <- beta_v2[, match(sample_sheet$Sentrix_ID, colnames(beta_v2))]
  message("Number of rows in beta: ", nrow(beta))

  ## openSesame(..., func = pOOBAH) doesn't accept collapseToPfx /
  ## collapseMethod in some sesame versions (they get forwarded to
  ## pOOBAH and error out as unused args). Instead we let pOOBAH run
  ## at the per-replicate granularity and then collapse manually to
  ## match beta's per-CpG space for EPICv2.
  detection_p <- sesame::openSesame(sample_sheet$Basename, func = pOOBAH,
                                    return.pval = TRUE, BPPARAM = bpparam)

  if (array_type == "EPICv2" &&
      any(grepl("_(BC|TC)\\d+$", utils::head(rownames(detection_p), 200)))) {
    message("Collapsing EPICv2 detection_p replicates -> one row per CpG (mean)...")
    dp_base <- sub("_.*$", "", rownames(detection_p))
    if (anyDuplicated(dp_base)) {
      # Per-CpG mean across replicates. rowsum() sums by group; divide
      # by per-group valid-count (non-NA) to get a correct mean.
      dp_mat <- as.matrix(detection_p)
      not_na <- !is.na(dp_mat)
      dp_mat[!not_na] <- 0
      sums   <- rowsum(dp_mat, dp_base)
      counts <- rowsum(matrix(as.numeric(not_na), nrow(dp_mat), ncol(dp_mat)),
                       dp_base)
      agg <- sums / counts
      agg[counts == 0] <- NA  # avoid NaN where every replicate was NA
      detection_p <- agg
    }
    message("Post-collapse detection_p shape: ",
            nrow(detection_p), " x ", ncol(detection_p))
  }

  detection_p_df = as.data.frame(detection_p)
  detection_p_df <- detection_p_df[, match(sample_sheet$Sentrix_ID, colnames(detection_p_df))]
  message("Number of rows in detection_p: ", nrow(detection_p)) 
  
  # Add predicted sex
  message("Predicting sex based on methylation patterns...")
  pred_sex <- NULL 
  mset_raw <- preprocessRaw(rgset)
  gset <- mapToGenome(mset_raw)
  pred_sex <- getSex(gset)

  # Extract predicted sex values from S4 object. Also keep the X/Y median
  # intensities (xMed/yMed, log2) — they drive minfi's call and, retained here,
  # let the karyotype step flag Loss-of-Y (single X by methylation but depleted
  # Y intensity, common in tumours) instead of discarding them.
  if (is(pred_sex, "DataFrame") || is(pred_sex, "data.frame")) {
    sample_info$pred_sex   <- as.character(pred_sex$predictedSex)
    sample_info$Minfi_xMed <- as.numeric(pred_sex$xMed)
    sample_info$Minfi_yMed <- as.numeric(pred_sex$yMed)
  } else {
    sample_info$pred_sex <- as.character(pred_sex)
  }

  fwrite(detection_p_df, file=file.path(output_dir, "detection_p.txt"), row.names=TRUE, sep="\t")
  m_values = BetaValueToMValue(beta)

  # Resolve the concrete array type once (the param may still be "auto").
  resolved_array_type <- if (array_type == "auto") determine_array_type(rgset) else array_type

  # GCT bisulfite-conversion control (sesame). Aligned to sample_info$Sample_ID
  # (= Sentrix_ID), the same per-sample identifier used elsewhere. Wrapped so a
  # GCT failure can never break preprocessing.
  gct <- tryCatch(
    compute_gct_scores(sample_sheet$Basename, sample_info$Sample_ID,
                       resolved_array_type, bpparam),
    error = function(e) {
      warning("compute_gct_scores() errored (", conditionMessage(e),
              "); skipping GCT table.")
      NULL
    }
  )

  message("Done Preprocessing!")
  write.table(sample_info, file=file.path(output_dir, "sample_info.txt"), row.names=TRUE, sep="\t")

  # Create result object
  result <- list(
    rgset = rgset,
    beta = beta,
    m_values = m_values,
    detection_p = detection_p_df,
    sample_info = sample_info,
    gct = gct,
    array_type = resolved_array_type,
    normalization = normalization
    )

  
  return(result)
}

#' Write beta values to file
#'
#' @param beta Beta values matrix
#' @param file_path Output file path
write_beta_values <- function(beta, file_path) {
  message("Writing beta values to file: ", file_path)
  beta_df <- as.data.frame(beta)
  beta_df <- cbind(data.frame(ProbeID = rownames(beta)), beta_df)
  
  # Write to file using data.table for better performance
  fwrite(beta_df, file = file_path, sep = "\t", quote = FALSE)
}

#' Read beta values from file
#'
#' @param file_path Input file path
#' @return Beta values matrix
read_beta_values <- function(file_path) {
  message("Reading beta values from file: ", file_path)
  beta_df <- fread(file_path, data.table = FALSE)
  rownames(beta_df) <- beta_df$ProbeID
  beta_df$ProbeID <- NULL
  return(as.matrix(beta_df))
}
