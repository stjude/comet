# CNV analysis module for methylation array analysis pipeline

#' Extract clean sample ID from full path
#'
#' @param sample_id Sample ID that may contain full path
#' @return Clean sample ID (basename only)
extract_clean_sample_id <- function(sample_id) {
  # Remove path information and extract just the basename
  # Handle both forward and backslashes
  clean_id <- basename(sample_id)
  
  # If it's still a path-like string, try to extract the last meaningful part
  if (grepl("/", clean_id)) {
    parts <- strsplit(clean_id, "/")[[1]]
    clean_id <- parts[length(parts)]
  }
  
  return(clean_id)
}

#' Fix column names that have been mangled by R when reading files
#'
#' @param col_names Vector of column names to fix
#' @return Vector of cleaned column names
fix_column_names <- function(col_names) {
  # Remove path prefixes and clean up column names
  cleaned_names <- sapply(col_names, function(name) {
    # Expected column names for CNV segments
    expected_cols <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "seg.median")
    
    # Check if any expected column name is at the end of the current name
    for (expected_col in expected_cols) {
      # Create pattern to match the expected column at the end
      pattern <- paste0("\\.", expected_col, "$")
      if (grepl(pattern, name)) {
        return(expected_col)
      }
    }
    
    # If no expected column found, check if the name ends with expected column (without preceding dot)
    for (expected_col in expected_cols) {
      if (name == expected_col) {
        return(expected_col)
      }
    }
    
    # Fallback: if the name contains path information, try to extract meaningful part
    if (grepl("\\.", name)) {
      parts <- strsplit(name, "\\.")[[1]]
      # Look for compound names like "seg.mean", "loc.start", etc.
      if (length(parts) >= 2) {
        # Check last two parts for compound names
        last_two <- paste(parts[(length(parts)-1):length(parts)], collapse = ".")
        if (last_two %in% expected_cols) {
          return(last_two)
        }
      }
      # Otherwise take the last part
      return(parts[length(parts)])
    }
    
    # Return the original name if no processing needed
    return(name)
  })
  
  # Return as vector without names
  return(as.vector(cleaned_names))
}

#' Run copy number variation analysis
#'
#' @param rgset RGChannelSet object
#' @param sample_info Sample information data frame
#' @param references Reference sample sheet or NULL to use internal references
#' @param method CNV calling method (conumee, ChAMP, cnAnalysis450k)
#' @param array_type Array type (450K, EPIC, EPICv2)
#' @param output_dir Output directory for CNV data (segments, RData)
#' @param plots_dir  Output directory for CNV plots. Defaults to
#'                   \code{file.path(output_dir, "plots")} for backward compat.
#' @param threads Number of CPU threads to use
#' @param gain_threshold seg.mean cutoff above which a segment is called a gain
#' @param loss_threshold seg.mean cutoff below which a segment is called a loss
#'                       (signed; e.g. -0.2)
#' @return List containing CNV results
run_cnv_analysis <- function(rgset, sample_info,
                           references = NULL,
                           method = "conumee",
                           array_type = "EPICv2",
                           output_dir = ".",
                           plots_dir  = NULL,
                           threads = 4,
                           gain_threshold =  0.18,
                           loss_threshold = -0.20) {

  message(paste("Running CNV analysis using", method, "method..."))

  # Normalise signs so callers can pass either form.
  gain_threshold <-  abs(gain_threshold)
  loss_threshold <- -abs(loss_threshold)

  if (is.null(plots_dir)) plots_dir <- file.path(output_dir, "plots")
  cnv_seg_dir <- file.path(output_dir, "segments")
  dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)
  dir.create(cnv_seg_dir, showWarnings = FALSE, recursive = TRUE)

  # Select method
  if (method == "conumee") {
    result <- run_conumee_cnv(rgset, sample_info, references, plots_dir, cnv_seg_dir,
                              threads, array_type,
                              gain_threshold = gain_threshold,
                              loss_threshold = loss_threshold)
  } else if (method == "ChAMP") {
    result <- run_champ_cnv(rgset, sample_info, references, plots_dir, cnv_seg_dir, threads, array_type)
  } else {
    stop("Unsupported CNV method: ", method)
  }
  
  return(result)
}

#' Run conumee CNV analysis
#'
#' @param rgset RGChannelSet object
#' @param sample_info Sample information data frame
#' @param references Reference sample sheet or NULL to use internal references
#' @param plots_dir Output directory for CNV plots
#' @param seg_dir Output directory for CNV segment files
#' @param threads Number of CPU threads to use
#' @return List containing CNV results
run_conumee_cnv <- function(rgset, sample_info, references, plots_dir, seg_dir,
                            threads, array_type = "EPICv2",
                            gain_threshold =  0.18,
                            loss_threshold = -0.20) {

  # conumee2 + yamapData are attached here (not at the top of
  # methylation_pipeline.R) so non-CNV steps stay runnable when these
  # packages are missing. Both are required for this step — fail loudly
  # with a message that points the user at setup.R rather than letting
  # `library()` raise the bare "no package called X" error.
  for (pkg in c("conumee2", "yamapData")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is required for conumee CNV analysis but is ",
        "not installed.\n  Re-run setup.R from the project root to provision ",
        "the R library:\n    Rscript setup.R"
      )
    }
  }
  suppressPackageStartupMessages({
    library(conumee2)
    library(yamapData)
  })

  # Normalize array_type to conumee2-accepted values (case-sensitive: 450k, EPIC, EPICv2, mouse)
  array_type <- switch(toupper(array_type),
    "450K"   = "450k",
    "EPIC"   = "EPIC",
    "EPICV2" = "EPICv2",
    "MOUSE"  = "mouse",
    array_type  # pass through unknown values and let conumee2 error
  )

  # Load reference samples if provided
  if (!is.null(references)) {
    message("Loading reference samples...")
    ref_rgset <- read.metharray.exp(base = NULL, targets = references, recursive = TRUE)
    #ref_rgset@annotation <- rgset@annotation
    
    # Preprocess reference samples
    ref_mset <- preprocessNoob(ref_rgset)
    ref_controls <- CNV.load(ref_mset)
  } else {
    message("Using internal reference samples...")
    # Load the yamapData internal reference dataset
    ref
    ref_mset <- preprocessNoob(ref)
    ref_controls <- CNV.load(ref_mset)
  }
  
  # Create conumee annotation
  message("Creating CNV annotation...")
  
  # Load required data from conumee2 package
  if (requireNamespace("conumee2", quietly = TRUE)) {
    # Load annotation data from conumee2
    data("exclude_regions")
    data("detail_regions")
      
    # For EPICv2 samples, create annotation using both EPIC and EPICv2 to get probe overlap
    # with the internal EPIC reference controls. For 450k, use only 450k.
    anno_array_type <- if (array_type == "450k") "450k" else c("EPIC", array_type)
    anno <- conumee2::CNV.create_anno(
      exclude_regions = exclude_regions,
      bin_minprobes = 15,
      bin_minsize = 50000,
      array_type = anno_array_type,
      detail_regions = detail_regions
    )
          
  }
  anno
  # Check if annotation was created successfully
  if (is.null(anno) || length(anno@probes) == 0) {
    stop("Failed to create CNV annotation or annotation contains 0 probes. ",
         "This may be due to incompatible array types or missing annotation data.")
  }
  
  message("Annotation object created with ", length(anno@probes), " probes")
  message("Annotation class: ", class(anno))

  # Get common probes between annotation and reference controls
  anno_probes <- rownames(as.data.frame(anno@probes))
  message("Annotation_probes: ", head(anno_probes))
  ref_probes <- rownames(ref_controls@intensity)
  message("Number of probes in ref_probes: ", head(ref_probes))
  common_probes <- intersect(anno_probes, ref_probes)
  message("Common probes between annotation and reference controls: ", length(common_probes))
  
  if (length(common_probes) == 0) {
    stop("No common probes found between annotation and reference controls")
  }
  
  # Subset reference controls to common probes (drop = FALSE preserves data.frame)
  ref_controls@intensity <- ref_controls@intensity[common_probes, , drop = FALSE]
  message("Reference controls subset to ", nrow(ref_controls@intensity), " common probes")

  # Subset annotation to the same common probes so that query, reference, and
  # annotation all share exactly the same probe set (required by CNV.fit)
  anno@probes <- anno@probes[names(anno@probes) %in% common_probes]
  message("Annotation subset to ", length(anno@probes), " common probes")
  
  # Process each sample
  message("Processing samples...")
  sample_ids <- sample_info$Sample_ID
  n_samples <- length(sample_ids)
  
  # Use serial processing to avoid rgset serialization issues
  message("Processing samples in serial mode to avoid serialization issues...")
  results <- lapply(seq_len(n_samples), function(i) {
    process_conumee_sample(rgset, sample_ids[i], ref_controls, anno, plots_dir, seg_dir,
                           gain_threshold = gain_threshold,
                           loss_threshold = loss_threshold)
  })
  
  # Combine segment files for IGV
  message("Combining segment files...")
  combined_segments <- tryCatch({
    # Extract segments from each result and add sample ID
    segments_list <- lapply(results, function(x) {
      segs <- x$segments
      if (is.null(segs) || nrow(segs) == 0) {
        message("Warning: No segments found for sample ", x$sample_id)
        return(NULL)
      }
      # Ensure it's a data frame
      if (!is.data.frame(segs)) {
        segs <- as.data.frame(segs)
      }
      
      # Fix column names that may have been mangled by R
      colnames(segs) <- fix_column_names(colnames(segs))
      
      # Standardize column names to match expected format
      # Common variations in column names from different sources
      col_mapping <- c(
        "chrom" = "chrom",
        "chr" = "chrom",
        "chromosome" = "chrom",
        "loc.start" = "loc.start", 
        "start" = "loc.start",
        "loc.end" = "loc.end",
        "end" = "loc.end",
        "num.mark" = "num.mark",
        "num.markers" = "num.mark",
        "nprobes" = "num.mark",
        "seg.mean" = "seg.mean",
        "mean" = "seg.mean",
        "log2" = "seg.mean"
      )
      
      # Rename columns based on mapping
      current_names <- colnames(segs)
      for (old_name in names(col_mapping)) {
        if (old_name %in% current_names) {
          colnames(segs)[colnames(segs) == old_name] <- col_mapping[old_name]
        }
      }
      
      # Add sample ID column (use clean sample ID)
      segs$ID <- extract_clean_sample_id(x$sample_id)
      
      # Ensure we have the required columns, add defaults if missing
      required_cols <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
      for (col in required_cols) {
        if (!col %in% colnames(segs)) {
          if (col == "num.mark") {
            segs[[col]] <- 1  # Default number of markers
          } else if (col == "seg.mean") {
            segs[[col]] <- 0  # Default log2 ratio
          } else {
            segs[[col]] <- NA  # Default for other columns
          }
        }
      }
      
      # Return only the required columns in the correct order
      return(segs[, required_cols])
    })
    
    # Remove NULL entries
    segments_list <- segments_list[!sapply(segments_list, is.null)]
    
    if (length(segments_list) == 0) {
      message("Warning: No valid segments found for any samples")
      return(data.frame())
    }
    
    # Combine all segments
    do.call(rbind, segments_list)
  }, error = function(e) {
    message("Error combining segments: ", e$message)
    return(data.frame())
  })
  
  # Segments should now have standardized column names and order
  if (nrow(combined_segments) > 0) {
    message("Successfully combined segments from ", length(results), " samples")
    message("Combined segments dimensions: ", nrow(combined_segments), " x ", ncol(combined_segments))
  } else {
    message("Warning: No segments to combine")
  }
  
  # Save combined segments
  combined_seg_file <- file.path(seg_dir, "combined_cnv_segments.seg")
  write.table(combined_segments, combined_seg_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Return results
  return(list(
    method = "conumee",
    segments = combined_segments,
    sample_results = results,
    combined_seg_file = combined_seg_file
  ))
}

#' Process a single sample with conumee
#'
#' @param rgset RGChannelSet object
#' @param sample_id Sample ID
#' @param ref_controls Reference controls
#' @param anno CNV annotation
#' @param plots_dir Output directory for CNV plots
#' @param seg_dir Output directory for CNV segment files
#' @return List containing sample CNV results
process_conumee_sample <- function(rgset, sample_id, ref_controls, anno, plots_dir, seg_dir,
                                   gain_threshold =  0.18,
                                   loss_threshold = -0.20) {
  message(paste("Processing sample:", sample_id))
  
  # Extract clean sample ID for file naming and output
  clean_sample_id <- extract_clean_sample_id(sample_id)
  message(paste("Clean sample ID:", clean_sample_id))
  
  # Subset RGSet for the current sample
  # 1. Exact match
  sample_idx <- which(colnames(rgset) == sample_id)

  # 2. Fallback: extract Sentrix barcode_position (e.g. 201465940014_R01C01)
  #    and match against rgset colnames. Handles cases where sample_info$Sample_ID
  #    carries extra prefixes (e.g. GSM3735546_GSM3735546_201465940014_R01C01).
  if (length(sample_idx) == 0) {
    sentrix_hit <- regmatches(sample_id,
                              regexpr("[0-9]{9,12}_R[0-9]{2}C[0-9]{2}", sample_id))
    if (length(sentrix_hit) == 1 && nchar(sentrix_hit) > 0) {
      sample_idx <- grep(sentrix_hit, colnames(rgset), fixed = TRUE)
      if (length(sample_idx) > 0)
        message(sprintf("  Matched '%s' -> rgset column '%s' via Sentrix ID",
                        sample_id, colnames(rgset)[sample_idx[1]]))
    }
  }

  if (length(sample_idx) == 0)
    stop("Sample '", sample_id, "' not found in rgset.\n",
         "  rgset colnames (first 5): ",
         paste(head(colnames(rgset), 5), collapse = ", "))

  if (length(sample_idx) > 1) {
    message("  Warning: multiple rgset columns matched '", sample_id,
            "'; using first: ", colnames(rgset)[sample_idx[1]])
    sample_idx <- sample_idx[1]
  }
  
  # Use proper subsetting for S4 objects
  sub_rgset <- rgset[, sample_idx]
  message("Number of probes in a sample ", nrow(sub_rgset))
  # Preprocess the sample - use preprocessSWAN to avoid getSex issues
  sub_mset <- tryCatch({
    preprocessSWAN(sub_rgset)
  }, error = function(e) {
    message("SWAN normalization failed, falling back to raw preprocessing: ", e$message)
    preprocessRaw(sub_rgset)
  })
  
  # Load CNV data
  cnv_data <- CNV.load(sub_mset)

  # EPICv2 samples come out of CNV.load with per-replicate rownames
  # (cg25324105_BC11, cg25324105_TC11, ...) while the CNV annotation
  # and the internal yamapData reference use plain EPIC cg-IDs.
  # Collapse replicates to one row per CpG (mean intensity across
  # replicates) so the three probe spaces align for CNV.fit.
  cnv_probes_raw <- rownames(cnv_data@intensity)
  if (any(grepl("_(BC|TC)\\d+$", utils::head(cnv_probes_raw, 200)))) {
    message("  Collapsing EPICv2 replicate probes in sample cnv_data (mean)...")
    base_ids <- sub("_.*$", "", cnv_probes_raw)
    int_mat  <- as.matrix(cnv_data@intensity)
    not_na   <- !is.na(int_mat)
    int_mat[!not_na] <- 0
    sums   <- rowsum(int_mat, base_ids)
    counts <- rowsum(matrix(as.numeric(not_na), nrow(int_mat), ncol(int_mat)),
                     base_ids)
    agg <- sums / counts
    agg[counts == 0] <- NA  # avoid NaN when every replicate was NA
    # CNV.data's @intensity slot requires a data.frame, not a matrix.
    cnv_data@intensity <- as.data.frame(agg)
    message("  After collapse: ", nrow(cnv_data@intensity), " probes")
  }

  # Subset query to exactly the probes present in the annotation (which has already
  # been aligned to the reference controls). CNV.fit requires query, reference, and
  # annotation to share the same probe set.
  anno_probes <- names(anno@probes)
  cnv_probes  <- rownames(cnv_data@intensity)
  common_probes <- intersect(anno_probes, cnv_probes)

  message("Sample ", sample_id, " - CNV data probes: ", length(cnv_probes))
  message("Sample ", sample_id, " - Common probes with annotation: ", length(common_probes))

  if (length(common_probes) == 0) {
    message("  rgset annotation: ",
            paste(rgset@annotation, collapse = " / "))
    message("  cnv_probes head: ",
            paste(utils::head(cnv_probes, 6), collapse = ", "))
    message("  anno_probes head: ",
            paste(utils::head(anno_probes, 6), collapse = ", "))
    message("  cnv_probes has '_' suffix in any of first 100: ",
            any(grepl("_", utils::head(cnv_probes, 100))))
    message("  anno_probes has '_' suffix in any of first 100: ",
            any(grepl("_", utils::head(anno_probes, 100))))
    stop("No common probes found between sample and annotation for sample ", sample_id)
  }

  # Subset cnv_data to annotation probes
  # Use drop = FALSE to preserve data.frame structure when cnv_data has one column
  cnv_data@intensity <- cnv_data@intensity[common_probes, , drop = FALSE]
  
  # Fit CNV data against reference controls
  cnv_fit <- CNV.fit(cnv_data, ref_controls, anno)
  
  # Perform binning
  cnv_bin <- CNV.bin(cnv_fit)
  
  # Add details to regions of interest
  cnv_detail <- CNV.detail(cnv_bin)
  
  # Segment the data
  cnv_segment <- CNV.segment(cnv_detail)
  
  # Plot genome-wide CNV profile.
  # `cols` is a 5-stop palette interpolated low→high (loss→neutral→gain):
  # dark teal → light teal → neutral grey → light brown → dark brown.
  # Matches the CNV heatmap's teal/brown diverging scheme.
  pdf_file <- file.path(plots_dir, paste0(clean_sample_id, "_cnv_profile.pdf"))
  pdf(pdf_file, width = 10, height = 5)
  CNV.genomeplot(
    cnv_segment,
    main = clean_sample_id,
    cols = c("#018571", "#80cdc1", "#f5f5f5", "#dfc27d", "#a6611a")
  )
  # Gain/loss threshold lines: neutral grey so they don't clash with the
  # new teal/brown palette. Asymmetric — gain at +gain_threshold,
  # loss at loss_threshold (already signed negative).
  abline(h = c(loss_threshold, gain_threshold), lty = 2, col = "grey50")
  dev.off()
  
  # Save segment results
  seg_file <- file.path(seg_dir, paste0(clean_sample_id, "_cnv_segments.seg"))
  CNV.write(cnv_segment, what = "segments", file = seg_file)
  
  # Extract segments as data frame for downstream processing
  segments_df <- tryCatch({
    # Try to extract segments from the S4 object
    if (class(cnv_segment)[1] == "CNV.analysis") {
      seg_data <- cnv_segment@seg
      # Ensure it's a data frame with proper column names
      if (!is.data.frame(seg_data)) {
        seg_data <- as.data.frame(seg_data)
      }
      seg_data
    } else {
      # Fallback: read the segments file we just wrote
      seg_data <- read.table(seg_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      # Fix column names that may have been mangled by R
      colnames(seg_data) <- fix_column_names(colnames(seg_data))
      seg_data
    }
  }, error = function(e) {
    message("Warning: Could not extract segments from CNV object, reading from file: ", e$message)
    # Read from the file we just wrote
    seg_data <- read.table(seg_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Fix column names that may have been mangled by R
    colnames(seg_data) <- fix_column_names(colnames(seg_data))
    seg_data
  })
  
  # Return results
  list(
    sample_id = clean_sample_id,  # Use clean sample ID
    original_sample_id = sample_id,  # Keep original for reference
    segments = segments_df,
    bins = cnv_bin,
    segments_file = seg_file,
    plot_file = pdf_file
  )
}

#' Run ChAMP CNV analysis
#'
#' @param rgset RGChannelSet object
#' @param sample_info Sample information data frame
#' @param references Reference sample sheet or NULL
#' @param plots_dir Output directory for CNV plots
#' @param seg_dir Output directory for CNV segment files
#' @param threads Number of CPU threads to use
#' @return List containing CNV results
run_champ_cnv <- function(rgset, sample_info, references, plots_dir, seg_dir, threads, array_type = "EPICv2") {
  
  # Check if ChAMP is available
  if (!requireNamespace("ChAMP", quietly = TRUE)) {
    stop("Package 'ChAMP' is required for ChAMP CNV analysis")
  }
  
  message("Preprocessing for ChAMP CNV analysis...")
  
  # Preprocess data for ChAMP
  mset <- preprocessFunnorm(rgset)
  beta <- getBeta(mset)
  
  # Run ChAMP CNV
  message("Running ChAMP CNA analysis...")
  champ_cnv <- ChAMP::champ.CNA(
    beta = beta,
    pheno = sample_info$Sample_Name,
    control = NULL,  # Use all samples as controls
    sampleCNA = 1:ncol(beta),  # Analyze all samples
    arraytype = gsub("IlluminaHuman", "", rgset@annotation["array"]),
    plotOutput = TRUE,
    resultsDir = plots_dir,
    suffix = ""
  )
  
  # Get segments
  segments <- do.call(rbind, lapply(names(champ_cnv), function(sample) {
    segs <- champ_cnv[[sample]]$output
    segs$ID <- sample
    return(segs)
  }))
  
  # Rename columns to match IGV format
  colnames(segments)[colnames(segments) == "chrom"] <- "chrom"
  colnames(segments)[colnames(segments) == "loc.start"] <- "loc.start"
  colnames(segments)[colnames(segments) == "loc.end"] <- "loc.end"
  colnames(segments)[colnames(segments) == "num.mark"] <- "num.mark"
  colnames(segments)[colnames(segments) == "seg.mean"] <- "seg.mean"
  
  # Ensure proper column order for IGV
  segments <- segments[, c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")]
  
  # Save combined segments
  combined_seg_file <- file.path(seg_dir, "champ_cnv_segments.seg")
  write.table(segments, combined_seg_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Return results
  return(list(
    method = "ChAMP",
    segments = segments,
    sample_results = champ_cnv,
    combined_seg_file = combined_seg_file
  ))
}


#' Generate CNV frequency plot
#'
#' @param segments Data frame with CNV segments
#' @param gain_threshold seg.mean cutoff above which a segment counts as a gain
#' @param loss_threshold seg.mean cutoff below which a segment counts as a loss
#'                       (signed; e.g. -0.2)
#' @param threshold Deprecated. Legacy single absolute cutoff applied
#'                  symmetrically; if supplied, overrides gain/loss defaults
#'                  with ±abs(threshold). Kept for back-compat with callers
#'                  that haven't been updated yet.
#' @param output_dir Output directory for plots
#' @return Path to frequency plot
generate_cnv_frequency_plot <- function(segments,
                                        gain_threshold =  0.18,
                                        loss_threshold = -0.20,
                                        threshold = NULL,
                                        output_dir = ".") {

  # Back-compat: if a caller passed the old single threshold, apply it
  # symmetrically and ignore the gain/loss defaults.
  if (!is.null(threshold)) {
    gain_threshold <-  abs(threshold)
    loss_threshold <- -abs(threshold)
  }
  gain_threshold <-  abs(gain_threshold)
  loss_threshold <- -abs(loss_threshold)
  
  message("Generating CNV frequency plot...")
  
  # Source the freqplot function from the provided code
  source_freqplot_functions()
  
  # Check if segments is empty
  if (is.null(segments) || nrow(segments) == 0) {
    message("No segments provided for frequency plot. Skipping plot generation.")
    return(NULL)
  }
  
  # Fix column names that may have been mangled by R
  colnames(segments) <- fix_column_names(colnames(segments))
  message("Column names after fixing: ", paste(colnames(segments), collapse = ", "))
  
  # Convert segments to format expected by freqplot
  # Need to ensure columns are in the right format
  segments$chrom <- as.numeric(gsub("chr", "", segments$chrom))
  segments$loc.start <- as.numeric(segments$loc.start)
  segments$loc.end <- as.numeric(segments$loc.end)
  segments$num.mark <- as.numeric(segments$num.mark)
  segments$seg.mean <- as.numeric(segments$seg.mean)
  
  # Remove rows with NA values in critical columns
  critical_cols <- c("chrom", "loc.start", "loc.end", "seg.mean")
  na_rows <- rowSums(is.na(segments[, critical_cols])) > 0
  
  if (any(na_rows)) {
    message("Removing ", sum(na_rows), " rows with NA values in critical columns")
    segments <- segments[!na_rows, ]
  }
  
  # Check if we still have segments after filtering
  if (nrow(segments) == 0) {
    message("No valid segments remaining after filtering NAs. Skipping plot generation.")
    return(NULL)
  }
  
  # Ensure chromosome values are valid (remove any remaining NAs or invalid values)
  segments <- segments[!is.na(segments$chrom) & segments$chrom >= 1 & segments$chrom <= 24, ]
  
  if (nrow(segments) == 0) {
    message("No valid segments with valid chromosome numbers. Skipping plot generation.")
    return(NULL)
  }
  
  # Generate plot
  pdf_path <- file.path(output_dir, "cnv_frequency_plot.pdf")
  pdf(pdf_path, width = 12, height = 5)
  
  # Call the freqplot function
  freqplot(segments,
           gain_threshold = gain_threshold,
           loss_threshold = loss_threshold,
           plot.title = "CNV Frequency Plot")
  
  dev.off()
  
  return(pdf_path)
}

#' Source freqplot functions
#'
#' This function sources the freqplot functions from the provided code
source_freqplot_functions <- function() {
  
  # Define the freqplot function and assign to global environment.
  # gain_threshold / loss_threshold are asymmetric: a segment is a gain iff
  # its seg.mean is strictly greater than gain_threshold, and a loss iff
  # strictly less than loss_threshold (loss_threshold is signed negative).
  freqplot <<- function(segs.data,
                        gain_threshold =  0.18,
                        loss_threshold = -0.20,
                        plot.title = "Copy Number Frequency") {
    gain_threshold <-  abs(gain_threshold)
    loss_threshold <- -abs(loss_threshold)

    segs.mtx <- segs.as.mtx(segs.data, anns=NULL)

    mean.pos.col <- which(colnames(segs.mtx) == "mean.pos")
    X <- as.matrix(segs.mtx[, -(1:mean.pos.col)])
    n.ints <- nrow(segs.mtx)
    int.size <- segs.mtx[, "loc.end"] - segs.mtx[, "loc.start"] + 1
    int.size <- int.size / max(int.size)
    g.start <- c(0, cumsum(int.size[-n.ints]))
    g.end <- cumsum(int.size)
    n.subj <- ncol(X)

    chr.num <- cumsum(c(1, segs.mtx$chrom[-1] != segs.mtx$chrom[-n.ints]))
    bg.col <- c("lightgray", "white")[chr.num %% 2 + 1]

    nloss <- rep(0, n.ints)
    ngain <- rep(0, n.ints)

    for(i in 1:n.subj) {
      seg.col <- bg.col
      vals <- X[, i]
      seg.gain <- !is.na(vals) & vals > gain_threshold
      seg.loss <- !is.na(vals) & vals < loss_threshold
      nloss <- nloss + seg.loss
      ngain <- ngain + seg.gain
    }
    
    max.chng <- max(ngain, nloss, na.rm=TRUE)
    
    par(xpd=TRUE)
    plot(c(0, g.end[n.ints]), c(-0.8, 0.8),  
         type="n", axes=FALSE,
         xlab="", ylab="")
    axis(side=2, at=c(-1, -0.5, 0, 0.5, 1),
         labels=c("100", "50", "0", "50", "100"), line=-1.8, las=1, cex=2.2)
    
    rect(g.start, -1, g.end, 1, col=bg.col, border=bg.col)
    
    rect(g.start, 0, g.end, 0 + (ngain/n.subj), col="#a6611a", border="#a6611a")

    rect(g.start, 0, g.end, 0 - (nloss/n.subj), col="#018571", border="#018571")
    
    mtext(plot.title, side=3, at=g.end[n.ints]/2, cex=1.3)
    mtext("% of Gain or Loss", side=2, at=0, cex=1.0)
    
    chr.start <- c(g.start[1], g.start[c(segs.mtx$chrom[-1] != segs.mtx$chrom[-n.ints])])
    if (length(chr.start) >= 22) {
      chr.start <- chr.start[-22]
    }
    
    chr.end <- c(g.end[c(segs.mtx$chrom[-1] != segs.mtx$chrom[-n.ints])], g.end[n.ints])
    if (length(chr.end) >= 22) {
      chr.end <- chr.end[-22]
    }
    
    chr.labels <- c(1:21, "X", "Y")
    chr.labels <- chr.labels[1:min(length(chr.start), length(chr.labels))]
    
    text((chr.start + chr.end)/2, -1.05, chr.labels, cex=0.6)
  }
  
  # Define supporting functions and assign to global environment
  segs.as.mtx <<- function(segs, value="seg.mean", threshold=NA, anns=NULL, max.per.row=25) {
    options(stringsAsFactors=FALSE)
    
    y <- segs[, "seg.mean"]
    if (!is.na(threshold)) y <- sign(y) * (abs(y) > abs(threshold))
    
    seg <- data.frame(
      samp.id = segs[, "ID"],
      chrom = segs[, "chrom"],
      loc.start = segs[, "loc.start"],
      loc.end = segs[, "loc.end"],
      value = y
    )
    colnames(seg)[1] <- "samp.id"
    colnames(seg)[5] <- "value"
    
    # Determine set of unique intervals defined by start and end locations of the segments
    uniq.ints <- define.unique.intervals(seg)   
    n.ints <- nrow(uniq.ints)
    
    # Get chr indices
    int.new.chr <- which(uniq.ints$chrom[-1] != uniq.ints$chrom[-n.ints])
    int.chr.ind <- cbind(c(1, int.new.chr + 1), c(int.new.chr, n.ints))
    
    ord <- order(seg$chrom, seg$samp.id, seg$loc.start)
    seg <- seg[ord, ]
    n.segs <- nrow(seg)
    
    seg.new.chr <- which(seg$chrom[-1] != seg$chrom[-n.segs])
    seg.chr.ind <- cbind(c(1, seg.new.chr + 1), c(seg.new.chr, n.segs))
    
    if (nrow(seg.chr.ind) != nrow(int.chr.ind)) stop("chr mismatch")
    if (any(seg$chrom[seg.chr.ind[, 1]] != uniq.ints$chrom[int.chr.ind[, 1]])) stop("Chr mismatch")
    if (any(seg$chrom[seg.chr.ind[, 2]] != uniq.ints$chrom[int.chr.ind[, 2]])) stop("Chr mismatch")  
    
    usamp <- unique(seg$samp.id)
    usamp.names <- make.names(usamp)
    nsamp <- length(usamp)
    X <- matrix(NA, n.ints, nsamp)
    colnames(X) <- usamp.names
    
    for (i in 1:nrow(int.chr.ind)) {
      int.chr.index <- (int.chr.ind[i, 1]:int.chr.ind[i, 2])
      chr.ints <- uniq.ints[int.chr.index, 2:3]
      chr.segs <- seg[seg.chr.ind[i, 1]:seg.chr.ind[i, 2], ]
      chr.ints.index <- 1:nrow(chr.ints)
      
      seg.int.index1 <- approx(chr.ints[, 2],
                              chr.ints.index,
                              xout=chr.segs$loc.start, f=1,
                              method="constant", rule=2:1)$y
      
      seg.int.index2 <- approx(chr.ints[, 1],
                              chr.ints.index,
                              xout=chr.segs$loc.end,
                              method="constant",
                              rule=1:2)$y
      
      nsegs.chr <- nrow(chr.segs)
      for (j in 1:nsegs.chr) {
        samp.id.index <- which(is.element(usamp, chr.segs$samp.id[j]))
        seg.int.index <- int.chr.index[seg.int.index1[j]:seg.int.index2[j]]
        X[seg.int.index, samp.id.index] <- chr.segs$value[j]
      }
    }
    
    # Remove rows with all NAs
    all.na <- (rowSums(is.na(X)) == nsamp)
    X <- X[!all.na, ]
    uniq.ints <- uniq.ints[!all.na, ]
    
    # Combine identical adjacent rows
    x.nints <- nrow(X)
    x.chng <- rowSums(abs(X[-1, ] - X[-x.nints, ])) > 0
    x.chng[is.na(x.chng)] <- TRUE
    chr.chng <- (uniq.ints[-1, 1] != uniq.ints[-x.nints, 1])
    new.row <- which(chr.chng | x.chng)
    
    new.chr <- uniq.ints[c(1, new.row + 1), 1]
    new.start <- uniq.ints[c(1, new.row + 1), 2]
    new.end <- uniq.ints[c(new.row, x.nints), 3]
    
    X <- X[c(1, new.row + 1), ]
    new.ints <- cbind.data.frame(
      chrom = new.chr,
      loc.start = new.start,
      loc.end = new.end
    )
    
    mean.val <- rowMeans(X, na.rm=TRUE)
    mean.abs <- rowMeans(abs(X), na.rm=TRUE)
    mean.neg <- rowMeans((X < 0) * X, na.rm=TRUE)
    mean.pos <- rowMeans((X > 0) * X, na.rm=TRUE)
    
    final.result <- cbind.data.frame(
      new.ints,
      mean.val,
      mean.abs,
      mean.neg,
      mean.pos,
      X
    )
    
    if (!is.null(anns)) {
      ann.result <- ext.ann.segs(final.result, anns, max.per.row)
      final.result <- ann.result
    }
    
    return(final.result)
  }
  
  define.unique.intervals <<- function(seg.data) {
    # Debug information
    message("Debug: define.unique.intervals called with ", nrow(seg.data), " segments")
    
    # Check for NAs in input data
    if (any(is.na(seg.data$chrom))) {
      message("Warning: Found ", sum(is.na(seg.data$chrom)), " NA values in chrom column")
    }
    if (any(is.na(seg.data$loc.start))) {
      message("Warning: Found ", sum(is.na(seg.data$loc.start)), " NA values in loc.start column")
    }
    if (any(is.na(seg.data$loc.end))) {
      message("Warning: Found ", sum(is.na(seg.data$loc.end)), " NA values in loc.end column")
    }
    
    ord <- order(seg.data$chrom, seg.data$loc.start)
    seg.data <- seg.data[ord, ]
    nsegs <- nrow(seg.data)
    
    chr <- seg.data$chrom
    loc.start <- seg.data$loc.start
    loc.end <- seg.data$loc.end
    
    new.chr.index <- which(chr[-1] != chr[-nsegs])
    chr.index <- cbind(
      c(1, new.chr.index + 1),
      c(new.chr.index, nsegs)
    )
    
    final.res <- NULL
    for (i in 1:nrow(chr.index)) {
      this.index <- (chr.index[i, 1]:chr.index[i, 2])
      this.start.end <- seg.data[this.index, c("loc.start", "loc.end")]
      this.res0 <- unique.start.end(this.start.end)
      this.res <- cbind.data.frame(
        chrom = seg.data$chrom[chr.index[i, 1]],
        loc.start = this.res0[, 1],
        loc.end = this.res0[, 2]
      )
      final.res <- rbind.data.frame(final.res, this.res)
    }
    
    return(final.res)
  }
  
  unique.start.end <<- function(st.end) {
    if (any(is.na(st.end))) {
      message("Debug: Found NAs in st.end matrix:")
      message("  Dimensions: ", paste(dim(st.end), collapse = " x "))
      message("  NA positions: ", paste(which(is.na(st.end), arr.ind = TRUE), collapse = ", "))
      message("  First few rows:")
      print(head(st.end))
      stop("NAs not allowed in st.end")
    }
    uniq.start <- unique(st.end[, 1])
    uniq.end <- unique(st.end[, 2])
    sort.ustart <- sort(uniq.start)
    sort.uend <- sort(uniq.end)
    n.end <- length(sort.uend)
    last.end <- sort.uend[n.end]
    new.start <- c(
      sort.ustart,
      sort.uend[-n.end] + 1
    )
    u.new.start <- unique(new.start)
    u.new.start <- sort(u.new.start)
    u.new.end <- c(
      u.new.start[-1] - 1,
      last.end
    )
    res <- cbind(u.new.start, u.new.end)
    return(res)
  }
  
  # Return invisible NULL to avoid printing
  invisible(NULL)
}
