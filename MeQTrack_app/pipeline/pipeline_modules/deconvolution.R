# deconvolution.R
# ---------------------------------------------------------------------------
# Reference-based cell-type deconvolution via the deconvMe package
# (omnideconv/deconvMe), which wraps EpiDISH / Houseman / methylCC /
# methylResolver / MethAtlas for Illumina array data.
#
# Standalone, OPT-IN step (--step deconvolution): never part of --step all,
# never gates Pass_QC. It reloads the preprocessed RGChannelSet and runs
# independently, exactly like reference_projection.
#
# Output: deconv/cell_fractions.csv — a tidy LONG table
#   (method, sample, celltype, value), including the "aggregated" method.
#
# Interpretation vs sesame Leukocyte_Fraction: deconvMe reports the immune
# COMPOSITION (which cell types, in what proportion); sesame's leukocyte
# fraction reports the AMOUNT of immune (a 2-component purity scalar). On
# tumours the blood references force-fit tumour methylation into blood cell
# types, so the summed-immune fraction over-estimates magnitude and is NOT
# comparable to sesame's number — compare composition, not magnitude.
# ---------------------------------------------------------------------------

#' Map our array_type to the deconvMe `array` argument ('450k' or 'EPIC').
#' EPICv2 isn't natively supported by deconvMe -> NA (skip).
.deconv_array <- function(array_type) {
  switch(toupper(as.character(array_type)),
         "450K"   = "450k",
         "EPIC"   = "EPIC",
         "EPICV2" = NA_character_,
         NA_character_)
}

#' Run cell-type deconvolution on a preprocessed run.
#'
#' @param rgset minfi RGChannelSet (from preprocessed_data.RData)
#' @param array_type "450k" / "EPIC" / "EPICv2"
#' @param methods deconvMe method names; default EpiDISH + Houseman (pure-R).
#'   methylcc / methatlas need a Python backend; methatlas is research-only.
#' @param output_dir the run's deconv/ directory
#' @return the tidy long data.frame (also written to deconv/cell_fractions.csv),
#'   or NULL if deconvolution was skipped (deconvMe absent, no rgset, EPICv2).
run_deconvolution <- function(rgset, array_type,
                              methods = c("epidish", "houseman"),
                              output_dir = ".") {
  if (!requireNamespace("deconvMe", quietly = TRUE)) {
    message("deconvMe is not installed - skipping deconvolution. ",
            "Install with: pak::pkg_install('omnideconv/deconvMe')")
    return(NULL)
  }
  if (is.null(rgset)) {
    message("No rgset available - deconvolution needs the preprocessed RGChannelSet.")
    return(NULL)
  }
  arr <- .deconv_array(array_type)
  if (is.na(arr)) {
    message("Cell-type deconvolution not supported for array type '", array_type,
            "' (deconvMe handles 450k/EPIC; EPICv2 would need EPIC conversion) ",
            "- skipping.")
    return(NULL)
  }

  message("Building MethylSet and running deconvMe (",
          paste(methods, collapse = ", "), ") on ", arr, "...")
  mset <- minfi::preprocessRaw(rgset)
  res <- deconvMe::deconvolute_combined(methyl_set = mset, array = arr,
                                        methods = methods)
  res <- as.data.frame(res)

  out_csv <- file.path(output_dir, "cell_fractions.csv")
  utils::write.csv(res, out_csv, row.names = FALSE)
  message(sprintf(
    "Deconvolution complete: %d rows (%d method-set(s) x samples x cell types) -> deconv/cell_fractions.csv",
    nrow(res), length(unique(res$method))))
  res
}
