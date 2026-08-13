#' Resolve a file under the vendored \code{Anno/} annotation dir
#'
#' The pipeline \code{setwd()}s into \code{pipeline/}, so \code{Anno/} sits one
#' level up. Probe a few candidate locations to stay robust to how the module
#' was sourced.
#'
#' @param rel_path Path relative to \code{Anno/} (e.g. "HM450/Clock_Horvath353.rds")
#' @return The first existing absolute/relative path, or NA_character_ if none.
find_anno_file <- function(rel_path) {
  candidates <- c(
    file.path("..", "Anno", rel_path),
    file.path("Anno", rel_path),
    file.path(getwd(), "..", "Anno", rel_path)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[1] else NA_character_
}

#' Locate the vendored Horvath353 epigenetic-clock model
#'
#' The model ships in the project's \code{Anno/} dir (from the Zhou Lab
#' InfiniumAnnotation repo) in sesame \code{predictAge()}-compatible form
#' (\code{intercept}/\code{param$slope}/\code{response2age}). We use the HM450
#' clock for every array type because its probe IDs are base \code{cg} form,
#' matching our (collapsed) beta matrices — including EPICv2, whose betas are
#' collapsed to base CpG IDs upstream.
#'
#' @return The loaded clock model list, or NULL if it cannot be found/read.
load_horvath_model <- function() {
  path <- find_anno_file(file.path("HM450", "Clock_Horvath353.rds"))
  if (is.na(path)) {
    warning("Horvath353 clock model not found under Anno/HM450/; age skipped.")
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) {
    warning("Could not read Horvath353 model: ", conditionMessage(e)); NULL
  })
}

#' Convert a collapsed EPICv2 beta vector into EPIC (v1) probe space
#'
#' Uses the Zhou Lab EPICv2->EPIC map (\code{Anno/EPICv2/EPICv2ToEPIC_map.tsv.gz},
#' a slim 3-column form of \code{EPICv2ToEPIC_conversion.tsv}). In that table the
#' EPICv2 base probe ID (suffix stripped) always equals the EPIC1 ID, and our
#' EPICv2 betas are already collapsed to base \code{cg} IDs — so conversion is a
#' reliability-filtered probe subset: keep probes whose EPIC1 target is
#' well-calibrated (\code{big_delta == FALSE}), drop the rest. This unlocks the
#' EPIC-only \code{estimateLeukocyte} for EPICv2 samples.
#'
#' @return Character vector of reliable EPIC probe IDs, or NULL if map missing.
#' Load the vendored EPICv2 Infinium-I C/T-extension probe lists for GCT.
#' Returns list(extC, extT) of EPICv2 probe IDs, or NULL if the file is absent.
#' Derived from the Zhou Lab EPICv2 manifest (type-I, nextBase R->extC, A->extT);
#' validated to reproduce sesame's native bisConversionControl on EPIC.
load_epicv2_ext_probes <- function() {
  path <- find_anno_file(file.path("EPICv2", "EPICv2.typeI.ext.rds"))
  if (is.na(path)) {
    warning("EPICv2 GCT ext-probe list not found under Anno/EPICv2/; ",
            "EPICv2 GCT skipped.")
    return(NULL)
  }
  ext <- tryCatch(readRDS(path), error = function(e) {
    warning("Could not read EPICv2 ext-probe list: ", conditionMessage(e)); NULL
  })
  if (is.null(ext) || !all(c("extC", "extT") %in% names(ext))) return(NULL)
  ext
}

load_epicv2_reliable_epic_probes <- function() {
  path <- find_anno_file(file.path("EPICv2", "EPICv2ToEPIC_map.tsv.gz"))
  if (is.na(path)) {
    warning("EPICv2->EPIC map not found under Anno/EPICv2/; ",
            "EPICv2 leukocyte estimation skipped.")
    return(NULL)
  }
  map <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE),
    error = function(e) { warning("Could not read EPICv2->EPIC map: ",
                                  conditionMessage(e)); NULL }
  )
  if (is.null(map)) return(NULL)
  unique(map$ID_EPIC1[map$big_delta == "FALSE"])
}

#' Convert a collapsed EPICv2 beta vector into EPIC (v1) probe space
#'
#' In the Zhou Lab EPICv2->EPIC map the EPICv2 base probe ID (suffix stripped)
#' always equals the EPIC1 ID, and our EPICv2 betas are already collapsed to
#' base \code{cg} IDs — so conversion is a reliability-filtered probe subset:
#' keep probes whose EPIC1 target is well-calibrated, drop the rest.
#'
#' @param betas Named numeric vector of collapsed EPICv2 betas (base cg IDs)
#' @param reliable Reliable EPIC probe IDs from load_epicv2_reliable_epic_probes()
#' @return Beta vector restricted to reliable EPIC probes, or NULL if none match
convert_epicv2_to_epic <- function(betas, reliable) {
  if (is.null(reliable)) return(NULL)
  keep <- intersect(names(betas), reliable)
  if (!length(keep)) return(NULL)
  betas[keep]
}

#' Base probe ID — strip the EPICv2 replicate suffix (cg00000029_TC21 -> cg00000029)
#' so probe sets match whether betas are collapsed or suffixed.
base_probe_id <- function(x) sub("_.*$", "", x)

#' X-linked (X-inactivation-informative) probe set for the platform, base cg IDs.
#' EPICv2 base IDs equal EPIC1 IDs, so EPIC's curated set applies to EPICv2 too.
load_xlinked_probes <- function(array_type) {
  ds <- switch(toupper(as.character(array_type)),
               "450K" = "HM450.probeInfo", "EPIC" = "EPIC.probeInfo",
               "EPICV2" = "EPIC.probeInfo", NULL)
  if (is.null(ds)) return(NULL)
  tryCatch(sesameDataGet(ds)$chrX.xlinked, error = function(e) {
    warning("X-linked probe set unavailable for karyotype: ", conditionMessage(e))
    NULL
  })
}

#' Infer a coarse sex karyotype from inferSex + X-inactivation heterozygosity.
#' Two active X chromosomes methylate X-linked CpGs to intermediate beta
#' (X-inactivation), so the fraction of X-linked probes with 0.3<beta<0.7
#' (X_het) cleanly separates two-X (~0.48) from one-X (~0.13) across EPIC/450k.
#' Combined with the MALE/FEMALE call this yields XX / XY and flags the two
#' realistically detectable aneuploidies as uncertain guesses (XXY? = male +
#' two-X; X0? = female + one-X).
#'
#' X-dosage (the call backbone) deliberately uses methylation, not Y-channel
#' intensity, which is noisy on degraded samples. Y intensity is used only for a
#' separate Loss-of-Y flag: when the sample looks male by methylation (one X)
#' but minfi's Y/X intensity gap (\code{y_minus_x}, log2 yMed-xMed) drops below
#' \code{loy_cutoff} (the same -2 minfi uses to call female), the Y is depleted
#' — common somatic LOY in tumours, and exactly the case where sesame says MALE
#' while minfi says FEMALE. Flagged as "XY (low Y - possible LOY)". When
#' y_minus_x is NA the flag simply never fires.
#'
#' Returns NA when sex is unknown or too few X-linked probes survive.
infer_karyotype <- function(beta_sample, sex, xlinked,
                            y_minus_x = NA_real_, loy_cutoff = -2) {
  if (is.null(xlinked) || is.na(sex))
    return(list(karyotype = NA_character_, x_het = NA_real_))
  bx <- beta_sample[base_probe_id(names(beta_sample)) %in% xlinked]
  bx <- bx[!is.na(bx)]
  if (length(bx) < 50) return(list(karyotype = NA_character_, x_het = NA_real_))
  x_het <- mean(bx > 0.3 & bx < 0.7)
  two_x <- x_het > 0.32
  one_x <- x_het < 0.25
  y_lost <- !is.na(y_minus_x) && y_minus_x < loy_cutoff
  kary <- if (sex == "FEMALE" && two_x) "XX"
          else if (sex == "MALE" && one_x && y_lost) "XY (low Y - possible LOY)"
          else if (sex == "MALE"   && one_x) "XY"
          else if (sex == "MALE"   && two_x) "XXY? (uncertain - verify)"
          else if (sex == "FEMALE" && one_x) "X0? (uncertain - verify)"
          else if (sex == "FEMALE") "XX"            # ambiguous X_het: trust sex
          else "XY"
  list(karyotype = kary, x_het = round(x_het, 3))
}

#' rs-SNP genotype fingerprint for sample-swap / identity matching. The ~59-65
#' Infinium rs probes genotype germline SNPs (beta ~0/0.5/1 -> AA/AB/BB), a
#' per-individual barcode: two arrays from the same person share ~all calls.
#' Replaces the removed sesame inferEthnicity (its model was dropped upstream)
#' with the sample-integrity signal actually wanted. Ordered by probe ID for
#' cross-sample comparability. Returns the genotype string and usable-SNP count.
#' Genotypes are encoded as letters A/H/B (AA / AB heterozygous / BB), NOT
#' digits: an all-digit string round-trips through read.csv as a number (a 60+
#' digit value collapses to scientific-notation garbage); letters stay character.
snp_fingerprint <- function(beta_sample) {
  rs <- beta_sample[grepl("^rs", base_probe_id(names(beta_sample)))]
  if (!length(rs)) return(list(fingerprint = NA_character_, n_snp = 0L))
  # Order by probe ID over the FULL rs set and keep NA probes as "." so position
  # i is the same SNP in every sample on a platform — required for the strings
  # to be comparable across samples (the sample-matching use case).
  rs <- rs[order(base_probe_id(names(rs)))]
  geno <- ifelse(is.na(rs), ".",
                 ifelse(rs < 1 / 3, "A", ifelse(rs > 2 / 3, "B", "H")))
  list(fingerprint = paste(geno, collapse = ""), n_snp = sum(!is.na(rs)))
}

#' Genotype concordance (%) between two position-aligned fingerprint strings,
#' over the SNPs both samples genotyped (ignoring "." positions). NA if fewer
#' than min_snp SNPs are jointly usable.
snp_concordance <- function(a, b, min_snp = 10) {
  if (is.na(a) || is.na(b)) return(NA_real_)
  ca <- strsplit(a, "")[[1]]; cb <- strsplit(b, "")[[1]]
  if (length(ca) != length(cb)) return(NA_real_)   # different platform/probe set
  ok <- ca != "." & cb != "."
  if (sum(ok) < min_snp) return(NA_real_)
  100 * mean(ca[ok] == cb[ok])
}

#' For each fingerprint, find its closest other sample in the batch by genotype
#' concordance. Turns the raw barcode into two interpretable numbers: the
#' best-matching Sample_ID and the % concordance. Same individual -> ~95-100%;
#' unrelated -> ~35-50%. Returns a data.frame(SNP_BestMatch, SNP_Match_Pct).
#' With <2 samples (nothing to compare against) both columns are NA.
snp_best_matches <- function(fingerprints, sample_ids, min_snp = 10) {
  n <- length(fingerprints)
  best_id  <- rep(NA_character_, n)
  best_pct <- rep(NA_real_, n)
  if (n >= 2) {
    for (i in seq_len(n)) {
      pcts <- vapply(seq_len(n), function(j) {
        if (j == i) return(NA_real_)
        snp_concordance(fingerprints[i], fingerprints[j], min_snp)
      }, numeric(1))
      if (any(!is.na(pcts))) {
        j <- which.max(pcts)
        best_id[i]  <- sample_ids[j]
        best_pct[i] <- round(pcts[j], 1)
      }
    }
  }
  data.frame(SNP_BestMatch = best_id, SNP_Match_Pct = best_pct,
             stringsAsFactors = FALSE)
}

#' Full pairwise genotype-concordance matrix (%) across the batch, for the
#' identity heatmap. Symmetric, diagonal = 100 (self), NA where two samples
#' share too few usable SNPs. Rows/cols named by Sample_ID.
snp_concordance_matrix <- function(fingerprints, sample_ids, min_snp = 10) {
  n <- length(sample_ids)
  m <- matrix(NA_real_, n, n, dimnames = list(sample_ids, sample_ids))
  for (i in seq_len(n)) {
    m[i, i] <- if (!is.na(fingerprints[i])) 100 else NA_real_
    if (i < n) for (j in (i + 1):n) {
      v <- snp_concordance(fingerprints[i], fingerprints[j], min_snp)
      m[i, j] <- v
      m[j, i] <- v
    }
  }
  round(m, 1)
}

#' Per-sample sesame sample-integrity inferences (sex, karyotype, age,
#' leukocyte fraction, SNP fingerprint)
#'
#' Valuable for detecting sample swaps / mislabelling and gauging tumour
#' purity. Sex uses sesame's curated X/Y probe model; karyotype refines it with
#' X-inactivation heterozygosity; age uses the Horvath353 clock via
#' sesame::predictAge(); leukocyte fraction uses sesame's two-component model
#' (EPIC/HM450, EPICv2 via the EPIC map); the SNP fingerprint barcodes the rs
#' probes. All operate on the beta matrix already computed upstream, each
#' wrapped so a failure yields NA rather than an error. All are informational
#' and never gate Pass_QC.
#'
#' @param beta_values Beta matrix (probes x samples)
#' @param array_type Array type ("450k","EPIC","EPICv2"); drives leukocyte platform
#' @param y_minus_x Optional named numeric vector (by Sample_ID) of minfi's
#'   log2 Y-minus-X median intensity (yMed - xMed) — enables the Loss-of-Y flag
#'   in the karyotype. NA/absent entries simply skip the flag.
#' @return data.frame(Sample_ID, Sesame_Sex, Karyotype, X_Het, Horvath_Age,
#'   Leukocyte_Fraction, SNP_Fingerprint, SNP_Count, SNP_BestMatch, SNP_Match_Pct).
#'   SNP_Fingerprint is computed for the concordance matrix but not surfaced in
#'   the report.
compute_sample_inferences <- function(beta_values, array_type = NULL,
                                      y_minus_x = NULL) {
  sample_ids <- colnames(beta_values)
  out <- data.frame(Sample_ID = sample_ids,
                    Sesame_Sex = NA_character_,
                    Horvath_Age = NA_real_,
                    Leukocyte_Fraction = NA_real_,
                    stringsAsFactors = FALSE)

  # Sex — sesame::inferSex auto-detects platform from probe names.
  out$Sesame_Sex <- vapply(sample_ids, function(sid) {
    tryCatch(as.character(sesame::inferSex(beta_values[, sid])),
             error = function(e) NA_character_)
  }, character(1))

  # Karyotype — inferSex + X-inactivation heterozygosity (informational).
  # y_minus_x (minfi yMed-xMed) adds the Loss-of-Y flag when available.
  xlinked <- load_xlinked_probes(array_type)
  yx <- function(sid) if (!is.null(y_minus_x) && sid %in% names(y_minus_x))
    as.numeric(y_minus_x[[sid]]) else NA_real_
  kary <- lapply(seq_along(sample_ids), function(i)
    infer_karyotype(beta_values[, sample_ids[i]], out$Sesame_Sex[i], xlinked,
                    y_minus_x = yx(sample_ids[i])))
  out$Karyotype <- vapply(kary, function(z) z$karyotype, character(1))
  out$X_Het     <- vapply(kary, function(z) z$x_het, numeric(1))

  # rs-SNP genotype fingerprint (sample-swap / identity), informational.
  # The raw barcode is kept for cross-run archival; the per-sample best-match
  # (closest other sample + % concordance) is the interpretable numeric summary.
  fp <- lapply(sample_ids, function(sid) snp_fingerprint(beta_values[, sid]))
  out$SNP_Fingerprint <- vapply(fp, function(z) z$fingerprint, character(1))
  out$SNP_Count       <- vapply(fp, function(z) z$n_snp, integer(1))
  bm <- snp_best_matches(out$SNP_Fingerprint, sample_ids)
  out$SNP_BestMatch  <- bm$SNP_BestMatch
  out$SNP_Match_Pct  <- bm$SNP_Match_Pct

  # Age — Horvath353 via the proper sesame::predictAge() on the vendored model.
  model <- load_horvath_model()
  if (!is.null(model)) {
    out$Horvath_Age <- round(vapply(sample_ids, function(sid) {
      tryCatch(as.numeric(sesame::predictAge(beta_values[, sid], model)),
               error = function(e) NA_real_)
    }, numeric(1)), 1)
  }

  # Leukocyte fraction — sesame::estimateLeukocyte (reference from sesameData).
  # estimateLeukocyte supports EPIC/HM450/HM27 only. EPICv2 has no native
  # leukocyte reference, so we first convert EPICv2 betas into EPIC space via
  # the Zhou Lab EPICv2->EPIC map, then estimate on the EPIC platform.
  at <- toupper(as.character(array_type))
  leuko_platform <- switch(at, "450K" = "HM450", "EPIC" = "EPIC",
                           "EPICV2" = "EPIC", NULL)
  convert_v2 <- identical(at, "EPICV2")
  # Load the EPICv2->EPIC reliable-probe set once (not per sample).
  reliable <- if (convert_v2) load_epicv2_reliable_epic_probes() else NULL
  if (!is.null(leuko_platform) && !(convert_v2 && is.null(reliable))) {
    out$Leukocyte_Fraction <- round(vapply(sample_ids, function(sid) {
      tryCatch({
        b <- beta_values[, sid]
        if (convert_v2) {
          b <- convert_epicv2_to_epic(b, reliable)
          if (is.null(b)) return(NA_real_)
        }
        as.numeric(sesame::estimateLeukocyte(b, platform = leuko_platform))
      }, error = function(e) NA_real_)
    }, numeric(1)), 4)
  }
  out
}

#' Perform quality control on methylation data
#'
#' @param rgset RGChannelSet object
#' @param beta_values Beta values matrix
#' @param sample_info Sample information data frame
#' @param detection_p_threshold Threshold for per-probe detection p-values
#' @param sample_detection_p_threshold Threshold for mean sample detection p-values
#' @param failed_probe_percent_threshold Max allowed percent of failed probes per sample
#' @param min_median_intensity Minimum acceptable median intensity (log2) for bisulfite check
#' @param gct Optional data frame of GCT bisulfite-conversion scores
#'   (columns Sample_ID, GCT_Score) from the preprocess step. When supplied,
#'   samples whose GCT exceeds \code{max_gct_score} fail QC. Samples with NA
#'   GCT (e.g. EPICv2, where GCT is not yet computed) are never failed on it.
#' @param max_gct_score GCT failure threshold. A score near 1.0 means complete
#'   bisulfite conversion; higher means more incomplete. Samples with
#'   GCT > max_gct_score fail QC.
#' @param array_type Array type ("450k","EPIC","EPICv2"); used to pick the
#'   leukocyte-fraction reference platform. NULL skips leukocyte estimation.
#' @param output_dir Output directory for QC data/report (CSV, RData)
#' @param plots_dir  Output directory for QC plots (PDF/HTML). Defaults to
#'                   \code{file.path(output_dir, "plots")} for backward compat.
#' @return List of QC results and plots
perform_qc <- function(rgset, beta_values, sample_info,
                       detection_p_threshold          = 0.01,
                       sample_detection_p_threshold   = 0.05,
                       failed_probe_percent_threshold = 25,
                       min_median_intensity           = 10.5,
                       gct                            = NULL,
                       max_gct_score                  = 1.3,
                       array_type                     = NULL,
                       output_dir = ".",
                       plots_dir  = NULL) {

  if (is.null(plots_dir)) plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ---------------------------------------------------------------------------
  # Detection p-values
  # ---------------------------------------------------------------------------
  message("Calculating detection p-values...")
  detP <- detectionP(rgset)

  mean_detP            <- colMeans(detP)
  failed_probes_count  <- colSums(detP > detection_p_threshold)
  failed_probes_percent <- failed_probes_count / nrow(detP) * 100

  # ---------------------------------------------------------------------------
  # Median intensity via minfi QC (raw, pre-normalisation)
  # ---------------------------------------------------------------------------
  message("Calculating median channel intensities (pre-normalisation)...")
  ms_raw        <- minfi::preprocessRaw(rgset)
  mfi_qc        <- minfi::getQC(ms_raw)
  median_meth   <- setNames(mfi_qc$mMed, rownames(mfi_qc))
  median_unmeth <- setNames(mfi_qc$uMed, rownames(mfi_qc))

  # ---------------------------------------------------------------------------
  # Build per-sample QC table
  # ---------------------------------------------------------------------------
  sample_ids <- colnames(detP)

  sample_qc <- data.frame(
    Sample_ID                    = sample_ids,
    Mean_Detection_P             = mean_detP[sample_ids],
    Failed_Probes_Count          = failed_probes_count[sample_ids],
    Failed_Probes_Percent        = round(failed_probes_percent[sample_ids], 3),
    Median_Meth_Intensity        = round(median_meth[sample_ids],   2),
    Median_Unmeth_Intensity      = round(median_unmeth[sample_ids], 2),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # ---------------------------------------------------------------------------
  # Sample-integrity inferences (sesame): predicted sex + Horvath epigenetic
  # age. Informational only — surfaced for sample-swap / mislabelling checks,
  # they do NOT contribute to Pass_QC. Aligned to sample_qc by Sample_ID.
  # ---------------------------------------------------------------------------
  message("Inferring sample sex, karyotype, epigenetic age, leukocyte fraction, and SNP fingerprint (sesame)...")
  # Minfi yMed-xMed per sample (from preprocess), keyed by Sample_ID, for the
  # karyotype Loss-of-Y flag. Absent on older preprocessed data -> NULL -> no flag.
  y_minus_x <- NULL
  if (!is.null(sample_info) &&
      all(c("Sample_ID", "Minfi_xMed", "Minfi_yMed") %in% names(sample_info))) {
    y_minus_x <- setNames(sample_info$Minfi_yMed - sample_info$Minfi_xMed,
                          as.character(sample_info$Sample_ID))
  }
  inferences <- tryCatch(
    compute_sample_inferences(beta_values, array_type = array_type,
                              y_minus_x = y_minus_x),
    error = function(e) {
      warning("Sample inference (sex/karyotype/age/leukocyte/SNP) failed: ",
              conditionMessage(e))
      NULL
    }
  )
  sample_qc$Sesame_Sex         <- NA_character_
  sample_qc$Karyotype          <- NA_character_
  sample_qc$X_Het              <- NA_real_
  sample_qc$Horvath_Age        <- NA_real_
  sample_qc$Leukocyte_Fraction <- NA_real_
  sample_qc$SNP_BestMatch      <- NA_character_
  sample_qc$SNP_Match_Pct      <- NA_real_
  sample_qc$SNP_Count          <- NA_integer_
  if (!is.null(inferences)) {
    lk <- function(col, cast) {
      v <- setNames(inferences[[col]], inferences$Sample_ID)[sample_qc$Sample_ID]
      cast(v)
    }
    sample_qc$Sesame_Sex         <- lk("Sesame_Sex",         as.character)
    sample_qc$Karyotype          <- lk("Karyotype",          as.character)
    sample_qc$X_Het              <- lk("X_Het",              as.numeric)
    sample_qc$Horvath_Age        <- lk("Horvath_Age",        as.numeric)
    sample_qc$Leukocyte_Fraction <- lk("Leukocyte_Fraction", as.numeric)
    sample_qc$SNP_BestMatch      <- lk("SNP_BestMatch",      as.character)
    sample_qc$SNP_Match_Pct      <- lk("SNP_Match_Pct",      as.numeric)
    sample_qc$SNP_Count          <- lk("SNP_Count",          as.integer)
    # The raw per-sample fingerprint is no longer surfaced in the report — the
    # pairwise concordance matrix below + SNP_BestMatch/Match_Pct replace it.
    # Write the full matrix for the in-app identity heatmap (>= 2 samples).
    if (nrow(inferences) >= 2 && "SNP_Fingerprint" %in% names(inferences)) {
      cm <- snp_concordance_matrix(inferences$SNP_Fingerprint, inferences$Sample_ID)
      utils::write.csv(cm, file.path(output_dir, "snp_concordance.csv"))
    }
  }

  # ---------------------------------------------------------------------------
  # Failure flags — NOTE: low intensity is INFORMATIONAL ONLY, not a failure.
  # Scanner gain settings legitimately vary across sites/batches; SWAN
  # normalisation may recover affected samples.
  # ---------------------------------------------------------------------------
  sample_qc$Flag_Mean_DetP     <- sample_qc$Mean_Detection_P      >= sample_detection_p_threshold
  sample_qc$Flag_Failed_Probes <- sample_qc$Failed_Probes_Percent >= failed_probe_percent_threshold

  # Informational intensity flag — does NOT contribute to Pass_QC
  sample_qc$Note_Low_Intensity <- sample_qc$Median_Meth_Intensity   < min_median_intensity |
                                  sample_qc$Median_Unmeth_Intensity  < min_median_intensity

  # ---------------------------------------------------------------------------
  # GCT bisulfite-conversion gate. Merge the per-sample GCT score (from the
  # preprocess step) and fail samples whose conversion is too incomplete.
  # NA GCT (e.g. EPICv2, not yet computed) is never a failure — Flag_GCT FALSE.
  # ---------------------------------------------------------------------------
  sample_qc$GCT_Score <- NA_real_
  if (!is.null(gct) && all(c("Sample_ID", "GCT_Score") %in% names(gct))) {
    gct_lookup <- setNames(gct$GCT_Score, as.character(gct$Sample_ID))
    sample_qc$GCT_Score <- as.numeric(gct_lookup[sample_qc$Sample_ID])
  }
  sample_qc$Flag_GCT <- !is.na(sample_qc$GCT_Score) &
                        sample_qc$GCT_Score > max_gct_score

  # Pass/fail based on detection p, failed probe rate, and GCT conversion.
  sample_qc$Pass_QC <- !(sample_qc$Flag_Mean_DetP |
                         sample_qc$Flag_Failed_Probes |
                         sample_qc$Flag_GCT)

  # ---------------------------------------------------------------------------
  # SWAN recovery check for low-intensity samples
  # Run SWAN normalisation and check whether median intensities improve above
  # the threshold — helps distinguish scanner gain artefacts from true failures.
  # ---------------------------------------------------------------------------
  low_int_ids <- sample_qc$Sample_ID[sample_qc$Note_Low_Intensity]

  sample_qc$SWAN_Median_Meth   <- NA_real_
  sample_qc$SWAN_Median_Unmeth <- NA_real_
  sample_qc$SWAN_Recoverable   <- NA

  if (length(low_int_ids) > 0) {
    message(sprintf(
      "%d sample(s) have low median intensity (< %.1f) — running SWAN normalisation check ...",
      length(low_int_ids), min_median_intensity
    ))
    tryCatch({
      ms_swan    <- minfi::preprocessSWAN(rgset, mSet = ms_raw, verbose = FALSE)
      qc_swan    <- minfi::getQC(ms_swan)
      swan_meth  <- setNames(qc_swan$mMed, rownames(qc_swan))
      swan_unmeth <- setNames(qc_swan$uMed, rownames(qc_swan))

      for (sid in low_int_ids) {
        if (sid %in% names(swan_meth)) {
          sm <- round(swan_meth[sid],   2)
          su <- round(swan_unmeth[sid], 2)
          sample_qc$SWAN_Median_Meth[sample_qc$Sample_ID   == sid] <- sm
          sample_qc$SWAN_Median_Unmeth[sample_qc$Sample_ID == sid] <- su
          # Recoverable = both channels exceed threshold after SWAN
          sample_qc$SWAN_Recoverable[sample_qc$Sample_ID   == sid] <-
            sm >= min_median_intensity & su >= min_median_intensity
        }
      }

      n_recovered <- sum(sample_qc$SWAN_Recoverable %in% TRUE)
      n_not       <- sum(sample_qc$SWAN_Recoverable %in% FALSE)
      message(sprintf(
        "  SWAN recovery: %d recoverable, %d not recoverable, %d not checked",
        n_recovered, n_not,
        sum(is.na(sample_qc$SWAN_Recoverable))
      ))
      if (n_not > 0) {
        not_rec <- sample_qc$Sample_ID[sample_qc$SWAN_Recoverable %in% FALSE]
        message("  Not recoverable after SWAN: ", paste(not_rec, collapse = ", "))
      }
    }, error = function(e) {
      warning("SWAN normalisation check failed: ", e$message)
    })
  }

  # ---------------------------------------------------------------------------
  # Human-readable failure reason (detection p / probe rate only)
  # and informational notes (intensity + SWAN result)
  # ---------------------------------------------------------------------------
  sample_qc$Failure_Reason <- apply(sample_qc, 1, function(r) {
    reasons <- character(0)
    if (isTRUE(as.logical(r["Flag_Mean_DetP"])))
      reasons <- c(reasons, sprintf("Mean detection p (%.4f) >= %.4f",
                                    as.numeric(r["Mean_Detection_P"]),
                                    sample_detection_p_threshold))
    if (isTRUE(as.logical(r["Flag_Failed_Probes"])))
      reasons <- c(reasons, sprintf("Failed probes (%.2f%%) >= %.1f%%",
                                    as.numeric(r["Failed_Probes_Percent"]),
                                    failed_probe_percent_threshold))
    if (isTRUE(as.logical(r["Flag_GCT"])))
      reasons <- c(reasons, sprintf("Incomplete bisulfite conversion (GCT %.3f > %.2f)",
                                    as.numeric(r["GCT_Score"]),
                                    max_gct_score))
    if (length(reasons) == 0) "PASS" else paste(reasons, collapse = "; ")
  })

  sample_qc$Notes <- apply(sample_qc, 1, function(r) {
    notes <- character(0)
    if (isTRUE(as.logical(r["Note_Low_Intensity"]))) {
      swan_rec <- r["SWAN_Recoverable"]
      intensity_note <- sprintf(
        "Low pre-norm intensity (Meth=%.1f, Unmeth=%.1f; threshold=%.1f)",
        as.numeric(r["Median_Meth_Intensity"]),
        as.numeric(r["Median_Unmeth_Intensity"]),
        min_median_intensity
      )
      swan_note <- if (is.na(swan_rec)) {
        "SWAN check not run"
      } else if (isTRUE(as.logical(swan_rec))) {
        sprintf("recoverable after SWAN (Meth=%.1f, Unmeth=%.1f)",
                as.numeric(r["SWAN_Median_Meth"]),
                as.numeric(r["SWAN_Median_Unmeth"]))
      } else {
        sprintf("NOT recoverable after SWAN (Meth=%.1f, Unmeth=%.1f)",
                as.numeric(r["SWAN_Median_Meth"]),
                as.numeric(r["SWAN_Median_Unmeth"]))
      }
      notes <- c(notes, paste0(intensity_note, " — ", swan_note))
    }
    if (length(notes) == 0) "" else paste(notes, collapse = "; ")
  })

  # Merge with sample info (keep all QC rows)
  sample_qc <- merge(sample_qc, sample_info, by = "Sample_ID", all.x = TRUE)

  # ---------------------------------------------------------------------------
  # Save CSV report (key QC columns first, then sample metadata columns)
  # ---------------------------------------------------------------------------
  key_cols  <- c("Sample_ID", "Pass_QC", "Failure_Reason", "Notes",
                 "Mean_Detection_P", "Failed_Probes_Count", "Failed_Probes_Percent",
                 "Median_Meth_Intensity", "Median_Unmeth_Intensity",
                 "GCT_Score", "Sesame_Sex", "Karyotype", "X_Het", "Horvath_Age",
                 "Leukocyte_Fraction", "SNP_BestMatch", "SNP_Match_Pct",
                 "SNP_Count",
                 "Flag_Mean_DetP", "Flag_Failed_Probes", "Flag_GCT", "Note_Low_Intensity",
                 "SWAN_Median_Meth", "SWAN_Median_Unmeth", "SWAN_Recoverable")
  extra_cols <- setdiff(names(sample_qc), key_cols)
  col_order  <- c(intersect(key_cols, names(sample_qc)), extra_cols)

  out_path <- file.path(output_dir, "sample_qc_report.csv")
  write.csv(sample_qc[, col_order],
            file      = out_path,
            row.names = FALSE)
  message("QC report saved: ", out_path)

  # ---------------------------------------------------------------------------
  # QC plots
  # ---------------------------------------------------------------------------
  qc_plots <- generate_qc_plots(rgset, detP, beta_values, sample_qc,
                                detection_p_threshold, sample_detection_p_threshold,
                                plots_dir)

  # minfi PDF QC report
  tryCatch(
    minfi::qcReport(rgset,
                    sampNames = if ("Sample_Name" %in% names(sample_qc))
                                  sample_qc$Sample_Name
                                else
                                  sample_qc$Sample_ID,
                    pdf = file.path(plots_dir, "minfi_qcReport.pdf")),
    error = function(e) warning("minfi::qcReport failed: ", e$message)
  )

  # ---------------------------------------------------------------------------
  # Summary message
  # ---------------------------------------------------------------------------
  n_pass     <- sum(sample_qc$Pass_QC)
  n_fail     <- nrow(sample_qc) - n_pass
  n_low_int  <- sum(sample_qc$Note_Low_Intensity)
  message(sprintf("QC complete: %d passed, %d failed (see sample_qc_report.csv)",
                  n_pass, n_fail))
  if (n_fail > 0) {
    message("Failed samples:")
    failed_df <- sample_qc[!sample_qc$Pass_QC, c("Sample_ID", "Failure_Reason")]
    for (i in seq_len(nrow(failed_df)))
      message(sprintf("  %s: %s", failed_df$Sample_ID[i], failed_df$Failure_Reason[i]))
  }
  if (n_low_int > 0) {
    message(sprintf(
      "%d sample(s) flagged for low pre-normalisation intensity (informational, not failed):",
      n_low_int
    ))
    low_df <- sample_qc[sample_qc$Note_Low_Intensity, c("Sample_ID", "Notes")]
    for (i in seq_len(nrow(low_df)))
      message(sprintf("  %s: %s", low_df$Sample_ID[i], low_df$Notes[i]))
  }
  n_gct_fail <- sum(sample_qc$Flag_GCT)
  if (n_gct_fail > 0) {
    message(sprintf(
      "%d sample(s) FAILED for incomplete bisulfite conversion (GCT > %.2f):",
      n_gct_fail, max_gct_score
    ))
    gct_df <- sample_qc[sample_qc$Flag_GCT, c("Sample_ID", "GCT_Score")]
    for (i in seq_len(nrow(gct_df)))
      message(sprintf("  %s: GCT %.3f", gct_df$Sample_ID[i], gct_df$GCT_Score[i]))
  }

  list(
    sample_qc                    = sample_qc,
    passed_samples               = sample_qc$Sample_ID[sample_qc$Pass_QC],
    failed_samples               = sample_qc$Sample_ID[!sample_qc$Pass_QC],
    low_intensity_samples        = sample_qc$Sample_ID[sample_qc$Note_Low_Intensity],
    swan_recoverable             = sample_qc$Sample_ID[sample_qc$SWAN_Recoverable %in% TRUE],
    swan_not_recoverable         = sample_qc$Sample_ID[sample_qc$SWAN_Recoverable %in% FALSE],
    gct_failed_samples           = sample_qc$Sample_ID[sample_qc$Flag_GCT],
    detection_p                  = detP,
    detection_p_threshold        = detection_p_threshold,
    sample_detection_p_threshold = sample_detection_p_threshold,
    max_gct_score                = max_gct_score,
    plots                        = qc_plots
  )
}


#' Generate QC plots
#'
#' @param rgset RGChannelSet object
#' @param detP Detection p-values matrix
#' @param beta_values Beta values matrix
#' @param sample_qc Sample QC metrics data frame
#' @param detection_p_threshold Threshold for detection p-values
#' @param sample_detection_p_threshold Threshold for mean sample detection p-values
#' @param output_dir Output directory for plots
#' @return List of plot paths
generate_qc_plots <- function(rgset, detP, beta_values, sample_qc, 
                              detection_p_threshold = 0.01,
                              sample_detection_p_threshold = 0.05,
                              output_dir = ".") {
  plots <- list()
  
  # 1. Mean detection p-value plot
  message("Generating mean detection p-value plot...")
  pdf_path <- file.path(output_dir, "mean_detection_pvalue.pdf")
  pdf(pdf_path, height = 8, width = 12)
  barplot(sample_qc$Mean_Detection_P, 
          names.arg = sample_qc$Sample_ID, 
          las = 2, 
          cex.names = 0.8,
          col = ifelse(sample_qc$Pass_QC, "forestgreen", "firebrick"),
          main = "Mean Detection P-value by Sample",
          ylab = "Mean Detection P-value")
  abline(h = sample_detection_p_threshold, col = "red", lty = 2)
  text(x = par("usr")[2] * 0.95, 
       y = sample_detection_p_threshold * 1.1, 
       labels = paste("Threshold (", sample_detection_p_threshold, ")", sep = ""), 
       cex = 0.8)
  dev.off()
  plots$mean_detection_pvalue <- pdf_path
  
  # 2. Sample density plot
  message("Generating sample density plot...")
  pdf_path <- file.path(output_dir, "beta_density.pdf")
  pdf(pdf_path, height = 8, width = 10)
  densityPlot(beta_values, 
              sampGroups = sample_qc$Pass_QC, 
              main = "Beta Value Density Plot",
              legend = FALSE)
  legend("topright", 
         legend = c("Pass QC", "Fail QC"), 
         col = c("black", "red"), 
         lty = 1)
  dev.off()
  plots$beta_density <- pdf_path
  
  # 3. Bean plot for beta distribution
  message("Generating bean plot...")
  pdf_path <- file.path(output_dir, "beta_bean_plot.pdf")
  pdf(pdf_path, height = 8, width = 10)
  densityBeanPlot(beta_values, 
                  sampGroups = sample_qc$Pass_QC,
                  main = "Beta Value Distribution")
  dev.off()
  plots$beta_bean <- pdf_path
  
 
  # 4. MDS plot if we have more than 3 samples
  if (ncol(beta_values) > 3) {
    message("Generating MDS plot...")
    pdf_path <- file.path(output_dir, "mds_plot.pdf")
    pdf(pdf_path, height = 8, width = 10)
    
    # Use tryCatch in case MDS calculation fails
    tryCatch({
      mds <- cmdscale(dist(t(beta_values)), k = 3)
      colnames(mds) <- c("PC1", "PC2", "PC3")
      par(mfrow = c(2, 1))
      plot(mds[, 1], mds[, 2], 
           col = ifelse(sample_qc$Pass_QC, "blue", "red"),
           pch = 19,
           main = "MDS Plot - PC1 vs PC2",
           xlab = "PC1", ylab = "PC2")
      text(mds[, 1], mds[, 2], labels = sample_qc$Sample_ID, pos = 3, cex = 0.8)
      
      plot(mds[, 1], mds[, 3], 
           col = ifelse(sample_qc$Pass_QC, "blue", "red"),
           pch = 19,
           main = "MDS Plot - PC1 vs PC3",
           xlab = "PC1", ylab = "PC3")
      text(mds[, 1], mds[, 3], labels = sample_qc$Sample_ID, pos = 3, cex = 0.8)
    }, error = function(e) {
      plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)
      text(1, 1, "MDS plot could not be generated.\nError: ")
      text(1, 0.8, e$message, col = "red")
      warning("Could not generate MDS plot: ", e$message)
    })
    
    dev.off()
    plots$mds <- pdf_path
  }
  
  # Create interactive plots using plotly if requested
  if (requireNamespace("plotly", quietly = TRUE) && 
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    
    # Try to create interactive plots with error handling
    tryCatch({
      # Interactive density plot
      message("Generating interactive density plot...")
      density_data <- lapply(1:ncol(beta_values), function(i) {
        dens <- density(beta_values[, i], na.rm = TRUE)
        data.frame(
          x = dens$x,
          y = dens$y,
          Sample = colnames(beta_values)[i],
          Pass_QC = sample_qc$Pass_QC[match(colnames(beta_values)[i], sample_qc$Sample_ID)]
        )
      })
      density_data <- do.call(rbind, density_data)
      
      p <- plotly::plot_ly()
      for (sample in unique(density_data$Sample)) {
        subset_data <- density_data[density_data$Sample == sample, ]
        pass_qc <- subset_data$Pass_QC[1]
        # Theme colors: teal for QC pass, red for QC fail.
        p <- p %>% plotly::add_lines(
          data = subset_data,
          x = ~x,
          y = ~y,
          name = sample,
          opacity = 0.9,
          line = list(color = ifelse(pass_qc, "#0d9488", "#ef4444"))
        )
      }
      p <- p %>% plotly::layout(
        title = "Beta Value Density Distribution",
        xaxis = list(title = "Beta Value"),
        yaxis = list(title = "Density")
      )
      
      html_path <- file.path(output_dir, "interactive_density_plot.html")
      htmlwidgets::saveWidget(plotly::as_widget(p), html_path)
      plots$interactive_density <- html_path
      
      # Interactive MDS plot if we have more than 3 samples
      if (ncol(beta_values) > 3 && exists("mds")) {
        message("Generating interactive MDS plot...")
        mds_data <- as.data.frame(mds)
        mds_data$Sample_ID <- rownames(mds_data)
        mds_data <- merge(mds_data, sample_qc, by = "Sample_ID")
        
        # Uniform marker styling: teal fill (#018571), brown border
        # (#a6611a). Pass/fail signaling is surfaced on the QC tab's
        # table; keeping the 3D MDS colors uniform makes the point cloud
        # read cleaner.
        p <- plotly::plot_ly(
          data = mds_data,
          x = ~PC1,
          y = ~PC2,
          z = ~PC3,
          text = ~Sample_ID,
          type = "scatter3d",
          mode = "markers",
          marker = list(
            size = 5,
            opacity = 0.9,
            color = "#018571",
            line = list(width = 0.8, color = "#a6611a")
          )
        ) %>%
          plotly::layout(
            title = "3D MDS Plot",
            scene = list(
              xaxis = list(title = "PC1"),
              yaxis = list(title = "PC2"),
              zaxis = list(title = "PC3")
            )
          )
        
        html_path <- file.path(output_dir, "interactive_mds_plot.html")
        htmlwidgets::saveWidget(plotly::as_widget(p), html_path)
        plots$interactive_mds <- html_path
      }
    }, error = function(e) {
      warning("Could not generate interactive plots: ", e$message)
    })
  }
  
  return(plots)
}