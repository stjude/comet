# dim_reduction.R
# Dimensionality reduction module: t-SNE, UMAP, and hierarchical clustering.
# Functions exported:
#   select_variable_probes(beta, n_probes, method)
#   run_tsne(beta, sample_info, perplexity, dimensions, output_dir)
#   run_umap(beta, sample_info, n_neighbors, output_dir)
#   run_hierarchical_clustering(beta, sample_info, method, distance, output_dir)

# ---------------------------------------------------------------------------
# Select top variable probes
# ---------------------------------------------------------------------------

#' Select the most variable probes from a beta matrix.
#'
#' @param beta        Numeric matrix, probes x samples.
#' @param n_probes    Number of top variable probes to retain.
#' @param method      Variance metric: "sd" (default) or "mad".
#' @return            Subset beta matrix (probes x samples).
select_variable_probes <- function(beta, n_probes = 10000, method = "sd") {
  stopifnot(is.matrix(beta) || is.data.frame(beta))
  beta <- as.matrix(beta)

  variability <- switch(method,
    sd  = apply(beta, 1, sd,  na.rm = TRUE),
    mad = apply(beta, 1, mad, na.rm = TRUE),
    stop("Unknown method '", method, "'. Use 'sd' or 'mad'.")
  )

  n_probes <- min(n_probes, nrow(beta))
  top_idx  <- order(variability, decreasing = TRUE)[seq_len(n_probes)]
  beta[top_idx, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# t-SNE
# ---------------------------------------------------------------------------

#' Run t-SNE on a beta matrix and save results.
#'
#' Duplicate samples (identical rows after transposition) are detected and
#' removed before passing data to Rtsne, then re-attached to the output.
#'
#' @param beta        Numeric matrix, probes x samples.
#' @param sample_info Data frame with sample metadata (rownames = sample IDs).
#' @param perplexity  t-SNE perplexity (capped at (n_samples - 1) / 3).
#' @param dimensions  Number of t-SNE dimensions (default 2).
#' @param output_dir  Directory for saving RData / coordinate files.
#' @param plots_dir   Directory for saving plot files. Defaults to \code{output_dir}.
#' @return            List with elements: coords, sample_info, duplicates.
run_tsne <- function(beta,
                     sample_info,
                     perplexity  = 5,
                     dimensions  = 2,
                     output_dir  = ".",
                     plots_dir   = NULL) {
  # Perplexity default of 5 matches default_config()$dim_reduction$tsne$perplexity
  # and is safe for small cohorts (Rtsne requires perplexity < (N-1)/3).
  # Users with larger N should override via config or the Settings UI.

  suppressPackageStartupMessages(library(Rtsne))
  suppressPackageStartupMessages(library(ggplot2))

  if (is.null(plots_dir)) plots_dir <- output_dir
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plots_dir))  dir.create(plots_dir,  recursive = TRUE)

  beta <- as.matrix(beta)
  # t-SNE operates on samples x probes
  t_beta <- t(beta)

  # ------------------------------------------------------------------
  # Duplicate removal
  # ------------------------------------------------------------------
  dup_mask   <- duplicated(t_beta)
  n_dups     <- sum(dup_mask)
  dup_names  <- rownames(t_beta)[dup_mask]

  if (n_dups > 0) {
    warning(sprintf(
      "run_tsne: %d duplicate sample(s) removed before t-SNE: %s",
      n_dups, paste(dup_names, collapse = ", ")
    ))
  }

  t_beta_uniq   <- t_beta[!dup_mask, , drop = FALSE]
  sample_info_u <- if (!is.null(sample_info)) {
    sample_info[rownames(t_beta_uniq), , drop = FALSE]
  } else {
    NULL
  }

  # ------------------------------------------------------------------
  # NA handling: keep only complete probes (no NA in any sample)
  # ------------------------------------------------------------------
  if (any(is.na(t_beta_uniq))) {
    complete_probes <- colSums(is.na(t_beta_uniq)) == 0
    n_dropped <- sum(!complete_probes)
    message(sprintf("run_tsne: dropping %d probes with missing values; %d complete probes retained.",
                    n_dropped, sum(complete_probes)))
    t_beta_uniq <- t_beta_uniq[, complete_probes, drop = FALSE]
  }

  n_samples  <- nrow(t_beta_uniq)
  perplexity <- min(perplexity, floor((n_samples - 1) / 3))

  if (perplexity < 1) {
    warning("Too few unique samples for t-SNE (need >= 4). Skipping.")
    return(list(coords = NULL, sample_info = sample_info_u,
                duplicates = dup_names))
  }

  message(sprintf("Running t-SNE: %d samples, %d probes, perplexity = %.0f",
                  n_samples, ncol(t_beta_uniq), perplexity))

  set.seed(42)
  tsne_out <- Rtsne::Rtsne(
    t_beta_uniq,
    dims             = dimensions,
    perplexity       = perplexity,
    check_duplicates = FALSE,   # duplicates already removed above
    verbose          = FALSE,
    max_iter         = 1000
  )

  coords <- as.data.frame(tsne_out$Y)
  colnames(coords) <- paste0("tSNE", seq_len(dimensions))
  rownames(coords) <- rownames(t_beta_uniq)

  if (!is.null(sample_info_u)) {
    coords <- cbind(coords, sample_info_u)
  }

  # Save
  tsne_results <- list(coords = coords, sample_info = sample_info_u,
                       duplicates = dup_names)
  save(tsne_results,
       file = file.path(output_dir, "tsne_results.RData"))

  # Plot (first two dims)
  p <- ggplot(coords, aes(x = tSNE1, y = tSNE2)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_bw(base_size = 12) +
    labs(title = "t-SNE", x = "t-SNE 1", y = "t-SNE 2")

  if (!is.null(sample_info_u) && "Sample_Group" %in% colnames(sample_info_u)) {
    p <- p + aes(colour = Sample_Group)
  }

  ggsave(file.path(plots_dir, "tsne_plot.pdf"), p,
         width = 7, height = 6)

  invisible(tsne_results)
}

# ---------------------------------------------------------------------------
# UMAP
# ---------------------------------------------------------------------------

#' Run UMAP on a beta matrix and save results.
#'
#' @param beta        Numeric matrix, probes x samples.
#' @param sample_info Data frame with sample metadata.
#' @param n_neighbors Number of nearest neighbours (default 15).
#' @param output_dir  Directory for saving RData / coordinate files.
#' @param plots_dir   Directory for saving plot files. Defaults to \code{output_dir}.
#' @return            List with elements: coords, sample_info.
run_umap <- function(beta,
                     sample_info,
                     n_neighbors = 15,
                     output_dir  = ".",
                     plots_dir   = NULL) {

  suppressPackageStartupMessages(library(umap))
  suppressPackageStartupMessages(library(ggplot2))

  if (is.null(plots_dir)) plots_dir <- output_dir
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plots_dir))  dir.create(plots_dir,  recursive = TRUE)

  beta   <- as.matrix(beta)
  t_beta <- t(beta)

  # Remove duplicates (UMAP is more tolerant but consistent with t-SNE)
  dup_mask <- duplicated(t_beta)
  if (any(dup_mask)) {
    warning(sprintf("run_umap: removing %d duplicate sample(s).", sum(dup_mask)))
    t_beta <- t_beta[!dup_mask, , drop = FALSE]
  }

  # NA handling: keep only complete probes (no NA in any sample)
  if (any(is.na(t_beta))) {
    complete_probes <- colSums(is.na(t_beta)) == 0
    n_dropped <- sum(!complete_probes)
    message(sprintf("run_umap: dropping %d probes with missing values; %d complete probes retained.",
                    n_dropped, sum(complete_probes)))
    t_beta <- t_beta[, complete_probes, drop = FALSE]
  }

  n_neighbors <- min(n_neighbors, nrow(t_beta) - 1L)

  umap_config             <- umap::umap.defaults
  umap_config$n_neighbors <- n_neighbors
  umap_config$random_state <- 42L

  message(sprintf("Running UMAP: %d samples, %d probes, n_neighbors = %d",
                  nrow(t_beta), ncol(t_beta), n_neighbors))

  umap_out <- umap::umap(t_beta, config = umap_config)

  coords <- as.data.frame(umap_out$layout)
  colnames(coords) <- c("UMAP1", "UMAP2")
  rownames(coords) <- rownames(t_beta)

  sample_info_u <- if (!is.null(sample_info)) {
    sample_info[rownames(t_beta), , drop = FALSE]
  } else {
    NULL
  }

  if (!is.null(sample_info_u)) {
    coords <- cbind(coords, sample_info_u)
  }

  umap_results <- list(coords = coords, sample_info = sample_info_u,
                       umap_object = umap_out)
  save(umap_results,
       file = file.path(output_dir, "umap_results.RData"))

  p <- ggplot(coords, aes(x = UMAP1, y = UMAP2)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_bw(base_size = 12) +
    labs(title = "UMAP", x = "UMAP 1", y = "UMAP 2")

  if (!is.null(sample_info_u) && "Sample_Group" %in% colnames(sample_info_u)) {
    p <- p + aes(colour = Sample_Group)
  }

  ggsave(file.path(plots_dir, "umap_plot.pdf"), p,
         width = 7, height = 6)

  invisible(umap_results)
}

# ---------------------------------------------------------------------------
# Hierarchical clustering
# ---------------------------------------------------------------------------

#' Run hierarchical clustering on a beta matrix and save a dendrogram.
#'
#' @param beta        Numeric matrix, probes x samples.
#' @param sample_info Data frame with sample metadata.
#' @param method      Agglomeration method passed to hclust() (default "complete").
#' @param distance    Distance metric: "pearson" (default), "spearman", or any
#'                    method accepted by dist().
#' @param output_dir  Directory for saving RData results.
#' @param plots_dir   Directory for saving plot files. Defaults to \code{output_dir}.
#' @return            hclust object.
run_hierarchical_clustering <- function(beta,
                                        sample_info,
                                        method   = "complete",
                                        distance = "pearson",
                                        output_dir = ".",
                                        plots_dir  = NULL) {

  suppressPackageStartupMessages(library(ggplot2))

  if (is.null(plots_dir)) plots_dir <- output_dir
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plots_dir))  dir.create(plots_dir,  recursive = TRUE)

  beta <- as.matrix(beta)

  # Compute distance matrix
  dist_mat <- if (distance %in% c("pearson", "spearman")) {
    cor_mat  <- cor(beta, method = distance, use = "pairwise.complete.obs")
    # Convert correlation to distance (0 = identical, 2 = opposite)
    as.dist(1 - cor_mat)
  } else {
    dist(t(beta), method = distance)
  }

  message(sprintf("Running hierarchical clustering: method = %s, distance = %s",
                  method, distance))

  hc <- hclust(dist_mat, method = method)

  # Annotate labels with group if available
  if (!is.null(sample_info) && "Sample_Group" %in% colnames(sample_info)) {
    grp <- sample_info[hc$labels, "Sample_Group"]
    hc$labels <- ifelse(is.na(grp), hc$labels,
                        paste0(hc$labels, " (", grp, ")"))
  }

  hclust_results <- list(hclust = hc, distance = distance, method = method)
  save(hclust_results,
       file = file.path(output_dir, "hclust_results.RData"))

  # Plot dendrogram
  pdf(file.path(plots_dir, "hclust_dendrogram.pdf"),
      width = max(10, ncol(beta) * 0.3), height = 8)
  par(mar = c(10, 4, 4, 2))
  plot(hc,
       main  = sprintf("Hierarchical Clustering (%s / %s)", method, distance),
       xlab  = "",
       ylab  = "Distance",
       cex   = min(1, 30 / ncol(beta)))
  dev.off()

  invisible(hclust_results)
}
