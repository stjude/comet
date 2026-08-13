# Configuration module for methylation array analysis pipeline

#' Recursively merge override values into a default config tree.
#'
#' Used by the driver so a user-supplied config file only needs to specify
#' the keys it wants to change — every other key falls back to default_config().
#' Lists are merged recursively (so e.g. config$dim_reduction$tsne$perplexity
#' can be overridden without losing the rest of dim_reduction); leaves are
#' replaced wholesale.
deep_merge <- function(default, override) {
  if (is.null(override)) return(default)
  for (key in names(override)) {
    if (is.list(default[[key]]) && is.list(override[[key]]) &&
        !is.null(names(override[[key]]))) {
      default[[key]] <- deep_merge(default[[key]], override[[key]])
    } else {
      default[[key]] <- override[[key]]
    }
  }
  default
}

#' Load configuration from file
#'
#' @param config_file Path to configuration file
#' @return Configuration list
load_config <- function(config_file) {
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }
  
  if (grepl("\\.R$", config_file)) {
    # Source R file
    source(config_file, local = TRUE)
    if (!exists("config")) {
      stop("Configuration file must define a 'config' variable")
    }
    return(config)
  } else if (grepl("\\.yaml$|\\.yml$", config_file)) {
    # Parse YAML file
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' is required for parsing YAML config files")
    }
    return(yaml::read_yaml(config_file))
  } else {
    stop("Unsupported configuration file format")
  }
}

#' Default configuration
#'
#' @return Default configuration list
default_config <- function() {
  # Structure matches the nested keys that run_pipeline() actually reads:
  #   config$preprocessing$normalization
  #   config$qc$<...>
  #   config$filtering$<...>
  #   config$dim_reduction$variable_probes / $tsne$... / $umap$... / $clustering$...
  #   config$cnv$<...>
  # Prior flat-key form was silently ignored — this restructure removes the
  # divergence. Flat-key fallbacks remain in run_pipeline (e.g. for legacy
  # `config$normalization`) so old config files still load.
  list(
    # Input (top-level — read directly on config$)
    sample_sheet      = "sample_sheet.csv",
    reference_samples = NULL,
    data_dir          = "./data",

    # Preprocessing
    preprocessing = list(
      normalization = "swan"   # Options: raw, illumina, functional, quantile, swan
    ),

    # QC
    qc = list(
      detection_p_threshold          = 0.01,
      sample_detection_p_threshold   = 0.05,
      failed_probe_percent_threshold = 25,
      min_median_intensity           = 10.5,
      max_gct_score                  = 1.3,   # GCT bisulfite-conversion fail cutoff (~1.0 = complete)
      filter_failed_samples          = TRUE
    ),

    # Probe filtering
    filtering = list(
      remove_sex_chromosomes  = TRUE,
      remove_snps             = TRUE,
      remove_cross_reactive   = TRUE,
      min_sample_success_rate = 0.75,
      keep_probe_list         = NULL,
      probe_list_column       = "x"
    ),

    # Dimensionality reduction
    dim_reduction = list(
      variable_probes = 10000,
      tsne = list(
        perplexity = 5,    # Safe for small cohorts (N < ~15).
        dimensions = 2
      ),
      umap = list(
        n_neighbors = 15
      ),
      clustering = list(
        method   = "complete",   # single / complete / average / ward.D2
        distance = "pearson"     # pearson / spearman / euclidean / manhattan
      )
    ),

    # CNV analysis. Asymmetric gain/loss thresholds:
    # a segment is called a gain when seg.mean > gain_threshold and a loss
    # when seg.mean < loss_threshold. The legacy single `threshold` key is
    # still honoured (treated as |threshold|, applied symmetrically) for
    # backward compatibility with run_manifest.json files written by older
    # versions of MeQTrack.
    cnv = list(
      method          = "conumee",   # Options: conumee, ChAMP, cnAnalysis450k
      gain_threshold  =  0.18,
      loss_threshold  = -0.20,
      frequency_plot  = TRUE
    ),

    # Reference projection (v2.0.0). Projects the user's samples onto a
    # pre-trained reference t-SNE embedding via snifter. `reference_dir`
    # is resolved relative to the pipeline/ directory when not absolute.
    reference_projection = list(
      enabled       = TRUE,
      dataset       = "COMET_1915",   # key into .REFERENCE_DATASETS
      reference_dir = "../reference",
      perplexity    = 5,
      knn_k         = 25              # neighbours for the nearest-class diagnostic
    ),

    # Output
    output = list(
      include_interactive_plots = TRUE,
      generate_report           = TRUE
    )
  )
}

#' Set up directory structure
#'
#' @param main_dir Main output directory
#' @return List of directory paths
setup_directories <- function(main_dir) {
  figures <- file.path(main_dir, "figures")
  dirs <- list(
    processed             = file.path(main_dir, "processed_data"),
    qc                    = file.path(main_dir, "qc"),
    dim_reduction         = file.path(main_dir, "dimensionality_reduction"),
    cnv                   = file.path(main_dir, "cnv"),
    reference_projection  = file.path(main_dir, "reference_projection"),
    figures               = figures,
    figures_qc            = file.path(figures, "qc"),
    figures_dim_reduction = file.path(figures, "dim_reduction"),
    figures_cnv           = file.path(figures, "cnv"),
    reports               = file.path(main_dir, "reports")
  )

  for (dir in dirs) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }

  return(dirs)
}

#' Log message to file and console
#'
#' @param message Message to log
#' @param log_file Log file path
log_message <- function(message, log_file = NULL) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  log_message <- paste(timestamp, message)
  
  # Print to console
  cat(log_message, "\n")
  
  # Write to log file if provided
  if (!is.null(log_file)) {
    write(log_message, file = log_file, append = TRUE)
  }
}