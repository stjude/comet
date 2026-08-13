# app/R/settings_module.R
# ---------------------------------------------------------------------------
# Wave 6 Theme 6d Phase 1 — user-tunable pipeline parameters.
#
# Numeric knobs across QC, dimensionality reduction, CNV and reference
# projection, plus the reference-dataset picker (which labelled cohort a
# query is projected onto). The reactive
# returned by settings_module_server() is consumed by run_controller, which
# passes the values into bridge_launch on each run. Per the v1.1 decision,
# these apply to the next launch only — there is no workspace-default
# persistence layer (past runs remember their own settings via the
# run_manifest.json parameters block, which the past-runs attach feeds back
# into these inputs).
#
# Public exports:
#   settings_module_ui(id)
#   settings_module_server(id, attach_run) -> reactive list of parameters
# ---------------------------------------------------------------------------

# Authoritative defaults — must match pipeline_modules/config.R$default_config.
# All numeric except refproj_dataset (a reference-dataset key); the reset and
# past-run-restore handlers branch on type so the selectInput is updated with
# updateSelectInput rather than updateNumericInput.
SETTINGS_DEFAULTS <- list(
  qc_detection_p          = 0.01,
  qc_failed_pct           = 25,
  qc_max_gct              = 1.3,
  dim_variable_probes     = 10000L,
  dim_tsne_perplexity     = 5L,
  dim_umap_neighbors      = 15L,
  cnv_gain_threshold      =  0.18,
  cnv_loss_threshold      = -0.20,
  refproj_dataset         = "COMET_1915",
  refproj_knn_k           = 25L,
  refproj_perplexity      = 5L
)

# Build the "<name> [?] (default X)" label for a numeric input. The info
# icon is wrapped in a bslib tooltip so hovering reveals the explanation.
.settings_label <- function(name, default, help_text) {
  shiny::tagList(
    name, " ",
    bslib::tooltip(
      shiny::icon("circle-info", class = "text-muted",
                  style = "cursor: help;"),
      help_text,
      placement = "right"
    ),
    " ",
    shiny::tags$small(class = "text-muted",
                      sprintf("(default %s)", default))
  )
}

settings_module_ui <- function(id) {
  ns <- shiny::NS(id)
  d  <- SETTINGS_DEFAULTS
  shiny::tagList(
    shiny::div(
      class = "row g-3",
      shiny::div(
        class = "col-md-6",
        shiny::h6(class = "text-muted text-uppercase small fw-semibold",
                  "Quality control"),
        shiny::numericInput(
          ns("qc_detection_p"),
          .settings_label(
            "Detection-p threshold", d$qc_detection_p,
            paste(
              "Per-probe cutoff on the detection p-value.",
              "Within a sample, a probe is counted as failed when its",
              "p-value exceeds this. Lower values = stricter calls."
            )
          ),
          value = d$qc_detection_p, min = 0, max = 1, step = 0.01
        ),
        shiny::numericInput(
          ns("qc_failed_pct"),
          .settings_label(
            "Failed-probe % per sample", d$qc_failed_pct,
            paste(
              "Maximum share of failed probes a sample is allowed to have",
              "before being flagged. Computed within each sample (failed",
              "probes ÷ total probes × 100), not across the cohort.",
              "A sample at or above this percent fails QC."
            )
          ),
          value = d$qc_failed_pct, min = 0, max = 100, step = 1
        ),
        shiny::numericInput(
          ns("qc_max_gct"),
          .settings_label(
            "Max GCT (bisulfite conversion)", d$qc_max_gct,
            paste(
              "GCT bisulfite-conversion control score (sesame, Zhou et al.",
              "2017). A value near 1.0 means complete conversion; higher",
              "means more incomplete. A sample whose GCT exceeds this value",
              "fails QC. Samples without a GCT score (e.g. EPICv2) are never",
              "failed on it."
            )
          ),
          value = d$qc_max_gct, min = 1, max = 5, step = 0.05
        )
      ),
      shiny::div(
        class = "col-md-6",
        shiny::h6(class = "text-muted text-uppercase small fw-semibold",
                  "Dimensionality reduction"),
        shiny::numericInput(
          ns("dim_variable_probes"),
          .settings_label(
            "# variable probes", d$dim_variable_probes,
            paste(
              "How many of the most variable CpG probes (by standard",
              "deviation across samples) to feed into t-SNE / UMAP /",
              "clustering. More probes = more signal but slower; 10,000",
              "is a common starting point."
            )
          ),
          value = d$dim_variable_probes, min = 100, max = 1e6, step = 1000
        ),
        shiny::numericInput(
          ns("dim_tsne_perplexity"),
          .settings_label(
            "t-SNE perplexity", d$dim_tsne_perplexity,
            paste(
              "t-SNE's neighborhood size — roughly the effective number",
              "of neighbors each point considers. Use small values (2–5)",
              "for tiny cohorts (<15 samples); 30 is typical for hundreds.",
              "Must be less than the number of samples."
            )
          ),
          value = d$dim_tsne_perplexity, min = 1, max = 100, step = 1
        ),
        shiny::numericInput(
          ns("dim_umap_neighbors"),
          .settings_label(
            "UMAP neighbors", d$dim_umap_neighbors,
            paste(
              "UMAP's local-vs-global trade-off. Smaller values preserve",
              "fine local structure (clusters); larger values preserve",
              "global topology at the cost of detail. 15 is the package",
              "default."
            )
          ),
          value = d$dim_umap_neighbors, min = 2, max = 200, step = 1
        )
      ),
      shiny::div(
        class = "col-md-6",
        shiny::h6(class = "text-muted text-uppercase small fw-semibold",
                  "Copy-number variation"),
        shiny::numericInput(
          ns("cnv_gain_threshold"),
          .settings_label(
            "CNV gain threshold", d$cnv_gain_threshold,
            paste(
              "seg.mean cutoff above which a segment is called a gain",
              "(seg.mean > this). 0.18 is the conumee default. Lower",
              "values = more sensitive."
            )
          ),
          value = d$cnv_gain_threshold, min = 0, max = 1, step = 0.01
        ),
        shiny::numericInput(
          ns("cnv_loss_threshold"),
          .settings_label(
            "CNV loss threshold", d$cnv_loss_threshold,
            paste(
              "seg.mean cutoff below which a segment is called a loss",
              "(seg.mean < this). Enter as a negative value, e.g. -0.2.",
              "Higher (closer to 0) = more sensitive."
            )
          ),
          value = d$cnv_loss_threshold, min = -1, max = 0, step = 0.01
        )
      ),
      shiny::div(
        class = "col-md-6",
        shiny::h6(class = "text-muted text-uppercase small fw-semibold",
                  "Reference projection"),
        shiny::selectInput(
          ns("refproj_dataset"),
          # No "(default ...)" suffix here — the dataset keys are long; the
          # tooltip carries the explanation instead.
          shiny::tagList(
            "Reference dataset", " ",
            bslib::tooltip(
              shiny::icon("circle-info", class = "text-muted",
                          style = "cursor: help;"),
              paste(
                "Which labelled reference cohort your samples are projected",
                "onto. Each dataset has its own t-SNE embedding, training",
                "beta matrix and tumour-class labels — the projection and",
                "the nearest-class hints are computed against the one",
                "picked here. COMET is a paediatric solid-tumour set;",
                "Capper is the CNS-tumour reference (GSE90496); GSE140686",
                "is the sarcoma reference."
              ),
              placement = "right"
            )
          ),
          choices  = reference_dataset_choices(),
          selected = d$refproj_dataset
        ),
        shiny::numericInput(
          ns("refproj_knn_k"),
          .settings_label(
            "Nearest-class k", d$refproj_knn_k,
            paste(
              "Number of nearest reference samples that vote on the tumour",
              "class for each projected sample. Larger = smoother but can",
              "blur small classes; 25 is the default."
            )
          ),
          value = d$refproj_knn_k, min = 1, max = 200, step = 1
        ),
        shiny::numericInput(
          ns("refproj_perplexity"),
          .settings_label(
            "Projection perplexity", d$refproj_perplexity,
            paste(
              "Neighbourhood size used when placing your samples onto the",
              "reference embedding (snifter's projection step). Distinct",
              "from the t-SNE perplexity above; 5 is the default."
            )
          ),
          value = d$refproj_perplexity, min = 1, max = 100, step = 1
        )
      )
    ),
    shiny::div(
      class = "mt-2",
      shiny::actionButton(
        ns("reset"), "Reset to defaults",
        icon = shiny::icon("arrow-rotate-left"),
        class = "btn-sm btn-outline-secondary"
      )
    )
  )
}

settings_module_server <- function(id, attach_run = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observeEvent(input$reset, {
      for (k in names(SETTINGS_DEFAULTS)) {
        v <- SETTINGS_DEFAULTS[[k]]
        if (is.character(v)) {
          shiny::updateSelectInput(session, k, selected = v)
        } else {
          shiny::updateNumericInput(session, k, value = v)
        }
      }
    })

    # When a past run is attached, populate inputs from its manifest if it
    # carries a parameters block. Missing keys leave the current value alone.
    shiny::observeEvent(attach_run(), {
      path <- attach_run()
      if (is.null(path)) return()
      params <- read_run_parameters(path)
      if (is.null(params)) return()
      apply_params_to_inputs(session, params)
    }, ignoreInit = TRUE)

    # Reactive bag consumed by run_controller. Each input falls back to its
    # default if the user blanks the field.
    shiny::reactive({
      list(
        qc.detection_p_threshold          = numeric_or_default(
          input$qc_detection_p, SETTINGS_DEFAULTS$qc_detection_p),
        qc.failed_probe_percent_threshold = numeric_or_default(
          input$qc_failed_pct, SETTINGS_DEFAULTS$qc_failed_pct),
        qc.max_gct_score                  = numeric_or_default(
          input$qc_max_gct, SETTINGS_DEFAULTS$qc_max_gct),
        dim.variable_probes               = integer_or_default(
          input$dim_variable_probes, SETTINGS_DEFAULTS$dim_variable_probes),
        dim.tsne_perplexity               = integer_or_default(
          input$dim_tsne_perplexity, SETTINGS_DEFAULTS$dim_tsne_perplexity),
        dim.umap_n_neighbors              = integer_or_default(
          input$dim_umap_neighbors, SETTINGS_DEFAULTS$dim_umap_neighbors),
        cnv.gain_threshold                = numeric_or_default(
          input$cnv_gain_threshold, SETTINGS_DEFAULTS$cnv_gain_threshold),
        cnv.loss_threshold                = numeric_or_default(
          input$cnv_loss_threshold, SETTINGS_DEFAULTS$cnv_loss_threshold),
        refproj.dataset                   = character_or_default(
          input$refproj_dataset, SETTINGS_DEFAULTS$refproj_dataset),
        refproj.knn_k                     = integer_or_default(
          input$refproj_knn_k, SETTINGS_DEFAULTS$refproj_knn_k),
        refproj.perplexity                = integer_or_default(
          input$refproj_perplexity, SETTINGS_DEFAULTS$refproj_perplexity)
      )
    })
  })
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

numeric_or_default <- function(x, default) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  as.numeric(x)
}

integer_or_default <- function(x, default) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  as.integer(x)
}

character_or_default <- function(x, default) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) return(default)
  as.character(x)
}

# Available reference datasets for the projection-dataset picker, read from
# the pipeline's authoritative registry (.REFERENCE_DATASETS in
# reference_projection.R) so the Settings panel and the pipeline never drift
# as datasets are added. Returns a named character vector (label -> key)
# suitable for selectInput(); falls back to the COMET default if the
# registry file cannot be located or parsed.
#
# reference_projection.R has no top-level side effects (it only defines the
# registry list and functions; library() calls live inside function bodies),
# so sourcing it into a throwaway environment is cheap and safe.
reference_dataset_choices <- function() {
  fallback <- c("Pour et al. COMET pediatric solid-tumor (1915 samples, GSE305405)" = "COMET_1915")
  reg_file <- tryCatch(
    file.path(pipeline_dir(), "pipeline_modules", "reference_projection.R"),
    error = function(e) NULL
  )
  if (is.null(reg_file) || !file.exists(reg_file)) return(fallback)
  reg <- tryCatch({
    e <- new.env()
    sys.source(reg_file, envir = e)
    get(".REFERENCE_DATASETS", envir = e)
  }, error = function(e) NULL)
  if (is.null(reg) || !length(reg)) return(fallback)
  stats::setNames(
    names(reg),
    vapply(reg, function(d) d$label %||% "(unnamed dataset)", character(1))
  )
}

# Read the parameters block from a run's manifest. Returns NULL if missing.
read_run_parameters <- function(run_dir) {
  if (is.null(run_dir) || !dir.exists(run_dir)) return(NULL)
  m_path <- file.path(run_dir, "run_manifest.json")
  if (!file.exists(m_path)) return(NULL)
  m <- tryCatch(jsonlite::read_json(m_path, simplifyVector = TRUE),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  m$parameters
}

# Update the input widgets from a saved parameters block.
# CNV thresholds: prefer the new gain/loss keys, but fall back to the
# legacy single `cnv.threshold` (older run_manifest.json files) by
# applying it symmetrically as ±abs(threshold).
apply_params_to_inputs <- function(session, params) {
  legacy_cnv <- params$cnv.threshold
  cnv_gain <- params$cnv.gain_threshold
  if (is.null(cnv_gain) && !is.null(legacy_cnv)) cnv_gain <-  abs(legacy_cnv)
  cnv_loss <- params$cnv.loss_threshold
  if (is.null(cnv_loss) && !is.null(legacy_cnv)) cnv_loss <- -abs(legacy_cnv)

  pairs <- list(
    qc_detection_p      = params$qc.detection_p_threshold,
    qc_failed_pct       = params$qc.failed_probe_percent_threshold,
    qc_max_gct          = params$qc.max_gct_score,
    dim_variable_probes = params$dim.variable_probes,
    dim_tsne_perplexity = params$dim.tsne_perplexity,
    dim_umap_neighbors  = params$dim.umap_n_neighbors,
    cnv_gain_threshold  = cnv_gain,
    cnv_loss_threshold  = cnv_loss,
    refproj_dataset     = params$refproj.dataset,
    refproj_knn_k       = params$refproj.knn_k,
    refproj_perplexity  = params$refproj.perplexity
  )
  for (input_id in names(pairs)) {
    v <- pairs[[input_id]]
    if (is.null(v) || length(v) != 1L || is.na(v)) next
    if (input_id == "refproj_dataset") {
      shiny::updateSelectInput(session, input_id, selected = as.character(v))
    } else {
      shiny::updateNumericInput(session, input_id, value = v)
    }
  }
}
