# app/R/help_module.R
# ---------------------------------------------------------------------------
# Wave 6 Theme 6c — In-app Getting Started page.
#
# Single-page reference rendered in the "Help" tab. Static content; the
# server-side scaffold exists so we can later swap in dynamic links (e.g.
# "open the workspace folder", "show the bundled example samplesheet path")
# without touching app.R.
# ---------------------------------------------------------------------------

help_module_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    class = "help-page",
    style = "max-width: 920px;",

    # --- Getting started -----------------------------------------------
    shiny::h3("Getting started"),
    shiny::p(
      class = "text-muted",
      "MeQTrack runs Illumina methylation array QC, dimensionality reduction, ",
      "and copy-number analysis on your laptop. No data leaves your machine."
    ),
    shiny::tags$ol(
      shiny::tags$li(
        shiny::tags$strong("Load a samplesheet."), " ",
        "Open the ", shiny::tags$em("Samplesheet"), " tab and pick a CSV. ",
        "The bundled example lives at ",
        shiny::tags$code("pipeline/data/example/samplesheet_epic.csv"), ". ",
        "Per-row validation runs immediately — green rows are good, ",
        "yellow/red rows have an issue (missing IDAT, duplicate ID, etc.)."
      ),
      shiny::tags$li(
        shiny::tags$strong("Optionally tweak Settings."), " ",
        "On the ", shiny::tags$em("Run"), " tab, the Settings card exposes ",
        "five tunable parameters across QC and dimensionality reduction. ",
        "Defaults work for typical small cohorts; hover the ",
        shiny::icon("circle-info"), " icons for guidance."
      ),
      shiny::tags$li(
        shiny::tags$strong("Click Run analysis."), " ",
        "The Stages panel shows live per-step progress, the log tail ",
        "streams the pipeline's output, and the result tabs ",
        "(QC, Dim. reduction, CNV, Report) populate as each stage finishes ",
        "— you don't have to wait for the whole pipeline to view QC."
      )
    ),

    shiny::hr(),

    # --- The tabs ------------------------------------------------------
    shiny::h3("The tabs at a glance"),
    shiny::tags$dl(
      class = "row",
      shiny::tags$dt(class = "col-sm-3", "Samplesheet"),
      shiny::tags$dd(class = "col-sm-9",
        "Pick a CSV, validate every row, see the auto-detected array type ",
        "(450K / EPIC / EPICv2). Gates the Run tab — invalid samplesheets ",
        "can't launch."
      ),
      shiny::tags$dt(class = "col-sm-3", "Run"),
      shiny::tags$dd(class = "col-sm-9",
        "Settings card on top, then Run / Cancel controls, the Stages panel ",
        "with per-step Run buttons (see ", shiny::tags$em("Per-step execution"),
        " below), the live log tail, and post-run actions ",
        "(open report, reveal in Finder)."
      ),
      shiny::tags$dt(class = "col-sm-3", "Past runs"),
      shiny::tags$dd(class = "col-sm-9",
        "Every prior run in your workspace, newest first. Click ",
        shiny::tags$strong("Open"), " on any row to attach it — the result ",
        "tabs render its artifacts and the per-step Run buttons start ",
        "operating against that run."
      ),
      shiny::tags$dt(class = "col-sm-3", "QC"),
      shiny::tags$dd(class = "col-sm-9",
        "Per-sample pass/fail, detection p, failed-probe %, intensity ",
        "medians, and the sesame GCT bisulfite-conversion score (gates ",
        "Pass_QC). Informational sample-integrity columns: predicted sex, ",
        "karyotype (incl. a Loss-of-Y flag), Horvath epigenetic age, ",
        "leukocyte fraction, and an rs-SNP identity match (SNP_Match_Pct). ",
        "Sub-tabs: Sample identity heatmap, per-sample Dye-bias QQ plots, ",
        "and embedded interactive density and MDS plots."
      ),
      shiny::tags$dt(class = "col-sm-3", "Dim. reduction"),
      shiny::tags$dd(class = "col-sm-9",
        "Interactive plotly t-SNE, UMAP, and dendrogram. Color points by ",
        "any metadata column (Sample_Group, Diagnosis, Batch...). QC-fail ",
        "samples render as triangles so they stay identifiable."
      ),
      shiny::tags$dt(class = "col-sm-3", "CNV"),
      shiny::tags$dd(class = "col-sm-9",
        "Per-sample genome-wide CNV profile (PDF embed), population ",
        "frequency plot, and an in-browser segment heatmap with a tunable ",
        "color-scale cap."
      ),
      shiny::tags$dt(class = "col-sm-3", "Report"),
      shiny::tags$dd(class = "col-sm-9",
        "In-app preview of the generated HTML report. Open in a new tab ",
        "to share by email — the report is fully self-contained."
      )
    ),

    shiny::hr(),

    # --- Pipeline stages -----------------------------------------------
    shiny::h3("The 6 pipeline stages"),
    shiny::tags$table(
      class = "table table-sm table-borderless",
      shiny::tags$thead(
        shiny::tags$tr(
          shiny::tags$th("Stage"),
          shiny::tags$th("What happens")
        )
      ),
      shiny::tags$tbody(
        shiny::tags$tr(
          shiny::tags$td(shiny::tags$strong("1. Preprocess")),
          shiny::tags$td(
            "Reads IDATs (auto-detects 450K / EPIC / EPICv2), applies SWAN ",
            "normalization (or sesame's QCDB for EPICv2), produces β-values."
          )
        ),
        shiny::tags$tr(
          shiny::tags$td(shiny::tags$strong("2. QC and probe filtering")),
          shiny::tags$td(
            "Per-sample detection p, failed-probe %, intensity medians, and ",
            "the sesame GCT bisulfite-conversion score; flags samples that ",
            "fall outside thresholds (GCT and detection failures gate ",
            "Pass_QC). Also derives informational sample-integrity signals ",
            "(sex/karyotype, age, leukocyte fraction, rs-SNP identity, ",
            "dye-bias QQ). Then removes sex-chromosome probes, SNP-affected ",
            "probes, cross-reactive probes, and applies array-specific keep-lists."
          )
        ),
        shiny::tags$tr(
          shiny::tags$td(shiny::tags$strong("3. Dimensionality reduction")),
          shiny::tags$td(
            "t-SNE, UMAP, and hierarchical clustering on the most variable ",
            "probes (default 10,000)."
          )
        ),
        shiny::tags$tr(
          shiny::tags$td(shiny::tags$strong("4. Reference projection")),
          shiny::tags$td(
            "Projects your samples onto a chosen reference t-SNE embedding ",
            "(COMET, Capper, or the GSE140686 sarcoma set) and assigns each ",
            "a nearest reference tumour group."
          )
        ),
        shiny::tags$tr(
          shiny::tags$td(shiny::tags$strong("5. Copy-number variation")),
          shiny::tags$td(
            "Per-sample CNV via conumee2, segment calls (.seg), genome-wide ",
            "frequency plot. Gain and loss are called against separate ",
            "thresholds (defaults: gain > 0.18, loss < -0.20)."
          )
        ),
        shiny::tags$tr(
          shiny::tags$td(shiny::tags$strong("6. Report")),
          shiny::tags$td(
            "Self-contained HTML report stitching every result together. ",
            "Open it in any browser, share by email."
          )
        )
      )
    ),

    shiny::hr(),

    # --- Per-step execution --------------------------------------------
    shiny::h3("Per-step execution"),
    shiny::p(
      "Each row in the Stages panel has a small ",
      shiny::tags$strong(shiny::icon("play"), " Run"),
      " button. Click it to re-run only that stage against the current ",
      "run directory — useful for tweaking a Setting and just re-doing ",
      "Dimensionality reduction, for example, without re-running ",
      "preprocessing."
    ),
    shiny::tags$ul(
      shiny::tags$li(
        "Buttons are disabled while a run is in progress, and disabled ",
        "until the upstream stage's outputs exist on disk. Hover over a ",
        "disabled button for the reason."
      ),
      shiny::tags$li(
        "Per-step runs reuse the current run directory; ",
        shiny::tags$strong("Run analysis"), " (the main button) creates ",
        "a fresh one."
      ),
      shiny::tags$li(
        "Combine with Past runs: open a prior run, then click Run on the ",
        "stages you want to redo with new Settings."
      )
    ),

    shiny::hr(),

    # --- Reference projection ------------------------------------------
    shiny::h3("Reference projection"),
    shiny::p(
      "The ", shiny::tags$em("Reference projection"), " sub-tab (under ",
      shiny::tags$em("Dim. reduction"), ") places your samples onto a ",
      "pre-built reference cohort and reports, per sample, the nearest ",
      "reference tumour group as a diagnostic hint. Choose which reference ",
      "to project onto with ", shiny::tags$em("Reference dataset"),
      " in the Run tab's Settings — the bundled COMET paediatric ",
      "solid-tumour set (~1,900 labelled methylomes), the Capper et al. ",
      "CNS-tumour reference (~2,800 methylomes, GSE90496), or the Koelsche ",
      "et al. sarcoma reference (1,077 methylomes, GSE140686)."
    ),
    shiny::tags$ul(
      shiny::tags$li(
        "The scatter shows your samples as dark diamonds on the reference ",
        "cloud; the table below lists every sample with its nearest class, ",
        "confidence, and ambiguity / distance flags."
      ),
      shiny::tags$li(
        "It runs as part of a full analysis, or on its own — the ",
        shiny::tags$em("Reference projection"), " stage has its own Run ",
        "button and needs only a samplesheet (no preprocessing first)."
      ),
      shiny::tags$li(
        shiny::tags$strong("First run is slower."), " The projection uses ",
        "a Python toolchain (openTSNE) that provisions a self-contained ",
        "environment on first use — a one-time few-minute download."
      )
    ),

    shiny::hr(),

    # --- Past runs -----------------------------------------------------
    shiny::h3("Past runs library"),
    shiny::p(
      "Every run writes to ",
      shiny::tags$code("~/MeQTrack/runs/<timestamp>_<samplesheet>/"),
      " with its own samplesheet, parameters, logs, and outputs. The ",
      shiny::tags$em("Past runs"), " tab lists them all with status, last ",
      "step, exit code, and sample count. ",
      shiny::tags$strong("Open"), " any run to attach it as the active ",
      "run — result tabs and per-step Run buttons start operating against ",
      "that directory, and the Settings card auto-populates from that ",
      "run's saved parameters."
    ),

    shiny::hr(),

    # --- Settings ------------------------------------------------------
    shiny::h3("Settings (tunable parameters)"),
    shiny::p(
      "Nine controls on the Run tab. Defaults match ",
      shiny::tags$code("pipeline_modules/config.R$default_config"),
      "; overrides apply to the next launch only and are persisted in the ",
      "run's ", shiny::tags$code("run_manifest.json"),
      " for reproducibility."
    ),
    shiny::tags$ul(
      shiny::tags$li(shiny::tags$strong("Detection-p threshold"),
        " (default 0.01) — per-probe cutoff. A probe is failed in a sample ",
        "if its p-value exceeds this."),
      shiny::tags$li(shiny::tags$strong("Failed-probe % per sample"),
        " (default 25) — a sample fails QC when more than this share of ",
        "its own probes failed."),
      shiny::tags$li(shiny::tags$strong("# variable probes"),
        " (default 10,000) — top variable CpGs (by SD) used for ",
        "dim. reduction."),
      shiny::tags$li(shiny::tags$strong("t-SNE perplexity"),
        " (default 5) — must be smaller than the sample count; small ",
        "cohorts often need 2–3."),
      shiny::tags$li(shiny::tags$strong("UMAP neighbors"),
        " (default 15) — local-vs-global trade-off."),
      shiny::tags$li(shiny::tags$strong("CNV gain threshold"),
        " (default 0.18) — seg.mean cutoff above which a segment is called a gain."),
      shiny::tags$li(shiny::tags$strong("CNV loss threshold"),
        " (default -0.20) — seg.mean cutoff below which a segment is called a loss; ",
        "enter as a negative value."),
      shiny::tags$li(shiny::tags$strong("Reference dataset"),
        " (default COMET) — which labelled reference cohort to project ",
        "onto: the COMET paediatric solid-tumour set, the Capper ",
        "CNS-tumour reference (GSE90496), or the GSE140686 sarcoma ",
        "reference."),
      shiny::tags$li(shiny::tags$strong("Nearest-class k"),
        " (default 25) — reference neighbours that vote on each projected ",
        "sample's tumour class."),
      shiny::tags$li(shiny::tags$strong("Projection perplexity"),
        " (default 5) — neighbourhood size for placing your samples onto ",
        "the reference embedding.")
    ),

    shiny::hr(),

    # --- Where data lives ----------------------------------------------
    shiny::h3("Where your data lives"),
    shiny::p(
      "On first launch the app creates a workspace at ",
      shiny::tags$code("~/MeQTrack/"), " with three subfolders:"
    ),
    shiny::tags$ul(
      shiny::tags$li(shiny::tags$code("samplesheets/"),
                     " — drop your CSVs here for quick access."),
      shiny::tags$li(shiny::tags$code("idats/"),
                     " — convenient location for IDAT files; the samplesheet ",
                     "can point anywhere on disk."),
      shiny::tags$li(shiny::tags$code("runs/<id>/"),
                     " — every run lands here. Open ",
                     shiny::tags$code("reports/methylation_analysis_report.html"),
                     " for the standalone report.")
    ),

    shiny::hr(),

    # --- Troubleshooting -----------------------------------------------
    shiny::h3("Troubleshooting"),
    shiny::tags$ul(
      shiny::tags$li(
        shiny::tags$strong("A run failed mid-stage."), " ",
        "Click ", shiny::tags$em("Open log"),
        " in the post-run actions to see the pipeline's full output. ",
        "Most failures point at a missing IDAT or a samplesheet column issue."
      ),
      shiny::tags$li(
        shiny::tags$strong("CNV step crashes with a BiocParallel error."), " ",
        "Usually a fork-related issue with sesame on macOS. Re-running ",
        "often clears it; if it persists, set the env var ",
        shiny::tags$code("MEQTRACK_SERIAL=1"), " to force serial execution."
      ),
      shiny::tags$li(
        shiny::tags$strong("Report is missing after a run."), " ",
        "Confirm pandoc is installed (",
        shiny::tags$code("pandoc --version"),
        "). The pipeline falls back to text-only output without it."
      ),
      shiny::tags$li(
        "More install / launch questions live in ",
        shiny::tags$code("QUICKSTART.md"),
        " in the project root."
      )
    )
  )
}

help_module_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Static content for now — no reactive logic. The scaffold is here so
    # future links (e.g. "show workspace folder", "open example samplesheet")
    # can be wired in without touching app.R.
    invisible(NULL)
  })
}
