# MeQTrack

**Me**thylation **Q**uality control and analysis **Track**ing — a
single-user desktop app for Illumina methylation arrays (450K, EPIC,
EPICv2). MeQTrack runs the analysis (QC, dimensionality reduction,
reference projection, copy-number variation, self-contained HTML report)
**and tracks every run**: each invocation writes a timestamped folder
with its samplesheet,
parameters, logs, and outputs, so any past run is reproducible and
re-openable.

The app runs entirely on your machine. No cloud, no upload, no account.
Your data never leaves your laptop.

## What it does

Drop a CSV samplesheet pointing at your IDAT files into the app. One
click runs a six-stage pipeline:

| Stage | What happens |
|---|---|
| **1. Preprocess** | Reads IDATs (auto-detects 450K / EPIC / EPICv2), applies SWAN normalization (or sesame's QCDB for EPICv2), produces β-values. |
| **2. QC and probe filtering** | Per-sample detection-p, failed-probe %, intensity medians, and the sesame **GCT bisulfite-conversion** score (gates Pass_QC). Adds informational sample-integrity signals — predicted **sex/karyotype** (incl. Loss-of-Y flag), **Horvath age**, **leukocyte fraction**, **rs-SNP identity** match, and per-sample **dye-bias** QQ. Then removes sex-chromosome probes, SNP-affected probes, cross-reactive probes, and applies array-specific keep-lists. |
| **3. Dimensionality reduction** | t-SNE, UMAP, and hierarchical clustering on the most variable probes (default 10,000). |
| **4. Reference projection** | Projects each sample onto a pre-built reference t-SNE embedding of thousands of labelled methylomes (COMET paediatric solid-tumour, Capper et al. CNS-tumour / GSE90496, or Koelsche et al. sarcoma / GSE140686), then assigns each sample a nearest reference tumour class by k-NN vote — a diagnostic hint with a confidence score and ambiguity / distant-from-reference flags. |
| **5. Copy-number variation** | Per-sample CNV via conumee2, genome-wide segment calls, frequency plot, and an in-app segment heatmap. |
| **6. Report** | A self-contained HTML report you can open in any browser or share by email. |

The Shiny UI surfaces every stage interactively:

- **Samplesheet tab** — pick a CSV, see per-row validation (missing
  IDATs, duplicate IDs, malformed rows) before running anything.
- **Run tab** — a Settings card with tunable parameters across QC,
  dim. reduction, reference projection, and CNV (detection-p threshold,
  failed-probe % per sample, # variable probes, t-SNE perplexity, UMAP
  neighbors, the reference-dataset picker, projection perplexity, and
  the k-NN neighbour count), Run / Cancel controls, the Stages panel
  with per-stage **▶ Run** buttons for re-running an individual stage,
  live log tail, and post-run actions.
- **Past runs tab** — every prior run in your workspace listed
  newest-first with status, last step, and sample count. Open any row
  to attach it: result tabs render its artifacts, per-step Run buttons
  start operating against that run, and the Settings card auto-populates
  from its saved parameters.
- **QC tab** — sortable per-sample metrics table (incl. GCT, sex/karyotype,
  age, leukocyte, SNP identity), a pairwise SNP-identity concordance heatmap,
  per-sample dye-bias QQ plots, and embedded interactive density and MDS
  plots; QC-fail samples are styled distinctly across every downstream view.
- **Dim. reduction tab** — interactive plotly t-SNE, UMAP, and
  dendrogram. Color points by any metadata column from your samplesheet
  (Sample_Group, Diagnosis, Batch, etc.). A **Reference projection**
  sub-tab places your samples (dark diamonds) onto the chosen reference
  cloud, with a per-sample table of nearest class, confidence, and
  ambiguity / distance flags.
- **CNV tab** — per-sample genome-wide CNV profile (PDF embed),
  population frequency plot, and an in-browser segment heatmap with a
  tunable color-scale cap.
- **Report tab** — in-app preview of the generated HTML report, plus
  one-click "Open in new tab" and "Show in Finder/Explorer".
- **Help tab** — single-page Getting Started reference; the in-app
  walkthrough.

Result tabs populate **as each stage finishes**, not only at end-of-run
— you can review QC the moment QC completes, without waiting for CNV
or the report. Combined with per-step Run buttons, you can iterate on a
single stage (e.g. retry dim. reduction with a different perplexity)
without re-running preprocessing.

## Installing & running

End-user install is **unzip + double-click**. Full instructions in
[QUICKSTART.md](QUICKSTART.md).

Prerequisites: **R ≥ 4.4** and **pandoc ≥ 2.x**. Everything else
(Bioconductor, conumee2, yamapData, sesame) installs automatically on
first launch.

Supported hosts: macOS 13+ and Windows 10+. Linux is best-effort.

## Workflow at a glance

```
Your IDATs  ─►  samplesheet.csv  ─►  Samplesheet tab (validate)
                                          │
                                          ▼
                                    Run tab (analyze)
                                          │
                                          ▼
    QC ── Dim. red. ── Ref. projection ── CNV ── Report (review in app)
                                          │
                                          ▼
                            self-contained HTML report
                            (shareable; opens offline)
```

The pipeline writes everything to your workspace folder
(`~/MeQTrack/runs/<timestamp>_<samplesheet>/` by default), so every run
is reproducible and you keep your raw data alongside its outputs.

## For developers / CLI users

- **CLI usage** of the underlying R pipeline (no UI required): see
  [`pipeline/README.md`](pipeline/README.md).
- **Working from a clone** (instead of the release zip): the
  "Developer setup" section of [QUICKSTART.md](QUICKSTART.md) walks
  through `Rscript setup.R` and the manual launch path.

## Limits & caveats

- One run at a time. The app doesn't queue concurrent runs.
- Sample-size envelope: 1–96 samples per run; tested up to 96 on a
  16 GB machine.
- The Settings card exposes tunable parameters across QC, dim.
  reduction, reference projection, and CNV. Other pipeline knobs —
  normalization method, filtering toggles — are still at their defaults;
  expanding the Settings surface is on the v1.2 roadmap.
- Reference projection uses a Python toolchain (openTSNE) that
  provisions a self-contained environment on first use — a one-time
  few-minute download, so the first projection run is slower.
- The app listens only on `127.0.0.1` (loopback). It is not
  network-accessible, by construction.
