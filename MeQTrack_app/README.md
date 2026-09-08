# MeQTrack — DNA methylation array analysis and run tracking

**Me**thylation **Q**uality control and analysis **Track**ing is a desktop
application and command-line pipeline for Illumina DNA methylation arrays
(HumanMethylation450, EPIC, EPICv2). It runs a complete analysis — preprocessing,
quality control, probe filtering, dimensionality reduction, projection onto
labelled reference cohorts, copy-number variation, and a self-contained HTML
report — **and tracks every run**: each invocation writes a timestamped folder
holding its samplesheet, parameters, logs, and outputs, so any past run is
reproducible and re-openable.

Everything runs locally. The app binds only to `127.0.0.1`; there is no cloud
service, no upload, and no account. Patient data never leaves the machine.

This repository contains the full source — the R analysis pipeline, the Shiny
desktop application that wraps it, the vendored annotation assets, and the
provisioning scripts — and is intended to let a reader understand the method,
audit the analysis logic, and run the whole thing end-to-end on the bundled
example data.

Current version: **3.0.0**.

## What's in this repository

Three layers, plus the assets they depend on:

| Component | Folder | What it is |
|---|---|---|
| Analysis pipeline | [`pipeline/`](pipeline/README.md) | The engine. A CLI driver (`methylation_pipeline.R`) plus seven modules that implement every analysis step. Runs standalone with no UI. |
| Desktop application | `app/` | A Shiny UI that wraps the pipeline — samplesheet validation, per-stage execution, live logs, interactive result views, and a past-runs browser. |
| Annotation assets | [`Anno/`](Anno/README.md) | Vendored sesame annotation files (Horvath clock models, the EPICv2→EPIC probe map, CNV control panels) committed so the pipeline runs offline. |
| Reference cohorts | `reference/` | Pre-built t-SNE embeddings and class labels for the three reference datasets a query sample can be projected onto. |
| Provisioning | `setup.R`, `scripts/` | One-shot dependency provisioning via `renv`, plus release-build, CNV-control-build, and HPC-submission helpers. |

The application is a **thin wrapper** over the pipeline — it shells out to the
same CLI driver a terminal user would call, and reads the same on-disk run
directory. See [`pipeline/README.md`](pipeline/README.md) for the analysis logic
itself (every module, every parameter, every output file) and
[`QUICKSTART.md`](QUICKSTART.md) for installation and day-to-day use.

## Repository layout

```
MeQTrack_app/
├── pipeline/             # Analysis engine: CLI driver + 7 analysis modules
│   ├── methylation_pipeline.R      # single entry point (Rscript ... --step)
│   ├── pipeline_modules/           # preprocess, qc, filtering, dim_reduction,
│   │                               #   reference_projection, cnv, visualization
│   └── data/                       # probe keep-lists + bundled example IDATs
├── app/                  # Shiny desktop UI (app.R + one module per result tab)
├── Anno/                 # Vendored sesame annotation assets (see Anno/README.md)
├── reference/            # Reference t-SNE embeddings + per-sample class labels
├── scripts/              # Release build, CNV control build, HPC submission
├── docs/decisions/       # Architecture decision records
├── setup.R               # One-shot renv provisioning (CRAN + Bioconductor)
├── meqtrack.command      # macOS double-click launcher
├── meqtrack.bat          # Windows double-click launcher
├── QUICKSTART.md         # Install, use, and developer-setup guide
├── CONTRIBUTING.md       # How to report issues and propose changes
├── CODE_OF_CONDUCT.md    # Contributor Covenant
├── SECURITY.md           # Vulnerability and data-exposure reporting
└── LICENSE               # GNU General Public License v3.0
```

Each of `pipeline/` and `Anno/` has its own `README.md` with full detail. Start
with [`pipeline/README.md`](pipeline/README.md) to understand the analysis, then
[Data and privacy](#data-and-privacy) below to understand what data backs it.

## How the codebase is structured

The project follows the same three-layer design throughout:

1. **Modules** (`pipeline/pipeline_modules/`) — pure analysis functions, one file
   per stage. Each takes a beta matrix (or an `RGChannelSet`) plus a config list
   and writes its outputs to a canonical subdirectory. No knowledge of the CLI,
   and none of the UI.
2. **CLI driver** (`pipeline/methylation_pipeline.R`) — parses flags, builds the
   config, creates the output tree, and calls the modules in order. `--step`
   selects one stage or `all`; because each stage persists its results as
   `.RData`, any later stage can be re-run without repeating the ones before it.
3. **Application** (`app/`) — one Shiny module per result tab
   (`qc_module.R`, `dimred_module.R`, `cnv_module.R`, …), a `pipeline_bridge.R`
   that invokes the CLI driver as a subprocess, and a `results_loader.R` that
   reads whatever artifacts the run directory currently holds.

Because the application only ever talks to the pipeline through its CLI and its
output directory, the two are independently usable: the pipeline runs headless on
an HPC cluster, and the app can attach to and browse a run produced there.

## Analysis stages

The pipeline runs seven stages, selectable individually with `--step`:

- **Preprocess** — reads IDATs and auto-detects the platform from probe count.
  Beta values and pOOBAH detection p-values come from sesame's `openSesame()`
  with `prep = "QCDB"` (channel inference, dye-bias correction, Noob background
  correction) [1]; SWAN normalization [2] is used on the recovery path for
  low-intensity samples. EPICv2 replicate probes are collapsed to base CpG IDs by
  prefix. Also computes the per-sample GCT bisulfite-conversion control score.
- **Quality control** — per-sample mean detection p, failed-probe percentage, and
  pre-normalization channel intensity medians via minfi [3]. A sample fails QC on
  detection p, failed-probe rate, or GCT score crossing their thresholds; low
  intensity is deliberately informational rather than a failure, because scanner
  gain varies legitimately across sites. Adds sample-integrity signals that never
  gate `Pass_QC`: predicted sex, Horvath 353-CpG epigenetic age [4], leukocyte
  fraction, rs-SNP identity concordance, and per-sample dye-bias QQ.
- **Filtering** — a detection-p filter, then an array-specific curated keep-list
  that removes sex-chromosome, SNP-affected, and cross-reactive probes.
- **Dimensionality reduction** — the top *N* most variable probes (default
  10,000) fed to t-SNE, UMAP [5], and hierarchical clustering on a correlation
  distance.
- **Reference projection** — projects each query sample onto a pre-built
  reference t-SNE embedding with `snifter`/openTSNE [6], harmonizing to the
  reference probe set on the query side only. A k-NN vote in the embedding space
  assigns each sample a nearest reference class with a confidence score and
  *ambiguous* / *distant-from-reference* flags. Three references ship: a
  paediatric solid-tumour cohort (GSE305405, default), Capper et al. CNS tumours
  (GSE90496) [7], and Koelsche et al. sarcoma (GSE140686) [8].
- **Copy-number variation** — per-sample CNV via conumee2 [9] against internal or
  user-supplied control references, genome-wide segment calls, a population
  gain/loss frequency plot, and a multi-sample segment heatmap.
- **Report** — a self-contained HTML report rendered with R Markdown, embedding
  every figure. Degrades to a plain-text report if pandoc is unavailable.

Class assignments from reference projection are a **diagnostic hint, not a
diagnosis** — they carry a confidence score and ambiguity flags precisely because
they are intended to be read alongside, not in place of, standard workup.

## The application

The Shiny UI surfaces each stage interactively. Result tabs populate **as each
stage finishes**, not only at end-of-run, and each stage has its own **▶ Run**
button, so you can iterate on one stage (say, re-running dimensionality reduction
at a different perplexity) without repeating preprocessing.

- **Samplesheet** — pick a CSV and see per-row validation (missing IDATs,
  duplicate IDs, malformed rows) before running anything.
- **Run** — a Settings card exposing tunable parameters across QC, dimensionality
  reduction, reference projection, and CNV; Run/Cancel controls; the per-stage
  Stages panel; and a live log tail.
- **Past runs** — every prior run in the workspace, newest first, with status,
  last step, and sample count. Open a row to attach it: result tabs render its
  artifacts and the Settings card repopulates from its saved parameters.
- **QC** — sortable per-sample metrics, a pairwise SNP-identity concordance
  heatmap, dye-bias QQ plots, and interactive density and MDS plots. QC-fail
  samples are styled distinctly across every downstream view.
- **Dim. reduction** — interactive t-SNE, UMAP, and dendrogram, colourable by any
  samplesheet metadata column, plus a Reference projection sub-tab placing your
  samples on the chosen reference cloud.
- **CNV** — per-sample genome-wide profile, population frequency plot, and an
  in-browser segment heatmap with a tunable colour-scale cap.
- **Report** and **Help** — in-app report preview with open-externally and
  show-in-folder actions, and a single-page getting-started reference.

## Data and privacy

MeQTrack processes patient-derived methylation data, so what this repository does
and does not ship matters. **Only public, de-identified data is distributed here.**

**What ships.** Probe keep-lists for each array type; vendored sesame annotation
assets (see [`Anno/README.md`](Anno/README.md)); pre-built reference t-SNE
embeddings with the per-sample class and colour labels that annotate them; and a
four-sample EPIC example drawn from public GEO series
[GSE130295](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130295).

**No raw patient arrays are distributed.** An Illumina methylation array carries
rs genotyping probes and is therefore inherently identifiable — an IDAT cannot be
de-identified by editing its samplesheet. Only the public GEO IDATs above are
included. The `samplesheet_450k.csv` and `samplesheet_epicv2.csv` examples carry
synthetic identifiers and no IDATs: they illustrate the samplesheet format for
those platforms and are not runnable. See
[`pipeline/data/example/README.md`](pipeline/data/example/README.md).

**Reference labels are minimal by construction.** Reference projection reads only
three fields per reference sample — array barcode, tumour group, and plotting
colour — so `reference/COMET_reference_labels.csv` carries exactly those three
columns for the 1,915 samples in the embedding and nothing more. Clinical exports
are excluded by `.gitignore` and must not be committed; regenerate the label file
from your own cohort export if the reference changes. The `GSE90496` and
`GSE140686` label files are the published supplementary tables from their
respective papers [7,8].

**What does not ship**, because of file size — each must be obtained separately
before the corresponding stage will run:

| Not committed | Why | How to obtain |
|---|---|---|
| `reference/beta_*.rds` | 77–201 MB each, over GitHub's per-file limit | Distributed separately; required for reference projection |
| `pipeline/data/yamapData_*.tar.gz` | 257 MB | Obtained separately, placed in `pipeline/data/` before `setup.R`; conumee2's internal control panel depends on it |
| `Anno/*/*_CNV_controls.rds` | ~40 MB each | Rebuild from control IDATs with `scripts/build_cnv_controls.R` |
| `Anno/EPICv2/EPICv2ToEPIC_conversion.tsv` | 115 MB upstream table | Only the slim 3-column map is committed and used |

**No run outputs are committed.** `runs/`, `dist/`, and `*.log` are gitignored;
the pipeline writes every result to a workspace folder outside the repository
(`~/MeQTrack/runs/<timestamp>_<samplesheet>/` by default). Cluster submission
templates in `pipeline/command_line.txt` and `scripts/submit_hpc_lsf.sh` use
placeholder paths rather than real site infrastructure.

## Getting started

Prerequisites: **R ≥ 4.4** and **pandoc ≥ 2.x**. Supported hosts are macOS 13+
and Windows 10+; Linux is best-effort.

For end users, installation is unzip-and-double-click — see
[`QUICKSTART.md`](QUICKSTART.md). From a clone:

```bash
Rscript setup.R
```

This provisions every CRAN and Bioconductor dependency into a project-local
`renv` library and snapshots it to `renv.lock`. It is idempotent; re-running skips
what is already installed.

As a smoke test, run the full pipeline on the four bundled EPIC samples:

```bash
Rscript pipeline/methylation_pipeline.R \
  --input      pipeline/data/example/samplesheet_epic.csv \
  --output     pipeline/results_example \
  --data_dir   pipeline/data \
  --array_type EPIC \
  --threads    4 \
  --step       all
```

This runs preprocessing through report generation and writes a complete run
directory to `pipeline/results_example/`, ending with a self-contained
`reports/methylation_analysis_report.html`. See
[`pipeline/README.md`](pipeline/README.md) for the full CLI reference, the
samplesheet format, and per-stage detail.

To launch the application instead, double-click `meqtrack.command` (macOS) or
`meqtrack.bat` (Windows), or follow the manual launch path in
[`QUICKSTART.md`](QUICKSTART.md).

## Requirements

R ≥ 4.4 with, from **CRAN**: `optparse`, `data.table`, `ggplot2`, `plotly`,
`Rtsne`, `umap`, `dendextend`, `circlize`, `htmlwidgets`, `rmarkdown`, `knitr`,
`DT`, `yaml`, `ggrepel`, plus `shiny`, `bslib`, and `callr` for the application;
and from **Bioconductor**: `minfi`, `sesame`, `limma`, `missMethyl`,
`matrixStats`, `snifter`, `DMRcate`, `conumee2`, `GenomicRanges`, `Gviz`, and the
450K/EPIC/EPICv2 manifest and annotation packages. `yamapData` installs from a
local source tarball. `setup.R` handles all of it; the exact pinned versions live
in `renv.lock`.

**External tools:** `pandoc` for HTML report rendering. Reference projection
provisions a self-contained Python environment (openTSNE) on first use — a
one-time few-minute download, so the first projection run is slower.

Everything runs on CPU. A reference point for sizing: 265 EPIC samples, 4
threads, Mac Studio (M2 Max, 32 GB) — 1 h 15 min end-to-end. Runtime scales
roughly with sample count.

## Limits and caveats

- One run at a time; concurrent runs are not queued.
- Tested up to 265 samples in a single run. Larger runs are untested.
- The Settings card exposes QC, dimensionality-reduction, reference-projection,
  and CNV parameters. Other knobs — normalization method, filtering toggles —
  remain at their defaults.
- The `remove_snps` and `remove_cross_reactive` CLI flags are effectively no-ops:
  SNP and cross-reactive probe removal is baked into the curated keep-lists
  rather than applied live.
- The GCT bisulfite-conversion score is unavailable on EPICv2 and records `NA`
  there; `NA` scores never fail QC.

## Contributing

Bug reports and small fixes are welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to
get started, and please follow our [Code of Conduct](CODE_OF_CONDUCT.md). Architecture
decisions are recorded in [`docs/decisions/`](docs/decisions/); if a change alters pipeline
behaviour or output layout, add an ADR alongside it.

When reporting an issue, **never attach IDATs, beta matrices, samplesheets, or QC reports from
real samples** — they carry identifiable data. Reproduce with the bundled public example where
you can. To report a security issue, or any file in this repository that you believe contains
identifiable data, see [`SECURITY.md`](SECURITY.md).

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE), matching the
other components of the COMET repository. GPL-3.0 is also required in practice because
MeQTrack depends on and loads GPL-licensed Bioconductor packages, including `minfi`,
`conumee2`, and `DMRcate`.

## References

[1] Zhou W, Triche TJ Jr, Laird PW, Shen H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. Nucleic Acids Research. 2018;46(20):e123.

[2] Maksimovic J, Gordon L, Oshlack A. SWAN: Subset-quantile within array normalization for Illumina Infinium HumanMethylation450 BeadChips. Genome Biology. 2012;13(6):R44.

[3] Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, Irizarry RA. Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. Bioinformatics. 2014;30(10):1363-9.

[4] Horvath S. DNA methylation age of human tissues and cell types. Genome Biology. 2013;14(10):R115.

[5] McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426. 2018.

[6] Poličar PG, Stražar M, Zupan B. openTSNE: A Modular Python Library for t-SNE Dimensionality Reduction and Embedding. Journal of Statistical Software. 2024;109(3).

[7] Capper D, Jones DTW, Sill M, et al. DNA methylation-based classification of central nervous system tumours. Nature. 2018;555(7697):469-474. (GEO: GSE90496)

[8] Koelsche C, Schrimpf D, Stichel D, et al. Sarcoma classification by DNA methylation profiling. Nature Communications. 2021;12(1):498. (GEO: GSE140686)

[9] Daenekas B, Pérez E, Boniolo F, et al. Conumee 2.0: enhanced copy-number variation analysis from DNA methylation arrays for humans and mice. Bioinformatics. 2024;40(2):btae029.

[10] Paediatric solid-tumour methylation reference cohort, GEO: [GSE305405](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305405).
