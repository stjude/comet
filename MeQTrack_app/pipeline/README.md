# MeQTrack R pipeline

A CLI-driven DNA methylation analysis pipeline for Illumina 450k, EPIC, and EPICv2 arrays. This folder contains the R source copied from `/Volumes/qtran/MeQTrack` and is the underlying engine that the MeQTrack desktop MVP will wrap.

## What it is

`methylation_pipeline.R` is the single entry point — a CLI driver (`Rscript methylation_pipeline.R …`) that runs a seven-step Illumina methylation array analysis end-to-end. It supports 450k, EPIC, and EPICv2 arrays, with `--array_type auto` inferring the platform from probe count. The seven steps are `preprocess → qc → filtering → dim_reduction → reference_projection → cnv → visualization`, and the `--step` flag lets you run all of them or just one. Each module is a sourced R file in `pipeline_modules/`.

A second, independent script, `cnv_heatmap.R`, is a standalone tool for producing a multi-sample CNV heatmap from an already-computed segmentation file — rows are samples, x-axis is genome coordinate, color is `seg.mean`, with optional metadata-driven side annotation bars. It's not called by the main pipeline; it's a downstream reporting utility.

## Quick start with the bundled example data

The `data/example/` folder ships with 4 EPIC (850k) samples from GEO series [GSE130295](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130295) (melanoma). Each sample contributes a paired `_Red.idat` / `_Grn.idat`, named `GSM3735546…GSM3735549_201465940014_R0{1..4}C01`. The matching samplesheet, `samplesheet_epic.csv`, sits next to the IDATs and references them by relative path.

To run the full pipeline on these 4 samples, from anywhere on disk:

```
Rscript pipeline/methylation_pipeline.R \
  --input      pipeline/data/example/samplesheet_epic.csv \
  --output     pipeline/results_example \
  --data_dir   pipeline/data \
  --array_type EPIC \
  --threads    4 \
  --step       all
```

This works regardless of where you invoke `Rscript` from because `methylation_pipeline.R` does `setwd(script_dir)` at startup, making the pipeline folder the effective working directory. The `Basename` values in `data/example/samplesheet_epic.csv` are therefore resolved as `data/example/GSM3735546_…R01C01` relative to that working directory.

After the run completes you should see:

```
pipeline/results_example/
├── processed_data/               beta_values.txt, preprocessed_data.RData, detection_p.txt, rgset.RData
├── qc/                           sample_qc_report.csv (8 rows), qc_results.RData
├── dimensionality_reduction/     tsne_results.RData, umap_results.RData, hclust_results.RData
├── reference_projection/         projected coords CSV, class-hints CSV, projection PDF
├── cnv/                          cnv_results.RData, segments/{GSM3735546…}.seg
├── figures/                      per-module PDFs + interactive HTML widgets
├── reports/                      methylation_analysis_report.html
└── pipeline_log.txt
```

To run a single step (e.g. re-generating only the CNV outputs from previously saved `preprocessed_data.RData`):

```
Rscript pipeline/methylation_pipeline.R \
  --input      pipeline/data/example/samplesheet_epic.csv \
  --output     pipeline/results_example \
  --data_dir   pipeline/data \
  --array_type EPIC \
  --step       cnv
```

## Samplesheet format

The `--input` CSV requires three columns:

| Column | Description |
|---|---|
| `Sentrix_ID` | Unique per-sample identifier, typically `<GSM>_<chip>_<RowCol>` |
| `Sample_Name` | Human-readable name used in plots and labels |
| `Basename` | Path **stem** of the IDAT pair (no `_Red.idat` / `_Grn.idat` suffix), absolute or relative to the pipeline folder |

Optional metadata columns (e.g. `Gender`, `Group`, `Batch`) are carried through to downstream plots and can be used for coloring — the bundled `samplesheet_epic.csv` includes `Gender` and a `path` column (the production HPC path) to show how additional metadata travels alongside the required three columns.

`command_line.txt` contains the exact production `bsub` commands used at St. Jude if you need to submit to LSF instead of running locally.

## Full CLI reference

```
Rscript methylation_pipeline.R \
  --input      <samplesheet.csv> \
  --output     <results_dir> \
  --array_type {450K|EPIC|EPICv2|auto} \
  --threads    <N> \
  --step       {all|preprocess|qc|filtering|dim_reduction|reference_projection|cnv|visualization} \
  [--config    <config.R|config.yaml>] \
  [--data_dir  <dir-containing-keep.probes.*.txt>] \
  [--hpc]       # generate HPC submission scripts and exit
```

## How a run is organized

The driver reads the input samplesheet, deduplicates on `Basename`, then calls `setup_directories()` (in `config.R`) to create a canonical output tree under `--output`:

```
<output>/
├── processed_data/               beta values, M values, rgset, sample_info, detection_p
├── qc/                           sample_qc_report.csv, qc_results.RData
├── dimensionality_reduction/     tsne_results.RData, umap_results.RData, hclust_results.RData
├── reference_projection/         projected coords CSV, class-hints CSV, projection PDF
├── cnv/                          cnv_results.RData, segments/
├── figures/{qc,dim_reduction,cnv}  PDFs and interactive HTML widgets
├── reports/                      methylation_analysis_report.html (or .txt fallback)
└── pipeline_log.txt
```

Each step persists its outputs as `.RData` so later steps can be re-run independently — e.g. `--step cnv` after a successful `--step preprocess` skips re-reading IDATs. If you run `--step all` against an existing directory, it overwrites in place.

## Step-by-step

**Preprocess (`pipeline_modules/preprocess.R`).** Reads the samplesheet, deduplicates on `Basename`, loads IDATs into an `RGChannelSet` via `minfi::read.metharray.exp()`, and runs sesame's `openSesame()` pipeline with `prep = "QCDB"` (QC → channel inference → dye bias → Noob background correction) to produce beta values and pOOBAH detection p-values. For EPICv2 arrays it additionally collapses replicate probes by prefix using the mean. It also predicts sex from raw methylation intensities (`minfi::getSex`). It computes the per-sample **GCT bisulfite-conversion control score** (sesame `bisConversionControl`, Zhou et al. 2017) and writes it to `qc/conversion_qc.csv` — a value near 1.0 means complete conversion, higher means more residual incomplete conversion. GCT auto-fetches its extension probes for EPIC/450k only; EPICv2 currently records `NA` with a note (Phase 2). The GCT score also feeds the QC step: samples with GCT above `config$qc$max_gct_score` (default 1.3) **fail** QC (NA scores never fail). There's an LSF-awareness hack in this module: it detects `LSB_JOBID` and forces `BiocParallel::SerialParam()` to avoid fork-based CPU-time-limit crashes on cluster jobs.

**QC (`pipeline_modules/qc.R`).** Computes per-sample detection p-values (`minfi::detectionP`), mean detection p, percent failed probes, and pre-normalization median methylated/unmethylated channel intensities (`minfi::getQC`). A sample fails QC when mean detection p crosses `sample_detection_p_threshold` (default 0.05), failed-probe rate crosses `failed_probe_percent_threshold` (default 25), or the GCT bisulfite-conversion score crosses `max_gct_score` (default 1.3, computed in the preprocess step; NA scores never fail). Low intensity is intentionally **not** a failure — it's treated as informational because scanner gain settings vary legitimately across sites. For low-intensity samples the module additionally runs SWAN normalization and reports whether the channels recover above threshold, so you can distinguish scanner gain artifacts from true failures. It also adds three informational sample-integrity inferences from the beta matrix (sesame): predicted sex (`Sesame_Sex`, sesame `inferSex`), epigenetic age (`Horvath_Age`, the Horvath 353-CpG clock via sesame `predictAge()` using the vendored model in `Anno/HM450/Clock_Horvath353.rds`), and leukocyte fraction (`Leukocyte_Fraction`, sesame `estimateLeukocyte`; EPICv2 betas are first converted to EPIC space via the Zhou Lab `Anno/EPICv2/EPICv2ToEPIC_map.tsv.gz` map, so it works across 450k/EPIC/EPICv2). These help catch sample swaps and gauge tumour purity; none affects `Pass_QC`. (sesame 1.30.0 provides no karyotype or ethnicity inference; its bundled `age.inference` data is incompatible with the installed `predictAge()`, which is why the clock model is vendored in `Anno/`.) Outputs: `sample_qc_report.csv` with `Pass_QC`, `Failure_Reason`, `Notes`, and SWAN columns; the standalone `conversion_qc.csv` GCT bisulfite-conversion table (written by the preprocess step, informational only); plus mean-detection-p bar, beta density, beanplot, MDS PDFs, the minfi canonical `qcReport.pdf`, and plotly HTML interactive density + 3D MDS widgets.

**Filtering (`pipeline_modules/filtering.R`).** Takes the beta matrix and applies two reduction steps. First, a detection-p filter: probes with detection p > 0.05 in more than `min_sample_success_rate` of samples are dropped. Second, an array-specific curated keep-list — this is what the `data/keep.probes.{450K,EPIC,EPICv2}.txt` files are for. For EPICv2 the match is done on the probe-name prefix (strips the `_*` suffix) because EPICv2 names carry a replicate suffix that the keep list doesn't. The CLI flags `remove_sex_chromosomes`, `remove_snps`, `remove_cross_reactive` are wired through the function signature, but the body here mostly defers cross-reactive/SNP removal to the curated keep list rather than calling `DMRcate::rmSNPandCH` live (the `rmSNPandCH` call is present but commented out). Output: `filtered_beta_values.txt` plus a `filtered_probes/filtering_summary.csv`.

**Dimensionality reduction (`pipeline_modules/dim_reduction.R`).** Selects the top N most variable probes by SD/MAD (default 10,000), then runs three techniques on the samples: t-SNE via `Rtsne` (auto-caps perplexity at `(n_samples-1)/3`, drops duplicate samples and NA-containing probes), UMAP via the `umap` package, and hierarchical clustering with a Pearson or Spearman correlation distance (`1 - cor`) and configurable linkage (default `complete`). Each emits coordinates + an `.RData` object + a PDF. If a `Sample_Group` column is present in `sample_info`, points/leaves are colored by it.

**Reference projection (`pipeline_modules/reference_projection.R`).** Self-contained — it reads IDATs and SWAN-normalizes the query itself (no preprocessing step required), so it also runs standalone via `--step reference_projection`. It projects each query sample onto a pre-built reference t-SNE embedding via `snifter::project()` (an R wrapper around `openTSNE`), harmonizing the query to the reference probe set on the query side only (EPICv2 probe-ID suffixes stripped, missing probes imputed with the reference per-probe mean). Three reference datasets ship in `reference/`: `COMET_1915` (default), `Capper_GSE90496` (CNS tumours), and `Sarcoma_GSE140686`; pick one with `config$reference_projection$dataset`. A k-NN vote in the embedding space (`config$reference_projection$knn_k`, default 25) assigns each sample a nearest reference tumour class with a confidence score and `ambiguous` / `distant_from_reference` flags. Outputs go to `reference_projection/`: projected coordinates CSV, class-hints CSV, and a plot PDF. Defaults: `enabled = TRUE`, `dataset = COMET_1915`, `perplexity = 5`, `knn_k = 25`.

**CNV (`pipeline_modules/cnv_analysis.R`).** Dispatches on `method` — `conumee` (using `conumee2`) or `ChAMP`. The conumee path is the primary one: per sample, it builds a `CNV.data` object, normalizes against either user-supplied reference IDATs (via `config$reference_samples`) or internal controls shipped with yamapData, runs segmentation, and writes per-sample CNV plots + `.seg` tables into `cnv/segments/`. Across samples, `generate_cnv_frequency_plot()` produces a genome-wide gain/loss frequency plot at threshold 0.18.

**Visualization (`pipeline_modules/visualization.R`).** Pulls `qc_results.RData`, `dim_reduction_results.RData`, `cnv_results.RData` from disk, builds an R Markdown template on the fly, and renders `methylation_analysis_report.html` via `rmarkdown::render()`. It gracefully degrades: if pandoc, `rmarkdown`, `DT`, or `knitr` is missing, it logs which one and falls back to writing `methylation_analysis_report.txt`. The Rmd embeds the per-module PDFs/HTMLs from `figures/` by relative path.

## The orthogonal pieces

**`config.R`.** Supplies `default_config()` (every tunable parameter with its default), `load_config()` (accepts `.R` or `.yaml`), `setup_directories()`, and `log_message()`. Nothing in the driver requires a config file — CLI flags plus `default_config()` are sufficient.

**`utils.R`.** Dependency management: `load_package`, `load_bioc_package`, `install_source_package`, and `install_dependencies()` — the last lists the full CRAN + Bioconductor dependency set and installs yamapData from `data/yamapData_0.0.3.tar.gz` (the 257 MB tarball not currently present in this folder). If you run `install_dependencies()` from a fresh R install, this is where it expects the tarball to live.

**`hpc.R`.** Not invoked during normal runs. When `methylation_pipeline.R` is called with `--hpc`, it generates LSF / SLURM / PBS submission scripts with per-step wrappers (queue, memory, time, threads) and exits without running analysis. `command_line.txt` shows the hand-written `bsub` incantations actually used in production at St. Jude (queue `rhel88_gpu`, 500 GB memory, LSF project `OrrLab`), suggesting the HPC generator is supplementary to these bespoke commands rather than the primary path.

## Folder contents

```
pipeline/
├── methylation_pipeline.R        main CLI driver (sources all modules)
├── cnv_heatmap.R                 standalone CNV heatmap tool
├── command_line.txt              example production bsub/Rscript invocations
├── .gitignore                    R-focused ignores (results*/, .Rhistory, ._*, etc.)
├── pipeline_modules/
│   ├── config.R                  default_config(), setup_directories(), log_message()
│   ├── utils.R                   package install helpers, yamapData bootstrap
│   ├── preprocess.R              read_sample_sheet, preprocess_methylation
│   ├── qc.R                      perform_qc (detection p, intensity, SWAN recovery)
│   ├── filtering.R               filter_probes, select_variable_probes, keep-list lookup
│   ├── dim_reduction.R           run_tsne, run_umap, run_hierarchical_clustering
│   ├── reference_projection.R    run_reference_projection (snifter/openTSNE projection + k-NN class hints)
│   ├── cnv_analysis.R            run_cnv_analysis (conumee2 / ChAMP backends)
│   ├── visualization.R           generate_report (rmarkdown HTML + text fallback)
│   └── hpc.R                     generate_hpc_scripts (LSF / SLURM / PBS)
└── data/
    ├── keep.probes.450K.txt      (7.4 MB)  curated 450k probe keep-list
    ├── keep.probes.EPIC.txt      (13.4 MB) curated EPIC  probe keep-list
    ├── keep.probes.EPICv2.txt    (12.3 MB) curated EPICv2 probe keep-list
    └── example/                  8-sample EPIC demo dataset (GSE130295)
        ├── samplesheet_epic.csv  samplesheet with relative Basenames + Gender/path metadata
        ├── GSM3735546_201465940014_R01C01_{Red,Grn}.idat
        ├── GSM3735547_201465940014_R02C01_{Red,Grn}.idat
        ├── GSM3735548_201465940014_R03C01_{Red,Grn}.idat
        ├── GSM3735549_201465940014_R04C01_{Red,Grn}.idat
        ├── GSM3735550_201465940014_R05C01_{Red,Grn}.idat
        ├── GSM3735551_201465940014_R06C01_{Red,Grn}.idat
        ├── GSM3735552_201465940014_R07C01_{Red,Grn}.idat
        └── GSM3735553_201465940014_R08C01_{Red,Grn}.idat
```

## Runtime dependencies

The declared dependency set is heavy.

**CRAN:** `optparse, data.table, ggplot2, plotly, Rtsne, umap, dendextend, circlize, htmlwidgets, rmarkdown, knitr, DT, parallel, yaml, ggrepel`.

**Bioconductor:** `minfi, limma, missMethyl, RColorBrewer, matrixStats, snifter, DMRcate, conumee2, GenomicRanges, IlluminaHumanMethylation450kanno.ilmn12.hg19, IlluminaHumanMethylationEPICanno.ilm10b4.hg19, IlluminaHumanMethylationEPICv2manifest, IlluminaHumanMethylationEPICv2anno.20a1.hg38, sesame, Gviz`. `conumee2` is provisioned by the top-level `setup.R` (Bioc first, GitHub fallback) — `install_dependencies()` here does not include it.

**Local source tarball:** `yamapData` from `data/yamapData_0.0.3.tar.gz` (not currently copied into this folder — `conumee2`'s internal reference panel depends on it).

**External tools:** `pandoc` is needed for HTML report rendering; without it the pipeline falls back to a plain-text report.

## Caveats worth knowing before productizing

A few things noticed while reading that the MVP effort should address.

The filtering module calls `rmSNPandCH` only via commented-out lines — the flags `remove_snps` and `remove_cross_reactive` in the API don't actually branch any behavior; all "real" SNP/cross-reactive removal is baked into the curated keep-list files, so those CLI options are effectively no-ops today.

`select_variable_probes` is defined in both `filtering.R` and `dim_reduction.R` — whichever is sourced last wins; it happens to be `dim_reduction.R`, and the two implementations differ (the `filtering.R` version also supports `iqr`).

The 257 MB `yamapData_0.0.3.tar.gz` was intentionally not copied from the source repo and still needs to land in `data/` before the pipeline can run end-to-end, since `conumee2`'s internal reference panel depends on it.
