# MeQTrack — Quickstart

End-user install path is **unzip + double-click**. Developer steps for
running from a clone are further down.

## Install (from the release zip)

You need **R ≥ 4.4** and **pandoc ≥ 2.x** on your machine before
launching:

- **R** — <https://cran.r-project.org/>
- **pandoc** — `brew install pandoc` on macOS, or
  <https://pandoc.org/installing.html> on Windows.

Then:

1. Download `MeQTrack_app-<version>.zip` and unzip it anywhere on disk
   (Desktop, `~/Applications`, an external drive — wherever you'd like).
2. **macOS:** double-click `meqtrack.command` inside the unzipped folder.
   **Windows:** double-click `meqtrack.bat`.
3. A Terminal / Command Prompt window opens. The **first launch** runs
   `setup.R` to install all R packages — this takes 5–15 minutes
   depending on network and CPU. Watch the window; the app starts
   automatically when setup finishes.
4. Subsequent launches are instant: the launcher detects the populated
   library cache and skips straight to starting the Shiny server, then
   opens the app in your default browser at `http://127.0.0.1:<port>`.
5. Close the Terminal / Command window to stop the app. Uninstall =
   delete the unzipped folder.

If R or pandoc isn't installed, the launcher prints a clear message
pointing at the install link instead of crashing.

## Using the app

The in-app **Help** tab walks through the workflow end-to-end. The
short version:

1. **Samplesheet tab** — pick a CSV. Per-row validation runs
   immediately. The bundled example lives at
   `pipeline/data/example/samplesheet_epic.csv`.
2. **Run tab** — optionally tweak the parameters in the Settings
   card (defaults work for typical small cohorts), then click
   **Run analysis**. The Stages panel shows live per-step progress.
3. **Result tabs** (QC / Dim. reduction / CNV / Report) populate as
   each pipeline stage finishes — you don't have to wait for the whole
   pipeline to view QC. (Reference projection appears as a sub-tab under
   Dim. reduction.)
4. **Past runs tab** — every prior run is listed; click **Open** on
   any row to attach it. Result tabs render its artifacts and the
   per-step Run buttons start operating against that run.

Each Stages-panel row has a small **▶ Run** button to re-run only that
stage against the current run directory — useful for tweaking a
parameter and just redoing dimensionality reduction without re-running
preprocessing.

The **Reference projection** sub-tab (under Dim. reduction) projects
your samples onto the bundled COMET reference cohort and reports, per
sample, the nearest reference tumour group. It runs in a full analysis,
or on its own via the Reference projection step's Run button.

## Where your data lives

On first launch the app creates `~/MeQTrack/` with three subfolders:

- `samplesheets/` — drop your CSVs here for quick access.
- `idats/` — convenient location for IDAT files; the samplesheet can
  point anywhere on disk.
- `runs/<timestamp>_<samplesheet>/` — every run lands here with this
  layout:

```
runs/<id>/
├── processed_data/             β-values, RGChannelSet, sample_info, etc.
├── qc/                         qc_results.RData, sample_qc_report.csv
├── dimensionality_reduction/   tsne / umap / hclust .RData
├── cnv/                        cnv_results.RData, segment .seg files
├── reference_projection/       projected coords, class hints, plot PDF
├── figures/                    interactive HTML plots, per-sample CNV PDFs
├── reports/
│   └── methylation_analysis_report.html   ← open this in a browser
├── logs/                       pipeline.log, bridge.out, bridge.err
├── run_manifest.json           samplesheet, parameters, status, timing
└── run_config.R                only present when Settings overrides used
```

`run_manifest.json` is the source of truth for the Past runs tab and
makes every run reproducible — re-attach a past run and the Settings
inputs auto-populate from its manifest.

---

## Developer setup (from a clone)

Skip this if you used the release zip above — these steps reproduce
the same environment manually. Useful for working on the code.

### Prerequisites

You need these installed before running anything in this doc:

1. **R ≥ 4.4** — check with `R --version`. Install from
   <https://cran.r-project.org/>.
2. **pandoc ≥ 2.x** — check with `pandoc --version`. Install via
   `brew install pandoc` on macOS, or
   <https://pandoc.org/installing.html> on Windows.
3. **Xcode Command Line Tools** on macOS (for compiling R packages
   from source) — `xcode-select --install` if not already installed.

All other packages (Bioconductor, `yamapData`, `conumee2`, `sesame`,
etc.) install automatically in step 1 below.

### Step 1 — Provision the project

From the project root run:

```sh
Rscript setup.R
```

This will:

- Install `renv` if you don't have it.
- Activate `renv` in this project (creates `.Rprofile` + `renv/`).
- Install CRAN + Bioconductor packages. Bioconductor's release
  schedule (twice a year, paired with R) doesn't fit a hardcoded map,
  so setup.R lets BiocManager pick whatever release is current for
  your running R — works on R 4.4, 4.5, 4.6, and any future minor R
  updates. First run takes 5–15 min depending on network and CPU
  because several Bioc packages compile from source.
- Install `yamapData` from the bundled
  `pipeline/data/yamapData_0.0.3.tar.gz`.
- Write `renv.lock` so the environment is reproducible.

Re-running `Rscript setup.R` is safe and cheap — installed packages
are skipped.

### Step 2 — Run the example pipeline end-to-end (CLI)

```sh
Rscript scripts/run_example.R
```

This invokes `pipeline/methylation_pipeline.R` with `--step all`
against `pipeline/data/example/samplesheet_epic.csv` and writes output
to `runs/example_<timestamp>/`. Useful for proving the pipeline works
without the Shiny UI.

### Step 3 — Launch the app (dev mode)

#### Via double-click (recommended)

**macOS:** in Finder, double-click `meqtrack.command`. A Terminal
window opens, runs the server, and your default browser opens to the
app.

**Windows:** in File Explorer, double-click `meqtrack.bat`. A Command
Prompt window opens, runs the server, and your default browser opens
to the app.

In both cases: **keep the terminal window open** while using the app.
Closing it stops the server.

#### Manual launch (fallback)

From the project root:

```sh
R -e 'shiny::runApp("app", launch.browser = TRUE, host = "127.0.0.1")'
```

## Troubleshooting

**`setup.R` fails installing a Bioconductor package.**
Most common cause on macOS is a missing system library (e.g.
`libxml2`, `libssl`). Check the error message; install via Homebrew
and re-run.

**Pipeline fails at the CNV step with a `BiocParallel` error.**
Usually a fork-related issue with sesame on macOS. Re-running often
clears it. If it persists, set the env var `MEQTRACK_SERIAL=1` before
launching to force serial execution.

**The report is missing.**
Confirm `pandoc` is installed (`pandoc --version`). The pipeline falls
back to a text-only output if pandoc is absent — you'll see a warning
in the pipeline log.

**Slow first run.**
Expected. Installing Bioconductor from source on a fresh machine is
the long pole. Subsequent runs reuse the `renv` library and start
fast.

**First reference-projection run is slow.**
The projection step uses a Python toolchain (openTSNE via `snifter`)
that provisions a self-contained Python environment on first use — a
one-time download of a few minutes. Later runs reuse it.

**Double-clicking `meqtrack.command` does nothing / opens a text
editor.**
On macOS, Finder may have lost the executable bit. From a terminal in
the project root:

```sh
chmod +x meqtrack.command
```

Then try again. If macOS Gatekeeper blocks it the first time
("cannot be opened because it is from an unidentified developer"),
right-click the file → Open → Open, then subsequent double-clicks
will work.

**Browser opens but the page fails to load.**
Shiny picks a random free port on 127.0.0.1. Wait a few seconds for
the server to finish starting — the terminal window will print
`Listening on http://127.0.0.1:<port>` when ready. Reload the browser
tab once that line appears.

**Per-step Run buttons are all disabled.**
You haven't started a run in this session yet, or you don't have a
past run attached. Either click **Run analysis** (full pipeline) or
go to the Past runs tab and Open a prior run, then per-step buttons
become available against that run directory.
