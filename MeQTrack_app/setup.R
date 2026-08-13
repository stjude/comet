# setup.R — one-time provisioning for MeQTrack_app.
#
# Run from the project root on the user's Mac, with R >= 4.4:
#   Rscript setup.R
#
# What it does:
#   1. Installs `renv` globally if missing.
#   2. Activates renv in this project (creates .Rprofile + renv/).
#   3. Installs BiocManager pinned to a Bioconductor release that matches R's minor version.
#   4. Installs CRAN + Bioconductor dependencies used by the pipeline.
#   5. Installs the bundled yamapData tarball from pipeline/data/.
#   6. Snapshots everything into renv.lock.
#
# Idempotent: safe to re-run. Packages already installed are skipped by renv.

# Use Posit Package Manager (PPM) as the default CRAN mirror on macOS.
# PPM hosts macOS arm64 binaries for current R versions more reliably than
# CRAN itself — including freshly-released package versions where CRAN
# hasn't built the binary yet. This is the main defence against the
# `<cmath> file not found` source-compile failure on macOS Tahoe SDKs.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))

# Force binary-only installs on macOS. Recent macOS SDKs (Tahoe / MacOSX26.x)
# break source compilation for R packages that use C++ — clang can't find
# <cmath> because libc++ headers moved within the SDK. CRAN and Bioconductor
# ship arm64 binaries for the supported R versions (4.4 / 4.5 / 4.6), so we
# never need to compile from source on this host. If a package has no binary
# available, fail fast rather than trying to compile — easier to see and fix.
if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Darwin") {
  options(pkgType = "binary")
  options(install.packages.compile.from.source = "never")
}

# Package type for the explicit install.packages() calls below. macOS pulls
# arm64 binaries from PPM (above); Linux/HPC builds from SOURCE. PPM's generic
# /cran/latest URL serves source on Linux, and forcing type = "binary" there
# yields broken/half-installed packages (e.g. askpass with no Meta/package.rds).
# HPC nodes have compilers, so source builds are the reliable path.
.pkg_type <- if (Sys.info()[["sysname"]] == "Darwin") "binary" else "source"

# Use all available cores for any compilation that does happen.
options(Ncpus = parallel::detectCores(logical = FALSE))

# ---------------------------------------------------------------------------
# 0. Sanity checks
# ---------------------------------------------------------------------------
r_ver <- getRversion()
message(sprintf("R version detected: %s", r_ver))

if (r_ver < "4.4.0") {
  stop(
    "MeQTrack_app targets R >= 4.4. Please upgrade R before running setup.R.\n",
    "  See https://cran.r-project.org/"
  )
}

# Let BiocManager pick the Bioconductor release for the running R. The
# pairing isn't 1:1 (Bioc has two releases per year against the same R
# minor — e.g. R 4.5 had both Bioc 3.21 and 3.22) so any hardcoded map
# we ship would go stale every six months. Passing version = NULL to
# BiocManager::install() means "use the current release for this R",
# which is what we want.
message(sprintf(
  "Will let BiocManager pick the Bioconductor release for R %s.", r_ver
))

# Verify we're running from the project root (pipeline/ must exist).
if (!dir.exists("pipeline") || !file.exists("pipeline/methylation_pipeline.R")) {
  stop(
    "setup.R must be run from the MeQTrack_app project root.\n",
    "  Current wd: ", getwd()
  )
}

# Verify yamapData tarball is present.
yamap_tarball <- file.path("pipeline", "data", "yamapData_0.0.3.tar.gz")
if (!file.exists(yamap_tarball)) {
  stop(
    "Missing yamapData tarball at ", yamap_tarball, ".\n",
    "  The CNV step cannot run without it. Place it in pipeline/data/ and re-run."
  )
}

# ---------------------------------------------------------------------------
# 1. renv bootstrap
# ---------------------------------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv...")
  install.packages("renv")
}

# Activate renv in the project (idempotent; creates .Rprofile + renv/ if missing).
if (!file.exists("renv/activate.R")) {
  message("Activating renv in the project...")
  renv::activate()
} else {
  source("renv/activate.R")
}

# Tell renv's implicit dependency scanner to ignore the optional CNV
# backends that the pipeline references via `::` but doesn't run in the
# default path. Without this, every snapshot will re-add ChAMP + ChAMPdata
# to the lockfile and demand they be installed.
renv::settings$ignored.packages(c("ChAMP", "ChAMPdata"))

# ---------------------------------------------------------------------------
# 2. BiocManager, pinned
# ---------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  renv::install("BiocManager")
}

# Initialize Bioconductor against whatever release pairs with the
# current R. No version arg means "current" — which is the right
# choice because Bioc's release schedule (twice a year, paired with R)
# doesn't fit a hardcoded map.
BiocManager::install(update = FALSE, ask = FALSE)

# ---------------------------------------------------------------------------
# 3. CRAN dependencies
# ---------------------------------------------------------------------------
cran_packages <- c(
  # Pipeline runtime
  "optparse", "data.table", "ggplot2", "plotly", "Rtsne", "umap",
  "dendextend", "circlize", "htmlwidgets", "rmarkdown", "knitr",
  "DT", "yaml", "ggrepel", "RColorBrewer",
  # Manages a pandoc binary for the HTML report on hosts without a system
  # pandoc (CRAN R ships none). Tiny package; the binary is fetched in
  # section 6b only when needed. See visualization.R for the render-time wiring.
  "pandoc",
  # App runtime
  "shiny", "bslib", "shinyFiles", "shinycssloaders",
  "promises", "future", "callr",
  "ggdendro", "jsonlite",
  # Needed for the conumee2 GitHub install (supports `subdir = ...`):
  #   remotes::install_github("hovestadtlab/conumee2", subdir = "conumee2")
  # remotes (not devtools) is the actual install_github workhorse: it has
  # almost no dependencies and a binary is always available. devtools merely
  # forwards to remotes and now loads it lazily via check_installed(), so
  # installing devtools alone leaves remotes missing and install_github fails
  # with: 'The package "remotes" is required.'
  "remotes"
)

message("Installing CRAN packages (may take a while on first run)...")
# Use install.packages directly instead of renv::install here. renv::install
# can silently no-op on packages that are already referenced in the lockfile
# but not actually installed (leaving the library incomplete). Going through
# install.packages() + the binary-only options above is more reliable on
# macOS Tahoe hosts.
#
# repos is a 2-element fallback chain: PPM first (faster, hosts the freshest
# arm64 binaries), then CRAN mainline. PPM occasionally lacks a binary for
# the newest R minor (e.g. devtools missing for R 4.6 / arm64-darwin) — in
# that case install.packages falls through to CRAN mainline rather than
# failing the whole setup.
.cran_repos <- c(
  PPM  = "https://packagemanager.posit.co/cran/latest",
  CRAN = "https://cran.r-project.org"
)
install.packages(cran_packages, type = .pkg_type, repos = .cran_repos)

# ---------------------------------------------------------------------------
# 4. Bioconductor dependencies
# ---------------------------------------------------------------------------
# NOTE: conumee2 is not in every Bioc release, so we handle it separately
# below with a GitHub fallback rather than including it in this list.
bioc_packages <- c(
  # core analysis
  "minfi", "limma", "missMethyl", "matrixStats", "snifter",
  "DMRcate", "GenomicRanges", "sesame", "Gviz",
  # conumee2 transitive deps. conumee2 is installed from GitHub below
  # and its source build needs methylumi to load (methylumi in turn
  # requires FDb.InfiniumMethylation.hg19). remotes::install_github
  # doesn't reliably resolve these via BiocManager — easier to install
  # them explicitly here so they're already in the renv library when
  # conumee2 is built.
  "methylumi", "FDb.InfiniumMethylation.hg19",
  # Bioc experiment-data companions. These are nominally pulled in as
  # dependencies of their parent packages, but BiocManager has been
  # observed to silently skip them on macOS (no binary, source build
  # fails to fetch the data payload, etc.) and the parent then fails to
  # load at runtime with "package X required by Y could not be found".
  # Listing them explicitly forces the install.
  "sesameData",
  # Illumina 450K / EPIC / EPICv2 manifests + annotations
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2manifest",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
)

# NOTE on ChAMP / ChAMPdata:
# The pipeline supports three CNV backends: conumee (default), ChAMP, and
# cnAnalysis450k. Only conumee runs in the default pipeline path and in the
# MVP example. ChAMP is optional — selected by setting
# config$cnv$method = "ChAMP". If a user wants that backend, they install
# ChAMP + ChAMPdata themselves:
#   BiocManager::install(c("ChAMP", "ChAMPdata"),
#                        update = FALSE, ask = FALSE)
# We skip it here to avoid a ~40 MB download from the Bioc 3.21 archive
# mirror (slow, frequently times out) for code that never runs in MVP.

message("Installing Bioconductor packages (this is the slow part)...")
BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)

# Data-only Bioc packages (Illumina manifests/annos, sesameData, FDb...)
# don't always ship a macOS binary on the Bioc mirror — particularly on a
# freshly-released Bioc minor (e.g. 3.23) where binaries lag the source
# release by days/weeks. The bulk install above silently skips them
# because of our global `install.packages.compile.from.source = "never"`
# setting (the macOS Tahoe SDK defense for *compiled* packages). For
# pure-R/data packages there's nothing to compile, so source install is
# safe and fast. Any missing one would otherwise cascade — e.g. missMethyl
# fails to load when its 450k anno dep didn't land.
data_only_pkgs <- c(
  # conumee2 transitive (methylumi → FDb)
  "FDb.InfiniumMethylation.hg19",
  # Illumina manifests + annotations consumed by minfi/missMethyl/sesame
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2manifest",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  # sesame's ExperimentHub-fed data companion
  "sesameData"
)
data_only_missing <- data_only_pkgs[
  !vapply(data_only_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(data_only_missing) > 0) {
  message("Installing data-only Bioc packages from source: ",
          paste(data_only_missing, collapse = ", "))
  BiocManager::install(data_only_missing, update = FALSE, ask = FALSE,
                       type = "source")
}

# conumee2: not in every Bioc release. Try Bioc first, fall back to GitHub.
# The GitHub repo hosts the package inside a `conumee2/` subdirectory, so
# we use remotes::install_github(subdir=) rather than renv::install(), which
# would look for DESCRIPTION at the repo root and fail.
if (!requireNamespace("conumee2", quietly = TRUE)) {
  message("Installing conumee2 (Bioc first, GitHub fallback)...")
  bioc_ok <- tryCatch({
    BiocManager::install("conumee2", update = FALSE, ask = FALSE)
    requireNamespace("conumee2", quietly = TRUE)
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (!isTRUE(bioc_ok)) {
    bioc_label <- tryCatch(
      as.character(BiocManager::version()),
      error = function(e) "(unknown)"
    )
    message("conumee2 not in Bioc ", bioc_label,
            "; installing from GitHub (hovestadtlab/conumee2, subdir=conumee2)...")
    # remotes is in the CRAN list above; this is a belt-and-suspenders
    # install in case the CRAN step skipped it for any reason. Same
    # PPM-then-CRAN fallback chain as the bulk install.
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", type = .pkg_type, repos = .cran_repos)
    }
    remotes::install_github("hovestadtlab/conumee2",
                            subdir = "conumee2",
                            upgrade = "never")
    # remotes::install_github can warn and return without actually installing
    # (e.g. when the host has no git, the GitHub API rate-limits, or the
    # tarball download fails partway). Verify the install landed; the final
    # post-install check at the end of setup.R catches this too, but failing
    # here gives a more localised error.
    if (!requireNamespace("conumee2", quietly = TRUE)) {
      stop(
        "conumee2 GitHub fallback did not deliver the package to the renv ",
        "library. Check the messages above for the underlying remotes/git ",
        "error (no network, missing git, GitHub rate limit, etc.) and re-run ",
        "setup.R after fixing it."
      )
    }
  }
}

# ---------------------------------------------------------------------------
# 5. yamapData from local tarball
# ---------------------------------------------------------------------------
if (!requireNamespace("yamapData", quietly = TRUE)) {
  message(sprintf("Installing yamapData from %s ...", yamap_tarball))
  install.packages(yamap_tarball, repos = NULL, type = "source")
  if (!requireNamespace("yamapData", quietly = TRUE)) {
    stop("yamapData installation failed.")
  }
} else {
  message("yamapData already installed.")
}

# ---------------------------------------------------------------------------
# 6. Verify every required package can actually be loaded.
# ---------------------------------------------------------------------------
# Several install paths above can succeed silently without delivering the
# package into the renv project library — most notably the conumee2 GitHub
# fallback (remotes::install_github can warn-and-bail without raising an
# error) and any Bioc install where the requested package isn't in the
# current Bioc release for the running R. If we let setup.R "complete" in
# that state, the user only finds out later when the pipeline crashes at
# `library(...)` time. Probe every package the pipeline attaches and stop
# loudly here instead.
required_pkgs <- c(
  # CRAN
  "data.table", "ggplot2", "plotly", "Rtsne", "umap", "dendextend",
  "stringr", "shiny", "bslib", "shinyFiles", "shinycssloaders",
  "promises", "future", "callr", "ggdendro", "jsonlite", "yaml",
  "DT", "ggrepel", "RColorBrewer", "rmarkdown", "knitr",
  # Bioc core
  "minfi", "limma", "missMethyl", "matrixStats", "snifter",
  "DMRcate", "GenomicRanges", "sesame", "sesameData", "Gviz",
  "methylumi", "FDb.InfiniumMethylation.hg19",
  "conumee2",
  # Bioc Illumina manifests + annotations. Listed in bioc_packages above
  # but BiocManager has been observed to silently skip individual data
  # packages (e.g. EPICv2anno.20a1.hg38) when their binary isn't
  # available; missMethyl then fails at runtime with
  # "package <X> required by missMethyl could not be found".
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2manifest",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  # Local tarball
  "yamapData"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "setup.R finished its install steps but the following packages are still ",
    "not loadable from the renv project library:\n  ",
    paste(missing_pkgs, collapse = ", "),
    "\n\nThis usually means an install warned-and-bailed without raising an ",
    "error (e.g. conumee2 from GitHub on a host without git/network access). ",
    "Re-run setup.R after fixing the underlying issue; the messages above ",
    "should point at the failing package."
  )
}

# ---------------------------------------------------------------------------
# 6b. Pandoc for the HTML report
# ---------------------------------------------------------------------------
# rmarkdown::render() (the report step) needs a pandoc binary. R installed
# from CRAN does not ship one — only RStudio bundles pandoc — so on a plain
# double-click install pandoc is absent and the report silently falls back
# to plain text. If no system pandoc is found, provision a managed copy via
# the `pandoc` package; visualization.R points RSTUDIO_PANDOC at it at render
# time. Non-fatal: the report degrades gracefully, so a failed pandoc
# download must not abort the whole setup.
if (!nzchar(rmarkdown::find_pandoc()$dir)) {
  message("No system pandoc found; provisioning a managed copy for the HTML report...")
  if (!requireNamespace("pandoc", quietly = TRUE)) {
    install.packages("pandoc", type = .pkg_type, repos = .cran_repos)
  }
  pandoc_ok <- tryCatch({
    if (!pandoc::pandoc_available()) pandoc::pandoc_install()
    isTRUE(pandoc::pandoc_available())
  }, error = function(e) {
    message("  pandoc provisioning failed: ", conditionMessage(e))
    FALSE
  })
  if (isTRUE(pandoc_ok)) {
    message("  Managed pandoc ready: ", pandoc::pandoc_bin())
  } else {
    message("  Could not provision pandoc. The HTML report will fall back to ",
            "plain text; install pandoc manually (https://pandoc.org/installing.html) ",
            "to enable it.")
  }
} else {
  message("Using system pandoc: ", rmarkdown::find_pandoc()$dir)
}

# ---------------------------------------------------------------------------
# 7. Pre-cache sesame data via ExperimentHub.
# ---------------------------------------------------------------------------
# sesame's openSesame() / pOOBAH() (called by the preprocess step) need
# data files (IDAT signatures, platform manifests, normalization tables)
# that aren't shipped with the sesame/sesameData packages — they're
# fetched lazily on first use via ExperimentHub. A fresh install therefore
# fails preprocess with `stopAndCache("idatSignature")` even though all
# packages are installed correctly. Pre-cache here so the (potentially
# multi-minute) download happens during setup, not mid-pipeline.
message("Caching sesame data via ExperimentHub (one-time, may take several minutes)...")
sesame_cache_ok <- tryCatch({
  sesameData::sesameDataCache()
  TRUE
}, error = function(e) {
  message("sesameDataCache() failed: ", conditionMessage(e))
  FALSE
})
if (!isTRUE(sesame_cache_ok)) {
  stop(
    "sesameDataCache() did not complete. The pipeline's preprocess step ",
    "needs sesame data files cached locally; without them sesame::openSesame() ",
    "fails with 'File idatSignature either not found or needs to be cached'.\n",
    "  Re-run setup.R after fixing network access to ExperimentHub, or run ",
    "`sesameData::sesameDataCache()` manually in R from the project root."
  )
}
# Sanity-check the specific data file the preprocess step needs first.
idat_sig_ok <- tryCatch({
  invisible(sesameData::sesameDataGet("idatSignature"))
  TRUE
}, error = function(e) {
  message("sesameDataGet('idatSignature') failed: ", conditionMessage(e))
  FALSE
})
if (!isTRUE(idat_sig_ok)) {
  stop(
    "sesameDataCache() reported success but 'idatSignature' is still not ",
    "retrievable. Check ExperimentHub connectivity and try again."
  )
}

# ---------------------------------------------------------------------------
# 8. Snapshot lockfile
# ---------------------------------------------------------------------------
message("Writing renv.lock...")
# type = "implicit" scans only packages actually referenced by project R
# scripts (setup.R, pipeline/*, app/*, scripts/*), so it ignores unrelated
# packages that may already live in the user's library.
# force = TRUE skips the interactive prompt when renv sees packages it
# can't classify (e.g. yamapData installed from a local tarball).
renv::snapshot(prompt = FALSE, type = "implicit", force = TRUE)

message("\nDone. Next: Rscript scripts/run_example.R")
