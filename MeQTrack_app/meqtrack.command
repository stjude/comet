#!/usr/bin/env bash
# meqtrack.command — macOS launcher for MeQTrack.
#
# Double-click this file in Finder to start the app. Terminal.app opens,
# this script runs, Shiny starts, and your default browser opens to the
# app on http://127.0.0.1:<port>.
#
# The terminal window stays open while the app is running so you can see
# server-level logs. Close it to stop the app.

set -e

# Resolve the script's own directory (robust to spaces, symlinks, and
# being run from anywhere). This lets the user put MeQTrack_app wherever
# they like — the launcher always cds to its sibling app/ and pipeline/.
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "${SCRIPT_DIR}"

# Silence renv's "project is out-of-sync" startup notice. setup.R installs
# current package versions and re-snapshots, but renv can't fully pin the
# GitHub (conumee2) / local-tarball (yamapData) packages, so it always flags
# a residual lock-vs-library difference. It's cosmetic — every package the
# app needs is installed — so we suppress the check in the launcher to avoid
# alarming end users.
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

echo "=============================================================="
echo "  MeQTrack — launching..."
echo "  Project: ${SCRIPT_DIR}"
echo "=============================================================="

# 1. Check R is on PATH. We don't try to install R ourselves.
if ! command -v Rscript >/dev/null 2>&1; then
  cat <<'EOS'

  ERROR: Rscript not found on your PATH.

  MeQTrack needs R >= 4.4 installed. Install from:
      https://cran.r-project.org/bin/macosx/
  and then re-run this launcher.

EOS
  echo "Press Return to close this window."
  read -r _
  exit 1
fi

# 2. Check pandoc (needed by the pipeline's report generator).
if ! command -v pandoc >/dev/null 2>&1; then
  cat <<'EOS'

  WARNING: pandoc not found on your PATH.

  The HTML report generator falls back to a plain-text report when pandoc
  is missing. Install via Homebrew to enable the full report:
      brew install pandoc

  Continuing anyway...

EOS
fi

# 3. First-launch provisioning. The distribution zip ships renv.lock and
#    renv/activate.R but NOT the renv/library/ cache (machine-specific
#    and would balloon the zip). On first launch we need to populate the
#    library — setup.R handles renv::restore + Bioconductor + yamapData
#    and is safe to re-run (no-op when already provisioned).
#
#    We probe by trying to load `shiny`, NOT by checking if renv/library/
#    exists. R keys its library by minor version (R-4.5 vs R-4.6 are
#    separate subdirs), so a user who upgrades R between launches has a
#    populated renv/library/ that still fails to load anything against
#    the new R. requireNamespace catches that case; setup.R then
#    populates the library against the current R.
NEED_SETUP=0
if ! Rscript -e 'quit(status = !requireNamespace("shiny", quietly = TRUE))' >/dev/null 2>&1; then
  NEED_SETUP=1
fi
if [ "$NEED_SETUP" = "1" ]; then
  echo
  echo "  First-time setup: installing R packages."
  echo "  This takes 5-15 minutes on first run; subsequent launches are instant."
  echo "  Watch this window for progress; the app starts automatically when done."
  echo
  if ! Rscript setup.R; then
    echo
    echo "  ERROR: setup.R failed. See the messages above."
    echo "  Press Return to close this window."
    read -r _
    exit 1
  fi
  echo
  echo "  Setup complete. Starting the app..."
  echo
fi

# 4. Start Shiny. runApp reads app/app.R; renv activates automatically via
#    the project's .Rprofile so the installed library is used.
echo "Starting Shiny server on http://127.0.0.1 ..."
echo "(Close this Terminal window to stop the app.)"
echo

exec Rscript -e 'shiny::runApp("app", launch.browser = TRUE, host = "127.0.0.1")'
