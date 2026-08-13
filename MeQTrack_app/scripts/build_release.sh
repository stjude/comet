#!/usr/bin/env bash
# scripts/build_release.sh — Build a distributable zip of MeQTrack_app.
#
# What goes in:
#   app/            full Shiny app
#   pipeline/       pipeline R code, example data, yamapData_*.tar.gz
#   reference/      COMET + Capper + GSE140686 sarcoma reference embeddings,
#                   metadata, and beta .rds files (the multi-hundred-MB
#                   source CSVs are excluded)
#   renv.lock       package lockfile
#   renv/           only activate.R, settings.json, .gitignore (NOT library/)
#   .Rprofile       activates renv on R startup
#   setup.R         first-launch package provisioning
#   meqtrack.command, meqtrack.bat   double-click launchers
#   QUICKSTART.md   end-user install guide
#
# What stays out:
#   renv/library/   machine-specific package cache (built on first launch)
#   runs/           local pipeline output
#   .git/, .DS_Store, *.Rhistory, .vscode/, .idea/, etc.
#
# Output: dist/MeQTrack_app-<git-sha>-<date>.zip
#
# Usage:
#   bash scripts/build_release.sh
#
# Run from any directory; the script resolves paths relative to itself.

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$( cd "${SCRIPT_DIR}/.." &> /dev/null && pwd )"

cd "${PROJECT_ROOT}"

echo "=============================================================="
echo "  MeQTrack release builder"
echo "  Project: ${PROJECT_ROOT}"
echo "=============================================================="

# ---------------------------------------------------------------------------
# 1. Preflight
# ---------------------------------------------------------------------------
yamap_tarball="pipeline/data/yamapData_0.0.3.tar.gz"
if [ ! -f "${yamap_tarball}" ]; then
  cat >&2 <<EOS
ERROR: ${yamap_tarball} is missing.

The yamapData tarball is gitignored because of its size (~257 MB), but
the release zip MUST include it — without it the CNV step has no
internal reference and setup.R will refuse to run.

Place the tarball at ${yamap_tarball} and re-run this script.
EOS
  exit 1
fi

# Reference beta matrices are gitignored (>100 MB each), but the release
# zip MUST include every registered dataset's .rds — the
# reference-projection step has nothing to project against otherwise.
# Each is the compact .rds that reference_projection.R::load_reference()
# reads. Add a line here when a new dataset joins .REFERENCE_DATASETS.
ref_betas=(
  "reference/beta_GSE305405_1915samples_top10K.rds"
  "reference/beta_GSE90496_top10K.rds"
  "reference/beta_GSE140686_1077Sarcoma_top10K.rds"
)
for ref_beta in "${ref_betas[@]}"; do
  if [ ! -f "${ref_beta}" ]; then
    cat >&2 <<EOS
ERROR: ${ref_beta} is missing.

Reference beta matrices are gitignored (>100 MB each), but the release
zip MUST include each registered dataset's .rds — the
reference-projection step has nothing to project against without it.

Place the .rds at ${ref_beta} and re-run this script.
EOS
    exit 1
  fi
done

if ! command -v rsync >/dev/null 2>&1; then
  echo "ERROR: rsync is required (preinstalled on macOS / Linux)." >&2
  exit 1
fi
if ! command -v zip >/dev/null 2>&1; then
  echo "ERROR: zip is required (preinstalled on macOS / Linux)." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 2. Compute release name (git-sha + date if in a repo, else just date).
# ---------------------------------------------------------------------------
release_date="$(date +%Y%m%d)"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  git_sha="$(git rev-parse --short HEAD)"
  release_tag="${git_sha}-${release_date}"
else
  release_tag="${release_date}"
fi
release_name="MeQTrack_app-${release_tag}"

dist_dir="${PROJECT_ROOT}/dist"
stage_dir="${dist_dir}/${release_name}"
zip_path="${dist_dir}/${release_name}.zip"

mkdir -p "${dist_dir}"
rm -rf "${stage_dir}" "${zip_path}"
mkdir -p "${stage_dir}"

# ---------------------------------------------------------------------------
# 3. Stage files via rsync. Each entry mirrored into stage_dir, with the
#    big-pile excludes applied to renv/.
# ---------------------------------------------------------------------------
echo
echo "Staging into ${stage_dir} ..."

rsync -a --exclude=".DS_Store" --exclude="*.Rhistory" \
  app/ "${stage_dir}/app/"

rsync -a --exclude=".DS_Store" --exclude="*.Rhistory" \
  pipeline/ "${stage_dir}/pipeline/"

# reference/: ship the embeddings, metadata, and the compact beta .rds.
# The multi-hundred-MB source CSVs (beta_*.csv) stay out — load_reference
# reads the .rds.
rsync -a --exclude=".DS_Store" --exclude="beta_*.csv" \
  reference/ "${stage_dir}/reference/"

# Anno/: small vendored sesame annotation assets (Horvath clock models). The
# QC step reads load_horvath_model() from here, so the zip MUST include it.
rsync -a --exclude=".DS_Store" \
  Anno/ "${stage_dir}/Anno/"

# renv/: ship activator, settings, and .gitignore only — never library/.
mkdir -p "${stage_dir}/renv"
for f in activate.R settings.json .gitignore; do
  if [ -f "renv/${f}" ]; then
    cp "renv/${f}" "${stage_dir}/renv/${f}"
  fi
done

cp renv.lock          "${stage_dir}/renv.lock"
cp .Rprofile          "${stage_dir}/.Rprofile"
cp setup.R            "${stage_dir}/setup.R"
cp meqtrack.command   "${stage_dir}/meqtrack.command"
cp meqtrack.bat       "${stage_dir}/meqtrack.bat"

# End-user documentation:
#   README.md                       — orientation: what MeQTrack is + workflow
#   QUICKSTART.md                   — install + run instructions
# Internal planning/requirements docs and docs/decisions/ stay out of the
# release — they're developer notes, not user-facing. The planning docs
# live in a separate private repo.
cp README.md                       "${stage_dir}/README.md"
cp QUICKSTART.md                   "${stage_dir}/QUICKSTART.md"

# Make sure the macOS launcher is executable inside the zip.
chmod +x "${stage_dir}/meqtrack.command"

# ---------------------------------------------------------------------------
# 4. Sanity-check key files are present in the staged tree.
# ---------------------------------------------------------------------------
required=(
  "app/app.R"
  "pipeline/methylation_pipeline.R"
  "pipeline/data/yamapData_0.0.3.tar.gz"
  "pipeline/data/keep.probes.EPIC.txt"
  "reference/beta_GSE305405_1915samples_top10K.rds"
  "reference/tSNE_embedding_GSE305405_top10K.RData"
  "reference/COMET_Labkey_August_12_2025.csv"
  "Anno/HM450/Clock_Horvath353.rds"
  "Anno/EPICv2/EPICv2.typeI.ext.rds"
  "reference/beta_GSE90496_top10K.rds"
  "reference/tSNE_embedding_GSE90496_top10K.RData"
  "reference/GSE90496_MC_MCF_color_labels_key.csv"
  "reference/beta_GSE140686_1077Sarcoma_top10K.rds"
  "reference/tSNE_embedding_GSE140686_top10K.RData"
  "reference/GSE140686_sarcoma_methylation_labels.csv"
  "renv.lock"
  "renv/activate.R"
  ".Rprofile"
  "setup.R"
  "meqtrack.command"
  "meqtrack.bat"
  "README.md"
  "QUICKSTART.md"
)
missing=()
for path in "${required[@]}"; do
  if [ ! -e "${stage_dir}/${path}" ]; then
    missing+=("${path}")
  fi
done
if [ "${#missing[@]}" -gt 0 ]; then
  echo "ERROR: required files missing from staged tree:" >&2
  printf '  %s\n' "${missing[@]}" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 5. Zip and report.
# ---------------------------------------------------------------------------
echo
echo "Zipping to ${zip_path} ..."
( cd "${dist_dir}" && zip -qr "${release_name}.zip" "${release_name}" )

# Tear down the staging dir; the zip is the artifact we keep.
rm -rf "${stage_dir}"

zip_size_h="$(du -h "${zip_path}" | awk '{print $1}')"
file_count="$(unzip -l "${zip_path}" | tail -1 | awk '{print $2}')"

cat <<EOS

==============================================================
  Release built.

  Path: ${zip_path}
  Size: ${zip_size_h}
  Files: ${file_count}

  Test instructions for a clean machine:
    1. unzip ${release_name}.zip
    2. cd ${release_name}
    3. macOS:   double-click meqtrack.command
       Windows: double-click meqtrack.bat
==============================================================
EOS
