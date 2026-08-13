#!/bin/bash
# =============================================================================
# MeQTrack — LSF submission script (full pipeline, one job)
#
# Sized for the 220-sample 450k Melanoma run
# (samplesheet_450k_hpc.csv, IDATs under /research/groups/orrgrp/...).
#
# Submit from the repo root on the cluster:
#     mkdir -p /research/rgs01/home/clusterHome/qtran/Melanoma/meqtrack  # once (parent for run dirs + LSF logs)
#     bsub < scripts/submit_hpc_lsf.sh
#
# Set the cluster-side REPO_DIR and SAMPLESHEET paths (marked # <-- EDIT) below,
# then submit. Project (melanoma), queue (rhel88_gpu), output (a timestamped
# dir under .../qtran/Melanoma/meqtrack/), modules (R/4.4.0, pandoc/3.8.3),
# and resources are already set.
# =============================================================================

#BSUB -J meqtrack_melanoma_450k
#BSUB -P melanoma                           # LSF project / account code
#BSUB -q rhel88_gpu                         # queue
##BSUB -gpu "num=1/host"                    # rhel88_gpu is a GPU queue but this pipeline is
                                            # CPU-only; uncomment ONLY if the queue rejects
                                            # CPU-only jobs (don't hold a GPU you won't use).
#BSUB -n 2                                  # few cores, big memory (keep == THREADS below)
#BSUB -R "span[hosts=1] rusage[mem=128GB]"  # 128GB RAM for 220-sample preprocess; raise if it OOMs
#BSUB -W 2880                               # walltime: 48h in minutes (220 samples + CNV + report)
#BSUB -o /research/rgs01/home/clusterHome/qtran/Melanoma/meqtrack/melanoma_450k_%J/lsf.%J.out
#BSUB -e /research/rgs01/home/clusterHome/qtran/Melanoma/meqtrack/melanoma_450k_%J/lsf.%J.err

set -euo pipefail

# ---- EDIT THESE PATHS (cluster-side, absolute) ------------------------------
REPO_DIR="/research/groups/orrgrp/projects/MeQTrack_app"                       # <-- EDIT: MeQTrack clone on the cluster
SAMPLESHEET="/research/groups/orrgrp/projects/Melanoma/samplesheet_450k_hpc.csv" # <-- EDIT: cluster path of the 220-sample sheet
OUTPUT_DIR="/research/rgs01/home/clusterHome/qtran/Melanoma/meqtrack/melanoma_450k_${LSB_JOBID}"  # run dir keyed to the LSF job id (matches %J in -o/-e above, so logs land here)
ARRAY_TYPE="450k"
THREADS=2                                   # keep equal to "#BSUB -n" above

# Pipeline steps to run, in order. Each runs as its own invocation; the pipeline
# reloads prior results from disk between steps. Note: the "qc" step already runs
# probe filtering, so there's no separate "filtering" entry. reference_projection
# is OMITTED here (its reference beta-matrices are gitignored / large and not
# melanoma-specific). To include it, insert reference_projection before cnv. To
# run the whole pipeline in one process instead, set: STEPS=(all)
STEPS=(preprocess qc dim_reduction cnv visualization)
# -----------------------------------------------------------------------------

# Cluster modules: R (>= 4.4) for the pipeline, pandoc for the HTML report.
module load R/4.4.0
module load pandoc/3.8.3

echo "Host:        $(hostname)"
echo "R:           $(command -v Rscript) ($(Rscript -e 'cat(R.version.string)' 2>/dev/null))"
echo "Repo:        $REPO_DIR"
echo "Samplesheet: $SAMPLESHEET"
echo "Output:      $OUTPUT_DIR"
echo "Array type:  $ARRAY_TYPE   Threads: $THREADS"

mkdir -p "$OUTPUT_DIR"

# Mirror everything (stdout + stderr, incl. pipeline errors) into the run dir so
# the full log lives with the results. LSF's own -o/-e land here too (same %J
# dir), but this guarantees the log regardless of LSF's -o dir-creation behavior.
exec > >(tee "$OUTPUT_DIR/pipeline.${LSB_JOBID}.log") 2>&1

# Run from the repo root so the project .Rprofile activates renv — this is what
# provides missMethyl, sesame, minfi, conumee2, etc. The pipeline itself
# setwd()s into pipeline/, so --input/--output/--data_dir are absolute.
cd "$REPO_DIR"

# One-time provisioning of the renv library on the cluster. Safe to leave
# uncommented (installed packages are skipped); comment out once provisioned.
# Rscript setup.R

for step in "${STEPS[@]}"; do
  echo "===== $(date '+%F %T')  step: ${step} ====="
  Rscript pipeline/methylation_pipeline.R \
    --input     "$SAMPLESHEET" \
    --output    "$OUTPUT_DIR" \
    --data_dir  "$REPO_DIR/pipeline/data" \
    --array_type "$ARRAY_TYPE" \
    --threads   "$THREADS" \
    --step      "$step"
done

echo "===== $(date '+%F %T')  pipeline finished -> $OUTPUT_DIR ====="

echo "Pipeline finished. Results in: $OUTPUT_DIR"
