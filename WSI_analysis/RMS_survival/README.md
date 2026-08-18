# `RMS_survival/`

**Task:** predict rhabdomyosarcoma (RMS) patient survival outcome — both a continuous time-to-event risk score and a binary 3-year-mortality proxy — from WSI, methylation, and integrative (combined) representations, using a penalized Cox proportional-hazards model.

This folder contains the cohort-splitting script and the per-modality wrapper scripts for this task. All of the actual survival-modeling logic (feature screening, PCA, Cox model fitting, risk-score evaluation) lives in [`../core_code/`](../core_code/README.md); this folder only builds the cross-validation splits and feeds them, one split at a time, to that core code.

## Why this folder exists separately from cohort cleaning

Building the analysis-ready RMS survival cohort involves several non-trivial cleaning steps (excluding a histologic subtype out of scope for this analysis, requiring sufficient follow-up time, collapsing risk categories, restricting to one slide per patient, and filtering out slides with too few usable tissue tiles to be reliable — see [`../data/README.md`](../data/README.md)). Rather than repeating that cleaning logic inside every survival script, it is done **once** and cached to [`../data/RMS_survival_metadata.p`](../data/README.md#rms_survival_metadatap). Everything in this folder starts from that already-cleaned cohort, which keeps the survival-modeling code itself focused purely on the modeling question.

## Files

### `build_rms_survival_splits.py`

Builds `rms_survival_cv_splits.p`: 200 stratified train/test splits (80/20, stratified on the cleaned cohort's binary 3-year-mortality label, `sklearn.model_selection.train_test_split` with a fixed base random seed) over the cohort defined by `../data/RMS_survival_metadata.p`. Each `split_i` entry is simply `{'train_idx', 'test_idx'}` — integer indices into that metadata pickle's parallel `Y`/`slides`/`sentrix`/`patient_ids` arrays.

This is a from-scratch re-implementation of the same splitting recipe used for the originally published RMS survival analysis (same test size, same stratification target, same fixed-seed approach), rather than a direct copy of the original splits file, because the WSI representation shared in this repository ([`../data/univ2_20x_sample_reps.p`](../data/README.md)) covers a very slightly different slide set than the one behind the original production run — so the two cohorts are not guaranteed to be index-for-index identical, and the splits have to be regenerated against the exact cohort that ships in this repository.

```bash
python build_rms_survival_splits.py
```

### `run_rms_survival_sampler.py`, `run_rms_survival_methylation.py`, `run_rms_survival_integrative.py`

Three wrappers, one per modality, all sharing the same structure:

1. Load `../data/RMS_survival_metadata.p` (the cleaned cohort) and `rms_survival_cv_splits.p` (the split indices).
2. For each split, use `train_idx`/`test_idx` to slice out that split's rows and write them to a CSV with columns `slide, sentrix, time, event, binary_3yr` (`time`/`event` are the true survival time and event indicator used to fit and evaluate the Cox model; `binary_3yr` is only used internally, for feature screening and as an auxiliary AUC metric — see [`../core_code/README.md`](../core_code/README.md#rms-survival-family)).
3. Hand the train/test CSVs to the corresponding core function:

| Wrapper | Core function |
|---|---|
| `run_rms_survival_sampler.py` | `core_code.rms_survival_sampler_core.classify_rms_survival_sampler` |
| `run_rms_survival_methylation.py` | `core_code.rms_survival_methylation_core.classify_rms_survival_methylation` |
| `run_rms_survival_integrative.py` | `core_code.rms_survival_integrative_core.classify_rms_survival_integrative` |

4. Aggregate `c_index` (Harrell's concordance index) and `auc` (3-year-mortality AUC of the fitted risk score) across splits.

```bash
python run_rms_survival_sampler.py \
    --metadata ../data/RMS_survival_metadata.p \
    --splits rms_survival_cv_splits.p \
    --sampler-reps ../data/univ2_20x_sample_reps.p \
    --out rms_survival_sampler_results.p

python run_rms_survival_methylation.py \
    --metadata ../data/RMS_survival_metadata.p \
    --splits rms_survival_cv_splits.p \
    --methylation-reps ../data/meth_top5000_sample_reps.p \
    --out rms_survival_methylation_results.p

python run_rms_survival_integrative.py \
    --metadata ../data/RMS_survival_metadata.p \
    --splits rms_survival_cv_splits.p \
    --sampler-reps ../data/univ2_20x_sample_reps.p \
    --methylation-reps ../data/meth_top5000_sample_reps.p \
    --out rms_survival_integrative_results.p

# add --n-splits 5 to any of the above for a quick test (200 splits is the default and can be slow)
```

## How this connects to the rest of the repo

- Depends on [`../data/RMS_survival_metadata.p`](../data/README.md) for the cleaned cohort (built by `../data/gen_rms_survival_metadata.py`).
- Depends on [`../data/univ2_20x_sample_reps.p`](../data/README.md) and/or [`../data/meth_top5000_sample_reps.p`](../data/README.md) for the actual feature representations, depending on modality.
- Depends on [`../core_code/rms_survival_sampler_core.py`, `rms_survival_methylation_core.py`, `rms_survival_integrative_core.py`](../core_code/README.md) for the Cox modeling pipeline.
- Unlike the classification tasks elsewhere in this repository (which reuse a shared classifier core across multiple diseases), the survival core code here is specific to this task — no other folder in this repository performs a time-to-event analysis.
