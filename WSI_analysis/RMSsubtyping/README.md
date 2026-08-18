# `RMSsubtyping/`

**Task:** classify a rhabdomyosarcoma (RMS) sample into one of its three major histologic subtypes — **ERMS** (embryonal), **ARMS** (alveolar), or **Spindle** cell/sclerosing — from WSI, methylation, integrative, and tile-level representations.

This folder contains **only wrapper scripts**. It intentionally has no core modeling code of its own: both classification approaches below reuse core code that already exists for other tasks in this repository, pointed at RMS-subtype-specific splits and label counts.

## Files

| File | Approach | Core code reused |
|---|---|---|
| `run_rms_subtype_classification.py` | Slide/sample-level classification (WSI / methylation / integrative) | [`../core_code/wsi_tumor_classification_core.py`](../core_code/README.md), `methylation_tumor_classification_core.py`, `integrative_tumor_classification_core.py`, via [`../tumor_classification/`](../tumor_classification/README.md)'s wrapper functions |
| `run_rms_subtype_mil.py` | Tile-level, multiple-instance learning | [`../core_code/mil_classification_core.py`](../core_code/README.md) |
| `rms_subtype_cv_splits.p` | — | cross-validation split definitions (see below) |

### `run_rms_subtype_classification.py`

Rather than re-implementing "load splits -> write CSVs -> call core classifier -> summarize" a third time, this script directly **imports and calls the `run_all_splits`/`summarize_results` functions already defined in [`../tumor_classification/`](../tumor_classification/README.md)** (`run_wsi_tumor_classification.py`, `run_methylation_tumor_classification.py`, `run_integrative_tumor_classification.py`), simply pointing them at `rms_subtype_cv_splits.p` instead of the disease-classification splits file. This is the clearest illustration in the repository of the "wrapper around wrapper" pattern: `tumor_classification/`'s wrappers are themselves thin wrappers around `core_code/`, and this script is a thin wrapper around *those* wrappers.

The one behavioral difference from the disease-classification task: RMS subtyping restricts the classifier grid to `['lda', 'lr', 'mlp']` (excluding QDA and Random Forest), matching which model families are actually used for this particular subtyping question.

```bash
python run_rms_subtype_classification.py --modality wsi \
    --sampler-reps ../data/univ2_20x_sample_reps.p

python run_rms_subtype_classification.py --modality methylation \
    --methylation-reps ../data/meth_top5000_sample_reps.p

python run_rms_subtype_classification.py --modality integrative \
    --sampler-reps ../data/univ2_20x_sample_reps.p \
    --methylation-reps ../data/meth_top5000_sample_reps.p

# add --n-splits 3 for a quick test; default splits file is rms_subtype_cv_splits.p
```

### `run_rms_subtype_mil.py`

Same wrapper pattern as `../NB_MYCN/run_nb_mycn_mil.py`: reads `rms_subtype_cv_splits.p`, writes per-split train/test CSVs, and calls `core_code.mil_classification_core.classify_mil` against the per-tile representations in `../data/rms_tiles/` (tumor-tile-filtered by default, using `../data/RMS_tumor_indicator.p`). Runs all 3 MIL architectures x 3 hyperparameter variants (see [`../core_code/README.md`](../core_code/README.md)) per split.

```bash
python run_rms_subtype_mil.py --tiles-dir ../data/rms_tiles --out results.p
python run_rms_subtype_mil.py --tiles-dir ../data/rms_tiles --n-splits 2   # quick test
python run_rms_subtype_mil.py --tiles-dir ../data/rms_tiles --device cuda # force GPU
```

### `rms_subtype_cv_splits.p`

100 stratified train/test splits (80/20, `sklearn.model_selection.train_test_split` with `stratify` on subtype) over the RMS-subtype cohort (445 slides: ERMS/ARMS/Spindle). Each `split_i` entry has the schema `{'slides_train', 'slides_lab_train', 'slides_test', 'slides_lab_test', 'patients_train', 'patients_test'}` — one representative slide per patient, integer labels `0`=ERMS, `1`=ARMS, `2`=Spindle (see the file's own `val_subtypes` key for the label ordering).

## How this connects to the rest of the repo

- Reuses the slide/sample-level classification wrappers from [`../tumor_classification/`](../tumor_classification/README.md) (which in turn call [`../core_code/`](../core_code/README.md)) for `run_rms_subtype_classification.py`.
- Calls [`../core_code/mil_classification_core.py`](../core_code/README.md) directly for `run_rms_subtype_mil.py`.
- Consumes representation files from [`../data/`](../data/README.md) (`univ2_20x_sample_reps.p`, `meth_top5000_sample_reps.p`, `rms_tiles/`, `RMS_tumor_indicator.p`).
- No code here is subtype-specific beyond the splits file and the `['lda','lr','mlp']` classifier restriction — everything else is identical to the disease-level and NB-MYCN classification tasks.
