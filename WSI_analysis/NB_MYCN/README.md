# `NB_MYCN/`

**Task:** predict *MYCN* gene amplification status (binary: amplified vs. not amplified) in neuroblastoma (NB), a key molecular biomarker associated with high-risk disease, from five different modalities/representations of the same underlying WSI+methylation cohort.

Like [`../RMSsubtyping/`](../RMSsubtyping/README.md), this folder contains **only wrapper scripts**. Every classifier family used here already exists in [`../core_code/`](../core_code/README.md); this folder's job is to point that code at the NB-MYCN cohort's splits and representation files, one wrapper per modality.

## Files

| File | Modality | Core code used |
|---|---|---|
| `run_nb_mycn_wsi_classification.py` | WSI (SAMPLER representation) | `core_code.wsi_tumor_classification_core.classify_wsi_sampler` |
| `run_nb_mycn_methylation_classification.py` | Methylation | `core_code.methylation_tumor_classification_core.classify_methylation` |
| `run_nb_mycn_integrative_classification.py` | WSI + methylation | `core_code.integrative_tumor_classification_core.classify_integrative` |
| `run_nb_mycn_tile_classification.py` | Tile-level (mean-softmax aggregation) | `core_code.tile_level_classification_core.classify_tile_level` |
| `run_nb_mycn_mil.py` | Tile-level (multiple-instance learning) | `core_code.mil_classification_core.classify_mil` |
| `nb_mycn_cv_splits.p` | — | cross-validation split definitions (see below) |

All five wrappers follow the exact same pattern used throughout this repository: read a train/test splits pickle, write per-split CSVs (`slide,label`), call one core classification function, and aggregate the returned metrics across splits.

- The first three (WSI / methylation / integrative) are functionally identical in structure to their same-named counterparts in [`../tumor_classification/`](../tumor_classification/README.md) — the only differences are the default splits file (`nb_mycn_cv_splits.p`) and the binary (vs. 16-way or 3-way) label space.
- The last two (tile-level and MIL) are the first place in this repository these two core modules are used for MYCN specifically; both use `../data/NB_tiles/` (rather than `../data/rms_tiles/`) as the tile-feature directory and `../data/NB_tumor_indicator.p` (rather than the RMS tumor indicator) for tumor-tile filtering.

### `nb_mycn_cv_splits.p`

100 slide-level stratified train/test splits (80/20) over a 319-slide NB cohort (81 MYCN-amplified, 238 non-amplified). Each `split_i` has the schema `{'train_idx', 'test_idx', 'slides_train', 'slides_test', 'slides_lab_train', 'slides_lab_test'}`; label `1` = MYCN-amplified, `0` = not amplified.

### Usage

```bash
# WSI
python run_nb_mycn_wsi_classification.py --sampler-reps ../data/univ2_20x_sample_reps.p --out results.p

# Methylation
python run_nb_mycn_methylation_classification.py --methylation-reps ../data/meth_top5000_sample_reps.p --out results.p

# Integrative
python run_nb_mycn_integrative_classification.py \
    --sampler-reps ../data/univ2_20x_sample_reps.p \
    --methylation-reps ../data/meth_top5000_sample_reps.p --out results.p

# Tile-level (mean-softmax aggregation)
python run_nb_mycn_tile_classification.py --tiles-dir ../data/NB_tiles --out results.p

# Tile-level (multiple-instance learning)
python run_nb_mycn_mil.py --tiles-dir ../data/NB_tiles --out results.p

# add --n-splits 5 to any of the above for a quick test; add --device cuda to the two
# tile-level wrappers to force GPU
```

## How this connects to the rest of the repo

- All five wrappers call into [`../core_code/`](../core_code/README.md); none of the actual model-fitting logic lives in this folder.
- Consumes representation files from [`../data/`](../data/README.md): `univ2_20x_sample_reps.p`, `meth_top5000_sample_reps.p`, `NB_tiles/`, `NB_tumor_indicator.p`.
- The three sample-level wrappers (WSI/methylation/integrative) mirror [`../tumor_classification/`](../tumor_classification/README.md) exactly, adapted to a binary MYCN label; the two tile-level wrappers mirror [`../RMSsubtyping/run_rms_subtype_mil.py`](../RMSsubtyping/README.md) and are the first callers of `tile_level_classification_core.py` in the repo.
