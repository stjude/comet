# `tumor_classification/`

**Task:** classify a sample into one of several tumor/disease types (a 16-way disease classification problem: 15 pediatric solid-tumor types plus a normal-tissue class — see `val_diseases` in `splitsD1_final_1norm_nomisc.p`), from three data modalities: whole-slide image (WSI), methylation, and a combined "integrative" model.

This folder contains only **wrapper scripts**. All of the actual model-fitting logic — PCA, LDA/QDA/Random Forest/Logistic Regression/MLP, hyperparameter grids — lives in [`../core_code/`](../core_code/README.md) and is not duplicated here. A wrapper's job is narrow and mechanical:

1. Load a train/test **splits pickle** (one entry per cross-validation split).
2. For each split, write the slide list + labels for that split's train/test sets to a temporary CSV.
3. Call the corresponding core classification function with those CSVs plus the relevant representation file(s) from [`../data/`](../data/README.md).
4. Collect the returned per-(classifier, hyperparameter) result rows across all splits, tag each with its split index, and summarize (mean/std of every metric) grouped by classifier/hyperparameter.
5. Save both the raw per-split results and the summary table to a single output pickle.

## Files

| File | Modality | Core function used |
|---|---|---|
| `run_wsi_tumor_classification.py` | WSI (SAMPLER representation) | `core_code.wsi_tumor_classification_core.classify_wsi_sampler` |
| `run_methylation_tumor_classification.py` | Methylation | `core_code.methylation_tumor_classification_core.classify_methylation` |
| `run_integrative_tumor_classification.py` | WSI + methylation | `core_code.integrative_tumor_classification_core.classify_integrative` |

All three share an identical structure (`run_all_splits`, `write_split_csvs`, `summarize_results`, a `main()` with an `argparse` CLI) — differing only in which core function they call and which representation file(s) they require as input. Every train/test split dictionary is expected to have the schema `{'slides_train': [...], 'slides_lab_train': [...], 'slides_test': [...], 'slides_lab_test': [...]}`.

### Splits file

`splitsD1_final_1norm_nomisc.p` is the disease-classification cross-validation splits file used for the results reported for this task: cell-line samples excluded, one representative slide per patient, stratified train/test splits (~80/20) repeated over many random seeds (`split_0`, `split_1`, ...). The pickle's `val_diseases` key gives the ordered list of the 16 disease-group class names that the integer labels (`0`-`15`) refer to.

### Usage

```bash
# WSI-only
python run_wsi_tumor_classification.py \
    --splits splitsD1_final_1norm_nomisc.p \
    --sampler-reps ../data/univ2_20x_sample_reps.p \
    --n-splits 5 --out wsi_tumor_class_results_5splits.p

# Methylation-only
python run_methylation_tumor_classification.py \
    --splits splitsD1_final_1norm_nomisc.p \
    --methylation-reps ../data/meth_top5000_sample_reps.p \
    --n-splits 5 --out methylation_tumor_class_results_5splits.p

# Integrative (WSI + methylation)
python run_integrative_tumor_classification.py \
    --splits splitsD1_final_1norm_nomisc.p \
    --sampler-reps ../data/univ2_20x_sample_reps.p \
    --methylation-reps ../data/meth_top5000_sample_reps.p \
    --n-splits 5 --out integrative_tumor_class_results_5splits.p
```

`--n-splits` limits how many of the splits in the pickle are actually run (useful for a fast smoke test); omitting it runs every split present in the file. Each script prints the summary table to stdout and writes a pickle containing both `all_results` (one row per split x classifier x hyperparameter combination) and `summary` (the same, aggregated to mean/std per classifier x hyperparameter, across splits).

## How this connects to the rest of the repo

- Depends on [`../core_code/`](../core_code/README.md) for all model-fitting logic.
- Depends on [`../data/`](../data/README.md) for the WSI SAMPLER representations and methylation representations.
- The same core classification functions used here are reused, unmodified, by [`../RMSsubtyping/`](../RMSsubtyping/README.md) (RMS subtype classification) and [`../NB_MYCN/`](../NB_MYCN/README.md) (MYCN amplification classification) — this folder is the simplest example of the wrapper pattern that those two folders also follow.
