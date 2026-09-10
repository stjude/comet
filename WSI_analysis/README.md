# COMET digital pathology / multi-omics analysis code

This repository contains the analysis code for a set of computational pathology and multi-omics classification/survival studies on pediatric solid tumors, built on top of whole-slide image (WSI) foundation-model features and DNA methylation array data. It accompanies the associated manuscript and is intended to let a reviewer or reader understand, reproduce the logic of, and run the core statistical/machine-learning analyses on a representative subset of the underlying data.

## What's in this repository

Five analysis tasks are implemented, each built from the same small set of reusable modeling building blocks:

| Task | Folder | Question | Modalities |
|---|---|---|---|
| Disease-level tumor classification | [`tumor_classification/`](tumor_classification/README.md) | Which of 16 tumor/normal classes does this sample belong to? | WSI, methylation, integrative |
| RMS histologic subtyping | [`RMSsubtyping/`](RMSsubtyping/README.md) | ERMS vs. ARMS vs. Spindle-cell rhabdomyosarcoma? | WSI, methylation, integrative, tile-level (MIL) |
| NB MYCN amplification | [`NB_MYCN/`](NB_MYCN/README.md) | Is this neuroblastoma *MYCN*-amplified? | WSI, methylation, integrative, tile-level, MIL |
| NB adrenergic vs. mesenchymal state | [`NB_ADRN_MES/`](NB_ADRN_MES/README.md) | Adrenergic (ADRN) or mesenchymal (MES) neuroblastoma cell state? | WSI |
| RMS survival prediction | [`RMS_survival/`](RMS_survival/README.md) | Time-to-event survival risk (Cox model) | WSI, methylation, integrative |

Every task is a thin **wrapper** over shared **core** modeling code — see [`core_code/README.md`](core_code/README.md) for the modeling logic itself (PCA/feature-selection strategies, every classifier and survival model, every hyperparameter grid) and each task folder's own README for how that task's wrapper scripts drive the core code end-to-end.

For a subset of these tasks, we also provide the final, fully-fitted classifiers/survival models themselves (not just the code to train them) so they can be applied directly to new data — see [`public_classifiers/`](public_classifiers/README.md).

## Repository layout

```
github_repo/
├── core_code/            # Shared, task-agnostic modeling code (classifiers + survival models)
├── data/                 # Minimal-example input data + the scripts that derived it from
│                         #   the internal cohort (see data/README.md for full disclosure)
├── tumor_classification/ # Task: disease-level classification (16 classes)
├── RMSsubtyping/         # Task: RMS histologic subtype classification
├── NB_MYCN/              # Task: NB MYCN amplification classification
├── NB_ADRN_MES/          # Task: NB adrenergic vs. mesenchymal classification
└── RMS_survival/           # Task: RMS survival prediction (Cox model)
```

Each folder above has its own `README.md` with full detail; start with [`core_code/README.md`](core_code/README.md) to understand the modeling approaches, then [`data/README.md`](data/README.md) to understand what data backs them, then whichever task folder is of interest.

## How a task is structured (the pattern used throughout this repo)

Every task folder follows the same three-layer design:

1. **Data** ([`data/`](data/README.md)) — pre-built, reviewer-shareable feature representations (WSI SAMPLER [1] vectors, per-tile UNIv2 [2] foundation-model features, methylation beta values) plus, for tasks with non-trivial cohort construction, a cached "cleaned cohort" pickle so that inclusion/exclusion logic is implemented exactly once and is fully auditable.
2. **Core code** ([`core_code/`](core_code/README.md)) — pure modeling functions. Each one takes a train CSV, a test CSV, and a path to a representation file, and returns a list of (classifier or model configuration, metric) result rows. No file I/O beyond reading its inputs, no task-specific assumptions, no knowledge of cross-validation.
3. **Task wrappers** (`tumor_classification/`, `RMSsubtyping/`, `NB_MYCN/`, `NB_ADRN_MES/`, `RMS_survival/`) — thin scripts that load a task-specific train/test splits file, loop over its splits, write the per-split CSVs the core code expects, call the appropriate core function once per split, and aggregate the results (mean/std of every metric) across splits.

Because the core code has no task-specific logic, most of it is **reused across multiple tasks unmodified** — e.g. the same WSI classification core powers disease classification, RMS subtyping, NB MYCN status, and NB ADRN/MES state. See the "quick reference" table at the bottom of [`core_code/README.md`](core_code/README.md) for exactly which core file backs which task.

## Modeling approaches used

- **Classification** (disease type, subtype, mutation status, cell state): PCA (or per-modality double-PCA for integrative data) on a training-fold-only fit, followed by one or more of: Linear/Quadratic Discriminant Analysis, Random Forest, Logistic Regression, and a small feed-forward neural network — swept over a fixed, fully-documented hyperparameter grid per modality (see [`core_code/README.md`](core_code/README.md) for exactly which grid applies to which modality).
- **Tile-level classification**: two complementary approaches operating directly on individual tile-level foundation-model features rather than a whole-slide summary — attention-based multiple-instance learning (ABMIL / CLAM-SB / DSMIL) and a simpler per-tile classifier with mean-softmax aggregation to the slide level.
- **Survival analysis**: a penalized (elastic-net) Cox proportional-hazards model, fit on PCA-reduced, t-test-screened features, evaluated by Harrell's concordance index and a 3-year-mortality AUC.

## Data and privacy

This repository ships **minimal-example, de-identified data** — representative subsets of the full study cohorts, sufficient to run every analysis end-to-end and inspect realistic results, but smaller than and not identical to the cohorts behind the associated publication's headline numbers. See [`data/README.md`](data/README.md) for a full, per-file account of:

- exactly what each shared data file contains and how it was derived,
- which upstream/raw sources (raw whole-slide-image caches, full methylation matrices, detailed clinical spreadsheets) were **not** included and why, and
- what cohort-cleaning/label-derivation logic was applied, even in cases where the underlying raw inputs to that logic are not shared.

Note: representations need to be downloaded from zenodo (https://zenodo.org/records/21997824) and copied to the data folder.

## Getting started

```bash
pip install -r requirements.txt
```

Each task folder's README gives copy-pasteable example commands. As a first smoke test, from `tumor_classification/`:

```bash
python run_wsi_tumor_classification.py \
    --splits splitsD1_final_1norm_nomisc.p \
    --sampler-reps ../data/univ2_20x_sample_reps.p \
    --n-splits 3 --out /tmp/results.p
```

This trains and evaluates every classifier/hyperparameter combination in [`core_code/wsi_tumor_classification_core.py`](core_code/README.md) on 3 cross-validation splits of the included disease-classification cohort and prints a summary table.

## Requirements

Python 3.9+ with `numpy`, `pandas`, `polars`, `scipy`, `scikit-learn`, `torch`, `lifelines`, `scikit-survival`, and `openpyxl` (see `requirements.txt`). All analyses run on CPU; a GPU is optional and only used by the tile-level/MIL neural-network models if `--device cuda` is requested.

## Contributing

Bug reports and small fixes are welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to get started, and please follow our [Code of Conduct](CODE_OF_CONDUCT.md). To report a security issue, see [`SECURITY.md`](SECURITY.md).

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE) (required because
it depends on and imports `scikit-survival`, which is GPL-3.0-licensed).

## References

[1] Mukashyaka P, Sheridan TB, Foroughi pour A, Chuang JH. SAMPLER: unsupervised representations for rapid analysis of whole slide tissue images. EBioMedicine. 2024 Jan 1;99.

[2] Chen RJ, Ding T, Lu MY, Williamson DF, Jaume G, Song AH, Chen B, Zhang A, Shao D, Shaban M, Williams M. Towards a general-purpose foundation model for computational pathology. Nature medicine. 2024 Mar;30(3):850-62.
