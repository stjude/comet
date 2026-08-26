# `NB_ADRN_MES/`

**Task:** classify a neuroblastoma (NB) WSI as **adrenergic (ADRN)** or **mesenchymal (MES)** — a binary cell-state distinction relevant to neuroblastoma biology and treatment resistance — from the WSI SAMPLER representation alone.

This is the simplest task folder in the repository: a single wrapper script and a single splits file, reusing the same WSI classification core code used everywhere else.

## Files

| File | Purpose |
|---|---|
| `run_nb_adrn_mes_classification.py` | Wrapper: reads the splits below, calls `core_code.wsi_tumor_classification_core.classify_wsi_sampler`, aggregates metrics across splits. |
| `nb_adrn_mes_cv_splits.p` | Train/test cross-validation splits over the ADRN/MES cohort. |

Structurally identical to `../tumor_classification/run_wsi_tumor_classification.py` / `../NB_MYCN/run_nb_mycn_wsi_classification.py` — the only difference is the default splits file and the binary ADRN/MES label.

```bash
python run_nb_adrn_mes_classification.py --sampler-reps ../data/univ2_20x_sample_reps.p --out results.p
python run_nb_adrn_mes_classification.py --sampler-reps ../data/univ2_20x_sample_reps.p --n-splits 5   # quick test
```

### `nb_adrn_mes_cv_splits.p`

100 slide-level stratified train/test splits (80/20) over a 332-slide cohort (90 MES, 242 ADRN; label `1`=MES, `0`=ADRN). Each `split_i` has the schema `{'train_idx', 'test_idx', 'slides_train', 'slides_test', 'slides_lab_train', 'slides_lab_test'}`.

### Cohort construction

The ADRN/MES label itself is **not** produced by anything in this folder or reproduced from scratch in this repository. It comes from a semi-supervised procedure (documented, but not re-executable, in `../data/gen_nb_adrn_mes_labels.py`): a small set of slides with clear ADRN/MES morphology on H&E was manually identified, that seed set was expanded via a nearest-neighbor similarity search in feature space, and any remaining slide was labeled from an independent methylation-derived mesenchymal-probability score instead. The final result — just the slide list and its binary label, with none of the seed slides, expansion mechanics, or probe-level probabilities exposed — is cached at [`../data/NB_ADRN_MES_slides.p`](../data/README.md#nb_adrn_mes_slidesp).

## How this connects to the rest of the repo

- Calls [`../core_code/wsi_tumor_classification_core.py`](../core_code/README.md) — the exact same function used by [`../tumor_classification/`](../tumor_classification/README.md) and [`../NB_MYCN/`](../NB_MYCN/README.md); this folder adds no new modeling code.
- Consumes [`../data/univ2_20x_sample_reps.p`](../data/README.md) (WSI representations) and [`../data/NB_ADRN_MES_slides.p`](../data/README.md) (cohort/label source, folded into `nb_adrn_mes_cv_splits.p`).
- Unlike the other classification tasks in this repository, only the WSI modality is analyzed here — methylation, integrative, and tile-level variants of this task are not part of this repository.
