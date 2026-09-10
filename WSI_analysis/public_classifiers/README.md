# Public classifiers

This folder contains a subset of the final, fully-fitted COMET classifiers/survival models
(built 2026-09-02) provided so that others can apply them directly to their own whole-slide-image
(WSI) foundation-model features, without retraining. These are not all of the classifiers described
in the manuscript — just the subset we consider most broadly useful to external users.

Each `.pkl` file is a Python dict of `{model_key: {..fitted sklearn/lifelines objects.., metadata}}`,
loadable with `pickle.load()`. Every entry documents its own required preprocessing order (e.g.
feature selection → PCA → classifier) — see each file's embedded `'__README__'` string for full
detail (`pickle.load(open(f, 'rb'))['__README__']`), summarized below. All were refit on the full
available COMET cohort for that task (not a single CV fold); reported accuracy/C-index figures are
CV estimates from held-out folds using the same configuration, not training-set performance.

| File | Task | Classes / target | Models included |
|---|---|---|---|
| `tumor_classification_public_classifiers.pkl` | 16-class disease classification | 16 tumor/normal types | 4 WSI foundation models (MELD-CCDI, MELD-COMET, UniV2, VirchowV2), Sampler aggregation |
| `rms_subtyping_public_classifiers.pkl` | RMS histologic subtyping | ERMS / ARMS / Spindle | 4 WSI foundation models (Sampler) + 4 slide-level foundation models (CHIEF, GigaPath-slide, PRISM, TITAN) |
| `mycn_amplification_public_classifiers.pkl` | NB MYCN status | Not amplified / Amplified | Same 8 foundation models as above |
| `adrn_mes_public_classifiers.pkl` | NB cell state | ADRN / MES | Same 8 foundation models as above |
| `rms_survival_public_classifiers.pkl` | RMS overall-survival (Cox) | Predicted hazard | Up to 4 WSI + 4 slide-level foundation models (some slide-FM models omitted where the feature-selection gate failed) |

## Notes

- **`tumor_classification_public_classifiers.pkl`** is a reduced, no-PCA-refit subset (4 of the full
  33 fitted classifiers) chosen to keep the file small enough to host directly in this repo. It is
  restricted to H&E-only, Sampler (all-tile) aggregation, for the 4 WSI foundation models listed
  above. The full 33-classifier release and other intermediate variants (all H&E-only classifiers,
  PCA-retained 4-model subset) are much larger (up to ~4.7 GB) and are not hosted in this repository;
  contact the authors if you need the full release.
- **`adrn_mes_public_classifiers.pkl`** / **`mycn_amplification_public_classifiers.pkl`**: each of
  the 8 models uses whichever classifier (MLDA or elastic-net logistic regression) had the best CV
  balanced accuracy for that specific model; see each entry's `'classifier_type'`.
- **`rms_subtyping_public_classifiers.pkl`**: all 8 models use MLDA (`LinearDiscriminantAnalysis`);
  PCA is kept (not dropped) for this task since dropping it cost 3-7pp balanced accuracy in testing.
- **`rms_survival_public_classifiers.pkl`**: fit on a small cohort (96 patients, 23 deaths) with a
  fixed pipeline (Welch's t-test filter → PCA(50) → elastic-net Cox PH model via
  `lifelines.CoxPHFitter`). CV C-index is modest (~0.58-0.66) with high fold-to-fold variance —
  treat these models as considerably less stable than the classification releases.

## How to apply to your own data

```python
import pickle

with open("tumor_classification_public_classifiers.pkl", "rb") as f:
    release = pickle.load(f)

print(release["__README__"])       # full documentation for this file
model = release["UniV2_Sampler"]    # one fitted model entry
label_classes = release["label_classes"]

# Typical per-entry preprocessing order (see the file's __README__ / per-entry 'description'
# for the exact steps that apply to your chosen model key):
x = your_feature_vector  # same foundation-model representation used to build this release
if model.get("feature_indices") is not None:
    x = x[model["feature_indices"]]
if model.get("pca") is not None:
    x = model["pca"].transform([x])
if model.get("scaler") is not None:
    x = model["scaler"].transform(x)
pred = model["classifier"].predict(x)
predicted_label = label_classes[pred[0]]
```

For the survival models (`rms_survival_public_classifiers.pkl`), apply `feature_indices` then `pca`
as above, then call `model["cox_model"].predict_partial_hazard(df)` on a DataFrame with columns
`f0..f<n_pca-1>` (higher output = higher predicted hazard).
