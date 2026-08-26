# `core_code/`

This folder contains the **model-fitting logic** for every analysis in the repository. Every file here is a self-contained Python module: it has no notion of "which task" it is being used for, no hardcoded splits, no hardcoded output paths. It just takes:

1. a **train CSV** and a **test CSV** describing which samples belong to each set (and their labels),
2. a **path (or paths)** to a pre-built representation pickle (see [`../data/README.md`](../data/README.md)), and
3. an optional **`params` dictionary** to override any default hyperparameter,

and returns a list of result dictionaries — one row per (classifier, hyperparameter setting) combination it tried.

Nothing in this folder reads a task-specific splits file, does cross-validation bookkeeping, or writes output to disk — that is the job of the **wrapper** scripts that live in each task folder (`tumor_classification/`, `RMSsubtyping/`, `NB_MYCN/`, `RMS_survival/`, `NB_ADRN_MES/`). See the root [`README.md`](../README.md) for how everything fits together. Each wrapper imports exactly one function from one file here, calls it once per cross-validation split, and aggregates the returned rows.

This separation is intentional: it means the exact same core code can be (and is) reused across multiple biological questions — e.g. `wsi_tumor_classification_core.py` powers the disease-level WSI classifier, the RMS-subtype WSI classifier, and the NB-MYCN WSI classifier, unmodified.

## Design conventions shared by every file in this folder

- **Module-level defaults, dictionary override.** Every hyperparameter (PCA dimensions, regularization grids, learning rates, model sizes, etc.) is declared as a plain module-level variable near the top of the file, with a comment. A `resolve_params(params)` function copies these into a dict and lets a caller override any subset of them by passing `params={...}`. Nothing is hidden inside a function body — you can see every knob just by reading the top of the file.
- **Train-fold-only preprocessing.** Any fitted preprocessing step (PCA, standardization, feature selection) is always fit on the **training split only** and then applied to both train and test. This avoids test-set leakage into the preprocessing.
- **Exhaustive small grid search, not internal parallelization.** Each core function loops over its full hyperparameter grid with plain Python `for` loops and returns one result row per combination. There is no multiprocessing/threading inside these files — the goal is to make the exact computation being performed obvious to a reader, not to be maximally fast. (Thread counts for the underlying numeric libraries are still capped via `torch.set_num_threads` / the `OMP_NUM_THREADS` environment variable, so that these scripts behave politely on a shared multi-core machine.)
- **Failures are recorded, not fatal.** Every classifier fit is wrapped in a `try/except`. If a particular (classifier, hyperparameter) combination fails for a given split (e.g. a covariance matrix is singular, or a feature-selection step finds too few significant features), the function records `ok=False` and a `reason` string instead of crashing — so a rare degenerate split doesn't invalidate the entire multi-split run. Downstream summaries simply average over the results that succeeded.
- **Same metrics vocabulary.** Every classification core reports: ROC-AUC (`ovr_auc`; one-vs-rest macro AUC for >2 classes, standard AUC of the positive class for binary problems), `accuracy`, `balanced_accuracy`, `f1_macro`, and a `confusion_matrix`. Every survival core reports Harrell's concordance index (`c_index`) and a companion `auc`. This lets every wrapper's summarization code be nearly identical.

---

## Tumor / disease classification family

These three files answer the same kind of question — "classify a slide/sample into one of several classes" — from three different data modalities. They are used directly for whole-cohort disease classification and are reused unmodified for RMS subtyping and for NB MYCN-amplification prediction (see those tasks' own READMEs).

### `wsi_tumor_classification_core.py` — H&E / SAMPLER-representation classifier

**Input representation:** a *SAMPLER representation* — for each slide, a small `(10, feature_dim)` array: 10 percentiles (5th, 15th, ..., 95th) computed **over all tiles of the slide**, for a single tile-level foundation-model embedding (see [`../data/README.md`](../data/README.md) for exactly how this is built). This compresses an entire whole-slide image, which may contain tens of thousands of tiles, into a single fixed-size vector per slide (the 10 percentile rows are flattened together), which is small enough to fit many slides in memory and to feed directly into classical ML models.

**Pipeline (`classify_wsi_sampler`):**
1. Load the SAMPLER pickle and the train/test CSVs (columns: `slide,label`).
2. Build `(X, Y)` matrices by looking up each slide's flattened SAMPLER vector; rows with a missing slide or a negative label are dropped.
3. For each candidate PCA dimensionality in `PC_dims`:
   - Fit a `PCA` on the training features only (dimensionality is capped to `min(target, n_train_samples, n_features)`), and transform both train and test.
   - Fit every requested classifier on the PCA-reduced training features and score it on the PCA-reduced test features:
     - **LDA** (`sklearn.discriminant_analysis.LinearDiscriminantAnalysis`, SVD solver)
     - **QDA** (`QuadraticDiscriminantAnalysis`), swept over a small grid of `reg_param` shrinkage values (helps when a class's covariance estimate is close to singular)
     - **Random Forest**, swept over `num_trees` x `rf_max_depths`
     - **Logistic Regression** (`sklearn.linear_model.LogisticRegression`, `class_weight='balanced'`), swept over a regularization grid `Cvec`
     - **2-layer MLP** (PyTorch): `Linear(in, hidden) -> ReLU -> Dropout -> Linear(hidden, n_classes)`, trained with Adam + cross-entropy + L1/L2 weight regularization, CPU-only, using the same architecture family as the WSI classifiers used elsewhere in the project.
4. Every (PCA dim, classifier, hyperparameter) combination is scored independently and returned as its own result row.

### `methylation_tumor_classification_core.py` — methylation classifier

**Input representation:** a methylation *beta-value* vector per sample (restricted to the most-variable CpG probes; see the data README), joined to a slide by that slide's Sentrix ID.

**Pipeline (`classify_methylation`):** structurally the same as the WSI core (PCA on train-only, then classify), but with a **narrower model set** — `LDA`, a single-layer torch Logistic Regression, and the 2-layer torch MLP — deliberately excluding QDA and Random Forest, mirroring which model families are actually used for methylation-based classification in this project (methylation feature vectors are very high-dimensional and low-sample-count, where QDA/RF are not well suited). The torch LR and torch MLP share one `build_torch_net`/`torch_predict_proba` implementation, differing only in whether a hidden layer is present.

### `integrative_tumor_classification_core.py` — combined WSI + methylation classifier

**Input representation:** the concatenation of a slide's WSI SAMPLER vector and its matched sample's methylation beta vector.

**Pipeline (`classify_integrative`):** same classifier family as the WSI core (LDA/QDA/RF/LR/MLP), but with **two alternative ways of reducing the concatenated feature vector**, controlled by the `pca_modes` hyperparameter (both are run, by default):
- `'concat_pca'` — fit a single PCA directly on the raw `[WSI-features || methylation-features]` concatenation. To keep the resulting feature count comparable to the alternative below, this mode uses **twice** the configured `PC_dims` value.
- `'double_pca'` — fit two *separate* PCAs, one on the WSI block and one on the methylation block (each reduced to `PC_dims` components), then concatenate the two reduced blocks. This lets each modality contribute a fixed, modality-balanced number of components regardless of its raw dimensionality, rather than letting whichever modality has larger raw variance dominate a single joint PCA.

Both PCA modes are always fit on the training fold only.

---

## RMS survival family

These three files fit the same **penalized Cox proportional-hazards** survival model on three different feature sets. All three share nearly identical logic (feature screening -> PCA -> Cox -> risk-score metrics); they are kept as separate files (rather than one file with a `modality` argument) to match this repo's convention of one self-contained, independently readable file per data modality.

**Common pipeline (`classify_rms_survival_sampler` / `_methylation` / `_integrative`):**
1. Build the feature matrix `X` for train/test (WSI SAMPLER vector, methylation beta vector, or their raw concatenation for the integrative version — same feature-construction logic as the corresponding tumor-classification core above).
2. Alongside `X`, three per-sample fields are carried through: a continuous **survival time**, a binary **event indicator** (died vs. censored/alive), and a binary **3-year mortality label** (died within 3 years vs. alive/censored with ≥3 years of follow-up). The 3-year label is *not* the survival outcome being modeled directly — it is used only as a fast, well-balanced group label for feature screening and for a companion AUC metric (see below).
3. **Feature screening:** an independent-samples t-test (Welch's, unequal variance) compares each feature's distribution between the two 3-year-label groups, computed on the training fold only. Features with `p < 0.05` are kept; if fewer than the target PCA dimensionality pass this filter, that configuration is skipped for that split (recorded as a failure with an explanatory reason, not silently dropped).
4. **PCA** is fit on the training fold's screened features (50 components by default) and applied to both train and test.
5. **Cox model:** a penalized Cox proportional-hazards model (`lifelines.CoxPHFitter`, elastic-net penalty — `penalizer=0.1`, `l1_ratio=0.5`) is fit on the training fold's PCA features plus (time, event). The fitted coefficients define a linear risk score, which is applied to the test fold's PCA features and min-max normalized to `[0, 1]`.
6. **Metrics:** Harrell's concordance index (`sksurv.metrics.concordance_index_censored`) between the test risk score and the true (time, event) pairs, plus an AUC of the same risk score against the 3-year mortality label (a fixed-horizon proxy that is easier to sanity-check than the concordance index alone).

### `rms_survival_sampler_core.py`
WSI SAMPLER representation only.

### `rms_survival_methylation_core.py`
Methylation beta-value representation only (joined by Sentrix ID, resolved once upstream — see [`../RMS_survival/README.md`](../RMS_survival/README.md)).

### `rms_survival_integrative_core.py`
The raw concatenation `[methylation || WSI SAMPLER]`, with a single t-test screen and a single PCA fit over the combined vector (no separate-PCA-then-concatenate option here — this matches how the survival analysis is actually done, unlike the tumor-classification integrative core above, which offers both).

---

## Tile-level family

Unlike the SAMPLER-based cores above (which compress an entire slide into one percentile vector before any modeling happens), the two files below work directly on **individual tile-level feature vectors**, one row per tile, and only combine tiles into a single slide-level prediction as the very last step.

### `mil_classification_core.py` — Multiple-Instance Learning (MIL)

**Input representation:** a per-slide "bag" of tile-level foundation-model feature vectors — a `(n_tiles, feature_dim)` array (see the data README for how these tile pickles are built). By default, bags are restricted to tiles flagged as tumor by a precomputed tumor/non-tumor indicator (see `tile_tum_detector.py` below and the data README); slides with too few tumor tiles remaining are dropped.

**Models (`ABMILClassifier`, `CLAMSB`, `DSMIL`, all PyTorch):**
- **ABMIL** — gated-attention multiple-instance learning (Ilse et al.): each tile is projected to a hidden representation, an attention network produces one scalar weight per tile, and the slide representation is the attention-weighted sum of tile representations, followed by a linear classification head.
- **CLAM-SB** — a single-branch variant of the CLAM architecture (Lu et al.): structurally similar gated attention pooling, trained with a plain bag-level cross-entropy loss only (no auxiliary instance-level clustering loss).
- **DSMIL** — dual-stream MIL (Li et al.): combines an attention-pooled bag prediction with a "critical instance" prediction (the single tile with the strongest per-tile logit), averaging the two branches' logits.

Each architecture is trained with three preset hyperparameter variants (`v1`/`v2`/`v3` — different hidden/attention widths, dropout, learning rate, and epoch budgets), for **3 models x 3 variants = 9 combinations** evaluated per split. Training runs bag-by-bag (one slide per optimizer step, no batching) with early stopping on training loss — deliberately simple so the same code runs unmodified on CPU or GPU without a custom batched-bag data loader.

### `tile_level_classification_core.py` — per-tile classifier with slide-level aggregation

**Input representation:** the same per-slide tile-feature pickles as the MIL core, tumor-tile-filtered the same way.

**Pipeline (`classify_tile_level`):** rather than treating a slide as a single "bag" processed by an attention network (as MIL does), this core:
1. Flattens **every kept tile from every training slide** into one large tile-level training set, where each tile inherits its parent slide's label.
2. Trains a simple **per-tile** classifier — logistic regression or a 2-layer MLP (`Linear -> ReLU -> Dropout -> Linear`), both implemented in PyTorch with an internal standardization layer — with Adam, cross-entropy, and optional class-balanced loss weighting.
3. At evaluation time, every tile of a test slide is scored independently, the resulting per-tile softmax probabilities are **averaged across all of that slide's tiles**, and the slide-level prediction is the arg-max of that averaged probability vector.

This "mean-softmax aggregation" is a much simpler alternative to attention-based MIL pooling, useful as a baseline for how much the learned attention mechanism in `mil_classification_core.py` is actually contributing.

### `tile_tum_detector.py` — pretrained tumor/non-tumor tile classifier

A small, self-contained wrapper around a **pre-fit** `sklearn.discriminant_analysis.LinearDiscriminantAnalysis` classifier (its weights are stored in `../data/`; see the data README) that scores a single tile as tumor vs. non-tumor from a tile-level foundation-model feature vector. This is not trained here — it is a fixed, previously-trained artifact — and is used only as a **preprocessing filter**: `../data/gen_*_tumor_indicator.py` scripts run it once over every tile in a cohort and cache the resulting boolean tumor mask to disk (see the data README), so that `mil_classification_core.py` and `tile_level_classification_core.py` can simply look up "is this tile tumor?" from a dictionary instead of re-running the detector at training/evaluation time.

---

## Quick reference: which core is used by which task

| Core file | Used directly by | Also reused (unmodified) by |
|---|---|---|
| `wsi_tumor_classification_core.py` | `tumor_classification/` | `RMSsubtyping/`, `NB_MYCN/`, `NB_ADRN_MES/` |
| `methylation_tumor_classification_core.py` | `tumor_classification/` | `RMSsubtyping/`, `NB_MYCN/` |
| `integrative_tumor_classification_core.py` | `tumor_classification/` | `RMSsubtyping/`, `NB_MYCN/` |
| `mil_classification_core.py` | `RMSsubtyping/` | `NB_MYCN/` |
| `tile_level_classification_core.py` | `NB_MYCN/` | — |
| `tile_tum_detector.py` | `../data/gen_*_tumor_indicator.py` (offline, once) | — |
| `rms_survival_sampler_core.py` / `_methylation_core.py` / `_integrative_core.py` | `RMS_survival/` | — |
