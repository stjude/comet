# MELD Documentation

This document is the full reference for the `meld` package: what MELD is, how it is
implemented, how to install and use it, and how to troubleshoot common issues. For a short
overview, see the top-level [README](../README.md); for a runnable walkthrough, see
[`notebooks/quickstart.ipynb`](../notebooks/quickstart.ipynb).

## Table of contents

- [1. Overview](#1-overview)
- [2. Installation](#2-installation)
- [3. Repository layout](#3-repository-layout)
- [4. Model card](#4-model-card)
- [5. Data format](#5-data-format)
- [6. API reference](#6-api-reference)
- [7. Weights: format, size, and Hugging Face](#7-weights-format-size-and-hugging-face)
- [8. Troubleshooting](#8-troubleshooting)
- [9. Status and scope](#9-status-and-scope)

## 1. Overview

MELD is a small transformer ("ViT-style") fusion network that combines tile-level embeddings
from two histopathology foundation models — **UNI-v2** (1536-d) and **Virchow2** (1280-d) —
into a single, compact 3000-d representation per tile. It is trained in a self-supervised
manner (as a fusion autoencoder) directly on cached tile embeddings, i.e. it does not consume
raw whole-slide images itself.

Two released checkpoints, both trained for 20 epochs:

| Checkpoint key | Name | Trained/distilled on |
|---|---|---|
| `"comet"` | MELD-COMET | COMET dataset |
| `"ccdi"`  | MELD-CCDI  | CCDI dataset |

Both checkpoints share the exact same architecture (`meld.create_architecture`) and only
differ in their weights.

## 2. Installation

```bash
git clone <this-repo-url>
cd MELD_github
pip install -e .
# or, without installing as a package:
pip install -r requirements.txt
```

### Requirements

- Python >= 3.9
- `tensorflow >= 2.16` (which bundles **Keras 3** by default)
- `numpy`

**Keras version matters.** The released weight files were produced with Keras 3 and use the
modern `.weights.h5` format. Loading them under the legacy `tf.keras` bundled with
TensorFlow <= 2.15 (Keras 2.x) will fail — see [Troubleshooting](#8-troubleshooting). This
package was developed/tested against `tensorflow==2.20.0` / `keras==3.13.0`.

## 3. Repository layout

```
MELD_github/
├── README.md                    # quickstart overview
├── docs/
│   └── DOCUMENTATION.md         # this file
├── pyproject.toml               # pip-installable package metadata
├── requirements.txt
├── meld/
│   ├── __init__.py              # public API re-exports
│   ├── architecture.py          # create_architecture(), cosine_similarity_loss(), constants
│   ├── model.py                 # load_model(), set_inference_mode(), extract_features()
│   └── train.py                 # train(), finetune()
├── weights/
│   ├── meld_comet.weights.h5    # MELD-COMET, 20 epochs, weights-only
│   └── meld_ccdi.weights.h5     # MELD-CCDI, 20 epochs, weights-only
└── notebooks/
    └── quickstart.ipynb         # end-to-end usage walkthrough
```

## 4. Model card

### Architecture summary

MELD is a functional Keras model with four self-attention ("ViT") blocks operating on a
`(100, 60)`-shaped reshaping of the input, plus a bottleneck `feats` layer and two auxiliary
projection heads used only during training:

1. **Input projection**: `Dense(6000)` → `Reshape(100, 60)`.
2. **ViT blocks 1–2**: `MultiHeadAttention(8 heads, key_dim=60)` → `LayerNormalization` →
   MLP (`Dense(360, gelu)` → `Dropout` → `Dense(60, gelu)` → `Dropout`) → residual `Add`.
3. **Bottleneck**: `Flatten` → `Dense(3000, name="feats")` → per-sample normalization
   (`anorm`, mean/std over the 3000 features) → `Reshape(100, 30)`. **This normalized
   3000-d `feats` output is the actual MELD embedding** returned by
   `meld.extract_features`.
4. **ViT blocks 3–4**: same pattern as blocks 1–2, but operating on the 30-d bottleneck
   representation (`key_dim=30` then `key_dim=60`).
5. **Reconstruction head**: `Flatten` → `Dense(INPUT_DIM, linear)` — reconstructs the raw
   8448-d input.
6. **Auxiliary heads** (training-only): two small MLPs project the bottleneck back to
   Virchow2-space (1280-d) and UNI-v2-space (1536-d) respectively, used by the training loss
   to encourage the bottleneck to preserve per-backbone structure.
7. **Output**: `Concatenate([reconstruction (8448), virchow_proj (1280), uni_proj (1536)])`
   → total output dimension **11264**. Only the `feats` layer matters at inference time; the
   rest of the graph exists to make the training loss well-defined.

Total parameters: **157,285,550** (identical for both checkpoints, since they share the same
architecture).

### Verification

The `create_architecture` function in this repo reconstructs the training-time graph layer
by layer, in the exact same creation order as the original training code, so that
Keras's auto-generated layer names (`dense`, `dense_1`, `multi_head_attention`, ...) line up
with the names baked into the released `.weights.h5` files. This was verified by loading the
original full training checkpoints (`.keras` files, with optimizer state) and confirming
that predictions from `create_architecture() + load_weights(...)` are **bit-for-bit
identical** (max abs difference = 0.0) to the original checkpoints, on both the full model
output and the `feats` layer output, for both MELD-COMET and MELD-CCDI.

### Training objective (informational)

MELD is trained with `meld.architecture.cosine_similarity_loss`, which combines:
- Cosine similarity between the reconstructed and true Virchow2/UNI-v2 blocks.
- L1 and max-abs reconstruction error.
- Pairwise-distance preservation losses (L1 and cosine) between samples in a batch, for both
  the reconstruction and the auxiliary projection heads.

This loss is what `meld.train` / `meld.finetune` compile the model with by default; you do
not need to understand its internals to use the package, but it is documented in the
[API reference](#6-api-reference) for completeness/reproducibility.

## 5. Data format

MELD consumes **precomputed tile-level embeddings**, not raw images. Each row corresponds to
one tile and has shape `(8448,)`:

```
[ virchow2_a (1280), virchow2_b (1280), virchow2_c (1280),
  uni_v2_a   (1536), uni_v2_b   (1536), uni_v2_c   (1536) ]
```

- `virchow2_*`: three 1280-d Virchow2 embeddings for the tile (e.g. from different
  crops/scales, matching how tiles were cached for COMET/CCDI/COG during training).
- `uni_v2_*`: three 1536-d UNI-v2 embeddings for the same tile, same convention.

A batch of tiles is therefore shape `(N_TILES, 8448)`.

If your own feature-extraction pipeline produces a different number of crops/scales, or only
a single embedding per backbone:
- Repeat/tile the single embedding 3x to match the expected layout as a quick approximation, or
- Fine-tune the released checkpoint (`meld.finetune`) on data laid out your way — the input
  dimension only has to stay at 8448 with the same [Virchow2 block, UNI-v2 block] ordering.

Tiling and feature-caching pipelines (i.e. going from whole-slide images to these cached
embeddings) are **not included** in this minimal release; see [Status and scope](#9-status-and-scope).

## 6. API reference

Everything below is available directly from the top-level `meld` package (e.g. `meld.load_model`).

### Constants (`meld.architecture`)

| Name | Value | Description |
|---|---|---|
| `VIRCHOW_DIM` | 1280 | Per-copy Virchow2 embedding dimension |
| `UNI_DIM` | 1536 | Per-copy UNI-v2 embedding dimension |
| `N_REPS` | 3 | Number of copies concatenated per backbone |
| `INPUT_DIM` | 8448 | Model input dimension (`3*1280 + 3*1536`) |
| `FEATURE_DIM` | 3000 | Dimension of the `feats` bottleneck / MELD embedding |
| `PRETRAINED_WEIGHTS` | dict | `{"comet": <path>, "ccdi": <path>}` resolved absolute paths |

### `meld.create_architecture(input_dim=INPUT_DIM) -> keras.Model`

Builds the MELD network with freshly initialized (random) weights. Use this to inspect the
graph, to train a model from scratch (§ `meld.train`), or as the target for `load_weights`
when loading a checkpoint manually.

```python
model = meld.create_architecture()
model.summary()
```

### `meld.load_model(weights="comet", input_dim=INPUT_DIM) -> keras.Model`

Creates the architecture and loads a checkpoint's weights into it.

- `weights`: `"comet"` or `"ccdi"` to use a bundled checkpoint, or a filesystem path to any
  `.weights.h5` file produced by `model.save_weights(...)` from this same architecture
  (e.g. your own fine-tuned model, or a checkpoint downloaded separately from Hugging Face).
- `input_dim`: only change if you trained a variant with a different input layout.
- Returns an **uncompiled** model. Raises `FileNotFoundError` with a clear message if the
  resolved weights path does not exist.

```python
model = meld.load_model("comet")
model = meld.load_model("ccdi")
model = meld.load_model("/path/to/my_finetuned.weights.h5")
```

### `meld.set_inference_mode(model) -> keras.Model`

Sets `model.trainable = False` (disables dropout/BN updates and gradient tracking) and
returns the same model object (mutated in place) for convenient chaining. Call this before
running inference if you plan to reuse the model object only for inference:

```python
model = meld.set_inference_mode(meld.load_model("comet"))
```

### `meld.extract_features(model, x, batch_size=512) -> np.ndarray`

Runs tiles through the model and returns the MELD embedding: the `feats` layer's output,
normalized per-sample (zero mean, unit std across the 3000 features), matching how the
embeddings were produced and consumed during training/evaluation.

- `x`: array-like of shape `(n_tiles, INPUT_DIM)`.
- Returns: `np.ndarray` of shape `(n_tiles, FEATURE_DIM)` i.e. `(n_tiles, 3000)`.

```python
features = meld.extract_features(model, tile_features)  # (n_tiles, 8448) -> (n_tiles, 3000)
```

### `meld.train(model, data, epochs=20, steps_per_epoch=None, learning_rate=1e-4, callbacks=None, optimizer=None) -> keras.Model`

Compiles `model` with `cosine_similarity_loss` and a default `Lamb` optimizer (cosine-decay
schedule with warmup, EMA), then calls `model.fit(data, epochs=epochs, callbacks=callbacks)`.
Use this both to train a freshly created model from scratch, and (equivalently, just with
different default hyperparameters — see `finetune` below) to continue training a loaded one.

- `data`: anything accepted by `keras.Model.fit`, yielding `(X, X)` pairs of shape
  `(batch, INPUT_DIM)` — e.g. a `tf.data.Dataset`, a `keras.utils.Sequence`, or a raw
  `(X, X)` numpy tuple.
- `steps_per_epoch`: only used to size the default LR schedule when `data` has no `__len__`
  (e.g. a plain generator); ignored once `optimizer` is supplied.
- `optimizer`: pass your own `keras.optimizers.Optimizer` to fully override the default.
- Returns the same `model` object, compiled and with updated weights.

```python
model = meld.create_architecture()
model = meld.train(model, train_dataset, epochs=20)
```

### `meld.finetune(model, data, epochs=5, steps_per_epoch=None, learning_rate=1e-5, callbacks=None, optimizer=None) -> keras.Model`

Identical implementation to `train` (in fact it calls `train` directly), but with defaults
suited to continuing a pretrained checkpoint: fewer epochs and a 10x smaller learning rate.

```python
model = meld.load_model("comet")
model = meld.finetune(model, my_dataset, epochs=5)
```

### `meld.cosine_similarity_loss(y_true, y_pred)`

The Keras loss function used by `train`/`finetune`. Exposed for users who want to `compile`
the model themselves with custom optimizer/callback logic instead of using `train`/`finetune`.
See [Training objective](#training-objective-informational) above for a summary of what it
computes.

## 7. Weights: format, size, and Hugging Face

Pretrained weights are distributed **weights-only**, using Keras's `.weights.h5` format
(no architecture graph, no optimizer state). This is why they are much smaller than the
original training checkpoints:

| File | This release (weights-only) | Original training checkpoint |
|---|---|---|
| `meld_comet.weights.h5` | ~601 MB | ~2.5 GB (`.keras`, includes Lamb optimizer state) |
| `meld_ccdi.weights.h5` | ~601 MB | ~2.5 GB (`.keras`, includes Lamb optimizer state) |

If you fine-tune a model and want to share just the resulting weights the same way:

```python
model.save_weights("my_model.weights.h5")
```

If this repository is cloned from a location that excludes the (still fairly large) weight
files, download them from the Hugging Face model repository (link to be added here once
published) and place them at `weights/meld_comet.weights.h5` and
`weights/meld_ccdi.weights.h5`, or point `meld.load_model` at wherever you saved them.

## 8. Troubleshooting

**`ValueError: Layer 'value' expected 2 variables, but received 0 variables during loading`**
(or similar "expected N variables" errors) when calling `meld.load_model(...)`.

This means you're running under the legacy Keras 2.x that ships with TensorFlow <= 2.15. The
released `.weights.h5` files require **Keras 3** (bundled by default with `tensorflow >= 2.16`).
Upgrade TensorFlow, e.g.:

```bash
pip install "tensorflow>=2.16"
```

**`FileNotFoundError: Could not find MELD weights at ...`**

Either pass `"comet"` / `"ccdi"` (bundled checkpoints, expected under `weights/`), or a valid
path to your own `.weights.h5` file. If you cloned the repo without the large weight files,
see [§7](#7-weights-format-size-and-hugging-face) for where to download them.

**Predictions look wrong / shapes don't match after editing `create_architecture`.**

`load_weights` matches variables by layer name, and layer names are auto-generated based on
*creation order* within the function. If you add, remove, or reorder any `Dense` /
`MultiHeadAttention` / `LayerNormalization` layer in `create_architecture`, existing
`.weights.h5` files will silently fail to load correctly (or raise a shape/count mismatch).
Keep the layer-creation order intact if you modify this file, or retrain from scratch.

## 9. Status and scope

This is a first, minimal public release intended to accompany a paper currently under
revision. It intentionally excludes:
- Whole-slide image tiling and feature-caching pipelines (Virchow2/UNI-v2 extraction) for
  COMET, CCDI, and COG.
- Downstream task heads (e.g. tumor detection, subtype classification) built on top of MELD
  embeddings.
- Extensive hyperparameter documentation for the original training runs.

A more comprehensive release, including these pieces and a full usage guide, is planned once
the paper is published.
