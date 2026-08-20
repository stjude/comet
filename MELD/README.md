# MELD

MELD is a fusion model that combines multi-scale tile-level **UNI-v2** (3*1536-d) and **Virchow-v2** (3*1280-d)
foundation-model embeddings into a compact per-tile representation for computational
pathology. This repository is a minimal, working release of the code and weights needed to
load. Note for the Virchow-v2 model, mean and std of input tiles should be normalized to 0.5 instead of default values.

Two checkpoints are included and can be downloaded from hugging face (https://huggingface.co/alifp04/MELD/tree/main)

| Name | Description |
|---|---|
| `comet` (MELD-COMET) | Trained/distilled on the COMET dataset |
| `ccdi` (MELD-CCDI) | Trained/distilled on the CCDI dataset |

📖 **Full documentation:** [`docs/DOCUMENTATION.md`](docs/DOCUMENTATION.md) (architecture
details, complete API reference, data format spec, troubleshooting).
📓 **Runnable walkthrough:** [`notebooks/quickstart.ipynb`](notebooks/quickstart.ipynb).

## Install

```bash
pip install -e .
# or: pip install -r requirements.txt
```

Requires `tensorflow >= 2.16` (Keras 3) — see [Troubleshooting](docs/DOCUMENTATION.md#8-troubleshooting)
if you hit weight-loading errors on an older TensorFlow/Keras.

## Quickstart

```python
import meld

model = meld.load_model("comet")          # or "ccdi", or a custom weights path
model = meld.set_inference_mode(model)

features = meld.extract_features(model, tile_features)  # (n_tiles, 8448) -> (n_tiles, 3000)
```

Other core functions: `meld.create_architecture()`, `meld.train(model, data)`,
`meld.finetune(model, data)`. See the [full API reference](docs/DOCUMENTATION.md#6-api-reference)
for details and signatures.

## Weights

Pretrained weights are stored **weights-only** (`.weights.h5`, no optimizer state) under
`weights/`, which keeps them a fraction of the size of the original training checkpoints
(~601 MB vs. ~2.5 GB each). See [§7 of the docs](docs/DOCUMENTATION.md#7-weights-format-size-and-hugging-face)
for details, including where to download them from Hugging Face if this repo was cloned
without the large files.

## Expected input format

MELD expects one row per tile of shape `(N_TILES, 8448)`:

```
[ virchow2_a (1280), virchow2_b (1280), virchow2_c (1280),
  uni_v2_a   (1536), uni_v2_b   (1536), uni_v2_c   (1536) ]
```

See [§5 of the docs](docs/DOCUMENTATION.md#5-data-format) for the full explanation and
guidance on adapting your own data layout.

## Package layout

```
meld/
  architecture.py   # create_architecture(), training loss
  model.py          # load_model(), set_inference_mode(), extract_features()
  train.py          # train(), finetune()
weights/
  meld_comet.weights.h5
  meld_ccdi.weights.h5
notebooks/
  quickstart.ipynb
docs/
  DOCUMENTATION.md
```

## Status

This is a first, minimal public release to accompany a paper currently under revision. The
API intentionally favors simplicity over completeness; expect a more comprehensive release
(with full tiling/caching utilities and documentation) once the paper is published.
