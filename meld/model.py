"""
Loading MELD models and running inference.

Pretrained weights are distributed weights-only (``.weights.h5``, no
optimizer state) to keep the release small. The full computation graph is
rebuilt from ``meld.architecture.create_architecture`` and the weights are
then loaded into it, so no large architecture file is needed.
"""

import os

import numpy as np
from tensorflow import keras

from .architecture import INPUT_DIM, anorm, create_architecture

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_WEIGHTS_DIR = os.path.join(os.path.dirname(_THIS_DIR), "weights")

# Convenience aliases for the two released checkpoints (see README for how to
# download these from Hugging Face if this ``weights/`` folder is empty).
PRETRAINED_WEIGHTS = {
    "comet": os.path.join(_WEIGHTS_DIR, "meld_comet.weights.h5"),   # MELD-COMET, 20 epochs
    "ccdi": os.path.join(_WEIGHTS_DIR, "meld_ccdi.weights.h5"),     # MELD-CCDI, 20 epochs
}


def load_model(weights="comet", input_dim=INPUT_DIM):
    """Create the MELD architecture and load pretrained (or custom) weights.

    Args:
        weights: one of "comet" / "ccdi" (bundled checkpoints), or a path to
            any ``.weights.h5`` file saved with ``model.save_weights(...)``
            from this same architecture (e.g. your own fine-tuned model).
        input_dim: input dimension of the fusion network (default 8448 =
            3x Virchow2 (1280) + 3x UNI-v2 (1536)). Only change this if you
            retrain from scratch with a different tile-feature layout.

    Returns:
        An uncompiled ``keras.Model``. Use ``set_inference_mode`` before
        inference, or ``meld.train`` to (further) train it.
    """
    weights_path = PRETRAINED_WEIGHTS.get(weights, weights)
    if not os.path.isfile(weights_path):
        raise FileNotFoundError(
            f"Could not find MELD weights at '{weights_path}'. "
            "Pass 'comet', 'ccdi', or a path to a .weights.h5 file. "
            "See the README for downloading the released checkpoints."
        )
    model = create_architecture(input_dim)
    model.load_weights(weights_path)
    return model


def set_inference_mode(model):
    """Freeze the model for inference (no dropout, no trainable updates).

    Use ``model.predict(x)`` or ``model(x, training=False)`` afterwards.
    """
    model.trainable = False
    return model


def extract_features(model, x, batch_size=512):
    """Run tiles through MELD and return the normalized MELD embedding.

    Args:
        model: a MELD model, e.g. from ``load_model``.
        x: array of shape (n_tiles, INPUT_DIM) with the concatenated
            Virchow2/UNI-v2 tile representations (see meld.architecture
            module docstring for the exact expected layout).
        batch_size: prediction batch size.

    Returns:
        np.ndarray of shape (n_tiles, FEATURE_DIM=3000): the per-tile MELD
        embedding (the model's 'feats' layer, normalized per-sample).
    """
    feats_model = keras.Model(inputs=model.input, outputs=model.get_layer("feats").output)
    raw = feats_model.predict(np.asarray(x, dtype=np.float32), batch_size=batch_size, verbose=0)
    return anorm(raw).numpy()
