"""
MELD architecture definition.

MELD is a small transformer ("ViT-style") encoder that fuses tile-level
foundation-model embeddings from UNI-v2 (1536-d) and Virchow2 (1280-d) into
a single compact representation per tile.

Expected input layout (one row per tile), concatenated along the last axis:

    [ virchow2_rep_0 (1280), virchow2_rep_1 (1280), virchow2_rep_2 (1280),
      uni_v2_rep_0   (1536), uni_v2_rep_1   (1536), uni_v2_rep_2   (1536) ]

    -> total input dimension = 3*1280 + 3*1536 = 8448

The "3" copies of each backbone's embedding correspond to how COMET/CCDI/COG
tiles were cached (e.g. multiple crops/scales per tile). If your pipeline only
has a single embedding per backbone, repeat/tile it 3x to match this layout,
or retrain/fine-tune (see meld.train) on your own data layout.

The model is trained as a masked-fusion autoencoder: it reconstructs the
concatenated input plus two auxiliary per-backbone projections, while an
intermediate "feats" layer (3000-d) serves as the actual MELD embedding used
downstream. This module only defines the network; see meld.model for loading
pretrained weights and meld.train for (re)training.
"""

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, regularizers

VIRCHOW_DIM = 1280
UNI_DIM = 1536
N_REPS = 3

INPUT_DIM = N_REPS * VIRCHOW_DIM + N_REPS * UNI_DIM  # 8448
FEATURE_DIM = 3000

# Split points of the 11264-d model output: [reconstruction | virchow-branch | uni-branch]
_DIM_V = N_REPS * VIRCHOW_DIM          # 3840
_DIM_RECON = INPUT_DIM                  # 8448
_DIM_V_BRANCH = VIRCHOW_DIM             # 1280
_DIM_U_BRANCH = UNI_DIM                 # 1536


def anorm(x):
    """Per-sample mean/std normalization, used both inside the graph (on the
    'feats' output) and as the final MELD embedding normalization."""
    xm = tf.math.reduce_mean(x, axis=1, keepdims=True)
    xs = tf.math.reduce_std(x, axis=1, keepdims=True)
    return (x - xm) / xs


def _mlp(x, hidden_units, dropout_rate, l1, l2):
    for units in hidden_units:
        x = layers.Dense(
            units,
            activation=keras.activations.gelu,
            kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2),
            kernel_initializer="he_normal",
        )(x)
        x = layers.Dropout(dropout_rate)(x)
    return x


def create_architecture(input_dim=INPUT_DIM):
    """Build the MELD encoder/decoder graph (uninitialized/random weights).

    NOTE: layers below are created in the exact same order used to train the
    released checkpoints. This matters because pretrained weights are shipped
    weights-only (see meld.model.load_model) and are matched to layers by
    their (auto-generated) names, which depend on creation order.
    """
    R = 1e-6
    S = 0.01
    l1, l2 = S * R, R

    inputs = keras.Input(shape=(input_dim,))
    x1 = layers.Dense(6000, activation="linear", kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2), kernel_initializer="he_normal")(inputs)
    x2 = layers.Reshape((100, 60))(x1)

    # ViT block 1
    x3 = layers.MultiHeadAttention(8, 60, dropout=0.1, kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2), kernel_initializer="he_normal")(x2, x2)
    x4 = layers.LayerNormalization(epsilon=1e-6)(x3)
    x5 = _mlp(x4, [360, 60], 0.1, l1, l2)
    x6 = layers.Add()([x4, x5])

    # ViT block 2
    x6 = layers.MultiHeadAttention(8, 60, dropout=0.1, kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2), kernel_initializer="he_normal")(x6, x6)
    x7 = layers.LayerNormalization(epsilon=1e-6)(x6)
    x8 = _mlp(x7, [360, 60], 0.1, l1, l2)
    x9 = layers.Add()([x7, x8])

    x9f = layers.Flatten()(x9)
    feats = layers.Dense(
        FEATURE_DIM,
        activation="linear",
        kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2),
        kernel_initializer="he_normal",
        name="feats",
    )(x9f)
    feats_norm = layers.Lambda(anorm)(feats)
    x9n = layers.Reshape((100, 30))(feats_norm)

    # ViT block 3
    x10 = layers.MultiHeadAttention(8, 30, dropout=0.1, kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2), kernel_initializer="he_normal")(x9n, x9n)
    x11 = layers.LayerNormalization(epsilon=1e-6)(x10)
    x12 = _mlp(x11, [360, 60], 0.1, l1, l2)

    # ViT block 4
    x14 = layers.MultiHeadAttention(8, 60, dropout=0.1, kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2), kernel_initializer="he_normal")(x12, x12)
    x15 = layers.LayerNormalization(epsilon=1e-6)(x14)
    x16 = _mlp(x15, [360, 60], 0.1, l1, l2)
    x17 = layers.Add()([x15, x16])

    # Reconstruction head
    xo1 = layers.Flatten()(x17)
    x_recon = layers.Dense(input_dim, activation="linear", kernel_regularizer=regularizers.L1L2(l1=l1, l2=l2), kernel_initializer="he_normal")(xo1)

    # Auxiliary per-backbone projection heads (used only for the training loss)
    xuv = layers.Flatten()(x9n)
    x_virchow = _mlp(xuv, [3 * VIRCHOW_DIM, VIRCHOW_DIM], 0.2, l1, l2)
    x_uni = _mlp(xuv, [3 * UNI_DIM, UNI_DIM], 0.2, l1, l2)

    output = layers.Concatenate()([x_recon, x_virchow, x_uni])
    return keras.Model(inputs, output)


def cosine_similarity_loss(y_true, y_pred):
    """Training/fine-tuning loss for MELD.

    y_true: (batch, INPUT_DIM) the raw concatenated Virchow2/UNI-v2 input.
    y_pred: (batch, INPUT_DIM + VIRCHOW_DIM + UNI_DIM) model output, i.e.
            [reconstruction, virchow_projection, uni_projection].
    """
    dim0, dim1 = _DIM_V, _DIM_RECON
    dim2 = dim1 + _DIM_V_BRANCH

    y_true_0, y_pred_0 = y_true[:, :dim0], y_pred[:, :dim0]
    y_true_1, y_pred_1 = y_true[:, dim0:], y_pred[:, dim0:dim1]

    def _cos_sim(a, b):
        a_n = tf.linalg.normalize(a, axis=-1)[0]
        b_n = tf.linalg.normalize(b, axis=-1)[0]
        return tf.reduce_sum(a_n * b_n, axis=-1)

    cosine_sim_0 = _cos_sim(y_true_0, y_pred_0)
    cosine_sim_1 = _cos_sim(y_true_1, y_pred_1)

    abserr = tf.abs(y_true - y_pred[:, :dim1])
    l1loss = tf.reduce_mean(abserr[:, :dim0], axis=1) + tf.reduce_mean(abserr[:, dim0:], axis=1)
    lmloss = tf.reduce_max(abserr[:, :dim0], axis=1) + tf.reduce_max(abserr[:, dim0:], axis=1)

    def _pairwise_l1_dist_loss(a, b):
        a_e0, a_e1 = tf.expand_dims(a, 0), tf.expand_dims(a, 1)
        a_dist = tf.reduce_sum(tf.abs(a_e0 - a_e1), axis=2)
        a_dist = a_dist / tf.expand_dims(tf.reduce_max(a_dist, axis=1), axis=1)

        b_e0, b_e1 = tf.expand_dims(b, 0), tf.expand_dims(b, 1)
        b_dist = tf.reduce_sum(tf.abs(b_e0 - b_e1), axis=2)
        b_dist = b_dist / tf.expand_dims(tf.reduce_max(b_dist, axis=1), axis=1)
        return tf.reduce_mean(tf.abs(a_dist - b_dist))

    def _pairwise_cos_sim_loss(a, b):
        a_n = tf.linalg.normalize(a, axis=1)[0]
        b_n = tf.linalg.normalize(b, axis=1)[0]
        a_cs = tf.matmul(a_n, a_n, transpose_b=True)
        b_cs = tf.matmul(b_n, b_n, transpose_b=True)
        return tf.reduce_mean(tf.abs(a_cs - b_cs))

    dist_loss_0 = _pairwise_l1_dist_loss(y_true_0, y_pred_0)
    dist_loss_1 = _pairwise_l1_dist_loss(y_true_1, y_pred_1)
    cos_dist_loss_0 = _pairwise_cos_sim_loss(y_true_0, y_pred_0)
    cos_dist_loss_1 = _pairwise_cos_sim_loss(y_true_1, y_pred_1)

    rpred_0 = y_pred[:, dim1:dim2]
    rpred_1 = y_pred[:, dim2:]
    rdist_loss_0 = _pairwise_l1_dist_loss(y_true_0, rpred_0)
    rdist_loss_1 = _pairwise_l1_dist_loss(y_true_1, rpred_1)
    cos_rdist_loss_0 = _pairwise_cos_sim_loss(y_true_0, rpred_0)
    cos_rdist_loss_1 = _pairwise_cos_sim_loss(y_true_1, rpred_1)

    loss = (
        0.025 * l1loss + 0.025 * lmloss
        + 2 - cosine_sim_0 - cosine_sim_1
        + 0.25 * dist_loss_0 + 0.25 * dist_loss_1
        + rdist_loss_0 + rdist_loss_1
        + 0.25 * cos_dist_loss_0 + 0.25 * cos_dist_loss_1
        + cos_rdist_loss_0 + cos_rdist_loss_1
    )
    return loss
