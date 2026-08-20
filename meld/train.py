"""
(Further) training a MELD model.

MELD is trained as a self-supervised fusion autoencoder: given a batch of
tile-level (Virchow2, UNI-v2) representations ``X`` (shape (batch,
INPUT_DIM)), it is trained with ``model.fit(data, ...)`` where ``data``
yields ``(X, X)`` pairs (input == target); the loss in
``meld.architecture.cosine_similarity_loss`` handles the rest.

Both ``train`` (fresh/full training) and ``finetune`` (continuing an already
loaded/pretrained model) are thin wrappers around the same compile+fit call;
they only differ in their default optimizer settings.
"""

from tensorflow import keras

from .architecture import cosine_similarity_loss


def _default_optimizer(learning_rate, total_steps):
    schedule = keras.optimizers.schedules.CosineDecay(
        learning_rate, total_steps, alpha=0.5, warmup_target=learning_rate * 10, warmup_steps=max(1, total_steps // 10)
    )
    return keras.optimizers.Lamb(learning_rate=schedule, use_ema=True, ema_momentum=0.99, clipvalue=3)


def train(model, data, epochs=20, steps_per_epoch=None, learning_rate=1e-4, callbacks=None, optimizer=None):
    """Compile and train a MELD model (from scratch or from loaded weights).

    Args:
        model: a keras.Model built with meld.architecture.create_architecture
            (optionally pre-loaded with weights via meld.model.load_model).
        data: anything accepted by keras.Model.fit as the first positional
            argument (a tf.data.Dataset, a keras.utils.Sequence/generator, or
            a (X, X) numpy tuple), yielding (input, target) pairs of shape
            (batch, INPUT_DIM).
        epochs: number of epochs to train.
        steps_per_epoch: needed for the default LR schedule if `data` doesn't
            expose `__len__` (e.g. a raw generator); ignored if `optimizer`
            is supplied explicitly.
        learning_rate: peak learning rate for the default optimizer.
        callbacks: optional list of keras callbacks (e.g. ModelCheckpoint).
        optimizer: optional custom optimizer; overrides `learning_rate`.

    Returns:
        The trained model (compiled, weights updated in place).
    """
    if optimizer is None:
        total_steps = steps_per_epoch or (len(data) if hasattr(data, "__len__") else 10000)
        optimizer = _default_optimizer(learning_rate, total_steps * epochs)
    model.compile(optimizer=optimizer, loss=cosine_similarity_loss)
    model.fit(data, epochs=epochs, callbacks=callbacks)
    return model


def finetune(model, data, epochs=5, steps_per_epoch=None, learning_rate=1e-5, callbacks=None, optimizer=None):
    """Continue training an already-loaded/pretrained MELD model.

    Identical to `train`, but defaults to fewer epochs and a smaller learning
    rate, appropriate for adapting a pretrained checkpoint to new data.
    """
    return train(
        model, data, epochs=epochs, steps_per_epoch=steps_per_epoch,
        learning_rate=learning_rate, callbacks=callbacks, optimizer=optimizer,
    )
