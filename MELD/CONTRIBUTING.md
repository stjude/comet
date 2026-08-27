# Contributing

Thanks for your interest in this project. This repository is a minimal public release
accompanying a paper currently under revision, so its primary purpose is transparency and
reproducibility rather than active feature development — but bug reports, clarifications, and
small fixes are welcome.

## Reporting issues

Please open a GitHub issue and include:
- Which function/module is affected (`meld.load_model`, `meld.extract_features`, etc.).
- The exact code you ran and the full error output.
- Your Python, TensorFlow/Keras, and package versions (`pip freeze`).
- Whether you hit the issue with the `comet` or `ccdi` checkpoint.

## Proposing changes

1. Fork the repository and create a branch for your change.
2. Keep pull requests focused and small — one logical change per PR.
3. Follow the existing code style in the file you're editing; see
   [`docs/DOCUMENTATION.md`](docs/DOCUMENTATION.md) for the architecture and API conventions
   used throughout this codebase.
4. Update the relevant section of [`docs/DOCUMENTATION.md`](docs/DOCUMENTATION.md) and/or
   `README.md` if your change affects behavior, inputs, outputs, or the package layout
   described there.
5. Open a pull request describing the motivation for the change and how you tested it.

## Running things locally

```bash
pip install -e .
# or: pip install -r requirements.txt
```

Requires `tensorflow >= 2.16` (Keras 3) — see
[Troubleshooting](docs/DOCUMENTATION.md#8-troubleshooting) if you hit weight-loading errors on
an older TensorFlow/Keras. Since this repository does not include a live test suite, the
recommended way to validate a change is to run
[`notebooks/quickstart.ipynb`](notebooks/quickstart.ipynb) end-to-end against both the `comet`
and `ccdi` checkpoints and confirm `extract_features` still returns the expected `(n_tiles,
3000)` output shape.

## Code of Conduct

Participation in this project is governed by our [Code of Conduct](CODE_OF_CONDUCT.md).
