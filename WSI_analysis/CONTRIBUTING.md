# Contributing

Thanks for your interest in this project. This repository accompanies a research
publication, so its primary purpose is transparency and reproducibility rather
than active feature development — but bug reports, clarifications, and small
fixes are welcome.

## Reporting issues

Please open a GitHub issue and include:
- Which script/folder is affected.
- The exact command you ran and the full error output.
- Your Python version and installed package versions (`pip freeze`).

## Proposing changes

1. Fork the repository and create a branch for your change.
2. Keep pull requests focused and small — one logical change per PR.
3. Follow the existing code style in the file you're editing (this codebase
   intentionally favors explicit, readable code over cleverness or heavy
   abstraction; see each folder's `README.md` for the conventions used there).
4. Update the relevant `README.md` if your change affects behavior, inputs,
   outputs, or file layout described there.
5. Open a pull request describing the motivation for the change and how you
   tested it.

## Running things locally

```bash
pip install -r requirements.txt
```

Every task folder's `README.md` (`tumor_classification/`, `RMSsubtyping/`,
`NB_MYCN/`, `NB_ADRN_MES/`, `RMS_survival/`) has copy-pasteable example commands
that exercise the code end-to-end against the data shipped in [`data/`](data/README.md).
Since this repository does not include a live test suite, the recommended way
to validate a change is to re-run the relevant task's example command(s) with
`--n-splits` set to a small number and confirm the output still has the expected
shape and reasonable metric values.

## Code of Conduct

Participation in this project is governed by our [Code of Conduct](CODE_OF_CONDUCT.md).
