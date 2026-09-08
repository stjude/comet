# Contributing

Thanks for your interest in this project. This repository is a public release accompanying the
COMET project, so its primary purpose is transparency and reproducibility rather than active
feature development — but bug reports, clarifications, and small fixes are welcome.

## Reporting issues

Please open a GitHub issue and include:
- Which pipeline step is affected (`preprocess`, `qc`, `filtering`, `dim_reduction`,
  `reference_projection`, `cnv`, `visualization`), or the app tab if it is a UI problem.
- The array type (450K, EPIC, or EPICv2) and the number of samples in the run.
- The exact command you ran, and `pipeline_log.txt` from the affected run directory.
- Your R version and `sessionInfo()`.

**Never attach IDATs, beta matrices, samplesheets, or QC reports from real samples to an
issue.** These carry identifiable data. Reproduce with the bundled public example where you
can, and redact sample identifiers otherwise. If you believe published files contain
identifiable data, follow [`SECURITY.md`](SECURITY.md) instead of opening an issue.

## Proposing changes

1. Fork the repository and create a branch for your change.
2. Keep pull requests focused and small — one logical change per PR.
3. Follow the existing style of the file you are editing. Analysis logic belongs in a module
   under `pipeline/pipeline_modules/`; the app should stay a thin wrapper that talks to the
   pipeline only through its CLI and its output directory.
4. Update [`README.md`](README.md) and [`pipeline/README.md`](pipeline/README.md) if your
   change affects behaviour, parameters, or the output layout they document.
5. If your change alters pipeline behaviour or the run-directory structure, add an ADR under
   [`docs/decisions/`](docs/decisions/).
6. Open a pull request describing the motivation and how you tested it.

## Running things locally

```bash
Rscript setup.R    # provisions CRAN + Bioconductor deps into a project-local renv library
```

Requires R >= 4.4 and pandoc >= 2.x. Since this repository has no automated test suite, the
recommended way to validate a change is to run the bundled four-sample public EPIC example
end-to-end and confirm it still produces a complete run directory:

```bash
Rscript pipeline/methylation_pipeline.R \
  --input      pipeline/data/example/samplesheet_epic.csv \
  --output     pipeline/results_example \
  --data_dir   pipeline/data \
  --array_type EPIC \
  --threads    4 \
  --step       all
```

Reference projection additionally needs the reference beta matrices, which are distributed
separately — see [Data and privacy](README.md#data-and-privacy).

## Code of Conduct

Participation in this project is governed by our [Code of Conduct](CODE_OF_CONDUCT.md).
