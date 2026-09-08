# Example data

## `samplesheet_epic.csv` — runnable

Four EPIC (850K) samples from GEO series
[GSE130295](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130295),
public data redistributed here as the project's smoke test. The IDATs are
present, so this samplesheet runs end-to-end — see the quick-start command in
[`../../README.md`](../../README.md).

## `samplesheet_450k.csv`, `samplesheet_epicv2.csv` — format illustration only

These carry **synthetic identifiers and no IDAT files**. They exist to show the
samplesheet layout for the other two supported platforms; running them will fail
at the IDAT-reading step, which is expected.

To use them, point `Basename` at your own IDAT path stems (no `_Grn.idat` /
`_Red.idat` suffix) and replace the `Sentrix_ID` values with your real chip
barcodes.

The 450K and EPICv2 arrays previously bundled here were internal patient samples
and were removed: an Illumina methylation array carries rs genotyping probes and
is therefore identifiable, so raw IDATs cannot be de-identified and are not
redistributable. Only public GEO data ships in this repository.

## Required columns

| Column | Required | Description |
|---|---|---|
| `Sentrix_ID` | yes | Unique per-sample identifier |
| `Sample_Name` | yes | Human-readable label used in plots |
| `Basename` | yes | IDAT path stem, absolute or relative to `pipeline/` |
| anything else | no | Carried through and available for colouring plots (e.g. `Sample_Group`) |
