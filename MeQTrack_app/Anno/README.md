# Anno — sesame annotation assets (vendored)

Annotation files used by the sesame-based QC steps. These are vendored
(committed) rather than downloaded at runtime so the pipeline runs offline. The
clock models are tiny (a few KB); the EPICv2→EPIC map is slimmed to 4.3 MB from
a 115 MB source (the raw source stays gitignored — see below).

## Source

Downloaded from the Zhou Lab Infinium annotation repository:
https://github.com/zhou-lab/InfiniumAnnotationV1/tree/main/Anno

Mirror the upstream directory layout (`Anno/<platform>/<file>`) so paths are
predictable.

## Contents

| File | Use |
|---|---|
| `HM450/Clock_Horvath353.rds` | Horvath 353-CpG epigenetic clock (Horvath 2013), `predictAge()`-compatible model. Used for the `Horvath_Age` QC column across all array types (probe IDs are base `cg` form, matching our collapsed beta matrices). |
| `EPICv2/Clock_Horvath353.EPICv2.345.rds` | EPICv2-native Horvath clock (suffixed probe IDs). Kept for a future EPICv2-specific age path; not used yet because our EPICv2 betas are collapsed to base CpG IDs. |
| `EPICv2/EPICv2ToEPIC_map.tsv.gz` | **Committed.** Slim 3-column map (`ID_EPIC1`, `ID_EPIC2`, `big_delta`) derived from the upstream `EPICv2ToEPIC_conversion.tsv`. Used to convert collapsed EPICv2 betas into EPIC space so the EPIC-only `estimateLeukocyte` works for EPICv2. |
| `EPICv2/EPICv2ToEPIC_conversion.tsv` | **Gitignored** (~115 MB). The full upstream benchmark table (adds six per-cell-line β-delta columns we don't need). Only the slim map above is committed; re-derive with `cut -f1,2,9 ... \| gzip`. |

## Notes

- The age-clock models are the structured `predictAge()` format
  (`intercept` / `param$slope` / `response2age`). The legacy coefficient table
  bundled in sesameData 1.30.0 (`age.inference`) is NOT compatible with the
  installed `predictAge()`, which is why these model files are vendored here.
- Leukocyte-fraction estimation (`estimateLeukocyte`) does NOT need a file from
  here — its reference (`leukocyte.betas`) ships in sesameData and is cached
  via ExperimentHub.

To refresh, re-download the same paths from the upstream repo `main` branch.
