# goals_v5 — v2.3.0 release scope: sesame sample/experiment QC

`goals_v4.md` was the v2.0.0 → v2.1.0 release scope (the reference-projection
marquee capability plus the Capper and sarcoma reference datasets). It shipped
and has since been superseded by v2.2.x (pandoc auto-provisioning, renv
notice suppression, reference-projection docs, COMET dataset relabelling).
Current state: `MEQTRACK_VERSION = "2.2.2"` in
`pipeline/methylation_pipeline.R`.

This file targets **v2.3.0** — the next *minor* release, reserved for
**sesame-based sample/experiment QC**:

1. **GCT bisulfite-conversion control score** — the anchor capability; it
   **gates `Pass_QC`** (configurable threshold, exposed in the Settings UI).
2. **Sample-integrity inferences** (informational, for sample-swap / purity
   checks): **predicted sex** (`inferSex`), **Horvath epigenetic age**
   (`predictAge` on a vendored clock model), and **leukocyte fraction**
   (`estimateLeukocyte`, with an EPICv2→EPIC conversion so it works on all
   array types).

All come from the sesame beta/SigDF path. A small vendored `Anno/` directory
(Zhou Lab annotation assets) backs the clock model and the EPICv2→EPIC map.
Resist folding unrelated work into this scope.

**Status — Phase 1 RELEASED-READY (P1-GATE passed); Phase 2 IN PROGRESS.** Work
is on branch `preprocess-sesame-migration`. Phase 1 (EPIC/450k GCT + all
sample-integrity inferences, including EPICv2 leukocyte) is verified end-to-end
and **P1-GATE passed** (live-app Conversion QC tab + pink fail-rows + banner
confirmed; a fail-row highlighting bug under bslib Bootstrap 5 was fixed).
**Phase 2** now covers more than the original EPICv2-GCT gap: at the user's
request it also adds **sex karyotype** and a **SNP genotype fingerprint**
(replacing the removed `inferEthnicity`). EPICv2 GCT is done and committed;
karyotype + fingerprint are implemented. **Phase 3** = packaging/release.

---

## Bisulfite-conversion QC

**Vision:** report a per-sample **GCT score** (Zhou et al. 2017) so users can
judge how complete bisulfite conversion was on each array. Infinium platforms
are intrinsically robust to *incomplete* conversion — non-converted probes fail
to hybridize — but residual incomplete conversion can be *quantified* from the
Infinium-I C/T-extension probes. A score near **1.0** means complete
conversion; **higher** values indicate more residual incomplete conversion.

This is the specific QC metric **sesame provides that minfi does not**, and the
reason preprocessing reads IDATs through sesame at all. The pipeline already
derived betas and pOOBAH detection-p via `sesame::openSesame()`; this capability
adds `sesame::bisConversionControl()` alongside them.

The metric **gates `Pass_QC`**: a sample whose GCT exceeds
`config$qc$max_gct_score` (default **1.3**) **fails** QC, with a
`Failure_Reason` and a dedicated `Flag_GCT` column, exactly like the
detection-p and failed-probe checks — so it also triggers sample removal when
`filter_failed_samples` is on. Samples with an **NA** GCT (e.g. EPICv2 in
Phase 1, where GCT is not yet computed) are **never** failed on it — an
unmeasured score is not a failing score.

## What this is NOT

This is **not** a normalization-method change. Investigation during planning
confirmed the `normalization` parameter (`raw`/`illumina`/`functional`/
`quantile`/`swan`/`sesame`) was already completely inert end-to-end: the app
never set it, the pipeline hardcoded `prep = "QCDB"` for betas regardless, and
the returned `result$normalization` was never read. Betas already come from
sesame. v2.3.0 leaves the beta and minfi paths untouched and *adds* the GCT
output. (The misleading "normalization is configurable" line in the README was
corrected as part of this work; the inert `normalization` *parameter* itself is
left for a later cleanup, not v2.3.0 scope.)

## API facts — validated against installed sesame 1.30.0

- `bisConversionControl(sdf, extR = NULL, extA = NULL, verbose = FALSE)` takes a
  sesame **SigDF** and returns one numeric GCT score per sample.
- Internally it reads `InfIR(sdf)` (Infinium-I probes), so the SigDF must have
  its Infinium-I channel inferred first → the minimal prep is **`prep = "C"`**
  (`inferInfiniumIChannel`). Heavier prep (noob/dyebias) would distort the raw
  extension-probe signal and is wrong here.
- For `platform %in% c("EPICplus", "EPIC", "HM450")` it **auto-fetches** the
  extension probes via `sesameDataGet(paste0(platform, ".probeInfo"))`.
- For any other platform (EPICv2, MSA, HM27) it hits
  `stopifnot(!is.null(extR) && !is.null(extA))` and **errors** unless extR/extA
  are supplied manually. This is the EPICv2 gap that splits Phase 1 from
  Phase 2.
- `openSesame(idats, prep = "C", func = bisConversionControl, BPPARAM = ...)`
  returns the per-sample score vector directly — reusing the same IDAT
  discovery and the LSF-safe `bpparam` already built in `preprocess.R`.
- GCT is **not** part of `sesameQC` / `sesameQC_calcStats`; it must be called
  directly.

## Outputs — standalone table + QC gate

GCT surfaces in two places:

1. **Standalone `qc/conversion_qc.csv`** (`Sample_ID, GCT_Score, Array_Type,
   Note`), written by the preprocess step. It is computed on a different
   (sesame SigDF) code path than the main QC table, and EPICv2 rows are
   explicitly `NA` + note, which reads cleanly in a dedicated table. The app
   surfaces it as a "Conversion QC" tab next to "Sample metrics".
2. **`sample_qc_report.csv` gate**: the preprocess GCT table is passed into
   `perform_qc()`, which merges `GCT_Score` and a `Flag_GCT` column into the
   main QC table and folds `Flag_GCT` into `Pass_QC` / `Failure_Reason`. The
   threshold is `config$qc$max_gct_score` (default 1.3, exposed in the Settings
   UI). NA scores never fail.

## Sample-integrity inferences (sex + age + leukocyte fraction)

Three informational columns in `sample_qc_report.csv`, computed in
`perform_qc()` from the beta matrix (so no preprocess threading), for
sample-swap detection and tumour-purity gauging:

- **`Sesame_Sex`** — sesame `inferSex(betas)` → MALE/FEMALE (curated X/Y probe
  model, auto-detects platform). Complements the existing minfi `getSex`
  prediction; a mismatch flags a likely swap.
- **`Horvath_Age`** — Horvath 353-CpG epigenetic age (Horvath 2013) via the
  proper sesame `predictAge(betas, model)`. The model is **vendored** at
  `Anno/HM450/Clock_Horvath353.rds` (from the Zhou Lab InfiniumAnnotation repo)
  in `predictAge()`-compatible form (`intercept`/`param$slope`/`response2age`).
  **Why vendored:** the `age.inference` data shipped in the pinned sesameData
  1.30.0 is the *legacy* coefficient table, incompatible with the installed
  `predictAge()` — so we supply the structured model file instead. The HM450
  clock (base `cg` IDs) is used for all array types, matching our collapsed
  beta matrices. `load_horvath_model()` resolves the path robustly from the
  pipeline working dir.
- **`Leukocyte_Fraction`** — sesame `estimateLeukocyte(betas, platform)` (0–1),
  its reference (`leukocyte.betas`) from sesameData (already cached). Works for
  450k/EPIC directly; for **EPICv2** the betas are first converted to EPIC space
  via the Zhou Lab `EPICv2ToEPIC_map.tsv.gz` (a reliability-filtered probe
  subset — see below), then estimated on the EPIC platform. So EPICv2 now gets
  a real leukocyte fraction instead of NA.

All three are **informational only** — none affects `Pass_QC`. **Not available
in sesame 1.30.0:** `inferSexKaryotypes` (XaY-style karyotype) and
`inferEthnicity`; they would require a sesame upgrade (deferred, see *Open
questions*).

### Vendored annotation — `Anno/`

`Anno/` mirrors the upstream `zhou-lab/InfiniumAnnotationV1/Anno` layout and
holds sesame annotation assets committed to the repo. Currently:
`HM450/Clock_Horvath353.rds` (~6 KB, used), `EPICv2/Clock_Horvath353.EPICv2.345.rds`
(kept for a future EPICv2-native age path — unused now because its suffixed
probe IDs don't match our collapsed betas), and
`EPICv2/EPICv2ToEPIC_map.tsv.gz` (4.3 MB, the slim 3-column EPICv2→EPIC map).
The raw 115 MB `EPICv2ToEPIC_conversion.tsv` is **gitignored** — only the slim
map (derived via `cut -f1,2,9 | gzip`) is committed. In that map the EPICv2
base ID always equals the EPIC1 ID, so converting our (already collapsed)
EPICv2 betas to EPIC space is a reliability-filtered probe subset (drop
`big_delta == TRUE`). `Anno/README.md` records provenance; `build_release.sh`
stages `Anno/` into the release zip.

## Progress

- **Phase 1 — IMPLEMENTED (EPIC/450k), unit-verified, not yet released.**
  - P1-T1. `compute_gct_scores()` helper added to
    `pipeline/pipeline_modules/preprocess.R`: guards on array type (EPIC/450k
    only), calls `openSesame(prep = "C", func = bisConversionControl)`, maps
    scores back to `Sample_ID`, returns the four-column data frame. EPICv2/other
    → `NA` + Phase-2 note. `tryCatch`-wrapped so GCT failure never aborts
    preprocessing. Threaded into `result$gct`. Resolves the concrete array type
    once (the param may be `"auto"`). Also fixed the stale `CDPB`→`QCDB`
    prep-code comment.
  - P1-T2. `methylation_pipeline.R` writes `result$gct` to
    `qc/conversion_qc.csv` in the preprocess step.
  - P1-T3. `app/R/results_loader.R` loads `conversion_qc.csv` into the results
    bundle (`RESULTS_PATHS$conversion_qc` + bundle slot + NULL-guard +
    contract doc).
  - P1-T4. `app/R/qc_module.R` adds the "Conversion QC" nav tab, a generic
    `DT::datatable` render reading `results()$conversion_qc`, and
    `QC_CONVERSION_TOOLTIPS` column help.
  - P1-T5. Docs — `pipeline/README.md` documents the GCT output and corrects
    the misleading "normalization is configurable" claim.
  - P1-T6. QC gating — `qc.R` `perform_qc()` gains `gct` / `max_gct_score`
    params, merges `GCT_Score` + `Flag_GCT` into `sample_qc`, folds `Flag_GCT`
    into `Pass_QC` and `Failure_Reason`, and exposes `gct_failed_samples`.
    `methylation_pipeline.R` threads `result$gct` (from both the fresh-preprocess
    and reload-from-RData branches) and reads `config$qc$max_gct_score`;
    `config.R` defaults it to 1.3. NA GCT never fails; old preprocessed data
    without `$gct` degrades gracefully (gating disabled).
  - P1-T7. Settings UI — `app/R/settings_module.R` adds a "Max GCT" numeric
    input (default 1.3) wired through the param bag, past-run restore, and
    `pipeline_bridge.R` `.write_run_config` (`config$qc$max_gct_score`).
  - P1-T8. Sample-integrity inferences — `qc.R` `compute_sample_inferences()`
    (sesame `inferSex` + Horvath age + leukocyte fraction), called in
    `perform_qc()` to add `Sesame_Sex` / `Horvath_Age` / `Leukocyte_Fraction`
    columns; `qc_module.R` tooltips. Informational, never gates `Pass_QC`.
  - P1-T9. Vendored clock model + proper `predictAge()` — committed
    `Anno/HM450/Clock_Horvath353.rds` (+ EPICv2 variant, unused) from the Zhou
    Lab repo; `load_horvath_model()` reads it and age now uses
    `sesame::predictAge()` (replacing the earlier manual Horvath math, which
    matched: 84.1y). Added leukocyte fraction via `estimateLeukocyte()`
    (EPIC/450k; NA EPICv2). `build_release.sh` bundles `Anno/`. Verified end-to-
    end from the pipeline working dir: model loads, age 84.1y, leukocyte 0.20.
  - P1-T10. EPICv2 leukocyte via conversion — committed the slim
    `Anno/EPICv2/EPICv2ToEPIC_map.tsv.gz` (raw 115 MB gitignored); `qc.R`
    adds `find_anno_file()`, `load_epicv2_reliable_epic_probes()`, and
    `convert_epicv2_to_epic()`; `compute_sample_inferences()` converts EPICv2
    betas → EPIC and fills `Leukocyte_Fraction` (no longer NA). Map loaded once
    per batch. Verified on an EPICv2 example: leukocyte 0.0427.
  - **Verified end-to-end** (`scripts/run_example.R`, EPIC example, 4 samples):
    pipeline completes; `conversion_qc.csv` written with GCT 1.10–1.12;
    `sample_qc_report.csv` carries `GCT_Score`, `Sesame_Sex` (FEMALE/FEMALE/
    MALE/MALE), `Horvath_Age` (3.5–26.3y, `predictAge` matching the prior manual
    84.1y on TCGA PAAD), and `Leukocyte_Fraction` (~0.016). EPICv2 example:
    leukocyte 0.0427 via conversion; GCT correctly `NA` + note, no error. GCT
    gating logic unit-checked (1.48 fails, 1.11 passes, NA passes, existing
    detection-p failure unaffected). All edited R files parse.
  - **P1-GATE — PARTIALLY DONE.** Pipeline run + CSV outputs verified above.
    Still open: (a) open a completed run in the live Shiny app and confirm the
    "Conversion QC" tab + new columns/tooltips render; (b) trigger a GCT failure
    (e.g. `max_gct_score = 1.11`) and confirm the red fail-row + reason + banner;
    (c) regression-confirm a pre-existing run's `Pass_QC` is unchanged.

- **Phase 2 — EPICv2 GCT + karyotype + SNP fingerprint (IN PROGRESS).**
  - P2-T1. **DONE.** Resolved the open question: `sesameData` 1.30.0 has **no**
    `EPICv2.probeInfo`, so EPICv2 ext probes are derived from the Zhou Lab EPICv2
    manifest. Rule (validated against sesame's curated EPIC set): type-I probes
    with manifest `nextBase == "R"` are ext-C, `nextBase == "A"` are ext-T. The
    `nextBase`-derived set reproduces sesame's *native* GCT exactly (0.00% diff)
    on the EPIC example. Vendored the slim lists as
    `Anno/EPICv2/EPICv2.typeI.ext.rds` (extC=47125, extT=16185).
  - P2-T2. **DONE.** `compute_gct_scores()` loads the vendored ext probes via
    `load_epicv2_ext_probes()` and passes `extR`/`extA` to
    `bisConversionControl` for EPICv2; EPIC/450k keep the native auto-fetch path.
    EPICv2 examples now score 1.19–1.92 (was `NA`). `build_release.sh`
    sanity-checks the vendored file. (committed `7f071b5`)
  - P2-T3. EPICv2 example run produces finite GCT in `conversion_qc.csv` and a
    real `GCT_Score` in `sample_qc_report.csv` (no Phase-2 note).
  - P2-T4. **Sex karyotype** (user-requested; `inferSexKaryotypes` is gone from
    modern sesame). Reimplemented on the modern API in
    `compute_sample_inferences()`: `Karyotype` (XX/XY, plus hedged `XXY?`/`X0?`)
    from `inferSex` + X-inactivation heterozygosity `X_Het` (two-X ~0.45 vs
    one-X ~0.13, clean on EPIC/450k). Deliberately avoids Y-channel intensity
    (unreliable on degraded samples → spurious Turner flags). Informational,
    never gates `Pass_QC`. Verified: EPIC → XX/XX/XY/XY, 450k → XY.
  - P2-T5. **SNP genotype fingerprint** (replaces the removed `inferEthnicity`,
    whose model was dropped from sesameData). `SNP_Fingerprint` + `SNP_Count`
    from the ~59–65 Infinium `rs` probes (beta→0/1/2, NA→`.`, position-aligned by
    probe ID for cross-sample comparability). Verified: same-individual pairs
    ~5% discordance vs ~60–65% across individuals — clean sample-swap signal.
  - **P2-GATE.** EPICv2 run produces finite GCT scores; karyotype + fingerprint
    columns populate in `sample_qc_report.csv` across 450k/EPIC/EPICv2; the
    Phase-2 GCT note disappears for v2.

- **Phase 3 — Packaging & release (NOT STARTED).**
  - P3-T1. In-app Help tab + QUICKSTART: explain the GCT score and the
    near-1.0 interpretation.
  - P3-T2. Optionally surface `conversion_qc.csv` in the HTML report
    (`visualization.R`) alongside the QC table.
  - P3-T3. Bump `MEQTRACK_VERSION` to `2.3.0`; build and smoke-test the release
    zip; merge `preprocess-sesame-migration` → `main`.
  - **P3-GATE.** Fresh-install smoke test passes with the Conversion QC tab
    working. **DONE 2026-06-19** — v2.3.0 fast-forward merged to `main`
    (`df92f70`); `--step all` smoke-tested end-to-end on the EPIC example.
    (Note: the "Conversion QC" tab was instead removed — `GCT_Score` lives in the
    main Sample-metrics table — and a "Sample identity" SNP-concordance heatmap +
    per-sample "Dye bias" QQ tab were added.)

- **Phase 4 — Cell-type deconvolution (deconvMe) — PLANNED (v2.4.0, post-v2.3.0).**
  A **standalone, opt-in pipeline step** (NOT part of QC) — reference-based
  immune cell-type deconvolution via the omnideconv **deconvMe** package
  (https://github.com/omnideconv/deconvMe), wrapping 5 methods for Illumina
  array data: **EpiDISH, Houseman, MethAtlas, methylCC, methylResolver**.
  Architected exactly like `reference_projection`: its own `--step
  deconvolution`, own `deconvolution.R` module, own `deconv/` output dir, own
  top-level "Deconvolution" app tab + Run-tab per-step button, reloading
  `preprocessed_data.RData` so it runs independently. It is an analysis, not a
  QC gate — never touches `Pass_QC`, and (decision) is **excluded from
  `--step all`** so a standard run never needs the deconv deps; the user runs
  it explicitly.
  - **Why / how it relates to sesame leukocyte (VERIFIED on the 450k example):**
    sesame `Leukocyte_Fraction` answers *how much* immune content (one scalar,
    2-component purity); deconvMe answers *which* cell types and in what
    proportion (B/CD4T/CD8T/NK/Mono/Neutrophil, EpiDISH adds an `other` =
    non-blood residual). **They are NOT interchangeable, especially on tumours.**
    On the bundled TCGA tumour 450k samples they disagree hard: sesame leukocyte
    ≈ **0.05**, but EpiDISH summed-immune ≈ **0.82–0.85** and Houseman ≈
    **1.11–1.14** (Houseman/estimateCellCounts2 isn't sum-to-1). Reason: the
    blood references force-fit tumour methylation into blood cell types — they
    estimate the *composition of the immune compartment*, not *how much* immune.
    → **sesame `Leukocyte_Fraction` stays the purity/contamination number;
    deconvMe is the relative immune-subtype profile** (lymphocyte-rich → TIL,
    neutrophil-high → handling). The summed-immune ≈ sesame cross-check only
    holds for blood-like samples; a large gap on tumours is expected, not a bug.
  - P4-T1. **Input:** the new step reloads `preprocessed_data.RData`, takes
    `result$rgset`, and builds a minfi MethylSet via `preprocessRaw(rgset)` (same
    call preprocess already makes for getSex). deconvMe takes `methyl_set` +
    `array` ('450k'/'EPIC'). **EPICv2 is not natively supported** → reuse the
    existing EPICv2→EPIC conversion (the leukocyte path) or skip with a note.
  - P4-T2. **Run** `deconvMe::deconvolute_combined(methyl_set, array = '450k',
    methods = c('epidish','houseman'))` in `deconvolution.R`, dispatched by
    `--step deconvolution` (guarded `step == "deconvolution"` only — NOT in
    "all"). Methods configurable in Settings + run_manifest. Default to the
    pure-R methods (EpiDISH + Houseman); **MethAtlas is research-license-only**
    and methylCC/methatlas need Python — gate those behind opt-in.
  - P4-T3. **Output (VERIFIED):** `deconvolute_combined()` returns a tidy LONG
    table `(method, sample, celltype, value)` with an extra `aggregated` method
    row-set — write it straight to **`deconv/cell_fractions.csv`** (no reshaping
    needed). EpiDISH names: B cell / T cell CD4+ / T cell CD8+ / NK cell /
    Monocyte / Neutrophil / other.
  - P4-T4. **App:** a new **top-level "Deconvolution" tab** (sibling of QC /
    Dim. reduction / Reference projection / CNV), loaded by `results_loader`
    from `deconv/cell_fractions.csv`, plus a Run-tab per-step "Deconvolution"
    button. **Caveat (VERIFIED):** `deconvMe::results_barplot()` works on a
    single-method `deconvolute()` result (ggplot) but ERRORS on the
    `deconvolute_combined()` result ("Column name `sample` must not be
    duplicated"). So either call `deconvolute()` per method for the barplot, or
    build our own stacked ggplot/plotly from the long CSV (preferred — method
    selector/facets). Include a note comparing summed-immune vs
    `Leukocyte_Fraction` (and the tumour caveat above so the gap doesn't read as
    a failure).
  - P4-T5. **Provisioning:** `pak::pkg_install("omnideconv/deconvMe")` (GitHub,
    not CRAN/Bioc) + the large `FlowSorted.Blood.450k` / `FlowSorted.Blood.EPIC`
    reference-data packages (yamapData-style) wired into `setup.R` and the
    release bundle. `methylcc`/`methatlas` pull a **Python backend**
    (`deconvMe::init_python()`, basilisk/reticulate-style) — EpiDISH + Houseman
    are pure-R, so default to those to avoid the Python provisioning. Note:
    `deconvolute()` takes no `array` arg (only `deconvolute_combined()` does).
    Informational only — never gates `Pass_QC`.
  - **Open Q:** which method(s) to make default; whether to surface the
    aggregated result or per-method; EPICv2 support depth; whether to ever fold
    it into `--step all` (default: no, keep opt-in).
  - **P4-GATE.** `--step deconvolution` on a 450k run writes
    `deconv/cell_fractions.csv` with sensible fractions; the top-level
    Deconvolution barplot tab renders; `--step all` is unaffected (does not run
    or require deconv); the tumour-vs-blood-ref interpretation note is shown.

---

## Scope discipline

The GCT bisulfite-conversion score is the anchor capability of v2.3.0 (it gates
`Pass_QC`); sesame sex, Horvath age, and leukocyte fraction ride along as
informational sample-integrity columns. Phase 1 ships EPIC/450k GCT plus *all*
sample-integrity inferences across array types (EPICv2 leukocyte works via the
EPICv2→EPIC map). Phase 2 closes the last EPICv2 gap — *GCT* for EPICv2 — and (added at the user's
request) reintroduces **sex karyotype** and an **identity SNP fingerprint**, both
reimplemented on the modern sesame API rather than via the removed functions.

**Scope change (2026-06-19):** karyotype and "ethnicity" were originally
out-of-scope on the belief they needed a sesame upgrade. Investigation showed
`inferSexKaryotypes` and `inferEthnicity` were **removed** from modern sesame
(old `sset` API; ethnicity model dropped from sesameData) — an upgrade would
*not* restore them, and a downgrade would break the SigDF-based Phase 1. So
instead: **karyotype** is reimplemented from `inferSex` + X-inactivation
heterozygosity (P2-T4), and **ethnicity** is replaced by an `rs`-SNP genotype
fingerprint for sample-swap/identity matching (P2-T5) — the sample-integrity
signal actually wanted, without the deprecated low-reliability ethnicity model.

**Explicitly out of scope:**
- **Ethnicity classification proper** — the deprecated, removed sesame model is
  not reproduced; we ship the identity SNP fingerprint instead.
- **Aneuploidy as definitive calls** — karyotype asserts only XX/XY; XXY?/X0?
  are hedged guesses (no aneuploid samples available to validate thresholds).
- **Sex karyotype/ethnicity via sesame upgrade** — not viable (functions removed
  upstream; upgrade churns `renv.lock` and risks minfi/basilisk interop).
- Removing/repairing the inert `normalization` parameter and its stale
  user-facing strings — a separate cleanup, not this capability.
- Any other sesame QC metric (e.g. `sesameQC_calcStats` intensity/detection
  stats) — the main QC table already covers detection-p and intensity via minfi.

## Open questions still to resolve

- **Threshold:** the GCT fail cutoff is `config$qc$max_gct_score`, default
  **1.3** (a pragmatic middle ground — ~1.0 is ideal, and the bundled HM450
  example at 1.48 fails while the EPIC example at 1.11 passes). Open: validate
  this default against a real cohort / literature and adjust if warranted.
- **EPICv2 ext probes:** RESOLVED — `sesameData` 1.30.0 has no
  `EPICv2.probeInfo` at all; derived ext-C/ext-T from the Zhou Lab EPICv2
  manifest (`nextBase` R/A), validated to reproduce native GCT exactly on EPIC.
- **Karyotype thresholds:** `X_Het` cutoffs (two-X >0.32, one-X <0.25) are
  calibrated on normal XX/XY examples only. Open: validate the aneuploidy
  (XXY?/X0?) band against real Klinefelter/Turner samples if available.
- **Report surfacing (P3-T2):** include the conversion table in the shareable
  HTML report, or keep it app-only? Leaning toward including it for parity with
  the main QC table.
- **sesame upgrade for karyotype/ethnicity:** a newer sesame/sesameData would
  restore `inferSexKaryotypes` (XaY/XaXi), `inferEthnicity`, and a `predictAge`
  compatible with the bundled clock models. Worth it only if those signals are
  wanted — weigh against `renv.lock` churn and pipeline-wide compatibility
  (minfi interop, the basilisk/snifter stack). Not v2.3.0.
