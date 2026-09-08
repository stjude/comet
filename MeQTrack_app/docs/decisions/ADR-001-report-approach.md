# ADR-001: Report approach for the MVP

Status: Accepted (Wave 1, 2026-04-22)

## Context

The Wave 1 plan asks for an explicit decision on how the app produces
its HTML report for a run. There are two viable options:

1. **Keep the existing `generate_report()` path** in
   `pipeline/pipeline_modules/visualization.R`. It builds an R Markdown document
   on the fly, knits to a self-contained HTML via `rmarkdown::render()` +
   `pandoc`, and already works against the pipeline's artifact tree.
2. **Replace it with a Shiny-rendered flat HTML.** The app would walk the
   pipeline outputs itself and produce an HTML from Shiny templates — no
   `pandoc`, no intermediate Rmd.

## Decision

Keep option 1 for the MVP: the app delegates report rendering to the existing
`generate_report()` pipeline module. The MVP does not fork or rewrite the
reporting logic.

## Why

- `generate_report()` is already integrated with the pipeline's output layout
  and file-naming conventions. Rewriting it would duplicate logic with no user
  benefit at MVP.
- `pandoc` is a documented host prerequisite (requirement S-INST-3);
  the user installs it once. In exchange, the app avoids carrying its own
  HTML-generation stack.
- The in-app views (t-SNE, UMAP, dendrogram, CNV) that Wave 4 adds are
  **separate** from the HTML report. The report is the shareable, portable
  artifact; the in-app views are for day-to-day inspection. Keeping them
  separate prevents report rigidity from being a blocker for UI iteration.
- If the Rmd approach proves too rigid when the UI evolves, the escape hatch
  is to introduce a Shiny-rendered HTML as a **side output** in Wave 4 rather
  than trying to evolve `generate_report()`. This is called out as a watch
  item under the plan's "Cross-cutting risks".

## Consequences

- Host must have `pandoc` >= 2.x (already required by `rmarkdown`).
- Any formatting change to the report that a user requests flows through
  edits to `pipeline_modules/visualization.R`, not through Shiny code. Document
  this in contributor notes when that time comes.
- The Wave 3 "Open report in browser" action simply opens the file produced
  by the pipeline; no UI-side rendering is involved.

## Revisit if

- A user requests interactive behaviour inside the report that the Rmd
  generator makes awkward (e.g. dynamic filtering across multiple samples).
- The Rmd generator accumulates conditional branches for parameter-specific
  sections to the point that maintaining it becomes a tax.
