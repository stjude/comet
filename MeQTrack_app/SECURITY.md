# Security Policy

## Reporting a Vulnerability

If you discover a security vulnerability in this repository (for example, a dependency with a
known exploit, or a way the code could be misused), please **do not** open a public GitHub
issue.

Instead, report it privately using GitHub's [private vulnerability reporting](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing/privately-reporting-a-security-vulnerability)
feature on this repository, or contact the repository maintainers directly.

Please include:
- A description of the vulnerability and its potential impact.
- Steps to reproduce it, if possible.
- Any suggested remediation.

We aim to acknowledge reports within 5 business days and to provide a plan for remediation as
soon as practical, consistent with St. Jude's internal vulnerability-response process.

## Reporting exposed data

This repository processes patient-derived methylation data, so **data exposure is treated as a
security issue**. If you believe any file here contains identifiable information — a clinical
export, a raw IDAT from a non-public sample, an internal identifier, or a sample-level
metadata column that should not have shipped — report it privately using the same channel
above rather than opening an issue.

Only public, de-identified data is intended to ship here; see the
[Data and privacy](README.md#data-and-privacy) section of the README for what that means in
practice and how the reference label file is derived.

## Scope

MeQTrack runs entirely on the user's own machine. The Shiny application binds only to
`127.0.0.1` and is not network-accessible by construction; there is no server component, no
account system, and no telemetry. Reports should focus on the code and its dependencies (for
example, a vulnerable CRAN or Bioconductor package pinned in `renv.lock`) rather than
infrastructure concerns.
