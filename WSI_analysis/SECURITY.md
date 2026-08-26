# Security Policy

## Reporting a Vulnerability

If you discover a security vulnerability in this repository (for example, a
dependency with a known exploit, or a way the code could be misused to expose
data), please **do not** open a public GitHub issue.

Instead, report it privately using GitHub's [private vulnerability reporting](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing/privately-reporting-a-security-vulnerability)
feature on this repository, or contact the repository maintainers directly.

Please include:
- A description of the vulnerability and its potential impact.
- Steps to reproduce it, if possible.
- Any suggested remediation.

We aim to acknowledge reports within 5 business days and to provide a plan for
remediation as soon as practical, consistent with St. Jude's internal
vulnerability-response process.

## Scope

This repository contains research analysis code and de-identified, minimal-example
data (see [`data/README.md`](data/README.md)). It does not run any network service
and does not process live patient data — reports should focus on the code and its
dependencies (e.g., vulnerable third-party packages) rather than infrastructure
concerns.
