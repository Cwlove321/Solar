# solar-code

Long-term maintained tools, workflows, and case-study code for solar physics research.

This repository is organized around three layers:

- `src/solarutils`: reusable Python utilities that should be stable across papers and case studies.
- `workflows`: repeatable analysis pipelines assembled from `solarutils`.
- `cases`: event-specific or project-specific notebooks, scripts, figures, and data manifests.

Large raw data and transient products should not be committed directly. Record them in a case-level
`data_manifest.md` instead.

