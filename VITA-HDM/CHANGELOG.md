# Repository preparation changes

## Numbered layout

The previous multi-directory code layout was replaced with one flat numbered
sequence:

```text
00                 shared non-scientific setup
01-06              simulation study
07-12              glioma case study
13-16              sensitivity analyses
90-91              environment and output checks
```

This follows the readability principle used by the TransHDM reference
repository while retaining project-relative paths and research-data safeguards.

## Organizational changes

- Consolidated the path and dependency helpers into `00_project_setup.R`.
- Removed the separate `R/`, `scripts/`, and `docs/` code hierarchy.
- Replaced internal source paths with explicit numbered filenames.
- Simplified the sensitivity summarizer so it loads the single numbered common
  file directly instead of searching for duplicate copies.
- Kept generated outputs separated by scientific workflow under `results/`.
- Retained environment-variable overrides for external data and results.

## Publication safeguards retained

- Raw CGGA/TCGA data and all RDS files remain ignored.
- SA6 coordinate tables containing `Sample_ID` remain under
  `private_participant_level/`.
- Generated results remain ignored unless deliberately reviewed and added.
- No package versions were invented.

## Scientific behavior retained

No intentional changes were made to model equations, parameters, random seeds,
data transformations, fitting order, result schemas, or figure-generation
logic. Numerical equivalence still requires a fresh-folder run in the verified
R environment.
