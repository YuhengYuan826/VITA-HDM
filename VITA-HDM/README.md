# VITA-HDM paper code

This repository contains the R code for the simulation study, glioma case
study, and sensitivity analyses.

The code files are kept in one numbered sequence so that their roles and
execution order are visible directly in the repository file list.

## Run order

Files `01`–`04`, `08`, and `13` contain functions and are sourced
automatically. Do not run them separately.

### Simulation study

```bash
Rscript 05_run_simulation_study.R
Rscript 06_summarise_simulation_results.R
```

### Glioma case study and sensitivity analyses

```bash
Rscript 07_prepare_glioma_data.R
Rscript 09_run_glioma_discovery.R
Rscript 10_run_glioma_fixed_signature_cv.R
Rscript 11_run_glioma_conditional_effects.R
Rscript 12_summarise_glioma_results.R
Rscript 14_run_glioma_selection_sensitivity.R
Rscript 15_run_glioma_prediction_sensitivity.R
Rscript 16_summarise_glioma_sensitivity.R
```

Files `10` and `11` both depend on file `09` and can run in either order.
File `14` expects their completed outputs for the SA8 audit. Files `14` and
`15` can run in either order after the main case-study analyses. File `16`
must run after both.

After all workflows finish, check the expected output set:

```bash
Rscript 91_verify_outputs.R
```

## Numbered file guide

| File | Run directly? | Purpose |
| --- | --- | --- |
| `00_project_setup.R` | No | Shared project-relative paths and dependency inventory |
| `01_simulation_utils.R` | No | Simulation covariance, scaling, censoring, BFDR, and C-index utilities |
| `02_simulation_data.R` | No | Simulation truth, source perturbation, platform shift, and AFT data generation |
| `03_simulation_aft_cavi.R` | No | Simulation CAVI engine |
| `04_simulation_vita_hdm.R` | No | Simulation design preparation and HDM/VITA-HDM wrapper |
| `05_run_simulation_study.R` | Yes | Main simulation and reduced structure-only sensitivity |
| `06_summarise_simulation_results.R` | Yes | Simulation summaries and figures |
| `07_prepare_glioma_data.R` | Yes | CGGA/TCGA loading, filtering, alignment, and local prepared object |
| `08_glioma_vita_functions.R` | No | Case-study design, model, transfer, selection, and prediction functions |
| `09_run_glioma_discovery.R` | Yes | Main glioma discovery and target-anchored reported set |
| `10_run_glioma_fixed_signature_cv.R` | Yes | Descriptive repeated CV for fixed signatures |
| `11_run_glioma_conditional_effects.R` | Yes | Conditional bootstrap effects for fixed mediator sets |
| `12_summarise_glioma_results.R` | Yes | Main case-study tables and text summary |
| `13_glioma_sensitivity_common.R` | No | Shared sensitivity inputs, baseline classifications, and helper functions |
| `14_run_glioma_selection_sensitivity.R` | Yes | SA1, SA2, SA3, SA7, and SA8 |
| `15_run_glioma_prediction_sensitivity.R` | Yes | SA4, SA5, and SA6 |
| `16_summarise_glioma_sensitivity.R` | Yes | Final sensitivity tables and figures |
| `90_check_environment.R` | Optional | Records R, package, and graphics-device information |
| `91_verify_outputs.R` | Yes, after completion | Checks the expected output set and SA6 privacy boundary |

## R environment

Do not upgrade all packages to their latest versions solely for publication.
Use the environment that produced the verified paper results and record it:

```bash
Rscript 90_check_environment.R
```

This writes the installed package versions, `sessionInfo()`, and graphics
capabilities under:

```text
results/run_metadata/
```

Required CRAN or recommended packages are:

```text
MASS
survival
dplyr
tidyr
ggplot2
cowplot
doParallel
foreach
data.table
stringr
scales
patchwork
```

Required Bioconductor packages are:

```text
AnnotationDbi
org.Hs.eg.db
```

`ggrepel` is optional. If it is absent, SA6 uses deterministic plain-text
labels. `patchwork` is required for the complete main sensitivity figure.

The exact working versions were not encoded in the submitted scripts. If
`renv` is adopted, create `renv.lock` from the already verified working
environment, not from a newly upgraded environment.

## Data placement

The repository does not contain participant-level CGGA or TCGA data. Place the
downloaded files as documented in `data/README.md`:

```text
data/
├── CGGA/
└── TCGA/
```

The prepared aligned object is generated locally under:

```text
results/case_study/data/
```

Raw data, expression matrices, mutation records, aligned RDS files, and fitted
objects must not be committed.

## Paths

No personal absolute paths are required. Every numbered runner locates
`00_project_setup.R` beside itself.

Default locations are:

```text
data/       input data
results/    generated output
```

External locations can be configured without editing the scripts:

```bash
export VITA_HDM_DATA_DIR=/path/to/local/data
export VITA_HDM_RESULTS_DIR=/path/to/local/results
```

`VITA_HDM_DATA_DIR` must contain the `CGGA/` and `TCGA/` subdirectories.
`VITA_HDM_RESULTS_DIR` is the parent directory under which `simulation/` and
`case_study/` are created.

## Output structure

```text
results/
├── simulation/
├── run_metadata/
└── case_study/
    ├── data/
    ├── discovery/
    ├── fixed_signature_cv/
    ├── conditional_effects/
    ├── summary/
    └── sensitivity/
        ├── tables/
        ├── figures/
        ├── private_participant_level/
        └── final_outputs/
```

Existing scientific output filenames, CSV schemas, and overwrite behavior are
retained.

SA6 PCA coordinate files contain `Sample_ID` and are written only to
`private_participant_level/`. Do not commit that directory. SA6 figures and
aggregate correlation, quadrant, and source-weight summaries remain available.

## Sensitivity-analysis scope

The sensitivity code is included because it supports reported robustness and
diagnostic claims:

- file `14`: noise knock-in, source dose, BFDR/stability grid, hyperparameter
  grid, and target-anchoring audit;
- file `15`: fixed-signature prediction, Cox-Snell diagnostics, and
  source-target transferability;
- file `16`: main and supplementary sensitivity outputs.

The primary SA1 setting adds 50 synthetic noise features. The optional
25/50/100 grid remains disabled, matching the submitted workflow.

## Random seeds

The submitted seed rules are retained. Important examples include:

- simulation truth and target data: `10000 + repetition`,
  `20000 + repetition`, and `30000 + repetition`;
- simulation source parameters and data:
  `40000 + 10 * repetition + source` and
  `50000 + 10 * repetition + source`;
- discovery bootstrap: `4100 + bootstrap`;
- fixed-signature CV: master seed 2030, with repeat-specific fold seeds;
- conditional effects: `5000 + bootstrap`;
- sensitivity-specific seed formulas remain next to their settings in files
  `14` and `15`.

No script changes `RNGkind()`. Retaining the verified R version is therefore
important for stochastic reproducibility.

## Scientific preservation and validation

The numbered reorganization does not intentionally change:

- model equations or algorithms;
- scientific parameters or thresholds;
- random seeds or seed formulas;
- preprocessing, cohort construction, or standardization;
- bootstrap, cross-validation, or aggregation rules;
- output values, schemas, or figure content.

The code was statically checked during repository preparation but was not run
against CGGA/TCGA data in the preparation environment. Before public release,
run it in a fresh folder using the verified R environment and compare all
generated scientific outputs with the original run.

Add the manuscript citation, author information, repository DOI, and selected
software license before public release.
