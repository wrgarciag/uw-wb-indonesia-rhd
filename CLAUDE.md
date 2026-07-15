# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Indonesia Rheumatic Heart Disease (RHD) Simulation Model — a Markov State-Transition population-level simulation projecting the future burden of RHD in Indonesia and evaluating prevention interventions. Built on a parent NCD Countdown Population-level Chronic disease Model (CCPM) framework.

## Running the Pipeline

The pipeline is a numbered sequence of R scripts under `scripts/`. Run them in order, or source the master orchestrator:

```r
source("scripts/00_run_all.R")   # single control panel: sets every parameter, sources 01→08
```

`00_run_all.R` is the SINGLE place to set every parameter (paths, year windows, ramp,
calibration settings, clinical/effect/coverage parameters, economic parameters, cores,
location list). Individual scripts can also be sourced standalone in R/RStudio — each reads
`00`'s globals when present and otherwise falls back to documented defaults via `here::here()`.

**Run order and purpose (each script reads the persisted outputs of the previous ones):**

| Script | Purpose | Writes |
|---|---|---|
| `01_prepare_inputs.R` | Load & filter GBD 2023 epidemiology CSVs | `data-raw/temp_baseline_rates_gbd.rds` |
| `02_build_demography.R` | WPP2024 population backbone; observed 1990–2024, projection 2025–2100 | `data/pop_observed_1990_2024.rds`, `data/pop_projection_2025_2100.rds`, `data/wpp/indonesia_rhd_demography.Rda` |
| `03_build_disease_model.R` | **Input builder only** (no run): data-fed rate arrays + clinical/effect/coverage params | `data/disease_model_inputs.rds` |
| `04_calibration_random_tp.R` | Random-search calibration of IR/CF, **full 0–95+**, window 2000–2019 | `data/adjusted_searo_part{1..10}.rds` + calibration CSVs |
| `05_build_baseline.R` | Assemble the initial state from 02+03+04 (no recompute) | `data/baseline_state.rds` |
| `06_run_prevention_model.R` | **Run** ref vs SAP scale-up; matrix WSD engine; parallel by location; no costs | `output/out_model/<location>.rds` (`$wsd` + `$tunnel`) |
| `07_make_outputs.R` | Compile the standard long table from `$wsd` | `output/tables/rhd_model_long.csv/.rds` |
| `08_economic_evaluation.R` | **All economics** (BCR, cost/death, cost/DALY) from `$tunnel`+`$wsd` | `output/tables/rhd_economic_summary.csv`, `rhd_budget_impact.csv` |

## Parameters (all set in `00_run_all.R`)

There are no longer any `run_*` trend/adjustment flags. `00_run_all.R` holds all parameters in
labelled blocks A–K: paths, locations, year/horizon windows, ramp window, incidence trend,
calibration settings (`run_calibration_par`, `SEARCH_HALFWIDTH`, `GRANULARITY`, `N_ITER`, …),
clinical params, intervention effect sizes + seed split, coverage targets, economic parameters
(block J — the ONLY monetary values in the pipeline), and parallel/core settings.

> **Provenance discipline**: scripts 01–07 contain NO monetary values (all economics live in
> `08`). Calibrated TPs from `04` are consumed AS-IS by 05/06 — no further adjustment is applied.

## Package Dependencies

Install CRAN packages:
```r
install.packages(c("here", "dplyr", "tidyr", "readr", "purrr", "tibble",
                   "ggplot2", "scales", "abind", "stringr", "readxl",
                   "cowplot", "data.table", "countrycode", "rlang",
                   "forecast", "RColorBrewer", "parallel", "doParallel",
                   "foreach", "gmodels"))
```

Install the WPP2024 demographic package (GitHub-only, required for `02_build_demography.R`):
```r
devtools::install_github("PPgp/wpp2024")
```

R >= 4.1.0 is required (uses the native pipe `|>`).

## Architecture

### Data Flow

```
data-raw/
  epidemiology/*.csv         ← GBD 2023 RHD + All-causes rates (not tracked by git)
  population/*.csv           ← GBD 2023 population (not tracked)
  temp_baseline_rates_gbd.rds← from 01
       ↓
data/
  pop_observed_1990_2024.rds       ← 02 (observed; 2024 = WPP jump-off)
  pop_projection_2025_2100.rds     ← 02 (medium-variant projection)
  wpp/indonesia_rhd_demography.Rda ← 02 (get.lt, locations, both tables)
  disease_model_inputs.rds         ← 03 (rate arrays + clinical/effect/coverage params)
  adjusted_searo_part{1..10}.rds   ← 04 (calibrated IR/CF, ages 0–95+, 2000–2019)
  baseline_state.rds               ← 05 (per-location initial state)
       ↓
output/
  out_model/<location>.rds  ← 06 (list: $wsd aggregate + $tunnel detail)  [gitignored]
  tables/                   ← 07 (rhd_model_long.*) + 08 (economics CSVs/rds)
```

### Core Model: Well-Sick-Dead Markov Recursion

Each `location × sex × cause` combo advances annual cohorts through three states:
- **Well** → via incidence rate (IR)
- **Sick** → via case-fatality rate (CF) or background mortality (BG.mx)
- **Dead**

Transition probabilities must satisfy: `IR + BG.mx ≤ 1` and `CF + BG.mx ≤ 1`.

### Calibration (`04_calibration_random_tp.R`)

It BUILDS an RHD-native baseline TP table (from GBD RHD rates + observed population) and then
calibrates by pure Monte Carlo random search (not hill-climbing). Purpose: improve the IR/CF
INPUTS only — it does NOT run the model or interventions. For each `location × sex × cause` combo:
- Draws `N_ITER` i.i.d. uniform multipliers in `[1 ± SEARCH_HALFWIDTH]` for IR and CF (candidate 0 = baseline, so calibrated is never worse than baseline)
- Full age range **0–95+** (paediatric groups `<1`, `12-23 months`, `2-4`, … included), window **2000–2019**
- Minimises weighted relative squared error against GBD RHD Prevalence + All-causes Deaths (background held fixed)
- Parallelises across combos; age-0 (births) inflow refreshes the youngest cohort each year

Key tunable constants (all controllable from `00`): `SEARCH_HALFWIDTH`, `GRANULARITY` ("combo" or "age_group"), `N_ITER`, `CONVERGE_TOL`, `SEED`.

### RHD Prevention Model (`06_run_prevention_model.R`)

Matrix-form (age × sex arrays) cohort state-transition model over the projection horizon
(2025–2100), whose "sick" state is resolved into three tunnel stocks:
- **Mild** (asymptomatic RHD): SAP reduces progression to severe
- **Severe** (heart failure): HF management reduces mortality; surgery moves to post-surgery
- **Post-surgery**: residual RHD mortality

Loads the initial state from `05` (no recompute), parallelises by location, and writes two tables
per location: `$wsd` (well-sick-dead aggregate, `sick = mild+severe+post`, ref+sap row-bound with a
`scenario` column) and `$tunnel` (the disaggregated sub-states + intervention volumes). NO costs
here — `07` emits the standard long table and `08` does all economics (BCR, deaths averted, DALYs).
`eff_ir`/`eff_cf` are reported RRR multipliers on [0,1] (1 = no reduction): `eff_ir = 1`
(secondary prevention doesn't cut incidence); `eff_cf = 1 − eff_sap_asymp × cov_sap`.

### Demography (`02_build_demography.R`)

Builds two tidy single-year population tables directly from the `wpp2024` R package: observed
1990–2024 (spline-smoothed, totals preserved; the 2024 slice is WPP's jump-off year, taken from
`popprojAge1dt`) and the 2025–2100 medium-variant projection. This closes the former 2024–2025 gap.
The former `get.par()` migration back-solve and 2020-anchored `sf.wpp` arrays are retired (the
forward engine they fed is absent from this repo). `get.lt()` (life-table fn) and the `locations`
lookup are saved alongside both tables in `data/wpp/indonesia_rhd_demography.Rda`.

### Shared Utilities

- `R/packages.R` — dependency checker; `source(here::here("R", "packages.R"))` at the top of any run script
- `R/utils.R` — `get.bp.prob()` (BP distribution shifts for HTN intervention), `calc_mortality_reduction()` (TFA policy), `create_age_groups()` (GBD age binning)

## File Path Conventions

- Scripts loaded via `00_run_all.R` use `wd`, `wd_raw`, `wd_data`, `wd_outp` globals (hardcoded absolute paths in `00_run_all.R`)
- Scripts meant to be sourced standalone use `here::here()` — requires the `.Rproj` file to be in the repo root
- `data-raw/` is not tracked by git; raw GBD and WPP CSVs must be obtained separately
