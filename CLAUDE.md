# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Indonesia Rheumatic Heart Disease (RHD) Simulation Model — a Markov State-Transition population-level simulation projecting the future burden of RHD in Indonesia and evaluating prevention interventions. Built on a parent NCD Countdown Population-level Chronic disease Model (CCPM) framework.

## Running the Pipeline

The pipeline is a numbered sequence of R scripts under `scripts/`. Run them in order, or source the master orchestrator:

```r
source("scripts/00_run_all.R")   # runs full pipeline (sets wd, global flags, sources all scripts)
```

Individual scripts can be sourced in R/RStudio using `here::here()` paths — they do not need to be run from within their own directory.

**Run order and purpose:**

| Script | Purpose | Run frequency |
|---|---|---|
| `01_prepare_inputs.R` | Loads GBD 2023 epidemiology/population CSVs from `data-raw/` | Each session |
| `02_build_demography.R` | Builds WPP2024 demographic backbone (`sf.wpp`); outputs `data/wpp/indonesia_ncd_demography.Rda` | Once per WPP version |
| `03_build_disease_model.R` | Loads and filters GBD cause-fraction data | Each run |
| `04_calibration_random_tp.R` | Random-search Monte Carlo calibration of Markov transition probabilities | When inputs change |
| `05_build_baseline.R` | Assembles baseline rates from calibrated TPs, WPP population, and COVID mortality | Each run |
| `06_run_prevention_model.R` | Runs reference vs. secondary-prevention scale-up scenarios | Each run |
| `07_make_outputs.R` | Compiles tables and figures from scenario results | After 06 completes |

## Key Global Flags (set in `00_run_all.R` before sourcing scripts)

```r
run_calibration_par  <- TRUE   # parallel calibration (recommended)
run_adjustment_model <- FALSE  # MUST be FALSE when using 04_calibration_random_tp.R
run_bgmx_trend       <- TRUE   # apply background mortality secular trends post-2019
run_CF_trend         <- TRUE   # apply case-fatality secular trends post-2019
run_CF_trend_80      <- TRUE   # use 80% of secular CF trend (default)
```

> **Critical**: `run_adjustment_model` must be `FALSE` when `04_calibration_random_tp.R` is used. Setting it `TRUE` re-applies the old 032 adjustment factors on top of already-calibrated TPs, double-adjusting IR and CF.

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
  epidemiology/*.csv   ← GBD 2023 cause-specific rates (not tracked by git)
  population/*.csv     ← GBD 2023 population (not tracked)
  GBD/                 ← Additional GBD files used by calibration
       ↓
data/
  wpp/indonesia_ncd_demography.Rda   ← sf.wpp from 02_build_demography.R
  tps_inpt_part*.rds                 ← baseline transition probabilities (input to calibration)
  adjusted_searo_part{1..10}.rds     ← calibrated TPs output by 04_calibration_random_tp.R
  model/scenarios/                   ← scenario results CSVs from 06
  model/baseline/                    ← baseline results from 05/06
       ↓
outputs/
  tables/   ← CSV tables from 07_make_outputs.R
  figures/  ← PNG figures from 07_make_outputs.R
```

### Core Model: Well-Sick-Dead Markov Recursion

Each `location × sex × cause` combo advances annual cohorts through three states:
- **Well** → via incidence rate (IR)
- **Sick** → via case-fatality rate (CF) or background mortality (BG.mx)
- **Dead**

Transition probabilities must satisfy: `IR + BG.mx ≤ 1` and `CF + BG.mx ≤ 1`.

### Calibration (`04_calibration_random_tp.R`)

Pure Monte Carlo random search (not hill-climbing). For each `location × sex × cause` combo:
- Draws `N_ITER = 400` i.i.d. uniform multipliers in `[1 ± SEARCH_HALFWIDTH]` for IR and CF
- Minimises weighted relative squared error against GBD 2009–2019 Deaths + Prevalence
- Parallelises across combos (not locations); Indonesia has one location so combo-level parallelism uses all available cores

Key tunable constants at the top of the file: `SEARCH_HALFWIDTH`, `GRANULARITY` ("combo" or "age_group"), `N_ITER`, `CONVERGE_TOL`, `SEED`.

### RHD Prevention Model (`06_run_prevention_model.R`)

Cohort state-transition model (2021–2030) with three disease stocks:
- **Mild** (asymptomatic RHD): SAP reduces progression to severe
- **Severe** (heart failure): HF management reduces mortality; surgery moves to post-surgery
- **Post-surgery**: residual RHD mortality

Two scenarios: `ref` (baseline) and `sap` (secondary-prevention scale-up). Outputs incremental costs, deaths averted, DALYs, and benefit-cost ratio.

### Demography (`02_build_demography.R`)

Builds `sf.wpp` from the `wpp2024` R package. `get.par()` back-solves net migration as the residual between WPP projected population and a forward CCPM projection. Output arrays are `[131 years × 2 sexes × 101 ages]` covering 2020–2150. `get.lt()` (life table function) and `locations` lookup are saved alongside `sf.wpp` in the `.Rda` file since the projection engine needs them at runtime.

### Shared Utilities

- `R/packages.R` — dependency checker; `source(here::here("R", "packages.R"))` at the top of any run script
- `R/utils.R` — `get.bp.prob()` (BP distribution shifts for HTN intervention), `calc_mortality_reduction()` (TFA policy), `create_age_groups()` (GBD age binning)

## File Path Conventions

- Scripts loaded via `00_run_all.R` use `wd`, `wd_raw`, `wd_data`, `wd_outp` globals (hardcoded absolute paths in `00_run_all.R`)
- Scripts meant to be sourced standalone use `here::here()` — requires the `.Rproj` file to be in the repo root
- `data-raw/` is not tracked by git; raw GBD and WPP CSVs must be obtained separately
