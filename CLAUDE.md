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
| `02_build_demography.R` | WPP2024 population backbone; observed `OBS_YEARS`, projection `PROJ_YEARS` (stable filenames, no year in the name) | `data/pop_observed.rds`, `data/pop_projection.rds`, `data/wpp/indonesia_rhd_demography.Rda` |
| `03_build_disease_model.R` | **Input builder only** (no run): data-fed rate arrays + A/B/C/D natural history + effects + surgery service + cascade params; **clips the model horizon to `ANALYSIS_YEARS`** | `data/disease_model_inputs.rds` (`meta$years` = analysis horizon) |
| `04_calibration_random_tp.R` | Random-search calibration of IR/CF, **full 0–95+**, window 2000–2019 | `data/adjusted_searo_part{1..10}.rds` + calibration CSVs |
| `05_build_baseline.R` | Assemble the initial state from 02+03+04 (no recompute) | `data/baseline_state.rds` |
| `06_run_prevention_model.R` | **Run** ref vs SAP scale-up; matrix A/B/C/D engine + surgery service; parallel by location; no costs | `output/out_model/<location>.rds` (`$wsd` + `$stages`) |
| `07_make_outputs.R` | Compile aggregate + stage + flow long tables | `output/tables/rhd_model_long.csv/.rds`, `rhd_stage_model_long.csv`, `rhd_annual_flows.csv` |
| `08_economic_evaluation.R` | **All economics** (BCR, cost/death, cost/DALY) from `$stages`+`$wsd` | `output/tables/rhd_economic_summary.csv`, `rhd_budget_impact.csv` |

## Parameters (all set in `00_run_all.R`)

There are no longer any `run_*` trend/adjustment flags. `00_run_all.R` holds all parameters in
labelled blocks A–L. Block C defines the **configurable analysis period** — set `ANALYSIS_YEAR_START`
/ `ANALYSIS_YEAR_END` (→ `ANALYSIS_YEARS`) and it propagates through demography, disease inputs,
baseline, model run, outputs, economics and the report. Block C keeps SEVEN periods conceptually
separate (analysis; observed-data `OBS_YEARS`; projection `PROJ_YEARS`; calibration `CAL_YEAR_*`;
GBD rate reference `RATE_BASE_YEAR`; GDP-pc reference `gdp_pc_base_year`; economic discount
`discount_base_year`, default = `ANALYSIS_YEAR_START`), plus the independent intervention ramp
(block D). Other blocks: paths, locations, incidence trend, calibration settings
(`run_calibration_par`, `SEARCH_HALFWIDTH`, `GRANULARITY`, `N_ITER`, …),
**A/B/C/D natural history** (transition + per-stage RHD-death probabilities, `rhd_d_fraction`),
**intervention effects** (`sap_rrr_rhd_death = 0.55`; surgery RRRs on C→D and D→death),
**care cascade + surgery service** (screening/diagnosis/optimal-treatment baseline+2050 targets;
surgery requirement fractions and coverage per arm), economic parameters (block J — the ONLY
monetary values in the pipeline), and parallel/core settings.

> **Provenance discipline**: scripts 01–07 contain NO monetary values (all economics live in
> `08`). Calibrated IR/CF from `04` are consumed AS-IS by 05/06 — no further adjustment is applied.
> A/B/C/D stage parameters come from `00`/`03` (stage calibration is pending local echo targets).

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
  pop_observed.rds                 ← 02 (observed OBS_YEARS; last year = WPP jump-off)
  pop_projection.rds               ← 02 (medium-variant projection PROJ_YEARS)
  wpp/indonesia_rhd_demography.Rda ← 02 (get.lt, locations, both tables)
  disease_model_inputs.rds         ← 03 (rate arrays + A/B/C/D + effects + surgery + cascade)
  calibrated_rhd_parameters.rds    ← 04 (single bundle: calibrated IR/CF + stage-calib interface)
  calibration_targets_stage_template.csv ← 04 (Layer-2 stage-target schema)
  baseline_state.rds               ← 05 (per-location initial state, seed A/B/C/D)
       ↓
output/
  out_model/<location>.rds  ← 06 (list: $wsd aggregate + $stages A/B/C/D + surgery trace)  [gitignored]
  tables/                   ← 07 (rhd_model_long.* + rhd_stage_model_long.* + rhd_annual_flows.*) + 08 (economics)
```

### Aggregate Well-Sick-Dead view (`$wsd`)

The `$wsd` output preserves the aggregate well-sick-dead accounting, where the "sick" compartment
is the sum of the four RHD stages:
- **Well** (No RHD) → incident RHD enters via calibrated incidence rate (IR)
- **Sick** = `A + B + C + D` (the four living RHD stages; see below)
- **Dead** (RHD-specific + competing other-cause)

### Disease structure: WHF RHD stages A/B/C/D

The living "sick" state is resolved into four stages (No RHD → A ↔ B ↔ C → D → RHD death, with
competing other-cause death from every stage; incident RHD enters **A**; no separate ARF state):
- **A** minimal/early RHD  · **B** mild established  · **C** advanced without complications
- **D** advanced WITH complications (heart failure, requiring surgery)

Per-cycle transitions (with optional adjacent regression) and per-stage RHD-death probabilities are
scalar `[CALIBRATE]` parameters broadcast over the age×sex matrices. **Surgery is a clinical
SERVICE, not a state** — see below.

### Calibration (`04_calibration_random_tp.R`) — two layers

**Layer 1 (runs):** BUILDS an RHD-native baseline IR/CF table (from GBD RHD rates + observed
population) and calibrates by pure Monte Carlo random search — improving the IR/CF INPUTS only
(does not run the model). For each `location × sex × cause` combo: draws `N_ITER` i.i.d. uniform
multipliers in `[1 ± SEARCH_HALFWIDTH]` for IR and CF (candidate 0 = baseline); full age range
**0–95+**, window **2000–2019**; minimises weighted relative squared error against GBD RHD
Prevalence + All-causes Deaths (background held fixed). Calibrated IR drives the inflow into stage
A; calibrated CF is a total-RHD-mortality anchor.

**Layer 2 (interface):** A/B/C/D stage calibration against LOCAL echocardiographic stage-prevalence
targets — target-table schema + weighted squared-log-error loss are defined, but run only if a
local target file is supplied. None exists yet, so it falls back to Layer 1, preserves the
uncalibrated stage params from `00`, writes `calibration_targets_stage_template.csv`, and prints a
prominent "stage calibration pending local echo targets" message. Output is a single self-describing
bundle `data/calibrated_rhd_parameters.rds` (`$tp` / `$factors` / `$diagnostics` /
`$stage_calibration` / `$meta`) — replacing the former ten `adjusted_searo_part*.rds` chunks.

Key tunable constants (all controllable from `00`): `SEARCH_HALFWIDTH`, `GRANULARITY` ("combo" or "age_group"), `N_ITER`, `CONVERGE_TOL`, `SEED`.

### RHD Prevention Model (`06_run_prevention_model.R`)

Matrix-form (age × sex arrays) cohort state-transition model over the analysis horizon
(`ANALYSIS_YEARS`, configurable in `00_run_all.R`; e.g. 2026–2050) with four living stocks
**A, B, C, D** advanced by elementwise matrix arithmetic
(`zero_mat`/`age_shift`/`melt_year`, terminal 100+ open group). Each cycle:
- new incident RHD `= no_rhd × IR` enters **A** (`eff_ir = 1`: no incidence effect);
- **SAP** cuts every stage's RHD-death probability by `1 − sap_rrr_rhd_death × effective_treatment_coverage`;
- **Surgery** (a SERVICE, not a state): a fraction of the C/D stock is operated each cycle; its only
  epidemiological effect is a risk reduction on `C→D` and `D→RHD-death` (population-average reach =
  requirement fraction × surgery coverage). Surgery never removes people from C or D;
- competing other-cause background mortality is data-fed and the same age×sex risk for every stage.

Only the care cascade differs between the reference and SAP arms (surgery held equal in both).
Writes two tables per location: `$wsd` (aggregate, `sick = A+B+C+D`, ref+sap row-bound with a
`scenario` column) and `$stages` (four stocks + all transition flows + per-stage deaths + cascade
volumes + surgery trace). NO costs here — `07` emits the long tables and `08` does all economics.
`eff_ir`/`eff_cf` are reported RRR multipliers on [0,1] (1 = no reduction): `eff_ir = 1`;
`eff_cf = 1 − sap_rrr_rhd_death × effective_treatment_coverage` (the SAP RHD-**mortality** multiplier).

### Demography (`02_build_demography.R`)

Builds two tidy single-year population tables directly from the `wpp2024` R package: observed
`OBS_YEARS` (default 1990–2024; spline-smoothed, totals preserved; the last slice is WPP's jump-off
year, taken from `popprojAge1dt`) and the `PROJ_YEARS` medium-variant projection (default starts the
year after observed ends and runs to `ANALYSIS_YEAR_END` — no fixed 2100 endpoint). Written to the
stable filenames `pop_observed.rds` / `pop_projection.rds`. This closes the former 2024–2025 gap.
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
