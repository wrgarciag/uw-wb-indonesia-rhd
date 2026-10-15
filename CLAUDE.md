# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Rheumatic Heart Disease (RHD) secondary-prevention investment case — a population-level Markov
state-transition model (WHF stages A/B/C/D) projecting RHD burden and evaluating echo screening +
diagnosis + secondary antibiotic prophylaxis (SAP) scale-up, with benefit-cost / cost-effectiveness.
Originally built for **Indonesia**; now **multi-country** (Indonesia, Uganda) via a single switch.
Descends from the NCD Countdown CCPM framework (legacy helpers remain in `R/utils.R`, unused).

## Running the Pipeline

```r
source("scripts/00_run_all.R")   # single control panel: sets every parameter, sources 01 → 08
```

`00_run_all.R` is the SINGLE place to set parameters. Scripts can be sourced standalone (each reads
`00`'s globals when present, else falls back to defaults via a `getp()` / `exists()` guard and
`here::here()`). **Gotcha:** standalone fallbacks point at the *root* `data/` and `output/`, not the
country-scoped folders — set `wd_data`/`wd_outp` to `data/<COUNTRY>/`, `output/<COUNTRY>/` first.

| Script | Purpose | Writes (under `data/<COUNTRY>/` or `output/<COUNTRY>/`) |
|---|---|---|
| `01_prepare_inputs.R` | Read every CSV in `data-raw/epidemiology/` + `population/`, keep RHD + All causes for `LOCATION` (retains GBD `upper`/`lower` UIs) | `temp_baseline_rates_gbd.rds`, `temp_population_gbd.rds` |
| `02_build_demography.R` | WPP2024 population backbone: observed `OBS_YEARS`, medium-variant `PROJ_YEARS` | `pop_observed.rds/.csv`, `pop_projection.rds/.csv`, `wpp/rhd_demography.Rda` |
| `03_build_disease_model.R` | **Input builder only**: rate arrays + A/B/C/D priors + effects + surgery + cascade; clips horizon to `ANALYSIS_YEARS` | `disease_model_inputs.rds` |
| `04_calibration_random_tp.R` | **Layer 1** aggregate calibration of IR + RHD CF (2000–2019) | `calibrated_rhd_parameters.rds` + `calibration_layer1_*.csv` |
| `04b_calibrate_structural.R` | **Layer 2** structural A/B/C/D calibration on the real engine | fills `$stage_calibration` in the bundle + `calibration_structural_*.csv` |
| `05_build_baseline.R` | Assemble initial state (uses 04b `best_parameters`); realise cascade trajectories | `baseline_state.rds` |
| `06_run_prevention_model.R` | **Run** ref vs SAP arms (parallel by location); no costs | `out_model/<location>.rds` (`$wsd`, `$stages`, `$meta`) |
| `07_make_outputs.R` | Long tables | `tables/rhd_model_long.*`, `rhd_stage_model_long.*`, `rhd_annual_flows.*` |
| `08_economic_evaluation.R` | **All economics** (BCR, cost/death, cost/DALY, budget impact) | `tables/rhd_economic_summary.csv`, `rhd_budget_impact.csv`, `rhd_economic_results.rds` |

### Tests and report

```r
Rscript tests/test_abcd_smoke.R          # deterministic engine test (mass balance, hazards, SAP, surgery)
Rscript tests/verify_abcd_pipeline.R     # 23 spec checks on persisted 03–08 outputs (run after pipeline)
rmarkdown::render("analysis/report.RMD", params = list(country = "Uganda"))   # reads outputs only
```

Both test scripts exit non-zero on failure. `verify_abcd_pipeline.R` defaults to root `data/`/`output/`
— point `wd_data`/`wd_outp` at the country folders. The report resolves the repo root from its own
location and never re-runs the model.

### Windows / CLI notes

- R is not on PATH; use R 4.5.1 (`C:/Program Files/R/R-4.5.1/bin/Rscript.exe`) — `wpp2024` segfaults under 4.3.
- Load `data.table` **before** `wpp2024` (OpenMP segfault otherwise); `00` and `02` already do.
- Pandoc is not available to CLI Rscript: `rmarkdown::render` works only from RStudio. To smoke-test the
  report from the CLI use `knitr::knit()` with `knitr::opts_knit$set(rmarkdown.pandoc.to = "html")`.

## Parameters (`00_run_all.R`, blocks A–L)

- **A** `COUNTRY` (`"Indonesia"` | `"Uganda"`) → `ISO3`, paths. `data/<COUNTRY>/` and `output/<COUNTRY>/`
  are country-scoped; `data-raw/` is shared. Most later parameters are `switch(COUNTRY, ...)`; adding a
  country means adding an arm to every switch.
- **B** `LOCATION`/`LOCATIONS`, per-country sanity bands (`POP2020_*`, `RHD_PREV_*`, `RHD_DEATH_*`,
  `ALLC_DEATH_*`) — order-of-magnitude gates in 02/03/05/06; don't loosen globally.
- **C** Year windows, kept conceptually separate: analysis `ANALYSIS_YEAR_START/END` (default
  **2026–2100**; propagates to everything downstream), observed `OBS_YEARS` (1990–2024), projection
  `PROJ_YEARS` (contiguous with observed, ends at analysis end), calibration `CAL_YEAR_*` (2000–2019),
  `RATE_BASE_YEAR` (2023), `gdp_pc_base_year`, `discount_base_year` (= analysis start). `stopifnot`
  checks enforce consistency.
- **D** Intervention ramp `ramp_start`/`ramp_end` (2026→2050), independent of the analysis period.
- **E** `incidence_trend` (annual multiplier).
- **F** Layer-1 calibration: weights `W_INC/W_PREV/W_DEATH`, `USE_IV_WEIGHTS`, `INCIDENCE_CALIBRATION_MODE`
  (`fixed`/`anchored`/`free`), `CALIB_AGE_BANDS`, priors `SIGMA_ALPHA/BETA`, `LAMBDA_*`, bounds
  `MULT_LO/HI`, `CALIB_OPTIMIZER` (`multistart` default; `nloptr`/`deoptim` optional), `N_STARTS`, `SEED`.
- **G** A/B/C/D transition and per-stage RHD-death **prior centers** + `rhd_d_fraction` (04b calibrates them).
- **G2** Layer-2 settings: `run_structural_calibration`, `STRUCT_MATCH_FROM`, `STRUCT_BURN_YEARS`,
  `SIGMA_MORT` (tight) / `SIGMA_TRANS` (loose), `RESID_MULT_LO/HI`, `MD_LO/HI`, `STRUCT_N_STARTS`.
- **H** Effects: `sap_rrr_rhd_death = 0.55`; `eff_surgery_C_to_D`, `eff_surgery_D_to_rhd_death`.
- **I** Care cascade (screen/diagnosis/treatment base + target) and surgery service (requirement
  fractions, coverage per arm); optional `screen_age_restrict`.
- **J** Economics — the ONLY monetary values in the pipeline (consumed by 08 only).
- **K** `run_model_par`, `MAX_CORES`. **L** sources 01 → 08.

Commented-out blocks in `00` are the superseded hand-tuned values, kept for reference.

> **Provenance discipline**: scripts 01–07 contain NO monetary values. Calibrated IR/CF (04) and stage
> parameters (04b) are baked in ONCE and consumed as-is by 05/06 — never re-adjust downstream. The former
> hand-tuned mortality scalars (1/2.5 Indonesia, 1/8 Uganda) were removed; do not reintroduce them.

## Architecture

### Directory layout

```
scripts/        pipeline 00–08 (+ 04b)
R/              abcd_engine.R (shared one-cycle engine), packages.R (dependency check), utils.R (legacy NCD helpers, unused)
tests/          test_abcd_smoke.R, verify_abcd_pipeline.R
analysis/       report.RMD (parameterised by country)
docs/           metadata_inputs.csv (input-data catalogue), WHF guidelines/papers, state-diagram PNGs + generator .py
dev/            refactor prompts / design notes
data-raw/       GBD 2023 epidemiology + population CSVs (not tracked; shared across countries)
data/<COUNTRY>/ processed inputs, calibration bundle + CSVs, baseline state
output/<COUNTRY>/ out_model/ (gitignored) + tables/
```

Files directly under `data/`, `output/tables/`, and `data-raw/*.rds`, plus year-suffixed
`pop_*_YYYY_YYYY.*` files, are stale pre-multi-country artifacts — the pipeline does not read them.

### Disease structure and engine

No RHD → A ↔ B ↔ C → D → RHD death, with competing other-cause death from every stage; incident RHD
enters **A** (no ARF state). A minimal/early · B mild · C advanced w/o complications · D advanced with
HF. **`R/abcd_engine.R`** (`abcd_one_cycle`, `abcd_compete`) is the single source of truth for one
annual cycle, used by both 06 and 04b. Competing risks combine on the **hazard scale**
(`p_any = 1 − exp(−Σh)`, split by hazard share) so no clipping is needed and stocks balance exactly.

- **SAP** multiplies every stage's RHD-death probability by `1 − sap_rrr_rhd_death × effective_treatment_coverage` (`eff_cf`); no incidence effect (`eff_ir = 1`).
- **Surgery is a SERVICE, not a state**: reach = requirement fraction × coverage; it reduces `C→D` and `D→death` only and never removes stock. Held equal in both arms.
- Background mortality = GBD all-cause − RHD, same age×sex risk for every stage. Ages 0–100 (100+ open).
- Only the care cascade differs between `ref` and `sap`.

06 writes `$wsd` (aggregate well-sick-dead, `sick = A+B+C+D`) and `$stages` (stocks, all flows,
per-stage deaths, cascade volumes, surgery trace); collapsing `$stages` reproduces `$wsd` exactly.

### Calibration — two layers

**Layer 1 (04):** per location × sex, one log-scale incidence multiplier (α) and one aggregate
RHD-mortality multiplier (β) per broad age band, for a well-sick-dead proxy with hazard competing risks.
Loss = log-scale squared error vs GBD RHD **incidence, prevalence and RHD-specific deaths** (optional
inverse-variance weights from GBD UIs) + Gaussian priors + 2nd-difference smoothness; bounded
multi-start L-BFGS-B, sequenced mortality → incidence → joint. All-cause deaths are validation only.
Bundle: `calibrated_rhd_parameters.rds` with `$tp`, `$factors`, `$diagnostics`, `$layer1`,
`$stage_calibration`, `$meta`.

**Layer 2 (04b):** runs the real engine (equilibrated by a burn-in) to fit transitions, severity-ordered
stage mortality and `rhd_d_fraction` to GBD RHD prevalence + RHD deaths. Mortality has tight clinical
priors, transitions loose priors, so excess deaths are resolved by moving the stock (notably
`p_A_to_no_rhd`), not by suppressing mortality. A residual mortality multiplier at its bound is a
validation FAILURE. Stage split is partially identified: near-optimal sets and an identifiability table
are saved. A local echo stage-prevalence target file, if supplied, adds a stage loss
(`calibration_targets_stage_template.csv` shows the schema).

### Demography (02)

Built directly from `wpp2024`: observed table spline-smoothed with totals preserved (2024 is taken from
`popprojAge1dt`, WPP's jump-off year), projection taken as-is. Stable filenames (no year in the name).

## Dependencies

R ≥ 4.1 (native pipe). `R/packages.R` checks: here, dplyr, tidyr, readr, purrr, tibble, ggplot2,
scales, abind, stringr, readxl, cowplot, data.table, countrycode, rlang. Also used: foreach,
doParallel, parallel (06); optional nloptr / DEoptim (04 alternative optimisers); knitr, rmarkdown,
kableExtra, DT (report). `wpp2024` is GitHub-only: `devtools::install_github("PPgp/wpp2024")`.

## Conventions

- `data.table` throughout pipeline scripts; parameters tagged `[PAPER]`, `[LIT]`, `[CALIBRATE]`,
  `[TODO: confirm for Uganda]`.
- Scripts communicate only via persisted files — never rely on in-memory objects from an earlier script.
- New country-specific values go in `00` as `switch(COUNTRY, ...)` arms, never hard-coded in 01–08.
