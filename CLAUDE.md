# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Rheumatic Heart Disease (RHD) secondary-prevention investment case — a population-level Markov
state-transition model (WHF stages A/B/C/D) projecting RHD burden and evaluating echo screening +
diagnosis + secondary antibiotic prophylaxis (SAP) scale-up, with benefit-cost / cost-effectiveness.
Originally built for **Indonesia**; now **multi-country** (Indonesia, Uganda) via a single switch.
Descends from the NCD Countdown CCPM framework (legacy helpers remain in `R/utils.R`, unused).

**Current primary specification** (set in `00`): `mortality_structure = "D_only"` (RHD death only
from stage D) and `sap_mechanism = "progression"` (SAP lowers progression hazards, no direct mortality
effect). `"all_stage"` and `"direct_mortality"` are labelled legacy comparators — never silent
fallbacks. Many effect sizes and priors are tagged **PROVISIONAL** pending clinical review.

## Running the Pipeline

```r
source("scripts/00_run_all.R")   # single control panel: sets every parameter, sources 01 → 08
```

Run once per `COUNTRY`. `04b` dominates runtime (~11 min per country solo; much slower with several
R processes in parallel). For partial runs, evaluate `00` without block L, then source the scripts you
need. Scripts can be sourced standalone (each falls back to defaults via a `getp()`/`exists()` guard
and `here::here()`). **Gotcha:** standalone fallbacks point at the *root* `data/` and `output/`, not
the country folders — set `wd_data`/`wd_outp` to `data/<COUNTRY>/`, `output/<COUNTRY>/` first.

| Script | Purpose | Writes (under `data/<COUNTRY>/` or `output/<COUNTRY>/`) |
|---|---|---|
| `01_prepare_inputs.R` | Read every CSV in `data-raw/epidemiology/` + `population/`, keep RHD + All causes for `LOCATION` (retains GBD `upper`/`lower` UIs) | `temp_baseline_rates_gbd.rds`, `temp_population_gbd.rds` |
| `02_build_demography.R` | WPP2024 population backbone: observed `OBS_YEARS`, medium-variant `PROJ_YEARS` | `pop_observed.*`, `pop_projection.*`, `wpp/rhd_demography.Rda` |
| `03_build_disease_model.R` | **Input builder only**: rate arrays + A/B/C/D priors + mortality structure + SAP/surgery effects + cascade; clips horizon to `ANALYSIS_YEARS` | `disease_model_inputs.rds` |
| `04_calibration_random_tp.R` | **Layer 1** aggregate IR/CF calibration (2000–2019) **and** builds the Layer-2 target/input grids | `calibrated_rhd_parameters.rds` (`$layer2_targets`, `$layer2_inputs`, …) + `calibration_layer1_*.csv`, `calibration_layer2_*.csv` |
| `04b_calibrate_structural.R` | **Layer 2** joint structural calibration on the production engine; bridges calibrated state to the analysis start | `$stage_calibration` in the bundle, `calibrated_rhd_parameters_last_valid.rds`, `calibration_structural_*.csv` |
| `05_build_baseline.R` | Assemble initial state from 04b (`best_parameters`, `initial_state`); realise cascade/adherence trajectories; spec-consistency checks | `baseline_state.rds` |
| `06_run_prevention_model.R` | **Run** ref vs SAP arms (parallel by location); no costs | `out_model/<location>.rds` (`$wsd`, `$stages`, `$diag`, `$meta`) |
| `07_make_outputs.R` | Long tables | `tables/rhd_model_long.*`, `rhd_stage_model_long.*`, `rhd_annual_flows.*` |
| `08_economic_evaluation.R` | **All economics** (BCR, cost/death, cost/DALY, budget impact) | `tables/rhd_economic_summary.csv`, `rhd_budget_impact.csv`, `rhd_economic_results.rds` |

### Tests

```r
Rscript tests/test_abcd_smoke.R          # deterministic engine test (mass balance, hazards, SAP, surgery)
Rscript tests/verify_abcd_pipeline.R     # 23 spec checks on persisted 03–08 outputs (run after pipeline)
```

Both exit non-zero on failure. `verify_abcd_pipeline.R` defaults to root `data/`/`output/` — set
`wd_data`/`wd_outp` to the country folders first. The smoke test exercises only the legacy
`direct_mortality` SAP path; the `progression` / `D_only` paths have no unit test.

### Report and deck (read saved outputs only; never re-run the model)

- `analysis/report.RMD` — HTML report, `params = list(country = "Indonesia" | "Uganda")`.
- `reports/slideck-technical.Rmd` — Beamer/XeLaTeX technical deck covering both countries; uses
  `reports/beamer_preamble.tex`, `reports/apa.csl`, `library/references.bib`. Keep deck-specific LaTeX
  in its YAML `header-includes`. Re-render after any pipeline re-run.
- Both honour env var `RHD_REPO_ROOT`; render copies outside the repo (copy Rmd + preamble + csl to
  `<tmp>/reports/`, bib to `<tmp>/library/`) to keep the repo clean. Integrity `check()` chunks stop the
  render if outputs are inconsistent. Use ASCII (not →, ↓, −) in ggplot text.

### Windows / CLI notes

- R is not on PATH; use R 4.5.1 (`C:/Program Files/R/R-4.5.1/bin/Rscript.exe`) — `wpp2024` segfaults under 4.3.
- Load `data.table` **before** `wpp2024` (OpenMP segfault otherwise); `00` and `02` already do.
- Don't pipe Rscript output through `head`/`tail` (masks exit code); redirect to a file instead.
- CLI rendering needs RStudio's pandoc:
  `Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")` before
  `rmarkdown::render()`. The deck needs XeLaTeX (TinyTeX).

## Parameters (`00_run_all.R`, blocks A–L)

- **A** `COUNTRY` (`"Indonesia"` | `"Uganda"`) → `ISO3`, paths. `data/<COUNTRY>/`, `output/<COUNTRY>/`
  are country-scoped; `data-raw/` is shared. Adding a country means adding an arm to every `switch(COUNTRY, ...)`.
- **B** `LOCATION`/`LOCATIONS`, per-country sanity bands (`POP2020_*`, `RHD_PREV_*`, `RHD_DEATH_*`,
  `ALLC_DEATH_*`) — order-of-magnitude gates; don't loosen globally.
- **C** Year windows, kept conceptually separate: analysis `ANALYSIS_YEAR_START/END` (default
  **2026–2100**), observed `OBS_YEARS` (1990–2024), projection `PROJ_YEARS` (contiguous, ends at analysis
  end), calibration `CAL_YEAR_*` (2000–2019), `RATE_BASE_YEAR` (2023), `gdp_pc_base_year`,
  `discount_base_year` (= analysis start). `stopifnot` checks enforce consistency.
- **D** Intervention ramp `ramp_start`/`ramp_end` (2026→2050). **E** `incidence_trend`.
- **F** Layer-1 calibration (weights, `USE_IV_WEIGHTS`, `INCIDENCE_CALIBRATION_MODE`, `CALIB_AGE_BANDS`,
  priors, `MULT_LO/HI`, `CALIB_OPTIMIZER` = `multistart` | `nloptr` | `deoptim`, `N_STARTS`, `SEED`).
- **G** Natural history: transition prior centers; `mortality_structure`; `p_rhd_death_D_only_prior`
  (provisional, from `dev/TP Stage D to death_RHDparameters.xlsx` col K); legacy `p_rhd_death_A..D`
  (`all_stage` only); `rhd_d_fraction`; age structure (`remission_adult_mult`, `remission_age_cutoff`,
  `STRUCT_CALIB_AGE_STRUCTURE` — two adult multipliers on progression and remission at ages ≥ 20).
- **G2** Layer-2 (04b): prevalence-led weights (`STRUCT_PREV_DEATH_RATIO` = 2.5, `W_INC_STRUCT` = 0),
  priors/λ (`LAMBDA_PRIOR_STRUCT` = 0.05, `LAMBDA_PRIOR_MORT`), `STRUCT_BURN_YEARS` (120, checked vs 2×),
  `STRUCT_PREV_OBS = "mid"`, `STRUCT_TARGET_BASIS = "rate_x_pop"`, `STRUCT_DEATH_TARGET_AGE_MIN` (20),
  fit tolerances, `MD_LO/HI`, `TRANS_BOUNDS`, `STRUCT_BRIDGE_TO_ANALYSIS_START`.
- **H** Effects: `sap_mechanism`; GOAL-trial-derived hazard ratios `sap_hr_A_to_B`, `sap_hr_B_to_C`
  (from `dev/Summary progression rates_Sep26Update.xlsx`), `sap_hr_C_to_D = 0.45` (Coates 2021
  assumption); adherence modifiers per arm (default 1); legacy `sap_rrr_rhd_death = 0.55`; surgery RRRs.
- **I** Care cascade (screen/diagnosis/treatment base + target) and surgery service (requirement
  fractions, coverage per arm); optional `screen_age_restrict`.
- **J** Economics — the ONLY monetary values (consumed by 08 only). **K** `run_model_par`, `MAX_CORES`.
  **L** sources 01 → 08.

Commented-out blocks in `00` are superseded values kept for reference.

> **Provenance discipline**: scripts 01–07 contain NO monetary values. Calibrated parameters are baked
> in ONCE and consumed as-is by 05/06 — never re-adjust downstream. A calibration fitted under one
> `mortality_structure`/`sap_mechanism` is never consumed under another (05 stops). The former hand-tuned
> mortality scalars (1/2.5 Indonesia, 1/8 Uganda) are gone; do not reintroduce them.

## Architecture

### Directory layout

```
scripts/          pipeline 00–08 (+ 04b)
R/                abcd_engine.R (shared engine), packages.R (dependency check), utils.R (legacy, unused)
tests/            test_abcd_smoke.R, verify_abcd_pipeline.R
analysis/         report.RMD (HTML, parameterised by country)
reports/          slideck-technical.Rmd/.pdf, beamer_preamble.tex, apa.csl
library/          references.bib + key papers (WHF 2023 guidelines, Coates 2021 + appendix)
docs/             metadata_inputs.csv (input-data catalogue), state-diagram PNGs + generator .py, notes
dev/              task prompts, parameter workbooks (.xlsx) feeding provisional priors, meeting notes
data-raw/         GBD 2023 epidemiology + population CSVs (not tracked; shared across countries)
data/<COUNTRY>/   processed inputs, calibration bundle + diagnostics CSVs, baseline state
output/<COUNTRY>/ out_model/ (gitignored) + tables/
```

Files directly under `data/`, `output/tables/`, and `data-raw/*.rds`, plus year-suffixed
`pop_*_YYYY_YYYY.*` files, are stale pre-multi-country artifacts — the pipeline does not read them.

### Disease structure and engine

No RHD → A ↔ B ↔ C → D → RHD death, with competing other-cause death from every living stage and the
No-RHD pool; incident RHD enters **A** (no ARF state). A minimal/early · B mild · C advanced w/o
complications · D advanced with HF. **`R/abcd_engine.R`** is the single source of truth, shared by 04b
and 06: `abcd_one_cycle` (one annual cycle), `abcd_compete` (hazard-scale competing risks:
`p_any = 1 − exp(−Σh)`, split by hazard share — no clipping, exact mass balance), `abcd_sap_effects`,
`abcd_age_shift`, `abcd_burn_in`, `abcd_project` (start stock → flows → end stock → ageing).
Transitions may be scalars or age×sex matrices.

- **SAP `progression`**: `h_ij(c) = h_ij(0) × [1 − c(1 − HR_ij)]` on A→B, B→C, C→D, with
  c = effective treatment coverage × adherence; `eff_cf = 1`. Legacy `direct_mortality`: every stage's
  RHD-death probability × `(1 − sap_rrr_rhd_death × c)`.
- **D_only**: `pd[A] = pd[B] = pd[C] = 0` exactly, so RHD deaths in A–C are exactly zero.
- **Surgery is a SERVICE, not a state**: reach = requirement fraction × coverage; reduces `C→D` and
  `D→death` only, never removes stock; held equal in both arms.
- Calibrated parameters are **untreated** natural history; reference care is applied inside the engine
  in both calibration and projection. All-cause deaths are reconstructed from flows.
- Only the care cascade (and adherence) differs between `ref` and `sap`.

06 writes `$wsd` (aggregate, `sick = A+B+C+D`), `$stages` (stocks, flows, per-stage deaths, cascade
volumes, surgery trace; collapses exactly to `$wsd`), and `$diag` (deaths averted, D start share,
well deaths, effective protection, realised D fatality, …).

### Calibration — two layers

**Layer 1 (04):** per location × sex, log-scale incidence (α) and aggregate RHD-mortality (β)
multipliers per broad age band for a well-sick-dead proxy, fitted to GBD RHD incidence, prevalence and
RHD deaths (priors + smoothness; bounded multi-start L-BFGS-B). Its α are **starting values** for 04b.
04 also builds the Layer-2 targets (0–100 grid, 95+ aggregated, missing cells masked; GBD Number and
Rate × WPP population side by side) and validates any supplied stage-target table.

**Layer 2 (04b):** jointly estimates incidence multipliers (band × sex), transitions, the stage-D
mortality level and two adult age multipliers against GBD **prevalence (principal target) and RHD
deaths at ages ≥ 20**; incidence is validation-only. Burn-in to equilibrium is checked, then the
end-of-window state is bridged to `ANALYSIS_YEAR_START` and saved as `initial_state` (06's seed — no
re-seeding from GBD prevalence). Separate structural / numerical / fit statuses; a candidate failing
structural validation is never written as calibrated (status `FAILED_*`, script stops,
`_last_valid.rds` kept). Stage split is partially identified: near-optimal sets and an
identifiability table are saved. Known limits: GBD child (1–19) RHD deaths are unreachable under
D-only; Indonesia's child/adult prevalence age pattern is only partly reproduced.

### Demography (02)

Built directly from `wpp2024`: observed table spline-smoothed with totals preserved (2024 is taken from
`popprojAge1dt`, WPP's jump-off year), projection taken as-is. Stable filenames (no year in the name).

## Dependencies

R ≥ 4.1 (native pipe). `R/packages.R` checks: here, dplyr, tidyr, readr, purrr, tibble, ggplot2,
scales, abind, stringr, readxl, cowplot, data.table, countrycode, rlang. Also used: foreach,
doParallel, parallel (06); optional nloptr / DEoptim (04 optimisers); knitr, rmarkdown, kableExtra,
DT (report); XeLaTeX/TinyTeX (deck). `wpp2024` is GitHub-only: `devtools::install_github("PPgp/wpp2024")`.

## Conventions

- `data.table` throughout pipeline scripts; parameters tagged `[PAPER]`, `[LIT]`, `[CALIBRATE]`,
  `[PROVISIONAL]`, `[TODO: confirm for Uganda]`. Keep a value's actual source next to it.
- Scripts communicate only via persisted files — never rely on in-memory objects from an earlier script.
- New country-specific values go in `00` as `switch(COUNTRY, ...)` arms, never hard-coded in 01–08.
- Never relabel a legacy or sensitivity comparator (`all_stage`, `direct_mortality`, all-ages death
  targets) as primary without recalibrating.
