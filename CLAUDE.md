# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Rheumatic Heart Disease (RHD) secondary-prevention investment case — a population-level Markov
state-transition model (WHF stages A/B/C/D) that separates **disease history from first confirmed
diagnosis**, projects RHD burden and evaluates echo screening + confirmation + secondary antibiotic
prophylaxis (SAP/BPG) scale-up, with benefit-cost / cost-effectiveness. Multi-country (Indonesia,
Uganda) via a single switch. Descends from the NCD Countdown CCPM framework (legacy helpers remain in
`R/utils.R`, unused).

**Current engine**: `ABCD_ENGINE_SPEC = "abcd_first_diagnosis_v1_adult_prog_mult"` (`R/abcd_engine.R`).
Fixed specification: `mortality_structure = "D_only"` and `sap_mechanism = "progression"` — 03 now
**rejects** any other value (the old `all_stage` / `direct_mortality` comparators are gone). Many inputs
are tagged **PROVISIONAL** pending clinical review.

**Status (Nov 2026)**: 04b's epidemiological-fit checkpoint FAILS for both countries (structural and
numerical validity pass; adult RHD deaths ~3–4× GBD in Indonesia). 05–08 run under an explicit user
override (`PIPELINE_RUN_THROUGH <- "08"`), 05 labels the baseline `INADEQUATE`, and all outputs are
**preliminary**. Details: `dev/checkpoint_first_diagnosis_2026-11-05.md`.

## Running the Pipeline

```r
source("scripts/00_run_all.R")   # single control panel: sets every parameter, sources 01 → 08
```

Run once per `COUNTRY`. `04b` dominates runtime (~15–20 min per country). `PIPELINE_RUN_THROUGH`
(`"04b"` | `"08"`) stops after the structural calibration or runs the full pipeline. For partial runs,
evaluate `00` without block L, then source the scripts you need. Scripts can be sourced standalone
(`getp()`/`exists()` fallbacks + `here::here()`). **Gotcha:** standalone fallbacks point at the *root*
`data/` and `output/` — set `wd_data`/`wd_outp` to `data/<COUNTRY>/`, `output/<COUNTRY>/` first.

| Script | Purpose | Writes (under `data/<COUNTRY>/` or `output/<COUNTRY>/`) |
|---|---|---|
| `01_prepare_inputs.R` | Read every CSV in `data-raw/epidemiology/` + `population/`, keep RHD + All causes for `LOCATION` (retains GBD UIs) | `temp_baseline_rates_gbd.rds`, `temp_population_gbd.rds` |
| `02_build_demography.R` | WPP2024 backbone: observed `OBS_YEARS`, medium-variant `PROJ_YEARS` | `pop_observed.*`, `pop_projection.*`, `wpp/rhd_demography.Rda` |
| `03_build_disease_model.R` | **Input builder only**: rate arrays, natural history, SAP HRs, care/screening inputs, surgery; requires every block G/H/I setting from 00 (no hidden defaults); stamps `engine_spec` | `disease_model_inputs.rds` |
| `04_calibration_random_tp.R` | **Layer 1** aggregate IR/CF proxy calibration (2000–2019) + builds Layer-2 target grids | `calibrated_rhd_parameters.rds`, `calibration_layer1_*.csv`, `calibration_layer2_*.csv` |
| `04b_calibrate_structural.R` | **Layer 2** joint calibration on the production engine; bridges the full engine state to the analysis start | `$stage_calibration` in the bundle, `calibrated_rhd_parameters_last_valid.rds`, `calibration_structural_*.csv` |
| `05_build_baseline.R` | Seed from 04b `initial_state$engine_state` (no reseeding); rejects bundles from another engine / SAP spec; fit-status label | `baseline_state.rds` |
| `06_run_prevention_model.R` | **Run** ref vs scale-up arms (parallel by location); no costs | `out_model/<location>.rds` (`$wsd`, `$stages`, `$diag`, `$meta`) |
| `07_make_outputs.R` | Long tables | `tables/rhd_model_long.*`, `rhd_stage_model_long.*`, `rhd_annual_flows.*` |
| `08_economic_evaluation.R` | **All economics** (BCR, cost/death, cost/DALY, budget impact) | `tables/rhd_economic_summary.csv`, `rhd_budget_impact.csv`, `rhd_economic_results.rds` |

Every stage checks `engine_spec`: 04b stops if `disease_model_inputs.rds` was built for another engine,
05 if the bundle was — re-run from 03 after any engine change.

### Tests

`tests/test_abcd_smoke.R` and `tests/verify_abcd_pipeline.R` are **stale**: written for the pre-November
engine signature (`abcd_one_cycle(A, B, C, D, ...)`) and old output columns; the smoke test currently
fails with "unused arguments". Validation lives inside the scripts (04b structural/accounting checks,
05/06/07/08 internal `stop()`s, report/deck `check()` chunks). Rewrite the tests before relying on them.

### Report and deck (read saved outputs only; never re-run the model)

- `analysis/report.RMD` — HTML report, `params = list(country = "Indonesia" | "Uganda")`.
- `reports/slideck-technical.Rmd` — Beamer/XeLaTeX deck covering both countries; uses
  `reports/beamer_preamble.tex`, `reports/apa.csl`, `library/references.bib`. Its `ENGINE` constant
  (≈ line 109) must equal `ABCD_ENGINE_SPEC` or the render stops. Keep deck-specific LaTeX in its YAML
  `header-includes`. Re-render after any pipeline re-run.
- Both honour env var `RHD_REPO_ROOT`; render copies outside the repo (Rmd + preamble + csl to
  `<tmp>/reports/`, bib to `<tmp>/library/`). Integrity `check()` chunks stop the render if outputs are
  inconsistent. Use ASCII (not →, ↓, −) in ggplot text.

### Windows / CLI notes

- R is not on PATH; use R 4.5.1 (`C:/Program Files/R/R-4.5.1/bin/Rscript.exe`) — `wpp2024` segfaults under 4.3.
- Load `data.table` **before** `wpp2024` (OpenMP segfault otherwise); `00` and `02` already do.
- Don't pipe Rscript output through `head`/`tail` (masks exit code); redirect to a file instead.
- CLI rendering needs RStudio's pandoc:
  `Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")` before
  `rmarkdown::render()`. The deck needs XeLaTeX (TinyTeX).
- The repo is on OneDrive: transient file locks can fail an `.rds` write — just re-run.

## Parameters (`00_run_all.R`, blocks A–L)

- **A** `COUNTRY` (`"Indonesia"` | `"Uganda"`) → `ISO3`, paths. `data/<COUNTRY>/`, `output/<COUNTRY>/`
  are country-scoped; `data-raw/` is shared. Adding a country means adding an arm to every `switch(COUNTRY, ...)`.
- **B** `LOCATION`/`LOCATIONS`, per-country sanity bands (`POP2020_*`, `RHD_PREV_*`, `RHD_DEATH_*`,
  `ALLC_DEATH_*`) — order-of-magnitude gates; don't loosen globally.
- **C** Year windows, kept separate: analysis `ANALYSIS_YEAR_START/END` (**2026–2100**), observed
  `OBS_YEARS` (1990–2024), projection `PROJ_YEARS` (contiguous, ends at analysis end), calibration
  `CAL_YEAR_*` (2000–2019), `RATE_BASE_YEAR` (2023); economic anchors in J. `stopifnot` enforces consistency.
- **D** Ramp `ramp_start`/`ramp_end` (2026→2050). **E** `incidence_trend`.
- **F** Layer-1 calibration (weights prevalence-led 1/2.5/1, `USE_IV_WEIGHTS`, `INCIDENCE_CALIBRATION_MODE`,
  `CALIB_AGE_BANDS`, priors, `MULT_LO/HI`, `CALIB_OPTIMIZER` = `multistart` | `nloptr` | `deoptim`).
- **G** Natural history (annual probabilities → hazards once): `nh_child_last_age = 20` (child rules at
  attained age ≤ 20, adult > 20); six calibrated childhood transitions `p_A_to_no_rhd`, `p_A_to_B`,
  `p_B_to_A`, `p_B_to_C`, `p_C_to_B`, `p_C_to_D` (prior centres); fixed `p_B_to_no_rhd_adult`,
  `entry_mix_child` (A) / `entry_mix_adult` (B), `p_rhd_death_D` (per country, fixed, not fitted),
  **`adult_prog_hazard_mult = 0.48`** (adult B→C, C→D: `p_adult = 1 − (1 − p_child)^m`; adopted 6 Nov 2026);
  diagnostic-only `dx_mix_benchmark_*` and `d_share_proxy_gbd_hf`.
- **G2** Layer-2 (04b): `STRUCT_MATCH_FROM = 2010` (first scored year), `STRUCT_BURN_YEARS = 120`,
  `STRUCT_BURNIN_RECONCILIATION = "per_capita"`, prevalence-led objective (`STRUCT_PREV_DEATH_RATIO` 2.5,
  `W_INC_STRUCT = 0`, 0–19/20+ split shares 0.5, 21+ diagnostic), `LAMBDA_PRIOR_STRUCT = 0.05`,
  `TRANS_BOUNDS`, fit tolerances (`STRUCT_FIT_TOL_PREV` 0.15, `_DEATH` 0.25 on |log ratio|),
  `STRUCT_N_STARTS`, `STRUCT_PARALLEL_STARTS`, `STRUCT_BRIDGE_TO_ANALYSIS_START`.
- **H** Effects: GOAL-trial child HRs `sap_hr_A_to_B_child`, `sap_hr_B_to_C_child` (computed from arm
  risks); adult B→C, prior-A, C→D and D-death HRs all 1 (no adult benefit, no mortality effect);
  surgery RRRs on C→D and D→death.
- **I** Screening + care: one screen per person across ages 5–39 among the never-screened undiagnosed
  (`screen_sensitivity`, `confirm_completion`, `confirm_sensitivity`, specificity fixed at 1);
  `care_inputs` (screen window prob, clinical detection, BPG initiation/retention, adherence) per arm as
  `c(baseline, target)`; initial care shares; `sap_eligibility = c("A", "B")`; surgery service fractions/coverage.
- **J** Economics — the ONLY monetary values (08 only): `disc_rate`, VSL = `vsl_mult` × `gni_pc_base`
  (WDI 2019) grown with `gdp_growth`, `dalys_per_death`, unit costs (screen, SAP-year, surgery).
  Confirmatory echoes have **no unit cost** — 08 flags them as uncosted.
- **K** `run_model_par`, `MAX_CORES` (also caps 04b parallel starts). **L** `PIPELINE_RUN_THROUGH` + sources.

Commented-out blocks in `00` are superseded values kept for reference.

> **Provenance discipline**: scripts 01–07 contain NO monetary values. Calibrated parameters are baked
> in ONCE and consumed as-is by 05/06 — never re-adjust downstream. Never reintroduce hand-tuned
> mortality scalars or lower `p_rhd_death_D` to manufacture fit.

## Architecture

### Directory layout

```
scripts/          pipeline 00–08 (+ 04b)
R/                abcd_engine.R (shared engine), packages.R (dependency check), utils.R (legacy, unused)
tests/            test_abcd_smoke.R, verify_abcd_pipeline.R (both stale, see above)
analysis/         report.RMD (HTML, parameterised by country)
reports/          slideck-technical.Rmd/.pdf, beamer_preamble.tex, apa.csl
library/          references.bib + key papers
docs/             metadata_inputs.csv (input catalogue), state-diagram PNGs + generator .py, notes
dev/              task prompts/specs, checkpoint reports, parameter workbooks (.xlsx), meeting notes
data-raw/         GBD 2023 epidemiology + population CSVs (not tracked; shared across countries)
data/<COUNTRY>/   processed inputs, calibration bundle + diagnostics CSVs, baseline state
output/<COUNTRY>/ out_model/ + tables/
```

Stale / not read by the pipeline: files directly under `data/` and `output/tables/`, `data-raw/*.rds`,
year-suffixed `pop_*_YYYY_YYYY.*`, and archived bundles
`data/<COUNTRY>/calibrated_rhd_parameters_last_valid_prev_engine_*.rds` (kept for comparison).

### Engine (`R/abcd_engine.R`, single source of truth for 04b and 06)

State per age × sex cell (person counts): never-onset `N_ns/N_sc`; undiagnosed classifiable disease
`U_<s>_ns/_sc` (by screening history); diagnosed off BPG `X_<s>`; on BPG `T_A`, `T_B`; adult prior-A
`P_ns/P_sc/P_dx`; normalised history `R_ns/R_sc/R_dx`. **Classifiable prevalence (GBD mapping) = U + X + T**;
prior-A and normalised are excluded.

- **Natural history by attained age**: ≤ 20: A↔B↔C→D, A→normalised; > 20: no A, prior-A→B re-entry
  (shared A→B), B→normalised (fixed), B→C and C→D scaled by the adult hazard multiplier. D→C = 0.
- **New disease** enters A (≤ 20) or B (> 20), undiagnosed; screening never creates disease.
- **Ageing 20→21**: surviving A becomes prior-A (BPG ends; reclassification, not remission).
- **D_only**: RHD death only from D (A/B/C exact zeros); competing other-cause death everywhere;
  hazard-scale competing risks (`abcd_compete`: `p_any = 1 − exp(−Σh)`, exact mass balance).
- **Care update** each cycle: screening → true positives → confirmation → first confirmed diagnosis
  (U→X, at most once per person); clinical detection; BPG initiation (eligible: A ≤ 20, B), retention,
  dropout. On-BPG protection `h × [1 − adh × (1 − HR)]`, `adh = good + rel × partial`.
- **Surgery is a SERVICE, not a state**: reduces C→D and D→death; held equal in both arms.
- Calibrated parameters are **untreated** natural history; reference care is applied inside the engine
  in burn-in, calibration, bridge and projection. Only care inputs differ between `ref` and scale-up.
- Key functions: `abcd_nh_build`, `abcd_care_schedule`, `abcd_one_cycle(state, popm, irm, oth, nhe, eff,
  surg, care, light)`, `abcd_age_state`, `abcd_burn_in`, `abcd_project`. Full and `light` modes share one path.

06 writes `$wsd` (15-column aggregate; `sick` = classifiable disease, `newcases` = new modelled disease,
not diagnoses), `$stages` (stocks by care state, first diagnoses by route/stage, confirmations, re-entry,
reclassification, BPG flows, D deaths by diagnosis status), `$diag` (deaths averted, `care_state_by_year`,
D fatality by rule age, …).

### Calibration — two layers

**Layer 1 (04):** per location × sex, log-scale incidence (α) and aggregate mortality (β) multipliers
per broad age band for a well-sick-dead proxy vs GBD incidence, prevalence, RHD deaths. Its α are
**starting values** for 04b; its CF is diagnostic only. Also builds Layer-2 targets (rate × WPP pop and
GBD Number side by side).

**Layer 2 (04b):** runs the real engine (zero-start per-capita burn-in on the 2000 cross-section, then
2000–2019) and jointly fits band × sex incidence multipliers + the six childhood transitions against GBD
**classifiable prevalence (principal) and RHD deaths**, scored 2010–2019 with 0–19 / 20+ components;
incidence is validation-only. Separate structural / numerical / epidemiological-fit statuses; a candidate
failing structural validation is never written as calibrated (`FAILED_*`, script stops, `_last_valid.rds`
kept), while a fit outside tolerance is recorded (`__fit_outside_diagnostic_tolerance`) and passed on.
The calibrated end state is bridged to `ANALYSIS_YEAR_START` as `initial_state$engine_state`. Exports
near-optimal sets, identifiability, loss decomposition, first-diagnosis stage mix, stage-by-age tables.
Known limits: shared progression vs Indonesia's high childhood case fatality (A→B, C→D at upper bounds);
prior-A exclusion creates a 20–24 prevalence dip; GBD case fatality falls over 2010–2019.

### Demography (02)

Built from `wpp2024`: observed table spline-smoothed with totals preserved (2024 from `popprojAge1dt`),
projection taken as-is. Stable filenames (no year in the name).

## Dependencies

R ≥ 4.1 (native pipe). `R/packages.R` checks: here, dplyr, tidyr, readr, purrr, tibble, ggplot2,
scales, abind, stringr, readxl, cowplot, data.table, countrycode, rlang. Also used: foreach,
doParallel, parallel (04b, 06); optional nloptr / DEoptim (04); knitr, rmarkdown, kableExtra, DT
(report); XeLaTeX/TinyTeX (deck). `wpp2024` is GitHub-only: `devtools::install_github("PPgp/wpp2024")`.

## Conventions

- `data.table` throughout pipeline scripts; parameters tagged `[PAPER]`, `[LIT]`, `[CALIBRATE]`,
  `[PROVISIONAL]`, `[ASSUMPTION]`, `[TODO: confirm for Uganda]`. Keep a value's actual source next to it.
- Scripts communicate only via persisted files — never rely on in-memory objects from an earlier script.
- New country-specific values go in `00` as `switch(COUNTRY, ...)` arms, never hard-coded in 01–08.
- Changing the engine means bumping `ABCD_ENGINE_SPEC`, re-running 03 → 08 for both countries, and
  updating the deck's `ENGINE` constant; archive the previous last-valid bundle as `*_prev_engine_<spec>.rds`.
- Never present preliminary / INADEQUATE-fit results as a credible calibrated baseline.
