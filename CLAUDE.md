# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Rheumatic Heart Disease (RHD) secondary-prevention investment case — a population-level Markov
state-transition model (WHF stages A/B/C/D) that separates **disease history from first confirmed
diagnosis**, projects RHD burden and evaluates echo screening + confirmation + secondary antibiotic
prophylaxis (SAP/BPG) scale-up, with benefit-cost / cost-effectiveness. Multi-country (Indonesia,
Uganda) via a single switch. Descends from the NCD Countdown CCPM framework (legacy helpers remain in
`R/utils.R`, unused).

**Current engine**: `ABCD_ENGINE_SPEC = "abcd_first_diagnosis_v1_adult_prog_mult_sap55_dmort_age"`
(`R/abcd_engine.R`; `_dmort_age` = fixed age-specific stage-D mortality, 9 Nov 2026). Fixed
specification: `mortality_structure = "D_only"` and `sap_mechanism = "progression"`; 03 **rejects** any
other value. Many inputs are tagged **PROVISIONAL** pending clinical review.

**Status (9 Nov 2026)**: 04b's epidemiological-fit checkpoint FAILS for both countries (structural and
numerical validity pass). Fit, 2010–2019 model/GBD: Indonesia prevalence ~0.73×, RHD deaths 20+
1.5–2.0×, 0–19 0.5–0.7×; Uganda prevalence ~0.6×, deaths 3.2–4.1× at all ages. 05–08 run under an
explicit user override (`PIPELINE_RUN_THROUGH <- "08"`), 05 labels the baseline `INADEQUATE`, and all
outputs are **preliminary**. Background: `dev/checkpoint_first_diagnosis_2026-11-05.md`.

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
| `03_build_disease_model.R` | **Input builder only**: rate arrays, natural history (incl. `[age × sex]` D mortality), SAP progression multiplier, care/screening inputs, surgery; requires every block G/H/I setting from 00 (no hidden defaults); stamps `engine_spec` | `disease_model_inputs.rds` |
| `04_calibration_random_tp.R` | **Layer 1** aggregate IR/CF proxy calibration (2000–2019) + builds Layer-2 target grids | `calibrated_rhd_parameters.rds`, `calibration_layer1_*.csv`, `calibration_layer2_*.csv` |
| `04b_calibrate_structural.R` | **Layer 2** joint calibration on the production engine; bridges the full engine state to the analysis start | `$stage_calibration` in the bundle, `calibrated_rhd_parameters_last_valid.rds`, `calibration_structural_*.csv` |
| `05_build_baseline.R` | Seed from 04b `initial_state$engine_state` (no reseeding); rejects bundles from another engine / SAP spec / D-mortality inputs; fit-status label | `baseline_state.rds` |
| `06_run_prevention_model.R` | **Run** ref vs scale-up arms (parallel by location); no costs | `out_model/<location>.rds` (`$wsd`, `$stages`, `$diag`, `$meta`) |
| `07_make_outputs.R` | Long tables | `tables/rhd_model_long.*`, `rhd_stage_model_long.*`, `rhd_annual_flows.*` |
| `08_economic_evaluation.R` | **All economics** (BCR disc/undisc, cost/death, cost/DALY, budget impact, per capita, INT$) | `tables/rhd_economic_summary.csv`, `rhd_budget_impact.csv`, `rhd_economic_results.rds`, `rhd_economic_annual_per_capita.csv`, `rhd_vsl_derivation.csv`, `rhd_price_adjustment_audit.csv`, `rhd_economic_parameter_register.csv`; once both countries exist, `output/rhd_country_comparison_2025.{csv,rds}` |

Every stage checks `engine_spec`: 04b stops if `disease_model_inputs.rds` was built for another engine,
05 if the bundle was — re-run from 03 after any engine change.

### Tests

`tests/test_abcd_smoke.R` and `tests/verify_abcd_pipeline.R` are **stale**: written for the pre-November
engine signature (`abcd_one_cycle(A, B, C, D, ...)`) and old output columns; the smoke test currently
fails with "unused arguments". Validation lives inside the scripts (04b structural/accounting checks,
05/06/07/08 internal `stop()`s, report/deck `check()` chunks). Rewrite the tests before relying on them.

### Reports, decks, poster (read saved outputs only; never re-run the model)

- `analysis/report.RMD` — HTML report, `params = list(country = "Indonesia" | "Uganda")`; rendered copies
  `analysis/report_Indonesia.html`, `analysis/report_Uganda.html`.
- `reports/slideck-technical.Rmd` — Beamer/XeLaTeX technical deck, both countries; uses
  `reports/beamer_preamble.tex`, `reports/apa.csl`, `library/references.bib`. Its hard-coded `ENGINE`
  constant (≈ line 109) must equal `ABCD_ENGINE_SPEC` or the render stops. Keep deck-specific LaTeX in its
  YAML `header-includes`.
- `reports/slideck-conference.Rmd` — 15-min conference deck, results recomputed for 2026–2050 only, no
  "preliminary" labels (user's choice). biblatex APA: pass options via YAML `biblatexoptions`, never a
  second `\usepackage{biblatex}`; any content before the first `##` creates a blank frame.
- `reports/poster.Rmd` — self-contained HTML of copyable content for `reports/WCRHD2026_POSTER.pptx`
  (2026–2050). Shows income-group labels only (Uganda → LIC, Indonesia → UMIC). Writes a sanitised copy
  of the bib to `tempdir()` (never edits `references.bib`); falls back to
  `calibrated_rhd_parameters_last_valid.rds` only if it matches `baseline_state.rds` exactly.
- Updated for `_dmort_age` and re-rendered 9 Nov 2026: they read `nh$p_rhd_death_D_child` / `_adult` /
  `d_death_adult_age_min`, never a single D value. Both decks allow a per-country child entry mix and
  adult multiplier (cross-country checks require only the same engine, years, discounting, SAP spec and
  adult entry mix). The conference deck's result-stating titles are `claim()`-guarded or data-derived.
- All honour env var `RHD_REPO_ROOT`; render copies outside the repo (Rmd + preamble + csl to
  `<tmp>/reports/`, bib to `<tmp>/library/`). Integrity `check()` chunks stop the render if outputs are
  inconsistent. Use ASCII (not →, ↓, −) in ggplot text. Re-render after any pipeline re-run.

### Windows / CLI notes

- R is not on PATH; use R 4.5.1 (`C:/Program Files/R/R-4.5.1/bin/Rscript.exe`) — `wpp2024` segfaults under 4.3.
- Load `data.table` **before** `wpp2024` (OpenMP segfault otherwise); `00` and `02` already do.
- Don't pipe Rscript output through `head`/`tail` (masks exit code); redirect to a file instead.
- CLI rendering needs RStudio's pandoc:
  `Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")` before
  `rmarkdown::render()`. The decks need XeLaTeX (TinyTeX).
- The repo is on OneDrive: transient file locks can fail an `.rds` write — just re-run.

## Parameters (`00_run_all.R`, blocks A–L)

- **A** `COUNTRY` (`"Indonesia"` | `"Uganda"`) → `ISO3`, paths. `data/<COUNTRY>/`, `output/<COUNTRY>/`
  are country-scoped; `data-raw/` is shared. Adding a country means adding an arm to every
  `switch(COUNTRY, ...)` and to the `*_BY_COUNTRY` lists (`TP_START_BY_COUNTRY`, `TRANS_BOUNDS_BY_COUNTRY`).
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
  `p_B_to_A`, `p_B_to_C`, `p_C_to_B`, `p_C_to_D` with **per-country** starting values / prior centres in
  `TP_START_BY_COUNTRY`; fixed `p_B_to_no_rhd_adult`, `entry_mix_child` (Indonesia A = 1; Uganda
  75% A / 25% B), `entry_mix_adult` (B); `adult_prog_hazard_mult` (Indonesia 0.48, Uganda 0.5; adult
  B→C, C→D: `p_adult = 1 − (1 − p_child)^m`). **Stage-D mortality is fixed and age-specific**:
  `p_rhd_death_D_child` at attained age < `d_death_adult_age_min` (= 20) and `p_rhd_death_D_adult` at ≥ 20
  (Indonesia 0.16 / 0.03, Uganda 0.06 / 0.03; provisional test values). The old scalar `p_rhd_death_D`
  is rejected by 03 and the engine. Diagnostic-only `dx_mix_benchmark_*` and `d_share_proxy_gbd_hf`.
- **G2** Layer-2 (04b): `STRUCT_MATCH_FROM = 2010` (first scored year), `STRUCT_BURN_YEARS = 120`,
  `STRUCT_BURNIN_RECONCILIATION = "per_capita"`, prevalence-led objective (`STRUCT_PREV_DEATH_RATIO` 2.5,
  `W_INC_STRUCT = 0`, 0–19/20+ split shares 0.5, 21+ diagnostic), `LAMBDA_PRIOR_STRUCT = 0.05`,
  per-country `TRANS_BOUNDS_BY_COUNTRY` (Uganda narrowed; starting values must lie inside or 00 stops),
  fit tolerances (`STRUCT_FIT_TOL_PREV` 0.15, `_DEATH` 0.25 on |log ratio|), `STRUCT_N_STARTS`,
  `STRUCT_PARALLEL_STARTS`, `STRUCT_BRIDGE_TO_ANALYSIS_START`.
- **H** Effects: "assumed 55% reduction in progression hazards (expert-informed scenario), informed by
  evidence on prevention of recurrent rheumatic fever" — `sap_progression_rrr = 0.55`,
  `sap_progression_hazard_mult` M = 0.45 (assumed full-adherence efficacy) on treated A→B and B→C only
  (`sap_effect_transitions`), at attained age ≤ `sap_protection_age_max` = 39 (independent of the age-20
  boundary and of eligibility; B on BPG at 40+ is costed, no benefit). Evidence is Manyemba et al. RR 0.45
  for recurrent RF (PMID 12137650) — **not** a hazard ratio or valve-progression estimate; never call M a
  "progression HR" or convert the RR with the model's B→C risk. Prior-A, C→D and D-death multipliers 1
  (no mortality effect). `sap_effect_scale = "hazard_ratio"` names the implementation scale only. GOAL-trial
  child HRs (0.0385 / 0.1998) remain commented out as a documented alternative. Surgery RRRs on C→D and D→death.
- **I** Screening + care: one screen per person across ages 5–39 among the never-screened undiagnosed
  (`screen_sensitivity`, `confirm_completion`, `confirm_sensitivity`, specificity fixed at 1);
  `care_inputs` (screen window prob, clinical detection, BPG initiation/retention, adherence) per arm as
  `c(baseline, target)`; per-country screening baseline, initial care shares and surgery service
  fractions/coverage; `sap_eligibility = c("A", "B")`.
- **J** Economics — the ONLY monetary values (08 only). Results in **constant 2025 US$** (no fallback to
  another price year); INT$ (`_int2025`) for final cross-country comparison only. `ECON_MACRO_OBS` holds
  every macro observation (WDI + BLS CPI-U, retrieved 2026-11-07): 08 makes no network calls and stops
  on a missing value. Per-country unit-cost bundles (screen, SAP-year, surgery), assumed 2019 US$,
  converted with `COST_PRICE_ADJUSTMENT_METHOD = "local_cpi_xr"` (own-country CPI + XR; never US CPI).
  VSL = reference case (Robinson et al. 2019): US $9.4m / $57,900 (2015) scaled by PPP GNI pc 2019,
  income and time elasticity 1.5, floor 20× GNI pc, grown with `gdp_growth`; no GDP add-on. `disc_rate`,
  `discount_base_year`, `dalys_per_death`, per-capita denominator N0 (WPP, first analysis year).
  `gni_pc_base`, `gdp_pc_base`, `vsl_mult` are **deprecated** (kept for report compatibility; unused by
  08). Confirmatory echoes have **no unit cost** (flagged). Sensitivity analysis: pending, not run.
- **K** `run_model_par`, `MAX_CORES` (also caps 04b parallel starts). **L** `PIPELINE_RUN_THROUGH` + sources.

Commented-out blocks in `00` are superseded values kept for reference.

> **Provenance discipline**: scripts 01–07 contain NO monetary values. Calibrated parameters are baked
> in ONCE and consumed as-is by 05/06 — never re-adjust downstream. Stage-D mortality is fixed input,
> never fitted: never reintroduce hand-tuned mortality scalars or lower D mortality to manufacture fit.

## Architecture

### Directory layout

```
scripts/          pipeline 00–08 (+ 04b)
R/                abcd_engine.R (shared engine), packages.R (dependency check), utils.R (legacy, unused)
tests/            test_abcd_smoke.R, verify_abcd_pipeline.R (both stale, see above)
analysis/         report.RMD (HTML, parameterised by country)
reports/          slideck-technical/-conference.Rmd (+ .pdf), poster.Rmd, WCRHD2026_POSTER.pptx,
                  beamer_preamble.tex, apa.csl
library/          references.bib + key papers
docs/             metadata_inputs.csv (input catalogue), state-diagram PNGs + generator .py, notes
dev/              task prompts/specs, checkpoint reports, parameter workbooks (.xlsx), meeting notes
data-raw/         GBD 2023 epidemiology + population CSVs (not tracked; shared across countries)
data/<COUNTRY>/   processed inputs, calibration bundle + diagnostics CSVs, baseline state
output/<COUNTRY>/ out_model/ + tables/
output/           rhd_country_comparison_2025.{csv,rds} (cross-country, from 08)
```

Stale / not read by the pipeline: files directly under `data/` and `output/tables/`, `data-raw/*.rds`,
year-suffixed `pop_*_YYYY_YYYY.*`, `output/<COUNTRY>/tables/_pre_vsl_update/`, and archived bundles
`data/<COUNTRY>/calibrated_rhd_parameters_last_valid_prev_engine_<spec>.rds` (kept for comparison).

### Engine (`R/abcd_engine.R`, single source of truth for 04b and 06)

State per age × sex cell (person counts): never-onset `N_ns/N_sc`; undiagnosed classifiable disease
`U_<s>_ns/_sc` (by screening history); diagnosed off BPG `X_<s>`; on BPG `T_A`, `T_B`; adult prior-A
`P_ns/P_sc/P_dx`; normalised history `R_ns/R_sc/R_dx`. **Classifiable prevalence (GBD mapping) = U + X + T**;
prior-A and normalised are excluded.

- **Natural history by attained age**: ≤ 20: A↔B↔C→D, A→normalised; > 20: no A, prior-A→B re-entry
  (shared A→B), B→normalised (fixed), B→C and C→D scaled by the adult hazard multiplier. D→C = 0.
- **New disease** enters by the entry mix (A at ≤ 20, B at > 20), undiagnosed; screening never creates disease.
- **Ageing 20→21**: surviving A becomes prior-A (BPG ends; reclassification, not remission).
- **D_only**: RHD death only from D (A/B/C exact zeros), base probability `nhe$pd` = `[age × sex]` matrix
  by current attained age (child < `d_death_adult_age_min`, adult ≥), same for every diagnosis history.
  **Two boundaries**: age 20 follows *child* progression rules (≤ 20) but *adult* D mortality (≥ 20).
  Competing other-cause death everywhere; hazard-scale competing risks (`abcd_compete`:
  `p_any = 1 − exp(−Σh)`, exact mass balance). `abcd_nh_fixed()` records the fixed inputs (incl. the
  three D-mortality fields) for stale-bundle checks 03 → 04b → 05.
- **Care update** each cycle: screening → true positives → confirmation → first confirmed diagnosis
  (U→X, at most once per person); clinical detection; BPG initiation (eligible: A ≤ 20, B), retention,
  dropout. On-BPG protection `h × [1 − adh × (1 − M)]` on A→B / B→C at age ≤ 39 (1 above),
  `adh = good + rel × partial` (reference 0.70 → 0.615). The SAP spec travels 03 → 04b (`spec$sap_hr`) →
  05 (`meta$sap_hr`) as a typed named list of 7 fields; 05 refuses incomplete / mismatched specs.
  `n_effectively_protected` (on BPG × adh) counts every age, including 40+ where there is no benefit.
- **Surgery is a SERVICE, not a state**: reduces C→D and D→death; held equal in both arms.
- Calibrated parameters are **untreated** natural history; reference care is applied inside the engine
  in burn-in, calibration, bridge and projection. Only care inputs differ between `ref` and scale-up.
- Key functions: `abcd_nh_build`, `abcd_nh_fixed`, `abcd_care_schedule`, `abcd_one_cycle(state, popm, irm,
  oth, nhe, eff, surg, care, light)`, `abcd_age_state`, `abcd_burn_in`, `abcd_project`. Full and `light`
  modes share one path.

06 writes `$wsd` (15-column aggregate; `sick` = classifiable disease, `newcases` = new modelled disease,
not diagnoses), `$stages` (stocks by care state, first diagnoses by route/stage, confirmations, re-entry,
reclassification, BPG flows, D deaths by diagnosis status), `$diag` (deaths averted, `care_state_by_year`,
D fatality by mortality age group 0–19 / 20+, …).

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
kept, `calibration_structural_failed_candidate.*` written), while a fit outside tolerance is recorded
(`__fit_outside_diagnostic_tolerance`) and passed on. The calibrated end state is bridged to
`ANALYSIS_YEAR_START` as `initial_state$engine_state`. Exports near-optimal sets, identifiability, loss
decomposition, first-diagnosis stage mix, stage-by-age tables, D fatality by mortality age group.
Known limits: transitions at bounds (Indonesia A→B, B→C, C→D; Uganda all six); prior-A exclusion creates a 20–24 prevalence
dip; GBD case fatality falls over 2010–2019; model survival at ages 90–99 exceeds WPP cohort survival,
so 95+ prevalence is far above GBD.

### Demography (02)

Built from `wpp2024`: observed table spline-smoothed with totals preserved (2024 from `popprojAge1dt`),
projection taken as-is. The open 100+ interval keeps WPP's raw count: the spline is fitted over 0–100,
then only 0–99 are rescaled to their source total (fitting 0–99 alone zeroes Uganda's male age 99).
Stable filenames (no year in the name).

## Dependencies

R ≥ 4.1 (native pipe). `R/packages.R` checks: here, dplyr, tidyr, readr, purrr, tibble, ggplot2,
scales, abind, stringr, readxl, cowplot, data.table, countrycode, rlang. Also used: foreach,
doParallel, parallel (04b, 06); optional nloptr / DEoptim (04); knitr, rmarkdown, kableExtra, DT
(report); XeLaTeX/TinyTeX + biber (decks). `wpp2024` is GitHub-only: `devtools::install_github("PPgp/wpp2024")`.

## Conventions

- `data.table` throughout pipeline scripts; parameters tagged `[PAPER]`, `[LIT]`, `[CALIBRATE]`,
  `[PROVISIONAL]`, `[ASSUMPTION]`, `[TODO: confirm for Uganda]`. Keep a value's actual source next to it.
- Scripts communicate only via persisted files — never rely on in-memory objects from an earlier script.
- New country-specific values go in `00` as `switch(COUNTRY, ...)` arms (or `*_BY_COUNTRY` lists), never
  hard-coded in 01–08.
- Changing the engine means bumping `ABCD_ENGINE_SPEC`, re-running 03 → 08 for both countries, and
  updating the report/decks (technical deck `ENGINE` constant); archive the previous last-valid bundle as
  `*_prev_engine_<spec>.rds`.
- Never present preliminary / INADEQUATE-fit results as a credible calibrated baseline.
