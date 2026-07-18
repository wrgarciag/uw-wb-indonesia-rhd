# Claude Code task: refactor the RHD A/B/C/D calibration pipeline (`scripts/00–08` + `analysis/report.RMD`)

You are refactoring a World Heart Federation A/B/C/D rheumatic-heart-disease (RHD) state-transition model. The pipeline is a country-switched, age–sex-structured annual cohort model. Read this whole brief before editing anything, then work in the staged order given in Part B. **Do not attempt the full rewrite in one pass.** There is a hard validation gate after Stage 1 that you must not cross until it passes.

Repo shape (all paths are country-scoped by the `COUNTRY` switch in `00_run_all.R`):
- `scripts/00_run_all.R` — master control panel; the ONLY place parameters are set. `COUNTRY <- "Uganda" | "Indonesia"`, `wd_data = data/<COUNTRY>/`, `wd_outp = output/<COUNTRY>/`.
- `scripts/01_prepare_inputs.R` → `data-raw/temp_baseline_rates_gbd.rds` (GBD 2023, RHD + All-causes, Deaths/Prevalence/Incidence, Number + Rate, by GBD age group × sex × year).
- `scripts/02_build_demography.R` → `pop_observed_1990_2024.rds`, `pop_projection_2025_2100.rds`.
- `scripts/03_build_disease_model.R` → `disease_model_inputs.rds` (`rates_by_year{ir_rhd, mort_rhd, mort_all, oth_mort, prev_seed}`, `transitions`, `p_rhd_death=c(A,B,C,D)`, `effects`, `surgery`, `stage_split`, `coverage`, `meta{…, rhd_d_fraction}`).
- `scripts/04_calibration_random_tp.R` → `calibrated_rhd_parameters.rds` (`$tp` with IR/CF/BG.mx, `$factors`, `$diagnostics`, `$stage_calibration`, `$meta`).
- `scripts/05_build_baseline.R` → `baseline_state.rds` (assembles 02+03+04; seeds A/B/C/D = `prev_seed × pop × stage_split`).
- `scripts/06_run_prevention_model.R` → `output/<COUNTRY>/out_model/<location>.rds` = `list($wsd, $stages, $diag, $meta)`.
- `scripts/07_make_outputs.R` → `rhd_model_long.csv/.rds`, `rhd_stage_model_long.csv`, `rhd_annual_flows.csv`.
- `scripts/08_economic_evaluation.R` → `rhd_economic_results.rds` (`$summary,$budget_impact,$incremental,$params`), `rhd_economic_summary.csv`, `rhd_budget_impact.csv`.
- `analysis/report.RMD` — reads `out_model/*.rds` (`$wsd,$stages,$meta,$diag`), `baseline_state.rds`, `rhd_economic_results.rds`, and GBD/pop for validation panels.

---

## A. Root cause you are fixing (verify this in the code before you touch anything)

There are **two disconnected representations of RHD mortality** in this pipeline, and they never meet:

1. `04` calibrates a **well-sick-dead proxy**: it fits a scalar case-fatality `CF` (sick→dead) and `IR` (well→sick) by pure random search (400 i.i.d. uniform multipliers on `[0.5,1.5]`), against a **weighted relative squared error** whose deaths term uses **GBD all-cause deaths**, not RHD deaths.
2. `06`, the **actual A/B/C/D production engine**, produces every RHD death from block-G `p_rhd_death[A..D] × sap_mult` (× surgery for D). It **never reads `CF`**. `05` even documents that calibrated `CF` "is carried as a total-RHD-mortality ANCHOR only" — i.e. it is a no-op in the engine.

Consequences, which reproduce the observed Uganda symptoms:
- The only thing reconciling `06`'s deaths to GBD is the hand-tuned scalar `rhd_mortality_calibration_mult <- switch(COUNTRY, Indonesia = 1/2.5, Uganda = 1/8)` in block G of `00_run_all.R`, applied via `scale_probability(p, mult) = 1 - (1-p)^mult` to all stage death probabilities. For Uganda this drives stage-C death from 1.0%→~0.126%/yr and stage-D from 7.5%→~0.97%/yr — **clinically implausible C/D→death** — and A/B death probs are hard-zeroed.
- Because deaths in `06` are `Σ_stage stock × p_death`, **total RHD deaths are driven by the stage stock distribution** (how much sits in C/D) as much as by the per-stage probabilities. If total prevalence is inflated or the C/D share is wrong, deaths inflate regardless of the per-stage values — so scaling all mortality down is the wrong lever and merely hides a stock/progression problem.
- `04`'s deaths term is all-cause, dominated by background mortality, so it barely constrains RHD mortality at all — `CF` is weakly identified and, in any case, unused downstream.
- With ~52 independent age-group multipliers fit by 400 random draws, calibrated IR is ragged across age → **implausible age-prevalence patterns**.

**The fix is Part 9 of the objective: calibrate the production A/B/C/D model directly against GBD RHD-specific deaths, prevalence, and incidence, replacing the well-sick-dead proxy and the `1/8` scalar with a low-dimensional, prior-anchored, hazard-based calibration.**

---

## B. Staged plan with a HARD GATE (do not reorder, do not skip the gate)

**Stage 1 — Improve Layer 1 (aggregate) without touching the engine structure yet.** New GBD targets (RHD-specific), log-scale loss, RHD-deaths mortality target, hazard competing risks, dimensionality reduction, DE + L-BFGS-B optimizer, log-scale multipliers, priors + smoothness. Validate against the acceptance criteria in Part H that apply to Layer 1.

> **GATE:** Do not begin Stage 2 until Stage 1 passes its validation checks and you have run the full `01→08` pipeline for Uganda **and** Indonesia and confirmed no downstream script breaks. If Stage 1 fails validation, stop and report; do not silently proceed.

**Stage 2 — Structural A/B/C/D calibration.** Put the actual `06` engine (or an extracted, side-effect-free copy of its one-cycle update) inside the calibration objective so the stage mortalities and transitions are fit to reproduce GBD RHD deaths and prevalence with plausible stage-mortality ordering. Remove the `1/8` scalar as the default. Add partial-identification handling (save near-optimal sets) when no local echo targets exist.

**Stage 3 — Diagnostics, report, regression.** Wire the new diagnostics into the bundle and `report.RMD`; run the before/after comparison; confirm the source of the former excess deaths.

Commit at the end of each stage with a message summarizing what changed and the validation result.

---

## C. Non-negotiable guardrails

1. Inspect the persisted object schemas (`readRDS(...)` then `str()`) before editing any consumer. Preserve the `01→08` contract and every existing output column: `$wsd` = `scenario, age, cause, sex, year, well, sick, newcases, dead, pop, all.mx, eff_ir, eff_cf, location, intervention`; `$stages` stocks+flows+surgery trace as consumed by `07`/`08`; `07`'s `OUT_COLS/STAGE_COLS/FLOW_COLS`; `08`'s reads of `$stages`/`$wsd`. Add new columns/objects; do not rename or drop existing ones.
2. Preserve `COUNTRY` switching. **All new settings live in `00_run_all.R`** behind `switch(COUNTRY, …)` or plain globals consumed via the existing `getp()` pattern. **No country-specific constants hard-coded in `03`–`08`.**
3. **No double-calibration.** A calibrated multiplier is baked in exactly once (into IR / stage mortality) and consumed AS-IS downstream. `eff_ir`/`eff_cf` in the engine remain intervention multipliers only.
4. Reproducibility: explicit seeds for every stochastic step (DE, multi-start). Report whether multi-start gives materially different optima.
5. Do NOT fabricate stage-specific prevalence. If no local echo target file exists, treat stage transitions as partially identified: keep clinical priors/bounds and emit a prominent partial-identification warning.
6. **The `1/8` (and `1/2.5`) scalar must not be the mechanism that reconciles deaths.** An optional *residual* mortality multiplier may remain, bounded to a modest configurable range (e.g. `[0.5, 2]`), but **a solution that pushes this residual to its bound is a VALIDATION FAILURE, not an accepted fit** — that is the guard that prevents `1/8` reappearing under a new name.

---

## D. Part 1 — GBD calibration targets (Stage 1)

In `01`/`04`, retain GBD by location × sex × year × GBD age group × measure, at full age–sex–year resolution including paediatric groups. Build RHD-specific targets as both **rate and number**: (1) RHD incidence, (2) RHD prevalence, (3) RHD deaths. Keep all-cause deaths **for validation only**. Do not use all-cause deaths in the calibration objective. Add safeguards for zero / near-zero target cells (ε floor; optionally drop cells with GBD number below a configurable threshold from the loss but keep them for validation).

If GBD uncertainty intervals are available in the prepared inputs, implement an **optional inverse-variance weighting** on the log-scale residuals and make it the default when intervals exist; otherwise fall back to unit weights. (This is the principled route to non-arbitrary weights — prefer it over hand-set `w_*` where data allow.)

## E. Part 2–3 — Anchored, low-dimensional incidence (Stage 1)

GBD incidence defines the baseline No-RHD→A inflow. Convert the total-population GBD rate to an at-risk probability/hazard: `p_{NoRHD→A} ≈ r_GBD / (1 − P_RHD)` with numerical safeguards; prefer a hazard parameterization. Add `INCIDENCE_CALIBRATION_MODE <- switch(COUNTRY, … "anchored")` in `00`, values `"fixed" | "anchored" | "free"`, default `"anchored"`:
`IR_model = IR_GBD × exp(f_s(a))`, `f_s(a)` a **low-dimensional smooth age adjustment centred at zero** (broad-band or spline coefficients), NOT an independent multiplier per GBD age group.

Reduce dimensionality: calibrate on **configurable broad age bands** (default `0–14, 15–24, 25–44, 45–64, 65+`) while keeping high-resolution targets. Estimate at most one incidence correction and one aggregate RHD-mortality correction per band × sex (or smooth spline coefficients if simpler/stabler). Expose the bands in `00`.

> **Guard from review:** the historical Uganda bug was GBD incidence running ~10× what the prevalent pool can sustain, because GBD RHD incidence includes transient/non-chronic disease. Do **not** set the anchor prior so tight that this re-enters. Let the incidence adjustment breathe (generous `σ_α`), and add a **mass-balance diagnostic**: prevalence ≈ incidence × mean dwell time, reported per age band. This is the relation that broke before; make it first-class, not implicit.

## F. Part 4–8 — Objective, priors, sequencing, hazards (Stage 1)

**Objective (log-scale, normalized).** For incidence, prevalence, RHD deaths:
`L_X = (1/N_X) Σ_j [log(X_j^model + ε) − log(X_j^GBD + ε)]²`, `L_data = w_I L_I + w_P L_P + w_D L_D`, weights configurable in `00`. All-cause deaths excluded from `L_data`, retained as validation: `D_all_model = D_RHD_model + D_other_model`.
*(Optional, preferred where UIs exist: replace the squared-log terms with a Poisson/negative-binomial count log-likelihood so the weighting follows from the variance structure rather than hand-set `w_*`.)*

**Priors + smoothness.** With `α_g` = incidence multiplier, `β_g` = aggregate RHD-mortality multiplier, both optimized on the log scale (`α_g = exp(η_g)`, `β_g = exp(ζ_g)`):
`L_prior = Σ_g (log α_g / σ_α)² + Σ_g (log β_g / σ_β)²`,
`L_smooth = Σ_g (Δ² log α_g)² + Σ_g (Δ² log β_g)²`,
`L_Layer1 = L_data + λ_prior L_prior + λ_smooth L_smooth`. Expose `σ_α, σ_β, λ_prior, λ_smooth` in `00`.

**Sequential then joint.** Step 1: hold incidence fixed/anchored, calibrate aggregate RHD mortality to GBD RHD deaths (never all-cause). Step 2: hold mortality, estimate small incidence corrections vs GBD incidence + prevalence. Step 3: joint refinement from the sequential solution under the full penalized objective. Save separate diagnostics for mortality-only, incidence-only, and joint.

**Hazards for competing risks (replaces the additive-probability + `pmax(...,0)` clipping currently in `06`).** Convert every annual probability to a hazard `h = −log(1−p)`; for people with RHD combine RHD and other-cause hazards:
`p_any = 1 − exp[−(h_RHD + h_other)]`, `p_RHD = p_any · h_RHD/(h_RHD+h_other)`, `p_other = p_any · h_other/(h_RHD+h_other)`, with explicit zero-total-hazard handling. This guarantees `p_RHD + p_other ≤ 1` without clipping. Apply the same hazard construction to the multi-way stage outflows in Stage 2 (progression + regression + death competing within a cycle) so stage `*_stay` is a proper residual, not a floored subtraction.

**Optimizer (replaces pure random search as default).** `DEoptim::DEoptim(...)` (bounded global) → `nlminb(...)` or `optim(..., method="L-BFGS-B")` (bounded local polish), on log-scale parameters, several reproducible seeds/starts, report multi-start spread. If `DEoptim` is an unacceptable dependency, provide a configurable fallback (e.g. `optim` multi-start) but **do not keep 400 i.i.d. random vectors as the production default**. Keep the "candidate/​start 0 = baseline (all multipliers = 1)" guarantee so the fit is never worse than baseline on the penalized objective.

## G. Part 9–12 — Structural A/B/C/D calibration (Stage 2)

Run the **actual production engine** (extract `06`'s single-cycle update into a pure function both `06` and `04`/a new `04b` can call, so there is one source of truth) inside the structural objective. Structure: `NoRHD → A ↔ B ↔ C → D → RHD death`, optional `D→C` if configured, background mortality on every living stage. For each age × sex × year, aggregated to GBD target age groups, enforce in the loss:
`A+B+C+D ≈ P_RHD^GBD` and `A·m_A + B·m_B + C·m_C + D·m_D ≈ D_RHD^GBD` (on counts).

Estimate only a limited, prior-bounded set: `p_A_to_no_rhd, p_A_to_B, p_B_to_A, p_B_to_C, p_C_to_B, p_C_to_D, optional p_D_to_C, rhd_d_fraction`, plus a **parsimonious stage-mortality parameterization**: one overall level `m_D` with severity ratios `m_C = ρ_C m_D`, `m_B = ρ_B m_D`, `m_A = ρ_A m_D`, subject to the clinical ordering `0 ≤ m_A ≤ m_B < m_C < m_D` and configurable bounds/priors on `m_D` and the ratios. **Remove the `1/8` scalar as the default mechanism** (keep only the bounded residual multiplier of Part C.6, whose bound-hit is a validation failure).

**Excess-death resolution rule (Part 10):** if the model cannot hit GBD RHD deaths with plausible stage mortality, the optimizer must adjust severe-stage prevalence, C/D progression, regression, or the initial D share — **not** suppress stage mortality. Encode this as: stage mortalities carry tight clinical priors; the stock/progression parameters carry looser priors; so the penalized objective preferentially moves the stock, not the mortality.

**Partial identification (Part 11):** keep the stage-target template. If a valid local echo target file exists, add `L_stage = Σ_j w_j [log(P̂_stage,j + ε) − log(P_target,j + ε)]²`. If not: do not fabricate; use priors/bounds; save **all near-optimal parameter sets** within `L(θ) ≤ L_min + δ` (δ configurable) for later uncertainty analysis, rather than presenting one matrix as identified. In the report and bundle, state explicitly **which parameters are data-identified vs prior-dominated** (e.g. compare fitted values to priors, or report objective flatness along each parameter direction).

**Surgery (Part 12):** unchanged in kind — a service flow conditional on prevalent C/D stocks, reducing C→D and D→RHD-death, never a state, never removing stock. Keep requirement/coverage/cost/effect independently switchable; do not scale surgery coverage inside the SAP/screening scenario unless a separate surgery scale-up scenario is explicitly enabled.

## H. Part 13–14 — Diagnostics + acceptance criteria (Stages 1–3)

Emit a diagnostic table by location × sex × age group × year with: GBD vs modelled RHD incidence, prevalence, RHD deaths, all-cause deaths; A/B/C/D counts and shares of prevalent RHD; every flow (`new_rhd_A, A_to_B, B_to_A, B_to_C, C_to_B, C_to_D, D_to_C`); RHD deaths by stage; other-cause deaths by state; surgery requirement and delivery for C/D. Report implied total case fatality `CF_implied = D_RHD_model / (A+B+C+D)` and compare it to the Layer-1 mortality anchor. Add an **excess-death attribution diagnostic** that classifies whether excess deaths trace to: excessive total prevalence, excessive C/D share, excessive progression into C/D, excessive stage mortality, or **duplicate mortality application** — and add a validation check that RHD is not killed twice in a cycle (stage RHD deaths must sum exactly to `$wsd` `dead`; living states + deaths reconcile with population accounting).

Acceptance criteria (the refactor is complete only if ALL hold): objective uses GBD RHD-specific deaths; all-cause deaths validation-only; GBD incidence retained as baseline age–sex pattern; incidence adjustments low-dimensional and anchored to 1; final penalized objective ≤ baseline; all probabilities and competing-risk sums valid; state counts finite and non-negative; living + deaths reconcile with population; stage RHD deaths sum exactly to total model RHD deaths; A/B/C/D sums exactly to total model RHD prevalence; mortality not applied twice; **Uganda no longer needs a default `1/8` multiplier and the bounded residual multiplier does not sit at its bound**; fitted Uganda stage-D mortality within a documented plausible range; if deaths remain excessive with plausible mortality, diagnostics name the responsible stock/progression parameter; Indonesia and Uganda separately configurable; `07`/`08`/`report.RMD` still run; surgery remains a service flow.

## I. Part 15 — Output bundle

Extend `calibrated_rhd_parameters.rds` to at least:
`$layer1{incidence_parameters, mortality_parameters, objective_components, optimizer_diagnostics, validation}`; `$stage_calibration{status, best_parameters, accepted_parameter_sets, objective_components, targets_file, priors, bounds}`; `$diagnostics{fit_by_age_sex_year, stage_stocks, stage_flows, deaths_by_stage, all_cause_validation}`; `$meta`. Keep writing human-readable CSV diagnostics. Do not remove `$tp`/`$factors` if any consumer still needs them — add alongside and migrate consumers explicitly.

## J. `analysis/report.RMD` updates (Stage 3)

- Replace any methods text describing the `1/2.5` / `1/8` mortality scalar and the "CF carried as an anchor" language with the new calibration description (RHD-specific deaths target, anchored incidence, hazard competing risks, severity-ratio stage mortality, priors/smoothness, DE+L-BFGS-B, partial identification).
- Add the exact mathematical objective actually implemented (Parts F/G) to the methods section.
- Add figures/tables from the new diagnostics: GBD-vs-model fit for incidence/prevalence/RHD deaths by age; A/B/C/D shares; `CF_implied` vs anchor; the excess-death attribution; and, where echo targets are absent, the near-optimal parameter-set spread with a partial-identification caveat.
- Keep the model schematic and existing narrative accurate to the refactored structure. Do not break existing chunks that read `out_model/*.rds`, `baseline_state.rds`, or `rhd_economic_results.rds`.

## K. Part 16 — Deliverables (report back at the end)

1. Summary of every modified file and why. 2. The new calibration sequence in prose. 3. Every new `00_run_all.R` setting, with default and meaning. 4. The exact implemented objective. 5. Assumptions forced by missing local A/B/C/D targets. 6. Full pipeline run for Uganda. 7. Before/after Uganda comparison: RHD prevalence, RHD deaths, C and D prevalence, implied total CF, stage-C mortality, stage-D mortality, deaths by stage, intervention deaths averted. 8. A definitive statement of whether the former excess-death problem arose from mortality probabilities, stage distribution, progression, initialization, or duplicate mortality accounting — backed by the attribution diagnostic. 9. A regression run for Indonesia confirming the second country still works.

Work incrementally, validate at each gate, and stop-and-report rather than proceeding on a failed check.
