# Claude Code prompt: implement the default RHD natural history and improve childhood-death calibration

Work in the existing `work-uw-rhd-modelling` repository. Implement the changes below, rather than only proposing a plan.

## Goal

Produce a simplified but plausible RHD natural-history model and a credible, calibrated baseline. Improve reproduction of GBD RHD-attributable deaths **below age 20**, which the current model underestimates, while retaining a reasonable prevalence age pattern and adult mortality fit. Use the default implementation described below, informed by the reference workbook in `dev/`.

The workbook is a conceptual reference, **not a production input interface**. All parameters needed for this task are specified here. Do not reference spreadsheet cell addresses, read parameter values from workbook cells, or add an Excel dependency to the production pipeline. Values described as provisional are assumptions, not validated clinical estimates.

## Scope and constraints

- Edit only necessary lines within `scripts/00*.R`–`scripts/08*.R`, including `04b`, and `R/abcd_engine.R`. This is an allowlist, not a requirement to edit every file. Preserve `01` and `02` unless an actual input-retention or demographic-interface defect requires a small change.
- Preserve existing functions, naming, data structures, dependencies, parallel execution, file locations and workflow as far as practical. Extend existing structures where additional history is essential; avoid a wholesale refactor.
- Put **all user-editable assumptions, input values, calibration settings, bounds, weights and age rules in `00_run_all.R`**. Downstream scripts consume the resolved configuration. Avoid competing silent defaults that override it.
- Implement one default specification. Do not add sensitivity scenarios, alternative structures, parameter selectors or comparison grids. Keep the existing reference-care and intervention-scale-up arms: they evaluate the same natural history under different care.
- Do not create unrelated folders, scripts, documentation systems or dependencies. Use existing diagnostic/output locations. Small helpers needed for shared dynamics belong in `R/abcd_engine.R`.
- Do not edit `analysis/report.Rmd`, `reports/slideck-technical.Rmd`, other report/deck files, or their supporting assets. Read them only to understand output contracts and existing assertions. Their methods text will be updated separately.
- Preserve the configured country selector. The numeric reference natural history below comes from the illustrative Indonesia workbook. Do not silently turn Indonesia-specific inputs into Uganda observations or overwrite Uganda's country-specific mortality, demographic, cost and care inputs.

## Inspect the current implementation first

Trace `00 → 03 → 04 → 04b → 05 → 06 → 07 → 08` and the shared engine. Before editing, identify where each required change enters the calculation and its affected callers and checks.

The supplied version has several relevant features:

1. Incident disease enters A only. With D-only mortality and one event per annual cycle, reaching D can take several cycles; this can suppress childhood deaths.
2. `STRUCT_DEATH_TARGET_AGE_MIN <- 20L` excludes childhood deaths from the Layer-2 objective. Removing this exclusion is necessary but does not by itself fix the natural history.
3. Current adult modifiers depend on attained age. They do not preserve early- versus later-onset history.
4. A, B and C RHD mortality are already exact zeros in the active D-only branch; retain this property.
5. Calibration and production already share `abcd_one_cycle`, `abcd_project`, `abcd_burn_in` and ageing helpers. Preserve that single source of truth.
6. Current screening is an annualized two-window encounter schedule. The requested workbook default instead uses one screen during ages 5–39 and a pooled diagnosis/BPG cascade. This is a deliberate change: update its connected quantities coherently rather than combining incompatible interpretations.

Do not claim that childhood underestimation has been quantified from these source files alone. Establish the actual discrepancy from the repository's GBD targets and saved or regenerated results.

## 1. Explicit default natural-history inputs in `00`

Use annual probability inputs. Convert them to cause-specific hazards once with `h = -log1p(-p)` before combining competing events. These are not realized competing-risk transition proportions and must not be annualized again.

### Pathways and stage-entry composition

Maintain two onset-history pathways for life:

- **Y:** first modeled onset at age **≤20**.
- **O:** first modeled onset at age **≥21**.

The use of age 20 in Y follows the workbook's modeling convention. It is separate from the outcome “deaths below 20,” which includes ages 0–19 only. Do not relabel this cutoff as a verified WHF diagnostic rule.

| Parameter | Default | Treatment in this implementation |
|---|---:|---|
| Last age for Y onset and childhood regression rules | 20 years | Fixed structural boundary, editable in `00` |
| First age for O onset and adult normalization | 21 years | Derived from boundary |
| Y incident entry mix, A/B/C/D | 0.50 / 0.30 / 0.15 / 0.05 | Fixed provisional mix |
| O incident entry mix, A/B/C/D | 0 / 0.60 / 0.30 / 0.10 | Fixed provisional mix |
| O:Y multiplier on B→C and C→D hazards | 0.48 | Fixed provisional multiplier; applies by onset history |
| Adult B→No RHD normalization input | 0.02 per year | Fixed for both Y and O at attained ages >20 |
| D→C | 0 | Exact structural zero |
| RHD mortality in A, B and C | 0 | Exact structural zeros |
| Residual RHD mortality multiplier | 1 | Fixed; no second mortality scale |

Use these Y annual transition values as starting values and prior centers:

| Transition | Y at attained age ≤20 | Y at attained age >20 | O |
|---|---:|---:|---:|
| A→No RHD | 0.311 | 0 | No A state |
| A→B | 0.054 | No A stock remains | No A state |
| B→A | 0.177 | 0 | 0 |
| B→C | 0.074 | Same Y parameter | Hazard-scaled Y value, multiplier 0.48 |
| C→B | 0.050 | 0 | 0 |
| C→D | 0.112 | Same Y parameter | Hazard-scaled Y value, multiplier 0.48 |
| D→C | 0 | 0 | 0 |
| B→No RHD | 0 | 0.020 | 0.020 |

For the O progression probabilities derive `p_O = 1 - (1 - p_Y)^0.48`. At the starting values, O B→C is approximately **0.0362302884**, and O C→D approximately **0.0554211359**. Derive these in code; do not hard-code rounded duplicates. Early-onset Y adults retain Y progression rates, rather than switching to O rates on their 21st birthday.

The childhood regression values, the entry mixes and the O multiplier are provisional inputs. In particular, 0.177 B→A and 0.050 C→B must not be presented as well-established estimates for clinically detected disease. Their calibration status and provenance must remain transparent.

### Retain high D mortality

For Indonesia fix the Y stage-D annual RHD-death input at **0.119** and the O:Y D-death hazard ratio at **2.0**, giving:

```r
p_D_O <- 1 - (1 - 0.119)^2  # 0.223839
```

These are provisional reference-model inputs, not fitted country estimates. Retain them as fixed in this default so childhood fit is not obtained by reducing D fatality. Competing other-cause mortality and the existing reference surgery service can make realized annual D fatality smaller than the input.

For Uganda preserve the existing country-specific Y D-death input, **0.0597**, rather than substituting 0.119. If running the shared default in Uganda, explicitly label use of the same O:Y hazard ratio and other shared workbook assumptions as extrapolations; derive O mortality from the selected country's Y value. Do not invent country-specific evidence.

Remove D mortality and the residual multiplier from the free parameter vector for this active default. Do not retain unused fitted mortality parameters or a fitted seed D fraction that cannot affect the forward-generated initial state.

### Parameter estimation versus fixed inputs

Keep the existing parsimonious incidence corrections and allow the six nonzero Y transition parameters to calibrate around the centers above, shared across sexes unless an existing evidence-supported sex distinction is already present. Keep their existing annual-input bounds:

```r
TRANS_BOUNDS <- list(
  p_A_to_no_rhd = c(0.005, 0.60),
  p_A_to_B     = c(0.001, 0.25),
  p_B_to_A     = c(0.002, 0.40),
  p_B_to_C     = c(0.001, 0.20),
  p_C_to_B     = c(0.001, 0.15),
  p_C_to_D     = c(0.001, 0.30),
  p_D_to_C     = c(0, 0)
)
```

Fixed zeros do not enter log/logit transforms or the optimizer. Fix entry mixes, adult normalization, the O progression multiplier and both D mortality inputs at their default values. Do not also fit the old attained-age adult progression/remission multipliers: they would conflict with the new pathway rules. Use identifiable parameters only, and retain uncertainty/identifiability diagnostics within the one fitted specification.

Retain the existing structural regularization settings unless an actual numerical defect requires a documented correction: `SIGMA_TRANS = 2.0`, `LAMBDA_PRIOR_STRUCT = 0.05`, `SIGMA_ALPHA_STRUCT = 1.0`, `LAMBDA_SMOOTH_STRUCT = 0.5`, incidence multiplier bounds **0.2–5.0**, and the existing reproducible optimizer settings. Report fitted versus fixed inputs and boundary hits.

## 2. Minimal history-aware states and coherent annual dynamics

Represent Y-A/B/C/D and O-B/C/D internally, retaining Y/O history when people age. Also retain separate normalized/recovered No-RHD stocks for Y and O. Aggregate these internal disease stocks back to the existing A/B/C/D interface for downstream consumers.

- Only never-onset people are eligible for first modeled onset. Normalized/recovered people remain alive and face background mortality but cannot have another onset in this simplified default. Explicitly record this no-recurrence assumption; normalization does not prove complete anatomical healing.
- A→No RHD in Y children and B→No RHD in Y/O adults leave prevalence, SAP eligibility and SAP costs. Do not return these people silently to the never-onset susceptible pool.
- At the age-20→21 boundary, transfer surviving Y-A to Y-B, including end-of-cycle incident Y-A at age 20. Record the transfer separately. It conserves living RHD stock and is neither new onset nor observed biological A→B progression.
- This A→B transfer is a **provisional modeling assignment**, not a WHF requirement or proof that the echo satisfies adult stage-B criteria. Do not use the adult reclassification proportion from longitudinal follow-up to set biological normalization.
- Adult Y-A must be exactly zero. Adult B→A and C→B are zero; adult B normalization remains active.
- Preserve one event per annual cycle. Start stocks experience competing events; new incident cases enter the selected stage mix at year end; then survivors age. A case entering D this year first faces D-specific death in the following cycle. Do not quietly add same-cycle onset-and-death or a separate acute-carditis state.
- Competing other-cause mortality applies to every living state. Retain real age/sex/country background mortality and the existing WPP demographic backbone. Do not substitute the workbook's illustrative mortality table.
- Reconcile the additional No-RHD history stocks with the exogenous population backbone using the existing demographic convention. Explicitly account for births or residual population inflows where needed, and reject impossible negative pools or stock-over-population states. Do not rescale disease stocks to hide a deficit or manufacture deaths.

Use exactly the same state transitions, ageing, observation convention, care and population reconciliation in burn-in, calibration, bridge and production. `light = TRUE` and the full engine must agree numerically on all shared outputs.

Do not reproduce the workbook's two closed 100,000-person examples, 1% initial prevalence, illustrative onset rates, 50% adult Y share, or fixed 30 DALYs per death as national inputs. Generate the national pathway/stage distribution from incidence, history and forward dynamics. Do not impose a 50% Y share on all adults.

## 3. Default screening, pooled care and intervention effects

Implement the reference workbook's simple pooled care interpretation without adding a full stage-by-detection-by-treatment model. All modeled disease stages A–D are eligible in this illustrative default; recognize this as a program assumption rather than a stage-specific clinical recommendation. Normalized/recovered people are excluded from SAP.

For Indonesia centralize these inputs in `00`:

| Input | Reference baseline → target | Scale-up baseline → target |
|---|---:|---:|
| Target probability of one screen across ages 5–39 inclusive | 0.05 → 0.05 | 0.05 → 0.80 |
| Annual clinical detection among undiagnosed | 0.02 → 0.02 | 0.02 → 0.02 |
| BPG initiation among diagnosed untreated | 0.50 → 0.50 | 0.50 → 0.90 |
| Annual BPG retention | 0.90 → 0.90 | 0.90 → 0.95 |
| Good-adherence share among those on BPG | 0.60 → 0.60 | 0.60 → 0.85 |
| Partial-adherence share among those on BPG | 0.20 → 0.20 | 0.20 → 0.10 |

Additional defaults: screening sensitivity **0.90**; relative protection under partial adherence **0.50**; first modeled observed-history diagnosed share **0.05** and on-BPG share **0.04**; initial ever-screened share **0** as a provisional initialization. Initialize these once at the beginning of the historical simulation, carry their history forward and do not reset them to these values in 2026. Use reference care throughout historical initialization/calibration and the bridge. Use the configured scale-up ramp, currently **2026–2050**, with targets held after 2050. Preserve other configured time horizons.

For a constant lifetime-window target `c`, the annual probability among never-screened eligible people is `q = 1 - (1 - c)^(1/35)`. In the ramp use the current year's target to derive q, document that cumulative realized cohort reach during a changing ramp need not equal the current target, and track the never-screened share by cohort. Encounters are living eligible people × never-screened share × q. Screening continues to include people without disease; it is not RHD prevalence × coverage.

For the pooled annual care update, let `d` and `b` be the start-of-cycle diagnosed and on-BPG shares among living RHD, `u` the never-screened share, `s` sensitivity, `k` clinical detection, `r` retention, and `i` initiation:

```text
d_after = d + (1-d) * [1 - (1-k)*(1-q*u*s)]
b_retained = b*r
b_after = min(d_after, b_retained + (d_after-b_retained)*i)
adherence_factor = good_share + 0.5*partial_share
effective_protection = b_after*adherence_factor
```

Use q=0 outside screening ages; clinical detection and ongoing SAP can continue outside those ages. Age and carry pooled care consistently, dilute diagnosed/treated shares when undiagnosed incident cases enter, and remove disease exits proportionally under this explicitly pooled approximation. Do not count the same diagnosis twice or turn normalized cases into treated RHD. Reject invalid shares; do not mask a broken state update by arbitrary clipping.

Use **on-BPG person-years**, not adherence-weighted protection, for SAP costs. Apply adherence once to health effects. Replace the current `min(screen, diagnosis, treatment)` coverage shortcut for this default; annual screening probability must not cap continuing SAP coverage. For other countries retain existing country-specific baseline care information where available and label any shared new cascade values as assumptions.

SAP acts through progression hazards only:

- A→B HR = `log(1-0.003)/log(1-0.075)` = approximately **0.0385383481**.
- B→C HR = `log(1-0.025)/log(1-0.119)` = approximately **0.1998285475**.
- C→D HR = **0.45**, provisional assumption.
- D→RHD death HR = **1**: no direct SAP mortality benefit.

For effective protection `c_eff`, apply `h_protected = h_untreated * [1 - c_eff*(1-HR)]`. The GOAL-related inputs are transcribed reference inputs; do not claim that extrapolation to adult pathways is a directly estimated adult effect. The adult B normalization input **0.02** is an extrapolation from limited longitudinal evidence, not a directly observed annual adult stage-B estimate. Do not recalibrate it from an age-based “classified normal” proportion.

Keep the existing production surgery service and its country-specific inputs, background mortality, economic valuation and discounting. Surgery must remain identical in the reference and scale-up arms unless the existing authorized specification says otherwise. Apply it once, in both calibration and production; do not add a new surgery state or remove surgery simply because it is absent from the illustrative workbook.

Retain existing country-specific unit costs. Indonesia defaults remain **US$110 per on-BPG person-year** and **US$1.10 per screening encounter**. Update `08` only where corrected service volumes require it. Do not replace existing age-specific health-loss calculations with the workbook's illustrative 30 DALYs per death.

## 4. Calibration: include childhood deaths without hiding adult misfit

Preserve the **existing calibration window 2000–2019**. Also preserve `STRUCT_MATCH_FROM = 2010`: in the supplied version Layer 2 scores 2010–2019, while earlier window years remain initialization/validation years. Do not silently extend or move either window. Preserve observed years **1990–2024**, the GBD rate reference year **2023**, analysis years **2026–2100** and the 2026–2050 scale-up ramp unless the repository configuration has explicitly changed them.

Set `STRUCT_DEATH_TARGET_AGE_MIN = 0L` in `00` and remove the active exclusion of deaths below age 20. Keep missing observations missing. Investigate the existing count-floor mask; retain `TARGET_MIN_COUNT = 1` for individual-cell loss stability, but report every available childhood target, including small cells, and show age-aggregated childhood fit so the floor does not conceal the problem.

Retain prevalence-led joint calibration on the actual production engine:

- Prevalence weight **2.5**, total death weight **1.0**, direct incidence weight **0**; incidence remains an anchored input and validation output, with joint incidence corrections informed by prevalence/deaths.
- Retain the existing broad incidence bands and country/sex structure initially. Do not add numerous age-specific free parameters to fit every GBD cell.
- Within the death term use a configurable childhood share **0.50**: `L_death = 0.5*L_death_age_0_19 + 0.5*L_death_age_20plus`. Each component is a normalized mean log-count error over its usable target cells. This is an operational weighting decision for this task, not a published statistical standard. Export both components and the counts/weights entering them; stop clearly if a required component has no usable data rather than replacing it with zero loss.
- Keep `STRUCT_EPS_COUNT = 1`, `STRUCT_PREV_OBS = "mid"`, `USE_IV_WEIGHTS_STRUCT = FALSE`, and existing usable stage-target handling.
- Retain `STRUCT_TARGET_BASIS = "rate_x_pop"`: construct targets from GBD Rate × the model's WPP population for the same age/sex/year cell. Also report published GBD Numbers and their difference from the population-harmonized target. Do not imply that matching one automatically matches the other.
- Use fixed entry mixes to permit existing advanced disease at modeled onset, instead of forcing every case to spend years passing through A. Interpret direct entry into C/D as a coarse annual representation of advanced disease at modeled entry, not proof of a specific biological bypass mechanism.
- Preserve the existing zero-start historical burn-in approach and its **120-cycle** default, updated for all new history/care states. Check convergence of pathway/stage age profiles and childhood deaths, not just total prevalence. Evaluate a doubled burn-in at the best solution to diagnose initialization dependence; do not repeat expensive searches unnecessarily. Do not force a no-screening-history equilibrium check on an absorbing screened-history stock.
- Carry the fitted state through bridge years **2020–2025** under reference care into 2026. Preserve all pathway, recovered, screening and care histories. Do not reset disease prevalence or stage mix to GBD in 2026.

The first calibration layer is a proxy, not the final natural history. Make only the edits needed to its target construction, metadata or starting values; the final baseline, incidence conversion, mortality and child-death claims must come from Layer 2 on the updated engine. Distinguish incidence per total population from first-onset hazard per never-onset susceptible population. Do not retain a conversion based only on `1-prevalence` when recovered people are excluded from susceptibility, or silently let a proxy case-fatality estimate overwrite fixed D-only mortality.

Run the updated calibration when the repository data and runtime are available. This request authorizes necessary changes to `04` and `04b` and regeneration of affected calibration outputs. A saved calibration from the old structure is not a valid final fit. Preserve the existing last-valid-output protection and add the minimal specification check needed to prevent stale bundles from being consumed.

Do not force a successful status. If fixed D fatality, entry composition and the simplified pathway structure cannot jointly reproduce GBD childhood deaths and adult prevalence/deaths, retain the discrepancy and report the limiting assumption. Do not add an acute pathway, lower D mortality, change the time window, suppress difficult ages or widen bounds silently to make the fit pass.

## 5. Preserve downstream output contracts truthfully

Maintain existing filenames and top-level RDS components, including calibration `$stage_calibration`, `best_parameters`, baseline state, per-location `$wsd`, `$stages`, metadata, long outputs and economic summaries. Existing scenario/location/age/sex/year keys must remain stable. Add history detail and diagnostics within the existing structures where practical.

In particular:

- Aggregate Y/O disease stocks to A/B/C/D. Aggregate `sick` consistently with the existing stock/observation timing. `well` includes both never-onset and recovered/normalized people; it must not be mislabeled as the susceptible stock.
- Keep `new_rhd_A` as **actual incident entries into A**. Add incident entries into B/C/D and an explicit total-onset field. Set `$wsd$newcases` to total onset across entry stages. Do not disguise total onset as A-only entry to satisfy an old assertion.
- Add adult B→No RHD flows, boundary transfers, recovered stocks and total disease exits. Correct stock-flow checks: end RHD = start RHD + total onset − A resolution − adult B normalization − RHD deaths − other deaths among RHD. Boundary transfers conserve the total.
- Retain stage-specific death fields and ensure RHD deaths outside D are exactly zero. Provide pathway-specific D mortality detail without falsely presenting a single pooled D input as the two-pathway mechanism.
- Retain `n_on_optimal_treatment` as delivered on-BPG person-years and `n_screened` as actual screening encounters. Cost each service once. Expose effective protection separately from treatment volume.
- Keep natural-history and reference-care metadata consistent between calibration and baseline; no stale inherited attained-age multipliers or contradictory all-stage interpretation.
- Preserve meaningful legacy fields when they can be calculated truthfully. Do not fabricate scalar parameters, annual encounter-rate vectors, calibration statuses or coverage trajectories just to make old report assertions pass.

Read the report/deck setup checks before finalizing compatibility. The new screening history, pooled care and pathway-specific mortality may require future changes to some report assertions or methods text. List these precise dependencies in the final response, but do not edit the reports or misrepresent their renderability. The numerical outputs should remain usable for that later reporting update.

## 6. Verification and diagnostics

Use focused checks in existing validation sections or temporary scratch checks; do not introduce an unrelated test framework. At minimum verify:

1. Nonnegative, finite stocks/flows; stage-wise competing-event probabilities plus staying probability sum to one; full within-cycle population conservation, with exogenous demographic adjustments accounted for separately.
2. Y/O identities survive ageing; adult Y-A is zero; the 20→21 transfer happens once and is separate from onset; normalized people do not re-enter susceptibility or SAP; no recurrent onset in recovered people.
3. D-only mortality is exact; an incident D case does not die in its entry cycle; fixed D inputs and the O:Y hazard ratio are preserved after calibration and baseline construction.
4. Annual adult B normalization is an input of 0.02 and its realized competing-risk flow is no larger than the isolated input risk, with no childhood B-normalization flow.
5. The shared engine's light and full modes agree. Replaying reference care from the calibrated state reproduces the corresponding calibration/bridge dynamics. Parallel workers receive every newly needed helper and configuration field.
6. Screening history never screens a person twice under the single-screen default; diagnosis and on-BPG shares are valid; on BPG ≤ diagnosed, and protection ≤ on BPG. Population changes, disease onset/exits and zero-RHD strata do not create impossible histories. SAP costs use treated volumes, not protected volumes.
7. Produce before/after fit diagnostics against identical target definitions and years. Report deaths at ages 0–19 and 20+, available constituent GBD childhood groups, sex, and year; show model/target counts, ratios, childhood share of all RHD deaths and residuals. Aggregate counts before taking ratios. Also show prevalence by age/sex, adult and all-age deaths, stage/pathway shares, incidence, objective decomposition and parameter boundary hits.
8. Use existing fit thresholds `STRUCT_FIT_TOL_PREV = 0.15` and `STRUCT_FIT_TOL_DEATH = 0.25` on the absolute log-ratio scale as diagnostic thresholds. Apply the death threshold separately to childhood and adult aggregate deaths as well as all ages. Report results across scored years, not just one favorable endpoint; zero targets need absolute-error diagnostics rather than undefined ratios. These are implementation tolerances, not clinical validation.
9. Reconcile `07` aggregates and `08` screening/SAP costs to engine service volumes. Preserve the existing age-specific health/economic calculations and sign conventions. Verify downstream objects and keys without editing or claiming to have updated report text.

Use comparable existing saved outputs for the before condition where available. If they are absent or stale, say so; do not invent a numerical improvement or run an additional legacy sensitivity model just to manufacture a comparison. If GBD data, packages or runtime are unavailable, finish the authorized code changes and feasible checks, then state exactly what could not be run.

## Final response

Summarize the files and substantive lines changed, the fixed versus calibrated parameters, the retained calibration/scoring windows, verification performed, and whether childhood deaths actually improved without unacceptable prevalence or adult-death deterioration. Give the exact rerun sequence and existing output paths. Identify residual epidemiological limitations, any extrapolated assumptions and report/deck compatibility work that remains. Do not label a numerically completed but poorly fitting run a credible calibrated baseline.
