# Claude Code task (v3): simplify the RHD model into entry-age pathways and build a credible calibrated baseline

**Repo:** `uw-wb-indonesia-rhd` · **Countries:** Indonesia (primary) and Uganda. Both must run through `source("scripts/00_run_all.R")` by changing only `COUNTRY`.
**Reference workbook:** `RHD_entry_age_pathway_model_Neema_updated*.xlsx` (Neema revision, 3 Nov 2026). Find its real path and do not assume it is in `dev/`. Do not edit it.

This is an **implementation task**. Inspect the repository, make the changes, run the checks, and run the revised calibration and projections when the data and runtime are available. Before editing, write a short plan into the diagnostics folder, then proceed without waiting for approval. Do not stop after the plan. Optimiser convergence alone does not make a baseline credible.

---

## 0. Scope, precedence and house rules

**Authorised source changes (only where the new specification requires them):** `scripts/00_run_all.R`, `01_prepare_inputs.R`, `02_build_demography.R`, `03_build_disease_model.R`, `04_calibration_random_tp.R`, `04b_calibrate_structural.R`, `05_build_baseline.R`, `06_run_prevention_model.R`, `07_make_outputs.R`, `08_economic_evaluation.R`, and `R/abcd_engine.R`.

- Uploaded copies may carry suffixes such as `00_run_all(10).R`. Resolve the actual repository files before editing.
- Do **not**:
  - add new source files, rename scripts, edit the workbook, or edit reports or slide decks;
  - refactor unrelated code;
  - add dependencies beyond base R and recommended packages (`splines`, `mgcv`, `stats`).
- Put validation code inside the authorised scripts or the engine.
- Read downstream report consumers read-only, to protect their input contracts.
- Generated products, diagnostics, the parameter registry and the methods note go into the existing COUNTRY-scoped `data/<COUNTRY>/` and `output/<COUNTRY>/` folders.

**Precedence:**
1. This task's decisions.
2. The workbook's primary inputs (Set 1).
3. Compatible existing code.

The workbook still has residual-state branches and nonzero later-onset regression switched on by default (`Inputs!B53`, `B113`–`B119`). This task **supersedes** them; do not reproduce them.

The workbook is illustrative and uncalibrated. Do **not** import any of these as country data:
- its two closed 100,000-person cohorts
- the 1% starting prevalence
- the background-mortality table
- the 50% adult Y share
- the incidence table

Keep the population-based country model: WPP backbone, GBD observations, both sexes, single ages 0–100, and the configurable analysis horizon.

**House style to keep:**
- `data.table` loaded before `wpp2024`
- `getp()` so each script still runs standalone, but standalone runs must load or validate the same registry and never revive legacy defaults
- fail-loud `stop()` before anything is persisted
- stable output filenames
- **no monetary values outside 08**

Give every new parameter a source tag (`[WORKBOOK Inputs!Bxx]`, `[LIT]`, `[ASSUMPTION – analyst chosen]`, `[CALIBRATED]`).

---

## 1. Model structure: two lifelong pathways defined by age at onset

### 1.1 Pathways and states

**Onset boundary:** `ONSET_AGE_MAX_EARLY <- 20L`.
- **Early onset (Y):** onset at attained age **≤ 20**.
- **Later onset (O):** onset at age **≥ 21**.

Use this inclusive rule identically in reconstruction, calibration, initialisation, projection, validation and outputs. GBD age groups that straddle the boundary are split on the single-age grid with WPP population weights.

States:
- **Disease-free:**
  - `N0`: never had onset; this is the only susceptible pool.
  - `R_Y`: recovered early-onset (Y A→No RHD).
  - Both age, die of background causes, and conserve people.
  - `R_Y` is **not susceptible**: there is no re-incidence and no first-ever O onset from recovered Y, matching the workbook's exclusion. No fitted recurrence.
- **Y (early onset):** `Y_A, Y_B, Y_C, Y_D`. Membership is lifelong. A Y case keeps Y transitions as it ages past 20, with no Y→O transfer on a birthday or at adult diagnosis.
- **O (later onset):** `O_B, O_C, O_D`.
  - There is no O-A, no residual compartment and no resolution branch; residual abnormalities sit inside O-B.
  - O cases stay put, progress B→C→D, or die.
  - **Exact structural zeros:** O B→No RHD, B→A, C→B and D→C are all 0, not small fitted values.
- **Absorbing:** RHD death and other-cause death.
- **Background mortality** competes in every living state, including `N0` and `R_Y`.

**Zero adult regression = zero regression for O only.** You selected this interpretation explicitly. Do **not** switch off Y regression when a Y case reaches adulthood.

D→C = 0 in both pathways.

### 1.2 Biological stocks vs diagnostic reporting

Keep **biological stocks** separate from **diagnostic reporting**. The engine moves people between biological stocks only. Diagnostic categories are **derived masks** applied to those stocks for reporting, screening and treatment eligibility. They never cause a flow.

**Adult Y-A is kept internally.**
- A Y case still in A at age 21 or older stays in `Y_A`. It keeps its early-onset history and its usual Y transitions: A→B progression and A→No RHD (`R_Y`) remission at the fixed Y rates.
- Ageing past 20 does **not** move adult Y-A to B or to No RHD. WHF limits stage A to ages ≤ 20, but leaving that diagnostic category is not biological recovery.
- Assert that no ageing step moves people from `Y_A` to any other stock.

**Two reported prevalence measures** (by age, sex, year; also aggregated to GBD age groups):

| Measure | Definition | Role |
|---|---|---|
| **All-stage model prevalence** | `Y_A + Y_B + Y_C + Y_D + O_B + O_C + O_D` at all ages, including adult Y-A | **Primary calibration target** against GBD RHD prevalence, through an explicit, **provisional** mapping (below) |
| **WHF-classifiable prevalence (approximation)** | All-stage prevalence **minus Y-A at age > 20**: `Y_A[age ≤ 20] + Y_B + Y_C + Y_D + O_B + O_C + O_D` | Reported. It is also the stock that screening can detect (section 5.2) |
| Adult unclassifiable A-like stock | `Y_A[age > 20]` | Reported diagnostic; equals all-stage minus WHF-classifiable exactly |

**Assumptions to document** in the methods note and the registry:
- Modelled B–D, in either pathway and at any age, are **assumed** to meet the relevant WHF 2023 diagnostic criteria, including the adult requirement of at least 2 morphological features for stage B. This is an approximation, not a verified mapping.
- The mapping of all-stage prevalence to GBD RHD prevalence is provisional. It is unknown whether GBD includes borderline or A-like disease, or adult A-like abnormalities.
- Do not hide a prevalence mismatch by dropping adult Y-A from the stock or declaring it resolved.
- Fitting GBD to WHF-classifiable prevalence instead is allowed only as a labelled mapping sensitivity (section 8.2), never the primary fit.

### 1.3 Cycle convention

Use one annual convention everywhere:

**start-of-year stocks → within-year events → end-of-year stocks → ageing → next start**

- New onset arrives at the end of the cycle. It cannot progress or die of RHD until the next cycle, but onset itself competes with background death in `N0`.
- State this convention in the engine header and validate conservation under it.
- Prevalence observation stays at the existing mid-year approximation `(start + end)/2`, used identically in reconstruction, loss and validation.

---

## 2. One natural-history registry (`rhd_natural_history` in 00 block G)

Replace the scattered `p_*` globals, `TRANS_BOUNDS`, `remission_adult_mult`, `STRUCT_CALIB_AGE_STRUCTURE`, `AGE_MULT_*`, `SIGMA_TRANS`, `rhd_d_fraction` and the `all_stage` branch with **one** named list.

**Fields for each entry:**
- `value`
- `scale` (input probability / hazard / multiplier / share)
- `pathway`
- `sens_lo`, `sens_hi` (admissible sensitivity range)
- `opt_lo`, `opt_hi` (fitted entries only)
- `status` (`fixed` / `fitted` / `derived` / `structural_zero`)
- `source`
- `workbook_cell`
- `note`

Fixed parameters must **never** enter the optimisation vector. A range with no evidence behind it is labelled "analyst chosen"; never attribute an invented bound to a paper. 03 exports the registry, with every transformation logged, as `natural_history_registry.csv`.

**Import rule:** check the workbook formulas before importing.
- The Set-1 values below are annual **input probabilities** that the workbook turns into hazards (`h = −ln(1−p)`) and combines through competing risks. They are not realised one-year transition probabilities.
- Do not annualise them a second time.
- Do not de-treat the already de-treated C→D input again.

| Parameter | Y | O | Status | Source |
|---|---|---|---|---|
| A→No RHD (→ `R_Y`) | 0.311 | — | fixed | `Inputs!B36`: WHF 50%/2y, competing-risk annualised jointly with GOAL A→B |
| A→B | 0.054 | — | fixed | `B37`: GOAL control 7.5%/2y |
| B→A | 0.177 | 0 (structural) | fixed, **provisional** | `B38`: far above the Cannon clinic-detected fit 0.0072 (`C38`). Keep this warning; 0.0072 is a labelled fixed sensitivity, never fitted |
| B→C | 0.074 | hazard × `m_O` (≈0.036 at 0.48) | fixed / derived | `B39`: GOAL control 11.9%/2y |
| C→B | 0.050 | 0 (structural) | fixed | `B40`: expert ≤5%/yr |
| C→D | 0.112 | hazard × `m_O` (≈0.055 at 0.48) | fixed / derived | `B41`: Cannon S1 refit, already de-treated |
| D→C | 0 | 0 | structural zero | `B42` |
| RHD death A/B/C | 0 | 0 | structural zero | D-only structure |
| RHD death D `p_D_Y` | start 0.119 IDN / 0.0597 UGA | `1−(1−p_D_Y)^2` | **fitted** (Y); O derived | `B57`/`B58`: provisional column-K proxy anchors, not clinical estimates |
| O/Y D-death hazard ratio | — | 2 | fixed | `B59`: Cannon S2/S1 2.1, an imperfect proxy for onset history (S2 = ages 15–24); sensitivity 1–4 |
| `m_O` (O progression hazard multiplier, shared by B→C and C→D, never applied to Y) | — | start 0.48 | **fitted** | `B45`: geometric mean of the Cannon S2/S1 HRs 0.63/0.37; workbook sensitivity 0.3–0.8 |
| Entry mix | A/B/C/D 0.50/0.30/0.15/0.05 | 0/0.60/0.30/0.10 | fixed, provisional | `B29:B32`, `C29:C32` (expert estimate, Oct 2026) |
| Onset boundary | ≤ 20 | ≥ 21 | fixed | `B52` |

Keep the Set-2 values (`Inputs!C36:C42`) as `value_set2` for labelled sensitivity runs (`NH_SET <- 1L`).

**Mortality rules (D-only):**
- RHD-attributable death comes **only** from stage D. RHD mortality in A, B and C is **exactly zero** in both pathways. Other-cause mortality stays active in every living state, including `N0` and `R_Y`.
- **How to read it:** in this simplified model, fatal disease must pass through D. D is the clinically complicated state, which WHF 2023 defines as established valve disease with overt complications such as heart failure, arrhythmia, stroke or the need for cardiac surgery.
  - A real death from C-level valve disease (sudden death, stroke, endocarditis) appears here as C→D progression followed by D mortality.
  - Under the annual cycle convention, someone entering D can die of RHD at the earliest in the next cycle.
  - So `p_D_Y` is a **D-stage RHD death probability for the complicated state**. It is not a case fatality for all RHD and not a clinical untreated-HF fatality.
- **No C-mortality sensitivity** and no direct C→death pathway.
- A single residual RHD-mortality multiplier is fixed at exactly 1.
- No severity ratios, no sex or age death multipliers.
- Never take the log of a structural zero.
- Write this reading into the registry note for `p_D_Y` and into the methods note.

**Optimiser bounds:**
- `p_D_Y`: the existing [0.03, 0.25]. Check this against the estimand: at 0.25, the derived O value is 0.44. Report any bound hit and never relax a bound silently.
- `m_O`: [0.2, 1.0], analyst chosen. Flag any fit outside the workbook's 0.3–0.8.

---

## 3. Onset: prevalence-derived, never GBD incidence

**GBD incidence is validation only.** Do not use it in 04 or 04b, or anywhere else, as:
- the healthy→disease inflow
- an onset curve, prior centre or starting value
- a target
- a fallback

This includes the retired Layer-1 alpha starting values. Setting its loss weight to zero is not enough; the dependency must be gone.

### 3.1 Deterministic reconstruction (04)

The reconstruction uses prevalence from **1990 to `CAL_YEAR_END`** only, so that 2020–2023 remains an independent holdout. If later years are ever used, state that the holdout is no longer independent.

1. **Prevalence surface.**
   - For each sex and year, turn GBD prevalence **rates** (per 100,000 → per person) by age group into single ages 0–94. Fit a smooth curve on the log-rate scale through the group midpoints and rescale it so each group's population-weighted mean is preserved.
   - Ages 95–100 take the 95+ rate.
   - `P[a,s,t] = rate × WPP Nx`: the same denominator as the forward model (`rate_x_pop`). Keep the published GBD Numbers for reconciliation.
   - Missing cells stay missing; a measured zero stays zero.
2. **Cohort balance, solved with the actual engine arithmetic.**
   - Equation: `P_end[a+1,t+1] = transported surviving prevalent stock + surviving new-onset contribution`.
   - Survival among prevalent cases uses the fixed natural history, GBD-proxy RHD mortality, background mortality and Y remission, all at the current iterate's stage/pathway composition.
   - Allocate onset by pathway (age at event) and entry mix. Solve for onset per **susceptible `N0`** person through the shared competing-risk split **exactly once**. Never divide a susceptible-based quantity by `(1 − prevalence)` again.
   - A cross-sectional age derivative or `P(t+1) − P(t)` is **not** acceptable.
   - Iterate as a fixed point: iteration 0 uses the entry mix as the stage composition; iterations 1–2 use the engine-implied composition from a forward run with the current onset.
   - Stop when total onset changes by < 1% (`BACKCALC_MAX_ITER <- 3L`).
3. **Negative residuals** signal noise, aggregation, migration, denominator differences or structural mismatch.
   - Smooth with a documented non-negative method; do not just floor cells to zero and call the result identified incidence.
   - Report the count, location and mass of negative cells (by age, sex and year).
4. **Fixed smoothing.** Produce a smooth non-negative reference hazard `h_ref(a,s,t)`, e.g. a quasi-Poisson fit of the non-negative balance on a **fixed** B-spline age basis by sex plus a linear year term, with an offset of log `N0`.
   - Report its degrees of freedom.
   - Keep it within the same parsimonious form used in calibration; no hidden high-dimensional fitting.
5. **Boundaries.**
   - Births enter age 0 disease-free.
   - The open 100+ group retains its stock.
   - Ages outside the observed support carry a documented extension.
   - Model predictions are always **aggregated back to the GBD age groups** for comparison; disaggregated cells are not independent observations.
6. **Exports:**
   - `onset_cohort_balance_raw.csv`: raw cohort-change residuals, including negatives.
   - `onset_reference_surface.csv`
   - `onset_backcalc_iterations.csv`
   - `onset_vs_gbd_incidence.csv`: labelled validation only.

Prevalence alone does not identify true incidence or entry stage. The reconstruction is **conditional** on the provisional mortality and stage assumptions, and the joint fit in 3.2 must be able to reconcile them.

### 3.2 Fitted onset correction (04b)

`h(a,s,t) = (h_ref(a,s,t) + h_floor(a)) × exp( Σ_k β_k B_k(a) + τ·(t − t_ref) [+ δ·1{male}] )`

- `B_k` is a **fixed** natural-spline age basis with **4–5 coefficients** (state the knots).
- `τ` is one period trend.
- `δ` is an optional sex offset, used only if diagnostics show it is needed. Default off.
- `h_floor(a)` is a small documented regularisation (e.g. 1% of the age-band mean of `h_ref`). It keeps zero reference cells fittable where the data support positive onset, with no GBD-incidence fallback.
- No independent parameter per age, sex or year, and no separately fitted age, period and cohort effects.

**Future onset:**
- The fitted age pattern at anchor `t_ref = CAL_YEAR_END`, continued with the fitted `τ`.
- Safeguards: `TAU_FUTURE_BOUNDS` and an optional linear taper of `τ` to 0 by `ONSET_TREND_TAPER_END`.
- **Retire `incidence_trend` from the active branch.** Never multiply it on top of `τ`, and never use future observed data.

---

## 4. One reduced calibration (04b) with the shared engine

**Free parameters, all of them:**
- `β_1..β_K` (K = 4–5)
- `τ`
- `log p_D_Y`
- `log m_O`
- optionally `δ`

That is about 7–8 in total. Assert that the parameter-vector names equal this allowed set.

Nothing else is fitted. Specifically, none of these: initial D fraction, remission or age multipliers, entry mixes, mortality ratios, CF/beta corrections, stage shares, care inputs.

**Repurposed sequence:**
- **04:** harmonised targets plus the deterministic onset reconstruction plus an initial bundle and diagnostics. Retire the well-sick-dead IR/CF Layer-1 calibration as a driver. Keep `$tp`/CF only as clearly marked diagnostics if a consumer needs them, and never populate fields for parameters that are no longer fitted.
- **04b:** a single joint fit, then the accepted historical trajectory, then the bridge seed.

**Objective:**
- `L = wP·MSLE(prevalence) + wD·MSLE(RHD deaths) + priors + λ·smoothness(β)`, keeping `wP:wD = 2.5:1` with the total held fixed.
- Prevalence here means **all-stage model prevalence** (section 1.2), mapped provisionally to GBD. The 3.1 reconstruction uses the same definition. WHF-classifiable prevalence is reported, not fitted, in the primary run.
- The MSLE is computed on GBD age groups × sex × year, on the `rate_x_pop` basis, with missing cells masked and never zero-filled.
- All-cause deaths, GBD incidence and the heart-failure (HF) proxy D share are **outside** the loss; report them as diagnostics.
- Do not treat correlated GBD UIs across years as independent (`USE_IV_WEIGHTS_STRUCT <- FALSE`).
- Priors:
  - `β` and `τ`: weak, centred on 0
  - `log p_D_Y`: log-sd 0.35 around the column-K anchor
  - `log m_O`: log-sd 0.5 around 0.48
- Export the data terms and the penalties separately.

**Death-target age floor:** `STRUCT_DEATH_TARGET_AGE_MIN <- 0L` (primary). Onset straight into D now makes child deaths possible. Report child cells separately: pooled agreement must not hide implausible child deaths. Run age floor 20 as a labelled sensitivity.

**Fit years and holdout:**
- Fit `STRUCT_MATCH_FROM` (2010) to `CAL_YEAR_END` (2019).
- 2020–2023 is the **holdout**: onset projected from the anchor with `τ`.

**Optimiser:**
- Seeded, bounded multi-start L-BFGS-B with a few reproducible starts, including the uncorrected reconstruction (`β = 0`, `τ = 0` relative to `h_ref`).
- Reuse the existing runtime controls; no large random search.
- Run a small profile or near-optimal-set analysis on `p_D_Y`, `m_O` and `τ` against the onset shape. Report weak identification honestly instead of adding levers.

**Historical initialisation:**
- Burn in with the shared engine under a documented pre-observation onset: the reconstructed 1991–1995 mean, held, with 1990 population.
- Check convergence at 2× burn-in, and check sensitivity to the prehistory length and to one alternative prehistory-onset assumption.
- Total prevalence does not identify the historical Y/O split. Show that uncertainty and add no fitted initial shares.
- If a prevalence-scaled seed is unavoidable, use it once, document it, and exclude the seeded year from any claim of fitted accuracy.
- **Never rescale to observed prevalence** during calibration, the bridge or projection.

**Bridge:**
- Forward-project the accepted state under **reference care** through 2020–2025 to `ANALYSIS_YEAR_START` (2026), on the same engine and rules.
- Check the bridge against the 2020–2023 GBD data.
- Never combine 2019 shares, 2023 prevalence and 2026 population.

**Population:**
- WPP stays exogenous. Demographic discrepancies and migration are absorbed by `N0`.
- Export the annual reconciliation residuals.
- Stock > population is a hard violation. It is never clipped or fixed by rescaling RHD.
- Newborn and open-age boundary flows must reconcile.

**Statuses:** report four separate statuses:
- structural validity
- optimiser convergence
- epidemiological fit
- identification

A numerically successful but poorly fitting candidate is **not** promoted to `best_parameters`. Save it as a failed candidate with its reasons. If the fixed structure cannot reproduce the data, show the age and year pattern of the failure and name the assumptions that need scientific review. Do not add hidden mortality scalars, remission or fitted transitions.

---

## 5. SAP, care cascade and surgery (fixed programme assumptions, never fitted)

### 5.1 SAP effect

SAP acts as a hazard multiplier `1 − c(1 − HR)` on A→B, B→C and C→D only, applied **once**, after the O multiplier on the untreated hazard. It has no effect on onset, regression or (by default) D death.

`c` is the effective protection of the **SAP-eligible** group the stock belongs to (5.2). For any stock that is ineligible under the active mask, `c = 0` and the multiplier is exactly 1. In the primary specification no A stock is eligible, so the A→B HR has no effect; it acts only in the "A ≤ 20 eligible" sensitivity.

| HR | Value | Source |
|---|---|---|
| A→B | `log(1−0.003)/log(1−0.075)` ≈ 0.039 | GOAL stage-specific (`Inputs!B64`); pooled 0.094 (`B66`) as a sensitivity |
| B→C | `log(1−0.025)/log(1−0.119)` ≈ 0.20 | `B65` |
| C→D | 0.45 at attained age < 40; 0.45 at ≥ 40 | `B70:B72`, Coates 2021 assumption. The engine must accept an [age × sex] multiplier |
| D death | 1 | `B73`; sensitivity 0.9 |

Stage mappings and the extrapolation to later onset and adults are provisional; keep their sources visible.

### 5.2 Pooled care cascade with WHF-consistent confirmation and an eligibility mask

This replaces the current `pmin(screen, diagnosis, treatment)` ceiling and the two-child/one-adult encounter schedule in the active branch. Do not apply both.

**SAP eligibility mask** (a derived mask on biological stocks; applied **before** any pooled coverage, volume or protection):

| Model group | Primary SAP eligibility | Switch |
|---|---|---|
| A, age ≤ 20 (Y-A) | **Off** | `SAP_ELIGIBLE_A_CHILD <- FALSE`; `TRUE` = sensitivity (WHF 2023: SAP "can be commenced" for stage A, class 2A) |
| Y-A, age > 20 | **Off** (not WHF-classifiable) | none: always off |
| B–D, either pathway, any age | **Eligible** | WHF 2023: SAP recommended for stage B; for C and D according to local recommendations |

- Diagnosed shares, on-BPG shares, delivered SAP person-years, effective protection and SAP costs are computed **only over the eligible stock**.
- Ineligible stocks get no SAP volume, cost or effect.
- Diagnosis of a stage-A child (≤ 20) is still recorded as a diagnostic volume in the primary run, even though it leads to no SAP.

**Screening and confirmation.** Confirmation is a **service flow, not a disease state**; add no disease or care-status states for it.
- **One routine screen per person between ages 5 and 39 inclusive** (`W = 35`). Annual probability among eligible **never-screened** people: `s = 1 − (1 − T)^(1/W)`.
  - When `T` ramps over calendar time, call this mapping an approximation and report the cohort coverage actually achieved.
  - Track the ever-screened share by age and sex through ageing and survival. Encounters are bounded by the never-screened stock.
  - Pre-2026 cohorts are assumed never screened (`Inputs!B17`).
- **The exam type depends on attained age at screening.** This follows WHF 2023: screening criteria apply only to ages ≤ 20, and confirmatory criteria apply to active case finding over age 20.
  - **Ages 5–20:** initial screen → confirmatory echo for **every** positive screen that completes referral → confirmed diagnosis.
    - Positive screens include false positives.
    - Confirmation demand = all positive screens × `confirm_completion` (fixed proportion).
    - A confirmed diagnosis = true positives × `confirm_completion`.
    - Assumption: confirmation is definitive. False positives are ruled out at confirmation and generate cost but no diagnosis or SAP. Uncompleted true positives stay undiagnosed and remain detectable clinically.
  - **Ages 21–39:** a single **comprehensive echo using adult confirmatory criteria**. It is costed directly per encounter, and its positives are confirmed diagnoses.
    - **No** separate confirmatory echo is charged, because the one exam already confirms. Never charge both.
- **Positivity, needed for confirmation demand.** Sensitivity alone cannot set the cost, because false positives drive demand. Implement both options and choose with `SCREEN_POSITIVITY_MODE`:
  - `"sens_spec"` (default): positives = `sens × detectable + (1 − spec) × non-detectable` among those screened.
  - `"observed_rate"`: positives = `observed_positivity_rate × screened`, with true positives capped at `sens × detectable`. Report the implied false-positive count, and stop if it is negative.
  - The **detectable** stock is the **WHF-classifiable** stock in the never-screened, undiagnosed population: `Y_A` at ≤ 20 plus B–D. Adult Y-A is non-detectable and counts in the specificity term.
  - Children (5–20) and adults (21–39) get separate sensitivity and specificity parameters (`screen_sens_child`, `screen_spec_child`, `exam_sens_adult`, `exam_spec_adult`). The child sensitivity defaults to the workbook's 0.9 (`B87`).
- **Clinical detection** applies at every age, to the undiagnosed **WHF-classifiable** stock.
- **Diagnosed and on-BPG shares.**
  - Pooled shares of the **eligible** stock, kept with `on-BPG ≤ diagnosed-eligible ≤ eligible stock` through new onset (which dilutes the shares), deaths, remission and age-related changes in eligibility.
  - Document the pooled approximation used when a diagnosed child's stock changes eligibility (e.g. Y-A ageing past 20, or A→B progression).
  - Treatment continues after people age out of screening, subject to retention.
- **Protection:** `c = on-BPG share × (good + φ·partial)`, applied once and only to eligible stocks.
- **Volumes** are separate quantities:
  - child screens
  - child positive screens (true and false)
  - confirmatory echoes completed
  - adult comprehensive echoes
  - new confirmed diagnoses (by age group)
  - diagnosed prevalent cases (eligible / ineligible)
  - SAP-eligible stock
  - on-BPG person-years
  - effectively protected person-years

  On-BPG person-years use the workbook's **mid-year** convention. Update the volume definition and every dependent check consistently.

**Care inputs** (linear 2026 → 2050 ramp, held afterwards; reference arm flat):

| Input | Ref | Scale-up 2026 → 2050 | Cell |
|---|---|---|---|
| Screening target T | 0.05 | 0.05 → 0.80 | `B77:E77` |
| Clinical detection | 0.02 | 0.02 → 0.02 | `B78:E78` |
| BPG initiation | 0.50 | 0.50 → 0.90 | `B79:E79` |
| Annual retention | 0.90 | 0.90 → 0.95 | `B80:E80` |
| Good adherence share | 0.60 | 0.60 → 0.85 | `B81:E81` |
| Partial adherence share | 0.20 | 0.20 → 0.10 | `B82:E82` |
| Partial relative efficacy φ | 0.5 | 0.5 | `B83` |
| Screening sensitivity, ages 5–20 | 0.9 | 0.9 | `B87` |
| Screening specificity, ages 5–20, or observed positivity rate | **required, no workbook value** | same | new |
| Adult comprehensive-echo sensitivity / specificity, ages 21–39 | **required, no workbook value** | same | new |
| Confirmatory completion proportion (positive screens, 5–20) | **required, fixed** | same | new |
| 2026 diagnosed share | IDN 0.05 / UGA 0.009 | — | `B88:B89` |
| 2026 on-BPG share | IDN 0.04 / UGA 0.0075 | — | `B90:B91` |

`screen_base`/`screen_target`, `diagnosis_*` and `treatment_*` in 00 block I are replaced by this table. Keep it as its own block, separate from natural history and economics.

**Rules for new required inputs:** specificity or positivity rate, adult exam accuracy, completion proportion, and the unit costs in 5.3.
- Search the repository and its documents for sourced, country-specific values first.
- If none exist, set a **provisional placeholder** labelled `[ASSUMPTION – analyst chosen, PROVISIONAL]`, with no literature attribution.
- Placeholders must appear in the registry export, the completion report's open issues, and a `provisional_inputs` field in `rhd_economic_results.rds`.
- Never present a placeholder as sourced.
- Run a sensitivity on specificity or positivity, because it drives confirmation cost.

**Historical care:** the same reference-care mechanism runs inside calibration and production. Care before 2026 is a documented **fixed** backward extrapolation of the 2026 reference inputs, not observed history.

**Reference arm:** its inputs are constant but its realised stock shares evolve. Retire the assertions that require them to be constant.

### 5.3 Echo unit costs (00 block J → consumed only by 08)

| Cost | Applies to | Note |
|---|---|---|
| `cost_screen_child` | each child screen, ages 5–20 | Existing `cost_screen_per_person` (IDN 1.10 / UGA 0.40); confirm that it means an initial screening encounter |
| `cost_confirm_echo` | each **completed** confirmatory echo after a positive child screen, true or false positive | New and required, country-specific |
| `cost_comprehensive_echo_adult` | each adult screening encounter, ages 21–39 | New and required. It replaces, never adds to, the screening cost for adults |
| `cost_sap_per_year` | each mid-year on-BPG person-year (eligible stocks only) | Existing |

No monetary values outside block J and 08.

### 5.4 Surgery

`surgery_enabled <- FALSE` by default. Off means:
- zero delivered surgeries
- zero surgery cost
- neutral surgery multipliers in calibration and in **both** arms

Equal nonzero coverage in both arms is not "off". Keep the code path and parameters so `TRUE` restores the current behaviour, and keep zero-valued legacy surgery fields for compatibility.

---

## 6. File responsibilities

| File | Responsibility |
|---|---|
| `R/abcd_engine.R` | The **one** implementation of pathway dynamics (`N0`, `R_Y`, Y A–D, O B–D), hazard competing risks, onset allocation, ageing, observation mapping (all-stage and WHF-classifiable masks), the SAP eligibility mask, and pooled care, screening and confirmation bookkeeping, used by reconstruction, calibration and production. Update the light and full modes together (identical disease outcomes). Validate invalid inputs instead of clipping them. Keep exact zeros. Update the header, including the D-only reading (1.2 and section 2 mortality rules). |
| `00` | Single configuration: `rhd_natural_history`, onset boundary, onset block (reconstruction, basis, `τ` safeguards), allowed fit vector, care block, `surgery_enabled`, seeds and runtime controls (`QUICK_CHECK <- FALSE`, bounded full-run settings). Remove contradictory legacy settings. Add a `stopifnot` block: entry mixes sum to 1 with O-A = 0, probabilities in [0,1], starts inside optimiser bounds, Y ordering A→B < B→C < C→D and A→NoRHD > B→A > C→B ≥ D→C (`Inputs!B47`). |
| `01` | Keep annual GBD prevalence, deaths and incidence, Rate + Number, UIs, units and age groups for 1990–2023. Validate keys and the required observations. Incidence is flagged validation-only. |
| `02` | Unchanged unless historical or bridge coverage or cohort alignment requires a change. |
| `03` | Pathway-specific fixed natural history as [age × sex] matrices, entry mixes, D-only mortality, care definitions, SAP HR set, surgery switch and registry export. No active GBD-incidence inflow (rename it `gbd_incidence_validation`). Remove the stage-split, `rhd_d_fraction` and severity-ordering checks. |
| `04` | Targets plus the deterministic Lexis onset reconstruction (3.1). Remove the Layer-1 optimiser as a driver. |
| `04b` | Reduced joint fit (section 4), accepted historical trajectory, versioned accepted bundle, bridge seed with the pathway dimension and the `R_Y` stock. |
| `05` | Consume the fitted onset, registry, pathway seed and care configuration as they are: no refitting, no second application of corrections, no replacing onset history with age. Remove the legacy seed fallbacks; stop if specs differ. |
| `06` | Run ref and SAP from identical accepted states. Keep Y/O through ageing. Trace every onset, entry-stage, transition, death, screening, detection and treatment flow. Update the parallel worker exports. |
| `07` | Keep the existing tables. Add `rhd_pathway_stage_long` and provenance/definition metadata. Reconcile totals with no double counting. |
| `08` | Cost each volume exactly once: child screens, completed confirmatory echoes (all positives × completion), adult comprehensive echoes (no extra confirmation charge) and delivered mid-year on-BPG person-years (eligible stocks only). Surgery cost = 0 when off. Add these cost components to `rhd_budget_impact.csv`. Keep the GNI-anchored VSL and discounting. Fix the null-effect checks (section 8). |

---

## 7. Output contracts and cache safety

**Keep these filenames:**
- `disease_model_inputs.rds`, `calibrated_rhd_parameters.rds`, `baseline_state.rds`
- the per-location RDS (`$wsd`, `$stages`, `$diag`, `$meta`)
- `rhd_model_long`, `rhd_stage_model_long`, `rhd_annual_flows`
- `rhd_budget_impact.csv`, `rhd_economic_summary.csv`, `rhd_economic_results.rds`

**Column rules:**
- Aggregate A/B/C/D = sums across pathways, with A = Y-A only.
- `new_rhd_A` must mean actual entry into A. Add `onset_total`, `onset_Y`, `onset_O` and entry-stage flows, and update consumers. Do not relabel all onset as A.
- Append `Y_A, Y_B, Y_C, Y_D, O_B, O_C, O_D, R_Y`, the screening, diagnosed and on-BPG stocks, and the volumes, as additive fields.
- Append `prev_all_stage`, `prev_whf_classifiable`, `Y_A_adult_unclassifiable`, `sap_eligible`, `sap_ineligible_diagnosed`, `screens_child`, `screen_pos_true_child`, `screen_pos_false_child`, `confirm_echo_completed`, `exam_adult_comprehensive` and `new_dx_confirmed`.
- The legacy `n_screened` = child screens + adult comprehensive exams (encounters), with the breakdown above alongside it.
- `sick` in `$wsd` stays the all-stage stock, used for population accounting.
- 07 adds `prevalence_definition` (`all_stage` / `whf_classifiable`) to the prevalence summary output.

**Spec fingerprint:** add a model-specification hash covering:
- pathways and boundary
- fixed registry values
- onset method and basis
- care assumptions (eligibility mask, screening/confirmation design, positivity mode)
- prevalence mapping
- input-file identities
- fitted window

Every bundle carries the hash. 05, 06, 07 and 08 reject a mismatch with a precise "regenerate with …" message, so old calibrated bundles and baseline caches are refused. Before writing the new spec, copy the previous valid outputs to a dated `*_prev_spec/` folder. Failed candidates are saved separately and never overwrite valid results.

**Full recalibration is required.** Provide `QUICK_CHECK` (synthetic or small grid plus a short real-data smoke run) and a bounded full run. Never launch the old high-dimensional calibration. If data or R packages are missing:
- finish all the work that can be implemented;
- state the exact blocker and the rerun commands;
- never substitute synthetic country data.

---

## 8. Validation and acceptance

### 8.1 Structural and numerical checks
Use small synthetic cases inside the scripts or engine, run in `QUICK_CHECK`; no large test framework.

1. **Conservation.** Finite, non-negative stocks and flows; hazard-split event and stay probabilities sum to 1; mass is conserved across onset, remission, deaths, ageing and population reconciliation, with stated absolute and relative tolerances.
2. **Exact O structural zeros and D-only mortality.** No O-A, O remission, O C→B or O D→C. RHD deaths in A, B and C are exactly 0 in both pathways (there is no C-mortality path or sensitivity), all RHD deaths come from D, and background death stays active in every state.
3. **Boundary ages.** Ages 19, 20 and 21: onset at 20 → Y and at 21 → O. A Y case crossing 20 keeps Y transitions and still regresses.
4. **Onset bookkeeping.** Entry-stage flows sum to total onset and use the `N0` denominator. Zero onset gives zero new cases. `R_Y` never produces O onset.
5. **Only allowed parameters move.** Fixed Y transitions are bitwise unchanged by optimisation. `m_O` acts on both O hazards and never on Y.
6. **Mandatory incidence-independence check.** Delete or radically perturb GBD incidence while holding prevalence and deaths fixed. Reconstructed onset, the fit and the baseline must not change; only the incidence-validation tables may.
7. **Light vs full engine.** Identical disease outcomes. Calibration and production reproduce the same trajectory from the same state, years, rates and reference care.
8. **Surgery off.** Zero delivery, cost and effect in every path.
9. **Screening and care bounds.** Screening only within ages 5–39. The never-screened stock reconciles. `on-BPG ≤ diagnosed ≤ eligible`. Treatment persists after age 39.
10. **Null effect.** Identical arms or HRs = 1 give identical disease trajectories. Costs follow delivered volumes. No positive benefit is demanded.
11. **Workbook mechanics replication** (a mechanics check, not calibration). Feed the engine a single closed cohort with the `Pop1_Ref` inputs, with workbook regression scenario 2 (zero O regression) and residual branches off to match this spec. Recompute the corresponding workbook values, or document each convention difference (start-of-year vs mid-year, residual handling). Compare RHD deaths and stage stocks in 2026, 2030, 2050 and 2085.

### 8.2 Epidemiological credibility
Write all of this into `output/<COUNTRY>/reports/baseline_calibration_report.md` plus CSVs and plots.

- **Fit tables.** Observed vs predicted prevalence and RHD deaths by sex, GBD age group and year: absolute error, relative error (when defined), fitted/held-out flag and target-inclusion status.
- **Fit summaries** for ≤ 20, > 20, 0–39, 40+ and all ages, with children inspected separately.
- **Declared provisional criteria.** Last fitted year, all-age totals, both sexes pooled, `rate_x_pop` basis:
  - prevalence within ±15% (`STRUCT_FIT_TOL_PREV`)
  - RHD deaths within ±25% (`STRUCT_FIT_TOL_DEATH`)
  - the same criteria reported for the 2020–2023 holdout
  - subgroup failures and bound hits reported separately
  - thresholds are never relaxed after seeing results, and difficult cells are never suppressed
- **Onset.** Onset rates and entry-stage allocation by age and period, compared with GBD incidence as **validation only**. Note that this compares two GBD outputs and is not clinical validation.
- **Stocks and mortality.** Y/O share by attained age, adult Y-A share, D share vs the HF proxy (3.21% IDN / 0.86% UGA), realised D-death risk, conservation and bridge continuity.
  - The workbook implies a D share of about 25% by 2050.
  - Expect tension and report it as a definition-sensitive diagnostic, not a target.
- **External plausibility.** Compare with Cannon/GOAL and the workbook's `Cannon_Validation`, noting the differences in stage, treatment and sampling. Cannon S2 (ages 15–24) is not a clean > 20-onset cohort. Do not force exact reproduction.
- **Sensitivities** (only after the primary run works, and fixed-assumption only):
  - Y B→A 0.0072
  - an alternative entry mix
  - historical prehistory / Y-O initialisation
  - O/Y death HR 1 and 4
  - death-target age floor 20

**Acceptance:** a baseline is labelled **credible** and passed through 05–08 as the primary result **only if** structural validity passes **and** the declared fit criteria pass. Otherwise, save the failed-fit diagnostics, state which criteria failed and what assumption is implicated, and do not call the baseline final.

### 8.3 Economics
- Reconcile the aggregate tables with the pathway outputs. 08 uses delivered volumes exactly once.
- Keep deaths averted and incremental costs **signed**. Remove the 06, 07 and 08 checks that fail when the SAP arm is not better; report the result instead.
- BCR = benefit / positive discounted incremental cost. Zero benefit with positive cost gives BCR = 0. Cost-neutral or cost-saving cases get an explicit undefined status, and undefined cost-per-outcome ratios are flagged, never fabricated.
- Do not import the workbook's illustrative costs or DALYs as epidemiological parameters.

---

## 9. Completion report (end of task)

1. File-by-file description of the changes.
2. `natural_history_registry.csv`: cells, sources, transformations, fixed/fitted flags, bounds, provisional labels.
3. A short generated methods note covering: cohort reconstruction (the equations and numerical method actually implemented), prevalence mapping, onset-stage allocation, lifelong pathways and `R_Y`, population reconciliation, reduced calibration, historical initialisation, the bridge, and baseline-care assumptions.
4. The exact free-parameter list and count, with start and final values, bounds, objective decomposition, seed, runtime, the four statuses and the fit diagnostics.
5. Accepted historical fit, 2026 baseline, and ref/SAP outputs through 08, **only if acceptance passes**. Otherwise, the failed-fit diagnostics.
6. The commands executed and their results (quick check and full run), the rerun instructions and any real blockers. Clearly separate code that was verified from calibration that was actually run.
7. Open issues, at least:
   - whether GBD prevalence includes A-like disease
   - the entry mix
   - Y B→A evidence (screen- vs clinic-detected disease)
   - O irreversibility
   - the D definition vs the HF proxy
   - historical care extrapolation
   - the surgery service when it is switched back on

Work autonomously within this scope and resolve routine details conservatively, documenting consequential choices. Simplify the model and expose its limits rather than adding fitted mechanisms nobody asked for.
