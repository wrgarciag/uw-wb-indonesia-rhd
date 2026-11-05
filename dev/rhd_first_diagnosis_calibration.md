# Simplify the RHD model and establish a credible calibrated baseline

Work in the existing `work-uw-rhd-modelling` codebase. Plan first, implement and test the upstream model and calibration, and pass the checkpoint below **before modifying scripts 05–08**. Implement one default specification. Do not build sensitivity alternatives.

## Goal

Produce a simplified but plausible RHD model that distinguishes disease history from first confirmed diagnosis, represents how screening actually finds existing disease, and reproduces GBD prevalence and RHD deaths more closely. The current problems are underestimated prevalence and overestimated deaths above age 20. Improve both jointly, giving prevalence greater calibration weight than deaths. A completed numerical run is not sufficient evidence of a credible baseline.

The conceptual reference is `dev/RHD_first_diagnosis_model.xlsx` (locate the matching workbook in that directory if its filename has a suffix). It illustrates first diagnosis and follow-up, not a complete national population model. Use it only as a reference. **Do not reference spreadsheet cells in your plan, code, comments, or final explanation; do not read the workbook at runtime.** The explicit inputs and their intended roles are specified below.

## Scope and preservation requirements

- Modify only the necessary lines in existing `scripts/00*.R` through `scripts/08*.R`, including `04b_calibrate_structural.R`, and `R/abcd_engine.R`. This is an allowed scope, not a requirement to edit every file. Leave 01 and 02 unchanged unless a demonstrated dependency requires a small correction.
- Keep the current architecture, country switch, naming conventions, parallel execution where valid, and output locations. Support Indonesia and Uganda with their existing country-specific inputs.
- Put **every user-editable assumption, input override, parameter centre, bound, calibration weight, threshold, and execution-stage control in `00*`**. Downstream scripts should consume the persisted settings and validate them. Avoid introducing competing hidden defaults.
- Do not refactor the whole pipeline, introduce another modelling framework, or add unrelated scripts, folders, dependencies, or scenario branches. Use existing diagnostic destinations and existing test infrastructure; temporary test harnesses may remain outside the repository.
- Preserve the current calendar settings: `CAL_YEAR_START = 2000`, `CAL_YEAR_END = 2019`, and `STRUCT_MATCH_FROM = 2010`. Thus 2000–2009 remains initialisation/validation and 2010–2019 remains the scored structural period. Preserve the existing observed, projection, analysis, and intervention-ramp windows as separately defined in 00. Do not widen the calibration period to improve fit.
- Outputs must continue feeding the existing `report.Rmd` and `slideck-technical.Rmd`. You may read them to inspect dependencies, but **do not edit or render these files yet**. Do not modify other source files outside the allowed scope.
- This task authorizes recalibration under the changed structure. Preserve the previous results through the existing last-valid/checkpoint mechanism; never relabel an old calibration as compatible with the new model.

## 1. Plan and diagnose before editing

Read the relevant scripts, engine, existing calibration diagnostics, and workbook. Explain the minimum changes needed and identify which parameters are fixed, calibrated, derived, or diagnostic only.

Trace the current causes of the adult mismatch, including:

1. Fixed stage mixes at biological/modelled incidence versus stage at first diagnosis.
2. Automatic transfer of surviving stage A to B at the age boundary.
3. Lifetime early/late-onset risk differences inferred from model entry age.
4. The adult-onset D-death hazard multiplier of 2, which currently changes Indonesia's base D probability of 0.119 into approximately 0.223839.
5. Progression into D, D occupancy, reference care, burn-in population balancing, and the bridge into the analysis period.
6. Population-denominator differences between GBD Numbers and GBD rates multiplied by WPP populations.

Do not assume that changing entry labels alone fixes mortality or prevalence. Show how the revised stocks and flows address the mechanism producing the mismatch.

## 2. One coherent default: population burden plus first diagnosis

Retain a national population model containing undiagnosed disease, diagnosed disease, normalised history, and people without prior modelled disease. A diagnosed-only screening cohort cannot reproduce total GBD prevalence/deaths or estimate the population effect of screening. Undiagnosed people must continue to progress, normalise where allowed, and die.

- **First confirmed diagnosis is entry into the diagnosed/care cohort.** It can occur in A, B, C, or D when age-eligible. It can result from screening with confirmation or clinical detection. Diagnosis age does not identify disease-onset age.
- **New disease and new diagnoses are different flows.** Screening detects existing disease. It must not create RHD, inflate biological incidence, or add the same person twice to prevalent disease.
- Carry sufficient stage-specific diagnosis and BPG counts to apply eligibility and protection correctly. Extend the existing state representation only as needed; keep one engine for calibration, bridge, preview, and production. Existing Y/O bookkeeping may be retained if this minimises changes, but its labels must retain an honest meaning and must no longer imply different hazards merely because diagnosis occurred later.
- Stage at diagnosis should follow the available undiagnosed stage stocks and the stated detection process. Do not reshuffle people's disease stages to force an illustrative diagnosis mix.
- Preserve clinical detection outside screening ages. Count first diagnoses through screening and clinical presentation without overlap. Confirmation failure leaves disease undiagnosed; it does not remove disease.
- Preserve the current simplified no-recurrence assumption after normalisation, with an explicit history stock. Do not interpret normalisation as proven anatomical healing.

### Disease entry and diagnosis composition

The workbook's illustrative first-diagnosis mixes are:

| Age at first diagnosis | A | B | C | D |
|---|---:|---:|---:|---:|
| Age ≤20 | 0.50 | 0.30 | 0.15 | 0.05 |
| Age >20 | 0 | 0.60 | 0.30 | 0.10 |

Record these in 00 as **provisional diagnosis-composition benchmarks**, not observed national stage targets, and not a mandatory fit term. Do not apply their 5%/10% D shares to every new biological disease entry.

For a simple explicit starting specification, use first **modelled classifiable disease** entry mixes `c(A=1, B=0, C=0, D=0)` at age ≤20 and `c(A=0, B=1, C=0, D=0)` at age >20. These are new simplifying assumptions, not estimates from the workbook or proof of the age of biological onset. Advanced first diagnoses then arise from progression of previously undiagnosed disease. Put these mixes in 00, keep them fixed initially, and assess whether they permit plausible calibration. If they cannot, diagnose the incompatibility and ask a focused question with alternatives rather than quietly adding direct C/D entry or new free parameters.

## 3. Explicit natural-history inputs

Use the following annual event probabilities as reference centres/starting values. They are provisional modelling inputs, not validated untreated WHF-stage transition estimates. Use attained age for the age-dependent rules. Initially calibrate the six nonzero childhood transition parameters within the existing bounds; use the fitted B→C and C→D values at adult ages as well, with no additional adult multiplier.

| Transition | Age ≤20 | Age >20 | Role |
|---|---:|---:|---|
| A → normalised history | 0.311 | 0 | Childhood parameter calibrated |
| A / adult prior A → B | 0.054 | 0.054 | Childhood parameter calibrated; adult prior-A uses the same fitted value as a provisional extrapolation |
| B → A | 0.177 | 0 | Childhood parameter calibrated; adult structural zero |
| B → C | 0.074 | 0.074 | One shared calibrated parameter |
| B → normalised history | 0 | 0.020 | Fixed adult extrapolation |
| C → B | 0.050 | 0 | Childhood parameter calibrated; adult structural zero |
| C → D | 0.112 | 0.112 | One shared calibrated parameter |
| D → C | 0 | 0 | Exact structural zero |
| A/B/C → RHD death | 0 | 0 | Exact structural zeros |

Keep these annual-probability calibration bounds in 00 unless the plan demonstrates a specific incompatibility requiring a decision:

```r
TRANS_BOUNDS <- list(
  p_A_to_no_rhd = c(0.005, 0.60),
  p_A_to_B      = c(0.001, 0.25),
  p_B_to_A      = c(0.002, 0.40),
  p_B_to_C      = c(0.001, 0.20),
  p_C_to_B      = c(0.001, 0.15),
  p_C_to_D      = c(0.001, 0.30),
  p_D_to_C      = c(0, 0)
)
```

### Mortality and age boundary

- Retain the existing country-specific base annual stage D RHD-death probability: **Indonesia 0.119; Uganda 0.0597**. Uganda's input must not be replaced by the illustrative workbook's Indonesia value.
- Apply the country's base D mortality to all disease/diagnosis histories. Remove the extra late-entry D-death hazard multiplier: its effective value is **1**, not 2. Keep `residual_rhd_mortality_mult = 1`. Do not lower base D mortality to manufacture fit.
- Remove the automatic 0.48 late-onset progression multiplier and any exploratory adult progression reductions. Their effective value is **1** in the single updated default. Do not infer slower progression solely from later diagnosis.
- Retain country- and age-specific competing non-RHD mortality from existing inputs. Do not substitute the workbook's illustrative background mortality table.
- For consistency with this reference, use age **≤20 versus >20** as the declared modelling convention, configured in 00. Retain the separate GBD reporting split of **0–19 versus 20+**; also report a 21+ diagnostic where needed to isolate the convention. Do not describe these two cutoffs as identical or assert that the convention independently verifies a clinical guideline.
- At the 20→21 ageing step, surviving A becomes a separate **adult prior-A lesion/history stock**. It is excluded from classifiable RHD prevalence and is ineligible for SAP. Do not force it into B, call reclassification remission, or drop people from the population. Its subsequent prior-A→B probability uses the shared A→B parameter, starts from the provisional centre 0.054, and is explicitly an extrapolation. Track that return into classifiable disease separately from first disease entry and first diagnosis; avoid repeat first diagnoses.
- D represents clinical complications. Keep D→C at zero; do not invent a valve-regression rate.

Convert annual probability inputs once to hazards `h = -log1p(-p)`. Apply care multipliers to the relevant hazards and combine competing events coherently. Preserve accounting under the existing annual cycle ordering; document when new disease entries can first face stage-specific death. Do not silently clip invalid transitions or repair population deficits by deleting disease.

## 4. Screening, confirmation, and care

Keep the current longitudinal programme settings from 00 rather than importing the workbook's synthetic cohort of 100,000 people, 1% prevalence, one screening round, 90% care coverage, or 25-cycle follow-up.

- Retain the existing single-screen programme window **ages 5–39**, inclusive, and its never-screened history: annual probability among never-screened people is `q = 1 - (1 - c)^(1/35)` for a constant window coverage target `c`. Preserve the existing interpretation during the ramp; do not silently substitute `2/15`, `1/20`, or the one-round workbook attendance assumption.
- Screening window coverage: Indonesia reference **0.05**, Uganda reference **0.01**; scale-up target **0.80**. Preserve the existing ramp.
- Screening sensitivity **0.90**. Add explicit confirmatory echo completion **0.90** and confirmation sensitivity **1.00**. Use perfect specificity as a declared simplifying assumption; count true positives and completed confirmation separately. Do not invent a confirmation unit cost. If no existing unit cost is available, flag its omission from total costs.
- Retain annual clinical detection **0.02** among undiagnosed disease. Define this as confirmed clinical diagnosis so the screening confirmation factor is not inadvertently applied twice.
- Retain reference/scale-up BPG initiation **0.50/0.90**, annual retention **0.90/0.95**, good adherence shares **0.60/0.85**, partial adherence shares **0.20/0.10**, and partial-adherence relative protection **0.50**. Coverage, retention, and adherence must have distinct denominators and roles.
- Retain the existing provisional initial diagnosis/on-BPG shares: Indonesia **0.05/0.04**; Uganda **0.009/0.0075**; initial ever-screened share **0**. Initialise these once, allocate consistently to eligible stage stocks, and carry them forward. Do not reset care histories in the analysis start year.
- SAP is optional in A at age ≤20 and modelled for confirmed B. Set default SAP eligibility to **A at age ≤20 and B**; C/D and adult prior-A have no SAP programme volume in this simplified default. Preserve specialist/surgical services separately.

Use these default SAP effects on progression hazards:

| Event | Age ≤20 | Age >20 |
|---|---:|---:|
| A → B | `log(1-0.003) / log(1-0.075)` ≈ 0.03853835 | 1 for adult prior-A |
| B → C | `log(1-0.025) / log(1-0.119)` ≈ 0.19982855 | 1 |
| C → D | 1 | 1 |
| D → RHD death | 1 | 1 |

These are explicit conservative choices for this reference-based default. The childhood stage mapping is provisional; the adult B→C effect is not extrapolated from the child trial. Do not annualise the derived hazard ratios again. Actual protection applies only to eligible diagnosed people receiving BPG, weighted by adherence. Screening alone has no treatment or mortality effect. Explain that adult B treatment has programme costs but no quantified progression benefit under the adult HR of 1.

Retain the current surgery service, country-specific need fractions/coverage, and effect inputs from 00; the workbook's omission of surgery does not require deleting the existing service. Keep it equal across arms unless current settings explicitly differ. Avoid applying any treatment benefit twice when interpreting clinically observed transition estimates.

## 5. Prevalence-led calibration using the production engine

- Preserve the existing Layer-1/Layer-2 handoff. Layer 1 supplies starts and diagnostics; it must not override the fixed stage D mortality or the final structural parameters.
- In 00, set Layer-1 weights **`W_PREV = 2.5`, `W_DEATH = 1`, `W_INC = 1`**, replacing the present death-led Layer-1 weighting.
- Retain structural prevalence-to-death weighting **2.5:1**: `STRUCT_PREV_DEATH_RATIO = 2.5`, `STRUCT_DATA_WEIGHT_TOTAL = 3.5`, yielding `W_PREV_STRUCT = 2.5`, `W_DEATH_STRUCT = 1`, and `W_INC_STRUCT = 0`. Incidence remains an anchored input and a validation output in Layer 2.
- Give adult prevalence an explicit contribution: split the prevalence term into equally weighted **normalised mean** errors for ages 0–19 and 20+, configured in 00. Keep `STRUCT_DEATH_CHILD_SHARE = 0.50` for equally weighted normalised childhood/adult death errors. Report raw and weighted components so a larger target grid or very small death cells cannot quietly dominate.
- Retain current log-count loss, masking, priors, smoothness, optimiser settings, incidence multiplier bounds **0.2–5.0**, `STRUCT_TARGET_BASIS = "rate_x_pop"`, and mid-cycle prevalence observation unless a demonstrated inconsistency requires a minimal correction. Record every consumed setting.
- Fit **total classifiable disease**, diagnosed plus undiagnosed, against the declared GBD prevalence mapping. Exclude adult prior-A and normalised history. Compare total RHD deaths from diagnosed and undiagnosed D against GBD RHD deaths. Keep all-cause mortality for validation.
- Preserve GBD Numbers alongside the rate×WPP targets. Explain denominator differences; do not claim exact agreement with published GBD counts from agreement only with rate×WPP targets. If the classifiable-stage mapping and GBD case definition differ, state the limitation explicitly.
- Keep the 120-year `per_capita` burn-in initialisation and its balancing trace initially. Test convergence against twice its duration. No rescaling of disease stocks in scored history, the bridge, or projection.
- Carry the calibrated end-of-window state through the reference-care bridge to the analysis start without reseeding stage shares, prevalence, diagnosis, or screening histories. Use exactly the same engine and reference care throughout.
- Do not use the illustrative diagnosis mixes or the GBD heart-failure proxy as validated population stage targets. Diagnose weak identification, prior dependence, optimiser convergence, and bound hits. An optimiser success code alone is not fit adequacy.

## 6. Required checkpoint before changing 05–08

Add only the minimal execution control in 00 needed to stop the initial run after 04b; sourcing 00 must not accidentally run 05–08. First implement and test the engine plus required changes in 00–04b. **Do not modify 05–08 yet.**

Run the upstream pipeline and structural calibration for Indonesia and Uganda using available data. Then use the shared engine directly from the calibrated bridged state for a small reference/scale-up preview through 2050 and the current analysis endpoint, without invoking or modifying 05–08. This preview is a check of anticipated future behaviour, not a new production pipeline.

Provide a checkpoint with:

1. Before/after fit for prevalence and RHD deaths, by country, sex, GBD age group and scored year, including all ages, 0–19, 20+, and the 21+ diagnostic. Show ratios against both the fitted targets and published GBD Numbers where available.
2. Unweighted and weighted losses; estimated versus fixed parameters; bound hits, convergence, and burn-in checks. Show D occupancy and realised D fatality, and why adult deaths changed without lowering base D mortality.
3. Checks for population and stage accounting, diagnosis/BPG eligibility, one-time diagnoses, prior-A reclassification, competing deaths, and the calibration→bridge→preview handoff.
4. Future trajectories of total and diagnosed prevalence, D occupancy, deaths, diagnoses, SAP person-years and screening encounters. Check continuity, age patterns, implausible accumulation, extinction, and unexpected jumps at age/year boundaries.
5. A clear **PASS or FAIL** for structural validity, numerical validity, epidemiological fit, and future plausibility, separately.

Use the existing absolute-log-ratio tolerances **0.15 for prevalence** and **0.25 for deaths**, configured in 00. Apply the checkpoint to all-age and adult prevalence and all-age/childhood/adult deaths across every scored year, with age/sex diagnostics reviewed for material hidden mismatches. These thresholds are operational criteria, not evidence of external validation. Do not loosen them or mask difficult target cells to obtain a pass.

**Proceed to 05–08 only if this checkpoint passes.** A pass permits continuation without another approval request. If it fails, iterate within the upstream allowed scope, explaining the limiting mechanism. If credible fit requires changing a fixed assumption, adding unsupported parameters, or choosing between materially different structures, stop with the evidence and a focused question offering 2–3 alternatives and your recommendation. If required data/runtime are unavailable, say which checks could not run and do not claim a pass or proceed on that basis.

## 7. After a pass: finish 05–08 with minimal edits

- Make the smallest handoff, state, and metadata changes needed in 05/06 to run the verified default. Bump the engine specification and reject stale inputs/calibration bundles. Store model assumptions, eligibility, observation definitions, country, target basis, weights, calendar windows, and fit status.
- Keep `$wsd` and `$stages` and existing filenames. Preserve the aggregate output column names/order: `scenario, age, cause, sex, year, well, sick, newcases, dead, pop, all.mx, intervention, location, eff_ir, eff_cf`. `sick` represents total classifiable disease; `newcases` remains new modelled disease, not diagnoses; `dead` remains RHD deaths. Preserve compatible stage/flow fields and append necessary diagnostics.
- Add honest outputs for undiagnosed/diagnosed disease, adult prior-A, confirmed first diagnoses, confirmation completions, and re-entry from prior-A. Update stock-flow identities to account explicitly for reclassification and prior-A→B. Preserve mass across these flows.
- Keep legacy Y/O columns only where their original history remains meaningful; do not fabricate historical data, silently redefine them as diagnosis histories, or fill unavailable quantities with zeros for compatibility. If a hard dependency in an unchanged report prevents honest compatibility, identify it and ask rather than editing the reports or concealing the conflict.
- In 08, cost screening encounters and eligible on-BPG person-years from the actual engine volumes. Preserve existing surgery costing and VSL valuation. No extra benefit for diagnosis itself, no direct SAP mortality multiplier, and no double-counting cases or deaths. Flag uncosted confirmation if its price is unavailable; do not invent a BCR improvement.
- Verify the report/deck input contracts without changing or rendering the Rmd files. If no update is needed in 07 or 08 beyond consuming compatible outputs, leave their implementation intact.

## 8. Focused verification and final handoff

Use meaningful checks for the changed mechanics: zero screening retains undiagnosed burden; detection alone conserves disease stocks; zero confirmation completion blocks screening-confirmed diagnoses; zero SAP benefit makes extra detection alone produce no health gain when all other care is held equal; equal arm inputs give equal outputs; stage A ageing conserves people without creating B; adult prior-A receives no SAP; no RHD death occurs outside D; and diagnosis/treatment counts remain nested within their eligible stocks. Check light/full engine equivalence and production/calibration equivalence at identical inputs.

Run the completed pipeline for both countries after the checkpoint and confirm that saved output schemas and accounting pass. Keep runtime reasonable by using the existing multistart/parallel controls and reusing compatible prepared inputs. Do not rerun expensive stages without a reason.

Conclude with changed files and why, the final parameter table with sources/assumption labels, a concise before/after calibration table, checkpoint results, commands run, production output locations, and any remaining limitations. Distinguish credible fit from numerical completion. Do not claim that any provisional transition, adult normalisation rate, diagnosis mix, or disease-entry assumption was clinically validated by matching GBD.
