# Minimal change: make the assumed 55% SAP progression-hazard reduction the default

Implement this change in the current production codebase for **Indonesia, then Uganda**, recalibrate both countries, and update/render the existing report and technical deck. Make only the edits required below. **No refactors, unrelated renaming, reformatting, line-ending conversion, new model layers, or changes to other assumptions.** Read the current files first; line numbers below are approximate. Preserve unrelated user edits.

## 1. Evidence framing and unchanged assumptions

Use this evidence statement in the relevant comments, metadata, and report/deck:

> Manyemba et al. (Cochrane, PMID 12137650) report RR = 0.45 for recurrent rheumatic fever (55% RRR), from one trial. This is not a hazard ratio and not a valve-progression estimate.

Use this exact label for the new default:

> assumed 55% reduction in progression hazards (expert-informed scenario), informed by evidence on prevention of recurrent rheumatic fever

The scenario is the **computational default**, not an empirically estimated progression effect. Applying the same assumed full-adherence multiplier, M = 0.45, to A→B and B→C is an explicit extrapolation. Never call it a “Manyemba progression HR.” Do not convert the recurrence RR with `log(1 - RR*p0)/log(1 - p0)` using the model's B→C risk: the required p0 would be the untreated **RF-recurrence** risk over the corresponding follow-up.

Keep these four decisions separate:

| Decision | Required default |
| --- | --- |
| Progression effect | M = 0.45 on treated A→B and B→C hazards only; same effect on both transitions is assumed |
| Protection age | Attained age ≤39, i.e. <40; independent of the natural-history boundary at age 20 |
| Eligibility and duration | Diagnosed A at age ≤20, and diagnosed B; no new age stop for stage B BPG. Existing initiation, retention, dropout, stage eligibility, and A→adult-prior-A rules remain unchanged. C/D remain ineligible |
| Adherence | Preserve `a = good + rel * partial` and `h_BPG = h_untreated * [1 - a*(1 - M)]`; M is assumed full-adherence efficacy, which the published RR does not establish |

Stage B patients aged ≥40 may remain on BPG and incur costs, but receive no quantified progression benefit in this model. The protection cutoff does not imply BPG stops at 40. “No age stop” does not mean lifelong guaranteed treatment: existing dropout and stage-exit rules still apply.

**Do not change** `adult_prog_hazard_mult = 0.48`: it modifies untreated adult B→C and C→D hazards, independently of treatment. Also preserve D mortality, calibration targets/window/bounds/tolerances/objective/starts/seeds, screening, care schedules, surgery, economic parameters, and all other model assumptions.

## 2. Preserve the old comparison before overwriting outputs

Before changing or rerunning anything, preserve the existing country-specific calibration/results/metadata and render inputs needed for an old-versus-new comparison, using a dated backup outside the live output paths. Record the old engine specification, SAP settings, country, analysis horizon, and economic conventions. Preserve existing work; do not clean or delete unrelated outputs.

Compare to the **actual pre-change outputs**, not to a reconstructed baseline presented as historical output. If an old output is missing, stale, or inconsistent, mark it unavailable or non-comparable and explain why; never invent a value. Do not rerun the old model unless necessary and clearly identified as a separately reproduced comparator.

## 3. Required edits by file

### `scripts/00_run_all.R`, block H

Replace the active GOAL-derived progression defaults with:

```r
sap_progression_rrr          <- 0.55
sap_progression_hazard_mult  <- 1 - sap_progression_rrr   # M = 0.45; ASSUMED full-adherence efficacy
sap_protection_age_max       <- 39L                      # attained age <=39 (<40); independent of nh_child_last_age
sap_effect_transitions       <- c("A_to_B", "B_to_C")     # same M on both transitions is an assumption
```

- Comment out the existing GOAL arm-risk inputs, conversion helper, and child HR assignments as a documented alternative, retaining their existing values and provisional stage mapping. Do not leave an active path using the old child/adult defaults or add a scenario-switch framework.
- Remove the active `sap_hr_B_to_C_adult` default; the new multiplier/cutoff replace the three child/adult HR fields.
- Keep `sap_hr_C_to_D`, `sap_hr_D_to_rhd_death`, and `sap_hr_prior_A_to_B_adult` equal to 1.
- Keep `sap_mechanism = "progression"` and the existing internal `sap_effect_scale = "hazard_ratio"` identifier to avoid unrelated schema changes. Explain that this identifier describes the implementation scale, not the published evidence.
- Add the evidence label and brief comments separating benefit age from treatment eligibility/duration. Keep block I's `sap_eligibility = c("A", "B")` and care inputs unchanged.
- Document the unchanged adherence equation and the full-adherence assumption. Current reference adherence weight is 0.70; do not change it.

### `R/abcd_engine.R`

- Set `ABCD_ENGINE_SPEC` to `"abcd_first_diagnosis_v1_adult_prog_mult_sap55"`, including its header mention, so old bundles fail existing compatibility checks.
- In `abcd_sap_hr()`, replace lookups for `sap_hr_A_to_B_child`, `sap_hr_B_to_C_child`, and `sap_hr_B_to_C_adult` with `eff$sap_progression_hazard_mult` and `eff$sap_protection_age_max`.
- Construct the protection mask from **attained ages**, replicated over sexes to match existing `[age × sex]` matrices. For both returned `A_to_B` and `B_to_C` matrices, use `prot * M + (1 - prot)`, where `prot = (ages <= sap_protection_age_max)`. Preserve return names, dimensions, and downstream interfaces. Do not use `child_m`/`adult_m` to determine protection age.
- Validate: M is one finite numeric scalar in `(0, 1]`; the cutoff is one finite integer-valued scalar within the modeled age range; `sap_progression_rrr` is in `[0, 1)` and agrees with `1 - M`; and `sap_effect_transitions` is exactly `c("A_to_B", "B_to_C")`. Missing/invalid settings must fail clearly. Keep the three unaffected HRs constrained to 1.
- Update only comments describing the changed SAP mechanism/evidence, including stale header text.
- **Leave `abcd_one_cycle()` arithmetic unchanged**, especially `mAB`/`mBC`: it already applies `1 - adh*(1 - HR)` once. Do not multiply probabilities directly, adjust adherence twice, or apply M a second time outside the engine.
- Preserve BPG person-year accounting. The current `n_effectively_protected` diagnostic is `on-BPG * adherence`; it is not proof of a progression benefit at age ≥40. Do not silently redefine this diagnostic or modify its arithmetic; describe it accurately wherever used in the report/deck.

### `scripts/03_build_disease_model.R`

- In `effects`, add `sap_progression_rrr`, `sap_progression_hazard_mult`, `sap_protection_age_max`, and `sap_effect_transitions`, all read from 00. Remove the three obsolete child/adult HR entries. Retain the unaffected HR fields.
- Rewrite `sap_effect_sources` with the exact evidence framing, assumed transition mapping, attained-age cutoff, and adherence interpretation. Correct the nearby stale comment suggesting a supported direct-mortality mode or C→D SAP effect.
- Update the numeric `prob_params` vector and SAP message to use the new fields. **Do not put the age cutoff or character transition list in a probability vector.** Validate those separately through the engine.
- Retain existing mechanism, eligibility, and scale validation. Messages must distinguish an assumed progression-hazard multiplier from a published recurrence RR.

### `scripts/04b_calibrate_structural.R`

- Search for `sap_hr`, `sap_effect`, and `spec$`; update every affected record of the SAP inputs used in reference care, including `fixed_tab` and `sc$spec$sap_hr`.
- Keep the existing `sap_hr` metadata container name to minimize changes, but store the new four fields plus the three unaffected HR fields as a **named list**, preserving numeric/integer/character types. Do not flatten mixed types with `unlist()`.
- In the numeric fixed-parameter table, record RRR, M, cutoff, and the unaffected HR values; record the transition mapping in an appropriate existing text/source or specification field. Keep parameter/role/value vector lengths aligned. Use the required evidence framing for SAP entries without relabeling unrelated natural-history parameters as treatment assumptions.
- Reference-care SAP must run through the same updated engine during burn-in, calibration, bridge, and projection. Do not add a separate calibration harness or post-hoc correction.

### `scripts/05_build_baseline.R`

- Update the 04b-versus-03 SAP consistency check to compare all new fields and the three unaffected HR fields, preserving their types, plus existing eligibility checks. Require the complete expected set of fields: missing fields must fail rather than allow a vacuous NULL-versus-NULL match.
- Fix `meta$sap_hr`, which currently requests nonexistent `sap_hr_A_to_B` / `sap_hr_B_to_C` fields. Keep the container name and record the same typed SAP specification used by 04b.
- Add `sap_progression_rrr`, `sap_progression_hazard_mult`, `sap_protection_age_max`, and `sap_effect_transitions` to `model_assumptions`; retain evidence/source metadata. Old specifications must require regeneration, not compatibility bypasses.

### `scripts/06_run_prevention_model.R`

Add **one sanity check** alongside the existing output checks: both `sap_prog_mult_A_to_B` and `sap_prog_mult_B_to_C` must be finite and equal to 1, within the existing numerical tolerance, at attained ages above `sap_protection_age_max`. Obtain the cutoff from the persisted model specification, not a hard-coded 39. No other edits, except directly stale SAP wording if present.

### `scripts/04_calibration_random_tp.R`, `07_make_outputs.R`, `08_economic_evaluation.R`

No functional changes are expected. Costing must continue using all on-BPG person-years, including stage B ages ≥40. Update comments only if they incorrectly describe the active default as GOAL-derived. Any strictly necessary compatibility edit discovered elsewhere must be minimal and explicitly itemized in the delivery.

## 4. Focused verification and production reruns

Use focused checks or temporary verification scripts; do not build new infrastructure or alter unrelated tests.

Verify before the long reruns:

- All modified R files parse; obsolete child/adult HR names remain only in the commented GOAL alternative or clearly historical material.
- Engine matrices have M at ages 20, 21, and 39 and 1 at age 40 and above, for both sexes and both transitions. Matrix values do not confer A eligibility above 20; the existing eligibility rules still govern treatment.
- At protected ages, `a = 0` gives 1, `a = 1` gives 0.45, `a = 0.90` gives 0.505, and reference `a = 0.70` gives 0.615. Above the cutoff, the effective multiplier is 1 regardless of adherence.
- Invalid/missing settings and old engine bundles are rejected; new SAP metadata survives the 03→04b→05 handoff without NULL fields or type coercion.
- Adult natural-history multiplier 0.48, C→D SAP multiplier 1, D-death SAP multiplier 1, prior-A multiplier 1, and BPG eligibility/accounting remain unchanged.

Then run **Indonesia first, Uganda second**, in separate fresh R sessions through the existing `00_run_all.R` entry point. The required sequence is **03 → 04 → 04b → 05 → 06 → 07 → 08**; 06 and 07 must not be skipped. The current 00 sources 01 and 02 as well: running that existing full pipeline is acceptable. Do not add skip/resume machinery just to avoid these steps. Use the existing country selector and country-specific paths; restore the selector to its original value after both runs.

Rebuild disease inputs and **force actual recalibration**, rather than reuse an old SAP bundle or stale cache. Preserve numerical calibration settings. Do not relax tolerances, change targets, or disable checks to obtain a pass. If the existing authorized workflow permits preliminary outputs despite poor epidemiological fit, retain that behavior and clearly report the failed fit status. A completed run is not a passed fit check. On a genuine execution error, retain logs and diagnose only within this scope; do not fabricate completion.

## 5. Update and render the report and technical deck

Find the existing files; use their actual paths and spelling. The deck may be named **`slideck-technical.Rmd`**, rather than `slidedeck-technical.Rmd`. Do not create a duplicate or rename it. Search both documents and directly used SAP tables/includes for `GOAL`, `hazard ratio`, `0.0385`, `0.1998`, `sap_hr`, and `not extrapolated to adults`.

- Update only affected SAP text, parameter tables, assumptions, checks, and result-linked content. Retain GOAL as a documented alternative/direct childhood progression evidence where relevant; do not describe it as the active effect.
- Show the evidence distinction: recurrence RR 0.45 = recurrence RRR 55%; progression M = 0.45 is **assumed**, not estimated. Use the exact expert-informed scenario label.
- Explain the four separate decisions in section 1, including treatment costs without modeled benefit at age ≥40 and the unchanged A age restriction/C-D ineligibility.
- Calculate the adherence examples from persisted M and care inputs. Reference `a = 0.70` yields 0.615; illustrative `a = 0.90` yields 0.505. Clearly distinguish illustrative adherence from the actual reference schedule.
- Pull all model parameters, fitted values, and results from current persisted inputs/outputs; derive displayed values and examples inline. The published evidence figures and a clearly labeled illustrative adherence value may be stated explicitly. Do not hard-code result numbers or falsely imply published recurrence evidence was estimated by the model.
- Explain existing adherence-weighted BPG exposure diagnostics accurately; they can include people outside the progression-benefit age range.
- Knit both files after both country runs. Retain checks against stale or incompatible country outputs. Report the actual rendered filenames and whether each render succeeded. If files/dependencies/data are absent, report the concrete blocker rather than create replacement reports or claim success.

## 6. Deliverables

1. **Unified diff per changed source file**, preserving original formatting and line endings.
2. **Rerun log summary** for each country: executed steps, recalibration/cache status, engine/SAP specification, completion status, fit-tolerance status, and render outcomes, with paths to full logs.
3. **Old-versus-new comparison per country**: calibrated transition probabilities (identify child/base versus derived adult quantities), fit metrics and existing tolerance flags, deaths averted, BCR, and cost per DALY. Use matching scenarios, horizon, currency/price year, discounting, and existing metric definitions; label units and distinguish discounted from undiscounted outcomes. Include absolute changes where comparable. State unavailable/non-comparable values explicitly and identify their source files.
4. **Scope audit**: list any changed line outside the requested scope, with the file/line and reason; state “None” if there are none. Distinguish expected regenerated outputs/rendered artifacts from source edits. Do not hide collateral changes in generated results or reformatting.

Proceed with the authorized edits and runs. Report actual outcomes and unresolved blockers precisely; do not treat this scenario as clinical guidance or as an empirically established valve-progression effect.

Evidence references supplied with this task:

- [Review discussing recurrent rheumatic fever evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC7017848/)
- [Manyemba et al., Cochrane, PMID 12137650](https://pubmed.ncbi.nlm.nih.gov/12137650/)
- [GOAL trial, PMID 34767321](https://pubmed.ncbi.nlm.nih.gov/34767321/)
