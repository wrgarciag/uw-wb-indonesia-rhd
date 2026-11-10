# Implement fixed age-specific stage D mortality and recalibrate

Update the existing RHD model to test a more realistic age pattern of stage D mortality and improve agreement with GBD, particularly the currently overestimated adult RHD deaths. Improvement is a hypothesis to evaluate, not a result to assume. Assess prevalence and childhood mortality alongside adult deaths.

## Scope

- Make minimal, task-specific edits only to the existing `scripts/00*.R` through `scripts/08*.R` files, including `04b_calibrate_structural.R`. Modify `R/abcd_engine.R` only where necessary to support the age-specific mortality consistently.
- Do not create new source files, helpers, tests, documentation, or folders. Do not edit any other existing files, including reports, slide decks, workbooks, or raw inputs. Normal pipeline-generated caches, calibration bundles, diagnostics, and outputs may be regenerated through the existing workflow.
- Avoid refactoring, renaming unrelated objects, formatting whole files, or changing unrelated assumptions. Preserve existing user edits. Inspect all stages for affected references, but edit only those that need a change.
- Declare all new parameter values, assumptions, and any relevant prior/starting-value settings centrally in `scripts/00_run_all.R`; downstream scripts and the engine must consume them without hidden constants or defaults.

## Required mortality specification

Replace the single country-specific stage D annual RHD-death probability with two fixed inputs:

| Country | Attained age <20 | Attained age ≥20 |
| --- | ---: | ---: |
| Indonesia | 0.16 | 0.03 |
| Uganda | 0.06 | 0.03 |

These are fixed annual probability assumptions, not hazard multipliers, estimated priors, or free calibration parameters. They must remain exactly fixed during random-TP and structural calibration; do not fit, sample, penalize, or optimize them. Retain calibration of the other currently fitted parameters under the existing objective, weights, priors, bounds, and calibration window.

- Select mortality by **current attained age**, regardless of age at onset, first diagnosis, screening, or treatment history. Apply the same base age-specific probabilities to all stage D histories.
- For mortality specifically, age 20 receives the adult value, 0.03. Declare this mortality boundary explicitly in `00`; keep the existing `nh_child_last_age = 20` convention for progression, entry, remission, and other non-mortality rules unchanged. Do not silently move those boundaries.
- Preserve D-only RHD mortality: direct RHD deaths from A/B/C remain exactly zero. Preserve the existing competing-risk calculation, probability-to-hazard conversion, background mortality, and surgery effects. The table specifies the base D probability before existing service effects and competing risks; realized annual fatality need not equal it.
- Keep the residual mortality multiplier fixed at 1 and preserve current country-specific adult progression multipliers, SAP effects, entry mixes, and all other settings.

## Propagate the change consistently

Trace the current scalar `p_rhd_death_D` / `nhe$pd` through `03`, `04`, `04b`, `05`, `06`, `07`, `08`, and the shared engine. Use the smallest compatible change that supplies the correct age-specific values in burn-in, calibration, bridge years, baseline, and intervention projections. Check matrix dimensions and age alignment; avoid R vector recycling.

Update affected validation, fixed-input signatures, bundle metadata, handoff checks, diagnostic labels, and exports so both probabilities and the mortality age boundary are recorded accurately. Reject incompatible old bundles/caches through the existing mechanism. Replace scalar assumptions in fatality guards and diagnostics with correctly aligned age-specific checks; do not weaken the guards. Keep legacy fields only where needed for compatibility and ensure they cannot silently reinstate one mortality probability for every age.

## Execute in this order

1. Briefly state the affected locations and implementation plan, then proceed with the minimal edits.
2. Parse the edited R files. Run an **Indonesia-only preliminary test through `00`–`04b`**, using the actual existing pipeline. If needed, temporarily reduce the existing run controls in `00` for this preliminary test, then restore full production settings. Do not build a separate calibration harness or treat a reduced run as final calibration.
3. Verify fixed base D mortality at ages 19, 20, and 21; applicability across D histories; zero direct A/B/C RHD mortality; valid probabilities, competing risks, and population accounting; and exclusion of the mortality inputs from the fitted parameter vector. Check that metadata and diagnostics use the new age specification. Fix task-related failures before proceeding.
4. Once the preliminary test passes, run **fresh, full production calibrations and the complete `00`–`08` pipeline for Indonesia and Uganda**, including `04b` and `05`–`08`. Reuse the completed Indonesia calibration only if it already used the full production settings and the final code. Do not reuse pre-change calibrations. Continue automatically after successful checks; no additional approval is needed.
5. Compare against available pre-change results without creating extra comparison files. Report modeled versus GBD deaths and prevalence, by <20 and 20+, over the existing calibration window and at its endpoint, plus the existing calibration loss. Aggregate both model and GBD into identical age groups; flag unavailable age detail rather than relabeling ≤20 / 21+ as <20 / 20+. State whether adult mortality fit improved and whether childhood mortality or prevalence fit worsened.

Finish with a concise account of files changed, checks passed, full runs completed, fixed values actually used, and the main calibration tradeoffs. Report failures or poor fit candidly; do not adjust other assumptions to manufacture improvement.
