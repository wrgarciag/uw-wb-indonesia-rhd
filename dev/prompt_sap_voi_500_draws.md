# Claude Code prompt: isolated SAP-focused VOI analysis

Implement and run the **smallest SAP-focused value-of-information analysis** for Indonesia and Uganda as a **separate analysis**.

**Do not modify the current production pipeline, its inputs, calibrated bundles, caches, outputs, reports, or slides.** Production scripts 00–08 and `R/abcd_engine.R` must remain unchanged. Do not source `00_run_all.R` or rerun production to perform this analysis.

## Scope

Compare only the two existing options:

1. Reference care (`ref`).
2. Current screening-and-SAP scale-up (`sap`).

Vary only these four SAP inputs:

- M_AB: full-adherence progression-hazard multiplier for A→B.
- M_BC: full-adherence progression-hazard multiplier for B→C.
- `sap_protection_age_max`.
- `partial_adherence_relative_protection`.

Keep C→D and direct D-mortality SAP multipliers fixed at **1**. Hold calibrated natural history, the initial state, care-cascade inputs, surgery assumptions, and unit costs fixed. Recalculate service volumes and total costs for every draw.

## Create these new files

| File | Responsibility |
|---|---|
| `scripts/voi/00_run_voi.R` | Independent configuration and entry point |
| `scripts/voi/R/abcd_engine_voi.R` | Isolated copy of the current engine, minimally extended to separate A→B and B→C SAP multipliers |
| `scripts/voi/01_sample_sap_parameters.R` | Joint SAP sampling and prior diagnostics |
| `scripts/voi/02_validate_voi.R` | Production-equivalence and analysis checks |
| `scripts/voi/03_run_psa.R` | Paired reference/scale-up projections and economics |
| `scripts/voi/04_analyse_and_plot.R` | Information-value estimates, summaries, and figures |

Also create a short `scripts/voi/README.md` documenting execution, assumptions, limitations, and deferred work.

Declare analysis settings and final prior specifications for this analysis in the new `00_run_voi.R`; implement sampling functions in `01_sample_sap_parameters.R`.

## Reference workbook and parameter specification

**Consult `RHD_EVPPI_parameter_register.xlsx` as a reference** (including any supplied copy with a filename suffix such as `RHD_EVPPI_parameter_register (2).xlsx`). It provides evidence notes and proposed distributions; it is **not an authoritative runtime input, user interface, or final specification**.

- Verify parameter names, definitions, units, age boundaries, hazard-versus-probability interpretation, adherence mapping, and current point values against the **current production pipeline**.
- Final parameters and distributions adopted for this separate analysis must be aligned with that pipeline's mechanisms and explicitly coded in the **new VOI scripts**. Do not change production to match the workbook.
- Where the workbook is stale or inconsistent, use the current pipeline for implemented behavior and document the resolution in the VOI README.
- Do not read the workbook at runtime to populate parameters.
- Distinguish implemented production values from proposed uncertainty distributions. Explicitly coding a distribution does not make it expert-agreed: retain its provisional status where appropriate.

## Implementation requirements

- Implement the analytical workflow in **R**, prioritizing **`data.table`** for tabular data processing, joins, grouped summaries, reshaping, and CSV input/output (`fread()` / `fwrite()`). Use idiomatic `data.table` operations, including `:=`, `.()`, and `rbindlist()` where appropriate. Preserve the engine's existing matrix/array representation for age-by-sex calculations, and use `copy()` before by-reference changes to shared inputs. Use suitable R packages for regression and plotting.
- Inspect the actual current code and saved outputs before implementing.
- Read existing country-specific `baseline_state.rds` and `rhd_economic_results.rds` without modifying them.
- Reuse the saved applied 2025 unit costs, annual VSL schedule, discount rate, discount base year, and analysis horizon. Do not repeat price conversion on already-adjusted costs.
- In the isolated engine, support separate `sap_hr_A_to_B` and `sap_hr_B_to_C`, with backward compatibility for the existing shared multiplier. Preserve parameter validation and prepare effect matrices correctly for every draw.
- Load the VOI engine in an isolated environment/process so production functions and constants are not overwritten.
- Use the same sampled SAP parameters in both arms. Rebuild the care schedule when partial-adherence protection changes.
- Write every new result, checkpoint, log, and figure under `output/voi/<run_id>/<country>/`. Never overwrite production files.
- Record input fingerprints, engine version, random seed, prior specifications, and the fixed calibration assumptions.

## Priors

Use the reviewed parameter register as a starting point, preserving a justified dependence structure between M_AB and M_BC. Clearly identify proposed priors as provisional and unagreed.

Resolve the evidence-to-parameter mapping explicitly: pooled GOAL progression is not automatically an A→B estimate, and an RF-recurrence RR is not automatically a progression hazard multiplier. Distinguish efficacy extrapolation from the separately modelled adherence response. Keep mixture weights fixed within each prior scenario.

## Validation and execution

Before the PSA, reproduce existing reference and scale-up projections and economic results using current point values, within documented numerical tolerances. Confirm that the engine extensions do not change those results.

Run a small pilot, then **500 paired parameter draws per country**: each draw is evaluated under both reference care and scale-up. Set `N_DRAWS <- 500L` in the new configuration script. **Do not automatically increase the number of draws now; we can increase it later.** Report convergence and Monte Carlo limitations at this initial sample size. Individual-parameter EVPPI estimates are required for all four SAP inputs; label estimates exploratory when the initial sample is insufficient for stable estimation.

Check parameter bounds, dependence, numerical failures, and information-value stability. Save compact annual summaries sufficient to recompute valuations without rerunning projections, and support extending the saved draws in a later run.

If required saved inputs are unavailable, report the specific blocker; do not fabricate results or regenerate production outputs.

## Analysis and deliverables

- Compute joint SAP information value directly from simulated net benefits. Because SAP parameters are the only uncertain inputs, joint SAP EVPPI equals **EVPI conditional on the fixed model assumptions**.
- Calculate and report **separate single-parameter EVPPI estimates for each of the four SAP inputs**: M_AB, M_BC, `sap_protection_age_max`, and `partial_adherence_relative_protection`. These are required outputs alongside the joint SAP EVPPI, to assess which uncertainties within this cluster carry the most decision value.
- For each individual EVPPI, learn only that parameter while integrating over the conditional joint distribution of the remaining SAP inputs. Preserve dependence, particularly between M_AB and M_BC. Do not replace this with one-way sensitivity analysis that fixes the other inputs at point values.
- Use regression-based conditional mean estimates for continuous parameters (for example, a suitably checked single-input GAM). If protection age is sampled categorically, estimate its EVPPI from category-specific mean net benefits and category probabilities. Report category counts and Monte Carlo uncertainty; this is a stratified estimate, not an exact calculation.
- Reuse the same 500 paired draws for all individual-parameter estimates. Bootstrap complete paired draw rows and refit the estimator to assess precision. Check fit, stability, and the theoretical bounds 0 <= individual EVPPI <= joint SAP EVPPI, allowing for estimation error. Investigate apparent violations; do not silently clip estimates to hide problems. If estimation is unstable, still report the exploratory estimate with diagnostics and flag that more draws are needed later. Do not automatically increase the draw count.
- Save a tidy comparison table and a comparison figure for each country, showing all four individual EVPPI estimates and the joint SAP value, with uncertainty intervals where estimable. Include clearly labelled threshold curves and a ranked comparison at explicit reference thresholds, since ranking can change with valuation.
- Explain that individual EVPPI values **are not additive contributions** and must not be summed or displayed as a stacked decomposition of joint SAP EVPPI. Correlations and interactions can create overlapping information value. Interpret the comparison as the value of learning each input separately, rather than a unique allocation of total VOI.
- Report expected incremental net benefit, probability of positive INB, expected opportunity loss, and Monte Carlo uncertainty.
- Produce country-specific curves across clearly labelled valuation thresholds. Distinguish VSL-based valuation from US$/DALY thresholds; identify the existing deaths-to-DALYs conversion as a proxy.
- If comparing with trial costs, describe EVPPI as an upper bound on research value. Do not conclude that a particular trial is worthwhile without EVSI.

## Pending: later extensions only

Flag the following as pending **in the VOI README and code comments only**, without adding pending placeholders to figures:

**Expanded SAP analysis:** add a third option—scale-up with BPG continued after B→C—and uncertainty in the C→D effect. This requires a treated-C stock, continuation/dropout rules, and corresponding service costs. **Do not implement it now.**

Also document that this first analysis conditions on the existing calibration and initial state; calibration sensitivity and trial-specific EVSI remain later extensions.

Proceed with implementation, validation, and the separate analysis. Finish by summarizing files created, checks completed, results, and remaining limitations.

