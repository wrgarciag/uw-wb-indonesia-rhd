# Claude Code prompt: revise RHD programme populations, costing, VSL benefits, and reporting

Implement the following focused revision in the existing RHD A/B/C/D investment-case repository. Complete the code changes and update the report and technical slide deck. Preserve the pipeline from calibration outputs through baseline construction, projection, economic evaluation, and reporting.

## Scope: strict edit boundaries

Only these source files may be modified, and only where required for this revision:

- scripts/00_run_all.R
- scripts/03_build_disease_model.R
- scripts/05_build_baseline.R
- scripts/06_run_prevention_model.R
- scripts/07_make_outputs.R
- scripts/08_economic_evaluation.R
- analysis/report.Rmd
- reports/slideck-technical.Rmd

Treat the uploaded filenames with numeric suffixes as reference snapshots; edit the canonical repository paths above. Read the current repository versions before changing anything.

Read scripts/01_prepare_inputs.R, scripts/02_build_demography.R, scripts/04_calibration_random_tp.R, scripts/04b_calibrate_structural.R, and R/abcd_engine.R as needed to understand dependencies, but DO NOT modify them. Do not modify the Beamer preamble, bibliography, other reports, or unrelated files. Regenerate the normal derived outputs when feasible; do not manually edit result files. Do not introduce unrelated refactoring, new model mechanisms, or a new pipeline.

## Mandatory: reuse saved calibration; do not run calibration

Do not modify any scripts/04*.R files, including 04_calibration_random_tp.R and 04b_calibrate_structural.R, or calibration-related settings and logic elsewhere. Do not run, source for execution, or indirectly trigger the calibration process: it is computationally expensive and outside this task. Use the existing saved calibration files as read-only inputs to the downstream model. Do not overwrite or regenerate them.

Before running downstream scripts, verify the saved calibration files exist and match the configured country, model specification, and reference care. Load settings without executing the full 00 pipeline if it would trigger calibration; use the existing resume/skip mechanism or execute only the required downstream stages. If saved calibration files are missing or incompatible, complete the independent edits and report the exact blocker. Do not recalibrate, bypass compatibility checks, or silently substitute uncalibrated parameters.

## Inspect first

Trace how screening, diagnosis, treatment coverage, adherence, programme volumes, unit costs, VSL, and economic summaries travel through the allowed files. Inspect both Rmd consumers and their existing references.

The attached snapshots show several specific issues to address; confirm their current equivalents:

- 00 sets screen_age_restrict = FALSE and an old upper age of 15.
- 03 carries the older screening mask and cumulative cascade interpretation.
- 05 computes eff_diagnosis = pmin(screen, diagnosis) and eff_treatment = pmin(screen, diagnosis, treatment).
- 06 currently computes n_screened = popm * screen_mask * c_screen, so the default costs annual screening of the entire population.
- 08 already costs n_screened and n_on_optimal_treatment from the projection trace, but computes VSL from GDP per capita and adds a separate GDP-per-capita benefit.
- 05 verifies that reference treatment, adherence, and surgery coverage match reference care embedded in calibration outputs. Preserve that safeguard.

### Task 1 — correct programme populations and screening costs

#### A. Annual screening encounters

Use these explicit scenario assumptions:

| Group | Inclusive age window | Screenings across the window | Window length | Annual encounter multiplier |
| --- | --- | --- | --- | --- |
| Children/adolescents | 5–19 | 2 | 15 years | 2/15 |
| Adults | 20–39 | 1 | 20 years | 1/20 = 0.05 |

The adult multiplier is 0.05, not 0.04. Derive multipliers from screening counts and inclusive window lengths instead of scattering rounded constants through scripts.

For scenario s and year t, with common screening coverage c_screen,s,t and unit cost u_screen:

N_screen,s,t = c_screen,s,t × [(2/15) P_5–19,t + (1/20) P_20–39,t]

C_screen,s,t = u_screen × N_screen,s,t

Implement this at the model's existing age × sex × location × year resolution, then aggregate. Under these assumptions, ages below 5 and above 39 generate no screening encounters. All sexes are included unless the existing scenario explicitly specifies otherwise. Use the existing population arrays; do not substitute RHD prevalence for the screening denominator.

Put the age windows and screening counts in 00 and pass them through the existing 03 → 05 → 06 structures with minimal additive changes. Prevent old mask defaults from silently overriding the new primary scenario. If retaining a legacy sensitivity switch, make its interaction explicit and unambiguous.

These are expected annual encounter volumes, potentially fractional in a deterministic model. They approximate a programme distributed across each age window. They do not model exact individual screening dates, probabilities of treatment, or unique people ever screened. Do not implement individual screening histories in this revision.

In 06, update n_screened to the corrected annual encounter volume. In 08, retain the single costing step:

c_screen = sum(n_screened) × cost_screen_per_person

Do not apply the encounter multipliers again in 08. Apply the same scheduling definition to reference and scale-up, with their respective coverage trajectories.

#### B. SAP denominator and coverage

For this provisional run, SAP treatment volume is:

N_SAP,s,t = living prevalent RHD_s,t × c_SAP,s,t

Here c_SAP means the proportion of ALL prevalent RHD receiving SAP, not the proportion of the total population, the fraction screened that year, or a conditional cascade probability. Use living_rhd_start, consistent with the existing start-of-cycle treatment denominator.

Continue SAP in years without screening and after people age out of the screening windows. Clinical diagnosis can also lead to SAP independently of screening. Do not restrict ongoing treatment to ages 5–39.

Keep annual encounter multipliers completely separate from SAP coverage. In particular, never cap SAP by c_screen × 2/15, c_screen × 1/20, or n_screened/pop.

For this immediate revision, retain the existing pmin reach constraints where they preserve current treatment trajectories and calibration compatibility, but explicitly describe them as a simplifying programme-reach/cumulative-coverage ceiling, not a mechanistic screening-only route to diagnosis or treatment. The screening reach parameter and annual encounter count have different denominators and meanings. Do not claim that this establishes a measured care cascade.

Do not add an explicit clinical-diagnosis pathway or change the coverage policy to implement one. Document that independent clinical diagnosis is subsumed by the provisional coverage assumption and that the reach ceiling is a modeling simplification. If an existing repository option already permits uncapped coverage, preserve it without changing the primary setting.

Maintain the existing SAP adherence convention: n_on_optimal_treatment is delivered treatment person-years, and the adherence modifier attenuates health protection once. Do not apply adherence twice or automatically lower programme costs by the adherence multiplier.

Do not implement guideline-specific eligibility, prophylaxis duration, initiation/discontinuation, or retention states. Explicitly list these as limitations of prevalence × coverage.

#### C. Preserve the calibrated health model

Do not change calibrated transition probabilities, mortality structure, mortality bounds, initial disease stocks, SAP effect sizes or mechanism, surgery effects, demography, or reference-care parameters as part of this revision. Preserve exact A/B/C RHD mortality zeros when D-only mortality is selected, and competing other-cause mortality in every stage.

The screening revision should change encounter volumes and screening costs while leaving health trajectories unchanged when SAP coverage, adherence, and other health inputs are unchanged. Do not multiply SAP health effects by screening frequency.

Maintain 05's reference-care compatibility check against calibration metadata. Do not bypass it, reinterpret stored calibration values, or silently use uncalibrated fallbacks. If the current repository contains a genuine incompatibility that cannot be resolved within these boundaries, report the precise conflict instead of modifying calibration code or outputs.

#### D. Label the evidence correctly

The WHO 2024 guideline supports population echocardiographic screening at ages 5–19 in populations with moderate/high RHD prevalence. The proposed adult window and both screening frequencies are scenario assumptions; do not present them as WHO-prescribed schedules.

Verified source: WHO guideline on the prevention and diagnosis of rheumatic fever and rheumatic heart disease, summary of recommendations:
https://www.ncbi.nlm.nih.gov/books/NBK609692/

Use an existing matching bibliography entry if present. Do not invent citation keys or edit the bibliography. If missing, provide a clearly identified source link/note in the allowed Rmd files.

### Task 2 — correct benefit monetization in 08

Use VSL alone to monetize deaths averted. Compute the baseline VSL using GNI PER CAPITA and project it forward using the GDP-per-capita growth path. Remove the additive GDP productivity benefit.

For country j, let b be the GNI/VSL anchor year:

VSL_j,b = m_j × GNIpc_j,b

VSL_j,t = VSL_j,b × (GDPpc_j,t / GDPpc_j,b)

With a constant real GDP-per-capita growth assumption g_j:

VSL_j,t = VSL_j,b × (1 + g_j)^(t − b)

This implements unit income elasticity over time for this revision. Retain the existing multiplier as a transparent provisional assumption, now applied to GNI per capita, unless an already documented repository specification establishes another GNI-based transfer method. Do not invent a new multiplier or silently change the established valuation methodology. Clearly identify its provenance or provisional status.

Requirements:

1. Add explicit GNI-per-capita anchor value/year parameters in 00 and consume them in 08. Keep income anchor years distinct from analysis, calibration, price, and discount reference years.
2. Find a documented country-specific GNI-per-capita value in existing inputs/references, or verify a value from an authoritative source. Record source, indicator definition, year, currency, and price basis. Do not relabel existing GDP-per-capita values as GNI or silently fall back to GDP. If GNI is unavailable, leave an explicit required parameter with an actionable validation error; complete other work and clearly report the blocked economic run.
3. Use the existing GDP-per-capita projection when available. If the repository provides only gdp_growth, retain it as an explicit constant-growth scenario assumption; do not call it an externally sourced forecast. Do not use aggregate GDP growth as per-capita growth without population adjustment.
4. Align any projected GDP series to the GNI/VSL anchor year before taking ratios. Keep VSL, unit costs, and discounting on a consistent documented monetary basis. Do not compound nominal growth into constant-price values. If source conversion is needed, make it explicit and reproducible.
5. Apply country-specific parameters for the currently supported countries. Avoid an Indonesia fallback for a Uganda run.

For each location/year:

deaths_averted_t = RHD_deaths_ref,t − RHD_deaths_sap,t

benefit_vsl_t = deaths_averted_t × VSL_t

benefit_total_t = benefit_vsl_t

DF_t = (1 + discount_rate)^−(t − discount_base_year)

PV_benefits = Σ_t DF_t × benefit_total_t

PV_incremental_costs = Σ_t DF_t × (C_sap,t − C_ref,t)

BCR = PV_benefits / PV_incremental_costs

Net_benefit = PV_benefits − PV_incremental_costs

C_s,t includes the existing screening, SAP, and surgery costs. Preserve current incremental rather than gross programme costing and preserve existing cost-effectiveness definitions outside this correction. Count each component once. Keep death differences signed; do not clamp adverse annual differences to zero.

Remove the GDP-per-capita add-on everywhere: calculations, console text, comments, tables, figures, equations, and report/deck narrative. GDP per capita remains a growth index, not a separately added benefit.

Preserve existing output schemas. If benefit_gdp is retained for compatibility in the incremental RDS table, set it to exactly zero and identify it as a retired compatibility field; ensure all consumers use the VSL-only total. Do not rename existing total_benefit, net_benefit, bcr, or other established result columns.

Check economic guards after the costing change. For positive incremental costs, compute BCR normally. For zero/negative incremental costs, do not fabricate a positive denominator or a finite BCR; report BCR as undefined and identify cost-neutral/cost-saving outcomes transparently. Do not fail solely because a legitimate revised scenario is cost-saving. Retain checks for missing/invalid inputs.

### Task 3 — update both Rmd consumers

Update analysis/report.Rmd and reports/slideck-technical.Rmd to describe the final implemented methods and regenerated results.

- Explain the two age windows, inclusive denominators, encounter formula, and scenario assumptions.
- Distinguish annual screening encounters, programme reach, prevalent RHD on SAP, and effective protection.
- Explain continued SAP between screenings and the provisional treatment-eligibility limitations.
- Replace GDP-based VSL and VSL + GDP equations with GNI-anchored, GDP-growth-projected VSL-only valuation.
- Make current results dynamic from existing output files and metadata: no hard-coded revised costs, deaths averted, BCRs, net benefits, or parameter values.
- Refresh relevant tables, figures, captions, messages, and appendix details while retaining the existing narrative and slide design.
- Explain that screening-only cost reductions improve BCR if benefits are held fixed and incremental costs remain positive. The combined revision also changes valuation, so the final BCR may rise or fall; report its actual direction from regenerated outputs. Do not imply greater biological effectiveness.
- Retain the existing DALY-per-death approximation and identify its limitations where discussed; do not introduce a new DALY model.
- If required outputs cannot be regenerated, explicitly mark the affected results as unavailable/stale. Do not present old results under the revised methods.

Use reports/beamer_preamble.tex and the existing bibliography without modifying either.

## Preserve output contracts

Maintain existing file paths, formats, list elements, table keys, scenario/location labels, horizons, and required column order, including:

- output/out_model/<location>.rds and its existing wsd/stages/meta structure.
- output/tables/rhd_model_long.csv and .rds.
- output/tables/rhd_stage_model_long.csv.
- output/tables/rhd_annual_flows.csv.
- output/tables/rhd_budget_impact.csv.
- output/tables/rhd_economic_summary.csv.
- output/tables/rhd_economic_results.rds, including summary, budget_impact, incremental, params, and discount_base_year.

Keep n_screened, now explicitly defined as annual screening encounters, and n_on_optimal_treatment as SAP person-years. Prefer small additive metadata/diagnostics over renaming fields. Keep 07 focused on compiling health/volume outputs; economic calculations remain in 08. Modify 07 only if necessary for corrected volumes, validation, or metadata.

## Validation and execution

Use focused checks in allowed scripts or temporary, untracked checks; do not create a new permanent test framework or modify disallowed files.

1. Before regeneration, preserve a temporary comparison of existing health and economic results when available.
2. Check inclusive boundaries: ages 4/40 contribute zero screening; ages 5/19 receive 2/15; ages 20/39 receive 1/20.
3. Hand check: P_5–19 = 1500, P_20–39 = 2000, c_screen = 0.8 gives 240 annual encounters; at u_screen = 1.10, screening cost is 264.
4. Verify per-age encounter counts sum to the formula for each location/scenario/year and frequencies are applied exactly once.
5. Verify SAP volume equals start-of-cycle prevalent RHD × realised SAP coverage, including outside screening ages. Changing screening frequency alone must not change treatment coverage, treatment volumes, or health outputs.
6. Confirm reference care still matches calibration metadata and calibrated parameters pass unchanged into projection.
7. Hand check: GNIpc_b = 2000 and m = 30 gives VSL_b = 60000; at g = 0.03, VSL_b+1 = 61800. These are synthetic validation values, not country inputs.
8. Confirm benefit_total = benefit_vsl = deaths_averted × VSL and any compatibility benefit_gdp is zero. Independently reconcile discounted benefits, incremental costs, net benefit, and BCR.
9. Check zero/negative incremental-cost handling and required output schemas. Compute pooled TOTAL BCR from pooled benefits/costs, not a mean of location ratios.
10. Parse changed R files, run the affected downstream pipeline using existing calibration outputs where possible, and render both Rmd files. Inspect report/slide rendering for broken tables, equations, or slide overflow.
11. Do not run calibration at all for this revision, directly or indirectly. Reuse compatible saved calibration files without modifying them. If 00 normally runs calibration, use the repository's existing resume/skip mechanism or run only downstream scripts with the required settings loaded. Do not alter calibration settings just to run the revision.
12. Use git diff/status to confirm that source edits are restricted to the eight allowed paths and no unrelated changes were introduced.

If data, GNI parameters, R dependencies, LaTeX, or outputs are missing, complete everything that can be verified, describe the exact remaining blocker, and distinguish executed checks from unexecuted ones.

## Final response

Summarize the changes file by file, the final screening/SAP denominators, GNI/VSL source and growth assumptions, preserved calibration/output contracts, validation and rendering results, and any blockers. Where comparable results exist, report before/after screening volumes, costs, benefits, and BCR from actual outputs, separating the screening-cost correction from the valuation correction. Do not invent numeric results.

