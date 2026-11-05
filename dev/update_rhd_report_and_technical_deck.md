# Update the RHD report and technical slide deck

You are working in the RHD prevention investment-case repository. Implement the updates; do not stop at a plan.

## Goal and strict edit scope

Update ONLY:

1. `analysis/report.Rmd`
2. `reports/slideck-technical.Rmd`

Bring both documents into agreement with the current implemented methods and the latest existing saved results. Every numerical result and result-dependent interpretation must be computed dynamically from persisted output files at render time.

Do not edit any other source, configuration, bibliography, preamble, model, input, output, or documentation file. Do not create shared helper scripts, regenerate model outputs, rerun calibration, or source the pipeline to obtain missing values. Keep any new reporting helpers inside these two Rmd files. Preserve unrelated user changes. Inspect `git status` before work and check the scoped diff afterwards. Render artifacts and intermediates only in a temporary directory outside the repository; disable caches or other render-time writes into the repository.

## Read-only grounding

Read both documents fully, then inspect these current sources without changing them:

- `scripts/00_run_all.R`
- `scripts/01_prepare_inputs.R`
- `scripts/02_build_demography.R`
- `scripts/03_build_disease_model.R`
- `scripts/04_calibration_random_tp.R`
- `scripts/04b_calibrate_structural.R`
- `scripts/05_build_baseline.R`
- `scripts/06_run_prevention_model.R`
- `R/abcd_engine.R`
- Relevant saved-output builders in scripts 07/08, solely to understand schemas, units, economics, and aggregation.
- `reports/beamer_preamble.tex`, existing model schematics, and `library/references.bib`.

If available, read the existing `dev/beautiful-deck-SKILL` for presentation guidance. Do not add dependencies or modify that guidance. Use only references already in `library/references.bib`; do not invent citation keys or assume a paper proves an implementation choice. Distinguish literature evidence, provisional mappings, and assumptions.

Treat executed code and saved artifact metadata as authoritative; some source comments still describe obsolete methods. Verify all field names and structures before using them. Current configuration alone does not establish which specification generated an existing output.

## Methods that must be represented correctly

### Two-layer calibration

- Layer 1 (`04`) is a prior-anchored aggregate well–sick–dead calibration, using broad age-band-by-sex incidence and mortality corrections, hazard-based competing risks, and bounded optimization with multiple starting points. Its filename does not mean the old simple random search is still used.
- Layer 2 (`04b`) jointly estimates incidence corrections and natural history using the actual production A/B/C/D engine. Layer-1 incidence multipliers are starting values, not fixed final incidence. Explain that the final Layer-2 incidence is transferred to production without applying correction factors twice.
- Show the actual loss, preferably as a compact equation: weighted mean squared log errors for prevalence, RHD deaths, and optional incidence, plus separately weighted mortality and structural priors, incidence smoothness, and optional usable stage targets.
- The current primary prevalence:death weight ratio is 2.5:1, total weight 3.5; direct incidence weight is zero. Incidence is still estimated, prior-anchored, and reported for validation. The current mortality prior weight is 1 and structural prior weight is 0.05. Obtain displayed settings from saved metadata when available, rather than inserting these numbers as reporting constants.
- Describe inverse-variance weighting separately for Layers 1 and 2; Layer 1 enables it by default and Layer 2 disables it by default. These are distinct settings, and GBD uncertainty intervals are not independent measurement errors across cells.
- Show the stored calibration window and scored years separately. The current configuration projects 2000–2019 and scores 2010 onward. Explain held-out years.
- Targets default to GBD rates multiplied by matching WPP age–sex–year populations, rather than mixing published GBD counts with WPP denominators. Show target reconciliation where available, label the selected basis, and do not compare a model/target ratio against the wrong denominator.
- Missing target cells and counts below the stored threshold are masked rather than treated as observed zeros.
- Explicitly disclose the current primary exclusion of death-target cells below age 20 from the loss. Retain childhood deaths in validation and distinguish scored adult deaths from all-age totals. Explain the unresolved acute/rapidly fatal childhood pathway issue as a structural limitation, without declaring exclusion clinically validated. Read the age floor from saved calibration metadata.

### Mortality and stage dynamics

- Primary `D_only`: RHD-attributable mortality in A/B/C is exactly zero; competing other-cause mortality remains in every living state. One stage-D mortality level is estimated and the residual multiplier is fixed at one. Do not plot a log of a structural zero or apply all-stage mortality-order checks to D-only results.
- Retain a correct explanation of `all_stage` only where relevant to saved results or as a labelled alternative. Remove the obsolete claim that a positive minimum stage-A mortality necessarily prevents fitting Uganda under the current D-only specification.
- Stage-D mortality prior centers come from provisional pooled workbook proxy estimates. They are not validated annual transitions or final calibrated estimates. Current D mortality is scalar across age/sex; age-stratified workbook values are not implemented. Show fitted values, prior centers, bounds, and their provenance distinctly, from persisted fields where available.
- Adult progression and A-to-No-RHD remission hazard multipliers provide parsimonious age variation at the stored age cutoff. They do not imply age-specific D mortality. Display the calibrated values and explain their extrapolation limits.
- Stage shares, including the D fraction, emerge from the calibrated dynamics; the current structural fit does not independently estimate a free D-share parameter. Remove claims based on obsolete comments to the contrary.
- Initialization uses zero-start burn-in, a convergence comparison with twice the burn-in, and forward projection of the calibrated final state through the bridge to analysis start under reference care. Distinguish the equilibrium initialization from the evolved historical state and the analysis-start seed. Do not imply the analysis-start seed was reset to a fixed Cannon split or exactly matched GBD prevalence.
- Shared annual sequence: entering-year stocks → competing within-year flows → year-end stocks → ageing. Prevalence observation is the stored convention, currently the mean of start/end stocks before ageing. Deaths and incidence are flows.
- Explain hazard-based competing risks with a short equation and proper residual survival. Surgery is a service affecting C→D and D mortality; it is not a separate postoperative stock.
- Population comes from WPP2024, whereas disease-model background mortality is GBD all-cause minus RHD. Do not describe WPP mortality or fertility as direct disease-engine inputs merely because WPP embeds demographic assumptions.

### Intervention mechanism and adherence

- Primary SAP mechanism is `progression`: pathway-specific hazard reductions on A→B, B→C, and C→D. No direct SAP mortality multiplier in this mode; mortality effects emerge through stage trajectories.
- Show `h_ij(c) = h_ij(0) [1 - c(1 - HR_ij)]`, where `c = effective treatment coverage × adherence`. This is a population-average hazard approximation, not an exact individual-level mixture of annual risks.
- Use the stored pathway effects. Explain the provisional GOAL-to-WHF stage mapping and the assumed C→D effect separately. The legacy direct mortality mechanism is an explicit comparator, not the primary interpretation.
- Adherence is implemented separately for each arm, with baseline/target trajectories. Defaults of one mean no additional attenuation of already effective treatment coverage; they do not mean universal treatment or an assumed adherence improvement. Show actual saved trajectories.
- Untreated natural-history parameters are fitted with reference treatment, adherence, and surgery applied inside the engine; production applies that care consistently once.
- There is no direct SAP effect on new-disease incidence, but do not force observed incident-case counts to be identical: changes in the susceptible pool can indirectly alter counts. Prevalence and stage composition can also change. Derive their direction from outputs.

## Dynamic results and artifact handling

Reuse the documents' existing read-only loaders and path conventions wherever practical. Discover and validate the actual schemas of:

- `data/<COUNTRY>/calibrated_rhd_parameters.rds`, especially `$layer1`, `$stage_calibration`, `$layer2_targets`, and `$meta`.
- `data/<COUNTRY>/baseline_state.rds` and `disease_model_inputs.rds` where needed.
- Existing `calibration_layer2_target_reconciliation.csv` and `calibration_structural_*.csv` diagnostics.
- Existing model outputs under `output/<COUNTRY>/out_model/` and economic outputs under `output/<COUNTRY>/tables/`.
- Existing GBD and population artifacts used by the documents.

Prefer the self-contained RDS bundles and use CSV diagnostics where they add information absent from those bundles. Do not assume every diagnostic is present or that a CSV and RDS from different runs are compatible.

1. Derive all years, country-specific parameters, coverage, adherence, counts, ratios, percentages, costs, discounting, BCRs, rankings, and stage shares from saved artifacts. Never paste a current result into prose, a slide title, speaker notes, captions, axes, or tables. Render methods settings dynamically from the saved run when supported; disclose genuinely unrecorded settings rather than substituting current globals as historical facts.
2. Make every data-dependent statement conditional or neutral. Replace assertions such as “Uganda cannot fit,” “both countries return the same ratio,” “prevalence increases,” or “deaths mostly occur outside D” with current-output-derived text. A result changing should update the narrative, not fail a guard that enforces the old conclusion.
3. Keep genuine integrity checks: location, scenario pairing, common years, dimensions, finiteness, units, aggregation identities, specification consistency, and model/economics agreement. Distinguish input consistency checks from hypotheses about the expected result. Handle zero denominators and nonpositive incremental costs explicitly without manufacturing ratios.
4. Check available metadata for agreement across calibration, baseline, model, and economics. Do not mix D-only/progression calibration with legacy projection or stale economic outputs. If a detected incompatibility affects a required result, stop with a specific diagnostic naming the incompatible files and upstream regeneration needed; do not regenerate them. Where provenance cannot be verified, state that limitation.
5. Missing optional diagnostics should produce a concise “not available in the saved outputs” note. Missing required artifacts should produce an actionable error. Do not silently substitute old outputs, priors, zeros, or a different country.
6. Keep the report country-parameterized and preserve the deck's existing Indonesia/Uganda comparison. Validate each country independently; do not assume identical fitted parameters or conclusions.
7. Compute cumulative deaths and incident cases as flows; present prevalence as a dated stock or an explicitly labelled person-year measure, never a sum described as unique prevalent people. Preserve annual versus cumulative and discounted versus undiscounted distinctions.
8. Reuse saved economics as authoritative. Cross-check any report-side summaries against those objects; do not reimplement economic evaluation or revise its assumptions. If the existing valuation is per RHD death averted, do not relabel it as life-years, QALYs, or lifetime welfare.

## Required presentation updates

Update the complete documents, including executive summaries, methods, equations, diagrams generated in chunks, figure captions, tables, conclusions, limitations, speaker notes, appendices, and reproducibility text. Static image assets are read-only: if an existing diagram is inconsistent, replace its inclusion with an accurate diagram generated within the Rmd, or explain its limitation; do not alter the asset.

Specifically replace the deck's old `sapm()` direct-mortality summaries and associated interpretation when saved results use progression. Replace old positive-mortality-floor explanations and all-stage-only residual-bound assumptions. Remove claim guards enforcing historical conclusions. Replace the report's “SAP bites on mortality,” “no progression effect,” and “deaths averted, not just delayed” interpretations. Report reductions within the saved horizon without asserting that every avoided death is permanently prevented.

Include dynamic calibration reporting for: prevalence/death fit by age and year; scored versus excluded observations; incidence as validation; target-basis reconciliation; fitted transitions and mortality; stage composition; burn-in/bridge diagnostics; and separate structural, optimizer, and epidemiological-fit statuses. State that a `calibrated` status alone does not certify satisfactory fit. Show usable stage-target counts and partial-identification wording conditionally. Near-optimal optimizer solutions and one-at-a-time flatness are not posterior intervals or joint identification proof.

Preserve the report's useful existing health/economic panels and the deck's fluid narrative, typography, preamble, and main/appendix balance. Put key findings and qualifications in the main slides; detailed equations, fit panels, parameter tables, and reproducibility diagnostics can go in backup. Titles should communicate a current supported message, with compact tables and readable figures. Avoid unrelated styling refactors.

## Validation and final handoff

- Check all code chunks and inline R for syntax, unique chunk labels, missing fields, unsupported assumptions, and stale wording.
- If required data and render dependencies exist, render the country-specific HTML reports and the Beamer deck with all generated files/intermediates directed outside the repository. Inspect the rendered report and slides for broken equations, clipping, missing plots, and readability.
- If rendering is blocked by missing outputs or dependencies, do not alter the pipeline or install dependencies; state the exact blocker and which checks were completed. Still complete all source updates possible within scope.
- Verify that the existing working-tree changes were preserved and that your only source modifications are the two authorized Rmd files.
- Return a concise summary of the methods/narrative updates, the existing artifacts consumed, validation performed, and any unresolved data/provenance problems. Do not claim calibration improved unless the saved diagnostics demonstrate it.
