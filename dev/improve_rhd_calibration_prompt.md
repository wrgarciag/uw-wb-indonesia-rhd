# Claude Code prompt: improve RHD calibration while preserving downstream outputs

You are working in the existing RHD investment-case repository. Implement focused improvements to calibration and the minimum supporting input, engine, and baseline code needed for the calibrated parameters to feed the production model correctly. Explicitly implement both exact D-only RHD mortality and an evidence-supported progression-based SAP mechanism in the allowed files. These are complementary requirements: SAP changes supported progression pathways, and its mortality benefit arises through the resulting stage trajectories and reduced entry into D.

Work through implementation and appropriate validation. Do not stop at an action plan.

## 1. Goal and scope

The final pipeline must remain:

prepared inputs -> calibration -> baseline state -> production model -> existing output tables -> existing report and technical slide deck.

The downstream consumers are:

- scripts/07_make_outputs.R
- scripts/08_economic_evaluation.R
- analysis/report.Rmd
- reports/slideck-technical.Rmd

Resolve their actual spelling and capitalization in the repository before reading them. Do not create duplicates if, for example, the report is analysis/Report.Rmd.

Strict scope: modify only the minimum necessary calibration-related lines in the seven existing files below. This is an allowlist, not a requirement to edit every listed file. Leave any listed file unchanged if the improvement does not require an edit there.

| Allowed file | Authorized purpose |
| --- | --- |
| scripts/00_run_all.R | Calibration settings and structural options only; leave unrelated control-panel blocks unchanged |
| scripts/03_build_disease_model.R | Necessary parameter and input definitions supporting calibration and the revised dynamics |
| scripts/04_calibration_random_tp.R | Target construction, calibration input preparation, and the existing joint-calibration interface |
| scripts/04b_calibrate_structural.R | Joint incidence/natural-history calibration, objective, and calibration diagnostics |
| scripts/05_build_baseline.R | Transfer final calibrated parameters and initialize a consistent model state |
| scripts/06_run_prevention_model.R | Minimum changes needed for consistent calibration/projection dynamics and output compatibility |
| R/abcd_engine.R | Shared disease-cycle/projection helpers needed by calibration and production |

Do not modify scripts/01_prepare_inputs.R or scripts/02_build_demography.R. Read them only as needed to understand their persisted outputs. Keep the existing GBD preparation and WPP population-backbone implementations unchanged. Implement calibration-specific aggregation, validation, and transformations in 03/04/04b/05 as appropriate. If an upstream defect cannot be addressed correctly within this scope, document the exact blocker without editing the protected script or concealing the problem.

All other source files are read-only, including scripts/07_make_outputs.R, scripts/08_economic_evaluation.R, R/packages.R, R/utils.R, report/deck sources, configuration, and dependency files. Do not add, delete, rename, or relocate source files. Put any essential new shared functions in R/abcd_engine.R rather than creating a new helper module. Temporary verification files and generated model/diagnostic outputs are permitted; they do not authorize new source modules or edits outside the allowlist.

Inspect the actual repository before assuming these exact files or functions exist. If an allowed filename is absent, report the mismatch rather than silently expanding the edit scope. Preserve sequential script entry points.

Read scripts/07_make_outputs.R, scripts/08_economic_evaluation.R, analysis/report.Rmd, and reports/slideck-technical.Rmd to establish the output contract. Leave these downstream consumers unchanged. Do not change economic assumptions, costs, report layout, slide narrative, or unrelated intervention features. The D-only and SAP progression/adherence changes specified below are explicitly in scope; do not interpret the restriction on unrelated interventions as permission to defer those requirements.

Do not reorganize directories, rename existing outputs/functions for style, rewrite working modules, upgrade dependencies, change the package stack, or perform repository-wide formatting. Use data.table and the existing project conventions.

Respect existing user edits. Inspect git status and applicable AGENTS.md instructions before editing.

## 2. Evidence and scientific direction

The agreed calibration direction is:

1. Make RHD prevalence, including its age pattern, the principal calibration target.
2. Fit RHD-attributable deaths jointly while protecting prevalence fit.
3. Remove incidence from the direct primary loss, but continue estimating and reporting incidence.
4. Allow incidence corrections to move within documented plausible constraints rather than freezing them after an incompatible aggregate calibration.
5. Implement exact D-only RHD mortality and the progression-based SAP mechanism described below; recalibrate under each selected mechanism instead of transferring parameters from an incompatible mechanism.
6. Diagnose conflicting inputs and observation definitions before interpreting poor fit as an optimization failure.

Use a prevalence-to-death weight ratio of 2.5:1 as the primary configurable setting, with 2:1 and 3:1 as sensitivity settings. These are operational choices for this implementation, not weights established by the meeting or a published likelihood.

Important qualifications:

- GBD prevalence is modeled evidence, not a direct national survey at every age/year.
- Evidence for young-age echo prevalence may be stronger than for older-age prevalence. Permit documented age/evidence weights, but do not invent country-specific quality scores or silently exclude older ages.
- A good incidence fit does not validate disease duration or stage progression.
- Fitting total prevalence and deaths does not uniquely identify every stage transition.
- An increase in prevalence after effective treatment can reflect improved survival; do not force intervention prevalence to decrease.

Read relevant local methods documents and references, especially the Coates et al. 2021 paper and supplementary appendix if present. Search likely repository locations by filename rather than assuming an attachment path. Use library/references.bib and local source documents; do not invent references or claim unread documents were checked.

### Coates precedent: adapt the reasoning, not the numerical constants

The supplied Coates appendix describes:

- Inconsistency between GBD 2017 prevalence, severe disease/heart failure, and deaths in Africa.
- A primary hybrid mortality pattern: DisMod age pattern below 45, constrained by the cardiovascular-death envelope at ages 5-19, and published CODEm RHD mortality above 45. This produced about 35,000 versus about 18,000 published RHD deaths in the AU.
- Untreated severe-RHD annual mortality assumed at 25%, with uncertainty roughly 15-40%.
- Back-calculated progression to severe disease, capped at 10% annually.
- Alternative published-GBD and Heart-of-Soweto-informed specifications.
- Reconstruction of starting populations using prevalence, survival, and incident flows.
- Remission assumptions informed mainly by evidence under age 20.
- A 55% secondary-prophylaxis effect on progression pathways, not a uniform direct mortality reduction in every stage.

These are historical, setting-specific assumptions. Do not automatically apply an African mortality multiplier, 25% stage-D fatality, a 10% progression cap, or the male incidence scalar to Indonesia or Uganda. Coates severe RHD/heart failure is not automatically identical to WHF stage D.

Keep unmodified published GBD 2023 targets as the default. Support an optional explicit target-override table only if it can be implemented simply. Require source, rationale, age/sex/year identifiers, and provenance; save original and adjusted targets side by side. Never let the optimizer change its own targets, construct undocumented mortality targets, or infer unavailable DisMod inputs. When no justified overrides exist, run the published-target analysis and report any remaining incompatibility.

## 3. Inspect and freeze the existing output contract first

Before changing code, inspect all readRDS/load/fread calls and named object/column accesses in the report, deck, scripts 05-08, and existing tests.

Record the required:

- Filenames, country-scoped paths, and file formats.
- List names, nested fields, table columns, key types, dimensions, dimnames, scenarios, and location/sex/age/year identifiers.
- Meaning and timing of legacy fields such as sick, dead, newcases, all.mx, seed, IR, and CF.

In the uploaded version, key persisted interfaces included:

- disease_model_inputs.rds
- calibrated_rhd_parameters.rds with tp, factors, diagnostics, stage_calibration, layer1, and meta
- stage_calibration$best_parameters, including transitions, p_rhd_death, stage_split, rhd_d_fraction, and stage_shares_by_age
- baseline_state.rds with locations, states, and meta
- Production location output with wsd, stages, diag, and meta
- Existing tables consumed by the two R Markdown files

Verify against the actual repository; this list is not exhaustive.

Preserve required legacy structures and their documented meanings. Append diagnostics through existing extension points or separate new calibration diagnostic files. Do not delete fields, silently repurpose them, rename scenarios, or require report/deck edits.

If an aggregate proxy or legacy CF is retained for compatibility, label its role truthfully in additive metadata. Do not claim it anchors the stage objective unless it actually does. Update stale explanatory comments only where touched by this work.

Do not satisfy compatibility by placing new, semantically different values under old field names without documentation.

## 4. Correct target construction and rate semantics

Audit and fix the following where present:

### Age and population alignment

- Use the same internal single-age grid for calibration and production, preferably 0-100 with an open 100+ group, matching the current population backbone.
- Aggregate all ages 95 and older when comparing with a GBD 95+ target. Never use the population of age 95 alone as the denominator for 95+.
- Retain existing reporting age groups. Distinguish internal age resolution from observation aggregation.
- Compare GBD Number with GBD Rate multiplied by the chosen population for matched age-sex-year cells.
- Keep population-source differences separate from model-dynamics residuals. Do not silently mix a WPP denominator and a GBD count target without a documented reconciliation.
- Audit the effect of any existing population smoothing on age-group denominators. Do not introduce additional smoothing as a calibration fix.
- Verify persons versus thousands and rates per 100,000 versus fractions.

### Completeness and definitions

- Assert unique input keys and complete required grids before widening or array assignment.
- Distinguish missing observations from genuine zero values. Missing target cells must be masked, not silently filled as zero.
- Treat missing required model inputs as errors with cell identifiers unless an existing documented imputation applies.
- Explicitly check both sexes and all needed years.
- Document the target case definition and its mapping to A/B/C/D. Do not assume GBD prevalence includes every early/borderline state without checking.
- If a definitive case-definition crosswalk is unavailable, retain a clearly labeled provisional mapping and expose the issue in diagnostics. Do not fabricate a conversion.

### Rates and probabilities

- Define units for each input: rate/hazard, net annual event probability, or observed competing-risk transition probability.
- Convert central mortality rates to hazards/probabilities consistently; passing a rate directly as a probability and then applying -log(1-p) is not a valid general conversion.
- Ensure calibration and production use identical incidence denominator and competition conventions.
- If GBD incidence is per total population, explicitly map it to the susceptible-person hazard needed by the engine.
- Preserve raw reference incidence and report realized incident counts/rates.
- Do not double-convert literature probabilities or treat observed competing transition probabilities as independent net risks without justification.

## 5. Use the production disease dynamics for calibration

The uploaded code calibrated incidence/CF in an aggregate well-sick-dead model without remission, then held incidence fixed while fitting a stage model with remission. Remove this inconsistency.

Use the actual A/B/C/D engine and the same essential projection wrapper in the objective and production. Share initialization/ageing/observation logic through minimal functions in the existing R/abcd_engine.R where needed; do not create another source file.

Preferred implementation:

- Keep scripts/04_calibration_random_tp.R as the existing entry point for assembling targets, inputs, starting values, and any retained proxy diagnostics.
- Let scripts/04b_calibrate_structural.R jointly estimate low-dimensional incidence corrections and the selected natural-history parameters against prevalence and RHD deaths.
- Use the aggregate model only for starting values or a labeled comparison, not to freeze incompatible incidence parameters.
- Persist the final incidence corrections and stage parameters consistently in the existing bundle.
- Apply each calibrated factor exactly once. Ensure 05 and 06 consume the final values without recalibration or duplicate scaling.

Reuse broad age bands if appropriate; avoid introducing one free parameter per age-year. Retain smoothness and documented weak incidence priors. Setting incidence loss to zero must not leave incidence effectively fixed through tight priors, bounds, or an unmodifiable upstream table.

Use a configurable direct incidence weight of zero in the primary objective. Incidence remains a validation output. A weak nonzero incidence-loss specification may be retained as an explicitly labeled sensitivity analysis.

## 6. Implement exact D-only mortality and avoid redundant parameters

Implement a configurable mortality specification in 00, 03, R/abcd_engine.R, 04b, and 05, with the minimum compatible propagation in 06:

- The revised primary structural experiment must impose p_rhd_death_A = p_rhd_death_B = p_rhd_death_C = 0 exactly. Retain competing non-RHD mortality in every living state, including No-RHD.
- Retain the previous all-stage pathway only as a labeled sensitivity/comparison option. Do not silently restore it if the D-only model fails to fit.
- Treat exclusive-D RHD mortality as an explicit modeling assumption requiring clinical review, not a fact established by Coates or WHF.
- Branch structural parameter packing/unpacking, priors, starting values, and validation by mortality structure. Remove the A-C mortality-ratio parameters and their log-prior terms from the D-only branch. Do not evaluate log(0) or approximate structural zeros with tiny positive values.
- Replace strict B<C<D mortality ordering checks with branch-appropriate validation: exact A=B=C=0 and valid D mortality in D-only mode. Keep appropriate legacy checks in the legacy branch.
- Prefer one D-mortality scale. Fix the residual mortality multiplier to 1 in the initial D-only specification rather than estimating two interchangeable scale parameters. Preserve required legacy fields as truthfully documented fixed/derived fields.
- Review D-mortality bounds against the selected estimand and clinical source. Do not truncate workbook estimates to existing bounds or widen bounds without documenting why. Do not replace country evidence automatically with Coates's 25%.
- If the selected evidence provides age-specific D mortality, support it through parameter inputs, calibration, serialization, baseline transfer, and production using a minimally invasive scalar/matrix representation. Validate dimensions, dimnames, country/sex handling, and age coverage. Do not invent sex differences when evidence is pooled.
- Keep untreated/reference-care D parameters separate from the realized competing-risk death probability. Export rhd_deaths_D / D_start as a diagnostic, using start-of-cycle stock, not end stock. For D_start=0, use an explicitly documented undefined-value convention and flag, with zero deaths required; keep undefined values outside legacy fields that require finite values.
- Preserve mathematical equivalence between full and light=TRUE shared-engine calculations for the stocks and deaths used in calibration.

Inspect TP Stage D to death_RHDparameters.xlsx if locally available. Column K is the previously requested provisional scenario, not a validated annual transition probability. Preserve numerator, denominator definition, age strata, year, and alternative severity/heart-failure definitions. Calibrate D stock and deaths consistently; do not mix a D population from one definition with fatality derived from another. If these workbook inputs are unavailable, expose a clear missing-evidence status and keep any existing-D-prior run explicitly exploratory, never labeled column-K calibrated.

Use transition-specific bounds/priors justified by existing evidence. Do not retain a universal 0.0001-0.60 range as though it were a literature-based interval for every transition.

Where justified, permit parsimonious age variation in progression/remission. Do not force high adult remission to absorb incompatible incidence or initialize the model. Make any age restrictions explicit and sensitivity-test them.

## 7. Make the objective interpretable

Use a transparent primary objective:

L = wP * normalized prevalence error
  + wD * normalized RHD-death error
  + prior penalty
  + smoothness penalty
  + optional observed stage-target error.

Requirements:

- Set wP/wD = 2.5 primarily; expose 2 and 3 as controlled sensitivity settings.
- Use comparable, explicitly normalized errors; weighted mean squared log error is acceptable.
- Define how zeros, small counts, missing cells, sex weights, age weights, and year weights enter.
- If uncertainty-based weights are used, implement them consistently in the structural objective, with a variance floor or weight cap.
- Do not interpret GBD uncertainty intervals as independent measurement errors across every age/year.
- Keep the overall data-weight scale controlled across weight-ratio comparisons so prior strength does not change accidentally.
- Export unweighted and weighted losses for each target and every prior/smoothness contribution.
- Assess whether summed prior penalties overwhelm averaged target errors. Calibrate penalty strength transparently; do not indiscriminately discard clinical evidence to force fit.
- Report age-specific fit alongside total fit so national totals cannot hide a poor age pattern.

A prevalence-only preliminary optimization may supply starting values. The final fit must evaluate prevalence and deaths together.

## 8. Align historical calibration, initialization, and future projection

Use one consistent time convention for start stock, within-year flows, end stock, and ageing.

- Preserve legacy output-field meanings; add explicit start/end/average prevalence diagnostics as needed.
- Compare GBD prevalence with the corresponding modeled observation, not a differently timed stock.
- Distinguish point prevalence, average annual prevalence/person-time, and the number of distinct people occupying a state during a year.
- Use a historical initial state or a documented initialization strategy informed by available observations.
- If retaining burn-in, assess convergence by age/stage and compare a longer burn-in. Do not label a fixed 60-cycle zero-start run as equilibrium without checking.
- Include births/new entrants and the handling of migration/population reconciliation explicitly.
- Avoid disease stock exceeding the age-sex population. Surface violations instead of hiding them with pmax or pmin.
- Forward-project the calibrated historical state to the analysis start through the intervening years with consistent rates and care assumptions.
- Avoid combining 2019 stage shares, a fresh 2023 prevalence seed, and 2026 population without an explicit reconciliation.
- Anchor incidence trends to an explicit calendar year. Changing the analysis start must not accidentally redefine the incidence trajectory.
- Do not introduce a new demographic model or modify scripts/02_build_demography.R. Consume its existing population outputs and make only necessary calibration/model-side alignment corrections within the allowed files. Document the population-backbone approximation.

Historical calibration must represent existing care. If the parameters are untreated natural history, apply estimated historical/reference treatment and surgery within calibration. If parameters are effective under existing care, do not apply the same baseline benefit again. Save the chosen interpretation and assumptions.

## 9. Implement progression-based SAP and coherent baseline care

SAP progression effects are part of this implementation because reference treatment influences calibration and future treatment changes the stage distribution. Implement the selectable mechanism and adherence interface in the allowed files; do not reduce this section to a recommendation for future work.

### Mechanism and effect scale

- Add transition-specific SAP effects only for supported pathways. A-to-B, B-to-C, and C-to-D are candidate pathways, not permission to apply one effect to all three. Do not assign a combined progression endpoint independently to several transitions.
- Verify source definitions in the available literature/progression workbook. Coates's 55% prophylaxis parameter is not evidence for a uniform direct mortality reduction in A-D, nor automatically a measured effect for every WHF stage transition.
- Record whether each effect is a risk ratio, hazard ratio, or another measure. Apply it on the corresponding scale or document the approximation. If stage mapping requires extrapolation, label and parameterize the assumption explicitly.
- For an assumed hazard-ratio effect, an annual population-average formulation is h_ij(c) = h_ij(0) * [1 - c * (1 - HR_ij)], where c is effective protection in the eligible source-stage population. Convert hazards into engine probabilities consistently afterward. This is an explicit modeling approximation, not automatically an exact mixture of annual risks.
- If fitting natural history under reference effective coverage c0, either estimate untreated parameters within calibration or normalize by the reference factor: h_ij(c)/h_ij(c0) = [1-c*(1-HR_ij)]/[1-c0*(1-HR_ij)]. Validate the denominator and do not apply baseline care twice. Use the corresponding convention when evidence is on a different scale.
- In progression mode, disable the existing direct SAP mortality multiplier. Keep the former direct-mortality mechanism as an explicitly selected legacy/provisional comparator only. Never retain the full 55% direct mortality effect automatically when adding progression benefits.
- A combined progression-plus-direct-mortality mechanism requires separately documented evidence for a residual direct benefit. Keep any HF-management mortality effect distinct from SAP and avoid attributing the same evidence twice.

### Cascade, adherence, and treatment volumes

- Preserve existing ref/sap interfaces. Define screening, diagnosis, SAP initiation, and sustained effective adherence denominators explicitly. pmin(screen, diagnosis, treatment) is not by itself a coherent cascade if screening is a fraction of total population and diagnosis/treatment are fractions of RHD cases.
- Multiply only genuinely conditional cascade proportions; do not multiply already cumulative coverage measures. Use the existing aggregate cascade where possible and clarify its interpretation instead of building a new screening/detection model.
- Add an explicit annual population-average adherence modifier defaulting to 1 in 00 and propagate it through 03/05/06 and calibration. Default 1 means no additional adherence attenuation of the chosen effect; it does not imply 100% treatment coverage. Support country-specific values and time trajectories through existing configuration structures where evidence is available; do not invent values.
- If existing optimal-treatment coverage already includes adherence, do not apply another absolute adherence probability. Document whether the new modifier is a relative scenario attenuation or a separate adherence component and preserve equivalence at modifier=1 within the selected mechanism. This does not require reproducing the old mortality-only results after changing mechanism.
- Separate delivered SAP/person-years from effective protection in the model's diagnostics. Incomplete adherence may still incur delivery costs. State whether optimal treatment includes services beyond antibiotics.
- Preserve the meaning of every legacy volume consumed by unchanged scripts 07/08. Do not silently repurpose n_on_optimal_treatment as delivered SAP or claim that economics costs imperfectly adherent treatment correctly if its current consumer does not. Append delivered/effectively protected volumes through existing diagnostic extension points or separate calibration diagnostic outputs. Identify any economic integration requiring a later change to protected scripts as outside this task.
- Keep annual population-average adherence for the first implementation. Do not claim to track person-level persistence or interruptions and do not add treatment-history states without a demonstrated need.

### Evidence gaps and exclusions

If pathway-specific evidence is unavailable or cannot be mapped defensibly, implement the mechanism interface, required-input validation, and clearly labeled synthetic verification. Complete independent D-only/calibration repairs. Mark the affected production progression run as blocked or exploratory with the exact missing inputs; do not silently fall back to direct mortality or invent efficacy. Use the existing clinically supported prior only with its actual source and provisional status.

Keep surgery need, coverage, effects, and service formulation fixed across the first comparisons. Do not add an incidence effect, postoperative state, adult incident-entry allocation, acute-carditis fatal pathway, new screening model, or economic model. Preserve incident entry through A as an explicit assumption; provide a sensitivity hook only when defensible inputs exist. Record stage-A age interpretation and possible fatal exceptions as clinical questions, not automatic new pathways.

Use the selected SAP mechanism, baseline care, and adherence consistently in both calibration and projection. A calibration fitted under one mechanism must not be consumed under another without recalibration or an explicitly justified conversion and compatibility check.

## 10. Fix stage-target handling and calibration status

The uploaded version detected a local echo-target file and could label the result calibrated_with_echo_targets without reading those targets into the objective.

Fix this:

- Read and validate any supplied stage-target table.
- Map stage definitions, ages, sex, year, denominators, and units explicitly.
- Evaluate model observations matching the target definition.
- Include the resulting loss when usable targets exist.
- Distinguish detected, usable, and actually used targets in metadata.
- Never label a fit as using echo targets if it did not.
- Do not fabricate local stage targets or stage distributions.

If only total prevalence and deaths are available, report stage parameters as partially identified/prior-informed. Do not label narrow one-parameter sweeps of a penalized objective as proof of data identification.

Save complete parameter sets for accepted solutions, including incidence and all fitted transitions. Optimizer endpoints alone are not a posterior or a complete uncertainty region. Preserve reproducible seeds, convergence codes, boundary diagnostics, and the evaluated baseline solution.

## 11. Diagnostics and validation

Run existing tests where available without modifying their source files. Perform any additional focused verification in temporary scratch files or nonpersistent R commands; do not add a new script or test module to the repository. Add runtime validation inside the seven allowed files only where necessary for the revised pipeline.

Required checks:

1. Complete/unique input grids; correct 95+ aggregation; unit and rate conversions.
2. Exact per-stage event accounting and finite, nonnegative flows.
3. RHD stock-flow identity before ageing:
   S_end = S_start + incident RHD - remission - RHD deaths - other deaths among RHD.
4. Correct full/fast shared-engine equivalence for the outputs used in calibration.
5. Identical calibration/production results under identical inputs, initial state, care, time convention, and horizon.
6. Exact zero RHD deaths in A-C under the D-only specification, with background deaths still present; total RHD deaths equal D deaths; D fatality uses D_start with an explicit zero-denominator convention.
7. Incidence multipliers applied once; removing direct incidence loss really allows incidence to move.
8. Consistent baseline-care interpretation; adherence modifier=1 preserves the selected mechanism's unattenuated effect, and zero efficacy or identical reference/scale-up protection gives identical health outcomes. At fixed start stocks, greater effective protection lowers the supported progression risks without activating a direct SAP mortality effect in progression mode. Verify that adherence is not counted twice and delivered/protected volumes preserve legacy consumer semantics.
9. Actual use of supplied stage targets, including a small synthetic target check.
10. Preservation of output filenames, required fields/columns, types, dimensions, identifiers, and downstream read access.

Reconstruct all-cause deaths from the modeled competing-risk flows, including those among people without RHD. Do not retain rhd_deaths + oth * population as a flow identity if it disagrees with the actual competing-risk accounting.

For data-enabled fits, export:

- Observed versus modeled prevalence, RHD deaths, and incidence by age, sex, and year, in counts and compatible rates.
- National fit summaries and loss decomposition.
- Stage stocks/shares, progression/regression flows, and deaths by stage.
- D_start, realized D death risk and its defined/undefined flag, mortality structure, SAP mechanism, effective adherence/protection, delivered SAP volumes, and pathway-specific effect modifiers in additive diagnostics.
- Parameter values, sources/priors/bounds, boundary hits, optimizer convergence.
- Initialization and age-tail sensitivity.
- Original versus optional adjusted targets.
- Baseline-to-analysis-start reconciliation and held-out fit where feasible.

Label prevalence/incidence as a cross-sectional ratio, not mean disease duration. Derive disease duration from cohort trajectories if needed; this is optional.

Use separate structural, numerical-convergence, and epidemiological-fit status fields. Define transparent, configurable fit tolerances; do not invent a universal acceptance threshold or certify poor fits because they pass probability bounds.

Only finalize a successful calibration bundle after validation. A failed candidate must not overwrite a valid production bundle or remain labeled as calibrated. Keep failed diagnostics separately and record the failure reason.

## 12. Execution and comparison

Before editing, give a short list of the specific allowed files/sections you will change, explain why each is necessary, and identify the compatibility contracts you found. Then implement without waiting for routine confirmation. Do not broaden the edit scope to address unrelated findings.

Run a controlled comparison where inputs/runtime permit:

- Existing specification captured before edits.
- Prevalence-led loss with direct incidence weight zero.
- Revised shared-engine joint calibration.
- D-only versus previous mortality structure using otherwise matched assumptions.
- D-only with the legacy direct-mortality mechanism versus D-only with the supported progression mechanism, using mechanism-consistent recalibration and otherwise matched assumptions. Label missing-evidence comparisons as not run rather than presenting synthetic results as empirical.

The revised mechanistic specification is D-only RHD mortality plus supported progression-based SAP. The provisional direct-mortality comparator remains available for attribution, not as an undocumented fallback. A failed or evidence-blocked mechanistic run must not overwrite the last valid production result.

Avoid a full factorial sensitivity exercise unless necessary to resolve a finding. Do not change several unrelated assumptions and then attribute the improvement to one weight.

Use synthetic checks when data are unavailable. If actual calibration cannot run, state exactly which inputs or packages are missing and which checks were completed. Do not claim improved empirical fit, successful report rendering, or identified parameters without evidence.

Check that the existing report and technical deck can still access their required outputs. Render them if the data and toolchain are available, without editing their sources. If rendering cannot run, verify the consumer contract and report the limitation.

## 13. Completion response

Provide:

1. Files changed and the calibration-related reason for each.
2. Final objective, exact D-only implementation, SAP progression pathways/effect scales, adherence convention, incidence parameterization, and baseline-care interpretation.
3. Compatibility checks showing how calibrated outputs still feed 05 -> 06 -> existing downstream files.
4. Validation/comparison results, including remaining age-specific discrepancies and prior dependence.
5. Scientific assumptions or missing evidence requiring review, including any unsupported stage mapping and any cost-volume integration deferred because scripts 07/08 are protected.
6. Exact commands to reproduce calibration and downstream output generation.

Before finishing, inspect the final diff and verify that every source change is confined to the seven-file allowlist, scripts 01 and 02 remain unchanged, and no unnecessary scripts or modules were modified or created. Distinguish pre-existing user changes from your own edits; do not revert user work. Report this scope check in the completion response.

Keep the diff focused. Do not claim that model fit improved unless the actual fitted outputs support that statement.

