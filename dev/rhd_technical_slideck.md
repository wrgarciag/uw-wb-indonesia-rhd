# Claude Code prompt: Indonesia RHD technical slide deck

Create `reports/slideck-technical.Rmd`, a Beamer presentation for a **40-minute technical presentation** of the Indonesia rheumatic heart disease (RHD) screening and secondary-prevention investment-case model. The audience is technical collaborators and health-economic/modeling researchers. Explain the model clearly, show the current high-level findings, and make the remaining methodological decisions concrete.

## Scope: only two editable repository files

You may create or modify **only**:

1. `reports/slideck-technical.Rmd`
2. `reports/beamer_preamble.tex`
3. `library/references.bib`

Use the existing preamble and preserve its house style. Make only the minimal, backward-compatible additions necessary for this deck; avoid changing the appearance or compilation behavior of other decks. Do not modify scripts, data, saved model outputs, the report, bibliography, documentation, existing decks, or any other repository file or folder. Do not create outline files, helper scripts, addons files, new bibliography files, or separate diagram assets. Put deck-specific logic in the Rmd, and draw any new schematic inline. Temporary compilation and visual-review artifacts must stay outside the repository. Do not rerun the production pipeline or overwrite its outputs.

These explicit scope instructions override conflicting procedural defaults in the design skill. Proceed through inspection, internal storyboarding, implementation, and verification without an outline-approval pause.

## Inspect the project before drafting

- Read `CLAUDE.md` and applicable repository instructions.
- Read and apply the complete `dev/beautiful-deck-SKILL` to **every slide**, including its standards for narrative, assertion titles, visual hierarchy, source verification, and clean rendering. Resolve the actual filename if it has a `.md` suffix. Use the supplied `beautiful-deck-SKILL.md` if that is the available copy. Adapt its workflow to the two-file limit above; do not install its additional toolchain files.
- Read all relevant `scripts/00-08*.R`, including intermediate calibration scripts such as `04b` where present. Follow calls to existing engine/helper files as needed to understand the implemented model. Do not infer model behavior solely from filenames or comments.
- Read `analysis/Report.Rmd`. If capitalization or location differs, locate the actual existing report (the supplied copy is named `report.RMD`) rather than creating or renaming one.
- Read `reports/beamer_preamble.tex`, inspect existing deck setup conventions, and inspect `library/references.bib`.
- Inspect `docs/rhd_abcd_state_transition.png` and include that existing image in the model-structure presentation.
- Read the supplied `Abstract WCRHD26 DW_NM.docx` if available. It describes hypothetical-country illustrations: **do not present its estimates, populations, prevalence assumptions, targets, or economic findings as current Indonesia results**. Current scripts and persisted Indonesia outputs determine the deck's numerical findings. Flag material discrepancies rather than silently combining different versions. Although in slides report eexplicity ugand and indonesia
- Inspect the saved inputs and result objects actually read by the report, including their metadata, units, scenario definitions, analysis years, and economic assumptions. Read only what is relevant.

## Narrative and slide design

Develop one coherent argument: the policy need → the decision question → how the model represents disease and care → what the comparison changes → projected health and economic consequences → what needs refinement before stronger investment conclusions.

Use roughly **30–35 main slides**, adjusting for complexity, with brief slide-level speaker notes and planned timing totaling 40 minutes. Budget approximately 35 minutes for presentation and 5 minutes for discussion. Appendices and references are outside that speaking budget. Do not pad the talk to meet a slide count.

- Open with the research question, evidence base, and a verified reason the investment decision matters; avoid an agenda or literature-review opening.
- Give every substantive slide a concise **message/assertion title**, with one inferential purpose. Read together, the titles should tell the story. Section labels may organize the source without producing a divider slide for every subsection.
- Introduce technical ideas through intuition, a diagram or example, and then essential mathematics. Keep sufficient modeling detail in the main presentation for technical scrutiny, while reserving full derivations for backup.
- Prefer one legible figure, compact table, equation, or comparison per slide. Split overloaded content. Never shrink fonts to force a fit.
- Reuse the existing World Bank–UW house style and consistent project scenario colors. Use a restrained, consistent color for provisional assumptions. Explain every abbreviation when first introduced.
- Main slides must exclude code listings, chunk names, filesystem paths, object names, software inventories, and reproducibility instructions. These belong in the appendix. Short scientific source citations and units remain on main figures/tables.
- Add concise speaker notes with the intended transition and any explanation that should be spoken instead of printed. Keep notes out of the audience-facing PDF.
- End the main narrative with the conclusion the evidence currently supports and a small set of substantive discussion questions; do not finish with a generic “Questions?” or “Thank you” slide.

## Suggested content and timing

Use the following sequence as an adaptable narrative, not a requirement to create an isolated divider for every item:

| Content | Approximate speaking time | Required emphasis |
| --- | ---: | --- |
| Background, motivation, objectives | 4 min | Indonesia decision context; why screening and secondary prevention merit an investment case; explicit study questions and scope. |
| Model methods and data | 11 min | Process schematic; state structure; population dynamics; inputs and provisional parameters; calibration and validation. |
| Analysis plan and scenarios | 4 min | Comparator, intervention cascade, coverage/adherence, implementation timing, intervention mechanisms. |
| Economic evaluation | 4 min | Costing, benefit valuation, discounting, decision metrics. |
| Current results | 7 min | Validation context; health effects; care/service implications; costs and economic value. |
| Discussion, work in progress, limitations | 5 min | Interpretation; parameter updates; strongest objections; next decisions. |
| Focused questions/discussion | 5 min | Questions that help collaborators resolve the remaining modeling and evidence choices. |

### Methods: explain what the code actually does

Include two distinct visuals:

1. **An inline process schematic** connecting epidemiological and demographic inputs, A–D parameterization/calibration, baseline construction, reference and scale-up simulations, health/service outputs, and economic evaluation. Use TikZ or another existing render-compatible mechanism without writing a new asset. Show logical dependencies, not a row of script filenames.
2. **The existing state-transition image**, `docs/rhd_abcd_state_transition.png`, displayed large enough to read. Explain state definitions and intervention levers from the actual code; identify any discrepancy between image and implementation rather than silently redrawing the disease process.

Explain, where implemented:

- Population scope, age/sex structure, cycle length, horizon, population aging/entry, incidence, background mortality, and RHD mortality.
- Meaning of stages A–D; permitted progression/regression transitions; how the engine handles competing transitions and deaths; relationship between disease stocks, transition flows, and service volumes.
- Screening → diagnosis → optimal treatment/adherence; eligible populations and denominators; distinction between nominal and realized cascade coverage.
- Surgery requirements, coverage, and effects; whether surgery coverage is held equal across the reference and scale-up arms. Verify this from current outputs/configuration.
- Calibration targets, free/fixed parameters, objective/constraints, structural calibration if present, validation evidence, and partial identification. Calibration fit does not establish that each transition parameter is uniquely identified.

Include a compact **data-input table**, split across slides if necessary. Show input/domain, source, model role or scale, and status. Distinguish observed inputs, literature-based values, calibrated quantities, and provisional assumptions. Include transition probabilities (TPs), stage-specific mortality, intervention effects, cascade coverage/adherence, demographics, costs, and valuation inputs as relevant. Put the full parameter table, numerical values/ranges, and source details in backup. Label provisional values clearly; do not imply they are validated estimates.

Use a few carefully selected equations in the main deck: the implemented state update, the intervention's effect on an important transition or mortality mechanism, and the central economic decision metric. Derive the notation from the actual implementation. A generic `n[t+1] = n[t] P[t]` is acceptable only if it represents the engine's ordering and population accounting; otherwise show the relevant entry, aging, competing-risk, or mortality terms. Define every symbol, scale, and time unit. Reserve complete matrices and derivations for backup.

### Scenarios and economic analysis

Distinguish the implemented analysis from planned extensions. The report currently describes `ref` as holding the care cascade at baseline and `sap` as scaling screening, diagnosis, and prophylaxis; confirm this against the current code and artifacts. Take targets, timing, screening ages, and adherence definitions from the active configuration and saved metadata. Do not import scenarios from unrelated projects.

Explain the economic perspective, currency and price year, costing components, unit-cost denominators, incremental costs, benefit valuation method, discount rates and base year. Show benefit-cost ratio and net economic benefit if implemented. Do not equate benefit-cost analysis with cost-effectiveness analysis, or monetize a health benefit twice. Only present cost-per-health-outcome metrics if supported by the current outputs; label exploratory deck-only calculations transparently.

### Results, interpretation, and ongoing updates

Select the strongest available results rather than copying the whole report. Candidate views include baseline fit, annual/cumulative RHD deaths averted, changes in disease-stage composition, service volumes, incremental costs/budget impact, and discounted economic benefits. Separate counts from rates, annual from cumulative outcomes, and discounted from undiscounted quantities. Make the comparator and time horizon explicit.

Describe deterministic results as point estimates. Present PSA/DSA, uncertainty intervals, or robustness findings only if they actually exist; otherwise put uncertainty analysis in the work plan. Explain model-dependent interpretation and avoid claiming empirical causal identification.

The work-in-progress material should distinguish **implemented**, **provisional**, and **planned**. Identify evidenced needs such as updating progression/regression TPs, stage-D mortality, intervention assumptions, calibration targets, costs, and structural assumptions only when supported by project materials. Do not implement those updates as part of creating this deck. Explain which outputs or conclusions each unresolved issue could affect and what evidence or decision is needed next.

Address the strongest limitations with their consequences: identifiability, transferability of parameters to Indonesia, aggregate calibration versus stage-specific natural history, mortality attribution, implementation/adherence, and economic valuation where applicable. Finish with specific questions collaborators can answer.

## Data provenance and references

- **Use only `library/references.bib` for citations and the generated bibliography.** Verify every citation key exists; never invent a key, source, author-year attribution, or reference entry. Do not modify the bibliography or add another one.
- Verify what a source supports using accessible project source text. A bibliography title alone does not verify an effect size or numerical claim. Omit unsupported claims or visibly identify the evidence gap. Do not browse for additional references.
- Use a compatible citation configuration with APA formatting where supported, without loading conflicting citation packages. Generate references from cited entries only; do not handwrite a reference list or cite the entire bibliography.
- Read persisted results; never source scripts that run or mutate the model. Reuse report calculations and aggregation definitions in read-only deck code as needed. Adapt HTML widgets/tables to static Beamer-compatible plots and tables.
- All reported model estimates, including quantitative titles, must be computed from current saved objects through R. Add lightweight checks for units, country, scenario pairing, years, finite values, and claims made in result titles. Do not hard-code abstract results or copy numbers manually.
- Use Indonesia-scoped paths and resolve the repository root portably. Do not carry over the report's machine-specific absolute fallback path.
- If a required artifact is missing, stop with an informative error identifying its producing script. Do not rerun the pipeline, substitute another country's results, fabricate numbers, or silently omit a core result. Optional absent analyses may be explicitly described as not yet available.

## Appendix and verification

After the main conclusion/discussion, add separately numbered backup slides and then references. Backup should include detailed model mathematics, full input/TP assumptions, calibration and identifiability diagnostics, detailed health/economic outputs, available sensitivity analyses, and reproducibility. In the appendix only, provide source script, artifact path, object/component, dimensions, units, and transformations for the major displayed results. Do not create a separate provenance file.

Render with `rmarkdown::beamer_presentation` and XeLaTeX, including `reports/beamer_preamble.tex` through the correct relative path. Keep R warnings visible to verification rather than suppressing them globally. If an existing compatible gate is available, run it without modifying or installing files. Otherwise perform equivalent citation, source, render-log, and visual checks within the permitted scope. Render to a temporary directory outside the repository and inspect every PDF page for clipping, overflow, unreadable text, broken citations, image quality, and narrative flow. Resolve issues within the two allowed source files.

Do not claim a clean render if dependencies or artifacts prevent it. Report the exact blocker and what was verified. Check the final repository diff: only the two authorized files may have changed; preserve all pre-existing unrelated changes.

Complete the work, then give a concise final response identifying the Rmd, main/backup slide counts, timing, preamble changes, verification status, and remaining evidence gaps. Do not create any additional repository deliverable.
