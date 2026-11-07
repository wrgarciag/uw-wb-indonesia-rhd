# Create a 15-minute RHD conference presentation

Create **`reports/slideck-conference.Rmd`** for a **15-minute conference presentation** on the health and economic consequences of scaling up rheumatic heart disease (RHD) screening and secondary antibiotic prophylaxis in **Indonesia and Uganda**.

**The only permitted project change is creating the new `reports/slideck-conference.Rmd`. Do not modify, overwrite, rename, move, or delete any existing file, script, or folder. Do not create any other project files or folders.** The existing `reports/` folder may receive this one new file only. If the target `.Rmd` already exists, stop and report the conflict instead of overwriting it. This restriction overrides any beautiful-deck instruction to write an outline, install helpers alongside the deck, or save rendered artifacts in the project.

**Report results only through 2050, inclusive. Do not show 2100 or any post-2050 results anywhere in the presentation, including backup slides and speaker notes.**

**Use the beautiful-deck skill, subject to the single-file restriction above.** Read the supplied `beautiful-deck-SKILL(1).md` (or its installed equivalent), the existing `beamer_preamble.tex`, the technical deck, the report, and the relevant code and saved artifacts before drafting. Plan the outline internally or in comments within the new `.Rmd`, then proceed directly through implementation and verification without waiting for outline approval.

## Purpose and framing

Build a conference talk for researchers in global health, health economics, and RHD. Assume an informed audience unfamiliar with this codebase. Establish credibility, explain why the question matters, and devote most of the presentation to the findings and their implications.

**Present the analysis as the current study results. Do not label anything “preliminary,” “provisional,” “draft,” “work in progress,” or “not ready for policy.”** Apply this throughout the conference deck, including the title, subtitle, footer, figures, notes, and backup slides. Do not reproduce the technical deck’s development-status framing, collaborator decisions, discrepancy logs, or outstanding-work lists.

This is a change in presentation and audience, not a change in evidence: describe projections as modelled estimates, state consequential assumptions plainly, and retain scientifically necessary qualifications in concise conference language. Do not claim successful validation, robust sensitivity results, observed intervention effects, or comprehensive costing unless the saved evidence supports those claims. Do not replace an unfavourable or small effect with a more attractive result. A substantive interpretation issue should be explained through its effect on the estimates, without a development-status label.

## Sources and scope

1. Locate and read `CLAUDE.md`/`AGENTS.md`, `reports/slideck-technical.Rmd`, the current report (`report.Rmd`/`report.RMD`), the project bibliography, `scripts/00_run_all.R`, `R/abcd_engine.R`, and the relevant pipeline scripts through `08_economic_evaluation.R`. The supplied copies are `slideck-technical(5).Rmd`, `report(4).RMD`, `beamer_preamble(1).tex`, and `references.bib`; use the repository’s current versions when available.
2. Use the **latest mutually consistent production artifacts for both countries**, with the current specification and saved metadata as the authority. Reuse the technical deck’s portable root detection, loaders, integrity checks, formatters, and result definitions where appropriate. Do not inherit stale hard-coded engine identifiers or parameter values.
3. Read saved model outputs, calibration and baseline bundles, and economic results as needed. The existing deck loads artifacts such as `output/<COUNTRY>/out_model/<COUNTRY>.rds`, `output/<COUNTRY>/tables/rhd_economic_results.rds`, `data/<COUNTRY>/baseline_state.rds`, and `data/<COUNTRY>/calibrated_rhd_parameters.rds`; inspect the actual schemas before using them.
4. **Do not rerun or alter the modelling pipeline, assumptions, calibration, inputs, or production outputs. Treat every existing project file and folder as read-only**, including scripts, reports, decks, bibliography, preambles, configuration, and internal status documentation. Create only `reports/slideck-conference.Rmd`; keep all new loaders, formatting functions, checks, plotting code, and LaTeX additions inside that file. Do not create a separate outline, support script, include file, image directory, or PDF in the project.
5. Rendering must read saved results only. Missing artifacts or incompatible runs must stop with an actionable error. Do not substitute earlier runs, scratch experiments, submitted-abstract numbers, guessed values, zeros, or placeholders. Report any blocking issue separately in your completion message.

## Required reporting horizon: through 2050 only

- Set a deck-specific reporting endpoint of **2050**. Filter all projected time series, figures, tables, and annual results to years at or before 2050. Historical inputs and calibration years may remain earlier than the projection start.
- Calculate cumulative deaths averted, other cumulative health effects, programme costs, monetised benefits, net benefits, and benefit–cost ratios over the **saved analysis start year through 2050 inclusive**, for both countries. Use 2050 for endpoint comparisons and final-year service volumes.
- Do not reuse full-horizon summaries that include years after 2050. Recompute the 2050 summaries within the deck from saved annual outcomes and economic components, retaining the pipeline's valuation, price basis, discount base year, rates, and timing conventions. Calculate ratios from the matching truncated totals; never prorate a 2100 total or average annual ratios.
- Reading annual components and calculating presentation summaries is permitted; rerunning the model or modifying production outputs is not. If the saved components do not support a valid through-2050 economic result, report the missing dependency separately instead of substituting a longer-horizon value.
- Use through-2050 results to choose all assertion titles, the central claim, and the closing message. Check that every displayed horizon and cumulative calculation respects this endpoint; restricting an axis alone is insufficient.

## Story and time budget

Plan the outline first, internally or in comments inside the new `.Rmd`, with one memorable central claim supported by the actual results. For every slide, identify its assertion title, single visual element, exact source object/fields, and speaking time. Do not create a separate outline file. Then implement it. Use approximately **13 main slides including the title and closing**, with no agenda or section-divider slides. The following is a suggested 15-minute arc; adjust the allocation to the evidence while keeping the total at 15 minutes:

| Slide | Role | Minutes |
| --- | --- | ---: |
| 1 | State the study’s central finding in the talk title | 0.5 |
| 2 | Establish the question, country settings, and evidence base | 1.0 |
| 3 | Show the RHD burden or treatment gap that motivates the comparison | 1.0 |
| 4 | Explain the model intuitively with one simplified schematic | 1.5 |
| 5 | Define reference care versus scale-up and the reporting horizon through 2050 | 1.0 |
| 6 | Explain the treatment mechanism and assumptions essential to interpretation | 1.0 |
| 7 | Present the main health effect in both countries | 1.5 |
| 8 | Show the care cascade or service volumes that explain that effect | 1.5 |
| 9 | Present incremental programme costs and their main drivers | 1.5 |
| 10 | Present the principal economic result and the country comparison | 1.0 |
| 11 | Address the strongest scientific objection and its implication | 1.0 |
| 12 | Draw the practical implications supported by the analysis | 1.0 |
| 13 | Close with one earned, memorable sentence | 0.5 |

Avoid a literature review and detailed calibration exposition in the main talk. Put essential supporting methods, fit diagnostics, parameter details, additional economic results, and completed sensitivity analyses in a short, separately numbered backup section. Do not present planned analyses as completed. Include concise speaker notes as R Markdown comments with the intended narration and time per main slide.

## Evidence and scientific wording

- Calculate **every analysis number directly from saved objects**, including numbers in titles, using inline R. Add claim guards or data-dependent titles so the presentation cannot silently contradict a later saved run.
- Clearly distinguish annual from cumulative outcomes, the scale-up endpoint from the reporting horizon ending in 2050, absolute from relative effects, and discounted from undiscounted values. Show currency/base year, valuation perspective, discounting, and denominators where relevant.
- Report deaths averted and other health outcomes only where available. Use the economic measures actually produced by the pipeline; distinguish monetised health benefits, benefit–cost ratios, programme costs, and budget impact. Do not call a benefit–cost ratio an ICER or equate monetised health gains with cash savings.
- Explain secondary prophylaxis, age eligibility, progression effects, and care uptake according to the saved specification. Distinguish assumed effects from trial estimates. Represent the model’s diagnosis and staging logic accurately without importing older model versions.
- Do not invent uncertainty intervals. If results are deterministic, identify them as scenario estimates; include uncertainty or sensitivity results only when corresponding saved analyses exist.
- Use only verified keys in the project bibliography and APA citations following the beautiful-deck workflow. Check the source behind literature claims. Generate references automatically after the talk/backup; do not edit the `.bib`.

## Design and technical requirements

- Use `rmarkdown::beamer_presentation`, XeLaTeX, 16:9, and the **existing project house template and palette**, including its World Bank/UW styling. Reference existing shared template files read-only; place all conference-specific footer text and LaTeX additions inside the new `.Rmd`, never in a new or modified include file.
- Use full-sentence assertion titles that carry the argument when read alone. One idea and one visual element per slide; no generic “Methods,” “Results,” or “Conclusions” headings. Derive the title-slide claim from the results, rather than assuming a favourable conclusion.
- Keep text sparse, usually under 25 words excluding citations and essential labels. Use legible plots, direct labels, consistent country/scenario colours, meaningful ordering, and ample white space. A compact comparison chart may show both countries when it proves one claim.
- Use a simplified model diagram suitable for someone seeing it once; leave stock accounting, optimizer details, code, and equations for backup. Include a concise source line under each figure/table.
- Never shrink text to fit. No `\tiny`, `\scriptsize`, slide shrinking, or scaled-down tables. Close with the take-home sentence, not “Thank you” or “Questions?”.

## Validation and delivery

Run the beautiful-deck source checks and full compile gate with **`--minutes=15`**, using the skill’s gate implementation, **only in a new temporary scratch directory outside the project**. Copy the new `.Rmd` and any required read-only dependencies into that scratch directory for testing; direct all intermediates, caches, figures, gate helpers, logs, PDFs, and page images there. Do not run the gate against the project copy if it creates or deletes files beside the deck. Preserve repository-root resolution and correct read-only artifact paths when rendering the scratch copy. Do not install or update project or system dependencies as part of this task.

Resolve warnings by editing only the new `.Rmd`. Render, rasterise, and inspect **every PDF page** in scratch, including backup and references, for layout, clipping, readable labels, and accurate citations. Verify displayed calculations against the saved objects, confirm that all result series and cumulative calculations end at 2050, and confirm that the title sequence tells one coherent story. Check the final project diff/status: the sole task-related change must be the new `.Rmd`; leave any pre-existing user changes untouched.

Deliver **only `reports/slideck-conference.Rmd`** as the project deliverable. Keep any verification PDF and other test artifacts in external scratch, without copying them back to the project. In your final response, identify the new file, main/backup slide counts, timing, central message, actual verification outcome, and confirmation that no existing project file or folder was modified. If the environment prevents rendering or checks, state that separately and accurately; do not insert a preliminary label into the presentation or claim a clean compile you did not run.
