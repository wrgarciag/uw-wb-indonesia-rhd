# Claude Code prompt: 15-minute conference presentation on RHD SAP value of information

## Goal

Produce and render a standalone, 15-minute conference Beamer presentation about the value of better evidence on secondary antibiotic prophylaxis (SAP) for rheumatic heart disease. Create **`reports/slideck-voi.Rmd`**, using **`reports/beamer_template.tex`** and the supplied **beautiful-deck** skill.

The talk should communicate the decision, the findings, and their practical meaning. It is a separate presentation from the modelling-pipeline deck: include only enough disease-model context to understand the VOI. Put technical detail in backup slides.

Use **2026–2050 inclusive** for all headline results, BCR calculations, VOI estimates, rankings, and normalization. Do not relabel 2026–2100 totals as 2026–2050.

## 1. Strict edit scope and workflow

Only these two implementation files may be created or modified, plus the append-only bibliography exception below:

1. Modify `scripts/voi/04_analyse_and_plot.R`.
2. Create `reports/slideck-voi.Rmd`.

Read any necessary project files, but leave the engine, scripts 00–03, production scripts, inputs, existing decks, template, preamble, skill, and other source files untouched. Preserve unrelated working-tree changes. Generated analysis tables, figures, provenance records, render intermediates, and the presentation PDF are allowed in a dedicated VOI conference-output directory; they are outputs, not additional implementation files.

Do not rerun calibration, production, or PSA. Reuse the existing paired draws and saved annual arm-level outputs. Keep the initial 500-draw design; do not silently extend it.

Read the beautiful-deck skill in full and follow its conference storytelling and design standards. The explicit constraints in this request override its defaults: do not create outline, addon, or gate source files; keep the storyboard and speaker notes as non-rendered comments inside the new Rmd. Use an existing compatible gate if available, otherwise perform equivalent checks without adding another project source file. Proceed through rendering without an outline-approval pause.

**Bibliography exception:** it is fine to add suggested additional references directly to `library/references.bib` when they support the presentation, VOI methods, SAP evidence, or normalization data. This is an explicitly authorized exception to the two-file implementation limit and to the beautiful-deck skill's default prohibition on editing the bibliography. Reuse existing entries first; verify new references against their original sources and append only missing entries with unique citation keys and accurate metadata. Preserve all existing entries and formatting, avoid duplicates, and cite the added entries from the deck where relevant. Do not create a temporary combined bibliography or leave needed references as proposed additions in comments. Report the citation keys added in the completion report. Never invent citation metadata.

## 2. Inspect the current implementation and available evidence

Read:

- Project `CLAUDE.md` / `AGENTS.md` instructions.
- All current `scripts/voi/` scripts and helpers, including the isolated engine, to establish what is actually estimated.
- The current run configuration, manifests, prior specification, validation checks, draw statuses, annual summaries, saved draw-level outputs, and prior/precision diagnostics for Indonesia and Uganda.
- `reports/beamer_template.tex`, any associated preamble, and an existing deck setup for the correct rendering convention. Determine whether the template is a full Pandoc template or a header include; use it accordingly.
- `library/references.bib`, and the actual sources supporting the methods and SAP assumptions.

Use current saved objects as authoritative. Do not infer numerical results from an earlier deck, this prompt, or a manuscript. Identify the run explicitly; do not accidentally combine different runs or priors.

## 3. Recompute the conference analysis for 2026–2050

Implement a self-contained conference analysis path inside script 04. Put helpers and configuration there rather than editing upstream scripts. Prefer **R data.table** for loading, joins, filtering, and summaries.

Use the saved annual paired arm-level results, filtered to 2026–2050, and join the corresponding draw parameters by country and draw ID. Recalculate:

- Discounted incremental cost, monetized mortality benefits, and incremental net benefit.
- Deaths averted and any reported DALY proxy.
- Discounted BCR.
- Decision probabilities, joint SAP EVPPI, individual EVPPIs, bootstrap intervals, and rankings.
- Relevant convergence and estimator diagnostics for this horizon.

Preserve the applied price year, VSL schedule, and discount base/rate from the saved economics. Do not repeat price conversion on costs already harmonized. All monetary numerators and denominators must be comparable. State that these are benefits and costs accrued through 2050; effects after 2050 are excluded.

Recalculate deterministic point-value comparisons from saved annual point results where available. Do not use full-horizon point totals as the comparator.

Preserve existing analysis outputs and schema where practical: write clearly named conference-specific outputs in a separate `conference_2026_2050/` subdirectory under the selected VOI run. Rendering the Rmd must only read completed results; it must not source the runner or execute the analysis.

Do not silently drop failed draws. Report requested, successful, and excluded counts in technical notes. If existing failures or missing paired years invalidate inference, stop the affected analysis with an actionable error outside the deck. Do not present biased complete-case estimates as valid VOI.

## 4. Include BCR and connect it to the decision

Use discounted **incremental** monetized benefits relative to reference care divided by discounted incremental programme costs. Identify this definition explicitly.

Primary BCR outputs:

- Median and 95% PSA uncertainty interval.
- Ratio of expected benefits to expected incremental costs.
- Probability BCR exceeds 1, with optional thresholds 2, 3, and 5.
- A simple probability-versus-BCR-hurdle curve, if it advances the narrative.

Distinguish `E[B/C]` from `E[B]/E[C]`; do not choose the intervention by the mean draw-level ratio. For positive incremental costs, BCR > 1 and INB > 0 are equivalent. Explain that the existing monetary VOI at the reference valuation already addresses that decision.

Treat zero/negative incremental-cost draws separately. Do not reverse inequalities, remove such draws silently, or compute unstable ratios without disclosure. Use net benefit to handle these configurations.

Keep higher BCR hurdles mainly descriptive. If VOI is presented for a hurdle r > 1, explicitly define the objective `G_r = B - r*C`. Explain in backup that this reweights costs and is a different decision objective; do not label its value as ordinary welfare gains. The connection to the existing VSL-multiplier analysis is `G_r = r * [(1/r)*B - C]`.

## 5. Standardize VOI relative to programme costs, GDP, and government health expenditure

For each country, express the joint SAP EVPPI and individual EVPPIs in:

1. Constant-price discounted US dollars.
2. Multiples of programme costs.
3. Multiples of GDP.
4. Multiples of government health expenditure.

Use the same fixed country-level denominator for all parameter EVPPIs at a given valuation scenario. Divide bootstrap estimates by that same denominator so intervals remain consistent.

**Programme-cost denominator:** use expected discounted total scale-up programme costs over 2026–2050, derived from annual `sap` arm costs across the same valid draws. Label this as gross programme costs. Report expected discounted incremental costs separately; if showing VOI/incremental-cost ratios in backup, label them distinctly and require a positive denominator. Do not conflate the two.

**GDP and health-expenditure denominators:** prefer a documented 2026 annual GDP and annual **domestic general government health expenditure (GGHE-D)**, converted to the same constant-price US-dollar basis. These are contextual scale comparisons: explicitly label them as “2026–2050 VOI / one year of GDP” and “2026–2050 VOI / one year of government health spending.” Do not present these as shares of cumulative 25-year resources.

If reliable 2026 values are unavailable, use the latest verifiable observed reference year, preferably a common year, and print that actual year in the label. Do not invent forecasts or assume GDP equals GNI. Do not use PPP international dollars against market-exchange-rate US-dollar VOI. Distinguish GGHE-D from total health expenditure and from a ministry-only budget.

Consult existing project inputs first; retrieve missing denominator data from authoritative sources such as national accounts/World Bank and WHO GHED. Save source, indicator code, country, reference year, original units, price conversion, and denominator definition in conference-output metadata and Rmd comments. Source these comparisons in the presentation.

If a denominator cannot be verified, omit only that unsupported comparison and explain the omission in the completion report; do not substitute a fabricated number or display a warning/placeholder on a slide.

For readability, multiples below one may also be displayed as percentages, with the underlying multiple retained in the saved table. These ratios communicate scale, **not fiscal savings, available research budgets, or research benefit-cost ratios**.

## 6. Be precise about the information being valued

Establish from current code whether these remain the four uncertain inputs:

- Full-adherence A→B progression-hazard multiplier.
- Full-adherence B→C progression-hazard multiplier.
- Attained-age protection cap.
- Relative protection under partial adherence.

State clearly, in plain language, that the analysis estimates the value of perfectly resolving these selected uncertainties **for choosing reference care versus the current screening + SAP scale-up package**, conditional on the fixed model, calibration, initial state, implementation cascade, costs, and valuation.

Joint SAP EVPPI equals EVPI only within this restricted uncertainty model. Individual EVPPIs integrate over the conditional joint distribution, are overlapping, and must not be summed or stacked.

If still true in current code, disclose that SAP effects on C→D and direct D-mortality are fixed at no effect. Do not imply that advanced-disease prophylaxis or every SAP uncertainty has been assessed.

Distinguish:

- PSA outcome intervals from Monte Carlo/bootstrap estimator intervals.
- Biological efficacy from achieved adherence and coverage.
- Perfect information (EVPPI) from a feasible study (EVSI).
- Value of information from the economic benefits of implementing SAP.
- A mortality-based DALY proxy from modelled DALYs, if the proxy is used.

Interpret zero VOI cautiously: it can mean the same option is preferred across sampled states at that valuation and horizon. It does not establish absence of disease-model uncertainty or absence of scientific value.

## 7. Conference narrative and country-specific slides

Plan approximately **12–14 main slides including title and closing**, paced to 15 minutes, followed by technical backup and references. Use assertion titles supported by the computed results, one message per slide, minimal prose, and no literature-review or pipeline walkthrough.

A workable arc is: why uncertainty matters; the investment decision and BCR; the selected uncertainty and what learning could change; country-specific information value; research priorities; the strongest limitation; and a defensible closing message.

Use the requested recurring label:

**Pillar 3 · What would better evidence be worth?**

Treat it as a small kicker/subtitle, not a duplicate title or an extra divider on every slide. Combining countries is fine for framing and synthesis, but include **dedicated country-specific slides for both Indonesia and Uganda** addressing:

- **“The value of knowing how well secondary prophylaxis works”**: show joint information value, the country-specific decision/BCR context, and its scale relative to programme costs, GDP, and government health expenditure. Use a compact graphic or a headline number with a concise scale comparison.
- **Which uncertainty is worth resolving**: show individual parameter EVPPIs with uncertainty and joint value separately. If necessary, place the detailed ranking in country-specific backup to preserve the 15-minute pace.

**Do not precommit to “Research on effectiveness is worth more than it costs.”** That claim requires a specified study, credible research costs, and an appropriate expected benefit from that study. EVPPI alone is an upper bound, and VOI/programme-cost ratios do not establish research cost-effectiveness. If study economics are absent, use an evidence-supported assertion such as “Perfect information about SAP effectiveness is worth up to [computed value] in [country]” or “Better effectiveness evidence could change the scale-up decision.” If no decision change is supported, say so in an appropriately restrained title.

End with a message justified by both countries' results, acknowledging differences rather than forcing a common positive conclusion.

## 8. Design, figures, and limitations

- Use the beautiful-deck conference standard and the existing template, without redesigning the theme.
- Prefer ggplot2 with the **RColorBrewer PuOr palette**: `scale_colour_brewer(palette = "PuOr")`, `scale_fill_brewer(palette = "PuOr")`, or a fixed named PuOr-derived vector. Keep parameter/country mappings consistent. Do not imply ordinal meaning through colour where none exists.
- Highlight the key result; use muted context and direct labels. No duplicate plot titles when the frame title already states the message.
- Avoid tiny text, crowded multipanel figures, dense parameter tables, and engine diagrams in the main talk.
- Use country-specific scales where required for legibility, clearly marked; do not let one country's monetary magnitude obscure the other.
- All analytical numbers, including title numbers, must come from saved results through inline R or explicit data bindings.
- Include one concise, substantive main-talk slide addressing the strongest limitation: selected provisional SAP uncertainties, fixed calibration/model assumptions, and an upper bound on research value. Move detailed priors, evidence mapping, C→D exclusions, diagnostics, equations, denominator definitions, and sensitivity analyses to backup.
- **Do not print runtime warnings, error traces, “pending,” TODOs, or verification placeholders in the rendered main text.** Keep technical issues in comments/logs/completion reporting. This does not permit hiding material scientific limitations or unresolved invalid results.
- Capture warnings in external render diagnostics and fix their causes; do not suppress genuine problems merely to achieve a clean-looking PDF. Resolve compatibility with any existing gate within the allowed files.

## 9. Validation and deliverables

Run script 04's conference analysis path only, then render the Rmd to Beamer PDF with the project's supported toolchain. Keep all outputs traceable to the selected run and 2026–2050 horizon.

Check:

- Paired draw completeness and annual horizon.
- Recomputed benefits, costs, BCR decisions, and positive/negative cost handling.
- Joint and individual VOI bounds, estimator diagnostics, and uncertainty labels.
- Correct decision direction when VOI is zero; never hardcode “all draws favour scale-up.”
- No undefined-ratio summary failures or misleading rank-first probabilities when estimates tie.
- Denominator definitions, source years, conversions, and matching normalization of bootstrap intervals.
- Country-specific claims, citations, 15-minute pacing, and separation of backup numbering.
- Clean compilation and visual inspection of every rendered page for overflow, clipping, missing glyphs, and unreadable text.
- A final diff confirms that only the two authorized implementation files changed, plus any verified new reference entries appended to `library/references.bib`.

If required results, template, or rendering dependencies are unavailable, deliver all feasible work and report the exact blocker outside the slides. Do not claim a render or validation succeeded unless it was run.

Return the new Rmd and rendered PDF paths, the script-04 command used, a brief account of the principal findings and what is being estimated, the denominator choices and sources, the citation keys of any references added to `library/references.bib`, and the actual validation status.

