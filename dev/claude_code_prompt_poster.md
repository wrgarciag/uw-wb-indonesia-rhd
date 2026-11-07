Create **`reports/poster.Rmd`** to render a **self-contained HTML collection of poster content** that my coauthors can copy into `WCRHD2026_POSTER.pptx`. The HTML does **not** need to reproduce the poster layout.

**Scope**
- Only create or edit `reports/poster.Rmd`. Do not modify any existing scripts, reports, decks, inputs, parameters, caches, or model outputs.
- Keep all new R code, CSS, captions, and diagram definitions inside this Rmd.
- Read existing saved outputs only. Do not run the modelling or calibration pipeline.
- Test rendering in a temporary directory, keeping generated test files outside the repository.

**Read first**
- The current `reports/report.Rmd` (attached reference: `report(6).Rmd`).
- `WCRHD2026_POSTER.pptx`, for its requested content and visual style.
- `library/references.bib`, for the bibliography and existing citation keys. Cite relevant background claims, guidelines, data sources, intervention effects, and economic valuation methods using those keys (for example, `[@citation_key]`). Configure `poster.Rmd` to use this bibliography and render in-text citations and a concise reference list in the HTML. Verify that each cited entry supports the associated claim. Do not invent citation keys or references, and do not modify `library/references.bib`; record any missing references only in non-rendered source comments, and omit unsupported claims from the rendered text.
- Relevant saved outputs and their metadata for **Indonesia and Uganda**, plus upstream code as needed to understand their definitions.

Use the report’s loading logic and relevant consistency checks, adapted within the new Rmd. Treat the poster’s numerical placeholders, hypothetical country descriptions, and suggested coverage targets as draft text. Actual results and assumptions must come from the saved runs.

**Content**

Provide clearly labelled, independently reusable elements:

1. **Poster text:** title, authors and affiliations from the template, concise background, objective, methods, results, conclusions, and a short take-home message. Preserve verified author details and record unresolved information only in non-rendered source comments.
2. **Figure 1: model schematic.** Illustrate the current implemented model, including WHF stages, diagnosis/screening, prophylaxis, progression, and mortality. Verify the diagram against the current implementation rather than reproducing an older conceptual pathway. Keep it simple enough for a poster, with a concise caption.
3. **Figure 2: principal results figure.** Prioritize a clear comparison of Indonesia and Uganda showing cumulative RHD deaths averted, incremental costs, monetized benefits, net benefits, and BCRs. Use aligned panels with appropriate units rather than combining incompatible quantities on one axis. Include a complete caption and essential interpretation caveats.
4. **Supporting results:** a compact, copyable country comparison table and short numerical statements covering deaths averted through 2050, percentage reduction in RHD mortality in 2050, average annual incremental cost per capita, discounted net benefit, and BCR. Include the GDP-per-capita comparison only if supported by available data with a clear denominator.
5. **Essential limitations and references:** concise text suitable for the poster, including current calibration limitations and omitted costs where applicable.

**Analysis rules**
- Restrict all projection figures and results to **2026–2050 inclusive**.
- Recalculate economic totals from the saved annual economic series for this window using the saved discount factors. Do not reuse full-horizon summary totals when they extend beyond 2050.
- Derive BCR as discounted benefits divided by discounted incremental costs, and net benefit as their difference. Preserve signed results and handle undefined ratios explicitly.
- Distinguish cumulative deaths averted from the percentage reduction in annual mortality in 2050. Label counts and rates accurately.
- Read intervention coverage, adherence, treatment effects, currency, price year, and discounting assumptions from saved metadata.
- Generate numerical prose dynamically. Do not invent results, uncertainty intervals, or favorable conclusions.
- Determine calibration status separately for each country from saved diagnostics. Preserve material caveats without copying potentially stale hard-coded warnings.
- If required artifacts are missing or inconsistent, stop with an actionable error. Do not substitute another country, older runs, or fabricated values.

**Rendered content and internal notes**
- Do not display anything as pending in the rendered HTML. Do not include TODO/TBD labels, placeholders, missing-item notices, editorial queries, work-in-progress banners, or task-status flags in visible text, figures, captions, tables, or appendices.
- Record pending items, missing references, unresolved author details, and follow-up tasks only in HTML comments (`<!-- ... -->`) or R comment lines inside chunks with `include=FALSE` or `echo=FALSE`. Ensure these notes do not print into the HTML.
- If optional content cannot be supported, omit it from the rendered document and explain the omission in source comments. Missing or inconsistent required numerical artifacts must still stop rendering with an actionable error rather than produce fabricated results.
- Keep material scientific limitations as concise factual prose where needed for interpretation, without presenting them as pending tasks or editorial warnings.

**Design and usability**
- Use a simple scrolling HTML with selectable text, copyable tables, hidden code, and large static figures.
- Match the poster’s palette approximately: navy **`#0C2A55`**, dark blue **`#002664`**, red **`#D61C26`**, orange accent **`#F3631B`**, pale background **`#F4F7FA`**, and white.
- Keep country/scenario colors consistent and use labels or line types alongside color.
- Make figures suitable for transfer to PowerPoint, with readable labels and high-resolution PNGs embedded in the HTML. Provide embedded SVG download links where feasible, without additional permanent source files.
- Keep technical provenance and verification details in a short final appendix.

Render and visually inspect the HTML in a temporary directory. Check numerical agreement, the 2050 cutoff, captions, legibility, and clipping. Finish by reporting the file created, the render command, and any unresolved data limitations.

