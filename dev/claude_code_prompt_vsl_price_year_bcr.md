# Claude Code prompt: VSL correction, 2025 USD, BCR and per-capita outputs

Implement a narrowly scoped update to the Indonesia/Uganda RHD economic evaluation.

## Scope

- **Only two source files may change:** `scripts/00_run_all.R` (block J only; leave blocks A–I, K and L byte-identical) and `scripts/08_economic_evaluation.R`.
- **Do not create or modify any other repository file.** That rules out new scripts, helpers, data-raw CSVs, R/ functions, tests, docs, changelog files, `reports/slideck*.Rmd`, `analysis/report.Rmd`, 01–07, the engine, epidemiology, calibration, intervention effects, care coverage and the analysis horizon. Files that 08 writes when it runs (under `output/<COUNTRY>/` and `output/`) are run outputs and are allowed.
- **Macro observations go in block J.** Retrieve the official WDI/BLS values interactively during this session, from the WDI API or the WDI R package, outside the scripts. Declare them in block J as an explicit named list: one entry per series × country × year, each with value, unit, series code, source, WDI last-updated date and retrieval date. 08 reads only these globals and makes no network calls. If a value cannot be retrieved, leave an `NA` placeholder with a `# TODO: supply` comment. 08 must then stop with a message naming that observation. Never type a value from memory.
- **Use the existing model outputs as they are. This is an economics-only update, so do not re-run the model, and above all do not re-run the expensive calibration (04/04b).** Source only block J's parameters and `08_economic_evaluation.R`, reading the saved `output/<COUNTRY>/out_model/*.rds` for each country. Never source the full `00_run_all.R` driver, because it launches 01–08. If an existing output is missing or incompatible, stop and report it instead of regenerating it.
- Keep the existing service-bundle costing. Do not introduce microcosting, new services, new unit-cost estimates, additional GDP/productivity benefits, morbidity valuation or revised DALY methods. Confirmation remains uncosted and flagged.
- Do not implement or run sensitivity scenarios. Explicitly flag sensitivity as pending and expose parameters for later variation.
- Do not commit or push.

First inspect the actual current parameters and output metadata. Do not assume calibration is inadequate: preserve the status recorded in the existing outputs, and report it as unknown if unavailable. The workbook is a reference only, not a runtime input or user interface.

## 0. Downstream consumers: keep the reports working (read-only)

08's outputs and block J's globals will later feed `reports/slideck*.Rmd` and `analysis/report.Rmd`. Do not edit those files, but keep them working.

**Before changing anything,** read them (read-only) to find what they use from 08 and block J. Search for:
- file names: `rhd_economic_summary`, `rhd_budget_impact`, `rhd_economic_results`;
- the columns and RDS list elements they reference;
- block-J globals such as `vsl_mult`, `gni_pc_base`, `gdp_pc_base`, `cost_*`, `disc_rate`, `dalys_per_death`;
- any hard-coded wording such as "30 x GNI", "base-year US$" or "2019 US$".

Record this list.

**Compatibility rules:**
- **Additive changes only.** Keep every existing output file name, CSV column name and column order. Append new columns at the end of each table. Keep every existing RDS list element under its current name and add new elements alongside. Existing headline columns (`total_cost`, `total_benefit`, `net_benefit`, `bcr`, `cost_per_death`, `cost_per_daly`, …) stay as aliases of the discounted 2025 USD results.
- **Keep existing global names.** If a report reads a global you are replacing (for example `vsl_mult`, `gni_pc_base`, `gdp_pc_base`), keep that name defined in block J with a `# DEPRECATED: retained for report compatibility; not used by 08` comment. Do not delete it.
- **Labels for reports.** Add a `report_labels` element to `rhd_economic_results.rds` so the Rmd files can print units without hard-coding them. Include:
  - `currency` ("USD"), `price_year` (2025), `currency_label` ("constant 2025 US$");
  - `comparison_currency_label` ("constant 2025 international $");
  - `discounting_label`, e.g. "discounted to 2026 at 3% per year";
  - `vsl_method_label`, e.g. "Reference-case VSL: PPP GNI, elasticity 1.5, floor 20 × GNI per capita";
  - `analysis_period_label`;
  - `sensitivity_status` ("PENDING; not run");
  - `calibration_status`, as recorded in the model outputs.
- **Definitions.** Add a `definitions` element: a named list giving the formula and unit for every new column, so the reports can cite them.
- **Report interface map.** In your final response, list:
  - every report dependency you found, and whether it still resolves unchanged;
  - every place where a report's *text or meaning* will now be out of date (for example a sentence saying "VSL = 30 × GNI per capita", or figures labelled in 2019 US$);
  - the new outputs and labels each report could use.

  Do not fix the reports. That is a later task.
- If a report would break (an object it reads disappears or changes type), treat that as a defect in your change and fix it **in 08 or block J**.

## 1. Reporting currency and price year

Use `ECON_CURRENCY = "USD"` and `PRICE_YEAR = 2025L`, meaning **constant 2025 US dollars**, for every monetary result. Do not fall back to 2024 or silently substitute another reporting year. Keep `discount_base_year = ANALYSIS_YEAR_START`; a result may therefore be in constant 2025 USD discounted to 2026.

Preserve the current bundled price amounts and record each original price year separately. For the six existing prices, retain the model's **assumed 2019 USD** basis unless the current code already contains documented 2025 values. Do not inflate a price twice.

Treat the existing Indonesia and Uganda prices as **local healthcare service bundles reported in USD**, not U.S. resource costs. Use a country-specific local-CPI/exchange-rate adjustment as the declared approximation for these aggregate bundles. Keep them bundled: no ingredients, imported/local splits or microcosting are required.

For a bundle in country c, originally reported in USD for year y0:

```r
cost_LCU_y0 = cost_original_USD * XR_c_y0
cost_LCU_2025 = cost_LCU_y0 * CPI_c_2025 / CPI_c_y0
cost_USD_2025 = cost_LCU_2025 / XR_c_2025
# Combined: cost_original_USD * XR_c_y0 *
#           (CPI_c_2025 / CPI_c_y0) / XR_c_2025
```

All exchange rates are **LCU per USD**. Recover the original local amount using the conversion rate reported by the source study when available. If unavailable, use the official annual-average exchange rate for y0 as an explicit proxy and flag that assumption. Convert to reporting USD at the 2025 annual-average market exchange rate. Indonesia uses its own CPI and IDR/USD rates; Uganda uses its own CPI and UGX/USD rates.

For a cost originally in LCU, omit the first USD-to-LCU conversion. For an already-2025 USD price, retain the value without applying inflation again. If a bundle actually comes from a different source country, identify that provenance rather than assuming that its display currency establishes its country of origin.

Retrieve official country CPI observations (for example WDI `FP.CPI.TOTL`) and exchange rates (`PA.NUS.FCRF`) for y0 and 2025, from a compatible source/vintage. Record the CPI definition and index base. **US CPI is not the automatic inflator for a local service merely because it is quoted in USD.** US inflation or updated international prices may be appropriate for an established internationally priced imported commodity, but do not assume an entire service bundle is such a commodity or construct a component split in this task.

For every bundle, save source setting, original value/currency, assumed/verified original year, historical conversion rate and its source/proxy status, local CPI values and series, 2025 exchange rate, adjustment factor, adjusted value, reporting currency/year and sources. Missing macro observations must be identified explicitly; do not substitute U.S. inflation or a 2024 reporting year silently.

Set `COST_PRICE_ADJUSTMENT_METHOD = "local_cpi_xr"` explicitly. This method is an aggregate local-service price adjustment, not a new estimate of the bundle's resource requirements. Do not also apply US CPI to its adjusted 2025 USD result.

## 2. Correct VSL using the guideline's preferred default

Apply the global-health BCA reference-case preferred default: PPP-based GNI per capita, cross-country elasticity **1.5**, temporal elasticity **1.5** for this default implementation, and a floor of **20 times real GNI per capita**. Do not implement the other valuation scenarios yet.

Replace the fixed `30 * GNIpc_Atlas` rule. **Atlas GNI is not a PPP-to-market-USD currency conversion.** Use explicit PPP and market exchange-rate conversion factors.

Separate these concepts:
- the guideline's U.S. reference observation and price year;
- the target-country income anchor year;
- the 2025 reporting price year;
- the model's analysis and discount years.

Use the published U.S. reference VSL **9.4 million in 2015 USD** and reference income **57,900 in 2015 international dollars** consistently. The often-cited ratio of 160 is rounded; do not impose it as an additional exact equality. Do not silently apply 160 to current U.S. income and call that the same reference estimate.

For a concrete implementation, retain **2019** as the target-country income anchor year Y. Retrieve official 2019 observations for IDN and UGA:

| Series | Unit |
|---|---|
| `NY.GNP.PCAP.PP.CD`: GNI per capita, PPP, current international dollars | 2019 international dollars/person |
| `PA.NUS.PPP`: GDP PPP conversion factor | 2019 LCU/international dollar |
| `FP.CPI.TOTL`: country CPI, annual observations | Same index series in 2019 and 2025 |
| `PA.NUS.FCRF`: annual-average official exchange rate | 2019 and 2025 LCU/USD |

Verify that the PPP conversion factor is compatible with the PPP income series and vintage. Declare the observations in block J (see Scope), with units, year, source, retrieval date and available release metadata. No network calls in 08. Missing or duplicate required observations must produce a specific error; do not invent data or substitute Atlas GNI.

### Align the reference price basis

Express both U.S. reference quantities on the income anchor's price basis before comparing with 2019 target income. **US CPI is appropriate here because these are U.S. reference values**, unlike the local service bundles above. Use BLS all-items CPI-U, U.S. city average, annual average, not seasonally adjusted: 2015 = 237.017 and 2019 = 255.657.

```r
reference_price_factor = CPI_US[Y] / CPI_US[2015]
VSL_US_ref_Y = 9.4e6 * reference_price_factor
income_US_ref_Y = 57900 * reference_price_factor
income_target_Y = GNIpc_PPP_target_Y
raw_vsl_anchor_intY =
  VSL_US_ref_Y * (income_target_Y / income_US_ref_Y)^vsl_income_elasticity
```

Here the reference income is the published reference observation expressed in Y-price units, not a new U.S. income observation. Document all intermediate units.

### Real income growth and the floor

Retain the existing country-specific real per-capita growth assumptions, currently 3% for Indonesia and 3.2% for Uganda. Explicitly label GDP-per-capita growth as a **proxy for real GNI-per-capita growth**, not an official GNI forecast.

```r
growth_index_t = (1 + income_growth)^(year - Y)
income_target_t_intY = income_target_Y * growth_index_t
raw_vsl_t_intY =
  raw_vsl_anchor_intY * growth_index_t^vsl_time_elasticity
floor_vsl_t_intY = vsl_floor_mult * income_target_t_intY
vsl_t_intY = max(raw_vsl_t_intY, floor_vsl_t_intY)
```

Apply the floor **in every year**. Preserve the unfloored anchor; do not grow an already floored anchor at elasticity 1.5 and prevent the floor from releasing correctly. Do not impose a VSL/income upper bound of 160: that is not a universal guideline cap.

### Convert to constant 2025 market-rate USD

```r
vsl_t_LCU_Y = vsl_t_intY * PPP_LCU_per_intDollar_Y
vsl_t_LCU_2025 = vsl_t_LCU_Y * CPI_country[2025] / CPI_country[Y]
vsl_t_USD2025 = vsl_t_LCU_2025 / XR_LCU_per_USD_2025
```

This expresses benefits in the same 2025 local-price basis and at the same 2025 reporting exchange rate as costs. Keep the income-anchor PPP factor and 2025 reporting exchange rate fixed across the future valuation path. Real income growth is handled separately above. Do not introduce nominal future exchange-rate forecasts or additionally multiply the result by US CPI.

Then retain:
`benefit_t = (RHD_deaths_ref_t - RHD_deaths_sap_t) * vsl_t_USD2025`.
Keep `benefit_gdp = 0`.

Store the reference values, inflation adjustment, income ratio, unfloored and floored VSL, floor flags, PPP/XR factors, annual VSL, growth assumptions and sources. If VSL growth exceeds discounting, report that numerical fact without changing the parameters or claiming that late benefits necessarily dominate the totals.

## 3. Discounted and undiscounted results in 08

Let `C_t = cost_sap_t - cost_ref_t`, `B_t = benefit_t`, and
`df_t = (1 + disc_rate)^-(year - discount_base_year)`.

Add explicit summary columns:

- `total_cost_undisc = sum(C_t)`
- `total_benefit_undisc = sum(B_t)`
- `net_benefit_undisc = sum(B_t - C_t)`
- `bcr_undisc = sum(B_t) / sum(C_t)`
- `total_cost_disc = sum(df_t * C_t)`
- `total_benefit_disc = sum(df_t * B_t)`
- `net_benefit_disc = sum(df_t * (B_t - C_t))`
- `bcr_disc = sum(df_t * B_t) / sum(df_t * C_t)`.

BCR is undefined when its cost denominator is zero or negative: return NA and the corresponding cost-neutral/cost-saving status. Allow zero or negative benefit results; validations must test arithmetic and data validity rather than require favorable results.

Retain existing headline columns as aliases of the discounted results for backward compatibility. Preserve the formulas for `cost_per_death`, `cost_per_daly`, `dalys_per_death` and other out-of-scope outcomes, while consistently updating their monetary units.

Carry currency, price year, analysis start/end, discount rate and discount reference year in the CSVs and results metadata. BCR itself is dimensionless.

## 4. Cost and benefit per capita

Use WPP population from existing model outputs, all ages and both sexes. Verify one population observation per location/age/sex/year/arm and agreement between arms. Do not add the two arms' populations together.

For each location, let `N0` be population in the **first analysis year**. Report cumulative results **per resident at program start**:

```r
cost_pc_undisc = total_cost_undisc / N0
benefit_pc_undisc = total_benefit_undisc / N0
net_benefit_pc_undisc = net_benefit_undisc / N0

cost_pc_disc = total_cost_disc / N0
benefit_pc_disc = total_benefit_disc / N0
net_benefit_pc_disc = net_benefit_disc / N0
```

These cost measures are incremental costs, not gross scale-up-arm spending. Use the **same undiscounted N0 denominator** for discounted and undiscounted cumulative results. Label the denominator, its year and the unit explicitly. Do not call these annual per-capita spending or divide by discounted population.

Also export an annual table with `N_t`, `C_t`, `B_t`, `df_t`, `C_t/N_t`, `B_t/N_t`, `df_t*C_t/N_t` and `df_t*B_t/N_t`. Label the latter as present-value contributions per resident of year t.

For TOTAL, sum location monetary totals and starting populations, then calculate ratios. Do not average location BCRs or per-capita results, and do not double-count a pre-existing TOTAL row.

## 4b. Secondary reporting in constant 2025 international dollars (comparison add-on)

The **primary results stay in constant 2025 USD**: every calculation, the budget impact and all validations run in USD. Add parallel headline columns in **constant 2025 international dollars (INT$)** so Indonesia and Uganda can be compared in purchasing-power terms. This is a conversion of final results only. Do not run any calculation in INT$.

### Conversion

Convert each final 2025 USD amount back to 2025 local currency at the same exchange rate used for costs and VSL. Then convert to international dollars with the 2025 PPP conversion factor:

```r
# XR_c_2025  : PA.NUS.FCRF, LCU per USD, 2025 annual average (the same observation used in sections 1-2)
# PPP_c_2025 : PA.NUS.PPP,  LCU per international dollar, 2025 (GDP PPP conversion factor)
usd_to_int_2025_c = XR_c_2025 / PPP_c_2025            # INT$ per USD
value_INT2025     = value_USD2025 * usd_to_int_2025_c
```

Use one fixed factor per country for every year of the horizon. All amounts are already constant 2025 prices, so do not apply future PPP or exchange-rate paths.

### Data

- Add `PA.NUS.PPP` **for 2025**, for IDN and UGA, to the block-J macro observations. Use the same retrieval, metadata and no-network-in-08 rules as the other series.
- Record the ICP benchmark round behind the 2025 value. WDI extrapolates PPP factors from the most recent ICP benchmark (currently ICP 2021). Record whether the value is a benchmark or an extrapolation.
- Confirm it belongs to the same WDI vintage as the 2019 `PA.NUS.PPP` value used for the VSL.
- **If the 2025 PPP factor is missing**, still produce the primary USD results. Set every INT$ column to NA and print a message naming the missing observation. Never substitute another year or source silently.
- An alternative source, such as the IMF WEO implied PPP conversion rate, may be used only if declared explicitly as a block-J parameter with its source and vintage recorded.

### What to convert

Add `_int2025` columns alongside the unchanged USD columns for:
- discounted and undiscounted total cost, total benefit and net benefit;
- the six cumulative per-capita results from section 4;
- `cost_per_death` and `cost_per_daly`;
- VSL in the first and last analysis years (in the VSL derivation table);
- the adjusted 2025 unit-cost bundles (in the price-adjustment audit).

Carry `currency = "INT$"`, `price_year = 2025` and `usd_to_int_2025_c` (with its XR and PPP inputs and sources) in the CSVs and the results RDS.

Do not duplicate unit-free outputs: BCRs and deaths averted are the same in either unit.

### Ratios and pooling

- **Per country:** BCR in INT$ equals BCR in USD exactly, because both numerator and denominator scale by the same factor. Validate this to 1e-12 for each country and each discounting basis.
- **Pooled across countries:** this invariance does NOT hold, because each country has its own factor. A pooled multi-country TOTAL BCR or per-capita value differs between USD and INT$. If a cross-country TOTAL is reported, compute it separately in each unit from that unit's pooled totals, label it, and never mix units.

### Country comparison table

Each country's outputs are already written to its own `output/<COUNTRY>/` folder. 08 runs once per COUNTRY, so it can only build the comparison after both countries' summaries exist on disk.

When both country summaries exist, have 08 write `output/rhd_country_comparison_2025.csv` with one row per country. It should contain:
- the USD and INT$ totals, discounted and undiscounted;
- the per-capita results in both units;
- the discounted and undiscounted BCRs;
- `N0` and its year;
- the conversion factors and their sources;
- the calibration status recorded in each country's outputs.

Build the table by reading the saved summaries. Do not rerun anything. If one country's summary is missing, skip the table and say which one is missing.

Include a definitions note in the table metadata or the RDS. INT$ values describe purchasing-power-equivalent magnitudes for comparison. They are not budget figures. Budget impact remains in 2025 USD.

## 5. Explicit parameters; sensitivity pending

Declare the calculation inputs once in 00 block J and pass a named parameter object into the new/refactored calculation functions in 08. Keep changes proportionate; do not build a new evaluation framework.

Expose at least:
- reporting currency and price year;
- each original bundle value, original currency/year and source;
- source setting, country CPI series/year/value inputs, historical/2025 exchange rates, historical-rate proxy flag and price-adjustment method;
- US CPI inputs used only to align the U.S. VSL/reference-income price basis;
- VSL reference value, reference income and their price/reference years;
- income anchor year and the target income/PPP/XR input series;
- `vsl_income_elasticity = 1.5`, `vsl_time_elasticity = 1.5`, `vsl_floor_mult = 20`, and `vsl_floor_enabled = TRUE`;
- country-specific income-growth rate and proxy interpretation;
- discount rate and reference year;
- per-capita denominator convention and year;
- secondary reporting currency `ECON_COMPARISON_CURRENCY = "INT$"`, its price year (2025), the 2025 PPP conversion factor series (`PA.NUS.PPP`), and the source to use if it is missing (default: none, so INT$ columns are NA);
- unchanged economic assumptions still used elsewhere, including `dalys_per_death`.

Export a parameter register with value, unit, role, source, year and future-sensitivity flag. Save effective numerical inputs, not only filenames or NULL defaults. Do not invent sensitivity ranges or implement alternative-method switches for this task.

Print and save: **Sensitivity analysis: PENDING; not run.**
Flag future work on VSL/reference assumptions, elasticity, floor, income growth, discounting and bundle prices/price years. Changes to epidemiological effects would require rerunning the health model and are outside this economics-only update.

## 6. Validation and delivery

- Snapshot existing economic outputs before overwriting them. Copy them to `output/<COUNTRY>/tables/_pre_vsl_update/`; that is a run output, not a repository source file. Do not alter saved epidemiological outputs.
- **File-scope check:** before finishing, list every repository file whose content changed (for example with `git status` / `git diff --stat`, without committing). Only `scripts/00_run_all.R` and `scripts/08_economic_evaluation.R` may appear, apart from generated run outputs under `output/`. In `00_run_all.R`, confirm the diff touches only block J. Revert anything else.
- **Report compatibility check:** confirm that every object, column and global the reports read (from section 0) still exists, with the same name and type. Do not knit or edit the reports.
- Verify all six 2019-to-2025 bundle conversions independently using each country's CPI and exchange rates, and verify already-2025 USD values remain unchanged. Confirm that a USD display currency never selects the US inflator automatically.
- Check reference price-unit alignment, PPP-to-LCU-to-USD conversion, yearly floor behavior, and positive finite VSL.
- Independently reconcile discounted/undiscounted costs, benefits, net benefits and BCRs. At discount rate zero, the two sets must agree. This is a validation check, not a reported sensitivity scenario.
- Check every INT$ value equals its USD value × `XR_c_2025 / PPP_c_2025`. Check per-country BCRs are identical in USD and INT$ (1e-12). Check no USD column changed when the INT$ add-on was added. Check a missing 2025 PPP factor produces NA INT$ columns and a named message, not an error in the USD results.
- Check each cumulative per-capita result times N0 reproduces its matching total. Check population denominators are not duplicated by arm and TOTAL is correctly pooled.
- Confirm the underlying deaths, populations, service volumes, scenarios and horizon are unchanged.
- Run 08 for Indonesia and Uganda from existing outputs. Do not source the full 00 driver in a way that launches 01–07. If an input is unavailable, complete the code and identify exactly which observation blocks the run; never manufacture a result.
- Export the updated existing economic files, annual per-capita table, VSL derivation, price-adjustment audit and parameter register. When both countries have run, also export the INT$/USD country comparison table. Include definitions in the results RDS.
- Give a short changelog **in your final response** (no changelog file). Cover:
  - the files changed and the formulas;
  - source dates;
  - before/after results;
  - price-year assumptions;
  - the actual calibration caveat;
  - the report interface map from section 0;
  - **Sensitivity analysis: PENDING**.

  The same information goes in the results RDS (`definitions`, `report_labels`, parameter register).

References:
- Robinson et al. (2019), Reference Case Guidelines for Benefit-Cost Analysis in Global Health and Development, Chapters 3–4:
  https://www.pseau.org/outils/ouvrages/bmgf_reference_case_guidelines_for_benefit_cost_analysis_in_global_health_and_development_2019.pdf
- Turner et al. (2019), Adjusting for Inflation and Currency Changes Within Health Economic Studies, Value in Health 22(9):1026-1032. Local-currency inflation followed by base-year currency conversion; already-USD local costs can first be converted back to local currency:
  https://doi.org/10.1016/j.jval.2019.03.021
- BLS annual CPI-U, for U.S. reference-value price alignment only:
  https://www.bls.gov/regions/mid-atlantic/data/ConsumerPriceIndexAnnualandSemiAnnual_Table.htm
- World Bank WDI metadata and official series/API for the listed observations, including `PA.NUS.PPP` (PPP conversion factor, GDP, LCU per international $) and its ICP benchmark/extrapolation notes. Atlas conversion uses a smoothed conversion factor and must not replace the specified PPP/XR conversion.

