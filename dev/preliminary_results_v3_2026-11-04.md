# Preliminary v3 results — Indonesia and Uganda (4 Nov 2026)

> **PRELIMINARY AND INADEQUATELY CALIBRATED — for discussion only.** Both calibrations are structurally
> valid and numerically converged but fall **outside the diagnostic fit tolerances in every scored year**
> (status `..._fit_outside_diagnostic_tolerance`; `baseline_state.rds` meta `calibration_adequacy` =
> "INADEQUATE"). Model RHD deaths are 3–4× GBD, so health-impact and especially economic results are
> **not credible estimates**. Economic results are **exploratory**.

**Run:** engine `abcd_pathway_YO_pooled_care_v3`, default settings in `scripts/00_run_all.R` (D-only, progression
SAP, `adult_progression_hr` = 1, six incidence bands, `STRUCT_BURNIN_RECONCILIATION = "per_capita"`, calibration
2000–2019 scored 2010–2019, bridge 2020–2025, analysis 2026–2100, scale-up ramp 2026–2050). Each country ran
independently, 03→08, using the existing 01/02 inputs. No assumptions or parameters were changed beyond the
existing calibration. Outputs: `data/<COUNTRY>/`, `output/<COUNTRY>/`. Ratios are model ÷ GBD (Rate × WPP
population), shown as 2019 [range over 2010–2019].

| Section | Metric | Indonesia | Uganda |
|---|---|---|---|
| **Status** | Pipeline 03→08 | Completed, all stages | Completed, all stages |
| | Calibration fit | **INADEQUATE** (outside tolerance) | **INADEQUATE** (outside tolerance) |
| | Optimiser: starts converged; objective L (prior-centre L) | 5/5; 1.504 (1.719) | 5/5; 4.717 (13.001) |
| | Error terms: prevalence / deaths 0–19 / deaths 20+ | 0.232 / 0.795 / 0.872 | 0.603 / 2.465 / 2.029 |
| | Parameters at bounds | childhood C→B (upper) | childhood C→B (upper); incidence multiplier 0–14 male (lower, 0.20) |
| | Integrity checks (structure, replays, bridge, seed, fixed inputs, 06–08 reconciliation) | 9/9 passed | 9/9 passed |
| **Fit vs GBD** | Prevalence, all ages | 0.73 [0.73–0.73] | 0.54 [0.50–0.54] |
| | Prevalence, ages 0–19 / 20+ | 0.86 / 0.72 | 0.54 / 0.53 |
| | RHD deaths, all ages | 3.05 [2.41–3.05] | 4.33 [3.51–4.33] |
| | RHD deaths, ages 0–19 | 1.28 [0.90–1.28] | 6.18 [5.16–6.18] |
| | RHD deaths, ages 20+ | 3.36 [2.74–3.36] | 4.06 [3.25–4.06] |
| | Scored years within tolerance: prevalence / deaths / 0–19 / 20+ | 0 / 0 / 9 / 0 of 10 | 0 / 0 / 0 / 0 of 10 |
| | Childhood share of RHD deaths (model vs GBD) | 6.2% vs 14.8% | 18.1% vs 12.7% |
| | Stage-D share vs GBD heart-failure proxy | 11.2% vs 3.2% | 5.1% vs 0.9% |
| | Incidence (validation only) | 0.80 | 0.24 |
| **2026–2050** | Living RHD entering 2026 (D share) | 985,872 (11.2%) | 465,134 (5.1%) |
| | Living RHD end-2050: reference / scale-up | 873,875 / 879,330 | 705,337 / 705,108 |
| | New RHD onsets (both arms) | 953,610 | 632,593 |
| | RHD deaths per year: 2026 / 2050 reference / 2050 scale-up | 17,714 / 16,172 / 15,290 | 1,915 / 2,886 / 2,814 |
| | Cumulative RHD deaths: reference / scale-up | 434,525 / 426,108 | 61,461 / 60,908 |
| | **RHD deaths averted** (% of reference); ages 0–19 | **8,417 (1.9%)**; 134 | **553 (0.9%)**; 16 |
| | Screening encounters: scale-up / reference | 58.6 M / 5.5 M | 18.3 M / 0.33 M |
| | On-BPG person-years: scale-up / reference | 7.37 M / 6.22 M | 4.38 M / 3.52 M |
| | 2050 on-BPG share / effective protection: scale-up vs reference | 40% / 36% vs 29% / 20% | 39% / 35% vs 27% / 19% |
| **Economics (EXPLORATORY)** — discounted 3% to 2026, US$ | 2026–2050 incremental cost / VSL benefit | 110.9 M / 1,248.5 M | 24.8 M / 17.4 M |
| | 2026–2050 BCR; cost per death averted; cost per DALY averted | 11.3; 13,179; 739 | 0.70; 44,807; 2,583 |
| | Full horizon 2026–2100 (08 output): BCR; deaths averted | 27.9; 53,453 | 3.68; 9,204 |
| | 2026–2050 scale-up cost shares: screening / BPG / surgery | 5% / 68% / 27% | 3% / 63% / 34% |

## Remaining limitations
- **Calibration inadequate in both countries.** Prevalence is too low (0.73× / 0.54×) and RHD deaths too high
  (3.1× / 4.3×). Deaths also drift upward relative to GBD over 2010–2019. The stage-D share is 3–6× the GBD
  heart-failure proxy.
- **Uganda.** Childhood deaths are 6.2× GBD even with childhood incidence multipliers at their 0.2 floor, and
  fitted incidence is 0.24× GBD. Uganda uses the Indonesia workbook's natural-history and care assumptions as
  extrapolations.
- **Missing adult exits.** GBD's own targets imply adult net outflows of about 2–3.6%/yr per prevalent case
  beyond RHD and background deaths (Indonesia cohort accounting). The model has only 2%/yr B normalisation, so
  it drains prevalence through D deaths. The mechanism behind GBD's outflow is unresolved: regression,
  deaths coded to other causes, or inconsistency within GBD.
- **Oldest ages.** Prevalence is overstated at 80+ (Indonesia 95+ ×6) because WPP populations shrink faster
  than GBD mortality implies during the observed years.
- **Provisional inputs.** The entry mixes are undocumented expert estimates and may describe stage at
  presentation rather than biological onset. Childhood regression rates, D mortality (column K) and the
  GOAL-derived hazard ratios extended to adults are all provisional.
- **Small modelled impact.** Scale-up raises effective protection from about 20% to 35% by 2050, averting
  1–2% of RHD deaths through 2050.
- **Economics are exploratory.** The VSL benefits scale with overstated death counts, DALYs use a flat
  30 per death, and the results are horizon-dependent (2050 vs 2100).
- **Not yet updated.** `tests/`, `analysis/report.RMD`, `reports/slideck-technical.Rmd` and `CLAUDE.md` still
  describe earlier output contracts.
- **Stale files.** `data/Uganda/calibration_structural_failed_candidate.*` are left over from an earlier v1
  attempt and are not current.

## Next steps
1. Decide, with clinical input, how to represent the adult net outflow implied by GBD (adult regression or
   normalisation versus non-RHD-coded excess mortality) before any further recalibration.
2. Confirm the basis of the entry mixes (onset versus presentation) and obtain Uganda-specific natural-history
   and care evidence.
3. Address the observed-year old-age population mismatch between WPP and GBD.
4. Once the structure is agreed, update the tests, report and deck to the v3 contracts.
