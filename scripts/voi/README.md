# SAP value-of-information analysis (separate from production)

The smallest SAP-focused value-of-information (VOI) analysis for Indonesia and Uganda. It compares
the two existing options, **reference care (`ref`)** and the **current screening + SAP scale-up
(`sap`)**, while four SAP inputs are uncertain. It is a **separate analysis**: production scripts
`00`–`08`, `R/abcd_engine.R`, calibrated bundles, baselines, outputs, reports and slides are not
modified, sourced for execution or re-run.

> **Status: exploratory and preliminary.** The priors are **provisional and not expert-agreed**.
> The production baseline it conditions on is labelled `INADEQUATE` (04b fit outside its GBD
> tolerances in both countries). EVPPI here is an **upper bound** on the value of further research,
> conditional on the fixed model. It is not EVSI, and it does not show that any particular trial
> is worthwhile.

## Run

From the repository root, with R 4.5.1:

```bash
"C:/Program Files/R/R-4.5.1/bin/Rscript.exe" scripts/voi/00_run_voi.R
```

or `source("scripts/voi/00_run_voi.R")` in R. Stages run in order:
`sample -> validate -> pilot -> psa -> analyse`. Set `VOI_STAGES` before sourcing to run a subset.
Each stage script can also be sourced on its own; it then loads the configuration only
(`VOI_CONFIG_ONLY`). Runtime is about 2–3 min for both countries on 12 workers: PSA ~27 s and
bootstrap ~50 s per country, ~0.27 s per paired draw serially.

**Extending the draws later.** Raise `N_DRAWS` in `00_run_voi.R` and keep `VOI_RUN_ID`, `VOI_SEED`
and the prior unchanged. Draw *i* always comes from the *i*-th L'Ecuyer-CMRG stream of the seed, so
draws 1..500 are reproduced exactly; 01 verifies this and stops on any mismatch. Chunk checkpoints
carry a signature (draw ids, parameter values, input/engine/helper md5), so only the new draws are
simulated. A changed prior or seed is refused for an existing run id: use a new `VOI_RUN_ID`. The
draw count is **not** increased automatically.

## Files

| File | Responsibility |
|---|---|
| `00_run_voi.R` | Independent configuration and entry point: paths, run id, `N_DRAWS <- 500L`, seed, tolerances, **final prior specification**, valuation thresholds, stage runner, production-integrity check before/after |
| `R/abcd_engine_voi.R` | Isolated copy of the production engine, minimally extended (search `VOI EXTENSION`): `abcd_sap_hr()` accepts the production shared multiplier **or** separate `sap_hr_A_to_B` / `sap_hr_B_to_C`, each in (0, 1], with optional per-transition RRR consistency checks; supplying both forms is refused |
| `R/voi_functions.R` | Shared helpers: write guard (`output/voi/<run_id>/` only), logging, manifest, isolated engine loading, saved-input loading, the paired per-draw runner, economics |
| `01_sample_sap_parameters.R` | Joint SAP sampling (per-draw RNG streams), prior diagnostics and prior figure |
| `02_validate_voi.R` | Production-equivalence and analysis checks (below) |
| `03_run_psa.R` | Pilot, then paired ref / scale-up projections and economics (parallel, resumable checkpoints) |
| `04_analyse_and_plot.R` | EVPI, single-parameter EVPPI, bootstrap, thresholds, summaries and figures |

Outputs go only to `output/voi/<run_id>/<country>/`:

- `draws/`: sampled inputs, prior spec
- `validation/`: checks, point run, pilot reproducibility, production integrity
- `pilot/`, `psa/`: compact annual summaries per draw × arm × year, per-draw totals, flags,
  checkpoints
- `tables/`, `figures/`, `logs/`
- `manifest.{rds,csv}`: input fingerprints, engine version, seed, prior and fixed assumptions

The current run is `sap4_n500_2026-11-10` (~51 MB, mostly checkpoints and the bootstrap array).

**Isolation.**
- **Engine loading:** the VOI and production engines are each loaded with `sys.source()` into
  their own environment, and PSA workers are separate processes.
- **Global environment:** validation confirms no `abcd_*` object is created there.
- **Write guard:** every write passes a guard that refuses paths outside `output/voi/<run_id>/`.
- **Production fingerprints:** md5 of 146 production files (R/, scripts/, data/<country>/,
  output/<country>/out_model and tables, the country comparison) is taken before and after each
  run and compared (`validation/production_integrity.csv`). It was unchanged after every run.

## Decision problem and what is fixed

- **Options:** `ref` and `sap` exactly as in production 06.
- **Outcome:** incremental net benefit (INB) of `sap` vs `ref`, in constant 2025 US$, discounted at
  3% to 2026, national population, 2026–2100.
- **Recomputed for every draw:** service volumes (screening encounters, on-BPG person-years,
  surgeries), costed with the **saved applied 2025-US$ unit costs** from `rhd_economic_results.rds`.
  No price conversion is repeated.
- **Valuation:** RHD deaths averted, valued with the **saved annual VSL schedule**.
- **Uncosted:** confirmatory echoes, as in production.
- **Fixed** (read from the saved `baseline_state.rds`):
  - calibrated natural history and incidence (Indonesia adult multiplier 0.1, D death 0.17 / 0.05;
    Uganda 0.5, 0.06 / 0.03);
  - the bridged 2026 engine state;
  - care-cascade inputs and ramps;
  - surgery service;
  - C→D and direct D-death SAP multipliers (= 1).
- **Same draw in both arms:** each draw's SAP inputs apply to both arms, including reference care's
  small on-BPG share. The care schedule is rebuilt for each draw's partial-adherence value, and
  `abcd_project()` re-prepares the effect matrices for every draw (checked, V10).

Because the four SAP inputs are the only uncertain quantities and the model is deterministic
given them, **joint SAP EVPPI = EVPI conditional on the fixed model assumptions**.

## Parameters and priors (provisional, not agreed)

| Input | Production point | Prior used here | Source / status |
|---|---|---|---|
| Regime indicator (shared by M_AB, M_BC) | — | Bernoulli, w_trial = **0.5 fixed** | workbook S3; 0.25 / 0.75 explored by re-weighting only |
| M_AB (full-adherence multiplier, treated A→B) | 0.45 (shared M) | Trial regime: LN, 95% 0.03–0.29 (median 0.093). Recurrence regime: 1 − s(1 − RR_rec) | see mapping below |
| M_BC (full-adherence multiplier, treated B→C) | 0.45 (shared M) | Trial regime: LN, 95% 0.068–0.59 (median 0.20), log-correlation 0.5 with M_AB. Recurrence regime: = M_AB | see mapping below |
| RR_rec (RF-recurrence RR) | — | LN, 95% 0.22–0.69 | Bray et al. 2024 (6 RCTs) |
| s (share of early-disease progression hazard driven by recurrent ARF) | — | Beta(8, 2): mean 0.80, 95% 0.52–0.97 | **new provisional judgement**, needs elicitation |
| `sap_protection_age_max` | 39 | 20 / 39 / 100 with P = 0.25 / 0.50 / 0.25 (100 = top model age = no cap) | workbook S8 |
| `partial_adherence_relative_protection` | 0.5 | Beta(4.537, 4.537): mean 0.5, 95% 0.20–0.80 | workbook S9; de Dassel et al. 2018 |

Truncation: M ∈ (0, 1]. A draw with M > 1 is redrawn, which happened 0 times in the 500 PSA draws
and 254 times in 200,000 prior Monte Carlo draws.

Realised prior:

| Input | Median | 95% interval | P(≥ 0.45) |
|---|---|---|---|
| M_AB | 0.30 | 0.036–0.73 | 0.36 |
| M_BC | 0.38 | 0.08–0.74 | 0.39 |

- **Dependence:** Spearman(M_AB, M_BC) = 0.89; M_AB = M_BC in the ~50% recurrence-regime draws.
- **Independence:** protection age and partial protection are independent of M.
- **Countries:** the same draws are used in both countries; the efficacy prior is not
  country-specific.

### Evidence-to-parameter mapping (resolved explicitly)

- **M is full-adherence efficacy.** The engine applies `h × [1 − a (1 − M)]`. Programme adherence
  `a = good + rel × partial` is modelled **separately** (reference 0.6 + 0.2 rel; scale-up target
  0.85 + 0.1 rel), and `rel` is one of the four uncertain inputs. The efficacy extrapolation (M) and
  the adherence response (a, rel) are therefore distinct, and trial evidence is mapped to M only.
  GOAL's ITT effect under high-adherence trial delivery is taken as approximately full-adherence
  efficacy. This is slightly conservative and not adjusted.
- **Pooled GOAL progression is not an A→B estimate.**
  - *What the endpoint is:* the GOAL / Bray 2024 early-RHD estimate (3/399 vs 33/400; RR 0.09,
    95% 0.03–0.29; constant-hazard HR 0.088) pools borderline (~A) and mild-definite (~B)
    children. Its endpoint covers both A→B-like and B→C-like progression.
  - *The B→C subgroup:* 00 block H records the mild-definite subgroup risks 0.119 vs 0.025
    (HR 0.20), which inform M_BC.
  - *The implied A→B subgroup:* 00 also records 0.075 vs 0.003 (HR 0.039, unverified) for the
    borderline subgroup.
  - *Why the pooled interval is used for A→B:* the A-subgroup effect looks stronger than pooled,
    and the subgroup evidence rests on very few events. The pooled interval contains the
    borderline-subgroup value near its lower end, so it is used for M_AB in the trial regime as a
    deliberately wide, verified-count-based **approximation**, not as an A→B estimate.
  - *Still to verify:* both subgroup figures against the NEJM supplement.
- **An RF-recurrence RR is not a progression multiplier.**
  - *Bridge:* if a share *s* of the untreated progression hazard of early disease is driven by
    recurrent ARF, and BPG multiplies the recurrence hazard by RR_rec, then
    `h_BPG = h [(1 − s) + s RR_rec]`, i.e. **M = 1 − s (1 − RR_rec)**.
  - *Both transitions:* the same for A→B and B→C (one mechanism), so M_AB = M_BC in this regime.
  - *Relation to the workbook:* s = 1 reproduces the workbook's recurrence component (M = RR_rec)
    and the logic behind production's 0.45 (Manyemba RR). s < 1 is the recurrence-mediated reading.
  - *Approximation:* RR ≈ HR at low recurrence risk.
- **Mixture weight** w_trial is fixed within the base prior scenario. Changing it is a different
  prior scenario (re-weighting table only).

### Workbook (`dev/RHD_EVPPI_parameter_register.xlsx`) — how it was used

Used as a reference only and **never read at run time**. Names, units, age boundaries, hazard
interpretation and adherence mapping were checked against the current engine and the saved
baseline. Resolutions:

1. **S1 (M_AB) GOAL component:** pooled RR used, but re-labelled as an approximation for A→B (see
   above).
2. **S1/S2 recurrence component:** the workbook used M = RR directly. Here the share *s* bridge is
   added (new provisional Beta(8, 2)), and M_AB = M_BC within that regime (the workbook drew them
   with log-correlation 0.5).
3. **Natural-history values quoted in the workbook / CLAUDE.md** differ from the saved Indonesia
   baseline. The saved baseline has D death 0.17 / 0.05, adult multiplier 0.1, entry mix
   0.75 / 0.25. The VOI analysis uses the saved baseline; none of these are sampled.
4. **Protection age "no cap" = 100**, the top model age (accepted by the engine's integer
   validation).
5. **Thresholds:** the workbook suggests 0.5× / 1× GDP per capita per DALY. The saved inputs
   carry only the 08 income path (WDI PPP GNI per capita grown at `gdp_growth`), so the
   DALY thresholds are multiples of **GNI per capita in 2026, constant 2025 US$, derived from that
   path**: Indonesia US$4,830, Uganda US$1,322. They are illustrative, not an endorsed
   willingness-to-pay.
6. **C→D effect (S4–S6), direct mortality (S7), waning (S10):** not implemented (pending, below).
7. **Draws:** the workbook proposes 3,000–5,000; this initial run uses 500 as instructed.

## Validation (stage 02; all hard checks pass in both countries)

| Check | Result (both countries) |
|---|---|
| V1 Production engine re-run on the saved baseline reproduces the saved 06 outputs (annual deaths, deaths 0–19, prevalence, screening, on-BPG PY, surgeries, confirmations, diagnoses, protected PY) | max rel. diff ≤ 1.4e-15 (tol 1e-9) |
| V2 VOI engine, shared multiplier = production engine (annual sums, every cycle matrix, final state) | **exactly 0** |
| V3 VOI engine with separate M_AB = M_BC = 0.45 = shared form; effect matrices `identical()` | **exactly 0** |
| V4 Economics from saved unit costs / VSL / discounting reproduce saved 08 totals, BCR, cost/death, cost/DALY and every year's arm costs and deaths | ≤ 3.6e-15 (tol 1e-9) |
| V5 Each multiplier scales only its own transition, at ages ≤ protection age | exact |
| V6 M_AB = M_BC = 1 ⇒ deaths averted ≈ 0 (≈1e-11), INB = −incremental cost | pass |
| V7 Protection-age mask 20 / 100 on attained age | exact |
| V8 Partial protection rebuilds the care schedule (ref 0.64; scale-up 0.64 → 0.87 at rel 0.2) | exact |
| V9 Engine refuses 11 invalid/ambiguous SAP specs (M ≤ 0, > 1, NA, one split field only, shared + split, RRR mismatch, non-integer / out-of-range age, C→D or D-death ≠ 1) | pass |
| V10 Stale cached effect matrices are overwritten for every draw | exact |
| V11 Engines isolated; global environment free of engine objects | pass |
| PSA | 500/500 draws per country, 0 numerical failures (deficits, care violations, onset cap, non-finite, negatives); pilot draws 1–20 reproduced exactly in the main run |

Informational result at the point values: the B→C channel alone gives 6,304 (IDN) / 6,858 (UGA)
deaths averted, and the A→B channel alone 835 / 933, against 6,784 / 7,386 with both.

## Methods

- **INB per draw:**
  - VSL family: `k × B − C`, where B = Σ df·VSL_t·deaths averted_t. **k = 1 is production's
    reference-case VSL.**
  - DALY family: `λ × D − C`, where D = Σ df·deaths averted_t × 30. This is the existing
    **deaths-to-DALYs proxy** from 08, not a DALY model.
  - C = discounted incremental cost in both families.
- **Reference thresholds (ranked comparison):** VSL × 0.05, × 0.10, × 1; and 0.5× and 1× GNI per
  capita per DALY. Curves cover k = 0.005–2 and λ = 0–3× GNI per capita.
- **EVPI = joint SAP EVPPI:** computed directly, `mean(max(0, INB)) − max(0, mean(INB))`.
- **Single-parameter EVPPI:** learn one input, with the other three **integrated over their
  conditional joint prior**. The M_AB/M_BC dependence is preserved; this is never a one-way
  analysis.
  - *Continuous inputs:* a single-input GAM (`mgcv`, thin-plate, REML) of each net-benefit
    component (B, D, C) on log M_AB, log M_BC or rel. Then `g = a·E[B or D | θ] − E[C | θ]`,
    which is exact by linearity, so one set of fits serves every threshold.
  - *Protection age:* a stratified estimate from category means and sample frequencies (counts
    137 / 222 / 141). This is a Monte Carlo estimate, not an exact calculation.
- **Precision:** 1,000 bootstrap resamples of **complete paired draw rows**, with every estimator
  refitted (fixed basis size).
- **Checks:**
  - GAM basis adequacy: k-index ≈ 1, all p > 0.24, no refits needed.
  - Equal-count 10-bin cross-check agrees within ~10%, except Indonesia at VSL × 0.1: binned 0
    vs GAM ~1 m, flagged exploratory.
  - Split-half agreement.
  - Running estimates at n = 100…500.
  - 0 ≤ EVPPI_i ≤ EVPI: no violation; no individual estimate exceeds the joint; nothing clipped.
- **Status labels:**
  - *stable*: bootstrap relative SE ≤ 0.20.
  - *negligible*: 97.5% bootstrap bound < 5% of the joint.
  - *EXPLORATORY*: otherwise, meaning more draws are needed.
  - *EVPI = 0*: every draw favours one option.

## Results (N = 500 paired draws per country; preliminary)

All amounts are US$ million, constant 2025, discounted to 2026, national population, 2026–2100.
Bootstrap 95% intervals are in brackets. `*` marks an exploratory estimate.

**Indonesia** — PSA means vs production point values:

| Outcome | PSA mean (95% interval) | Production point |
|---|---|---|
| Deaths averted | 8,314 (2,337–18,660) | 6,784 |
| BCR | 24.2 (6.9–51.8) | 20.1 |
| Cost per DALY (US$) | 4,499 (1,750–11,280) | 4,244 |
| Incremental cost (US$ m) | 222–247 | — |

| Threshold | E[INB] (MC SE) | P(INB>0) | Joint SAP = EVPI | M_AB | M_BC | Age cap | Partial adh. |
|---|---|---|---|---|---|---|---|
| VSL × 1 (reference) | 5,399 (126) | 1.000 (500/500) | 0 | 0 | 0 | 0 | 0 |
| VSL × 0.10 | 333 (12) | 0.912 | 4.5 [3.0–6.0] | 1.1 [0.3–2.5]* | 1.3 [0.4–2.7]* | 0 | 0 |
| VSL × 0.05 | 52 (6) | 0.588 | 28.4 [24.5–32.3] | 15.5 [11.8–19.7] | 15.8 [12.1–19.8] | **18.5** [15.0–22.1] | 0 |
| 1× GNI pc (US$4,830/DALY) | 85 (6) | 0.708 | 19.5 [16.2–22.7] | 10.4 [7.5–13.7] | **10.8** [7.9–14.0] | 8.9 [5.6–12.1] | 0 |
| 0.5× GNI pc (US$2,415/DALY) | −72 (3) | 0.128 | 6.6 [4.7–8.5] | 0.0 [0–0.6]* | 0.0 [0–1.6]* | 0 | 0 |

The joint value peaks at ~40 m around k ≈ 0.04, and at ~38 m around US$3,600/DALY.

**Uganda** — PSA means vs production point values:

| Outcome | PSA mean (95% interval) | Production point |
|---|---|---|
| Deaths averted | 8,682 (1,791–21,860) | 7,386 |
| BCR | 16.7 (3.5–40.9) | 14.2 |
| Cost per DALY (US$) | 1,633 (457–5,045) | 1,240 |

| Threshold | E[INB] (MC SE) | P(INB>0) | Joint SAP = EVPI | M_AB | M_BC | Age cap | Partial adh. |
|---|---|---|---|---|---|---|---|
| VSL × 1 (reference) | 1,006 (29) | 1.000 (500/500) | 0 | 0 | 0 | 0 | 0 |
| VSL × 0.10 | 43 (3) | 0.660 | 7.1 [6.0–8.2] | 1.5 [0.9–2.5]* | 1.7 [1.0–2.6]* | **6.0** [5.0–7.1] | 0 |
| VSL × 0.05 | −10 (1.4) | 0.366 | 8.5 [7.0–10.1] | 2.0 [0.5–4.2]* | 2.5 [1.0–4.7]* | **5.9** [4.4–7.6] | 0 |
| 1× GNI pc (US$1,322/DALY) | 15 (2) | 0.572 | 11.1 [9.8–12.4] | 4.5 [3.4–5.9] | 4.6 [3.5–6.0] | **8.8** [7.4–10.1] | 0 |
| 0.5× GNI pc (US$661/DALY) | −24 (1) | 0.120 | 2.1 [1.5–2.7] | 0.0 [0–0.1]* | 0.0 [0–0.3]* | 0.0 [0–0.01] | 0 |

The joint value peaks at ~15 m around k ≈ 0.06, and at ~14 m around US$1,060/DALY.

Full tables: `tables/inb_summary.csv`, `evppi_comparison.csv` (tidy, ranked, with bootstrap rank
probabilities), `evppi_curves.csv`, `decision_curves.csv`. Figures: `figures/evppi_comparison.png`
(threshold curves + ranked comparison), `decision_uncertainty.png`, `gam_fits.png`,
`convergence.png`, `prior_draws.png`.

### Interpretation

- **At the production valuation (VSL × 1) there is no decision uncertainty.** INB > 0 in every
  draw in both countries, so EVPI and every EVPPI are 0. Under these priors, resolving SAP efficacy
  would not change the ref-vs-scale-up decision.
- **Information value appears only at much lower valuations.** These are roughly k ≈ 0.02–0.2 of
  the VSL, or per-DALY thresholds near 0.5–1.5× GNI per capita, where the decision becomes
  uncertain. **The ranking depends on the threshold**:
  - *Uganda:* protection age ranks first at every reference threshold with positive value
    (P(rank 1) = 1.00).
  - *Indonesia:* protection age ranks first at VSL × 0.05, M_BC (just ahead of M_AB) at
    1× GNI per capita and VSL × 0.10.
- **Why the protection-age cap matters:** it changes both arms' costs and benefits. A later cap
  keeps treated B in B for longer, which means more BPG person-years (costed at every age) and more
  deaths averted. Its effect is close to a step change across three categories.
- **M_AB and M_BC have almost equal EVPPI** even though the B→C channel produces ~93% of the deaths
  averted. This is because they are strongly dependent under the prior: learning M_AB largely
  reveals M_BC. Their EVPPIs overlap; they are not separate contributions.
- **Partial-adherence protection carries negligible value** everywhere (bootstrap upper bound
  ≤ 0.06 m). The partial-adherence share is small (0.2 → 0.1), so `rel` moves the effective
  multiplier little.
- **Not additive.** Individual EVPPIs **must not be summed or stacked**. Correlation creates
  overlap, and interactions make the reverse also possible. At 0.5× GNI per capita, for example,
  the joint value is 6.6 m (Indonesia) and 2.1 m (Uganda), while every single-input value is ~0:
  only *combinations* (strong efficacy *and* a late age cap) reverse the decision. Read each value
  as the value of learning that input alone.
- **Prior-scenario sensitivity** (re-weighting the regime weight; joint quantities only; effective
  sample size ~400): at VSL × 0.05 in Indonesia, EVPI is 37.9 / 28.4 / 18.6 m for
  w_trial = 0.25 / 0.5 / 0.75. The values depend materially on the unagreed extrapolation judgement.
- **Upper bounds.** These EVPPIs bound what research on these inputs could be worth for this
  decision, under this model. **No trial is shown to be worthwhile:** that needs EVSI with a
  specific design, sample size and cost.

### Monte Carlo limitations at N = 500

- **Exploratory estimates:** 4 of 20 (Indonesia) and 6 of 20 (Uganda) reference-threshold
  individual estimates are flagged EXPLORATORY. These are mainly M_AB / M_BC at VSL × 0.05–0.10 in
  Uganda, and near-zero values at 0.5× GNI per capita.
- **Stable estimates:** the joint EVPI and the larger individual estimates are stable (relative
  SE ≤ 0.2, split halves agree).
- **E[INB]:** changes by < 5% between n = 400 and 500.
- **Running estimates:** at Indonesia VSL × 0.05 they drift between n = 100 and 300
  (age cap 24.6 → 18.3 m) and then settle.
- **Protection age:** each category has ≥ 137 draws.
- **More draws:** raising `N_DRAWS` (e.g. to 2,000–5,000, as the workbook proposes) would tighten
  the exploratory estimates. This was not done automatically.

## Limitations

- **Calibration held fixed.** The analysis conditions on the saved 04b calibration and the bridged
  2026 initial state. These were produced with the production point SAP values applied to reference
  care's on-BPG share (~4% Indonesia, < 1% Uganda); re-sampling SAP inputs without recalibrating is
  an approximation.
- **Preliminary baseline.** The baseline fit is `INADEQUATE` (Indonesia prevalence ~0.77×; Uganda
  RHD deaths 3.2–4.1× GBD).
- **Fixed inputs.** Natural history, care cascade, unit costs, VSL and surgery are fixed. EVPI is
  conditional on them; it is not the total decision uncertainty.
- **DALY proxy.** DALYs use the 30-per-death proxy (no age-specific YLL, no YLD).
- **Thresholds.** The GNI-per-capita thresholds are illustrative.
- **Evidence still to verify.** GOAL subgroup risks (00 block H) and the GOAL-Post durability data
  have not been checked against full text. *s* is a new, unelicited judgement.
- **Structural assumptions in the prior.** One shared regime for both transitions; protection age
  independent of efficacy.

## Pending — later extensions (not implemented)

- **Expanded SAP analysis (third option):** scale-up with **BPG continued after B→C**, plus
  uncertainty in the **C→D** SAP effect. This needs:
  - an on-BPG C stock (`T_C`) in the engine;
  - continuation / dropout rules for C;
  - the corresponding SAP service costs;
  - a C→D multiplier < 1 (relaxing the `sap_hr_C_to_D == 1` check).

  With only the current two options, EVPPI(θ_CD) = 0 by construction. **Not implemented now.**
- **Calibration sensitivity:** recalibrate (04b) at alternative SAP values, or sample calibrated
  parameter sets with their bridged states.
- **EVSI** for specific trial designs (e.g. extended GOAL-type follow-up, age-extension studies).
- **More draws** under the same run id; **prior scenarios** with fully re-run EVPPIs (beyond
  joint-quantity re-weighting); effect waning; direct-mortality scenario.
