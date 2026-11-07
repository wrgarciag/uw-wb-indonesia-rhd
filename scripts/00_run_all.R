################################################################################
# INDONESIA RHD SECONDARY-PREVENTION MODEL — MASTER CONTROL PANEL
# scripts/00_run_all.R
# ------------------------------------------------------------------------------
# This is the SINGLE place to set every parameter needed to run scripts 01-08.
# Edit the blocks below, then source this file to run the whole pipeline:
#
#     source("scripts/00_run_all.R")
#
# Pipeline (each script reads the PERSISTED outputs of the previous ones):
#   01 prepare inputs  -> data-raw/temp_baseline_rates_gbd.rds
#   02 demography      -> data/pop_observed.rds, pop_projection.rds   (stable names)
#   03 disease inputs  -> data/disease_model_inputs.rds            (builder, not runner)
#   04 calibration     -> data/calibrated_rhd_parameters.rds       (calibrated IR/CF, 2000-2019)
#   05 initial state   -> data/baseline_state.rds                  (assembles 02+03+04)
#   06 run model       -> output/out_model/<location>.rds          (parallel by location)
#   07 make outputs    -> output/tables/rhd_model_long.csv/.rds    (+ stage + flow tables)
#   08 economics       -> output/tables/rhd_economic_summary.csv, rhd_budget_impact.csv
#
# DISEASE STRUCTURE (DEFAULT, engine "abcd_first_diagnosis_v1_adult_prog_mult"): World Heart
#   Federation RHD stages A/B/C/D in a national population model that separates
#   disease history from FIRST CONFIRMED DIAGNOSIS: undiagnosed and diagnosed (off /
#   on BPG) stocks by stage, adult prior-A lesion history (excluded from classifiable
#   prevalence), normalised history (no recurrence) and never-onset people, with
#   screening history. Rules depend on ATTAINED AGE only (<= 20 / > 20). New modelled
#   disease enters A (age <= 20) or B (age > 20); screening finds existing disease and
#   never creates it. Adult (> 20) B->C and C->D hazards = 0.48 x the childhood hazards
#   (block G, adult_prog_hazard_mult). RHD death occurs only from D (A/B/C exact zeros);
#   competing other-cause death from every living state. See R/abcd_engine.R.
#   SURGERY is a clinical SERVICE flow (fraction of C/D receiving surgery each
#   cycle), NOT a health state: it modifies C->D and D->RHD-death risk only.
#   Secondary antibiotic prophylaxis (SAP) lowers PROGRESSION hazards only (block H)
#   for confirmed, SAP-eligible people on BPG (block I). It does not reduce incidence
#   and has no direct mortality effect.
#
# PARAMETER TAGS: [PAPER] article main text; [LIT] literature value; [CALIBRATE]
#   tune to setting/appendix. All economic/monetary values live ONLY in block J
#   (consumed by 08); scripts 01-07 contain no monetary values.
################################################################################

rm(list = ls())

# NOTE: load data.table BEFORE wpp2024 (loaded inside 02) to avoid the Windows
# OpenMP-runtime segfault. Individual scripts also load what they need.
suppressPackageStartupMessages({
  library(data.table)   # <- first
  library(dplyr)
  library(tidyr)
  library(foreach)
  library(doParallel)
  library(parallel)
})

#===============================================================================
# A. COUNTRY + PATHS
#    COUNTRY is the SINGLE switch (Indonesia | Uganda). data/ and output/ are
#    COUNTRY-scoped so runs never overwrite each other; data-raw/ is SHARED (the
#    raw GBD/WPP CSVs contain every country and 01 filters to COUNTRY).
#===============================================================================
COUNTRY <- "Indonesia"           # <<< the ONE line to change to run another country
ISO3    <- switch(COUNTRY,
                  Indonesia = "IDN",
                  Uganda    = "UGA",
                  stop("Unknown COUNTRY '", COUNTRY,
                       "'. Add its ISO3 here and its block-B/J values below."))

wd      <- "C:/Users/wrgar/OneDrive - UW/02Work/WorldBank-Indonesia/uw-wb-indonesia-rhd/"
wd_code <- paste0(wd, "scripts/")
wd_raw  <- paste0(wd, "data-raw/")                 # raw GBD/WPP inputs (SHARED; not on GitHub)
wd_data <- paste0(wd, "data/",   COUNTRY, "/")     # processed inputs + calibrated TPs + baseline state
wd_outp <- paste0(wd, "output/", COUNTRY, "/")     # model outputs + tables
dir.create(wd_data, recursive = TRUE, showWarnings = FALSE)
dir.create(wd_outp, recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# B. LOCATIONS + COUNTRY-SPECIFIC SANITY BANDS
#===============================================================================
LOCATION  <- COUNTRY            # single location used by 01/03/04 (rate/TP build)
LOCATIONS <- c(COUNTRY)         # location list run by 05/06 (parallelised in 06)

# Expected national population ~2020 — 02's population sanity band (do NOT loosen;
# set per country). Indonesia ~275M, Uganda ~45M (UN WPP2024).
POP2020_LO <- switch(COUNTRY, Indonesia = 2.0e8, Uganda = 3.5e7)
POP2020_HI <- switch(COUNTRY, Indonesia = 3.5e8, Uganda = 5.5e7)

# Base-year RHD input-sanity bands (order-of-magnitude gates in 03 & 05), set per
# COUNTRY (do NOT loosen globally). Prevalence is high in both settings, but the
# absolute RHD DEATH count differs sharply: GBD 2023 reports ~393 RHD deaths for
# Uganda (young population; the crude base-year count is ~430), well below
# Indonesia's ~thousands — so Uganda's death-band LOWER bound is country-specific.
RHD_PREV_LO  <- switch(COUNTRY, Indonesia = 1e5, Uganda = 1e5)
RHD_PREV_HI  <- switch(COUNTRY, Indonesia = 1e7, Uganda = 1e7)
RHD_DEATH_LO <- switch(COUNTRY, Indonesia = 1e3, Uganda = 1e2)
RHD_DEATH_HI <- switch(COUNTRY, Indonesia = 1e5, Uganda = 1e5)
# Base-year ALL-CAUSE death band (06's order-of-magnitude gate). Indonesia ~2.5M;
# Uganda ~322k (GBD 2023 all-cause; smaller, younger population).
ALLC_DEATH_LO <- switch(COUNTRY, Indonesia = 5e5, Uganda = 2e5)
ALLC_DEATH_HI <- switch(COUNTRY, Indonesia = 4e6, Uganda = 6e5)

#===============================================================================
# C. HORIZON / YEAR WINDOWS
#    The pipeline distinguishes SEVEN conceptually-separate time references. Only
#    the ANALYSIS period (1) is the model run / reporting horizon; the others are
#    independent inputs and may differ from it. Change the analysis period HERE
#    and it propagates through demography, disease inputs, baseline, model run,
#    outputs, economics, plots, tables, text and report — nothing downstream
#    hard-codes a horizon or a 2100 endpoint.
#      (1) ANALYSIS period          — model + every output/report object (this block)
#      (2) OBSERVED-data period     — WPP2024 estimates (feeds calibration)   [OBS_YEARS]
#      (3) PROJECTION-demography    — medium-variant population backbone       [PROJ_YEARS]
#      (4) CALIBRATION window       — in-sample fit years (04/04b)             [CAL_YEAR_*]
#      (5) GBD rate reference year  — latest GBD rate year, held forward (03)  [RATE_BASE_YEAR]
#      (6) GDP-per-capita ref year  — year gdp_pc_base refers to (block J)     [gdp_pc_base_year]
#          GNI / VSL income anchor  — year gni_pc_base refers to (block J)     [gni_pc_base_year]
#      (7) discount reference year  — economic discounting anchor (block J)    [discount_base_year]
#    (The intervention scale-up ramp, block D, is an 8th, also independent.)
#-------------------------------------------------------------------------------
# (1) ANALYSIS PERIOD — the model run / reporting horizon. EDIT THESE TWO LINES.
#     Every model, output, economic and report object covers EXACTLY these years.
ANALYSIS_YEAR_START <- 2026L
ANALYSIS_YEAR_END   <- 2100L
ANALYSIS_YEARS      <- ANALYSIS_YEAR_START:ANALYSIS_YEAR_END

# (2) OBSERVED-data period (02 demography ESTIMATES; also feeds the 04/04b
#     calibration). Independent of the analysis window: WPP2024 reports estimates
#     through its jump-off year (2024) and must span the calibration window below.
OBS_YEARS <- 1990:2024

# (3) PROJECTION-demography period (02 medium-variant backbone the model draws
#     population from). CONTIGUOUS with OBS_YEARS (starts the year after observed
#     ends, or 02's join-continuity check fails) and need only reach the analysis
#     end — NO fixed 2100 endpoint. Widen PROJ_YEAR_END only to keep population
#     available beyond the analysis horizon.
PROJ_YEAR_START <- max(OBS_YEARS) + 1L
PROJ_YEAR_END   <- ANALYSIS_YEAR_END
PROJ_YEARS      <- PROJ_YEAR_START:PROJ_YEAR_END

# (4) CALIBRATION in-sample window (04/04b) — a fit period, NOT the analysis
#     period. Must fall within the observed-data period.
CAL_YEAR_START <- 2000L
CAL_YEAR_END   <- 2019L

# (5) GBD rate reference year — latest GBD rate year, held forward by 03.
RATE_BASE_YEAR <- 2023L

# --- consistency checks on the year windows (fail fast, before sourcing 01-08) --
stopifnot(
  "ANALYSIS_YEAR_START <= ANALYSIS_YEAR_END" = ANALYSIS_YEAR_START <= ANALYSIS_YEAR_END,
  "OBS_YEARS ascending"                       = min(OBS_YEARS)  <  max(OBS_YEARS),
  "PROJ_YEARS ascending"                      = min(PROJ_YEARS) <= max(PROJ_YEARS),
  "observed & projection contiguous"          = min(PROJ_YEARS) == max(OBS_YEARS) + 1L,
  "calibration window within observed period" =
    CAL_YEAR_START >= min(OBS_YEARS) && CAL_YEAR_END <= max(OBS_YEARS) &&
    CAL_YEAR_START <= CAL_YEAR_END,
  "analysis period within projection backbone" = all(ANALYSIS_YEARS %in% PROJ_YEARS)
)
message(sprintf(
  "Year windows | analysis %d-%d | observed %d-%d | projection %d-%d | calibration %d-%d | GBD rate ref %d",
  ANALYSIS_YEAR_START, ANALYSIS_YEAR_END, min(OBS_YEARS), max(OBS_YEARS),
  min(PROJ_YEARS), max(PROJ_YEARS), CAL_YEAR_START, CAL_YEAR_END, RATE_BASE_YEAR))

#===============================================================================
# D. INTERVENTION SCALE-UP RAMP WINDOW  (03/05) — independent of the analysis
#    period (a policy assumption). Defaults coincide with the analysis bounds but
#    may differ (e.g. a ramp that finishes before or after the analysis horizon).
#===============================================================================
ramp_start <- 2026L              # first year cascade coverage begins rising (baseline year)
ramp_end   <- 2050L              # year target cascade coverage is reached (held thereafter)

#===============================================================================
# E. INCIDENCE SECULAR TREND  (03/05)
#===============================================================================
incidence_trend <- switch(COUNTRY,
                          Indonesia = 0.985,
                          Uganda    = 0.985)  # [CALIBRATE] ~15%/decade RHD incidence decline
                                              # [TODO: confirm for Uganda] trend may differ from Indonesia's

#===============================================================================
# F. CALIBRATION SETTINGS  (04) — Layer-1 aggregate, hazard-based, prior-anchored
#    (REPLACES the former pure random search on IR/CF vs GBD all-cause deaths.)
#    Layer 1 now calibrates, per location x sex, a low-dimensional set of
#    log-scale multipliers — one INCIDENCE correction (alpha) and one aggregate
#    RHD-MORTALITY correction (beta) per broad AGE BAND — for a well-sick-dead
#    proxy with hazard competing risks, against GBD RHD-SPECIFIC incidence,
#    prevalence and DEATHS (all-cause deaths are retained for VALIDATION only).
#    The loss is log-scale + priors + 2nd-difference smoothness; the optimiser is
#    a bounded multi-start local search (see CALIB_OPTIMIZER). Layer 1 is a PROXY:
#    its alpha are STARTING VALUES for 04b and its CF is diagnostic only (it never
#    overrides the fixed D-only mortality); the final incidence, natural history,
#    mortality and initial state come from Layer 2 (04b) on the production engine.
#===============================================================================
run_calibration_par <- switch(COUNTRY, Indonesia = TRUE, Uganda = TRUE)
SEED                <- switch(COUNTRY, Indonesia = 42L,  Uganda = 42L)  # master RNG

## --- objective weights (log-scale, normalized loss terms) -------------------
# Prevalence-led (5 Nov 2026; replaces the former death-led 1/1/2 weighting).
W_INC   <- switch(COUNTRY, Indonesia = 1,   Uganda = 1)     # weight on RHD INCIDENCE fit
W_PREV  <- switch(COUNTRY, Indonesia = 2.5, Uganda = 2.5)   # weight on RHD PREVALENCE fit
W_DEATH <- switch(COUNTRY, Indonesia = 1,   Uganda = 1)     # weight on RHD-SPECIFIC DEATHS fit
# Inverse-variance weight the log-residuals using GBD 95% UIs when present
# (01 now retains upper/lower); FALSE => unit weights.
USE_IV_WEIGHTS <- switch(COUNTRY, Indonesia = TRUE, Uganda = TRUE)
# Drop target cells whose GBD Number is below this floor from the LOSS (they are
# still reported for validation); guards divide-by-tiny in the relative fit.
TARGET_MIN_COUNT <- switch(COUNTRY, Indonesia = 1, Uganda = 1)

## --- incidence calibration mode + broad age bands ---------------------------
# "fixed"    : IR held exactly at the GBD age-sex pattern (alpha = 1, not estimated)
# "anchored" : IR = IR_GBD x exp(f_band(a)), f low-dim & prior-anchored to 0 (DEFAULT)
# "free"     : same parameterisation but weak prior (alpha free to move)
INCIDENCE_CALIBRATION_MODE <- switch(COUNTRY, Indonesia = "anchored", Uganda = "anchored")
# Lower edges of the broad calibration age bands (last band is open-ended).
# 0-14 / 15-24 / 25-44 / 45-64 / 65-79 / 80+ — one alpha & one beta per band x sex
# (65+ split at 80; used by the dynamic band logic of 04 and 04b).
CALIB_AGE_BANDS <- c(0L, 15L, 25L, 45L, 65L, 80L)

## --- priors + smoothness (all on the log-multiplier scale) ------------------
SIGMA_ALPHA   <- switch(COUNTRY, Indonesia = 0.75, Uganda = 0.75)  # incidence prior sd (GENEROUS:
                                                                   #   let incidence breathe so the
                                                                   #   GBD-incidence >> prevalence
                                                                   #   tension does NOT re-enter)
SIGMA_BETA    <- switch(COUNTRY, Indonesia = 0.50, Uganda = 0.50)  # RHD-mortality prior sd
LAMBDA_PRIOR  <- switch(COUNTRY, Indonesia = 1.0,  Uganda = 1.0)   # weight on the prior penalty
LAMBDA_SMOOTH <- switch(COUNTRY, Indonesia = 1.0,  Uganda = 1.0)   # weight on 2nd-diff smoothness
MULT_LO       <- 0.2    # multiplier lower bound (log bounds = log(MULT_LO/HI))
MULT_HI       <- 5.0    # multiplier upper bound

## --- optimiser --------------------------------------------------------------
# "multistart" : bounded multi-start L-BFGS-B (base R optim; DEFAULT, no extra deps)
# "nloptr"     : global CRS2 -> L-BFGS-B polish (needs the installed nloptr)
# "deoptim"    : DEoptim global -> L-BFGS-B polish (only if DEoptim is installed)
CALIB_OPTIMIZER <- "multistart"
N_STARTS        <- 8L     # reproducible random starts (+ the baseline start-0), multistart/deoptim
CONVERGE_TOL    <- 1e-8   # optimiser convergence tolerance

#===============================================================================
# G. DEFAULT A/B/C/D NATURAL HISTORY — POPULATION BURDEN + FIRST DIAGNOSIS (03/04b/05/06)
#    Stages: A minimal/early, B mild, C advanced w/o complications, D advanced WITH
#    clinical complications (HF). Values marked [PROVISIONAL] are reference modelling
#    inputs, NOT validated untreated WHF-stage transition estimates. Nothing is read
#    from any workbook at run time.
#
#    Rules depend on ATTAINED AGE only (R/abcd_engine.R), never on when diagnosis
#    occurred: child rules at age <= nh_child_last_age (20), adult rules at > 20. This
#    is a declared MODELLING CONVENTION; it is separate from the GBD reporting split
#    0-19 / 20+ used for the calibration components (a 21+ diagnostic is reported) and
#    is not asserted to verify a clinical guideline.
#    At the 20 -> 21 ageing step surviving A becomes ADULT PRIOR-A lesion history:
#    excluded from classifiable RHD prevalence, not SAP-eligible, people conserved (not
#    remission, not dropped). Prior-A -> B (shared A->B parameter; an extrapolation)
#    is tracked as RE-ENTRY into classifiable disease, separately from first entry and
#    from first diagnosis. Normalised people keep an explicit history stock, face
#    background mortality only and have no recurrence (simplified assumption;
#    normalisation is not proven anatomical healing).
#    All inputs are ANNUAL PROBABILITIES, converted ONCE to hazards h = -log1p(-p).
#===============================================================================
mortality_structure <- "D_only"   # the implemented default: RHD death only from D (03 rejects anything else)

nh_child_last_age <- 20L          # [MODELLING CONVENTION] child rules at attained age <= 20; adult rules at > 20

# Annual transition PROBABILITIES = 04b STARTING VALUES and PRIOR CENTRES. The six
# non-zero ones are calibrated (shared across sexes) within TRANS_BOUNDS (block G2).
p_A_to_no_rhd <- 0.311   # [PROVISIONAL] A -> normalised history, age <= 20 (adult: 0)             [calibrated]
p_A_to_B      <- 0.054   # [PROVISIONAL] A -> B, age <= 20; the SAME fitted value drives adult
                         #   prior-A -> B (provisional extrapolation)                                [calibrated]
p_B_to_A      <- 0.177   # [PROVISIONAL] B -> A, age <= 20 (adult: structural zero)                [calibrated]
p_B_to_C      <- 0.074   # [PROVISIONAL] B -> C, age <= 20; adult (> 20) = hazard x adult_prog_hazard_mult [calibrated]
p_C_to_B      <- 0.050   # [PROVISIONAL] C -> B, age <= 20 (adult: structural zero)                [calibrated]
p_C_to_D      <- 0.112   # [PROVISIONAL] C -> D, age <= 20; adult (> 20) = hazard x adult_prog_hazard_mult [calibrated]
p_D_to_C      <- 0       # EXACT structural zero (D = clinical complications; no valve-regression rate)

# Fixed inputs (not calibrated):
p_B_to_no_rhd_adult <- 0.02   # [PROVISIONAL extrapolation] adult (> 20) B -> normalised history (child: 0)
# New modelled CLASSIFIABLE disease (first entry from the never-onset pool) enters:
entry_mix_child <- c(A = 1, B = 0, C = 0, D = 0)   # [NEW SIMPLIFYING ASSUMPTION] age <= 20 (fixed)
entry_mix_adult <- c(A = 0, B = 1, C = 0, D = 0)   # [NEW SIMPLIFYING ASSUMPTION] age  > 20 (fixed; no adult A)
#   Not estimates of biological onset age. Advanced first diagnoses arise from the
#   progression of previously undiagnosed disease, not from direct C/D entry.
# PROVISIONAL DIAGNOSIS-COMPOSITION BENCHMARKS (stage at FIRST CONFIRMED DIAGNOSIS by
# age at diagnosis): illustrative reference values, compared as a DIAGNOSTIC only —
# not observed national stage targets, not a fit term, never applied to new disease.
dx_mix_benchmark_child <- c(A = 0.50, B = 0.30, C = 0.15, D = 0.05)   # age <= 20 at first diagnosis
dx_mix_benchmark_adult <- c(A = 0,    B = 0.60, C = 0.30, D = 0.10)   # age  > 20 at first diagnosis

# Stage-D annual RHD-death probability (FIXED; applied to EVERY disease / diagnosis
# history — the former late-entry D-death hazard multiplier of 2 is removed, i.e.
# effective value 1; not fitted and never lowered to manufacture fit).
p_rhd_death_D <- switch(COUNTRY,
                        Indonesia = 0.119,   # [PROVISIONAL reference input; dev/TP Stage D to death_RHDparameters.xlsx col K]
                        Uganda    = 0.0597)  # [PROVISIONAL] Uganda-specific input (not replaced by Indonesia's value)
residual_rhd_mortality_mult <- 1     # fixed; no second mortality scale
# COMMON ADULT PROGRESSION HAZARD MULTIPLIER (fixed; not fitted). At attained age > 20
# (nh_child_last_age; age 20 itself uses the childhood probabilities) adult B->C and C->D
# use p_adult = 1 - (1 - p_child)^m, i.e. the adult HAZARD is m x the childhood hazard
# (never p x m). Applies to undiagnosed, diagnosed and on-BPG disease regardless of age at
# onset or diagnosis, in calibration (04b: burn-in, window, bridge) and projection (06).
# The calibrated (and prior-centred / bounded) parameters remain the CHILDHOOD values;
# adult prior-A -> B (= p_A_to_B), adult B normalisation and D mortality are NOT scaled.
# Adopted 6 Nov 2026 as the production reference for both countries after the Indonesia
# scratch test dev/experiments/2026-11_adult_prog_mult_048_pipeline (REPORT.md): total
# loss -15% vs m = 1, but adult RHD deaths remain ~3-4x GBD (a partial improvement).
adult_prog_hazard_mult <- switch(COUNTRY,
                                 Indonesia = 0.48,   # [PROVISIONAL] reference value (formerly the late-onset 0.48 progression multiplier)
                                 Uganda    = 0.48)   # [PROVISIONAL] same reference value; [TODO: confirm for Uganda]
stopifnot("adult_prog_hazard_mult must be one value > 0" =
            is.numeric(adult_prog_hazard_mult) && length(adult_prog_hazard_mult) == 1L && adult_prog_hazard_mult > 0)
# Realised annual D fatality can be smaller than the input
# because of competing other-cause mortality and the reference surgery service.
# For Uganda every shared input above except p_rhd_death_D is an EXTRAPOLATION of the
# illustrative reference values, not a Uganda-specific estimate.

# Diagnostic only: GBD "HF due to RHD" proxy D share (not a validated stage target).
d_share_proxy_gbd_hf <- switch(COUNTRY, Indonesia = 0.0321, Uganda = 0.00863)

#===============================================================================
# G2. STRUCTURAL (LAYER-2) CALIBRATION  (04b)
#     Runs the REAL engine (R/abcd_engine.R) inside the objective — the same burn-in,
#     care update, ageing, prior-A reclassification and observation as production —
#     and jointly estimates the band x sex incidence multipliers and the six non-zero
#     transition probabilities against GBD RHD prevalence (TOTAL classifiable disease,
#     diagnosed + undiagnosed; prior-A and normalised history excluded) + RHD deaths
#     (diagnosed + undiagnosed D). Entry mixes, adult normalisation, the adult
#     progression hazard multiplier and D mortality are FIXED (block G).
#===============================================================================
run_structural_calibration <- switch(COUNTRY, Indonesia = TRUE, Uganda = TRUE)
STRUCT_MATCH_FROM  <- 2010L          # first GBD year SCORED (2000-2009 = initialisation / validation years)
# Zero-start burn-in under first-window-year inputs and the provisional initial care
# shares (block I) generates prevalence from incidence + natural history; 04b CHECKS
# convergence of stage / prior-A / normalised totals, age profiles and deaths vs a 2x run.
STRUCT_BURN_YEARS  <- 120L
# Burn-in INITIALISATION reconciliation (R/abcd_engine.R abcd_burn_in; 04b, recorded in the bundle and
# checked by 05). The burn-in holds the first-window-year (2000) cross-section constant for 120 cycles:
#   "per_capita"  (DEFAULT): after each burn-in ageing step every population count stock (classifiable
#                 disease by care state, prior-A, normalised, never-onset) is rescaled to that population, giving an
#                 equilibrium in per-capita rates. The adjustments are INITIALISATION BALANCING (recorded
#                 separately, never deaths or onset). Adopted 4 Nov 2026: the previous rule built older cohorts'
#                 RHD counts from today's larger young cohorts (x1.5 at age 40, x2.7 at 60, x7 at 80 for Indonesia
#                 2000) and the calibration compensated with fast progression to D; a scratch recalibration
#                 lowered L 2.52 -> 1.50 with all three data components improved, but remained OUTSIDE the
#                 diagnostic tolerances (prevalence 0.73x, adult deaths 2.7-3.4x).
#   "no_rhd_only" : the previous rule (only No-RHD pools absorb the residual during the burn-in), kept to
#                 reproduce engine-v2 results.
# Observed years, the bridge and production ALWAYS use the No-RHD-only convention (disease stocks never rescaled).
STRUCT_BURNIN_RECONCILIATION <- "per_capita"
# Count-scale epsilon inside log() for the structural MSLE terms: ONE case
# (log(x + 1)); cells with GBD count < TARGET_MIN_COUNT are masked from the cell loss
# (still reported, and aggregated childhood fit is always shown).
STRUCT_EPS_COUNT   <- 1

# priors: logit-sd on the six calibrated Y transitions (literature RANGES in
# TRANS_BOUNDS carry the clinical constraint; the data pick the value within them).
SIGMA_TRANS  <- 2.00
# PRIOR STRENGTH: data terms are MEANS over ~400 target cells while prior penalties
# are SUMS over the parameters; 0.05 ~ treats the target grid as ~20 effective
# independent observations (04b exports the prior/data ratio).
LAMBDA_PRIOR_STRUCT <- 0.05

# --- PREVALENCE-LED JOINT OBJECTIVE (04b) -------------------------------------
#   L = wP * L_prev + wD * L_death + wI * MSLE(incidence)
#       + prior penalty + smoothness penalty + optional stage-target error
#   L_prev  = p * MSLE(prevalence, ages 0-19) + (1 - p) * MSLE(prevalence, ages 20+)
#   L_death = s * MSLE(deaths, ages 0-19)     + (1 - s) * MSLE(deaths, ages 20+)
# Each MSLE is a NORMALISED MEAN over its own usable cells, so a larger target grid or
# very small death cells cannot dominate; raw and weighted components are exported.
# Prevalence (incl. its age pattern) is the PRINCIPAL target; RHD deaths are fit
# jointly. Operational choices for this implementation, not published weights.
STRUCT_PREV_DEATH_RATIO  <- 2.5       # wP / wD
STRUCT_DATA_WEIGHT_TOTAL <- 3.5       # wP + wD (2.5:1 -> wP = 2.5, wD = 1)
W_PREV_STRUCT  <- STRUCT_DATA_WEIGHT_TOTAL * STRUCT_PREV_DEATH_RATIO / (1 + STRUCT_PREV_DEATH_RATIO)
W_DEATH_STRUCT <- STRUCT_DATA_WEIGHT_TOTAL / (1 + STRUCT_PREV_DEATH_RATIO)
W_INC_STRUCT   <- 0                   # incidence = anchored input + VALIDATION output (not fit directly)
# Childhood shares p (prevalence) and s (deaths): each component is a normalised mean
# log-count error over its own usable target cells (04b stops if any has no usable cell).
STRUCT_PREV_CHILD_SHARE    <- 0.50    # equal weight on prevalence ages 0-19 and 20+ (explicit adult contribution)
STRUCT_DEATH_CHILD_SHARE   <- 0.50
STRUCT_DEATH_CHILD_AGE_MAX <- 19L     # childhood component = GBD groups wholly below age 20 (ages 0-19; GBD reporting split)
# 21+ DIAGNOSTIC (reported, not fitted): isolates the <= 20 / > 20 modelling convention
# from the 0-19 / 20+ GBD split; target at single age 20 = GBD 20-24 rate x WPP age-20 population.
STRUCT_DIAG_ADULT_CONVENTION_AGE_MIN <- 21L
# Incidence corrections estimated JOINTLY: one log-multiplier per broad band
# (CALIB_AGE_BANDS) x sex on the GBD incidence pattern, weak prior + smoothness.
SIGMA_ALPHA_STRUCT  <- 1.0
LAMBDA_SMOOTH_STRUCT <- 0.5
STRUCT_ALPHA_LO <- 0.2; STRUCT_ALPHA_HI <- 5.0
USE_IV_WEIGHTS_STRUCT <- FALSE        # GBD UIs are not independent measurement errors across age/year
STRUCT_IV_WEIGHT_CAP  <- 5
STRUCT_AGE_WEIGHTS <- NULL            # documented band weights (e.g. c("0-14" = 1.5)); NULL = equal
# Target basis (built in 04): "rate_x_pop" = GBD Rate x the model's WPP population of
# the same cell (published GBD Numbers are reported side by side); or "number".
STRUCT_TARGET_BASIS <- "rate_x_pop"
STRUCT_PREV_OBS <- "mid"              # prevalence observation: mean of start and end (pre-ageing) stock
STRUCT_PREV_FIRST <- TRUE             # prevalence-only preliminary fit supplies a START only
# Fit tolerances on |log(model/target)| across EVERY scored year (operational checkpoint
# criteria, not external validation): prevalence tolerance on all-age and adult (20+)
# prevalence; death tolerance on all-age, childhood (0-19) and adult (20+) deaths.
STRUCT_FIT_TOL_PREV  <- 0.15
STRUCT_FIT_TOL_DEATH <- 0.25
# Age floor on the RHD-DEATH target cells: 0 = every published cell enters the loss
# (childhood deaths are no longer excluded).
STRUCT_DEATH_TARGET_AGE_MIN <- 0L

# Annual-probability bounds of the calibrated transitions (literature ranges, widened
# modestly). D->C is a fixed structural zero and never enters a transform or the optimiser.
TRANS_BOUNDS <- list(
  p_A_to_no_rhd = c(0.005, 0.60),
  p_A_to_B      = c(0.001, 0.25),
  p_B_to_A      = c(0.002, 0.40),
  p_B_to_C      = c(0.001, 0.20),
  p_C_to_B      = c(0.001, 0.15),
  p_C_to_D      = c(0.001, 0.30),
  p_D_to_C      = c(0, 0))

STRUCT_N_STARTS      <- 2L           # reproducible random starts (+ prior-centre, Layer-1 and prevalence-first starts)
STRUCT_PARALLEL_STARTS <- TRUE       # compute only: run the independent starts on parallel workers (<= MAX_CORES; identical results)
STRUCT_NEAROPT_DELTA <- 0.05         # near-optimal set: L <= L_min * (1 + delta)
# Forward-project the calibrated end-of-window state (all stage, diagnosis, BPG,
# prior-A, normalised and screening histories) through CAL_YEAR_END+1 ..
# ANALYSIS_YEAR_START-1 under REFERENCE care; 05 seeds 06 with it (no reseeding).
STRUCT_BRIDGE_TO_ANALYSIS_START <- TRUE

#===============================================================================
# H. INTERVENTION EFFECTS  (03)
#===============================================================================
# SAP acts through PROGRESSION hazards only, for confirmed SAP-ELIGIBLE people on BPG
# (block I): h = h(0) x [1 - adh x (1 - HR)], adh = good + rel x partial adherence.
# Screening alone has no treatment or mortality effect. No direct mortality effect.
sap_mechanism    <- "progression"     # the implemented default (03 rejects anything else)
sap_effect_scale <- "hazard_ratio"
# Child HRs from the GOAL trial (Uganda, children 5-17 with latent RHD, 24-month
# follow-up; control vs prophylaxis 2-year risks): A->B 0.075 vs 0.003; B->C 0.119 vs
# 0.025 -> HR = log(1 - p_sap) / log(1 - p_ctrl) (a hazard ratio; NOT annualised again).
# PROVISIONAL stage mapping (borderline ~ A, mild definite ~ B). Explicit conservative
# choices for this default: the child trial effect is NOT extrapolated to adults.
goal_ctrl_A_to_B <- 0.075; goal_sap_A_to_B <- 0.003
goal_ctrl_B_to_C <- 0.119; goal_sap_B_to_C <- 0.025
hr_from_arm_risks <- function(p_sap, p_ctrl) log1p(-p_sap) / log1p(-p_ctrl)
sap_hr_A_to_B_child <- hr_from_arm_risks(goal_sap_A_to_B, goal_ctrl_A_to_B)   # ~0.03853835 (age <= 20)
sap_hr_B_to_C_child <- hr_from_arm_risks(goal_sap_B_to_C, goal_ctrl_B_to_C)   # ~0.19982855 (age <= 20)
sap_hr_B_to_C_adult <- 1              # age > 20: adult B on BPG has programme cost but NO quantified progression benefit
sap_hr_prior_A_to_B_adult <- 1        # adult prior-A: not SAP-eligible
sap_hr_C_to_D <- 1                    # C is not SAP-eligible in this simplified default
sap_hr_D_to_rhd_death <- 1            # no direct SAP mortality benefit (03 enforces 1)
# Calibrated transitions are UNTREATED natural history; reference care (block I) is
# applied inside the engine in calibration, bridge and projection (applied once).
baseline_care_interpretation <- "untreated_parameters_reference_care_applied_in_engine"

# Surgery effects (C and D only), scaled by surgery reach = fraction requiring x coverage.
eff_surgery_C_to_D        <- 0.85  # [LIT]/[CALIBRATE] surgery RRR on C -> D progression
eff_surgery_D_to_rhd_death <- 0.85 # [LIT]/[CALIBRATE] surgery RRR on D -> RHD death

#===============================================================================
# I. SCREENING, CONFIRMATION, CARE + SURGERY SERVICE  (03/05)
#    ONE echo screen per person across ages 5-39 inclusive among people WITHOUT a
#    prior confirmed diagnosis (never-screened history tracked by the engine). True
#    screen positives -> confirmatory echo -> FIRST CONFIRMED DIAGNOSIS (entry into
#    the diagnosed / care cohort, any stage A-D). Screening detects existing disease;
#    confirmation failure leaves disease undiagnosed. Clinical detection (confirmed
#    clinical diagnosis) continues at all ages. BPG for confirmed SAP-ELIGIBLE people
#    (A at age <= 20, B); C/D and adult prior-A have no SAP programme volume.
#    Per arm: baseline -> target, reference held flat, scale-up ramped linearly over
#    ramp_start..ramp_end (block D) and held after. Historical initialisation,
#    calibration and the bridge use the REFERENCE values throughout.
#    Screening target c = probability of one screen across the window; annual
#    probability among never-screened q = 1 - (1 - c)^(1/35), derived from the CURRENT
#    year's target during the ramp (cumulative realised cohort reach during a changing
#    ramp need not equal the current target).
#    Uganda: baseline screening/diagnosed/on-BPG values retain the previous Uganda-
#    specific cascade inputs; every other care value is a shared [ASSUMPTION].
#===============================================================================
screen_age_lo <- 5L; screen_age_hi <- 39L        # one screen across ages 5-39 inclusive (35 years)
screen_sensitivity  <- 0.90                      # [PROVISIONAL] screening sensitivity (true positives among screened disease)
confirm_completion  <- 0.90                      # [PROVISIONAL, new] confirmatory echo completed among true screen positives
confirm_sensitivity <- 1.00                      # [SIMPLIFYING] confirmation sensitivity
screen_specificity  <- 1                         # [SIMPLIFYING, declared] perfect specificity: no false positives modelled
#   No confirmation unit cost exists in block J: confirmation is UNCOSTED (flagged in 08).
sap_eligibility <- c("A", "B")                   # confirmed A (age <= 20) and B; the engine implements exactly this
partial_adherence_relative_protection <- 0.50    # relative protection under partial adherence
# baseline / target per arm: c(baseline, target)
screen_base_c <- switch(COUNTRY, Indonesia = 0.05, Uganda = 0.01)   # Uganda: previous country-specific screen_base
care_inputs <- list(
  screen_window_prob = list(ref = c(screen_base_c, screen_base_c), up = c(screen_base_c, 0.80)),
  clinical_detection = list(ref = c(0.02, 0.02), up = c(0.02, 0.02)),   # annual confirmed clinical diagnosis, among undiagnosed
  bpg_initiation     = list(ref = c(0.50, 0.50), up = c(0.50, 0.90)),   # among diagnosed untreated
  bpg_retention      = list(ref = c(0.90, 0.90), up = c(0.90, 0.95)),   # annual
  adherence_good     = list(ref = c(0.60, 0.60), up = c(0.60, 0.85)),   # share of those on BPG
  adherence_partial  = list(ref = c(0.20, 0.20), up = c(0.20, 0.10)))   # share of those on BPG
# Provisional initial care state, set ONCE at the first observed-history year
# (CAL_YEAR_START) and carried forward (never reset in 2026); also the static
# pre-history care labels held during the burn-in. Diagnosed share applies to every
# stage stock; the on-BPG share is allocated to the SAP-ELIGIBLE stage stocks only
# (A at age <= 20, B), so the realised on-BPG share of all classifiable disease is lower.
care_init_diagnosed_share     <- switch(COUNTRY, Indonesia = 0.05, Uganda = 0.009)   # Uganda: previous diagnosis_base
care_init_on_bpg_share        <- switch(COUNTRY, Indonesia = 0.04, Uganda = 0.0075)  # Uganda: previous treatment_base
care_init_ever_screened_share <- 0

# SURGERY as a clinical service (NOT a health state), held equal in both arms.
frac_C_requiring_surgery <- switch(COUNTRY,
                                   Indonesia = 0.03,
                                   Uganda    = 0.1)
frac_D_requiring_surgery <- switch(COUNTRY,
                                   Indonesia = 0.20,
                                   Uganda    = 0.4)
cov_surgery_ref <- switch(COUNTRY,
                          Indonesia = 0.05,
                          Uganda    = 0.03)
cov_surgery_scale_up <- switch(COUNTRY,
                               Indonesia = 0.05,
                               Uganda    = 0.03)

#===============================================================================
# J. ECONOMIC PARAMETERS  (08 ONLY) — the only monetary values in the pipeline
#===============================================================================
# disc_rate        <- 0.03         # [PAPER] annual discount rate (costs & benefits)
# gdp_pc_base      <- switch(COUNTRY,   # [LIT] GDP per capita, base-year US$ (World Bank,
#                            Indonesia = 4150,   #   current US$, 2019: IDN ~4,107; UGA ~822)
#                            Uganda    = 822)    #   [TODO: confirm] if a non-WB source is preferred
# gdp_pc_base_year <- 2019L        # [LIT] year gdp_pc_base refers to (same basis for both)
# gdp_growth       <- 0.03         # [CALIBRATE] real per-capita GDP growth
# vsl_mult         <- 30           # [PAPER] value of a statistical life = 30 x GDP per capita
# dalys_per_death  <- 30           # [CALIBRATE] undiscounted DALYs per RHD death
# discount_base_year <- ANALYSIS_YEAR_START  # discount to this year (default = first analysis year)
# 
# # unit costs (base-year US$) -- [PAPER]/[LIT]
# cost_screen_per_person <- 1.10   # [PAPER] echo screening cost per person screened
# cost_sap_per_year      <- 110    # [PAPER] annual optimal secondary-prevention cost per person
# cost_surgery           <- 9000   # [LIT] valve surgery + first-year post-op cost per surgery

disc_rate <- switch(COUNTRY,
                    Indonesia = 0.03,
                    Uganda    = 0.03)

# BENEFIT VALUATION (08): VSL ALONE monetises RHD deaths averted (the former additive
# one-year GDP-per-capita "productivity" benefit is RETIRED). VSL is anchored on GNI
# PER CAPITA and projected with the GDP-per-capita growth path (unit income elasticity
# over time):
#     VSL_b = vsl_mult x gni_pc_base ;  VSL_t = VSL_b x GDPpc_t / GDPpc_b = VSL_b x (1 + gdp_growth)^(t - b)
# with b = gni_pc_base_year (the GNI/VSL income anchor; distinct from the analysis,
# calibration, GBD-rate and discount reference years).
# GNI per capita: World Bank WDI indicator NY.GNP.PCAP.CD, "GNI per capita, Atlas
# method (current US$)"; values for 2019 (WDI last updated 2026-10-08; retrieved
# 2026-10-22 from api.worldbank.org). Current US$ of 2019 = the same 2019-US$ basis the
# economic block already used (gdp_pc_base_year), so VSL is in constant 2019 US$ and is
# grown with REAL per-capita growth (no nominal compounding). [TODO: confirm] unit
# costs below carry no recorded price year; they are assumed to be on this basis.
# (For reference, WDI 2024 values are 4,920 Indonesia / 1,010 Uganda; moving the anchor
# to 2024 would require re-pricing the unit costs to 2024 US$ as well.)
gni_pc_base <- switch(COUNTRY,
                      Indonesia = 4020,     # [WDI NY.GNP.PCAP.CD 2019] US$, Atlas method
                      Uganda    = 810)      # [WDI NY.GNP.PCAP.CD 2019] US$, Atlas method
gni_pc_base_year <- switch(COUNTRY,
                           Indonesia = 2019L,
                           Uganda    = 2019L)
gni_pc_source <- "World Bank WDI NY.GNP.PCAP.CD (GNI per capita, Atlas method, current US$), 2019; WDI updated 2026-10-08"

# GDP-per-capita PATH: used ONLY as the growth index that carries VSL forward from the
# GNI anchor year (08 aligns it to gni_pc_base_year before taking the ratio, so its
# LEVEL cancels). gdp_growth is a constant REAL per-capita growth SCENARIO ASSUMPTION,
# not an externally sourced forecast. [TODO: confirm] Uganda's 1340 does not match WDI
# 2019 GDP per capita (822 current US$); harmless for valuation because the level cancels.
gdp_pc_base <- switch(COUNTRY,
                      Indonesia = 4150,
                      Uganda    = 1340)
gdp_pc_base_year <- switch(COUNTRY,
                           Indonesia = 2019L,
                           Uganda    = 2019L)
gdp_growth <- switch(COUNTRY,
                     Indonesia = 0.03,      # [ASSUMPTION] constant real GDP-per-capita growth
                     Uganda    = 0.032)     # [ASSUMPTION] constant real GDP-per-capita growth
# VSL multiple of GNI per capita. [PROVISIONAL] the 30x multiple is RETAINED from the
# previous 30 x GDP-per-capita specification (Coates 2021 used VSL ~30x income) and is
# now applied to GNI per capita; it is not a newly estimated income-elasticity transfer.
vsl_mult <- switch(COUNTRY,
                   Indonesia = 30,
                   Uganda    = 30)
dalys_per_death <- switch(COUNTRY,
                          Indonesia = 30,
                          Uganda    = 30)
# (7) Economic discount reference year — costs & benefits are discounted TO this
#     year. Default = the first analysis year (kept separate from (5) the GBD rate
#     reference year, (6) gdp_pc_base_year and the GNI/VSL anchor gni_pc_base_year,
#     which are their own anchors).
discount_base_year <- ANALYSIS_YEAR_START

cost_screen_per_person <- switch(COUNTRY,
                                 Indonesia = 1.10,
                                 Uganda    = 0.40)
cost_sap_per_year <- switch(COUNTRY,
                            Indonesia = 110,
                            Uganda    = 40)
cost_surgery <- switch(COUNTRY,
                       Indonesia = 9000,
                       Uganda    = 5000)

#===============================================================================
# K. PARALLEL / CORE SETTINGS
#===============================================================================
run_model_par <- TRUE            # parallelise 06 across locations
MAX_CORES     <- 14L             # cap on worker cores (04 and 06)

#===============================================================================
# L. RUN THE PIPELINE  (01 -> 08, in order; each reads the prior persisted output)
#    PIPELINE_RUN_THROUGH = last script sourced: "04b" stops after the structural
#    calibration (checkpoint stage: 05-08 are NOT run); "08" runs the full pipeline.
#    Set to "08" on 5 Nov 2026 by explicit user override although the epidemiological
#    fit checkpoint FAILED: outputs are PRELIMINARY (05 labels the baseline INADEQUATE).
#===============================================================================
PIPELINE_RUN_THROUGH <- "08"
stopifnot("PIPELINE_RUN_THROUGH must be '04b' or '08'" = PIPELINE_RUN_THROUGH %in% c("04b", "08"))
setwd(wd_code)   # scripts use absolute wd_* globals for I/O and here() for R/packages.R

source(paste0(wd_code, "01_prepare_inputs.R"))        # GBD inputs
source(paste0(wd_code, "02_build_demography.R"))      # WPP2024 population backbone
source(paste0(wd_code, "03_build_disease_model.R"))   # disease-model INPUTS (builder)
source(paste0(wd_code, "04_calibration_random_tp.R")) # Layer 1: calibrate IR + CF anchor (2000-2019)
source(paste0(wd_code, "04b_calibrate_structural.R")) # Layer 2: structural A/B/C/D stage calibration
if (PIPELINE_RUN_THROUGH == "08") {
  source(paste0(wd_code, "05_build_baseline.R"))        # assemble initial state (uses 04b stage params)
  source(paste0(wd_code, "06_run_prevention_model.R"))  # run ref + SAP (parallel by location)
  source(paste0(wd_code, "07_make_outputs.R"))          # standard long table
  source(paste0(wd_code, "08_economic_evaluation.R"))   # benefit-cost & cost-effectiveness
  message("\n================ PIPELINE COMPLETE (01 -> 08) ================")
} else message("\n================ PIPELINE STOPPED AFTER 04b (PIPELINE_RUN_THROUGH = '04b'; 05-08 not run) ================")

