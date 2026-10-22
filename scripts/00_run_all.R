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
# DISEASE STRUCTURE: World Heart Federation RHD stages A/B/C/D (living states)
#   No RHD -> A <-> B <-> C -> D -> RHD death, with competing other-cause death
#   from every living state. Incident RHD enters stage A (no separate ARF state).
#   SURGERY is a clinical SERVICE flow (fraction of C/D receiving surgery each
#   cycle), NOT a health state: it modifies C->D and D->RHD-death risk only.
#   Secondary antibiotic prophylaxis (SAP): selectable mechanism (block H). PRIMARY
#   = pathway-specific hazard reductions on supported PROGRESSION transitions (its
#   mortality benefit arises through reduced entry into D); LEGACY comparator =
#   the former 55% direct RHD-mortality RRR. It does not reduce incidence here.
#   RHD mortality structure (block G): PRIMARY = exact D-only; LEGACY = all-stage.
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
#    a bounded multi-start local search (see CALIB_OPTIMIZER). Calibrated IR feeds
#    the inflow into stage A (06); calibrated aggregate CF is the RHD-mortality
#    anchor consumed by the Stage-2 structural stage calibration.
#===============================================================================
run_calibration_par <- switch(COUNTRY, Indonesia = TRUE, Uganda = TRUE)
SEED                <- switch(COUNTRY, Indonesia = 42L,  Uganda = 42L)  # master RNG

## --- objective weights (log-scale, normalized loss terms) -------------------
W_INC   <- switch(COUNTRY, Indonesia = 1, Uganda = 1)   # weight on RHD INCIDENCE fit
W_PREV  <- switch(COUNTRY, Indonesia = 1, Uganda = 1)   # weight on RHD PREVALENCE fit
W_DEATH <- switch(COUNTRY, Indonesia = 2, Uganda = 2)   # weight on RHD-SPECIFIC DEATHS fit
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
# Default 0-14 / 15-24 / 25-44 / 45-64 / 65+ — one alpha & one beta per band x sex.
CALIB_AGE_BANDS <- c(0L, 15L, 25L, 45L, 65L)

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
# G. A/B/C/D NATURAL HISTORY  (03) — annual transition probabilities
#    Stages: A minimal/early, B mild established, C advanced w/o complications,
#    D advanced WITH complications (HF, requiring surgery). Adjacent regression
#    (A->No-RHD, B->A, C->B, D->C) is permitted; set any to 0 for an
#    irreversible specification. Competing OTHER-cause mortality is data-fed
#    (age x sex GBD background, 03/05) and NOT set here. Per-stage RHD-specific
#    death probabilities ARE the lever SAP (and, for D, surgery) act on.
#    Starting values from the standalone A/B/C/D prototype; all [CALIBRATE].
#===============================================================================
# p_A_to_no_rhd <- 0.005   # [CALIBRATE] A regression to No RHD /yr
# p_A_to_B      <- 0.020   # [CALIBRATE] A -> B progression /yr
# p_B_to_A      <- 0.010   # [CALIBRATE] B -> A regression /yr
# p_B_to_C      <- 0.030   # [CALIBRATE] B -> C progression /yr
# p_C_to_B      <- 0.005   # [CALIBRATE] C -> B regression /yr
# p_C_to_D      <- 0.060   # [CALIBRATE] C -> D progression /yr (surgery lowers this)
# p_D_to_C      <- 0.000   # [CALIBRATE] D -> C regression /yr (default 0)

# # annual UNTREATED RHD-specific death probability by stage -- [CALIBRATE]
# p_rhd_death_A <- 0.0005  # [CALIBRATE]
# p_rhd_death_B <- 0.0020  # [CALIBRATE]
# p_rhd_death_C <- 0.0200  # [CALIBRATE]
# p_rhd_death_D <- 0.0800  # [CALIBRATE]


# Temporal hard-coded fix scaling the per-stage RHD-specific death probabilities
# to match observed GBD 2023 mortality: 1/2.5 scales the untreated death probs
# down by a factor of 2.5 to better align with the observed data.
# [TODO: confirm for Uganda] — this multiplier was HAND-TUNED to INDONESIA's
# GBD-2023 RHD mortality and almost certainly needs re-tuning for other countries.


# rhd_mortality_calibration_mult <- 1 / 2.5
# 
# scale_probability <- function(p, multiplier) {
#   1 - (1 - p)^multiplier
# }
# 
# p_rhd_death_A <- scale_probability(0.0005, rhd_mortality_calibration_mult)
# p_rhd_death_B <- scale_probability(0.0020, rhd_mortality_calibration_mult)
# p_rhd_death_C <- scale_probability(0.0200, rhd_mortality_calibration_mult)
# p_rhd_death_D <- scale_probability(0.0800, rhd_mortality_calibration_mult)
# 
# 
# # Stage-D share of prevalent RHD at the seed year. GBD gives only TOTAL RHD
# # prevalence; A/B/C come from Cannon et al (56.5/27.2/16.2, normalised to the
# # non-D share in 03) and D is this RHD-with-heart-failure fraction.
# rhd_d_fraction <- 0.10   # [CALIBRATE] complicated-RHD / RHD-with-HF fraction
#                          # [TODO: confirm for Uganda] may be higher (later presentation)

#===============================================================================
# G. A/B/C/D NATURAL HISTORY  (03) — annual transition probabilities
#===============================================================================
# p_A_to_no_rhd <- switch(COUNTRY,
#                         Indonesia = 0.005,
#                         Uganda    = 0.055) # https://pubmed.ncbi.nlm.nih.gov/28972003/
# p_A_to_B <- switch(COUNTRY,
#                    Indonesia = 0.020,
#                    Uganda    = 0.020)
# p_B_to_A <- switch(COUNTRY,
#                    Indonesia = 0.010,
#                    Uganda    = 0.010)
# p_B_to_C <- switch(COUNTRY,
#                    Indonesia = 0.030,
#                    Uganda    = 0.030)
# p_C_to_B <- switch(COUNTRY,
#                    Indonesia = 0.005,
#                    Uganda    = 0.005)
# p_C_to_D <- switch(COUNTRY,
#                    Indonesia = 0.060,
#                    Uganda    = 0.060)
# p_D_to_C <- switch(COUNTRY,
#                    Indonesia = 0.000,
#                    Uganda    = 0.000)

# A/B/C/D annual transition PRIOR CENTERS (shared clinical priors; 04b calibrates
# them per country against GBD RHD prevalence + RHD deaths using the real engine).
# These are NOT final values and are NOT the former hand-tuned per-country hacks.
# 04b gives transitions LOOSER priors than mortality, so the penalized objective
# reconciles GBD deaths by moving the STOCK/progression, not by suppressing
# mortality (Part 10 excess-death rule). p_A_to_no_rhd (regression) is the lever
# that reconciles high GBD incidence with observed prevalence (mass balance).
p_A_to_no_rhd <- 0.05     # regression A -> No RHD  (loose prior; data-identified)
p_A_to_B      <- 0.020    # A -> B progression
p_B_to_A      <- 0.010    # B -> A regression
p_B_to_C      <- 0.030    # B -> C progression
p_C_to_B      <- 0.005    # C -> B regression
p_C_to_D      <- 0.060    # C -> D progression (surgery lowers this in the engine)
p_D_to_C      <- 0.000    # D -> C regression (fixed at 0 unless calibrate_p_D_to_C = TRUE)

# Per-stage UNTREATED RHD-death PRIOR CENTERS (annual), clinical and severity-
# ordered (A <= B < C < D). The former hand-tuned mortality scalars
# (rhd_mortality_calibration_mult = 1/2.5 Indonesia, 1/8 Uganda) and
# scale_probability() are REMOVED: they hid a stock/progression problem behind
# implausible stage mortality. 04b now calibrates the stock/progression to
# reproduce GBD RHD deaths while these mortalities stay near their clinical priors
# (TIGHT priors). An optional BOUNDED residual multiplier (RESID_MULT_LO/HI) may
# move mortality modestly, but a solution sitting AT that bound is a VALIDATION
# FAILURE (the 1/8 must not reappear under a new name).
p_rhd_death_A <- 0.001    # stage A (minimal/early)      ~0.1%/yr  [clinical prior]  (all_stage structure only)
p_rhd_death_B <- 0.005    # stage B (mild established)    ~0.5%/yr  [clinical prior]  (all_stage structure only)
p_rhd_death_C <- 0.020    # stage C (advanced, no compl.) ~2%/yr    [clinical prior]  (all_stage structure only)
p_rhd_death_D <- 0.080    # stage D (advanced w/ HF)      ~10%/yr   [clinical prior]  (all_stage structure)

# --- MORTALITY STRUCTURE (03/04b/05/06) ---------------------------------------
# "D_only"    (PRIMARY structural experiment): RHD-attributable death occurs ONLY
#             from stage D; p_rhd_death_A = _B = _C = 0 EXACTLY (structural zeros,
#             not tiny values). Competing non-RHD mortality remains in every living
#             state incl. No-RHD. ONE D-mortality scale is estimated (the residual
#             multiplier is FIXED at 1, reported as a documented fixed field).
#             This is an explicit MODELLING ASSUMPTION requiring clinical review —
#             it is not established by Coates 2021 or the WHF staging.
# "all_stage" (LEGACY comparator): severity-ordered A<=B<C<D mortality + bounded
#             residual multiplier (the previous specification). Selected only
#             explicitly; never a silent fallback if D-only fails to fit.
mortality_structure <- "D_only"

# D-only prior CENTER for the annual stage-D RHD-death probability. PROVISIONAL:
# column K of dev/TP Stage D to death_RHDparameters.xlsx = GBD 2023 RHD deaths /
# GBD "HF due to RHD" proxy D stock (all ages, both sexes pooled; Indonesia 11.9%,
# Uganda 5.97%). Column K is the previously requested provisional scenario, NOT a
# validated annual transition probability; its age strata (0-39 / 40+: 40% and 73%
# for 40+) are NOT used because the proxy D denominators are tiny (377 cases) and
# the implied annual probabilities are implausible. No sex difference is imposed
# (evidence pooled). Column J (literature 16.2% severe proxy -> 2.4% / 0.3%) is
# the alternative definition and is reported as a diagnostic, not mixed in.
# Coates's 25% untreated severe-RHD fatality is NOT applied automatically.
p_rhd_death_D_only_prior <- switch(COUNTRY,
                                   Indonesia = 0.119,   # [PROVISIONAL workbook col K]
                                   Uganda    = 0.0597)  # [PROVISIONAL workbook col K]
# D definition used with that prior (kept consistent in diagnostics: model D share
# is compared with the GBD HF-due-to-RHD proxy share, col G of the workbook).
d_share_proxy_gbd_hf <- switch(COUNTRY, Indonesia = 0.0321, Uganda = 0.00863)

# Stage-D share of prevalent RHD at the seed year (PRIOR CENTER; 04b calibrates it).
rhd_d_fraction <- 0.10

# Optional parsimonious AGE variation in A -> No-RHD regression: evidence for
# regression of early/borderline disease is mainly in people under ~20 (Coates 2021
# appendix; GOAL trial ages 5-17). Multiplier applied to p_A_to_no_rhd at ages
# >= remission_age_cutoff. Default 1 = NO age restriction (explicit; sensitivity-
# test with e.g. 0.5 or 0). Both calibration (04b) and projection (05/06) use it.
remission_adult_mult <- 1.0
remission_age_cutoff <- 20L
# PARSIMONIOUS AGE STRUCTURE calibrated in 04b (two parameters, log scale, weak
# priors centred on 1 and remission_adult_mult): a multiplier on the three
# PROGRESSION hazards (A->B, B->C, C->D) and one on A->No-RHD REMISSION at ages
# >= remission_age_cutoff. Justification (diagnosed, not assumed): progression /
# regression evidence comes from child and adolescent cohorts, while GBD's age
# pattern (deaths at ages 5-19 alongside a large adult prevalence with low adult
# case fatality) cannot be reproduced by one all-age set of transitions under
# D-only mortality. FALSE = fixed at remission_adult_mult / 1 (sensitivity).
STRUCT_CALIB_AGE_STRUCTURE <- TRUE
SIGMA_AGE_MULT <- 1.0                  # log-sd prior on the two adult multipliers
AGE_MULT_LO <- 0.05; AGE_MULT_HI <- 3.0

#===============================================================================
# G2. STRUCTURAL (LAYER-2) A/B/C/D CALIBRATION  (04b)
#     Runs the REAL engine (R/abcd_engine.R) inside the objective and calibrates
#     transitions + severity-ratio stage mortality + rhd_d_fraction against GBD
#     RHD prevalence + RHD-SPECIFIC deaths (2000-2019). Mortality carries TIGHT
#     clinical priors, stock/progression LOOSER priors, so excess deaths are
#     resolved by moving the stock, not by suppressing mortality.
#===============================================================================
run_structural_calibration <- switch(COUNTRY, Indonesia = TRUE, Uganda = TRUE)
calibrate_p_D_to_C <- FALSE          # keep D->C fixed at its prior (0) unless TRUE
STRUCT_MATCH_FROM  <- 2010L          # first GBD year matched (skip the seed transient)
# The structural projection is EQUILIBRATED by a zero-start burn-in using the
# first-window-year rates, so the scored stage distribution (and hence the 06
# seed) is the transition-implied EQUILIBRIUM — not the Cannon clinical-cohort
# split, whose slow-draining C/D would otherwise contaminate the fit. Prevalence
# is therefore GENERATED from incidence + natural history and fit to GBD (rather
# than seeded to match), so the regression rate that reconciles the mass balance
# is genuinely data-identified.
STRUCT_BURN_YEARS  <- 120L           # burn-in cycles to reach the stage equilibrium (04b CHECKS
                                     #   convergence against a 2x burn-in; 60 was not converged
                                     #   at adult ages where exit rates are small)
# Count-scale epsilon inside log() for the structural MSLE terms: ONE case
# (log(x + 1)), so a near-zero model cell against a small GBD count is penalised
# strongly but boundedly; cells with GBD count < TARGET_MIN_COUNT are masked.
STRUCT_EPS_COUNT   <- 1

# priors: mortality TIGHT (log-sd), transitions LOOSE (logit-sd) so the objective
# reconciles GBD deaths by moving the STOCK/progression, not the mortality.
SIGMA_MORT   <- 0.35                  # log-sd on stage mortalities -> stay near clinical
SIGMA_TRANS  <- 2.00                  # logit-sd on transitions      -> data-driven stock
SIGMA_RESID  <- 0.25                  # log-sd on the residual mortality multiplier (all_stage only)
# PRIOR STRENGTH (transparent): the data terms are MEANS over ~400 target cells
# while the prior penalties are SUMS over ~20 parameters, so at LAMBDA = 1 a
# one-sigma move of a single parameter costs as much as a 2.7x misfit of EVERY
# target cell — the priors, not the data, were setting the age structure (04b
# exports the prior/data ratio). The MORTALITY prior keeps full strength (it is
# the provisional column-K evidence on the one D-mortality scale, bounded by
# MD_LO/MD_HI); the STRUCTURAL priors (transitions, incidence multipliers, adult
# age multipliers) are down-weighted so the literature RANGES (TRANS_BOUNDS) carry
# the clinical constraint and the data pick the value within them. 0.05 ~ treats
# the target grid as ~20 effective independent observations.
LAMBDA_PRIOR_MORT   <- 1.0
LAMBDA_PRIOR_STRUCT <- 0.05

# --- PREVALENCE-LED JOINT OBJECTIVE (04b) -------------------------------------
#   L = wP * MSLE(prevalence) + wD * MSLE(RHD deaths) + wI * MSLE(incidence)
#       + prior penalty + smoothness penalty + optional stage-target error
# Prevalence (incl. its age pattern) is the PRINCIPAL target; RHD deaths are fit
# jointly. The prevalence:death weight RATIO is the primary configurable setting
# (2.5:1 primary; 2:1 and 3:1 sensitivity). The TOTAL data weight is held fixed
# across ratio comparisons so prior strength does not change accidentally. These
# ratios are operational choices for this implementation, not published weights.
STRUCT_PREV_DEATH_RATIO <- 2.5        # wP / wD  (sensitivity: 2, 3)
STRUCT_DATA_WEIGHT_TOTAL <- 3.5       # wP + wD  (fixed scale; 2.5:1 -> wP=2.5, wD=1)
W_PREV_STRUCT  <- STRUCT_DATA_WEIGHT_TOTAL * STRUCT_PREV_DEATH_RATIO / (1 + STRUCT_PREV_DEATH_RATIO)
W_DEATH_STRUCT <- STRUCT_DATA_WEIGHT_TOTAL / (1 + STRUCT_PREV_DEATH_RATIO)
# Direct incidence weight: 0 in the PRIMARY objective (incidence is estimated and
# REPORTED as a validation output, not fit directly). A weak non-zero value is an
# explicitly labelled sensitivity analysis.
W_INC_STRUCT   <- 0
# Incidence corrections estimated JOINTLY in 04b: one log-multiplier per broad
# age band (CALIB_AGE_BANDS) x sex on the GBD incidence pattern, with a WEAK prior
# centred on the GBD pattern (log-sd below) and a 2nd-difference smoothness
# penalty across bands. Layer-1 alpha is used ONLY as a starting value.
SIGMA_ALPHA_STRUCT  <- 1.0            # weak log-sd prior on the structural incidence multipliers
LAMBDA_SMOOTH_STRUCT <- 0.5           # 2nd-difference smoothness across bands (log scale)
STRUCT_ALPHA_LO <- 0.2; STRUCT_ALPHA_HI <- 5.0   # incidence-multiplier bounds (same as Layer 1)
# Optional inverse-variance weights from GBD UIs in the STRUCTURAL objective.
# Default FALSE: GBD UIs are not independent measurement errors across age/year.
# When TRUE, weights are normalised to mean 1 and capped at STRUCT_IV_WEIGHT_CAP.
USE_IV_WEIGHTS_STRUCT <- FALSE
STRUCT_IV_WEIGHT_CAP  <- 5
# Optional documented AGE weights on the prevalence/death targets (named by band
# label, e.g. c("0-14"=1.5, "15-24"=1.25)); NULL = equal weights. Older ages are
# never excluded — only re-weighted when a documented reason is supplied here.
STRUCT_AGE_WEIGHTS <- NULL
# Target construction basis for Layer 2 (built in 04, consumed by 04b):
#   "rate_x_pop" (DEFAULT): GBD Rate x WPP population of the SAME age-sex-year cell
#                 (the model's own denominator; avoids mixing a WPP denominator with
#                 a GBD count — the two population sources differ by up to 5x at 95+)
#   "number"    : published GBD Number directly.
# Both are saved side by side with their ratio (target reconciliation table).
STRUCT_TARGET_BASIS <- "rate_x_pop"
# Prevalence OBSERVATION convention vs GBD annual prevalence: "mid" = mean of
# start-of-year and end-of-year (pre-ageing) stocks (mid-year point prevalence
# approximation); "start" / "end" alternatives. Deaths/incidence are within-year flows.
STRUCT_PREV_OBS <- "mid"
# Prevalence-only preliminary optimisation to supply starting values (then the
# final fit ALWAYS evaluates prevalence and deaths together).
STRUCT_PREV_FIRST <- TRUE
# Fit-status TOLERANCES (|log(model/GBD)| at the last matched year, total): used to
# LABEL epidemiological fit status (within/outside tolerance), never to certify or
# fail a run. Configurable, not a universal acceptance threshold.
STRUCT_FIT_TOL_PREV  <- 0.15          # ~ +/-16%
STRUCT_FIT_TOL_DEATH <- 0.25          # ~ +/-28%
# EXPLICIT (never silent) age floor for the RHD-DEATH target cells: GBD age groups
# whose lower age is below this are MASKED from the deaths loss (still reported as
# validation). Diagnosed conflict: GBD attributes ~15% of RHD deaths to ages 1-19,
# which a D-only natural history entered through stage A cannot produce (it would
# require acute-carditis / rapidly fatal childhood pathways that are NOT modelled
# — a clinical question, not an automatic new pathway). Including them forces
# faster progression everywhere and inflates adult deaths 1.4-1.6x (see the
# labelled all-ages comparison). 20 = PRIMARY (documented observation-definition
# decision, flagged for clinical review); 0 = all published cells in the loss.
STRUCT_DEATH_TARGET_AGE_MIN <- 20L

# bounds (documented plausible ranges)
MD_LO <- 0.03; MD_HI <- 0.25         # stage-D untreated annual mortality bounds: col-K all-ages
                                     # values (0.119 / 0.0597) lie inside; Coates's 25% is the upper
                                     # plausible untreated severe-RHD fatality. NOT widened for the
                                     # implausible age-stratified col-K values (see block G).
RESID_MULT_LO <- 0.5; RESID_MULT_HI <- 2.0   # bounded residual mult (all_stage only; bound-hit = FAILURE)
# Transition-specific bounds (annual probabilities), informed by the ANNUALISED
# literature ranges in dev/Summary progression rates_Sep26Update.xlsx (Table 2 /
# "Range TPS, annualized"), widened modestly; replaces the former universal
# 0.0001-0.60 range. lo/hi per transition:
TRANS_BOUNDS <- list(
  p_A_to_no_rhd = c(0.005, 0.60),   # lit 0.007-0.588 (children; see remission_adult_mult)
  p_A_to_B      = c(0.001, 0.25),   # lit 0.0015-0.21
  p_B_to_A      = c(0.002, 0.40),   # lit 0.004-0.20
  p_B_to_C      = c(0.001, 0.20),   # lit 0.001-0.125
  p_C_to_B      = c(0.001, 0.15),   # lit 0.012-0.081
  p_C_to_D      = c(0.001, 0.30),   # lit 0.006-0.25
  p_D_to_C      = c(0.000, 0.10))   # lit 0.004-0.071 (fixed at prior unless calibrate_p_D_to_C)

STRUCT_N_STARTS      <- 2L           # reproducible multi-starts (+ prior-center + Layer-1 starts)
STRUCT_NEAROPT_DELTA <- 0.05         # near-optimal set: keep solutions with L <= L_min*(1+delta)
# Forward-project the calibrated end-of-window state to the ANALYSIS start (bridge
# years CAL_YEAR_END+1 .. ANALYSIS_YEAR_START-1) with the calibrated incidence x
# trend (anchored at CAL_YEAR_END), the held background mortality and REFERENCE
# care, so 05 seeds 06 with a dynamically consistent state (no 2019-shares x 2023-
# prevalence x 2026-population mix). TRUE = primary; FALSE = legacy seed.
STRUCT_BRIDGE_TO_ANALYSIS_START <- TRUE

#===============================================================================
# H. INTERVENTION EFFECTS  (03) — relative risk reductions
#===============================================================================
# --- SAP MECHANISM (R/abcd_engine.R; 03/04b/05/06) -----------------------------
# "progression"      (PRIMARY): SAP lowers the HAZARD of the SUPPORTED progression
#                    pathways; its mortality benefit arises through the stage
#                    trajectories (less entry into D). NO direct mortality effect.
# "direct_mortality" (LEGACY comparator, provisional): the former 55% RRR applied
#                    to every stage's RHD-death probability. Coates 2021's 55% is an
#                    ARF-recurrence risk ratio (Padmavati 1973; RR 0.45) that Coates
#                    applied to the non-HF -> HF progression — NOT evidence for a
#                    uniform direct mortality reduction. Kept for attribution only.
# A calibration fitted under one mechanism is never consumed under another: 04b
# applies the selected mechanism + REFERENCE care in the objective (05 checks).
sap_mechanism <- "progression"
sap_rrr_rhd_death <- 0.55   # [LEGACY/PROVISIONAL] used ONLY when sap_mechanism == "direct_mortality"

# Pathway-specific SAP effects, HAZARD-RATIO scale (sap_effect_scale). Applied as
# h_ij(c) = h_ij(0) x [1 - c (1 - HR_ij)], c = effective protection.
sap_effect_scale <- "hazard_ratio"
# A->B and B->C: GOAL trial (Uganda, children 5-17 with latent RHD, 24-month
# follow-up; control vs prophylaxis arm, near-perfect adherence) as curated in
# dev/Summary progression rates_Sep26Update.xlsx, Table 1 columns I (control) and J
# (trial arm): A->B 0.075 vs 0.003; B->C 0.119 vs 0.025 (2-year risks). Converted to
# a hazard ratio over the common follow-up: HR = log(1-p_sap)/log(1-p_ctrl).
# PROVISIONAL stage mapping: GOAL "borderline" ~ WHF stage A, "mild definite" ~ B.
# [TODO: confirm against Beaton et al. 2022 supplement; not in library/references.bib]
goal_ctrl_A_to_B <- 0.075; goal_sap_A_to_B <- 0.003
goal_ctrl_B_to_C <- 0.119; goal_sap_B_to_C <- 0.025
hr_from_arm_risks <- function(p_sap, p_ctrl) log1p(-p_sap) / log1p(-p_ctrl)
sap_hr_A_to_B <- hr_from_arm_risks(goal_sap_A_to_B, goal_ctrl_A_to_B)   # ~0.039
sap_hr_B_to_C <- hr_from_arm_risks(goal_sap_B_to_C, goal_ctrl_B_to_C)   # ~0.20
# C->D: NO pathway-specific trial evidence. Coates 2021 (appendix Table 6, 2b)
# ASSUMED the ARF-recurrence RR 0.45 (55% reduction) for "RHD without HF -> RHD
# with HF", which maps most closely to C->D here. PROVISIONAL assumption, kept with
# its actual source; set to 1 for the no-C->D-effect sensitivity.
sap_hr_C_to_D <- 0.45   # [Coates2021 assumption, provisional]

# Annual population-average ADHERENCE modifier (default 1 = no additional
# attenuation of the chosen effect; it does NOT mean 100% treatment coverage).
# Effective protection c = effective treatment coverage x adherence. It is a
# RELATIVE scenario attenuation on top of the cascade's "optimal treatment"
# coverage (which is interpreted as effective, adherent treatment), so 1 preserves
# the unattenuated effect. Country-specific values / trajectories may be set per
# arm (baseline -> target over the ramp) when evidence exists; none is invented here.
sap_adherence_ref_baseline <- 1.0; sap_adherence_ref_target <- 1.0
sap_adherence_up_baseline  <- 1.0; sap_adherence_up_target  <- 1.0
# Baseline-care interpretation: calibrated transitions/mortality are UNTREATED
# natural-history parameters; REFERENCE care (treatment_base x adherence, surgery
# cov_surgery_ref) is applied inside the engine in BOTH calibration and projection
# (applied once; no renormalisation).
baseline_care_interpretation <- "untreated_parameters_reference_care_applied_in_engine"

# Surgery effects (applied to C and D only; parameterised separately so either
# can be set to 0). Scaled in the engine by effective surgery reach = fraction
# requiring surgery x surgery coverage.
eff_surgery_C_to_D        <- 0.85  # [LIT]/[CALIBRATE] surgery RRR on C -> D progression
eff_surgery_D_to_rhd_death <- 0.85 # [LIT]/[CALIBRATE] surgery RRR on D -> RHD death

#===============================================================================
# I. CARE CASCADE + SURGERY SERVICE  (03/05)
#    Three cascade metrics (screening / diagnosis / optimal treatment) as
#    CUMULATIVE population coverages. Reference holds baselines flat; scale-up
#    ramps each linearly from the ramp_start baseline to the ramp_end target
#    (held at target through the end of the analysis horizon).
#    Effective diagnosis/treatment are capped by the earlier cascade stages (05).
#    [TODO: confirm for Uganda] the BASELINE cascade coverages and the surgery
#    requirement/coverage fractions below are health-system-specific and were set
#    for Indonesia; Uganda's baseline coverage and cardiac-surgery capacity likely
#    differ. Targets (0.80/0.65) are policy-scenario assumptions held common.
#===============================================================================
# screen_base    <- 0.05   # [PAPER] 2026 baseline screened
# diagnosis_base <- 0.05   # [PAPER] 2026 baseline diagnosed
# treatment_base <- 0.04   # [PAPER] 2026 baseline on optimal treatment
# 
# screen_target    <- 0.80 # [PAPER] 2050 screening target
# diagnosis_target <- 0.80 # [PAPER] 2050 diagnosis target
# treatment_target <- 0.65 # [PAPER] 2050 optimal-treatment target
# 
# # SURGERY as a clinical service (NOT a health state). A fixed fraction of the
# # prevalent C / D stock requires surgery each cycle; surgery coverage is the
# # share actually delivered. Held equal in both arms (tertiary care fixed) so
# # surgery is a background cost, not a driver of incremental program cost.
# frac_C_requiring_surgery <- 0.03  # [CALIBRATE] share of stage C needing surgery /yr
# frac_D_requiring_surgery <- 0.20  # [CALIBRATE] share of stage D needing surgery /yr
# cov_surgery_ref      <- 0.05      # [CALIBRATE] surgery coverage, reference arm
# cov_surgery_scale_up <- 0.05      # [CALIBRATE] surgery coverage, scale-up arm
# 
# # Optional school-age-only screening mask (RETIRED from the primary spec:
# # screening cost now applies to the TOTAL population screened). Kept OFF by
# # default for a sensitivity analysis only.
# screen_age_restrict <- FALSE      # [CALIBRATE] TRUE => restrict screening to [lo,hi]
# screen_age_lo       <- 5L         # [LIT] school-based echo screening window (low)
# screen_age_hi       <- 15L        # [LIT] school-based echo screening window (high)

screen_base <- switch(COUNTRY,
                      Indonesia = 0.05,
                      Uganda    = 0.01)
diagnosis_base <- switch(COUNTRY,
                         Indonesia = 0.05,
                         Uganda    = 0.009)
treatment_base <- switch(COUNTRY,
                         Indonesia = 0.04,
                         Uganda    = 0.0075)

screen_target <- switch(COUNTRY,
                        Indonesia = 0.80,
                        Uganda    = 0.80)
diagnosis_target <- switch(COUNTRY,
                           Indonesia = 0.80,
                           Uganda    = 0.80)
treatment_target <- switch(COUNTRY,
                           Indonesia = 0.65,
                           Uganda    = 0.65)

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

screen_age_restrict <- switch(COUNTRY,
                              Indonesia = FALSE,
                              Uganda    = FALSE)
screen_age_lo <- switch(COUNTRY,
                        Indonesia = 5L,
                        Uganda    = 5L)
screen_age_hi <- switch(COUNTRY,
                        Indonesia = 15L,
                        Uganda    = 15L)

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
gdp_pc_base <- switch(COUNTRY,
                      Indonesia = 4150,
                      Uganda    = 1340)
gdp_pc_base_year <- switch(COUNTRY,
                           Indonesia = 2019L,
                           Uganda    = 2019L)
gdp_growth <- switch(COUNTRY,
                     Indonesia = 0.03,
                     Uganda    = 0.032)
vsl_mult <- switch(COUNTRY,
                   Indonesia = 30,
                   Uganda    = 30)
dalys_per_death <- switch(COUNTRY,
                          Indonesia = 30,
                          Uganda    = 30)
# (7) Economic discount reference year — costs & benefits are discounted TO this
#     year. Default = the first analysis year (kept separate from (5) the GBD rate
#     reference year and (6) gdp_pc_base_year, which are their own anchors).
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
#===============================================================================
setwd(wd_code)   # scripts use absolute wd_* globals for I/O and here() for R/packages.R

source(paste0(wd_code, "01_prepare_inputs.R"))        # GBD inputs
source(paste0(wd_code, "02_build_demography.R"))      # WPP2024 population backbone
source(paste0(wd_code, "03_build_disease_model.R"))   # disease-model INPUTS (builder)
source(paste0(wd_code, "04_calibration_random_tp.R")) # Layer 1: calibrate IR + CF anchor (2000-2019)
source(paste0(wd_code, "04b_calibrate_structural.R")) # Layer 2: structural A/B/C/D stage calibration
source(paste0(wd_code, "05_build_baseline.R"))        # assemble initial state (uses 04b stage params)
source(paste0(wd_code, "06_run_prevention_model.R"))  # run ref + SAP (parallel by location)
source(paste0(wd_code, "07_make_outputs.R"))          # standard long table
source(paste0(wd_code, "08_economic_evaluation.R"))   # benefit-cost & cost-effectiveness

message("\n================ PIPELINE COMPLETE (01 -> 08) ================")

