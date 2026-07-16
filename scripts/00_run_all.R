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
#   02 demography      -> data/pop_observed_1990_2024.rds, pop_projection_2025_2100.rds
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
#   Secondary antibiotic prophylaxis (SAP) reduces RHD-specific MORTALITY (55%
#   RRR) among people on optimal treatment; it does not reduce incidence here.
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
# A. PATHS
#===============================================================================
wd      <- "C:/Users/wrgar/OneDrive - UW/02Work/WorldBank-Indonesia/uw-wb-indonesia-rhd/"
wd_code <- paste0(wd, "scripts/")
wd_raw  <- paste0(wd, "data-raw/")   # raw GBD/WPP inputs (not on GitHub)
wd_data <- paste0(wd, "data/")       # processed inputs + calibrated TPs + baseline state
wd_outp <- paste0(wd, "output/")     # model outputs + tables

#===============================================================================
# B. LOCATIONS
#===============================================================================
LOCATION  <- "Indonesia"        # single location used by 03/04 (rate/TP build)
LOCATIONS <- c("Indonesia")     # location list run by 05/06 (parallelised in 06)

#===============================================================================
# C. HORIZON / YEAR WINDOWS
#===============================================================================
# Demography (02): observed 1990-2024, medium-variant projection 2025-2100.
# These close the former 2024-2025 gap; keep them contiguous (obs ends the year
# before proj begins) or 02's join-continuity check will fail.
OBS_YEARS      <- 1990:2024
PROJ_YEARS     <- 2025:2100
RATE_BASE_YEAR <- 2023L          # latest GBD rate year (held forward by 03)

# Calibration in-sample window (04)
CAL_YEAR_START <- 2000L
CAL_YEAR_END   <- 2019L

#===============================================================================
# D. INTERVENTION SCALE-UP RAMP WINDOW  (03/05)
#===============================================================================
ramp_start <- 2026L              # first year cascade coverage begins rising (baseline year)
ramp_end   <- 2050L              # year target cascade coverage is reached (then held to 2100)

#===============================================================================
# E. INCIDENCE SECULAR TREND  (03/05)
#===============================================================================
incidence_trend <- 0.985         # [CALIBRATE] ~15%/decade RHD incidence decline

#===============================================================================
# F. CALIBRATION SETTINGS  (04) — pure random search
#===============================================================================
run_calibration_par <- TRUE      # parallelise calibration across combos
SEARCH_HALFWIDTH    <- 0.50      # IR/CF multipliers sampled in [1-hw, 1+hw]
GRANULARITY         <- "age_group"  # "combo" | "age_group"
N_ITER              <- 400       # i.i.d. candidates per combo (candidate 0 = baseline)
CONVERGE_TOL        <- 1e-4      # early stop when best weighted error < tol
SEED                <- 42        # master RNG seed (reproducible)
W_DEATHS            <- 2         # objective weight on (all-cause) deaths
W_PREV              <- 1         # objective weight on RHD prevalence
# (calibration age range is the full 0-95+; see AGE_LO/AGE_HI defaults in 04)

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
p_A_to_no_rhd <- 0.005   # [CALIBRATE] A regression to No RHD /yr
p_A_to_B      <- 0.020   # [CALIBRATE] A -> B progression /yr
p_B_to_A      <- 0.010   # [CALIBRATE] B -> A regression /yr
p_B_to_C      <- 0.030   # [CALIBRATE] B -> C progression /yr
p_C_to_B      <- 0.005   # [CALIBRATE] C -> B regression /yr
p_C_to_D      <- 0.060   # [CALIBRATE] C -> D progression /yr (surgery lowers this)
p_D_to_C      <- 0.000   # [CALIBRATE] D -> C regression /yr (default 0)

# annual UNTREATED RHD-specific death probability by stage -- [CALIBRATE]
p_rhd_death_A <- 0.0005  # [CALIBRATE]
p_rhd_death_B <- 0.0020  # [CALIBRATE]
p_rhd_death_C <- 0.0200  # [CALIBRATE]
p_rhd_death_D <- 0.0800  # [CALIBRATE]

# Stage-D share of prevalent RHD at the seed year. GBD gives only TOTAL RHD
# prevalence; A/B/C come from Cannon et al (56.5/27.2/16.2, normalised to the
# non-D share in 03) and D is this Indonesia RHD-with-heart-failure fraction.
rhd_d_fraction <- 0.10   # [CALIBRATE] complicated-RHD / RHD-with-HF fraction

#===============================================================================
# H. INTERVENTION EFFECTS  (03) — relative risk reductions
#===============================================================================
# Secondary prevention (SAP / optimal treatment): 55% RRR on the transition
# from EACH living RHD stage to RHD-specific death. It does NOT reduce incidence.
sap_rrr_rhd_death <- 0.55   # [PAPER] SAP RRR on RHD-specific mortality (55%)

# Surgery effects (applied to C and D only; parameterised separately so either
# can be set to 0). Scaled in the engine by effective surgery reach = fraction
# requiring surgery x surgery coverage.
eff_surgery_C_to_D        <- 0.85  # [LIT]/[CALIBRATE] surgery RRR on C -> D progression
eff_surgery_D_to_rhd_death <- 0.85 # [LIT]/[CALIBRATE] surgery RRR on D -> RHD death

#===============================================================================
# I. CARE CASCADE + SURGERY SERVICE  (03/05)
#    Three cascade metrics (screening / diagnosis / optimal treatment) as
#    CUMULATIVE population coverages. Reference holds baselines flat; scale-up
#    ramps each linearly from the 2026 baseline to the 2050 target (held to 2100).
#    Effective diagnosis/treatment are capped by the earlier cascade stages (05).
#===============================================================================
screen_base    <- 0.05   # [PAPER] 2026 baseline screened
diagnosis_base <- 0.05   # [PAPER] 2026 baseline diagnosed
treatment_base <- 0.04   # [PAPER] 2026 baseline on optimal treatment

screen_target    <- 0.80 # [PAPER] 2050 screening target
diagnosis_target <- 0.80 # [PAPER] 2050 diagnosis target
treatment_target <- 0.65 # [PAPER] 2050 optimal-treatment target

# SURGERY as a clinical service (NOT a health state). A fixed fraction of the
# prevalent C / D stock requires surgery each cycle; surgery coverage is the
# share actually delivered. Held equal in both arms (tertiary care fixed) so
# surgery is a background cost, not a driver of incremental program cost.
frac_C_requiring_surgery <- 0.03  # [CALIBRATE] share of stage C needing surgery /yr
frac_D_requiring_surgery <- 0.20  # [CALIBRATE] share of stage D needing surgery /yr
cov_surgery_ref      <- 0.05      # [CALIBRATE] surgery coverage, reference arm
cov_surgery_scale_up <- 0.05      # [CALIBRATE] surgery coverage, scale-up arm

# Optional school-age-only screening mask (RETIRED from the primary spec:
# screening cost now applies to the TOTAL population screened). Kept OFF by
# default for a sensitivity analysis only.
screen_age_restrict <- FALSE      # [CALIBRATE] TRUE => restrict screening to [lo,hi]
screen_age_lo       <- 5L         # [LIT] school-based echo screening window (low)
screen_age_hi       <- 15L        # [LIT] school-based echo screening window (high)

#===============================================================================
# J. ECONOMIC PARAMETERS  (08 ONLY) — the only monetary values in the pipeline
#===============================================================================
disc_rate        <- 0.03         # [PAPER] annual discount rate (costs & benefits)
gdp_pc_base      <- 4150         # [LIT] Indonesia GDP per capita, base-year US$
gdp_pc_base_year <- 2019L        # [LIT] year gdp_pc_base refers to
gdp_growth       <- 0.03         # [CALIBRATE] real per-capita GDP growth
vsl_mult         <- 30           # [PAPER] value of a statistical life = 30 x GDP per capita
dalys_per_death  <- 30           # [CALIBRATE] undiscounted DALYs per RHD death
discount_base_year <- 2025L      # discount to this year (default = first model year)

# unit costs (base-year US$) -- [PAPER]/[LIT]
cost_screen_per_person <- 1.10   # [PAPER] echo screening cost per person screened
cost_sap_per_year      <- 110    # [PAPER] annual optimal secondary-prevention cost per person
cost_surgery           <- 9000   # [LIT] valve surgery + first-year post-op cost per surgery

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
source(paste0(wd_code, "04_calibration_random_tp.R")) # calibrate IR/CF (2000-2019)
source(paste0(wd_code, "05_build_baseline.R"))        # assemble initial state
source(paste0(wd_code, "06_run_prevention_model.R"))  # run ref + SAP (parallel by location)
source(paste0(wd_code, "07_make_outputs.R"))          # standard long table
source(paste0(wd_code, "08_economic_evaluation.R"))   # benefit-cost & cost-effectiveness

message("\n================ PIPELINE COMPLETE (01 -> 08) ================")
