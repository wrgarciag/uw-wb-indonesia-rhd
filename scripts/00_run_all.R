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
#   04 calibration     -> data/adjusted_searo_part*.rds            (calibrated IR/CF, 2000-2019)
#   05 initial state   -> data/baseline_state.rds                  (assembles 02+03+04)
#   06 run model       -> output/out_model/<location>.rds          (parallel by location)
#   07 make outputs    -> output/tables/rhd_model_long.csv/.rds    (standard long table)
#   08 economics       -> output/tables/rhd_economic_summary.csv, rhd_budget_impact.csv
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
ramp_start <- 2025L              # first year coverage begins rising
ramp_end   <- 2030L              # year target coverage is reached (then held)

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
# G. DISEASE / CLINICAL PARAMETERS  (03) — annual transition probabilities
#    RHD-specific progression/mortality that the interventions act on.
#===============================================================================
p_mild_to_severe      <- 0.010   # [CALIBRATE] asymptomatic -> heart failure /yr
p_severe_death        <- 0.09    # [LIT] untreated severe RHD (HF) mortality /yr
p_surg_op_mortality   <- 0.03    # [PAPER] operative mortality
p_post_death_rhd      <- 0.020   # [LIT] residual RHD mortality after valve surgery
frac_severe_surg_elig <- 0.50    # [CALIBRATE] share of severe RHD surgery-eligible

#===============================================================================
# H. INTERVENTION EFFECT SIZES + SEED SPLIT  (03) — relative risk reductions
#===============================================================================
eff_sap_asymp <- 0.55            # [PAPER] SAP in asymptomatic RHD (55%)
eff_hf_mgmt   <- 0.60            # [PAPER] HF management (60%)
eff_surgery   <- 0.85            # [PAPER] surgery (85%)

# prevalent RHD severity split at the seed year (must sum to 1) -- [LIT]/[CALIBRATE]
seed_split <- c(mild = 0.96, severe = 0.03, post = 0.01)

#===============================================================================
# I. COVERAGE TARGETS  (03) — only SAP differs between reference and scale-up
#===============================================================================
cov_sap_ref    <- 0.05           # [PAPER] baseline SAP-in-asymptomatic coverage (reference)
cov_sap_target <- 0.40           # [PAPER] SAP coverage under scale-up (reached by ramp_end)
cov_hf         <- 0.08           # [PAPER] HF management coverage (held in both arms)
cov_surg       <- 0.05           # [PAPER] surgery coverage (held in both arms)
screen_age_lo  <- 5L             # [LIT] school-based echo screening age window (low)
screen_age_hi  <- 15L            # [LIT] school-based echo screening age window (high)

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

# unit costs (base-year US$) -- [LIT]/[CALIBRATE]
cost_screen_per_person <- 12     # echo screening cost per person
cost_sap_per_year      <- 45     # annual SAP cost per person
cost_hf_per_year       <- 120    # annual HF management cost per person
cost_surgery           <- 9000   # valve surgery + first-year post-op cost per person

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
