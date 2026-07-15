# ==============================================================================
# RHD secondary-prevention investment case: DISEASE-MODEL INPUT BUILDER
# scripts/03_build_disease_model.R
#
# Structure after: Coates et al., Lancet Glob Health 2021 (PMC9087136).
# Focus: scale-up of SECONDARY PREVENTION =
#   (a) echocardiographic screening to detect asymptomatic (mild) RHD, and
#   (b) secondary antibiotic prophylaxis (SAP) for screen-detected mild RHD.
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT (repositioned in this refactor)
# ------------------------------------------------------------------------------
# This script ONLY BUILDS AND VALIDATES the disease-model INPUTS. It does NOT run
# the baseline projection or the intervention scenarios — that is done by the
# initial-state assembler (05_build_baseline.R) and the model runner
# (06_run_prevention_model.R). The former hard-coded Australia scalars, the
# `run_scenario()` engine, the `ref`/`sap` runs, and ALL economic / cost values
# have been removed from here (the engine now lives in 06; economics in 08).
#
# It assembles data-fed age x sex ( x year) arrays plus the tagged clinical /
# intervention / coverage parameter blocks, validates them, and PERSISTS a single
# self-contained input bundle to disk so 05/06 can consume it without recomputation.
#
#   INPUTS
#     data-raw/temp_baseline_rates_gbd.rds   (from 01_prepare_inputs.R)
#         GBD 2023: RHD + All-causes; Deaths / Prevalence / Incidence; Number+Rate;
#         1990-2023; 22 GBD age groups. Used (metric = Rate, per capita):
#           RHD Incidence  -> incident inflow of new asymptomatic RHD
#           RHD Prevalence -> prevalent RHD seed at the first model year
#           RHD Deaths     -> RHD cause-specific mortality anchor
#           All-cause Deaths - RHD Deaths -> background (other-cause) mortality
#     data/pop_projection_2025_2100.rds      (from 02_build_demography.R)
#         single-year population by age x sex x year (persons), Indonesia. Read here
#         ONLY to establish the projection horizon (year grid) and a base-year
#         order-of-magnitude anchor. 05 re-reads the full population itself.
#
#   OUTPUT (written to wd_data = data/):
#     disease_model_inputs.rds  — a named list (see section 8) holding:
#       rates_by_year : ir/mort_rhd/mort_all/oth_mort as [age x sex x year] arrays,
#                       plus prev_seed [age x sex] (single-year prevalent seed);
#       clinical      : progression / case-fatality clinical parameters (tagged);
#       effects       : intervention effect sizes (tagged);
#       seed_split    : prevalent-pool severity split mild/severe/post (tagged);
#       coverage      : coverage-ramp STRUCTURE (baselines, targets, ramp window,
#                       screening age window) — realised into trajectories by 05;
#       meta          : AGES, SEXES, years, LOCATION, RATE_BASE_YEAR, incidence_trend.
#
# GBD gives TOTAL RHD prevalence (not split by severity) and rates only to 2023:
#   * the prevalent pool is split into mild/severe/post by a tagged [LIT] vector
#     (asymptomatic RHD dominates prevalence);
#   * for the projection horizon the GBD age-sex RATE pattern is held at its 2023
#     level, with an explicit incidence secular trend ([CALIBRATE]) baked into the
#     incidence array (anchored at the first projection year). Both are flagged.
#
# PARAMETER TAGS: [PAPER] article main text; [LIT] literature value (source in
#   comment); [CALIBRATE] tune to setting/appendix. Edit only the PARAMETER /
#   COVERAGE blocks; the demographic & epidemiological inputs come from the tables.
#   NO monetary values live here — unit costs / VSL / DALY weights are in 08 only.
# ==============================================================================

library(data.table)

# ------------------------------------------------------------------------------
# 0. PATHS + DIMENSIONS
# ------------------------------------------------------------------------------
if (!exists("wd_raw"))  wd_raw  <- paste0(here::here("data-raw"), "/")
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

# getp(): use the 00_run_all.R global when set, else the documented standalone
# default below. Every user-controllable parameter is funnelled through this so
# 00_run_all.R can act as the single control panel while 03 still runs standalone.
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

LOCATION <- getp("LOCATION", "Indonesia")   # engine runs per location
AGES     <- 0:100                # single-year ages (0..100, 100 = 100+ open group)
SEXES    <- c("Female", "Male")
n_age    <- length(AGES)

RATE_BASE_YEAR <- getp("RATE_BASE_YEAR", 2023L)  # latest GBD year
POP_PROJ_FILE  <- paste0(wd_data, "pop_projection_2025_2100.rds")

# ------------------------------------------------------------------------------
# 1. RAMP WINDOW (structure only; no discounting/economics here)
# ------------------------------------------------------------------------------
# Coverage ramp window (honour globals from 00_run_all.R when present).
ramp_start <- as.integer(getp("ramp_start", 2025L))   # first year coverage begins rising
ramp_end   <- as.integer(getp("ramp_end",   2030L))   # year target coverage is reached

# ------------------------------------------------------------------------------
# 2. DISEASE / CLINICAL PARAMETERS  (annual transition probabilities; tagged)
#    RHD-specific mortality & progression stay parameters because the SAP / HF /
#    surgery levers act on them; competing (other-cause) mortality is data-fed.
#    All values are user-controllable via 00_run_all.R (getp fallbacks shown).
# ------------------------------------------------------------------------------
incidence_trend <- getp("incidence_trend", 0.985)
                                 # [CALIBRATE] ~15%/decade incidence decline (Table 1 footnote)

clinical <- list(
  p_mild_to_severe      = getp("p_mild_to_severe",     0.010),  # [CALIBRATE] asymptomatic -> HF /yr
                                 #   (tuned so uncalibrated baseline RHD deaths sit
                                 #    near the GBD 2023 rate-implied level; 04 refines)
  p_severe_death        = getp("p_severe_death",       0.09),   # [LIT] REMEDY ~17%/2yr -> ~9%/yr untreated
  p_surg_op_mortality   = getp("p_surg_op_mortality",  0.03),   # [PAPER] 3% operative mortality
  p_post_death_rhd      = getp("p_post_death_rhd",     0.020),  # [LIT] residual RHD mortality post-surgery
  frac_severe_surg_elig = getp("frac_severe_surg_elig", 0.50)   # [CALIBRATE] share surgery-eligible
)

effects <- list(
  eff_sap_asymp        = getp("eff_sap_asymp", 0.55),  # [PAPER] Table 1 #2b: SAP in asymptomatic RHD, 55% (7-78)
  eff_hf_mgmt          = getp("eff_hf_mgmt",   0.60),  # [PAPER] Table 1 #3: HF management, 60% (30-80)
  eff_surgery          = getp("eff_surgery",   0.85)   # [PAPER] Table 1 #4: surgery, 85% (70-92)
)

# prevalent RHD severity split (GBD gives only TOTAL prevalence) -- [LIT]/[CALIBRATE]
# Echo-detected RHD is overwhelmingly subclinical/asymptomatic; symptomatic HF
# (severe) and post-surgical stocks are small shares.
seed_split <- getp("seed_split", c(mild = 0.96, severe = 0.03, post = 0.01))  # must sum to 1
stopifnot(abs(sum(seed_split) - 1) < 1e-9)

# ------------------------------------------------------------------------------
# 3. COVERAGE-RAMP STRUCTURE  (baselines/targets; realised into trajectories by 05)
#    Only SAP coverage differs between reference and scale-up; HF and surgery are
#    held at baseline in both arms to isolate the SECONDARY-prevention effect.
#    Screening age window is STRUCTURAL (who gets screened) — not a monetary value.
# ------------------------------------------------------------------------------
coverage <- list(
  sap_ref_baseline = getp("cov_sap_ref",    0.05), # [PAPER] Table 1 #2b: 5.0% flat (reference)
  sap_ref_target   = getp("cov_sap_ref",    0.05),
  sap_up_baseline  = getp("cov_sap_ref",    0.05), # [PAPER] Table 1 #2b: 5.0% baseline ->
  sap_up_target    = getp("cov_sap_target", 0.40), #                       40% under scale-up
  hf_baseline   = getp("cov_hf",   0.08), hf_target   = getp("cov_hf",   0.08),  # [PAPER] #3, held both arms
  surg_baseline = getp("cov_surg", 0.05), surg_target = getp("cov_surg", 0.05),  # [PAPER] #4, held both arms
  ramp_start = ramp_start,
  ramp_end   = ramp_end,
  screen_age_lo = as.integer(getp("screen_age_lo",  5L)),  # [LIT] school-based echo screening window
  screen_age_hi = as.integer(getp("screen_age_hi", 15L))
)
stopifnot(coverage$ramp_end > coverage$ramp_start)

# ==============================================================================
# 4. BUILD DATA-FED RATE INPUTS  (single age x sex arrays for LOCATION)
# ==============================================================================

## 4a. GBD 2023 base-year rates (per capita) at single-year age -----------------
gbd <- as.data.table(readRDS(paste0(wd_raw, "temp_baseline_rates_gbd.rds")))
gbd <- gbd[location_name == LOCATION & metric_name == "Rate" & year == RATE_BASE_YEAR,
           .(sex = sex_name, age_group = age_name,
             cause = cause_name, measure = measure_name, rate = val / 1e5)]

# single-year age -> GBD age-group label (matches 01's 22 groups incl. <5 split).
# This SAME mapping is mirrored in 04_calibration_random_tp.R.
age_to_gbd_group <- function(a) {
  fcase(
    a < 1,  "<1 year",
    a < 2,  "12-23 months",
    a < 5,  "2-4 years",
    a < 10, "5-9 years",
    a < 15, "10-14 years", a < 20, "15-19 years", a < 25, "20-24 years",
    a < 30, "25-29 years", a < 35, "30-34 years", a < 40, "35-39 years",
    a < 45, "40-44 years", a < 50, "45-49 years", a < 55, "50-54 years",
    a < 60, "55-59 years", a < 65, "60-64 years", a < 70, "65-69 years",
    a < 75, "70-74 years", a < 80, "75-79 years", a < 85, "80-84 years",
    a < 90, "85-89 years", a < 95, "90-94 years",
    default = "95+ years"
  )
}
age_map <- data.table(age = AGES, age_group = age_to_gbd_group(AGES))

# helper: pull one (cause, measure) rate -> [age x sex] matrix, missing -> 0
rate_matrix <- function(cause_name, measure_name) {
  d <- gbd[cause == cause_name & measure == measure_name,
           .(sex, age_group, rate)]
  d <- age_map[d, on = "age_group", allow.cartesian = TRUE]   # expand to single age
  m <- matrix(0, n_age, length(SEXES), dimnames = list(AGES, SEXES))
  for (s in SEXES) {
    ds <- d[sex == s]
    if (nrow(ds)) m[as.character(ds$age), s] <- ds$rate
  }
  m
}

ir_rhd0   <- rate_matrix("Rheumatic heart disease", "Incidence")   # incidence rate
prev_rhd0 <- rate_matrix("Rheumatic heart disease", "Prevalence")  # prevalent fraction
mort_rhd0 <- rate_matrix("Rheumatic heart disease", "Deaths")      # RHD death rate
mort_all0 <- rate_matrix("All causes",              "Deaths")      # all-cause death rate
oth_mort0 <- pmax(mort_all0 - mort_rhd0, 0)   # background (non-RHD) competing mortality

## 4b. establish the projection horizon (year grid) from 02's projection table --
if (!file.exists(POP_PROJ_FILE))
  stop("Missing population projection input:\n  ", POP_PROJ_FILE,
       "\n  Run 02_build_demography.R first.", call. = FALSE)
pop_proj <- as.data.table(readRDS(POP_PROJ_FILE))
pop_proj <- pop_proj[location == LOCATION & age %in% AGES & sex %in% SEXES]
years    <- sort(unique(pop_proj$year))     # 2025..2100  (HORIZON driven by input)
n_years  <- length(years)

## 4c. broadcast base-year [age x sex] matrices to [age x sex x year] arrays -----
##     incidence carries the secular trend (anchored at the first projection year);
##     the mortality envelope is held at the base-year pattern across the horizon.
broadcast_years <- function(mat, mult = rep(1, n_years)) {
  arr <- array(0, dim = c(n_age, length(SEXES), n_years),
               dimnames = list(AGES, SEXES, years))
  for (iy in seq_len(n_years)) arr[, , iy] <- mat * mult[iy]
  arr
}
trend_mult <- incidence_trend^(seq_len(n_years) - 1L)   # 1.0 at first projection year

ir_rhd_yr   <- broadcast_years(ir_rhd0, trend_mult)   # incidence, trend baked in
mort_rhd_yr <- broadcast_years(mort_rhd0)             # RHD death rate (held)
mort_all_yr <- broadcast_years(mort_all0)             # all-cause death rate (held)
oth_mort_yr <- broadcast_years(oth_mort0)             # background mortality (held)

# base-year population [age x sex] for the order-of-magnitude anchor below
py1  <- pop_proj[year == years[1]]
pop1 <- matrix(0, n_age, length(SEXES), dimnames = list(AGES, SEXES))
for (s in SEXES) {
  ps <- py1[sex == s]
  pop1[as.character(ps$age), s] <- ps$Nx
}

# ------------------------------------------------------------------------------
# 5. INPUT VALIDATION  (fail loudly BEFORE persisting)
# ------------------------------------------------------------------------------
chk_rate <- function(x, nm, prob = TRUE) {
  if (any(is.na(x)))  stop("Input '", nm, "' contains NA.", call. = FALSE)
  if (any(x < 0))     stop("Input '", nm, "' contains negative values.", call. = FALSE)
  if (prob && any(x > 1))
    stop("Input '", nm, "' has values > 1 (should be a per-capita probability/rate).",
         call. = FALSE)
}
chk_rate(ir_rhd_yr,   "RHD incidence rate")
chk_rate(prev_rhd0,   "RHD prevalence rate")
chk_rate(mort_rhd_yr, "RHD death rate")
chk_rate(mort_all_yr, "all-cause death rate")
chk_rate(oth_mort_yr, "background (other-cause) mortality")

# arrays complete & correctly shaped
for (nm in c("ir_rhd_yr", "mort_rhd_yr", "mort_all_yr", "oth_mort_yr")) {
  a <- get(nm)
  if (!identical(dim(a), c(n_age, length(SEXES), n_years)))
    stop("Array '", nm, "' has wrong dimensions.", call. = FALSE)
}
if (any(is.na(pop1)) || any(pop1 < 0) || sum(pop1) <= 0)
  stop("Base-year population slice is empty/invalid.", call. = FALSE)
if (!all(SEXES %in% pop_proj$sex))
  stop("Population projection is missing a sex.", call. = FALSE)

# GBD rates must actually carry RHD signal (not all zero)
if (sum(prev_rhd0) == 0 || sum(ir_rhd0) == 0)
  stop("GBD RHD prevalence/incidence rates are all zero -- check inputs.", call. = FALSE)

# competing-risk headroom at the base year: IR + background <= 1, CF proxy + bg <= 1
if (any(ir_rhd0 + oth_mort0 > 1 + 1e-9))
  stop("IR + background mortality exceeds 1 at some age-sex cell.", call. = FALSE)

# parameter-block sanity: probabilities/effects/coverage all in [0,1]
prob_params <- c(unlist(clinical), unlist(effects),
                 coverage$sap_ref_baseline, coverage$sap_ref_target,
                 coverage$sap_up_baseline,  coverage$sap_up_target,
                 coverage$hf_baseline, coverage$hf_target,
                 coverage$surg_baseline, coverage$surg_target)
if (any(prob_params < 0 | prob_params > 1))
  stop("A clinical/effect/coverage parameter is outside [0,1].", call. = FALSE)
if (coverage$screen_age_lo < min(AGES) || coverage$screen_age_hi > max(AGES) ||
    coverage$screen_age_lo > coverage$screen_age_hi)
  stop("Screening age window is invalid.", call. = FALSE)

# order-of-magnitude anchors vs GBD (base-year rate x base-year population)
rhd_prev_count  <- sum(prev_rhd0 * pop1)
rhd_death_count <- sum(mort_rhd0 * pop1)
all_death_count <- sum(mort_all0 * pop1)
message(sprintf(
  "Data-fed inputs OK | horizon %d-%d | ages %d-%d\n  base-year(%d): RHD prevalence ~ %s | RHD deaths ~ %s | all-cause deaths ~ %s",
  as.integer(min(years)), as.integer(max(years)),
  as.integer(min(AGES)), as.integer(max(AGES)), as.integer(years[1]),
  formatC(round(rhd_prev_count),  format = "d", big.mark = ","),
  formatC(round(rhd_death_count), format = "d", big.mark = ","),
  formatC(round(all_death_count), format = "d", big.mark = ",")))
# sane bands for Indonesia RHD (order of magnitude): prevalence 1e5-1e7, deaths 1e3-1e5
if (rhd_prev_count < 1e5 || rhd_prev_count > 1e7)
  stop("Base-year RHD prevalence count ", round(rhd_prev_count),
       " is outside the sane band 1e5-1e7.", call. = FALSE)
if (rhd_death_count < 1e3 || rhd_death_count > 1e5)
  stop("Base-year RHD death count ", round(rhd_death_count),
       " is outside the sane band 1e3-1e5.", call. = FALSE)

# ==============================================================================
# 6. ASSEMBLE + PERSIST THE DISEASE-MODEL INPUT BUNDLE
# ==============================================================================
disease_model_inputs <- list(
  rates_by_year = list(
    ir_rhd   = ir_rhd_yr,     # [age x sex x year] incidence (secular trend baked in)
    mort_rhd = mort_rhd_yr,   # [age x sex x year] RHD cause-specific death rate (held)
    mort_all = mort_all_yr,   # [age x sex x year] all-cause death rate (held)
    oth_mort = oth_mort_yr,   # [age x sex x year] background (non-RHD) mortality (held)
    prev_seed = prev_rhd0     # [age x sex] prevalent RHD fraction — single-year seed
  ),
  clinical   = clinical,      # progression / case-fatality clinical params [tagged]
  effects    = effects,       # intervention effect sizes [PAPER]
  seed_split = seed_split,    # prevalent-pool severity split mild/severe/post [tagged]
  coverage   = coverage,      # coverage-ramp STRUCTURE (baselines/targets/window)
  meta = list(
    LOCATION       = LOCATION,
    AGES           = AGES,
    SEXES          = SEXES,
    years          = years,
    RATE_BASE_YEAR = as.integer(RATE_BASE_YEAR),
    incidence_trend = incidence_trend,
    built_from = c("data-raw/temp_baseline_rates_gbd.rds",
                   basename(POP_PROJ_FILE))
  )
)

OUT_FILE <- paste0(wd_data, "disease_model_inputs.rds")
saveRDS(disease_model_inputs, file = OUT_FILE)

message(sprintf(
  "03_build_disease_model.R complete | wrote %s\n  arrays [age=%d x sex=%d x year=%d]; params: clinical(%d) effects(%d) coverage-struct(%d).",
  basename(OUT_FILE), n_age, length(SEXES), n_years,
  length(clinical), length(effects), length(coverage)))
message("  Next: 04_calibration_random_tp.R (calibrate IR/CF), then 05_build_baseline.R.")
