#===============================================================================
# 04_calibration_random_tp.R
#-------------------------------------------------------------------------------
# RANDOM-SEARCH TP CALIBRATION FOR THE INDONESIA RHD WELL-SICK-DEAD MARKOV MODEL
#
# PURPOSE (reaffirmed): this script does NOT run the model and does NOT run the
# interventions. Its SOLE objective is to CALIBRATE / IMPROVE the transition-
# probability INPUTS (IR incidence and CF case-fatality) within the in-sample
# window 2000-2019, so the well-sick-dead model reproduces the GBD targets BETTER
# than the raw GBD-rate-implied values. The OUTPUT is the calibrated, improved
# transition probabilities (adjusted_searo_part*.rds). The forward projection and
# the reference/SAP scenarios are run downstream (05 prepares the state, 06 runs
# the model). Background mortality is held fixed; only IR and CF are calibrated.
#
# WHAT THIS SCRIPT DOES
# ---------------------
# 1. BUILDS an RHD-native baseline transition-probability (TP) table from
#    upstream, rather than reading the parent NCD model's tps_inpt_part*.rds
#    (which do not exist in this repo). For the FULL 0-95+ age range (paediatric
#    and adult), both sexes, Indonesia, years 2000-2019 it derives, per single
#    age x sex x year:
#        IR        RHD incidence probability (well -> sick)      = GBD RHD Incidence Rate/1e5
#        CF        RHD case-fatality (sick -> dead)              = RHD Deaths / RHD Prevalence
#        PREVt0    prevalent RHD fraction (seed)                 = GBD RHD Prevalence Rate/1e5
#        DIS.mx.t0 RHD death rate per capita                     = GBD RHD Deaths Rate/1e5
#        ALL.mx    all-cause death rate per capita               = GBD All-causes Deaths Rate/1e5
#        BG.mx     background (non-RHD) mortality of the sick     = ALL.mx - DIS.mx.t0  (>=0)
#        BG.mx.all background (non-RHD) mortality of the pool     = ALL.mx - DIS.mx.t0  (>=0)
#        Nx        population                                     = 02's pop_observed table
#    INPUTS:  data-raw/temp_baseline_rates_gbd.rds  (from 01_prepare_inputs.R)
#             data/pop_observed_1990_2024.rds        (from 02_build_demography.R)
#
# 2. CALIBRATES the structural IR and CF (per location-sex-cause, at GRANULARITY)
#    by PURE RANDOM SEARCH over multiplicative adjustment factors, projecting the
#    well-sick-dead cohort 2000 -> 2019 and minimising a weighted RELATIVE squared
#    error against GBD counts. Candidate 0 = baseline (all multipliers = 1), so the
#    calibrated fit is never worse than baseline. The argmin over i.i.d. uniform
#    draws is kept (no Gaussian step / restart -- those are hill-climb constructs).
#
# CALIBRATION TARGETS  (GBD 2023 "Number", aggregated to GBD age groups: the
#                        paediatric groups <1/12-23mo/2-4/5-9/10-14/15-19 plus the
#                        5-year groups 20-24..90-94 and 95+)
# ------------------------------------------------------------------------
#   * PREVALENCE : GBD RHD Prevalence  vs the model sick stock          (RHD-specific)
#   * DEATHS     : GBD ALL-CAUSES Deaths vs the model's DECOMPOSED
#                  all-cause deaths = (RHD case-fatality deaths, sick->dead via CF)
#                                   + (background mortality of the pool, BG.mx.all).
#     i.e. model all-cause deaths = model RHD deaths + background deaths, and THAT
#     sum is compared to the GBD All-causes target (not RHD deaths in isolation).
#
#   BACKGROUND MORTALITY IS HELD FIXED. Only IR and CF are calibrated. BG.mx /
#   BG.mx.all are exogenous (= observed all-cause minus observed RHD deaths) and
#   are NOT free parameters, consistent with enforce_tp_constraints() which
#   PRESERVES BG.mx and treats IR/CF/BG.mx as competing risks summing to <= 1.
#   Consequence (documented): because background is pinned to (all-cause - RHD),
#   the all-cause death term chiefly validates the mortality envelope while the
#   RHD PREVALENCE term drives the IR/CF fit. This is the "RHD deaths + fixed BG
#   must sum to observed all-cause" option.
#
#   error = sum_{year, age.group} [
#               W_DEATHS * ((AllDeaths_model - AllDeaths_gbd)/(AllDeaths_gbd + EPS))^2
#             + W_PREV   * ((Prev_model      - Prev_gbd     )/(Prev_gbd      + EPS))^2 ]
#
# OUTPUT CONTRACT (unchanged)
# ---------------------------
#   * adjusted_searo_part{1..10}.rds in wd_data  -- calibrated TP rows, input schema.
#   * calibration_factors_random_tp.csv          -- per-combo(-age.group) IR/CF multipliers.
#   * calibration_diagnostics_random_tp.csv      -- baseline-vs-calibrated RMSE
#       (RMSE_deaths_* now measure ALL-CAUSE deaths), weighted error, % improvement.
#
# AVOID DOUBLE-CALIBRATION: this script bakes the calibrated multipliers directly
#   into IR/CF. Downstream scripts (05 prepares state, 06 runs the model) must
#   consume these calibrated TPs AS-IS and must NOT re-apply any further
#   adjustment/multiplier on top.
#
# Source AFTER 01_prepare_inputs.R and 02_build_demography.R.
#===============================================================================

library(data.table)
library(foreach)
library(doParallel)
library(parallel)

#===============================================================================
# 0. TUNABLE PARAMETERS  (all honour 00_run_all.R globals via getp(); else default)
#===============================================================================
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

## --- search SPACE -----------------------------------------------------------
SEARCH_HALFWIDTH <- getp("SEARCH_HALFWIDTH", 0.50)  # IR/CF sampled in [1-hw, 1+hw] = [0.5, 1.5]
GRANULARITY <- getp("GRANULARITY", "age_group")     # "combo" (primary) | "age_group" (sensitivity)

## --- search ALGORITHM (PURE RANDOM SEARCH) ----------------------------------
N_ITER       <- getp("N_ITER",       400)           # i.i.d. uniform candidates/combo (baseline = cand 0)
CONVERGE_TOL <- getp("CONVERGE_TOL", 1e-4)          # early stop if best weighted error < this
SEED         <- getp("SEED",         42)            # master seed; per-combo seed = SEED + ci*10000

## --- objective WEIGHTS / numerics -------------------------------------------
W_DEATHS <- getp("W_DEATHS", 2)                     # fatal weight
W_PREV   <- getp("W_PREV",   1)                     # non-fatal weight
EPS_REL  <- 1e-6                                    # denominator floor for relative error

## --- probability-constraint numerics ----------------------------------------
TP_EPS <- 0.005                                     # buffer kept below 1 when capping/renormalising

## --- calibration target window (run the cohort 2000 -> 2019) ----------------
CAL_YEAR_START <- as.integer(getp("CAL_YEAR_START", 2000))
CAL_YEAR_END   <- as.integer(getp("CAL_YEAR_END",   2019))

## --- RHD baseline TP build ---------------------------------------------------
LOCATION  <- getp("LOCATION", "Indonesia")
AGE_LO    <- as.integer(getp("AGE_LO", 0L))  # FULL cohort age range: 0..95+ (paediatric + adult).
AGE_HI    <- as.integer(getp("AGE_HI", 95L)) # AGE_LO is the birth/entrant age (see project_combo).
RHD_CAUSE <- "Rheumatic heart disease"
ALL_CAUSE <- "All causes"

## --- execution --------------------------------------------------------------
RUN_PAR   <- if (exists("run_calibration_par")) isTRUE(run_calibration_par) else TRUE
MAX_CORES <- as.integer(getp("MAX_CORES", 14L))
N_OUT_CHUNKS <- 10

## --- paths (honour 00_run_all.R globals; else here()) -----------------------
if (!exists("wd_raw"))  wd_raw  <- paste0(here::here("data-raw"), "/")
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

#===============================================================================
# 1. BUILD RHD-NATIVE BASELINE TP TABLE  (replaces tps_inpt_part*.rds)
#===============================================================================

## single-year age (0-95+) -> GBD age-group label (matches make_age_match).
## MIRRORS the mapping in 03_build_disease_model.R, INCLUDING the paediatric /
## young groups (<1, 12-23mo, 2-4, 5-9, 10-14, 15-19) so the full 0-95+ range is
## calibrated. (Previously this started at "20-24 years", which silently binned
## every age < 25 into the 20-24 group.)
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

build_rhd_tps <- function(gbd, pop, loc, y0, y1) {
  ages <- AGE_LO:AGE_HI
  am   <- data.table(age = ages, age_group = age_to_gbd_group(ages))

  ## GBD per-capita rates by group, RHD + All-causes, calibration window
  r <- gbd[location_name == loc & metric_name == "Rate" &
             year >= y0 & year <= y1 &
             cause_name %in% c(RHD_CAUSE, ALL_CAUSE),
           .(sex = sex_name, age_group = age_name, year,
             cause = cause_name, measure = measure_name, rate = val / 1e5)]
  r <- r[age_group %in% am$age_group]
  r[, cm := paste0(fifelse(cause == RHD_CAUSE, "rhd", "all"), "_",
                   fcase(measure == "Incidence", "inc",
                         measure == "Prevalence", "prev",
                         measure == "Deaths",     "dth",
                         default = "oth"))]
  w <- dcast(r, sex + age_group + year ~ cm, value.var = "rate")
  for (col in c("rhd_inc", "rhd_prev", "rhd_dth", "all_dth"))
    if (!col %in% names(w)) w[, (col) := 0]
  w <- am[w, on = "age_group", allow.cartesian = TRUE]     # expand groups -> single age

  ## population (single age x sex x year) from 02's observed table
  pp <- pop[location == loc & age %in% ages & year >= y0 & year <= y1,
            .(location, sex, age, year, Nx)]

  dt <- merge(w, pp, by = c("sex", "age", "year"), all.x = TRUE)
  ## fill any missing rate/pop cells with 0 (rare young/old gaps)
  for (col in c("rhd_inc", "rhd_prev", "rhd_dth", "all_dth"))
    dt[is.na(get(col)), (col) := 0]

  dt[, `:=`(
    cause     = RHD_CAUSE,
    PREVt0    = rhd_prev,
    DIS.mx.t0 = rhd_dth,
    ALL.mx    = all_dth,
    IR        = rhd_inc,
    CF        = fifelse(rhd_prev > 0, rhd_dth / rhd_prev, 0),
    BG.mx     = pmax(all_dth - rhd_dth, 0),
    BG.mx.all = pmax(all_dth - rhd_dth, 0)
  )]

  out <- dt[!is.na(Nx),
            .(age, sex, location, year, cause,
              BG.mx.all, ALL.mx, BG.mx, PREVt0, DIS.mx.t0, Nx, IR, CF)]
  setkey(out, location, sex, cause, year, age)
  out[]
}

gbd_raw <- as.data.table(readRDS(paste0(wd_raw, "temp_baseline_rates_gbd.rds")))
pop_obs <- as.data.table(readRDS(paste0(wd_data, "pop_observed_1990_2024.rds")))

## the input population MUST actually cover the 2000-onward calibration window
pop_years <- range(pop_obs[location == LOCATION, year])
if (pop_years[1] > CAL_YEAR_START)
  stop(sprintf("Population input starts in %d but calibration needs %d onward. ",
               pop_years[1], CAL_YEAR_START),
       "Re-run 02_build_demography.R with an earlier OBS_YEARS start.", call. = FALSE)
gbd_years <- range(gbd_raw[cause_name == RHD_CAUSE, year])
if (gbd_years[1] > CAL_YEAR_START)
  stop(sprintf("GBD RHD data starts in %d but calibration needs %d onward.",
               gbd_years[1], CAL_YEAR_START), call. = FALSE)

b_rates <- build_rhd_tps(gbd_raw, pop_obs, LOCATION, CAL_YEAR_START, CAL_YEAR_END)
locs    <- unique(b_rates$location)

## defensive clamps (as the parent 031 applied before calibrating)
b_rates[CF >= 1, CF := 0.99]
b_rates[IR >= 1, IR := 0.99]
b_rates[CF < 0,  CF := 0]
b_rates[IR < 0,  IR := 0]

## frozen copy of the INPUT for end-of-run schema / row-count validation
tps_input_cols <- copy(names(b_rates))
tps_input_nrow <- nrow(b_rates)

cat(sprintf("Built RHD baseline TPs: %d rows | ages %d-%d | years %d-%d | sexes %s\n",
            nrow(b_rates), min(b_rates$age), max(b_rates$age),
            min(b_rates$year), max(b_rates$year),
            paste(unique(b_rates$sex), collapse = ", ")))

## incoming AGE_LO (=age 0, births) population each calibration year, from 02's
## observed table. With AGE_LO = 0 these are the birth-cohort entrants that refresh
## the youngest age each year (previously age-20 entrants when the range began at 20).
pop_ent <- pop_obs[location %in% locs & age == AGE_LO &
                     year >= CAL_YEAR_START & year <= CAL_YEAR_END,
                   .(location, sex, age, year, Nx_ent = Nx)]

#===============================================================================
# 2. GBD CALIBRATION TARGETS
#    Prevalence target = GBD RHD Prevalence (Number).
#    Deaths     target = GBD ALL-CAUSES Deaths (Number)  <-- decomposed match.
#    Both by GBD age group x sex x year, calibration window, LOCATION.
#===============================================================================

gbd <- copy(gbd_raw)
setnames(gbd,
         c("sex_name", "age_name", "cause_name", "measure_name", "metric_name", "location_name"),
         c("sex",      "age",      "cause",      "measure",      "metric",      "location"))
gbd <- gbd[metric == "Number" & location %in% locs &
             year >= CAL_YEAR_START & year <= CAL_YEAR_END]

age_groups_keep <- unique(age_to_gbd_group(AGE_LO:AGE_HI))

prev_t <- gbd[cause == RHD_CAUSE & measure == "Prevalence" & age %in% age_groups_keep,
              .(location, sex, age, year, gbdPrev = val)]
alld_t <- gbd[cause == ALL_CAUSE & measure == "Deaths" & age %in% age_groups_keep,
              .(location, sex, age, year, gbdAllDeaths = val)]

targets <- merge(prev_t, alld_t, by = c("location", "sex", "age", "year"), all = TRUE)
targets[is.na(gbdPrev),      gbdPrev := 0]
targets[is.na(gbdAllDeaths), gbdAllDeaths := 0]
targets[, cause := RHD_CAUSE]                 # single disease cause in this model
setkey(targets, location, sex, cause)

if (nrow(targets) == 0)
  stop("No GBD calibration targets built -- check dx_include / age groups / years.", call. = FALSE)

#===============================================================================
# 3. HELPER FUNCTIONS
#===============================================================================

## GBD age-group labels for the full 0-95+ range (matches age_to_gbd_group).
make_age_match <- function() {
  am <- data.table(age = AGE_LO:AGE_HI)
  am[, age.group := age_to_gbd_group(age)]
  am
}

## draw ONE i.i.d. uniform multiplier vector (heart of the random search).
perturb_cvd_combo_random <- function(n_age_groups, granularity, lo, hi, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_par <- if (granularity == "age_group") 2L * n_age_groups else 2L
  runif(n_par, lo, hi)
}

## enforce probability / competing-risk constraints (in place). PRESERVES BG.mx
## (only capping the disease TP) except in the rare fallback where BG.mx alone
## leaves no room; those rows are flagged in bg_modified.
enforce_tp_constraints <- function(dt, tp_eps = TP_EPS) {
  dt[is.na(IR),    IR := 0]
  dt[is.na(CF),    CF := 0]
  dt[is.na(BG.mx), BG.mx := 0]
  dt[IR < 0, IR := 0]; dt[IR > 1, IR := 1]
  dt[CF < 0, CF := 0]; dt[CF > 1, CF := 1]
  dt[BG.mx < 0, BG.mx := 0]

  if (!("bg_modified" %in% names(dt))) dt[, bg_modified := 0L]

  dt[, headroom := 1 - BG.mx - tp_eps]
  dt[headroom >= 0 & IR > headroom, IR := headroom]
  dt[headroom >= 0 & CF > headroom, CF := headroom]

  dt[headroom < 0 & (IR + BG.mx) > 1, `:=`(
    IR_new   = IR    / (IR + BG.mx) - tp_eps,
    BGmx_new = BG.mx / (IR + BG.mx) - tp_eps,
    bg_modified = 1L)]
  dt[!is.na(IR_new), `:=`(IR = pmax(IR_new, 0), BG.mx = pmax(BGmx_new, 0))]
  dt[, c("IR_new", "BGmx_new") := NULL]

  dt[(CF + BG.mx) > 1, `:=`(
    CF_new   = CF    / (CF + BG.mx) - tp_eps,
    BGmx_new = BG.mx / (CF + BG.mx) - tp_eps,
    bg_modified = 1L)]
  dt[!is.na(CF_new), `:=`(CF = pmax(CF_new, 0), BG.mx = pmax(BGmx_new, 0))]
  dt[, c("CF_new", "BGmx_new") := NULL]

  dt[, headroom := NULL]
  dt[]
}

## build multiplier table from a flat parameter vector.
build_mtab <- function(par, age_groups, granularity) {
  if (granularity == "age_group") {
    n <- length(age_groups)
    data.table(age.group = age_groups, ir_mult = par[1:n], cf_mult = par[(n + 1):(2 * n)])
  } else {
    data.table(age.group = age_groups, ir_mult = par[1], cf_mult = par[2])
  }
}

## apply multipliers to a combo's TP rows, then enforce constraints (new DT).
apply_multipliers <- function(combo_rates, mtab, age_match) {
  cr <- copy(combo_rates)
  cr <- merge(cr, age_match, by = "age",       all.x = TRUE)
  cr <- merge(cr, mtab,      by = "age.group", all.x = TRUE)
  cr[is.na(ir_mult), ir_mult := 1]
  cr[is.na(cf_mult), cf_mult := 1]
  cr[, IR := IR * ir_mult]
  cr[, CF := CF * cf_mult]
  cr[, c("age.group", "ir_mult", "cf_mult") := NULL]
  enforce_tp_constraints(cr)
  cr
}

## project ONE combo through the well-sick-dead recursion, y0 -> y1.
## Returns sick (model prevalence), dead (model RHD deaths) and all.mx
## (model ALL-CAUSE deaths = RHD CF deaths + background of the pool).
project_combo <- function(cr, pop_combo, y0 = CAL_YEAR_START, y1 = CAL_YEAR_END) {
  br <- merge(cr, pop_combo, by = c("year", "location", "sex", "age"), all.x = TRUE)
  br[age == AGE_LO & year > y0, Nx := Nx_ent]    # refresh age-0 (birth) entrants each year
  br[, Nx_ent := NULL]

  br[year == y0 | age == AGE_LO, sick   := Nx * PREVt0]
  br[year == y0 | age == AGE_LO, dead   := Nx * DIS.mx.t0]
  br[year == y0 | age == AGE_LO, well   := Nx * (1 - (PREVt0 + ALL.mx))]
  br[year == y0 | age == AGE_LO, pop    := Nx]
  br[year == y0 | age == AGE_LO, all.mx := Nx * ALL.mx]

  br[CF > 0.9, CF := 0.9]
  br[IR > 0.9, IR := 0.9]

  n_steps <- y1 - y0
  for (s in 1:n_steps) {
    yr <- y0 + s
    b2 <- br[year <= yr & year >= yr - 1]
    setorder(b2, sex, location, cause, age, year)
    b2[, age2 := age + 1]

    b2[, sick2 := shift(sick) * (1 - (CF + BG.mx)) + shift(well) * IR,
       by = .(sex, location, cause, age)]
    b2[sick2 < 0, sick2 := 0]
    b2[, dead2 := shift(sick) * CF, by = .(sex, location, cause, age)]
    b2[dead2 < 0, dead2 := 0]
    b2[, pop2 := shift(pop) - shift(all.mx), by = .(sex, location, cause, age)]
    b2[pop2 < 0, pop2 := 0]
    ## all-cause deaths = disease (CF) deaths this combo + background mortality of pool
    b2[, all.mx2 := sum(dead2), by = .(sex, location, year, age)]
    b2[, all.mx2 := all.mx2 + (pop2 * BG.mx.all)]
    b2[all.mx2 < 0, all.mx2 := 0]
    b2[, well2 := pop2 - all.mx2 - sick2]
    b2[well2 < 0, well2 := 0]

    upd <- b2[year == yr & age2 <= AGE_HI,
              .(age = age2, year, sick2, dead2, well2, pop2, all.mx2)]
    br[upd, on = .(year, age), `:=`(
      sick = i.sick2, dead = i.dead2, well = i.well2,
      pop = i.pop2, all.mx = i.all.mx2)]
  }

  br[year >= y0, .(location, sex, cause, year, age, sick, dead, all.mx)]
}

## aggregate single-age projection to GBD age groups and join GBD targets.
## model Prevalence = sum(sick); model AllDeaths = sum(all.mx) (decomposed).
proj_vs_targets <- function(proj, combo_targets, age_match) {
  m  <- merge(proj, age_match, by = "age", all.x = TRUE)
  ms <- m[, .(Prevalence = sum(sick,   na.rm = TRUE),
              AllDeaths  = sum(all.mx, na.rm = TRUE)),
          by = .(location, sex, cause, year, age.group)]
  setnames(ms, "age.group", "age")
  j <- merge(combo_targets, ms,
             by = c("location", "sex", "cause", "year", "age"), all.x = TRUE)
  j[is.na(AllDeaths),  AllDeaths := 0]
  j[is.na(Prevalence), Prevalence := 0]
  j
}

## weighted RELATIVE squared error (search objective). All-cause deaths weighted 2x.
combo_error <- function(proj, combo_targets, age_match,
                        w_deaths = W_DEATHS, w_prev = W_PREV, eps = EPS_REL) {
  j <- proj_vs_targets(proj, combo_targets, age_match)
  j[, sum(
    w_deaths * ((AllDeaths  - gbdAllDeaths) / (gbdAllDeaths + eps))^2 +
    w_prev   * ((Prevalence - gbdPrev)      / (gbdPrev      + eps))^2,
    na.rm = TRUE)]
}

## absolute RMSE diagnostics per location-sex-cause-age.group.
## RMSE_deaths now measures ALL-CAUSE deaths (the decomposed target).
combo_diag <- function(proj, combo_targets, age_match) {
  j <- proj_vs_targets(proj, combo_targets, age_match)
  j[, .(
    RMSE_deaths = sqrt(mean((AllDeaths  - gbdAllDeaths)^2, na.rm = TRUE)),
    RMSE_prev   = sqrt(mean((Prevalence - gbdPrev)^2,      na.rm = TRUE))
  ), by = .(location, sex, cause, age)]
}

## calibrate ONE combo: PURE RANDOM SEARCH; candidate 0 = baseline.
calibrate_one_combo_random <- function(combo_rates, pop_combo, combo_targets, age_match,
                                       granularity, hw, n_iter, converge_tol,
                                       w_deaths, w_prev, eps, seed) {
  lo <- 1 - hw; hi <- 1 + hw
  age_groups <- sort(unique(combo_targets$age))
  n_g   <- length(age_groups)

  eval_par <- function(par) {
    mtab <- build_mtab(par, age_groups, granularity)
    cr2  <- apply_multipliers(combo_rates, mtab, age_match)
    proj <- project_combo(cr2, pop_combo)
    combo_error(proj, combo_targets, age_match, w_deaths, w_prev, eps)
  }

  n_par    <- if (granularity == "age_group") 2L * n_g else 2L
  best_par <- rep(1, n_par)
  best_err <- eval_par(best_par)
  base_err <- best_err
  n_eval   <- 1L

  for (it in 1:n_iter) {
    cand <- perturb_cvd_combo_random(n_g, granularity, lo, hi, seed = seed + it)
    err  <- eval_par(cand); n_eval <- n_eval + 1L
    if (err < best_err) { best_err <- err; best_par <- cand }
    if (best_err < converge_tol) break
  }

  list(mtab = build_mtab(best_par, age_groups, granularity),
       best_err = best_err, base_err = base_err,
       n_eval = n_eval, n_par = n_par,
       hit_bound = any(best_par <= lo + 1e-9 | best_par >= hi - 1e-9))
}

## drive one combo end-to-end.
run_combo <- function(ci, combos, b_rates, pop_ent, targets, age_match) {
  loc <- combos$location[ci]; sx <- combos$sex[ci]; cse <- combos$cause[ci]

  cr <- b_rates[location == loc & sex == sx & cause == cse]
  pc <- pop_ent[location == loc & sex == sx]
  ct <- targets[location == loc & sex == sx & cause == cse]

  base_rows <- enforce_tp_constraints(copy(cr))

  if (nrow(ct) == 0) {
    base_rows[, bg_modified := NULL]
    return(list(
      rows = base_rows,
      factors = data.table(location = loc, sex = sx, cause = cse,
                           age.group = NA_character_, ir_mult = 1, cf_mult = 1,
                           granularity = GRANULARITY),
      diag = data.table(location = loc, sex = sx, cause = cse, age = NA_character_,
                        RMSE_deaths_base = NA_real_, RMSE_prev_base = NA_real_,
                        RMSE_deaths_cal = NA_real_,  RMSE_prev_cal = NA_real_),
      err = data.table(location = loc, sex = sx, cause = cse,
                       base_err = NA_real_, cal_err = NA_real_,
                       n_eval = 0L, n_par = 0L, hit_bound = FALSE,
                       bg_modified_rows = 0L)))
  }

  fit <- calibrate_one_combo_random(cr, pc, ct, age_match,
                                    GRANULARITY, SEARCH_HALFWIDTH, N_ITER,
                                    CONVERGE_TOL, W_DEATHS, W_PREV, EPS_REL,
                                    SEED + ci * 10000L)

  cal_rows  <- apply_multipliers(cr, fit$mtab, age_match)
  bg_mod_n  <- sum(cal_rows$bg_modified)
  cal_rows[, bg_modified := NULL]

  base_diag <- combo_diag(project_combo(base_rows, pc), ct, age_match)
  cal_diag  <- combo_diag(project_combo(cal_rows,  pc), ct, age_match)
  setnames(base_diag, c("RMSE_deaths", "RMSE_prev"), c("RMSE_deaths_base", "RMSE_prev_base"))
  setnames(cal_diag,  c("RMSE_deaths", "RMSE_prev"), c("RMSE_deaths_cal",  "RMSE_prev_cal"))
  diag <- merge(base_diag, cal_diag, by = c("location", "sex", "cause", "age"), all = TRUE)

  factors <- copy(fit$mtab)
  factors[, `:=`(location = loc, sex = sx, cause = cse, granularity = GRANULARITY)]
  setcolorder(factors, c("location", "sex", "cause", "age.group",
                         "ir_mult", "cf_mult", "granularity"))

  err <- data.table(location = loc, sex = sx, cause = cse,
                    base_err = fit$base_err, cal_err = fit$best_err,
                    n_eval = fit$n_eval, n_par = fit$n_par,
                    hit_bound = fit$hit_bound, bg_modified_rows = bg_mod_n)

  list(rows = cal_rows, factors = factors, diag = diag, err = err)
}

#===============================================================================
# 4. RUN CALIBRATION OVER ALL location-sex-cause COMBOS
#===============================================================================

age_match <- make_age_match()
combos    <- unique(b_rates[, .(location, sex, cause)])
n_combos  <- nrow(combos)

cat(sprintf("Random-search TP calibration: %d combos | granularity = %s | window %d-%d\n",
            n_combos, GRANULARITY, CAL_YEAR_START, CAL_YEAR_END))
cat(sprintf("Search range per multiplier: [%.2f, %.2f] | %d i.i.d. candidates/combo\n",
            1 - SEARCH_HALFWIDTH, 1 + SEARCH_HALFWIDTH, N_ITER))

worker_exports <- c(
  "combos", "b_rates", "pop_ent", "targets", "age_match",
  "make_age_match", "perturb_cvd_combo_random", "enforce_tp_constraints",
  "build_mtab", "apply_multipliers", "project_combo", "proj_vs_targets",
  "combo_error", "combo_diag", "calibrate_one_combo_random", "run_combo",
  "age_to_gbd_group", "GRANULARITY", "SEARCH_HALFWIDTH", "N_ITER", "CONVERGE_TOL",
  "W_DEATHS", "W_PREV", "EPS_REL", "TP_EPS",
  "CAL_YEAR_START", "CAL_YEAR_END", "AGE_LO", "AGE_HI", "SEED"
)

if (RUN_PAR && n_combos > 1) {
  n_cores <- max(1L, min(MAX_CORES, n_combos, parallel::detectCores() - 1L))
  cat(sprintf("Running in parallel on %d cores...\n", n_cores))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  results <- foreach(ci = seq_len(n_combos),
                     .packages = c("data.table"),
                     .export   = worker_exports) %dopar% {
    setDTthreads(1)
    run_combo(ci, combos, b_rates, pop_ent, targets, age_match)
  }
  stopCluster(cl)
} else {
  cat("Running sequentially...\n")
  results <- lapply(seq_len(n_combos), function(ci) {
    res <- run_combo(ci, combos, b_rates, pop_ent, targets, age_match)
    cat(sprintf("  [%d/%d] %s | %s | %s : err %.3g -> %.3g\n",
                ci, n_combos, combos$location[ci], combos$sex[ci], combos$cause[ci],
                res$err$base_err, res$err$cal_err))
    res
  })
}

## --- consolidate ------------------------------------------------------------
calibrated  <- rbindlist(lapply(results, `[[`, "rows"),    use.names = TRUE, fill = TRUE)
factors_out <- rbindlist(lapply(results, `[[`, "factors"), use.names = TRUE, fill = TRUE)
diag_out    <- rbindlist(lapply(results, `[[`, "diag"),    use.names = TRUE, fill = TRUE)
err_out     <- rbindlist(lapply(results, `[[`, "err"),     use.names = TRUE, fill = TRUE)

setcolorder(calibrated, intersect(tps_input_cols, names(calibrated)))

#===============================================================================
# 5. WRITE OUTPUTS  (adjusted_searo_part{1..10}.rds in wd_data)
#===============================================================================

n     <- nrow(calibrated)
chunk <- ceiling(n / N_OUT_CHUNKS)
for (i in 1:N_OUT_CHUNKS) {
  start <- (i - 1) * chunk + 1
  end   <- min(i * chunk, n)
  if (start > n) {
    saveRDS(calibrated[0], file = paste0(wd_data, "adjusted_searo_part", i, ".rds"))
    next
  }
  saveRDS(calibrated[start:end], file = paste0(wd_data, "adjusted_searo_part", i, ".rds"))
}

fwrite(factors_out, paste0(wd_data, "calibration_factors_random_tp.csv"))

err_pct <- copy(err_out[, .(location, sex, cause, base_err, cal_err,
                            n_eval, n_par, hit_bound, bg_modified_rows)])
err_pct[, pct_improvement := 100 * (base_err - cal_err) / pmax(base_err, EPS_REL)]
diag_full <- merge(diag_out, err_pct, by = c("location", "sex", "cause"), all.x = TRUE)
fwrite(diag_full, paste0(wd_data, "calibration_diagnostics_random_tp.csv"))

#===============================================================================
# 6. VALIDATION  (hard failures; identify the offending combo/age on violation)
#===============================================================================

cat("\n", strrep("=", 70), "\nVALIDATION\n", strrep("=", 70), "\n", sep = "")

## locate any offending rows and stop() naming them (stronger than a bare stopifnot)
report_bad <- function(dt, bad_idx, msg) {
  bad_idx[is.na(bad_idx)] <- TRUE                 # NA in the test itself counts as bad
  if (any(bad_idx)) {
    bad <- dt[bad_idx]
    cat("OFFENDING ROWS (", msg, "):\n", sep = "")
    print(utils::head(bad[, .(location, sex, cause, age, year, IR, CF, BG.mx)], 10))
    stop(msg, " -- ", sum(bad_idx), " offending row(s); see above.", call. = FALSE)
  }
}
report_bad(calibrated, calibrated[, is.na(IR)],                "IR is NA")
report_bad(calibrated, calibrated[, is.na(CF)],                "CF is NA")
report_bad(calibrated, calibrated[, is.na(BG.mx)],             "BG.mx is NA")
report_bad(calibrated, calibrated[, IR < 0 | IR > 1],          "IR outside [0,1]")
report_bad(calibrated, calibrated[, CF < 0 | CF > 1],          "CF outside [0,1]")
report_bad(calibrated, calibrated[, BG.mx < 0 | BG.mx > 1],    "BG.mx outside [0,1]")
report_bad(calibrated, calibrated[, IR + BG.mx > 1 + 1e-9],    "IR + BG.mx > 1 (competing risk)")
report_bad(calibrated, calibrated[, CF + BG.mx > 1 + 1e-9],    "CF + BG.mx > 1 (competing risk)")

stopifnot(
  "row count != input" = nrow(calibrated) == tps_input_nrow,
  "schema != input"    = setequal(names(calibrated), tps_input_cols)
)
cat("All probability/row constraints satisfied (IR, CF, BG.mx in [0,1]; competing sums <= 1).\n")
cat(sprintf("Rows: %d (matches input: %d). Schema matches input: TRUE.\n",
            nrow(calibrated), tps_input_nrow))

bg_rows_total <- sum(err_out$bg_modified_rows, na.rm = TRUE)
if (bg_rows_total > 0) {
  cat(sprintf("NOTE: BG.mx renormalised (fallback) on %d rows where BG.mx alone left ",
              bg_rows_total), "no room for the disease TP.\n", sep = "")
} else {
  cat("BG.mx preserved on ALL rows (no fallback renormalisation).\n")
}

tot_base_err <- sum(err_out$base_err, na.rm = TRUE)
tot_cal_err  <- sum(err_out$cal_err,  na.rm = TRUE)
cat(sprintf("\nWeighted relative error (search objective), summed over combos:\n"))
cat(sprintf("  baseline = %.4g   calibrated = %.4g   reduction = %.1f%%\n",
            tot_base_err, tot_cal_err,
            100 * (tot_base_err - tot_cal_err) / max(tot_base_err, EPS_REL)))

abs_summary <- diag_out[, .(
  RMSE_deaths_base = mean(RMSE_deaths_base, na.rm = TRUE),
  RMSE_deaths_cal  = mean(RMSE_deaths_cal,  na.rm = TRUE),
  RMSE_prev_base   = mean(RMSE_prev_base,   na.rm = TRUE),
  RMSE_prev_cal    = mean(RMSE_prev_cal,    na.rm = TRUE))]
cat("\nMean absolute RMSE across combo x age.group cells (counts):\n")
cat(sprintf("  All-cause deaths : baseline = %.1f -> calibrated = %.1f\n",
            abs_summary$RMSE_deaths_base, abs_summary$RMSE_deaths_cal))
cat(sprintf("  RHD prevalence   : baseline = %.1f -> calibrated = %.1f\n",
            abs_summary$RMSE_prev_base, abs_summary$RMSE_prev_cal))

cat("\nPer-combo weighted error (baseline -> calibrated):\n")
print(err_out[order(cause, sex),
              .(location, sex, cause,
                base_err = round(base_err, 3), cal_err = round(cal_err, 3),
                n_eval, hit_bound, bg_modified_rows)])

if (any(err_out$hit_bound, na.rm = TRUE))
  cat("\nWARNING: some combos hit the search bound -- consider widening ",
      "SEARCH_HALFWIDTH.\n", sep = "")

cat("\nWrote:\n")
cat(sprintf("  %sadjusted_searo_part{1..%d}.rds\n", wd_data, N_OUT_CHUNKS))
cat(sprintf("  %scalibration_factors_random_tp.csv\n", wd_data))
cat(sprintf("  %scalibration_diagnostics_random_tp.csv\n", wd_data))
cat("\nReminder: 05/06 consume these calibrated TPs AS-IS (no further adjustment).\n")
