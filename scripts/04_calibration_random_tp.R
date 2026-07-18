#===============================================================================
# 04_calibration_random_tp.R
#-------------------------------------------------------------------------------
# LAYER-1 (AGGREGATE) CALIBRATION FOR THE RHD A/B/C/D STATE-TRANSITION MODEL
#
# This script calibrates the aggregate well-sick-dead PROXY that feeds the
# A/B/C/D engine (06): it fixes the incidence age-sex pattern (drives the inflow
# into stage A) and an aggregate RHD case-fatality anchor (consumed by the
# Stage-2 structural stage calibration). It REPLACES the former pure random
# search (400 i.i.d. uniform multipliers) against GBD ALL-CAUSE deaths with a
# low-dimensional, prior-anchored, hazard-based calibration against GBD
# RHD-SPECIFIC incidence, prevalence and deaths.
#
# WHAT CHANGED vs the old random search (root cause it fixes)
# -----------------------------------------------------------
#   * TARGET   : GBD RHD-SPECIFIC deaths (not all-cause). All-cause deaths are
#                retained for VALIDATION only, never in the loss.
#   * LOSS     : normalized log-scale squared error on incidence + prevalence +
#                RHD deaths, with OPTIONAL inverse-variance weights (GBD 95% UIs),
#                plus Gaussian priors and 2nd-difference smoothness on the
#                (log-scale) multipliers.
#   * RISK     : competing risks are combined on the HAZARD scale
#                (p_any = 1 - exp[-(h_a+h_b)]; split by hazard share), so
#                probabilities never need clipping and survival is a proper
#                residual.
#   * DIM      : instead of ~52 independent per-age-group multipliers found by
#                random draws, ONE incidence correction (alpha) and ONE aggregate
#                RHD-mortality correction (beta) per BROAD AGE BAND x sex, so the
#                calibrated IR is smooth (not ragged).
#   * ANCHOR   : incidence is anchored to the GBD age-sex pattern,
#                IR_model = IR_GBD x exp(f_band(a)), f prior-anchored to 0.
#   * OPTIMISER: bounded multi-start L-BFGS-B (base R; DEfault), or nloptr/DEoptim
#                if selected/installed — NOT 400 i.i.d. random vectors. The
#                baseline (all multipliers = 1) is always evaluated so the
#                penalized objective is never worse than baseline.
#   * SEQUENCE : (1) mortality-only vs GBD RHD deaths -> (2) incidence-only vs GBD
#                incidence + prevalence -> (3) joint refinement, each with saved
#                diagnostics.
#
# LAYER 2 (A/B/C/D STRUCTURAL / STAGE calibration) is the Stage-2 concern and
# lives in 04b_calibrate_structural.R; this script preserves its interface
# ($stage_calibration) so the pipeline runs whether or not 04b has been wired.
#
# OUTPUT CONTRACT (single self-describing bundle — schema PRESERVED for 05/tests)
# ------------------------------------------------------------------------------
#   data/<COUNTRY>/calibrated_rhd_parameters.rds  — named list:
#     $tp                : calibrated TP data.table (IR, CF, BG.mx, ... by
#                          location x sex x cause x year x age) — same schema as
#                          before; IR = IR_GBD x alpha, CF = CF_GBD x beta.
#     $factors           : per-combo x band alpha/beta multipliers (+ age.group map).
#     $diagnostics       : baseline-vs-calibrated loss + RMSE per combo.
#     $stage_calibration : Layer-2 interface (unchanged; Stage-2 fills it).
#     $layer1            : NEW — incidence_parameters, mortality_parameters,
#                          objective_components, optimizer_diagnostics, validation,
#                          mass_balance, fit_by_group_year.
#     $meta              : window, calib_last_year, bands, weights, optimiser, ...
#   plus human-readable CSVs (factors / diagnostics / fit / mass-balance).
#
# AVOID DOUBLE-CALIBRATION: the calibrated multipliers are baked ONCE into IR/CF.
#   05 prepares the state and 06 runs the model consuming calibrated IR/CF AS-IS.
#
# Source AFTER 01_prepare_inputs.R and 02_build_demography.R.
#===============================================================================

library(data.table)

#===============================================================================
# 0. TUNABLE PARAMETERS  (all honour 00_run_all.R globals via getp(); else default)
#===============================================================================
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

## --- objective weights + numerics -------------------------------------------
W_INC   <- getp("W_INC",   1)      # weight on RHD incidence fit
W_PREV  <- getp("W_PREV",  1)      # weight on RHD prevalence fit
W_DEATH <- getp("W_DEATH", 2)      # weight on RHD-SPECIFIC deaths fit (all-cause is validation-only)
USE_IV_WEIGHTS   <- isTRUE(getp("USE_IV_WEIGHTS", TRUE))
TARGET_MIN_COUNT <- getp("TARGET_MIN_COUNT", 1)   # drop cells below this GBD Number from the LOSS
EPS_LOG <- 1e-3                    # floor inside log() for the relative log-loss

## --- incidence calibration mode + broad bands -------------------------------
INCIDENCE_CALIBRATION_MODE <- getp("INCIDENCE_CALIBRATION_MODE", "anchored")  # fixed|anchored|free
CALIB_AGE_BANDS <- as.integer(getp("CALIB_AGE_BANDS", c(0L, 15L, 25L, 45L, 65L)))

## --- priors + smoothness (log-multiplier scale) -----------------------------
SIGMA_ALPHA   <- getp("SIGMA_ALPHA",   0.75)
SIGMA_BETA    <- getp("SIGMA_BETA",    0.50)
LAMBDA_PRIOR  <- getp("LAMBDA_PRIOR",  1.0)
LAMBDA_SMOOTH <- getp("LAMBDA_SMOOTH", 1.0)
MULT_LO       <- getp("MULT_LO", 0.2)
MULT_HI       <- getp("MULT_HI", 5.0)
# "free" mode loosens the incidence prior (x8) so alpha can move; "fixed" pins alpha=1.
SIGMA_ALPHA_EFF <- if (INCIDENCE_CALIBRATION_MODE == "free") SIGMA_ALPHA * 8 else SIGMA_ALPHA

## --- optimiser --------------------------------------------------------------
CALIB_OPTIMIZER <- getp("CALIB_OPTIMIZER", "multistart")  # multistart|nloptr|deoptim
N_STARTS        <- as.integer(getp("N_STARTS", 8L))
CONVERGE_TOL    <- getp("CONVERGE_TOL", 1e-8)
SEED            <- as.integer(getp("SEED", 42L))

## --- calibration window + age range -----------------------------------------
CAL_YEAR_START <- as.integer(getp("CAL_YEAR_START", 2000))
CAL_YEAR_END   <- as.integer(getp("CAL_YEAR_END",   2019))
LOCATION  <- getp("LOCATION", "Indonesia")
AGE_LO    <- as.integer(getp("AGE_LO", 0L))
AGE_HI    <- as.integer(getp("AGE_HI", 95L))
RHD_CAUSE <- "Rheumatic heart disease"
ALL_CAUSE <- "All causes"
TP_EPS    <- 0.005                 # buffer kept below 1 when capping IR/CF vs background

## --- paths (honour 00_run_all.R globals; else here()) -----------------------
if (!exists("wd_raw"))  wd_raw  <- paste0(here::here("data-raw"), "/")
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

OUT_BUNDLE  <- paste0(wd_data, "calibrated_rhd_parameters.rds")
STAGE_TMPL  <- paste0(wd_data, "calibration_targets_stage_template.csv")
STAGE_TARGET_CANDIDATES <- c(
  paste0(wd_data, "calibration_targets_stage.csv"),
  paste0(wd_raw,  "calibration_targets_stage.csv")
)
DISEASE_INPUTS_FILE <- paste0(wd_data, "disease_model_inputs.rds")

#===============================================================================
# 1. AGE-GROUP + BAND MAPPINGS
#===============================================================================
## single-year age (0-95+) -> GBD age-group label (MIRRORS 03/report age_to_grp).
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

## single-year age -> broad calibration band index (1..n_band) + labels.
band_edges <- sort(unique(as.integer(CALIB_AGE_BANDS)))
n_band     <- length(band_edges)
band_labels <- vapply(seq_len(n_band), function(i) {
  lo <- band_edges[i]
  hi <- if (i < n_band) band_edges[i + 1] - 1L else NA_integer_
  if (is.na(hi)) sprintf("%d+", lo) else sprintf("%d-%d", lo, hi)
}, character(1))
band_of <- function(a) findInterval(a, band_edges)   # 1..n_band

#===============================================================================
# 2. BUILD RHD-NATIVE BASELINE RATE TABLE  (per single age x sex x year)
#===============================================================================
build_rhd_tps <- function(gbd, pop, loc, y0, y1) {
  ages <- AGE_LO:AGE_HI
  am   <- data.table(age = ages, age_group = age_to_gbd_group(ages))

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
  w <- am[w, on = "age_group", allow.cartesian = TRUE]

  pp <- pop[location == loc & age %in% ages & year >= y0 & year <= y1,
            .(location, sex, age, year, Nx)]
  dt <- merge(w, pp, by = c("sex", "age", "year"), all.x = TRUE)
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

gbd_raw <- as.data.table(readRDS(paste0(wd_data, "temp_baseline_rates_gbd.rds")))
pop_obs <- as.data.table(readRDS(paste0(wd_data, "pop_observed_1990_2024.rds")))

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

## defensive clamps (as the parent applied before calibrating)
b_rates[CF >= 1, CF := 0.99]; b_rates[IR >= 1, IR := 0.99]
b_rates[CF < 0,  CF := 0];    b_rates[IR < 0,  IR := 0]

tps_input_cols <- copy(names(b_rates))
tps_input_nrow <- nrow(b_rates)

cat(sprintf("Built RHD baseline rates: %d rows | ages %d-%d | years %d-%d | sexes %s\n",
            nrow(b_rates), min(b_rates$age), max(b_rates$age),
            min(b_rates$year), max(b_rates$year),
            paste(unique(b_rates$sex), collapse = ", ")))
cat(sprintf("Calibration bands (%d): %s | incidence mode = %s | optimiser = %s\n",
            n_band, paste(band_labels, collapse = ", "),
            INCIDENCE_CALIBRATION_MODE, CALIB_OPTIMIZER))

#===============================================================================
# 3. GBD CALIBRATION TARGETS  (Number: RHD incidence, prevalence, RHD deaths;
#    all-cause deaths kept for VALIDATION only). Optional inverse-variance
#    weights from GBD 95% uncertainty intervals (upper/lower), if present.
#===============================================================================
gbd <- copy(gbd_raw)
setnames(gbd,
         c("sex_name", "age_name", "cause_name", "measure_name", "metric_name", "location_name"),
         c("sex",      "age",      "cause",      "measure",      "metric",      "location"))
gbd <- gbd[metric == "Number" & location %in% locs &
             year >= CAL_YEAR_START & year <= CAL_YEAR_END]

age_groups_keep <- unique(age_to_gbd_group(AGE_LO:AGE_HI))
has_ui <- all(c("upper", "lower") %in% names(gbd))

## one tidy target table per (measure) with value + optional UI-based log-variance.
grab_target <- function(cause_nm, measure_nm, valname) {
  cols <- c("location", "sex", "age", "year", "val")
  if (has_ui) cols <- c(cols, "upper", "lower")
  d <- gbd[cause == cause_nm & measure == measure_nm & age %in% age_groups_keep, ..cols]
  setnames(d, "val", valname)
  if (has_ui) {
    # var(log X) ~ ((log upper - log lower)/(2*1.96))^2 ; unit weight where UI missing/degenerate
    d[, logvar := ((log(pmax(upper, EPS_LOG)) - log(pmax(lower, EPS_LOG))) / (2 * 1.959964))^2]
    d[!is.finite(logvar) | logvar <= 0, logvar := NA_real_]
    d[, c("upper", "lower") := NULL]
    setnames(d, "logvar", paste0(valname, "_logvar"))
  }
  d
}
inc_t  <- grab_target(RHD_CAUSE, "Incidence",  "gbdInc")
prev_t <- grab_target(RHD_CAUSE, "Prevalence", "gbdPrev")
dth_t  <- grab_target(RHD_CAUSE, "Deaths",     "gbdDeath")
alld_t <- grab_target(ALL_CAUSE, "Deaths",     "gbdAllDeath")   # validation only

keyc <- c("location", "sex", "age", "year")
targets <- Reduce(function(x, y) merge(x, y, by = keyc, all = TRUE),
                  list(inc_t, prev_t, dth_t, alld_t))
for (cc in c("gbdInc", "gbdPrev", "gbdDeath", "gbdAllDeath"))
  targets[is.na(get(cc)), (cc) := 0]
if (nrow(targets) == 0)
  stop("No GBD calibration targets built -- check ages / years / measures.", call. = FALSE)

#===============================================================================
# 4. HAZARD HELPERS + MATRIX PROJECTOR  (aggregate well-sick-dead proxy)
#===============================================================================
# enforce probability / competing-risk constraints IN PLACE (adds no permanent
# columns). PRESERVES BG.mx and caps IR/CF into the remaining headroom so
# IR + BG.mx <= 1 and CF + BG.mx <= 1 without renormalising background.
enforce_tp_constraints <- function(dt, tp_eps = TP_EPS) {
  dt[is.na(IR), IR := 0]; dt[is.na(CF), CF := 0]; dt[is.na(BG.mx), BG.mx := 0]
  dt[IR < 0, IR := 0]; dt[IR > 1, IR := 1]
  dt[CF < 0, CF := 0]; dt[CF > 1, CF := 1]
  dt[BG.mx < 0, BG.mx := 0]
  dt[, headroom := pmax(1 - BG.mx - tp_eps, 0)]
  dt[IR > headroom, IR := headroom]
  dt[CF > headroom, CF := headroom]
  dt[, headroom := NULL]
  dt[]
}

# competing-risk split of two annual probabilities on the hazard scale.
# returns p_event1, p_event2, p_survive (sum == 1 exactly; zero-hazard safe).
compete2 <- function(p1, p2) {
  p1 <- pmin(pmax(p1, 0), 0.999999); p2 <- pmin(pmax(p2, 0), 0.999999)
  h1 <- -log1p(-p1); h2 <- -log1p(-p2); ht <- h1 + h2
  pany <- 1 - exp(-ht)
  share1 <- ifelse(ht > 0, h1 / ht, 0)
  list(e1 = pany * share1, e2 = pany * (1 - share1), surv = exp(-ht))
}

# Build, for one combo, the [n_age x n_yr] input matrices + the target arrays.
combo_data <- function(loc, sx) {
  ages  <- AGE_LO:AGE_HI; n_age <- length(ages)
  yrs   <- CAL_YEAR_START:CAL_YEAR_END; n_yr <- length(yrs)
  cr    <- b_rates[location == loc & sex == sx]
  to_mat <- function(col) {
    m <- matrix(0, n_age, n_yr, dimnames = list(ages, as.character(yrs)))
    idx <- cbind(match(cr$age, ages), match(cr$year, yrs))
    m[idx] <- cr[[col]]; m
  }
  IR0 <- to_mat("IR"); PREV0 <- to_mat("PREVt0"); RHDMX0 <- to_mat("DIS.mx.t0")
  BG0 <- to_mat("BG.mx"); NX <- to_mat("Nx")
  CFBASE <- ifelse(PREV0 > 0, RHDMX0 / PREV0, 0)

  # age -> group integer (1..G) in GBD-group order; group labels present.
  grp_lab_age <- age_to_gbd_group(ages)
  grp_present <- unique(grp_lab_age)                 # in ascending-age order
  grp_int     <- match(grp_lab_age, grp_present)     # 1..G per age
  G <- length(grp_present)

  # target matrices [G x n_yr] aligned to (grp_present x yrs), for this sex.
  tg <- targets[location == loc & sex == sx]
  tmat <- function(col) {
    m <- matrix(0, G, n_yr, dimnames = list(grp_present, as.character(yrs)))
    ri <- match(tg$age, grp_present); ci <- match(tg$year, yrs)
    ok <- !is.na(ri) & !is.na(ci)
    m[cbind(ri[ok], ci[ok])] <- tg[[col]][ok]; m
  }
  TGT <- list(inc = tmat("gbdInc"), prev = tmat("gbdPrev"),
              death = tmat("gbdDeath"), alldeath = tmat("gbdAllDeath"))
  # inverse-variance weights [G x n_yr] per measure (normalized to mean 1), else unit.
  wmat <- function(varcol, tgtmat) {
    W <- matrix(1, G, n_yr)
    if (has_ui && USE_IV_WEIGHTS && varcol %in% names(tg)) {
      m <- matrix(NA_real_, G, n_yr)
      ri <- match(tg$age, grp_present); ci <- match(tg$year, yrs)
      ok <- !is.na(ri) & !is.na(ci)
      m[cbind(ri[ok], ci[ok])] <- tg[[varcol]][ok]
      w <- 1 / m; w[!is.finite(w)] <- NA_real_
      if (any(is.finite(w))) { w[!is.finite(w)] <- mean(w[is.finite(w)]); W <- w / mean(w) }
    }
    W
  }
  WT <- list(inc = wmat("gbdInc_logvar", TGT$inc), prev = wmat("gbdPrev_logvar", TGT$prev),
             death = wmat("gbdDeath_logvar", TGT$death))
  # loss mask: keep cells with GBD Number >= TARGET_MIN_COUNT.
  MSK <- list(inc = TGT$inc >= TARGET_MIN_COUNT, prev = TGT$prev >= TARGET_MIN_COUNT,
              death = TGT$death >= TARGET_MIN_COUNT)

  list(ages = ages, yrs = yrs, n_age = n_age, n_yr = n_yr,
       IR0 = IR0, PREV0 = PREV0, CFBASE = CFBASE, BG0 = BG0, NX = NX,
       band = band_of(ages), grp_int = grp_int, grp_present = grp_present, G = G,
       TGT = TGT, WT = WT, MSK = MSK)
}

# Aggregate cohort projection with hazard competing risks, given band multipliers.
# alpha_band/beta_band are length-n_band multiplicative factors on IR / CF.
# Returns model [G x n_yr] matrices for incidence, prevalence, RHD deaths, all deaths.
project_aggregate <- function(cd, alpha_band, beta_band) {
  n_age <- cd$n_age; n_yr <- cd$n_yr
  aM <- alpha_band[cd$band]; bM <- beta_band[cd$band]     # per-age multipliers
  # per (age,year) annual probabilities
  p_inc_raw <- pmin(cd$IR0 * aM / pmax(1 - cd$PREV0, 0.5), 0.99)   # at-risk incidence prob
  cf_raw    <- pmin(cd$CFBASE * bM, 0.99)
  bg        <- pmin(cd$BG0, 0.99)
  cw <- compete2(p_inc_raw, bg)          # well: incidence vs background
  cs <- compete2(cf_raw,    bg)          # sick: RHD death vs background
  p_newsick <- cw$e1                     # [n_age x n_yr]
  p_rhddeath <- cs$e1
  surv_sick  <- cs$surv

  S <- matrix(0, n_age, n_yr); newsick <- matrix(0, n_age, n_yr)
  S[, 1] <- pmin(pmax(cd$NX[, 1] * cd$PREV0[, 1], 0), cd$NX[, 1])
  W1 <- cd$NX[, 1] - S[, 1]; newsick[, 1] <- W1 * p_newsick[, 1]
  for (t in 2:n_yr) {
    prev_surv <- S[, t - 1] * surv_sick[, t - 1] + newsick[, t - 1]   # cohort ages a-1 -> a
    S[2:n_age, t] <- prev_surv[1:(n_age - 1)]
    S[n_age, t]   <- S[n_age, t] + prev_surv[n_age]                   # open terminal group retains
    S[1, t]       <- cd$NX[1, t] * cd$PREV0[1, t]                     # age-0 births reseeded (GBD)
    S[, t] <- pmin(pmax(S[, t], 0), cd$NX[, t])
    Wt <- cd$NX[, t] - S[, t]
    newsick[, t] <- Wt * p_newsick[, t]
  }
  rhd_deaths <- S * p_rhddeath
  all_deaths <- rhd_deaths + cd$NX * bg          # RHD + background of whole pop (validation)
  agg <- function(M) rowsum(M, cd$grp_int)       # [G x n_yr]
  list(inc = agg(newsick), prev = agg(S), death = agg(rhd_deaths),
       alldeath = agg(all_deaths))
}

# normalized weighted log-scale squared error on one measure (masked cells only).
loss_term <- function(model, tgt, wt, msk) {
  if (!any(msk)) return(0)
  r <- (log(model[msk] + EPS_LOG) - log(tgt[msk] + EPS_LOG))^2 * wt[msk]
  sum(r) / sum(wt[msk])
}
# 2nd-difference smoothness penalty over ordered bands (0 if < 3 bands).
smooth_pen <- function(v) if (length(v) < 3) 0 else sum(diff(v, differences = 2)^2)

# Full penalized objective from a packed parameter vector eta||zeta (log-mults).
# `which` selects sub-objective for the sequential steps: "mortality","incidence","joint".
make_objective <- function(cd, which = "joint") {
  function(par) {
    eta  <- par[seq_len(n_band)]
    zeta <- par[n_band + seq_len(n_band)]
    pr <- project_aggregate(cd, exp(eta), exp(zeta))
    L_I <- loss_term(pr$inc,   cd$TGT$inc,   cd$WT$inc,   cd$MSK$inc)
    L_P <- loss_term(pr$prev,  cd$TGT$prev,  cd$WT$prev,  cd$MSK$prev)
    L_D <- loss_term(pr$death, cd$TGT$death, cd$WT$death, cd$MSK$death)
    pri_a <- sum((eta  / SIGMA_ALPHA_EFF)^2)
    pri_b <- sum((zeta / SIGMA_BETA)^2)
    smo   <- smooth_pen(eta) + smooth_pen(zeta)
    data_term <- switch(which,
      mortality = W_DEATH * L_D,
      incidence = W_INC * L_I + W_PREV * L_P,
      joint     = W_INC * L_I + W_PREV * L_P + W_DEATH * L_D)
    prior_term <- switch(which,
      mortality = LAMBDA_PRIOR * pri_b + LAMBDA_SMOOTH * smooth_pen(zeta),
      incidence = LAMBDA_PRIOR * pri_a + LAMBDA_SMOOTH * smooth_pen(eta),
      joint     = LAMBDA_PRIOR * (pri_a + pri_b) + LAMBDA_SMOOTH * smo)
    data_term + prior_term
  }
}

#===============================================================================
# 5. BOUNDED OPTIMISER  (multi-start L-BFGS-B default; nloptr / DEoptim optional)
#    Always evaluates the baseline (log-mults = 0) so the penalized objective is
#    never worse than baseline. Reports the multi-start spread.
#===============================================================================
optimise_block <- function(fn, npar, lower, upper, base_par, seed, n_starts = N_STARTS) {
  vals <- numeric(0); best_par <- base_par; best_val <- fn(base_par)
  base_val <- best_val
  set.seed(seed)
  # candidate start set: baseline + reproducible uniform draws within bounds
  starts <- c(list(base_par),
              lapply(seq_len(n_starts), function(i) runif(npar, lower, upper)))

  run_lbfgs <- function(p0) tryCatch(
    optim(p0, fn, method = "L-BFGS-B", lower = lower, upper = upper,
          control = list(factr = 1e7, pgtol = CONVERGE_TOL, maxit = 300)),
    error = function(e) NULL)

  if (CALIB_OPTIMIZER == "deoptim" && requireNamespace("DEoptim", quietly = TRUE)) {
    ctrl <- DEoptim::DEoptim.control(NP = 10 * npar, itermax = 120, trace = FALSE,
                                     reltol = CONVERGE_TOL, steptol = 30)
    de <- DEoptim::DEoptim(fn, lower, upper, control = ctrl)
    dp <- as.numeric(de$optim$bestmem)
    pol <- run_lbfgs(dp); if (!is.null(pol) && pol$value < fn(dp)) dp <- pol$par
    v <- fn(dp); vals <- c(vals, v); if (v < best_val) { best_val <- v; best_par <- dp }
  } else if (CALIB_OPTIMIZER == "nloptr" && requireNamespace("nloptr", quietly = TRUE)) {
    g <- tryCatch(nloptr::nloptr(base_par, fn, lb = lower, ub = upper,
             opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval = 4000, xtol_rel = 1e-6)),
             error = function(e) NULL)
    gp <- if (!is.null(g)) g$solution else base_par
    pol <- run_lbfgs(gp); if (!is.null(pol)) gp <- pol$par
    v <- fn(gp); vals <- c(vals, v); if (v < best_val) { best_val <- v; best_par <- gp }
  } else {
    for (p0 in starts) {
      r <- run_lbfgs(p0)
      if (!is.null(r)) { vals <- c(vals, r$value); if (r$value < best_val) { best_val <- r$value; best_par <- r$par } }
    }
  }
  list(par = best_par, value = best_val, base_value = base_val,
       start_values = vals,
       spread = if (length(vals)) c(min = min(vals), median = stats::median(vals), max = max(vals)) else
                c(min = base_val, median = base_val, max = base_val))
}

#===============================================================================
# 6. CALIBRATE EACH combo (location x sex): sequential -> joint
#===============================================================================
lower <- rep(log(MULT_LO), 2 * n_band)
upper <- rep(log(MULT_HI), 2 * n_band)
base0 <- rep(0, 2 * n_band)          # all multipliers = 1 (baseline)

combos <- unique(b_rates[, .(location, sex, cause)])
n_combos <- nrow(combos)
cat(sprintf("Layer-1 calibration: %d combo(s) | %d params/combo (%d bands x alpha,beta)\n",
            n_combos, 2 * n_band, n_band))

calibrate_combo <- function(ci) {
  loc <- combos$location[ci]; sx <- combos$sex[ci]; cse <- combos$cause[ci]
  cd  <- combo_data(loc, sx)
  seed_ci <- SEED + ci * 10000L

  ## Step 1 — mortality only (beta), incidence held at anchor (alpha = 0 on log scale)
  fn_m <- function(zeta) make_objective(cd, "mortality")(c(rep(0, n_band), zeta))
  o1 <- optimise_block(fn_m, n_band, lower[seq_len(n_band)], upper[seq_len(n_band)],
                       rep(0, n_band), seed_ci + 1L)
  zeta1 <- o1$par

  ## Step 2 — incidence only (alpha), mortality held at step-1 beta
  if (INCIDENCE_CALIBRATION_MODE == "fixed") {
    eta2 <- rep(0, n_band); o2 <- list(value = NA_real_, base_value = NA_real_, spread = c(min=NA,median=NA,max=NA))
  } else {
    fn_i <- function(eta) make_objective(cd, "incidence")(c(eta, zeta1))
    o2 <- optimise_block(fn_i, n_band, lower[seq_len(n_band)], upper[seq_len(n_band)],
                         rep(0, n_band), seed_ci + 2L)
    eta2 <- o2$par
  }

  ## Step 3 — joint refinement from the sequential solution (full penalized objective)
  fn_j <- make_objective(cd, "joint")
  start_joint <- c(eta2, zeta1)
  oj <- optimise_block(fn_j, 2 * n_band, lower, upper, base0, seed_ci + 3L)
  # ensure joint is seeded from the sequential solution too (compare, keep best)
  vj_seq <- fn_j(start_joint)
  if (is.finite(vj_seq) && vj_seq < oj$value) { oj$par <- start_joint; oj$value <- vj_seq }
  if (INCIDENCE_CALIBRATION_MODE == "fixed") oj$par[seq_len(n_band)] <- 0   # pin alpha

  eta <- oj$par[seq_len(n_band)]; zeta <- oj$par[n_band + seq_len(n_band)]
  alpha <- exp(eta); beta <- exp(zeta)

  ## baseline vs calibrated fit components (for diagnostics + acceptance)
  comp <- function(par) {
    e <- par[seq_len(n_band)]; z <- par[n_band + seq_len(n_band)]
    pr <- project_aggregate(cd, exp(e), exp(z))
    list(L_I = loss_term(pr$inc, cd$TGT$inc, cd$WT$inc, cd$MSK$inc),
         L_P = loss_term(pr$prev, cd$TGT$prev, cd$WT$prev, cd$MSK$prev),
         L_D = loss_term(pr$death, cd$TGT$death, cd$WT$death, cd$MSK$death),
         pr = pr)
  }
  cb <- comp(base0); cc <- comp(oj$par)
  pen <- function(cl, par) {
    e <- par[seq_len(n_band)]; z <- par[n_band + seq_len(n_band)]
    W_INC*cl$L_I + W_PREV*cl$L_P + W_DEATH*cl$L_D +
      LAMBDA_PRIOR*(sum((e/SIGMA_ALPHA_EFF)^2) + sum((z/SIGMA_BETA)^2)) +
      LAMBDA_SMOOTH*(smooth_pen(e) + smooth_pen(z))
  }
  pen_base <- pen(cb, base0); pen_cal <- pen(cc, oj$par)

  ## calibrated IR / CF written back to the tp rows (IR = IR_GBD x alpha; CF = CF_GBD x beta)
  cr <- copy(b_rates[location == loc & sex == sx])
  cr[, bnd := band_of(age)]
  cr[, IR := IR * alpha[bnd]]
  cr[, CF := CF * beta[bnd]]
  cr[, bnd := NULL]
  cr <- enforce_tp_constraints(cr)

  ## fit-by-group-year long table (model vs GBD, incl. all-cause validation)
  yrs <- cd$yrs; grp <- cd$grp_present
  mk_long <- function(mat, meas, src) data.table(
    location = loc, sex = sx, age = rep(grp, times = length(yrs)),
    year = rep(yrs, each = length(grp)),
    measure = meas, source = src, value = as.vector(mat))
  fit_long <- rbindlist(list(
    mk_long(cc$pr$inc,      "Incidence",  "model"),
    mk_long(cd$TGT$inc,     "Incidence",  "gbd"),
    mk_long(cc$pr$prev,     "Prevalence", "model"),
    mk_long(cd$TGT$prev,    "Prevalence", "gbd"),
    mk_long(cc$pr$death,    "Deaths",     "model"),
    mk_long(cd$TGT$death,   "Deaths",     "gbd"),
    mk_long(cc$pr$alldeath, "AllDeaths",  "model"),
    mk_long(cd$TGT$alldeath,"AllDeaths",  "gbd")))

  ## (mass-balance is built once at consolidation from pr_cal + TGT; see build_massbalance)

  factors <- data.table(location = loc, sex = sx, cause = cse,
                        band = band_labels, band_lo = band_edges,
                        alpha = alpha, beta = beta)

  diag <- data.table(location = loc, sex = sx, cause = cse,
    L_I_base = cb$L_I, L_P_base = cb$L_P, L_D_base = cb$L_D,
    L_I_cal = cc$L_I, L_P_cal = cc$L_P, L_D_cal = cc$L_D,
    penalized_base = pen_base, penalized_cal = pen_cal,
    step1_mortality_value = o1$value, step2_incidence_value = o2$value,
    step3_joint_value = oj$value,
    joint_spread_min = oj$spread[["min"]], joint_spread_med = oj$spread[["median"]],
    joint_spread_max = oj$spread[["max"]])

  list(rows = cr, factors = factors, diag = diag, fit_long = fit_long,
       cd = cd, alpha = alpha, beta = beta, pr_cal = cc$pr, pr_base = cb$pr)
}

results <- lapply(seq_len(n_combos), function(ci) {
  res <- calibrate_combo(ci)
  cat(sprintf("  [%d/%d] %s | %s : penalized %.4g -> %.4g | alpha[%s] beta[%s]\n",
              ci, n_combos, combos$location[ci], combos$sex[ci],
              res$diag$penalized_base, res$diag$penalized_cal,
              paste(sprintf("%.2f", res$alpha), collapse = ","),
              paste(sprintf("%.2f", res$beta),  collapse = ",")))
  res
})

#===============================================================================
# 7. CONSOLIDATE + MASS-BALANCE + BAKE INTO $tp
#===============================================================================
calibrated  <- rbindlist(lapply(results, `[[`, "rows"),     use.names = TRUE, fill = TRUE)
factors_out <- rbindlist(lapply(results, `[[`, "factors"),  use.names = TRUE, fill = TRUE)
diag_out    <- rbindlist(lapply(results, `[[`, "diag"),     use.names = TRUE, fill = TRUE)
fit_out     <- rbindlist(lapply(results, `[[`, "fit_long"), use.names = TRUE, fill = TRUE)
setcolorder(calibrated, intersect(tps_input_cols, names(calibrated)))

## mass-balance table per combo x band x year: model prevalence vs incidence x dwell,
## and the same for GBD. dwell (mean years in the RHD state) = prevalence / incidence.
build_massbalance <- function() {
  ages <- AGE_LO:AGE_HI
  out <- list()
  for (r in results) {
    cd <- r$cd; grp <- cd$grp_present; yrs <- cd$yrs
    # map each GBD group to a band via its lowest age
    grp_lo_age <- vapply(grp, function(g) min(ages[age_to_gbd_group(ages) == g]), numeric(1))
    grp_band   <- band_of(grp_lo_age)
    for (bi in seq_len(n_band)) {
      gi <- which(grp_band == bi); if (!length(gi)) next
      m_prev <- colSums(r$pr_cal$prev[gi, , drop = FALSE])
      m_inc  <- colSums(r$pr_cal$inc [gi, , drop = FALSE])
      g_prev <- colSums(cd$TGT$prev [gi, , drop = FALSE])
      g_inc  <- colSums(cd$TGT$inc  [gi, , drop = FALSE])
      out[[length(out) + 1L]] <- data.table(
        location = r$factors$location[1], sex = r$factors$sex[1],
        band = band_labels[bi], year = yrs,
        model_prev = m_prev, model_inc = m_inc,
        model_dwell = ifelse(m_inc > 0, m_prev / m_inc, NA_real_),
        gbd_prev = g_prev, gbd_inc = g_inc,
        gbd_dwell = ifelse(g_inc > 0, g_prev / g_inc, NA_real_))
    }
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}
massbalance <- build_massbalance()

## write human-readable diagnostics
fwrite(factors_out, paste0(wd_data, "calibration_layer1_factors.csv"))
fwrite(diag_out,    paste0(wd_data, "calibration_layer1_diagnostics.csv"))
fwrite(fit_out,     paste0(wd_data, "calibration_layer1_fit_by_group_year.csv"))
fwrite(massbalance, paste0(wd_data, "calibration_layer1_massbalance.csv"))
# keep the legacy filenames as aliases (nothing reads them, but avoid surprising a grep)
fwrite(factors_out, paste0(wd_data, "calibration_factors_random_tp.csv"))
fwrite(diag_out,    paste0(wd_data, "calibration_diagnostics_random_tp.csv"))

#===============================================================================
# 8. LAYER-2 INTERFACE  (unchanged; Stage-2 04b fills it against local echo targets)
#===============================================================================
calibration_targets_stage_template <- data.table(
  location = character(), year = integer(), sex = character(),
  age_lo = integer(), age_hi = integer(), stage = character(),
  target_prevalence = numeric(), target_type = character(), weight = numeric()
)
fwrite(calibration_targets_stage_template, STAGE_TMPL)

stage_params_uncalibrated <- if (file.exists(DISEASE_INPUTS_FILE)) {
  dmi_ <- readRDS(DISEASE_INPUTS_FILE)
  list(transitions = dmi_$transitions, p_rhd_death = dmi_$p_rhd_death,
       stage_split = dmi_$stage_split, rhd_d_fraction = dmi_$meta$rhd_d_fraction)
} else NULL

stage_target_file <- STAGE_TARGET_CANDIDATES[file.exists(STAGE_TARGET_CANDIDATES)][1]
stage_targets_present <- !is.na(stage_target_file)
stage_status <- if (stage_targets_present) "targets_present_not_yet_optimised" else
                "pending_local_echo_targets"
if (!stage_targets_present) {
  cat(strrep("!", 70), "\n", sep = "")
  cat("LAYER 2 (A/B/C/D STAGE CALIBRATION) IS PENDING LOCAL ECHO TARGETS.\n")
  cat("  No local stage-prevalence target file at either candidate path.\n")
  cat("  Stage transitions/mortality are PARTIALLY IDENTIFIED: Stage-2 (04b) uses\n")
  cat("  clinical priors/bounds + the Layer-1 RHD-mortality anchor, and saves the\n")
  cat("  near-optimal set rather than presenting one matrix as identified.\n")
  cat(strrep("!", 70), "\n", sep = "")
}

stage_calibration <- list(
  status = stage_status,
  targets_file = if (stage_targets_present) stage_target_file else NA_character_,
  targets_template = calibration_targets_stage_template,
  template_file = STAGE_TMPL,
  calibratable_params = c("p_A_to_B","p_B_to_C","p_C_to_D","p_A_to_no_rhd",
                          "p_B_to_A","p_C_to_B","rhd_d_fraction"),
  stage_params = stage_params_uncalibrated,
  loss_note = "Weighted squared-log-error on stage prevalence; Stage-2 in 04b."
)

#===============================================================================
# 9. ASSEMBLE + WRITE THE CALIBRATED-PARAMETER BUNDLE  (schema preserved + $layer1)
#===============================================================================
layer1 <- list(
  incidence_parameters = factors_out[, .(location, sex, cause, band, band_lo, alpha)],
  mortality_parameters = factors_out[, .(location, sex, cause, band, band_lo, beta)],
  objective_components = diag_out,
  optimizer_diagnostics = list(
    optimizer = CALIB_OPTIMIZER, n_starts = N_STARTS,
    incidence_mode = INCIDENCE_CALIBRATION_MODE,
    bands = band_labels, band_edges = band_edges,
    weights = c(W_INC = W_INC, W_PREV = W_PREV, W_DEATH = W_DEATH),
    priors = c(sigma_alpha = SIGMA_ALPHA_EFF, sigma_beta = SIGMA_BETA,
               lambda_prior = LAMBDA_PRIOR, lambda_smooth = LAMBDA_SMOOTH),
    used_iv_weights = has_ui && USE_IV_WEIGHTS),
  fit_by_group_year = fit_out,
  mass_balance = massbalance,
  validation = list(
    penalized_improved = all(diag_out$penalized_cal <= diag_out$penalized_base + 1e-8),
    total_penalized_base = sum(diag_out$penalized_base),
    total_penalized_cal  = sum(diag_out$penalized_cal))
)

calibrated_rhd_parameters <- list(
  tp = calibrated, factors = factors_out, diagnostics = diag_out,
  stage_calibration = stage_calibration, layer1 = layer1,
  meta = list(
    location = locs, calib_year_start = CAL_YEAR_START, calib_year_end = CAL_YEAR_END,
    calib_last_year = max(calibrated$year),
    granularity = "band", bands = band_labels, band_edges = band_edges,
    incidence_mode = INCIDENCE_CALIBRATION_MODE, optimizer = CALIB_OPTIMIZER,
    tp_schema = tps_input_cols,
    built_from = c("temp_baseline_rates_gbd.rds", "pop_observed_1990_2024.rds",
                   basename(DISEASE_INPUTS_FILE)))
)
saveRDS(calibrated_rhd_parameters, file = OUT_BUNDLE)

old_chunks <- list.files(wd_data, pattern = "^adjusted_searo_part[0-9]+\\.rds$", full.names = TRUE)
if (length(old_chunks)) { file.remove(old_chunks)
  cat(sprintf("Removed %d obsolete adjusted_searo_part*.rds chunk file(s).\n", length(old_chunks))) }

#===============================================================================
# 10. VALIDATION  (hard failures; identify the offending combo/age on violation)
#===============================================================================
cat("\n", strrep("=", 70), "\nLAYER-1 VALIDATION\n", strrep("=", 70), "\n", sep = "")
report_bad <- function(dt, bad_idx, msg) {
  bad_idx[is.na(bad_idx)] <- TRUE
  if (any(bad_idx)) {
    cat("OFFENDING ROWS (", msg, "):\n", sep = "")
    print(utils::head(dt[bad_idx, .(location, sex, cause, age, year, IR, CF, BG.mx)], 10))
    stop(msg, " -- ", sum(bad_idx), " offending row(s); see above.", call. = FALSE)
  }
}
report_bad(calibrated, calibrated[, is.na(IR)],             "IR is NA")
report_bad(calibrated, calibrated[, is.na(CF)],             "CF is NA")
report_bad(calibrated, calibrated[, is.na(BG.mx)],          "BG.mx is NA")
report_bad(calibrated, calibrated[, IR < 0 | IR > 1],       "IR outside [0,1]")
report_bad(calibrated, calibrated[, CF < 0 | CF > 1],       "CF outside [0,1]")
report_bad(calibrated, calibrated[, BG.mx < 0 | BG.mx > 1], "BG.mx outside [0,1]")
report_bad(calibrated, calibrated[, IR + BG.mx > 1 + 1e-9], "IR + BG.mx > 1 (competing risk)")
report_bad(calibrated, calibrated[, CF + BG.mx > 1 + 1e-9], "CF + BG.mx > 1 (competing risk)")
stopifnot("row count != input" = nrow(calibrated) == tps_input_nrow,
          "schema != input"    = setequal(names(calibrated), tps_input_cols))

if (!layer1$validation$penalized_improved)
  stop("Calibrated penalized objective exceeds baseline for some combo -- ",
       "the baseline start-0 guarantee failed; investigate.", call. = FALSE)

cat("All probability/row constraints satisfied (IR, CF, BG.mx in [0,1]; sums <= 1).\n")
cat(sprintf("Rows: %d (matches input). Schema matches input: TRUE.\n", nrow(calibrated)))
cat(sprintf("\nPenalized objective (sum over combos): baseline = %.4g -> calibrated = %.4g (%.1f%% lower)\n",
            layer1$validation$total_penalized_base, layer1$validation$total_penalized_cal,
            100 * (layer1$validation$total_penalized_base - layer1$validation$total_penalized_cal) /
              max(layer1$validation$total_penalized_base, 1e-9)))
cat("\nPer-combo loss (baseline -> calibrated):\n")
print(diag_out[, .(location, sex,
                   L_I = sprintf("%.3g->%.3g", L_I_base, L_I_cal),
                   L_P = sprintf("%.3g->%.3g", L_P_base, L_P_cal),
                   L_D = sprintf("%.3g->%.3g", L_D_base, L_D_cal),
                   joint_spread = sprintf("%.3g/%.3g/%.3g", joint_spread_min, joint_spread_med, joint_spread_max))])

## mass-balance headline (last calibration year): model vs GBD implied dwell by band
mb_last <- massbalance[year == max(year)]
cat(sprintf("\nMass balance (prevalence ~ incidence x dwell), %d, dwell in years:\n", max(massbalance$year)))
print(mb_last[, .(location, sex, band,
                  model_dwell = round(model_dwell, 1), gbd_dwell = round(gbd_dwell, 1))])

## all-cause deaths retained for VALIDATION only — report model vs GBD envelope
val_all <- fit_out[measure == "AllDeaths", .(v = sum(value)), by = source]
cat(sprintf("\nAll-cause deaths (VALIDATION only; NOT in the loss), summed over window: model=%.0f | GBD=%.0f\n",
            val_all[source == "model", v], val_all[source == "gbd", v]))

cat("\nWrote:\n")
cat(sprintf("  %scalibrated_rhd_parameters.rds  ($tp/$factors/$diagnostics/$stage_calibration/$layer1/$meta)\n", wd_data))
cat(sprintf("  %scalibration_layer1_{factors,diagnostics,fit_by_group_year,massbalance}.csv\n", wd_data))
cat(sprintf("\nLayer-2 stage calibration status: %s\n", stage_status))
cat("Reminder: 05/06 consume the calibrated IR/CF AS-IS; stage params come from 03 until Stage-2 (04b).\n")
