# ==============================================================================
# RHD secondary-prevention investment case: INITIAL-STATE ASSEMBLER
# scripts/05_build_baseline.R
#
# Structure after: Coates et al., Lancet Glob Health 2021 (PMC9087136).
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT (rewritten in this refactor)
# ------------------------------------------------------------------------------
# The previous version was parent-NCD leftover (blood pressure / salt / HTN / COVID
# / multi-country) and referenced many inputs that do not exist in this repo. It
# has been fully replaced.
#
# This script assembles the INITIAL STATE and every input needed to RUN the
# secondary-prevention model in 06 WITHOUT re-running 01-04. It reads only the
# PERSISTED outputs of the upstream refactored scripts:
#
#   INPUTS
#     data/pop_projection_2025_2100.rds  (02) population [age x sex x year], horizon
#     data/disease_model_inputs.rds      (03) rate arrays + clinical/effect params +
#                                             severity split + coverage-ramp structure
#     data/adjusted_searo_part{1..10}.rds(04) CALIBRATED transition probabilities
#                                             (IR incidence, CF case-fatality), 2000-2019
#
#   OUTPUT (written to wd_data = data/):
#     baseline_state.rds  — a single self-contained "baseline state" object that 06
#       loads directly. It is a per-location bundle (location-general; Indonesia only
#       for now) holding, per location:
#         pop      [age x sex x year]  population over the horizon (from 02)
#         ir       [age x sex x year]  incidence probability — CALIBRATED (04, last
#                                      calibration year) with the secular trend applied
#         cf       [age x sex]         CALIBRATED RHD case-fatality (04) — base-year
#                                      RHD-mortality ANCHOR (see NOTE below)
#         oth_mort [age x sex]         background (non-RHD) competing mortality (03)
#         seed     list(mild,severe,post) [age x sex] seeded prevalent pool at year 1
#         clinical / effects           tagged clinical + intervention parameters (03)
#         coverage                     REALISED coverage trajectories per scenario
#                                      (ref, sap) over the horizon + screening window
#     plus a `meta` block (AGES, SEXES, years, base_year, scenarios, ramp window, ...).
#
# NOTE on the calibrated TPs (documented modelling decision):
#   * CALIBRATED IR drives the incidence inflow of new (mild) RHD directly.
#   * CALIBRATED CF is carried in as the base-year RHD case-fatality ANCHOR. The
#     06 engine is the Coates mild->severe->post tunnel model whose within-sick RHD
#     mortality comes from the clinical parameters; 06 checks that its base-year
#     aggregate RHD death rate is within an order-of-magnitude band of CF x prevalence.
#   * Background mortality is held FIXED (from 03's GBD all-cause minus RHD). This is
#     consistent with 04, which calibrated only IR and CF and held background fixed.
#
# NO monetary / cost values live here — economics are in 08 only.
# ==============================================================================

library(data.table)
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

# ------------------------------------------------------------------------------
# 0. CONFIG (honour globals from 00_run_all.R; else standalone defaults)
# ------------------------------------------------------------------------------
LOCATIONS <- if (exists("LOCATIONS")) LOCATIONS else "Indonesia"
SCENARIOS <- c("ref", "sap")           # reference vs secondary-prevention scale-up

IN_POP     <- paste0(wd_data, "pop_projection_2025_2100.rds")
IN_DISEASE <- paste0(wd_data, "disease_model_inputs.rds")
OUT_FILE   <- paste0(wd_data, "baseline_state.rds")

for (f in c(IN_POP, IN_DISEASE))
  if (!file.exists(f))
    stop("Missing required input:\n  ", f,
         "\n  Run the upstream script (02 / 03) first.", call. = FALSE)

# ------------------------------------------------------------------------------
# 1. LOAD PERSISTED UPSTREAM OUTPUTS
# ------------------------------------------------------------------------------
message("── 05_build_baseline.R : assembling initial state ─────────")

dmi      <- readRDS(IN_DISEASE)            # 03 disease-model inputs
pop_long <- as.data.table(readRDS(IN_POP)) # 02 population (tidy long)

AGES  <- dmi$meta$AGES
SEXES <- dmi$meta$SEXES
years <- dmi$meta$years                    # projection horizon (2025..2100)
n_age <- length(AGES); n_sex <- length(SEXES); n_years <- length(years)
base_year <- min(years)

# calibrated TPs (04) — read + row-bind all chunks
tp_files <- list.files(wd_data, pattern = "^adjusted_searo_part[0-9]+\\.rds$",
                       full.names = TRUE)
if (length(tp_files) == 0)
  stop("No calibrated TP files (adjusted_searo_part*.rds) found in ", wd_data,
       ".\n  Run 04_calibration_random_tp.R first.", call. = FALSE)
calib <- rbindlist(lapply(tp_files, readRDS), use.names = TRUE, fill = TRUE)
calib_last_year <- max(calib$year)         # last calibration year (2019)

message(sprintf("  Loaded: pop(%d rows) | disease inputs(horizon %d-%d) | calibrated TPs(%d rows, %d-%d)",
                nrow(pop_long), min(years), max(years),
                nrow(calib), min(calib$year), max(calib$year)))

# ------------------------------------------------------------------------------
# 2. HELPERS  (tidy long -> [age x sex] matrix / [age x sex x year] array)
# ------------------------------------------------------------------------------
empty_mat <- function() matrix(0, n_age, n_sex, dimnames = list(AGES, SEXES))

# build an [age x sex] matrix for one location from a long DT with cols age,sex,<val>.
# ages beyond the calibrated top age (AGE_HI) inherit the top-age value (95+ open
# group); ages absent entirely stay 0.
mat_from_long <- function(dt, valcol) {
  m <- empty_mat()
  for (s in SEXES) {
    ds <- dt[sex == s]
    if (!nrow(ds)) next
    idx <- match(ds$age, AGES)
    keep <- !is.na(idx)
    m[idx[keep], s] <- ds[[valcol]][keep]
    # forward-fill ages above the max supplied age with that top value (95+ open)
    top_age <- max(ds$age)
    if (top_age < max(AGES)) {
      fill_rows <- which(AGES > top_age)
      m[fill_rows, s] <- m[as.character(top_age), s]
    }
  }
  m
}

# coverage ramp: baseline before window, linear across [start,end], target after
ramp_traj <- function(baseline, target, start, end) {
  frac <- (years - start) / (end - start)
  frac <- pmin(pmax(frac, 0), 1)
  baseline + (target - baseline) * frac
}

# ------------------------------------------------------------------------------
# 3. BUILD PER-LOCATION INITIAL STATE
# ------------------------------------------------------------------------------
cov <- dmi$coverage
build_location_state <- function(loc) {

  ## 3a. population array [age x sex x year] from 02 --------------------------
  pl <- pop_long[location == loc & age %in% AGES & sex %in% SEXES]
  if (!nrow(pl)) stop("No population rows for location '", loc, "'.", call. = FALSE)
  pop_arr <- array(0, dim = c(n_age, n_sex, n_years),
                   dimnames = list(AGES, SEXES, years))
  for (iy in seq_len(n_years)) {
    py <- pl[year == years[iy]]
    m  <- empty_mat()
    for (s in SEXES) {
      ps <- py[sex == s]
      m[match(ps$age, AGES), s] <- ps$Nx
    }
    pop_arr[, , iy] <- m
  }

  ## 3b. calibrated IR & CF at the last calibration year, [age x sex] --------
  cb <- calib[location == loc & year == calib_last_year]
  if (!nrow(cb))
    stop("No calibrated TP rows for location '", loc, "' at year ",
         calib_last_year, ".", call. = FALSE)
  ir_base <- mat_from_long(cb, "IR")   # incidence probability, calibrated
  cf_base <- mat_from_long(cb, "CF")   # RHD case-fatality, calibrated (anchor)

  ## 3c. incidence array over horizon: calibrated base pattern x secular trend
  ##     trend anchored at the last calibration year (continuous decline forward).
  trend <- dmi$meta$incidence_trend
  ir_arr <- array(0, dim = c(n_age, n_sex, n_years),
                  dimnames = list(AGES, SEXES, years))
  for (iy in seq_len(n_years))
    ir_arr[, , iy] <- ir_base * trend^(years[iy] - calib_last_year)

  ## 3d. background (non-RHD) mortality + prevalence seed from 03 -------------
  oth_mort <- dmi$rates_by_year$oth_mort[, , 1]   # base-year pattern (held)
  prev_seed <- dmi$rates_by_year$prev_seed        # [age x sex] prevalence fraction

  ## 3e. seed the prevalent pool at year 1 (split by severity) ---------------
  prev_pool <- prev_seed * pop_arr[, , 1]
  seed <- list(
    mild   = prev_pool * dmi$seed_split[["mild"]],
    severe = prev_pool * dmi$seed_split[["severe"]],
    post   = prev_pool * dmi$seed_split[["post"]]
  )

  ## 3f. realised coverage trajectories per scenario -------------------------
  #  Only SAP differs between arms; HF & surgery held at baseline in both.
  coverage <- list(
    ref = list(
      sap  = ramp_traj(cov$sap_ref_baseline, cov$sap_ref_target, cov$ramp_start, cov$ramp_end),
      hf   = ramp_traj(cov$hf_baseline,      cov$hf_target,      cov$ramp_start, cov$ramp_end),
      surg = ramp_traj(cov$surg_baseline,    cov$surg_target,    cov$ramp_start, cov$ramp_end)
    ),
    sap = list(
      sap  = ramp_traj(cov$sap_up_baseline,  cov$sap_up_target,  cov$ramp_start, cov$ramp_end),
      hf   = ramp_traj(cov$hf_baseline,      cov$hf_target,      cov$ramp_start, cov$ramp_end),
      surg = ramp_traj(cov$surg_baseline,    cov$surg_target,    cov$ramp_start, cov$ramp_end)
    ),
    screen_age_lo = cov$screen_age_lo,
    screen_age_hi = cov$screen_age_hi
  )

  list(pop = pop_arr, ir = ir_arr, cf = cf_base, oth_mort = oth_mort,
       seed = seed, clinical = dmi$clinical, effects = dmi$effects,
       coverage = coverage)
}

states <- setNames(lapply(LOCATIONS, build_location_state), LOCATIONS)

# ------------------------------------------------------------------------------
# 4. VALIDATION  (fail loudly BEFORE writing)
# ------------------------------------------------------------------------------
message("── Validation ─────────────────────────────────")

chk <- function(x, nm) {
  if (any(is.na(x)))  stop(nm, ": contains NA.", call. = FALSE)
  if (any(x < 0))     stop(nm, ": contains negative values.", call. = FALSE)
}
chk_prob <- function(x, nm) {
  chk(x, nm)
  if (any(x > 1)) stop(nm, ": contains values > 1 (must be a probability/rate).", call. = FALSE)
}

for (loc in LOCATIONS) {
  st <- states[[loc]]

  # shapes
  if (!identical(dim(st$pop), c(n_age, n_sex, n_years)))
    stop(loc, ": population array has wrong dimensions.", call. = FALSE)
  if (!identical(dim(st$ir), c(n_age, n_sex, n_years)))
    stop(loc, ": incidence array has wrong dimensions.", call. = FALSE)
  if (!identical(dim(st$cf), c(n_age, n_sex)))
    stop(loc, ": CF matrix has wrong dimensions.", call. = FALSE)

  # non-negativity / probability ranges
  chk(st$pop, paste0(loc, " pop"))
  chk_prob(st$ir, paste0(loc, " ir"))
  chk_prob(st$cf, paste0(loc, " cf"))
  chk_prob(st$oth_mort, paste0(loc, " oth_mort"))

  # competing risks at base year: IR + background <= 1 ; CF + background <= 1
  ir1 <- st$ir[, , 1]
  if (any(ir1 + st$oth_mort > 1 + 1e-9))
    stop(loc, ": IR + background mortality exceeds 1 at some age-sex cell.", call. = FALSE)
  if (any(st$cf + st$oth_mort > 1 + 1e-9))
    stop(loc, ": CF + background mortality exceeds 1 at some age-sex cell.", call. = FALSE)

  # seeded prevalent pool: non-negative and NOT exceeding the year-1 population
  sick1 <- st$seed$mild + st$seed$severe + st$seed$post
  chk(sick1, paste0(loc, " seeded sick pool"))
  if (any(sick1 > st$pop[, , 1] + 1e-6))
    stop(loc, ": seeded sick pool exceeds year-1 population somewhere.", call. = FALSE)

  # coverage trajectories in [0,1] for both arms
  for (sc in SCENARIOS) for (k in c("sap", "hf", "surg")) {
    v <- st$coverage[[sc]][[k]]
    if (length(v) != n_years) stop(loc, ": coverage ", sc, "/", k, " wrong length.", call. = FALSE)
    chk_prob(v, paste0(loc, " coverage ", sc, "/", k))
  }

  # order-of-magnitude anchor: seeded RHD prevalence vs a sane band
  prev_cnt <- sum(sick1)
  message(sprintf("  %s | seeded RHD prevalence (year %d) ~ %s | max sick/pop = %.4f%%",
                  loc, base_year, formatC(round(prev_cnt), format = "d", big.mark = ","),
                  100 * max(sick1 / pmax(st$pop[, , 1], 1))))
  if (prev_cnt < 1e5 || prev_cnt > 1e7)
    stop(loc, ": seeded RHD prevalence ", round(prev_cnt),
         " outside the sane band 1e5-1e7.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 5. ASSEMBLE + PERSIST THE BASELINE STATE
# ------------------------------------------------------------------------------
baseline_state <- list(
  locations = LOCATIONS,
  states    = states,
  meta = list(
    AGES = AGES, SEXES = SEXES, years = years, base_year = base_year,
    scenarios = SCENARIOS,
    ramp_start = cov$ramp_start, ramp_end = cov$ramp_end,
    incidence_trend = dmi$meta$incidence_trend,
    calib_last_year = calib_last_year,
    RATE_BASE_YEAR  = dmi$meta$RATE_BASE_YEAR,
    intervention_labels = c(ref = "none", sap = "echo_screening_plus_SAP"),
    built_from = c(basename(IN_POP), basename(IN_DISEASE),
                   "adjusted_searo_part*.rds")
  )
)

saveRDS(baseline_state, file = OUT_FILE)

message("── Saved ──────────────────────────────────────")
message(sprintf("  %s  | locations: %s | horizon %d-%d | scenarios: %s",
                basename(OUT_FILE), paste(LOCATIONS, collapse = ", "),
                min(years), max(years), paste(SCENARIOS, collapse = ", ")))
message("── 05_build_baseline.R complete ───────────────────────")
message("  Next: 06_run_prevention_model.R")
