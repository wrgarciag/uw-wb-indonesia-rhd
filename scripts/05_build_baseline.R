# ==============================================================================
# RHD secondary-prevention investment case: INITIAL-STATE ASSEMBLER (A/B/C/D)
# scripts/05_build_baseline.R
#
# Disease structure: WHF RHD stages A/B/C/D (No RHD -> A <-> B <-> C -> D ->
# RHD death, competing other-cause death from every stage). Surgery is a SERVICE
# (fraction of C/D treated each cycle), NOT a health state — there is no
# post-surgery stock here.
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT
# ------------------------------------------------------------------------------
# Assembles the INITIAL STATE and every input needed to RUN the model in 06
# WITHOUT re-running 01-04. It reads only the PERSISTED upstream outputs:
#
#   INPUTS
#     data/pop_projection_2025_2100.rds       (02) population [age x sex x year]
#     data/disease_model_inputs.rds           (03) rate arrays + A/B/C/D natural
#                                                  history + effects + surgery +
#                                                  stage_split + cascade structure
#     data/calibrated_rhd_parameters.rds       (04) CALIBRATED IR/CF ($tp), 2000-2019
#
#   OUTPUT (written to wd_data = data/):
#     baseline_state.rds — a single self-contained per-location bundle holding:
#       pop        [age x sex x year] population over the horizon (from 02)
#       ir         [age x sex x year] incidence probability — CALIBRATED (04) with
#                                     the secular trend applied (drives inflow to A)
#       cf         [age x sex]        CALIBRATED total-RHD case-fatality (04) — an
#                                     ANCHOR only (A/B/C/D deaths come from p_rhd_death)
#       oth_mort   [age x sex]        background (non-RHD) competing mortality (03)
#       seed       list(A,B,C,D) [age x sex] prevalent pool at year 1 (prev x split)
#       transitions / p_rhd_death / effects / surgery   A/B/C/D params (03)
#       coverage   REALISED cascade + surgery trajectories per scenario (ref, sap)
#     plus a `meta` block (AGES, SEXES, years, base_year, scenarios, ramp window).
#
# NOTE on the calibrated inputs (documented modelling decision):
#   * CALIBRATED IR drives the incidence inflow of new RHD into stage A directly.
#   * CALIBRATED CF is carried as a total-RHD-mortality ANCHOR only; the A/B/C/D
#     engine derives RHD deaths from the per-stage p_rhd_death probabilities.
#   * Background mortality is data-fed (GBD all-cause minus RHD) and held FIXED;
#     it is the SAME age x sex competing risk for every living stage.
#
# NO monetary / cost values live here — economics are in 08 only.
# ==============================================================================

library(data.table)
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

# getp(): honour a 00_run_all.R global when set, else the standalone default.
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
RHD_PREV_LO <- getp("RHD_PREV_LO", 1e5)   # seeded-prevalence sanity band (COUNTRY-settable;
RHD_PREV_HI <- getp("RHD_PREV_HI", 1e7)   #   wide default covers Indonesia and Uganda)

# ------------------------------------------------------------------------------
# 0. CONFIG (honour globals from 00_run_all.R; else standalone defaults)
# ------------------------------------------------------------------------------
LOCATIONS <- if (exists("LOCATIONS")) LOCATIONS else "Indonesia"
SCENARIOS <- c("ref", "sap")           # reference vs secondary-prevention scale-up

IN_POP     <- paste0(wd_data, "pop_projection_2025_2100.rds")
IN_DISEASE <- paste0(wd_data, "disease_model_inputs.rds")
IN_CALIB   <- paste0(wd_data, "calibrated_rhd_parameters.rds")
OUT_FILE   <- paste0(wd_data, "baseline_state.rds")

for (f in c(IN_POP, IN_DISEASE, IN_CALIB))
  if (!file.exists(f))
    stop("Missing required input:\n  ", f,
         "\n  Run the upstream script (02 / 03 / 04) first.", call. = FALSE)

# ------------------------------------------------------------------------------
# 1. LOAD PERSISTED UPSTREAM OUTPUTS
# ------------------------------------------------------------------------------
message("── 05_build_baseline.R : assembling A/B/C/D initial state ─────────")

dmi      <- readRDS(IN_DISEASE)             # 03 disease-model inputs
pop_long <- as.data.table(readRDS(IN_POP))  # 02 population (tidy long)
calib    <- readRDS(IN_CALIB)               # 04 calibrated-parameter bundle

AGES  <- dmi$meta$AGES
SEXES <- dmi$meta$SEXES
years <- dmi$meta$years                    # projection horizon (2025..2100)
n_age <- length(AGES); n_sex <- length(SEXES); n_years <- length(years)
base_year <- min(years)

# calibrated IR/CF table (Layer 1) — a single self-describing bundle now.
calib_tp <- as.data.table(calib$tp)
calib_last_year <- calib$meta$calib_last_year         # last calibration year (2019)

message(sprintf("  Loaded: pop(%d rows) | disease inputs(horizon %d-%d) | calibrated IR/CF(%d rows, %d-%d) | stage-calib %s",
                nrow(pop_long), min(years), max(years),
                nrow(calib_tp), min(calib_tp$year), max(calib_tp$year),
                calib$stage_calibration$status))

# ------------------------------------------------------------------------------
# 2. HELPERS  (tidy long -> [age x sex] matrix / coverage ramp)
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
cov  <- dmi$coverage
surg <- dmi$surgery

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
  cb <- calib_tp[location == loc & year == calib_last_year]
  if (!nrow(cb))
    stop("No calibrated IR/CF rows for location '", loc, "' at year ",
         calib_last_year, ".", call. = FALSE)
  ir_base <- mat_from_long(cb, "IR")   # incidence probability, calibrated
  cf_base <- mat_from_long(cb, "CF")   # total-RHD case-fatality, calibrated (anchor)

  ## 3c. incidence array over horizon: calibrated base pattern x secular trend
  trend <- dmi$meta$incidence_trend
  ir_arr <- array(0, dim = c(n_age, n_sex, n_years),
                  dimnames = list(AGES, SEXES, years))
  for (iy in seq_len(n_years))
    ir_arr[, , iy] <- ir_base * trend^(years[iy] - calib_last_year)

  ## 3d. background (non-RHD) mortality + prevalence seed from 03 -------------
  oth_mort  <- dmi$rates_by_year$oth_mort[, , 1]   # base-year pattern (held)
  prev_seed <- dmi$rates_by_year$prev_seed         # [age x sex] prevalence fraction

  ## 3e. seed the prevalent pool at year 1, split across A/B/C/D --------------
  prev_pool <- prev_seed * pop_arr[, , 1]
  ss <- dmi$stage_split
  seed <- list(
    A = prev_pool * ss[["A"]],
    B = prev_pool * ss[["B"]],
    C = prev_pool * ss[["C"]],
    D = prev_pool * ss[["D"]]
  )

  ## 3f. realised cascade + surgery coverage trajectories per scenario --------
  #  Reference holds baselines flat; scale-up ramps to the 2050 targets.
  #  Effective diagnosis/treatment are capped by the earlier cascade stages.
  cascade_for <- function(arm) {                          # arm = "ref" | "up"
    screen    <- ramp_traj(cov[[paste0("screen_",    arm, "_baseline")]],
                           cov[[paste0("screen_",    arm, "_target")]],   cov$ramp_start, cov$ramp_end)
    diagnosis <- ramp_traj(cov[[paste0("diagnosis_", arm, "_baseline")]],
                           cov[[paste0("diagnosis_", arm, "_target")]],   cov$ramp_start, cov$ramp_end)
    treatment <- ramp_traj(cov[[paste0("treatment_", arm, "_baseline")]],
                           cov[[paste0("treatment_", arm, "_target")]],   cov$ramp_start, cov$ramp_end)
    surgery   <- ramp_traj(surg[[paste0("coverage_", arm, "_baseline")]],
                           surg[[paste0("coverage_", arm, "_target")]],   cov$ramp_start, cov$ramp_end)
    # cumulative-coverage logic: effective downstream <= upstream cascade stages
    eff_diag <- pmin(screen, diagnosis)
    eff_trt  <- pmin(screen, diagnosis, treatment)
    list(screen = screen, diagnosis = diagnosis, treatment = treatment,
         eff_diagnosis = eff_diag, eff_treatment = eff_trt, surgery = surgery)
  }
  coverage <- list(
    ref = cascade_for("ref"),
    sap = cascade_for("up"),
    screen_age_restrict = isTRUE(cov$screen_age_restrict),
    screen_age_lo = cov$screen_age_lo,
    screen_age_hi = cov$screen_age_hi
  )

  list(pop = pop_arr, ir = ir_arr, cf = cf_base, oth_mort = oth_mort,
       seed = seed,
       transitions = dmi$transitions, p_rhd_death = dmi$p_rhd_death,
       effects = dmi$effects, surgery = surg,
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

  # competing risks at base year: IR + background <= 1
  ir1 <- st$ir[, , 1]
  if (any(ir1 + st$oth_mort > 1 + 1e-9))
    stop(loc, ": IR + background mortality exceeds 1 at some age-sex cell.", call. = FALSE)

  # seeded prevalent pool A+B+C+D: non-negative and NOT exceeding year-1 population
  sick1 <- st$seed$A + st$seed$B + st$seed$C + st$seed$D
  chk(sick1, paste0(loc, " seeded sick pool"))
  if (any(sick1 > st$pop[, , 1] + 1e-6))
    stop(loc, ": seeded sick pool exceeds year-1 population somewhere.", call. = FALSE)

  # coverage trajectories in [0,1]; effective <= stated for both arms
  for (sc in SCENARIOS) {
    cc <- st$coverage[[sc]]
    for (k in c("screen", "diagnosis", "treatment", "eff_diagnosis", "eff_treatment", "surgery")) {
      v <- cc[[k]]
      if (length(v) != n_years) stop(loc, ": coverage ", sc, "/", k, " wrong length.", call. = FALSE)
      chk_prob(v, paste0(loc, " coverage ", sc, "/", k))
    }
    if (any(cc$eff_treatment > cc$diagnosis + 1e-9) ||
        any(cc$eff_treatment > cc$screen + 1e-9) ||
        any(cc$eff_diagnosis > cc$screen + 1e-9))
      stop(loc, ": effective coverage exceeds an upstream cascade stage (", sc, ").", call. = FALSE)
  }
  # reference cascade must be constant over the horizon
  for (k in c("screen", "diagnosis", "treatment")) {
    v <- st$coverage$ref[[k]]
    if (diff(range(v)) > 1e-12)
      stop(loc, ": reference cascade '", k, "' is not constant over the horizon.", call. = FALSE)
  }

  # order-of-magnitude anchor: seeded RHD prevalence vs a sane band
  prev_cnt <- sum(sick1)
  message(sprintf("  %s | seeded RHD prevalence (year %d) ~ %s | A/B/C/D = %s | max sick/pop = %.4f%%",
                  loc, base_year, formatC(round(prev_cnt), format = "d", big.mark = ","),
                  paste(sprintf("%.0f%%", 100 * c(sum(st$seed$A), sum(st$seed$B),
                                                  sum(st$seed$C), sum(st$seed$D)) / prev_cnt),
                        collapse = "/"),
                  100 * max(sick1 / pmax(st$pop[, , 1], 1))))
  if (prev_cnt < RHD_PREV_LO || prev_cnt > RHD_PREV_HI)
    stop(loc, ": seeded RHD prevalence ", round(prev_cnt),
         " outside the sane band ", RHD_PREV_LO, "-", RHD_PREV_HI, ".", call. = FALSE)
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
    stages = c("A", "B", "C", "D"),
    ramp_start = cov$ramp_start, ramp_end = cov$ramp_end,
    incidence_trend = dmi$meta$incidence_trend,
    rhd_d_fraction  = dmi$meta$rhd_d_fraction,
    stage_split     = dmi$stage_split,
    calib_last_year = calib_last_year,
    stage_calibration_status = calib$stage_calibration$status,
    RATE_BASE_YEAR  = dmi$meta$RATE_BASE_YEAR,
    # base-year death sanity bands consumed by 06 (COUNTRY-specific; 00 via getp)
    rhd_death_lo  = getp("RHD_DEATH_LO",  1e3), rhd_death_hi  = getp("RHD_DEATH_HI",  1e5),
    allc_death_lo = getp("ALLC_DEATH_LO", 5e5), allc_death_hi = getp("ALLC_DEATH_HI", 4e6),
    intervention_labels = c(ref = "reference_cascade_held_at_baseline",
                            sap = "echo_screening_diagnosis_SAP_scale_up"),
    built_from = c(basename(IN_POP), basename(IN_DISEASE), basename(IN_CALIB))
  )
)

saveRDS(baseline_state, file = OUT_FILE)

message("── Saved ──────────────────────────────────────")
message(sprintf("  %s  | locations: %s | horizon %d-%d | scenarios: %s | stages: A/B/C/D",
                basename(OUT_FILE), paste(LOCATIONS, collapse = ", "),
                min(years), max(years), paste(SCENARIOS, collapse = ", ")))
message("── 05_build_baseline.R complete ───────────────────────")
message("  Next: 06_run_prevention_model.R")
