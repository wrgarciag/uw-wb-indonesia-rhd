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
#     data/pop_projection.rds                  (02) population [age x sex x year]
#     data/disease_model_inputs.rds           (03) rate arrays + A/B/C/D natural
#                                                  history + effects + surgery +
#                                                  stage_split + cascade structure;
#                                                  meta$years = the ANALYSIS horizon
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

IN_POP     <- paste0(wd_data, "pop_projection.rds")            # STABLE name (02)
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
years <- dmi$meta$years                    # the ANALYSIS horizon (set by 03)
n_age <- length(AGES); n_sex <- length(SEXES); n_years <- length(years)
base_year <- min(years)

# analysis-period validation: when 00_run_all.R is driving, the horizon inherited
# from 03 must equal ANALYSIS_YEARS exactly (guards against a stale disease-inputs
# bundle built under a different window).
if (exists("ANALYSIS_YEARS", inherits = TRUE)) {
  ay <- sort(as.integer(get("ANALYSIS_YEARS", inherits = TRUE)))
  if (!identical(sort(as.integer(years)), ay))
    stop("Baseline horizon (", min(years), "-", max(years),
         ") != ANALYSIS_YEARS (", min(ay), "-", max(ay),
         "). Re-run 03_build_disease_model.R under the current window.", call. = FALSE)
}

# calibrated IR/CF table (Layer 1) — a single self-describing bundle now.
calib_tp <- as.data.table(calib$tp)
calib_last_year <- calib$meta$calib_last_year         # last calibration year (2019)

message(sprintf("  Loaded: pop(%d rows) | disease inputs(horizon %d-%d) | calibrated IR/CF(%d rows, %d-%d) | stage-calib %s",
                nrow(pop_long), min(years), max(years),
                nrow(calib_tp), min(calib_tp$year), max(calib_tp$year),
                calib$stage_calibration$status))

# ------------------------------------------------------------------------------
# 1b. CHOOSE STAGE PARAMETERS: 04b CALIBRATED (baked once) or 03 priors (fallback)
#     04b (structural) writes $stage_calibration$best_parameters when it runs; we
#     consume them AS-IS (no further scaling — the former 1/2.5//1/8 mortality
#     scalars are gone). If absent, fall back to 03's clinical priors.
# ------------------------------------------------------------------------------
sc <- calib$stage_calibration
if (isTRUE(grepl("^FAILED", as.character(sc$status))))
  stop("04b structural calibration FAILED validation (status '", sc$status, "'): ",
       paste(sc$failure_reason, collapse = "; "),
       "\n  A failed candidate is never consumed. Fix the calibration or re-run 04b.", call. = FALSE)
use_calibrated_stages <- !is.null(sc$best_parameters) &&
  isTRUE(grepl("^calibrated", as.character(sc$status)))
if (use_calibrated_stages) {
  bp <- sc$best_parameters
  # MECHANISM COMPATIBILITY: a calibration fitted under one mortality structure /
  # SAP mechanism is never consumed under another (05 would otherwise mix specs).
  if (!is.null(bp$mortality_structure) && !identical(bp$mortality_structure, dmi$meta$mortality_structure))
    stop("Calibrated bundle was fitted under mortality_structure '", bp$mortality_structure,
         "' but disease inputs use '", dmi$meta$mortality_structure, "'. Re-run 03 -> 04 -> 04b under ONE spec.", call. = FALSE)
  if (!is.null(bp$sap_mechanism) && !identical(bp$sap_mechanism, dmi$effects$sap_mechanism))
    stop("Calibrated bundle was fitted under sap_mechanism '", bp$sap_mechanism,
         "' but disease inputs use '", dmi$effects$sap_mechanism, "'. Recalibrate under the selected mechanism.", call. = FALSE)
  if (identical(dmi$meta$mortality_structure, "D_only") && any(bp$p_rhd_death[c("A", "B", "C")] != 0))
    stop("D_only spec but calibrated p_rhd_death A/B/C are not exactly zero.", call. = FALSE)
  stage_transitions <- bp$transitions
  stage_p_rhd_death <- bp$p_rhd_death
  stage_split_use   <- bp$stage_split
  rhd_d_fraction_use <- bp$rhd_d_fraction
  message(sprintf("  Stage params: CALIBRATED by 04b (status %s; %s mortality; SAP %s) | p_rhd_death A/B/C/D = %s",
                  sc$status, if (is.null(bp$mortality_structure)) "legacy" else bp$mortality_structure,
                  if (is.null(bp$sap_mechanism)) "legacy" else bp$sap_mechanism,
                  paste(sprintf("%.4f", stage_p_rhd_death[c("A","B","C","D")]), collapse="/")))
} else {
  stage_transitions <- dmi$transitions
  stage_p_rhd_death <- dmi$p_rhd_death
  stage_split_use   <- dmi$stage_split
  rhd_d_fraction_use <- dmi$meta$rhd_d_fraction
  message(sprintf("  Stage params: UNCALIBRATED from 03 (status %s) -- 04b did not bake stage params.",
                  sc$status))
}

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
  cf_base <- mat_from_long(cb, "CF")   # Layer-1 aggregate case-fatality: LEGACY/DIAGNOSTIC anchor only
  # FINAL incidence: prefer 04b's jointly estimated per-SUSCEPTIBLE probability
  # (alpha applied ONCE in 04b); fall back to the Layer-1 tp IR (per total pop).
  use_l2_ir <- use_calibrated_stages && !is.null(sc$best_parameters$ir_calibrated)
  if (use_l2_ir) {
    ir_base <- sc$best_parameters$ir_calibrated
    if (!identical(dim(ir_base), c(n_age, n_sex)))
      stop("04b ir_calibrated has wrong dimensions for the production grid.", call. = FALSE)
    ir_base <- ir_base[as.character(AGES), SEXES]
    ir_anchor_year <- as.integer(sc$best_parameters$ir_anchor_year)
    incidence_source <- "layer2_joint (04b; per susceptible; alpha applied once)"
  } else {
    ir_base <- mat_from_long(cb, "IR")   # Layer-1 incidence probability (legacy path)
    ir_anchor_year <- as.integer(calib_last_year)
    incidence_source <- "layer1_tp (legacy; per total population)"
  }

  ## 3c. incidence array over horizon: calibrated base pattern x secular trend,
  ##     ANCHORED at an explicit calendar year (ir_anchor_year = CAL_YEAR_END), so
  ##     changing the analysis start never redefines the incidence trajectory.
  trend <- dmi$meta$incidence_trend
  ir_arr <- array(0, dim = c(n_age, n_sex, n_years),
                  dimnames = list(AGES, SEXES, years))
  for (iy in seq_len(n_years))
    ir_arr[, , iy] <- ir_base * trend^(years[iy] - ir_anchor_year)

  ## 3d. background (non-RHD) mortality + prevalence seed from 03 -------------
  oth_mort  <- dmi$rates_by_year$oth_mort[, , 1]   # base-year pattern (held)
  prev_seed <- dmi$rates_by_year$prev_seed         # [age x sex] prevalence fraction

  ## 3e. seed the prevalent pool at year 1 ------------------------------------
  ##     PRIMARY: 04b's analysis-start INITIAL STATE (calibrated end-of-window
  ##     state forward-projected through the bridge years with reference care) —
  ##     a dynamically consistent seed; no 2019-shares x 2023-prevalence x 2026-
  ##     population mixing. LEGACY fallback: GBD prevalence x population split by
  ##     04b's per-age stage shares (or 03's uniform split).
  prev_pool <- prev_seed * pop_arr[, , 1]
  ist <- if (use_calibrated_stages) sc$best_parameters$initial_state else NULL
  ssa <- if (use_calibrated_stages) sc$best_parameters$stage_shares_by_age else NULL
  if (!is.null(ist)) {
    if (!identical(as.integer(ist$year), as.integer(base_year)))
      stop("04b initial_state is for year ", ist$year, " but the analysis starts in ", base_year,
           ". Re-run 04b under the current ANALYSIS_YEAR_START.", call. = FALSE)
    seed <- lapply(ist[c("A", "B", "C", "D")], function(m) {
      if (!identical(dim(m), c(n_age, n_sex))) stop("04b initial_state has wrong dimensions.", call. = FALSE)
      m <- m[as.character(AGES), SEXES]; dimnames(m) <- list(AGES, SEXES); m })
    seed_source <- "layer2_initial_state (04b bridge)"
    message(sprintf("  Seed: 04b analysis-start state (%s) | total %.0f vs legacy GBD-prev x pop seed %.0f (ratio %.2f)",
                    ist$method, sum(seed$A + seed$B + seed$C + seed$D), sum(prev_pool),
                    sum(seed$A + seed$B + seed$C + seed$D) / sum(prev_pool)))
  } else if (!is.null(ssa)) {
    seed_source <- "legacy: GBD prevalence x pop x 04b per-age stage shares"
    ssa <- as.data.table(ssa)
    cal_age <- pmin(AGES, max(ssa$age))                 # ages beyond 95 inherit the 95+ fractions
    fmat <- function(col) {
      m <- empty_mat()
      for (s in SEXES) { ds <- ssa[sex == s]; m[, s] <- ds[[col]][match(cal_age, ds$age)] }
      m[is.na(m)] <- 0; m
    }
    seed <- list(A = prev_pool * fmat("fracA"), B = prev_pool * fmat("fracB"),
                 C = prev_pool * fmat("fracC"), D = prev_pool * fmat("fracD"))
  } else {
    seed_source <- "legacy: GBD prevalence x pop x 03 uniform stage split"
    ss <- stage_split_use
    seed <- list(A = prev_pool * ss[["A"]], B = prev_pool * ss[["B"]],
                 C = prev_pool * ss[["C"]], D = prev_pool * ss[["D"]])
  }

  ## 3e'. engine-side transitions: scalar list (contract) + optional explicit age
  ##      modifier on A -> No-RHD regression (remission_adult_mult at ages >= cutoff),
  ##      identical to the 04b calibration convention.
  rem_mult <- if (is.null(dmi$meta$remission_adult_mult)) 1 else dmi$meta$remission_adult_mult
  rem_cut  <- if (is.null(dmi$meta$remission_age_cutoff)) 20L else dmi$meta$remission_age_cutoff
  if (use_calibrated_stages && !is.null(bp$remission_adult_mult) &&
      (bp$remission_adult_mult != rem_mult || bp$remission_age_cutoff != rem_cut))
    stop("remission_adult_mult / cutoff differ between 04b calibration and 03 inputs; re-run 03 -> 04b.", call. = FALSE)
  # age structure: 04b's calibrated adult multipliers (progression + remission,
  # hazard scale) when present; else the fixed remission hook from 03
  ag <- if (use_calibrated_stages) bp$age_structure else NULL
  age_mults <- if (!is.null(ag)) c(prog = ag$adult_prog_mult, rem = ag$adult_remission_mult) else c(prog = 1, rem = rem_mult)
  age_cut   <- if (!is.null(ag)) ag$adult_age_cutoff else rem_cut
  transitions_engine <- stage_transitions
  mk_age <- function(p, m) { M <- empty_mat() + p; M[AGES >= age_cut, ] <- 1 - (1 - p)^m; M }
  if (age_mults[["rem"]] != 1) transitions_engine$p_A_to_no_rhd <- mk_age(stage_transitions$p_A_to_no_rhd, age_mults[["rem"]])
  if (age_mults[["prog"]] != 1) {
    transitions_engine$p_A_to_B <- mk_age(stage_transitions$p_A_to_B, age_mults[["prog"]])
    transitions_engine$p_B_to_C <- mk_age(stage_transitions$p_B_to_C, age_mults[["prog"]])
    transitions_engine$p_C_to_D <- mk_age(stage_transitions$p_C_to_D, age_mults[["prog"]])
  }

  ## 3f. realised cascade + surgery coverage trajectories per scenario --------
  #  Reference holds baselines flat; scale-up ramps to targets over [ramp_start,
  #  ramp_end] (from 03's coverage$ramp_*), held at target thereafter. Effective
  #  diagnosis/treatment are capped by the earlier cascade stages.
  cascade_for <- function(arm) {                          # arm = "ref" | "up"
    screen    <- ramp_traj(cov[[paste0("screen_",    arm, "_baseline")]],
                           cov[[paste0("screen_",    arm, "_target")]],   cov$ramp_start, cov$ramp_end)
    diagnosis <- ramp_traj(cov[[paste0("diagnosis_", arm, "_baseline")]],
                           cov[[paste0("diagnosis_", arm, "_target")]],   cov$ramp_start, cov$ramp_end)
    treatment <- ramp_traj(cov[[paste0("treatment_", arm, "_baseline")]],
                           cov[[paste0("treatment_", arm, "_target")]],   cov$ramp_start, cov$ramp_end)
    surgery   <- ramp_traj(surg[[paste0("coverage_", arm, "_baseline")]],
                           surg[[paste0("coverage_", arm, "_target")]],   cov$ramp_start, cov$ramp_end)
    # annual population-average SAP adherence modifier (default 1; 00 block H)
    adh_b <- cov[[paste0("adherence_", arm, "_baseline")]]; adh_t <- cov[[paste0("adherence_", arm, "_target")]]
    adherence <- ramp_traj(if (is.null(adh_b)) 1 else adh_b, if (is.null(adh_t)) 1 else adh_t,
                           cov$ramp_start, cov$ramp_end)
    # cumulative-coverage logic: effective downstream <= upstream cascade stages.
    # INTERPRETATION (documented): screen/diagnosis/treatment are CUMULATIVE
    # population-level coverages of the eligible (living RHD) pool, so taking the
    # minimum caps each downstream stage by its upstream stage; they are NOT
    # multiplied as conditional proportions. eff_treatment = share of living RHD on
    # effective (adherent) optimal treatment = DELIVERED SAP; effective PROTECTION
    # = eff_treatment x adherence (adherence is a relative attenuation, default 1).
    eff_diag <- pmin(screen, diagnosis)
    eff_trt  <- pmin(screen, diagnosis, treatment)
    list(screen = screen, diagnosis = diagnosis, treatment = treatment,
         eff_diagnosis = eff_diag, eff_treatment = eff_trt, surgery = surgery,
         adherence = adherence, eff_protection = eff_trt * adherence)
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
       transitions = stage_transitions, p_rhd_death = stage_p_rhd_death,
       effects = dmi$effects, surgery = surg,
       coverage = coverage,
       # additive fields (06 consumes; legacy fields above unchanged in meaning)
       transitions_engine = transitions_engine,       # scalar or [age x sex] (adult age multipliers)
       age_structure = list(adult_prog_mult = unname(age_mults[["prog"]]), adult_remission_mult = unname(age_mults[["rem"]]),
                            adult_age_cutoff = age_cut, source = if (!is.null(ag)) "04b calibrated" else "03 fixed hook"),
       seed_source = seed_source,
       incidence_source = incidence_source, ir_anchor_year = ir_anchor_year,
       cf_role = "Layer-1 aggregate case-fatality anchor: LEGACY/DIAGNOSTIC only; A/B/C/D deaths come from p_rhd_death")
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
    for (k in c("screen", "diagnosis", "treatment", "eff_diagnosis", "eff_treatment", "surgery",
                "adherence", "eff_protection")) {
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

  # mortality structure of the values 06 will run
  if (identical(dmi$meta$mortality_structure, "D_only") && any(st$p_rhd_death[c("A", "B", "C")] != 0))
    stop(loc, ": D_only spec but p_rhd_death A/B/C not exactly zero.", call. = FALSE)
  # reference care used in 04b must equal the reference arm run by 06 (applied once)
  if (use_calibrated_stages && !is.null(bp$reference_care_in_calibration)) {
    rc <- bp$reference_care_in_calibration
    if (abs(rc$treatment_effective - st$coverage$ref$eff_treatment[1]) > 1e-9 ||
        abs(rc$adherence - st$coverage$ref$adherence[1]) > 1e-9 ||
        abs(rc$surgery - st$coverage$ref$surgery[1]) > 1e-9)
      stop(loc, ": reference care in 04b (treatment ", rc$treatment_effective, ", adherence ", rc$adherence,
           ", surgery ", rc$surgery, ") differs from the 06 reference arm; recalibrate.", call. = FALSE)
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
    rhd_d_fraction  = rhd_d_fraction_use,
    stage_split     = stage_split_use,
    calib_last_year = calib_last_year,
    stage_calibration_status = calib$stage_calibration$status,
    RATE_BASE_YEAR  = dmi$meta$RATE_BASE_YEAR,
    # spec carried to 06 and its outputs (additive)
    mortality_structure = dmi$meta$mortality_structure,
    sap_mechanism = dmi$effects$sap_mechanism,
    sap_hr = unlist(dmi$effects[c("sap_hr_A_to_B", "sap_hr_B_to_C", "sap_hr_C_to_D")]),
    sap_effect_scale = dmi$effects$sap_effect_scale,
    sap_effect_sources = dmi$effects$sap_effect_sources,
    baseline_care_interpretation = dmi$meta$baseline_care_interpretation,
    adherence_convention = "annual population-average relative attenuation (default 1 = unattenuated effect; optimal-treatment coverage already denotes effective, adherent treatment)",
    incidence_trend_anchor_year = states[[1]]$ir_anchor_year,
    incidence_source = states[[1]]$incidence_source,
    seed_source = states[[1]]$seed_source,
    remission_adult_mult = dmi$meta$remission_adult_mult, remission_age_cutoff = dmi$meta$remission_age_cutoff,
    age_structure = states[[1]]$age_structure,
    transitions_note = "$states[[loc]]$transitions are the child/reference-age scalar probabilities; 06 runs $transitions_engine (adult hazard multipliers at ages >= adult_age_cutoff)",
    population_backbone_approximation = paste(
      "Population is EXOGENOUS (WPP2024 by age x sex x year); RHD stocks age by cohort; the No-RHD pool is",
      "population minus RHD stock each year (it absorbs migration and population-source differences);",
      "births enter age 0 with zero RHD; stock > population is surfaced as a violation, never clipped silently."),
    eff_cf_role = if (identical(dmi$effects$sap_mechanism, "progression"))
      "progression mode: no direct SAP mortality effect, eff_cf == 1 by construction (legacy column kept)" else
      "direct_mortality mode: eff_cf = 1 - sap_rrr_rhd_death x effective protection",
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
