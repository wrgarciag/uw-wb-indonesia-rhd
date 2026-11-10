# ==============================================================================
# RHD secondary-prevention investment case: INITIAL-STATE ASSEMBLER (A/B/C/D)
# scripts/05_build_baseline.R
#
# Default natural history (engine abcd_first_diagnosis_v1_adult_prog_mult, R/abcd_engine.R): WHF
# stages A/B/C/D held as undiagnosed / confirmed-diagnosis (off / on BPG) stocks,
# adult prior-A and normalised history, attained-age rules (<= 20 / > 20) and
# D-only mortality. Surgery is a SERVICE, not a state.
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT
# ------------------------------------------------------------------------------
# Assembles the INITIAL STATE and every input needed to RUN the model in 06
# WITHOUT re-running 01-04b. It reads only the PERSISTED upstream outputs:
#
#   INPUTS
#     data/pop_projection.rds                  (02) population [age x sex x year]
#     data/disease_model_inputs.rds           (03) rate arrays, natural history, effects,
#                                                  surgery, pooled care STRUCTURE; meta$years
#     data/calibrated_rhd_parameters.rds       (04/04b) $stage_calibration$best_parameters:
#                                                  calibrated natural history, FINAL incidence
#                                                  and the analysis-start engine state
#
#   OUTPUT (written to wd_data = data/):
#     baseline_state.rds — per-location bundle holding:
#       pop        [age x sex x year] population over the horizon (from 02)
#       ir         [age x sex x year] incidence per TOTAL population (04b, x secular trend);
#                                     the engine converts it per never-onset susceptible
#       cf         [age x sex]        Layer-1 aggregate case fatality: LEGACY/DIAGNOSTIC only
#       oth_mort   [age x sex]        background (non-RHD) competing mortality (03)
#       seed       list(A,B,C,D)      aggregate classifiable stocks ENTERING the first analysis year
#       seed_state engine state       full engine state (undiagnosed / diagnosed / on-BPG stage
#                                     stocks, adult prior-A, normalised, never-onset, screening
#                                     histories) = 06 seed
#       transitions / natural_history / p_rhd_death / effects / surgery
#       coverage   REALISED per-year pooled care trajectories per arm (ref, sap)
#     plus a `meta` block (AGES, SEXES, years, scenarios, ramp window, spec, definitions).
#
# The calibrated natural history and incidence are consumed AS-IS (baked once in
# 04b). Prevalence and stage mix are NOT reset to GBD at the analysis start: the seed
# is the calibrated state carried through the bridge years under reference care.
# NO monetary / cost values live here — economics are in 08 only.
# ==============================================================================

library(data.table)
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")
eng <- if (exists("wd")) paste0(wd, "R/abcd_engine.R") else here::here("R", "abcd_engine.R")
source(eng)

# getp(): honour a 00_run_all.R global when set, else the standalone default.
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
RHD_PREV_LO <- getp("RHD_PREV_LO", 1e5)   # seeded-prevalence sanity band (COUNTRY-settable)
RHD_PREV_HI <- getp("RHD_PREV_HI", 1e7)

# ------------------------------------------------------------------------------
# 0. CONFIG
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
         "\n  Run the upstream script (02 / 03 / 04 / 04b) first.", call. = FALSE)

# ------------------------------------------------------------------------------
# 1. LOAD PERSISTED UPSTREAM OUTPUTS
# ------------------------------------------------------------------------------
message("── 05_build_baseline.R : assembling the analysis-start state ─────────")

dmi      <- readRDS(IN_DISEASE)
pop_long <- as.data.table(readRDS(IN_POP))
calib    <- readRDS(IN_CALIB)

AGES  <- dmi$meta$AGES
SEXES <- dmi$meta$SEXES
years <- dmi$meta$years                    # the ANALYSIS horizon (set by 03)
n_age <- length(AGES); n_sex <- length(SEXES); n_years <- length(years)
base_year <- min(years)

if (exists("ANALYSIS_YEARS", inherits = TRUE)) {
  ay <- sort(as.integer(get("ANALYSIS_YEARS", inherits = TRUE)))
  if (!identical(sort(as.integer(years)), ay))
    stop("Baseline horizon (", min(years), "-", max(years),
         ") != ANALYSIS_YEARS (", min(ay), "-", max(ay),
         "). Re-run 03_build_disease_model.R under the current window.", call. = FALSE)
}

calib_tp <- as.data.table(calib$tp)
calib_last_year <- calib$meta$calib_last_year

# ------------------------------------------------------------------------------
# 1b. CALIBRATED NATURAL HISTORY (04b) — required; consumed AS-IS
#     Stale-bundle protection: the bundle must have been fitted by THIS engine under
#     the SAME fixed natural-history inputs and the same reference care as 03 now
#     describes; a FAILED or uncalibrated bundle is never consumed.
# ------------------------------------------------------------------------------
sc <- calib$stage_calibration
if (isTRUE(grepl("^FAILED", as.character(sc$status))))
  stop("04b structural calibration FAILED validation (status '", sc$status, "'): ",
       paste(sc$failure_reason, collapse = "; "),
       "\n  A failed candidate is never consumed. Fix the calibration or re-run 04b.", call. = FALSE)
if (is.null(sc$best_parameters) || !isTRUE(grepl("^calibrated", as.character(sc$status))))
  stop("No calibrated natural history in ", basename(IN_CALIB), " (status '", sc$status, "'). The default ",
       "model has no valid analysis-start state without 04b; run 04b_calibrate_structural.R.", call. = FALSE)
bp <- sc$best_parameters
if (!identical(bp$engine_spec, ABCD_ENGINE_SPEC) || !identical(dmi$meta$engine_spec, ABCD_ENGINE_SPEC))
  stop("Calibrated bundle engine '", if (is.null(bp$engine_spec)) "legacy (none recorded)" else bp$engine_spec,
       "' / disease inputs '", dmi$meta$engine_spec, "' differ from the engine '", ABCD_ENGINE_SPEC,
       "'. A calibration from another structure is never consumed: re-run 03 -> 04 -> 04b.", call. = FALSE)
if (!identical(bp$natural_history_fixed, dmi$meta$natural_history_fixed))
  stop("Fixed natural-history inputs changed since the calibration (00 block G); re-run 03 -> 04 -> 04b.", call. = FALSE)
# the natural history 06 will actually run must carry the same fixed inputs (incl. the adult
# progression hazard multiplier) and 00's current multiplier
if (!identical(abcd_nh_fixed(bp$natural_history), dmi$meta$natural_history_fixed) ||
    (exists("adult_prog_hazard_mult") && !identical(bp$natural_history$adult_prog_hazard_mult, adult_prog_hazard_mult)))
  stop("Calibrated natural history (adult_prog_hazard_mult = ", format(bp$natural_history$adult_prog_hazard_mult),
       ") differs from the disease inputs / 00 block G; re-run 03 -> 04 -> 04b.", call. = FALSE)
# fixed AGE-SPECIFIC stage-D mortality: the calibrated natural history must carry 00's current
# child / adult probabilities and mortality boundary (no single all-age value is accepted)
D_MORT_FIELDS <- c("p_rhd_death_D_child", "p_rhd_death_D_adult", "d_death_adult_age_min")
d_mort_bad <- !is.null(bp$natural_history$p_rhd_death_D) ||
  any(vapply(D_MORT_FIELDS, function(nm) is.null(bp$natural_history[[nm]]) ||
               (exists(nm) && !identical(as.numeric(bp$natural_history[[nm]]), as.numeric(get(nm)))), logical(1)))
if (d_mort_bad)
  stop("Calibrated stage-D mortality (", paste(sprintf("%s = %s", D_MORT_FIELDS, vapply(D_MORT_FIELDS, function(nm) { v <- bp$natural_history[[nm]]; if (is.null(v)) "missing" else format(v) }, "")), collapse = ", "),
       ") is missing, a legacy single value, or differs from 00 block G; re-run 03 -> 04 -> 04b.", call. = FALSE)
if (!identical(bp$mortality_structure, dmi$meta$mortality_structure) || !identical(bp$sap_mechanism, dmi$effects$sap_mechanism))
  stop("Calibrated bundle spec (", bp$mortality_structure, "/", bp$sap_mechanism, ") differs from the disease inputs (",
       dmi$meta$mortality_structure, "/", dmi$effects$sap_mechanism, "); recalibrate under one spec.", call. = FALSE)
if (any(bp$p_rhd_death[c("A", "B", "C")] != 0))
  stop("D_only spec but calibrated p_rhd_death A/B/C are not exactly zero.", call. = FALSE)
burnin_now <- getp("STRUCT_BURNIN_RECONCILIATION", "per_capita")
if (is.null(bp$initialisation$burnin_reconciliation) || !identical(bp$initialisation$burnin_reconciliation, burnin_now))
  stop("Calibrated bundle burn-in reconciliation '", if (is.null(bp$initialisation$burnin_reconciliation)) "not recorded" else bp$initialisation$burnin_reconciliation,
       "' differs from STRUCT_BURNIN_RECONCILIATION = '", burnin_now, "' (00 block G2); re-run 04b.", call. = FALSE)
if (is.null(bp$initial_state) || is.null(bp$initial_state$engine_state))
  stop("04b best_parameters has no analysis-start engine state (bridge disabled, or a bundle from another engine). Set STRUCT_BRIDGE_TO_ANALYSIS_START = TRUE and re-run 04b.", call. = FALSE)
SAP_SPEC_FIELDS <- c("sap_progression_rrr", "sap_progression_hazard_mult", "sap_protection_age_max",
                     "sap_effect_transitions", "sap_hr_prior_A_to_B_adult", "sap_hr_C_to_D", "sap_hr_D_to_rhd_death")
sap_missing <- c(setdiff(SAP_SPEC_FIELDS, names(sc$spec$sap_hr)[!vapply(sc$spec$sap_hr, is.null, logical(1))]),
                 setdiff(SAP_SPEC_FIELDS, names(dmi$effects)[!vapply(dmi$effects, is.null, logical(1))]))
if (length(sap_missing))
  stop("SAP specification incomplete in the 04b bundle or the disease inputs (missing: ", paste(unique(sap_missing), collapse = ", "),
       "); old specifications are not accepted -- re-run 03 -> 04 -> 04b.", call. = FALSE)
if (!identical(sc$spec$sap_eligibility, dmi$effects$sap_eligibility) ||
    !identical(sc$spec$sap_hr[SAP_SPEC_FIELDS], dmi$effects[SAP_SPEC_FIELDS]))
  stop("SAP eligibility / progression-multiplier specification used as reference care in 04b differs from the disease inputs; re-run 04b.", call. = FALSE)
message(sprintf("  natural history: CALIBRATED by 04b (status %s; engine %s; burn-in %s) | transitions %s | D death %.4f at age < %d / %.4f at age >= %d (fixed, all histories)",
                sc$status, bp$engine_spec, bp$initialisation$burnin_reconciliation,
                paste(sprintf("%s=%.4f", names(bp$transitions), unlist(bp$transitions)), collapse = ", "),
                bp$natural_history$p_rhd_death_D_child, bp$natural_history$d_death_adult_age_min,
                bp$natural_history$p_rhd_death_D_adult, bp$natural_history$d_death_adult_age_min))
if (!isTRUE(sc$statuses$fit_within_tolerance))
  message("  NOTE: the calibration is OUTSIDE its diagnostic fit tolerances (", sc$statuses$epidemiological_fit,
          ") -- results are a numerically complete run, not a credible calibrated baseline.")

# ------------------------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------------------------
empty_mat <- function() matrix(0, n_age, n_sex, dimnames = list(AGES, SEXES))
mat_from_long <- function(dt, valcol) {
  m <- empty_mat()
  for (s in SEXES) {
    ds <- dt[sex == s]
    if (!nrow(ds)) next
    idx <- match(ds$age, AGES); keep <- !is.na(idx)
    m[idx[keep], s] <- ds[[valcol]][keep]
    top_age <- max(ds$age)
    if (top_age < max(AGES)) m[which(AGES > top_age), s] <- m[as.character(top_age), s]
  }
  m
}
# baseline before window, linear across [start,end], target after
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
care_vec_names <- c("screen_window_prob", "clinical_detection", "bpg_initiation", "bpg_retention",
                    "adherence_good", "adherence_partial")

build_location_state <- function(loc) {

  ## 3a. population array [age x sex x year] from 02
  pl <- pop_long[location == loc & age %in% AGES & sex %in% SEXES]
  if (!nrow(pl)) stop("No population rows for location '", loc, "'.", call. = FALSE)
  pop_arr <- array(0, dim = c(n_age, n_sex, n_years), dimnames = list(AGES, SEXES, years))
  for (iy in seq_len(n_years)) {
    py <- pl[year == years[iy]]; m <- empty_mat()
    for (s in SEXES) { ps <- py[sex == s]; m[match(ps$age, AGES), s] <- ps$Nx }
    pop_arr[, , iy] <- m
  }

  ## 3b. FINAL incidence (04b; per TOTAL population, alpha applied once) x secular
  ##     trend anchored at an explicit calendar year (ir_anchor_year = CAL_YEAR_END)
  ir_base <- bp$ir_calibrated
  if (!identical(dim(ir_base), c(n_age, n_sex)))
    stop("04b ir_calibrated has wrong dimensions for the production grid.", call. = FALSE)
  ir_base <- ir_base[as.character(AGES), SEXES]
  ir_anchor_year <- as.integer(bp$ir_anchor_year)
  trend <- dmi$meta$incidence_trend
  ir_arr <- array(0, dim = c(n_age, n_sex, n_years), dimnames = list(AGES, SEXES, years))
  for (iy in seq_len(n_years)) ir_arr[, , iy] <- ir_base * trend^(years[iy] - ir_anchor_year)

  ## 3c. Layer-1 CF (legacy diagnostic only) + background mortality (03, held)
  cb <- calib_tp[location == loc & year == calib_last_year]
  if (!nrow(cb)) stop("No Layer-1 rows for location '", loc, "' at year ", calib_last_year, ".", call. = FALSE)
  cf_base  <- mat_from_long(cb, "CF")
  oth_mort <- dmi$rates_by_year$oth_mort[, , 1]

  ## 3d. analysis-start seed = 04b bridged state (no GBD reseeding)
  ist <- bp$initial_state
  if (!identical(as.integer(ist$year), as.integer(base_year)))
    stop("04b initial_state is for year ", ist$year, " but the analysis starts in ", base_year,
         ". Re-run 04b under the current ANALYSIS_YEAR_START.", call. = FALSE)
  if (!identical(sort(names(ist$engine_state)), sort(ABCD_STATE_FIELDS)))
    stop("04b initial_state engine fields differ from engine '", ABCD_ENGINE_SPEC, "'; re-run 04b.", call. = FALSE)
  seed_state <- lapply(ist$engine_state[ABCD_STATE_FIELDS], function(m) {
    if (!identical(dim(m), c(n_age, n_sex))) stop("04b initial_state has wrong dimensions.", call. = FALSE)
    m <- m[as.character(AGES), SEXES]; dimnames(m) <- list(AGES, SEXES); m })
  seed <- lapply(abcd_stage_aggregate(seed_state), function(m) { dimnames(m) <- list(AGES, SEXES); m })
  for (s in c("A", "B", "C", "D"))
    if (max(abs(seed[[s]] - ist[[s]][as.character(AGES), SEXES])) > 1e-6)
      stop("04b initial_state aggregate ", s, " does not equal its engine-state detail.", call. = FALSE)
  prev_pool <- dmi$rates_by_year$prev_seed * pop_arr[, , 1]          # legacy comparison only
  message(sprintf("  Seed: 04b analysis-start state (%s) | classifiable RHD %.0f (legacy GBD-prev x pop %.0f, ratio %.2f) | prior-A %.0f | normalised %.0f",
                  ist$method, sum(abcd_rhd_total(seed_state)), sum(prev_pool),
                  sum(abcd_rhd_total(seed_state)) / sum(prev_pool), sum(abcd_prior_A_total(seed_state)),
                  sum(abcd_normalised_total(seed_state))))

  ## 3e. realised pooled-care trajectories per arm (reference flat; scale-up ramped
  ##     over [ramp_start, ramp_end], held after). q is derived from the CURRENT
  ##     year's window target; realised cohort reach during the ramp is tracked by the
  ##     engine's never-screened share (it need not equal the current target).
  W <- cov$screen_window_years
  care_for <- function(arm) {
    out <- lapply(care_vec_names, function(nm) {
      v <- cov$care_inputs[[nm]][[arm]]
      ramp_traj(v[["baseline"]], v[["target"]], cov$ramp_start, cov$ramp_end) })
    names(out) <- care_vec_names
    out$screen_annual_prob <- abcd_screen_annual_prob(out$screen_window_prob, W)
    out$adherence_factor <- out$adherence_good + cov$partial_adherence_relative_protection * out$adherence_partial
    out$surgery <- ramp_traj(surg[[paste0("coverage_", arm, "_baseline")]], surg[[paste0("coverage_", arm, "_target")]],
                             cov$ramp_start, cov$ramp_end)
    out
  }
  coverage <- list(
    ref = care_for("ref"),
    sap = care_for("up"),
    constants = list(sensitivity = cov$screen_sensitivity, confirm_completion = cov$confirm_completion,
                     confirm_sensitivity = cov$confirm_sensitivity, partial_relative_protection = cov$partial_adherence_relative_protection,
                     screen_age_lo = cov$screen_age_lo, screen_age_hi = cov$screen_age_hi, screen_window_years = W),
    init = cov$init,
    eligibility = cov$eligibility)

  list(pop = pop_arr, ir = ir_arr, cf = cf_base, oth_mort = oth_mort,
       seed = seed, seed_state = seed_state,
       transitions = bp$transitions, natural_history = bp$natural_history,
       p_rhd_death = bp$p_rhd_death,
       effects = dmi$effects, surgery = surg,
       coverage = coverage,
       seed_source = "layer2_initial_state (04b bridge; full stage x care state, prior-A, normalised, never-onset and screening histories)",
       incidence_source = "layer2_joint (04b; per total population; alpha applied once; engine converts per never-onset susceptible)",
       ir_anchor_year = ir_anchor_year,
       cf_role = "Layer-1 aggregate case-fatality: LEGACY/DIAGNOSTIC only; never overrides the fixed D-only mortality")
}

states <- setNames(lapply(LOCATIONS, build_location_state), LOCATIONS)

# ------------------------------------------------------------------------------
# 4. VALIDATION  (fail loudly BEFORE writing)
# ------------------------------------------------------------------------------
message("── Validation ─────────────────────────────────")
chk <- function(x, nm) {
  if (any(!is.finite(x)))  stop(nm, ": contains NA / non-finite values.", call. = FALSE)
  if (any(x < 0))          stop(nm, ": contains negative values.", call. = FALSE)
}
chk_prob <- function(x, nm) { chk(x, nm); if (any(x > 1)) stop(nm, ": contains values > 1.", call. = FALSE) }

for (loc in LOCATIONS) {
  st <- states[[loc]]
  if (!identical(dim(st$pop), c(n_age, n_sex, n_years))) stop(loc, ": population array has wrong dimensions.", call. = FALSE)
  if (!identical(dim(st$ir), c(n_age, n_sex, n_years)))  stop(loc, ": incidence array has wrong dimensions.", call. = FALSE)
  chk(st$pop, paste0(loc, " pop")); chk_prob(st$ir, paste0(loc, " ir")); chk_prob(st$cf, paste0(loc, " cf"))
  chk_prob(st$oth_mort, paste0(loc, " oth_mort"))
  for (k in ABCD_STATE_FIELDS) chk(st$seed_state[[k]], paste0(loc, " seed_state$", k))
  # classifiable RHD + adult prior-A <= year-1 population (never-onset / normalised are reconciled by the engine)
  rhd1 <- abcd_rhd_total(st$seed_state)
  if (any(rhd1 + abcd_prior_A_total(st$seed_state) > st$pop[, , 1] + 1e-6))
    stop(loc, ": seeded classifiable RHD + prior-A exceeds year-1 population somewhere.", call. = FALSE)
  adult <- AGES > st$natural_history$child_last_age
  if (any(abcd_stage_total(st$seed_state, "A")[adult, ] != 0) || any(st$seed_state$T_A[adult, ] != 0))
    stop(loc, ": seeded adult classifiable A (or A on BPG) is not exactly zero (adult A is prior-A history).", call. = FALSE)
  if (any(abcd_prior_A_total(st$seed_state)[!adult, ] != 0))
    stop(loc, ": seeded prior-A history at child ages.", call. = FALSE)
  # care trajectories in [0,1], correct length; reference flat
  for (sc_ in SCENARIOS) for (k in c(care_vec_names, "screen_annual_prob", "adherence_factor", "surgery")) {
    v <- st$coverage[[sc_]][[k]]
    if (length(v) != n_years) stop(loc, ": coverage ", sc_, "/", k, " wrong length.", call. = FALSE)
    chk_prob(v, paste0(loc, " coverage ", sc_, "/", k))
  }
  for (k in c(care_vec_names, "surgery"))
    if (diff(range(st$coverage$ref[[k]])) > 1e-12) stop(loc, ": reference care '", k, "' is not constant over the horizon.", call. = FALSE)
  # the engine accepts both arms' schedules (validates good + partial <= 1 etc.)
  for (sc_ in SCENARIOS) invisible(abcd_care_schedule(c(st$coverage[[sc_]][c(care_vec_names, "surgery")], st$coverage$constants),
                                                     n_years, AGES, SEXES))
  # reference care used in 04b == the 06 reference arm (applied once, one care spec)
  rc <- bp$reference_care_in_calibration
  ref_now <- c(lapply(st$coverage$ref[c(care_vec_names, "surgery")], `[`, 1),
               st$coverage$constants[c("sensitivity", "confirm_completion", "confirm_sensitivity", "partial_relative_protection",
                                       "screen_age_lo", "screen_age_hi")])
  rc_cmp <- rc[names(ref_now)]
  if (!isTRUE(all.equal(lapply(rc_cmp, as.numeric), lapply(ref_now, as.numeric))) || !isTRUE(all.equal(rc$init[c("diagnosed_share", "on_bpg_share", "ever_screened_share")], cov$init[c("diagnosed_share", "on_bpg_share", "ever_screened_share")])))
    stop(loc, ": reference care / initial care shares used in 04b differ from the 06 reference arm; recalibrate.", call. = FALSE)
  prev_cnt <- sum(rhd1)
  message(sprintf("  %s | seeded classifiable RHD (entering %d) ~ %s | A/B/C/D = %s | diagnosed %.1f%% | on BPG %.1f%% | adult prior-A %s",
                  loc, base_year, formatC(round(prev_cnt), format = "d", big.mark = ","),
                  paste(sprintf("%.1f%%", 100 * vapply(st$seed, sum, numeric(1)) / prev_cnt), collapse = "/"),
                  100 * sum(abcd_diagnosed_total(st$seed_state)) / prev_cnt,
                  100 * sum(st$seed_state$T_A + st$seed_state$T_B) / prev_cnt,
                  formatC(round(sum(abcd_prior_A_total(st$seed_state))), format = "d", big.mark = ",")))
  if (prev_cnt < RHD_PREV_LO || prev_cnt > RHD_PREV_HI)
    stop(loc, ": seeded RHD prevalence ", round(prev_cnt), " outside the sane band ", RHD_PREV_LO, "-", RHD_PREV_HI, ".", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 5. ASSEMBLE + PERSIST THE BASELINE STATE
# ------------------------------------------------------------------------------
s1 <- states[[1]]
baseline_state <- list(
  locations = LOCATIONS,
  states    = states,
  meta = list(
    AGES = AGES, SEXES = SEXES, years = years, base_year = base_year,
    scenarios = SCENARIOS,
    stages = c("A", "B", "C", "D"), state_fields = ABCD_STATE_FIELDS, disease_stocks = ABCD_DISEASE_STOCKS,
    prior_A_stocks = ABCD_PRIOR_A_STOCKS, normalised_stocks = ABCD_NORMALISED_STOCKS,
    engine_spec = ABCD_ENGINE_SPEC, country = LOCATIONS[1],
    # model assumptions + observation definitions + calibration record (consumed as-is; never re-adjusted)
    model_assumptions = list(
      child_last_age = dmi$natural_history$child_last_age, entry_mix_child = dmi$natural_history$entry_mix_child,
      entry_mix_adult = dmi$natural_history$entry_mix_adult, p_B_to_no_rhd_adult = dmi$natural_history$p_B_to_no_rhd_adult,
      p_rhd_death_D_child = dmi$natural_history$p_rhd_death_D_child, p_rhd_death_D_adult = dmi$natural_history$p_rhd_death_D_adult,
      d_death_adult_age_min = dmi$natural_history$d_death_adult_age_min, d_mortality = dmi$natural_history$d_mortality,
      recurrence = dmi$natural_history$recurrence_after_normalisation,
      adult_prog_hazard_mult = dmi$natural_history$adult_prog_hazard_mult,
      adult_progression = dmi$natural_history$adult_progression,
      effective_adult_transitions = unlist(bp$natural_history$derived[c("p_B_to_C_adult", "p_C_to_D_adult", "p_prior_A_to_B_adult")]),
      sap_eligibility = dmi$effects$sap_eligibility, sap_effect_sources = dmi$effects$sap_effect_sources,
      sap_progression_rrr = dmi$effects$sap_progression_rrr, sap_progression_hazard_mult = dmi$effects$sap_progression_hazard_mult,
      sap_protection_age_max = dmi$effects$sap_protection_age_max, sap_effect_transitions = dmi$effects$sap_effect_transitions,
      dx_mix_benchmark = dmi$dx_mix_benchmark, screen_specificity = cov$screen_specificity),
    observation_definitions = list(
      sick = "total CLASSIFIABLE RHD (undiagnosed + diagnosed, stages A-D); adult prior-A and normalised history excluded",
      well = "population - sick (never-onset, normalised history and adult prior-A)",
      newcases = "new modelled classifiable disease (first entry from the never-onset pool), NOT diagnoses",
      dead = "RHD deaths (diagnosed + undiagnosed D only)",
      prevalence_mapping = sc$spec$prevalence_mapping),
    calibration_record = list(target_basis = sc$spec$target_basis, weights = sc$objective_components$weights,
                              calendar = sc$spec$calendar, fit_tolerances = sc$spec$fit_tolerances,
                              statuses = sc$statuses, fit_tolerance_summary = sc$fit_tolerance_summary),
    ramp_start = cov$ramp_start, ramp_end = cov$ramp_end,
    incidence_trend = dmi$meta$incidence_trend,
    seed_stage_shares = vapply(s1$seed, sum, numeric(1)) / sum(abcd_rhd_total(s1$seed_state)),   # derived (model), not an input
    calib_last_year = calib_last_year,
    stage_calibration_status = sc$status,
    fit_within_tolerance = isTRUE(sc$statuses$fit_within_tolerance),
    calibration_adequacy = if (isTRUE(sc$statuses$fit_within_tolerance)) "within diagnostic tolerances" else
      paste("INADEQUATE (not a credible calibrated baseline): outside diagnostic tolerances --", sc$statuses$epidemiological_fit),
    burnin_reconciliation = bp$initialisation$burnin_reconciliation,
    RATE_BASE_YEAR  = dmi$meta$RATE_BASE_YEAR,
    mortality_structure = dmi$meta$mortality_structure,
    sap_mechanism = dmi$effects$sap_mechanism,
    sap_hr = sc$spec$sap_hr[SAP_SPEC_FIELDS],   # typed SAP specification used by 04b (== disease inputs, checked above)
    sap_effect_scale = dmi$effects$sap_effect_scale,
    sap_effect_sources = dmi$effects$sap_effect_sources,
    baseline_care_interpretation = dmi$meta$baseline_care_interpretation,
    natural_history_fixed = dmi$meta$natural_history_fixed,
    care_definitions = list(
      eligibility = cov$eligibility, screening = cov$screening_definition, clinical_detection = cov$clinical_detection_definition,
      sap_volume = cov$sap_volume_definition, initial_state = cov$init$note,
      care_update = paste("start of cycle: never-screened people without a prior confirmed diagnosis screened with probability q;",
                          "first confirmed diagnosis = q x sens x confirmation completion x confirmation sensitivity among never-screened",
                          "undiagnosed disease, then clinical detection k among all remaining undiagnosed (no overlap); on BPG retained r,",
                          "dropouts rejoin diagnosed-off-BPG; initiation i among eligible diagnosed off BPG (A <= 20, B);",
                          "protection on BPG: h x [1 - adh (1 - M)] on A->B and B->C at attained age <= sap_protection_age_max (M = assumed full-adherence multiplier)")),
    screening_volume_definition = "n_screened = screening ENCOUNTERS = q x never-screened people without a prior confirmed diagnosis at ages 5-39 (one screen per person; people without RHD included)",
    sap_volume_definition = "n_on_optimal_treatment = on-BPG persons (confirmed A at age <= 20 and B; start of cycle, after the care update) = SAP person-years costed in 08; protection-weighted volume is n_effectively_protected",
    confirmation_volume_definition = "n_confirm_completed = confirmatory echoes completed among true screen positives (UNCOSTED: no unit cost in 00 block J)",
    incidence_trend_anchor_year = s1$ir_anchor_year,
    incidence_source = s1$incidence_source,
    seed_source = s1$seed_source,
    transitions_note = "$states[[loc]]$transitions are the calibrated CHILDHOOD (age <= 20) annual probabilities; adult (> 20) B->C and C->D = 1 - (1 - p)^adult_prog_hazard_mult (effective values in $natural_history$derived); A->B also drives adult prior-A -> B unscaled; 06 runs the full $natural_history through abcd_nh_build",
    population_backbone_approximation = paste(
      "Population is EXOGENOUS (WPP2024 by age x sex x year); every stock ages by cohort; each cycle the No-RHD",
      "population (population - classifiable RHD - prior-A) is reconciled with the carried never-onset + normalised",
      "survivors (births at age 0; backbone residual shared proportionally); disease and prior-A stocks are never",
      "rescaled; living disease above population is surfaced as a violation, never clipped silently."),
    eff_cf_role = "progression mode: no direct SAP mortality effect, eff_cf == 1 by construction (legacy column kept)",
    rhd_death_lo  = getp("RHD_DEATH_LO",  1e3), rhd_death_hi  = getp("RHD_DEATH_HI",  1e5),
    allc_death_lo = getp("ALLC_DEATH_LO", 5e5), allc_death_hi = getp("ALLC_DEATH_HI", 4e6),
    intervention_labels = c(ref = "reference_care_held_at_baseline",
                            sap = "echo_screening_diagnosis_BPG_scale_up"),
    built_from = c(basename(IN_POP), basename(IN_DISEASE), basename(IN_CALIB))
  )
)

saveRDS(baseline_state, file = OUT_FILE)

message("── Saved ──────────────────────────────────────")
message(sprintf("  %s  | locations: %s | horizon %d-%d | scenarios: %s | engine %s",
                basename(OUT_FILE), paste(LOCATIONS, collapse = ", "),
                min(years), max(years), paste(SCENARIOS, collapse = ", "), ABCD_ENGINE_SPEC))
message("── 05_build_baseline.R complete ───────────────────────")
message("  Next: 06_run_prevention_model.R")
