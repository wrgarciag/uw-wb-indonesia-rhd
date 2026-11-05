# ==============================================================================
# RHD secondary-prevention investment case: MODEL RUNNER (A/B/C/D)
# scripts/06_run_prevention_model.R
#
# Default natural history (engine abcd_first_diagnosis_v1, R/abcd_engine.R): WHF
# stages A/B/C/D held as undiagnosed / confirmed-diagnosis (off / on BPG) stocks,
# adult prior-A and normalised history (no recurrence), attained-age rules, D-only
# RHD mortality. Intervention package: one echo screen across ages 5-39 among people
# without a prior confirmed diagnosis + confirmation + BPG for eligible confirmed
# disease (A at age <= 20, B); SAP lowers PROGRESSION hazards only.
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT (the actual model runner)
# ------------------------------------------------------------------------------
# Runs the reference and scale-up arms from the analysis-start state assembled by 05
# (the 04b bridged state: no reseeding). Contains NO monetary values (08 only).
#
#   INPUT : data/baseline_state.rds                      (from 05)
#   OUTPUT: output/out_model/<location>.rds              (one RDS per location)
#
# THE ENGINE is the shared abcd_project -> abcd_one_cycle -> abcd_age_state path, the
# SAME code and time convention 04b used for burn-in, calibration and bridge:
#   start state -> care update -> competing events -> new disease entries at year end ->
#   ageing (+ A -> adult prior-A reclassification at 20 -> 21).
# Only the care cascade (screening target, initiation, retention, adherence) differs
# between arms; surgery coverage/effects are held equal (background service).
#
# OUTPUT TABLES PER LOCATION
#   (1) $wsd    — WELL-SICK-DEAD aggregate (ref + sap, `scenario` col). sick = total
#                 CLASSIFIABLE RHD A+B+C+D, undiagnosed + diagnosed (end-of-cycle stock);
#                 well = population - sick (never-onset, normalised and adult prior-A);
#                 newcases = new modelled classifiable disease (NOT diagnoses); dead =
#                 RHD deaths; eff_ir = 1; eff_cf = 1 (no direct effect).
#   (2) $stages — aggregate A/B/C/D stocks and flows (new disease by entry stage, A
#                 resolution, adult B normalisation, prior-A re-entry, A -> prior-A
#                 reclassification, per-stage deaths), undiagnosed / diagnosed / on-BPG /
#                 prior-A / normalised / never-onset stocks, care volumes (screening
#                 encounters, true positives, confirmations, first confirmed diagnoses by
#                 route and stage, on-BPG persons, protection-weighted volume), surgery
#                 trace and care trajectories. A+B+C+D reproduces $wsd `sick` exactly.
#
# PARALLELISM: by LOCATION (foreach/doParallel). Workers receive every engine helper
# (ABCD_ENGINE_FUNS) and call setDTthreads(1).
# ==============================================================================

library(data.table)
library(foreach)

if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")
if (!exists("wd_outp")) wd_outp <- paste0(here::here("output"), "/")

# Shared engine (single source of truth; also used by 03/04b/05) — always (re)source
eng <- if (exists("wd")) paste0(wd, "R/abcd_engine.R") else here::here("R", "abcd_engine.R")
source(eng)

RUN_PAR   <- if (exists("run_model_par")) isTRUE(run_model_par) else TRUE
MAX_CORES <- if (exists("MAX_CORES")) MAX_CORES else 14L

IN_FILE  <- paste0(wd_data, "baseline_state.rds")
OUT_DIR  <- paste0(wd_outp, "out_model/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(IN_FILE))
  stop("Missing baseline state:\n  ", IN_FILE, "\n  Run 05_build_baseline.R first.", call. = FALSE)

baseline_state <- readRDS(IN_FILE)
LOCATIONS <- baseline_state$locations
SCENARIOS <- baseline_state$meta$scenarios
if (!identical(baseline_state$meta$engine_spec, ABCD_ENGINE_SPEC))
  stop("baseline_state was built for engine '", baseline_state$meta$engine_spec, "' but R/abcd_engine.R is '",
       ABCD_ENGINE_SPEC, "'. Re-run 03 -> 05.", call. = FALSE)

if (exists("ANALYSIS_YEARS", inherits = TRUE)) {
  ay <- sort(as.integer(get("ANALYSIS_YEARS", inherits = TRUE)))
  my <- sort(as.integer(baseline_state$meta$years))
  if (!identical(my, ay))
    stop("baseline_state horizon (", min(my), "-", max(my), ") != ANALYSIS_YEARS (",
         min(ay), "-", max(ay), "). Re-run 03-05 under the current window.", call. = FALSE)
}

message(sprintf("── 06_run_prevention_model.R (%s) : %d location(s), scenarios %s ──",
                ABCD_ENGINE_SPEC, length(LOCATIONS), paste(SCENARIOS, collapse = ", ")))

# ==============================================================================
# ENGINE — self-contained per-location runner (safe to ship to a parallel worker)
# ==============================================================================
run_location <- function(loc, baseline_state) {
  setDTthreads(1)
  meta  <- baseline_state$meta
  st    <- baseline_state$states[[loc]]
  AGES  <- meta$AGES; SEXES <- meta$SEXES; years <- meta$years
  n_age <- length(AGES); n_sex <- length(SEXES); n_years <- length(years)
  SCEN  <- meta$scenarios
  ilab  <- meta$intervention_labels
  CAUSE <- "Rheumatic heart disease"

  nhe  <- abcd_nh_build(st$natural_history, AGES, SEXES)    # same builder as 04b
  eff  <- st$effects; surg <- st$surgery; oth <- st$oth_mort
  care_vec <- c("screen_window_prob", "clinical_detection", "bpg_initiation", "bpg_retention",
                "adherence_good", "adherence_partial", "surgery")
  win_lo <- st$coverage$constants$screen_age_lo; win_hi <- st$coverage$constants$screen_age_hi
  y_last <- nhe$boundary_age; p_Bn <- st$natural_history$p_B_to_no_rhd_adult   # y_last = last child-rule age (20)

  melt_year <- function(iy, mats) {
    dt <- data.table(age = rep(AGES, times = n_sex), sex = rep(SEXES, each = n_age), year = years[iy])
    for (nm in names(mats)) set(dt, j = nm, value = as.vector(mats[[nm]]))
    dt
  }

  run_one <- function(scenario) {
    covs <- st$coverage[[scenario]]
    care_year <- abcd_care_schedule(c(covs[care_vec], st$coverage$constants), n_years, AGES, SEXES)
    proj <- abcd_project(init = st$seed_state, pop = st$pop, ir = st$ir, oth = oth, nhe = nhe,
                         eff = eff, surg = surg, care_year = care_year, years = years, light = FALSE)
    wsd_l <- vector("list", n_years); stg_l <- vector("list", n_years)
    for (iy in seq_len(n_years)) {
      popm <- st$pop[, , iy]; cyc <- proj$cycles[[iy]]; s0 <- proj$start[[iy]]; s1 <- cyc$state
      wsd_l[[iy]] <- melt_year(iy, list(
        well = cyc$well_end, sick = cyc$sick_end, newcases = cyc$new_rhd_total,
        dead = cyc$rhd_deaths, pop = popm, all.mx = cyc$all_mx))[
          , `:=`(eff_ir = cyc$eff_ir, eff_cf = cyc$sap_mult)]
      rc <- proj$reclass[[iy]]
      stg_l[[iy]] <- melt_year(iy, list(
        A = cyc$A, B = cyc$B, C = cyc$C, D = cyc$D, pop = popm,
        living_rhd_start = cyc$rhd_start, B_start = abcd_stage_total(s0, "B"),
        # NEW MODELLED DISEASE (first classifiable entry; never diagnoses) by entry stage
        new_rhd_A = cyc$new_rhd_A, new_rhd_B = cyc$new_rhd_B, new_rhd_C = cyc$new_rhd_C, new_rhd_D = cyc$new_rhd_D,
        new_rhd_total = cyc$new_rhd_total,
        A_to_no_rhd = cyc$A_to_no_rhd, B_to_no_rhd = cyc$B_to_no_rhd,
        A_to_B = cyc$A_to_B, B_to_A = cyc$B_to_A, B_to_C = cyc$B_to_C,
        C_to_B = cyc$C_to_B, C_to_D = cyc$C_to_D, D_to_C = cyc$D_to_C,
        # adult prior-A: re-entry into classifiable B (not first entry, not first diagnosis) and the
        # A -> prior-A reclassification at the 20 -> 21 ageing step (people conserved; recorded at age 20)
        reentry_prior_A_to_B = cyc$reentry_prior_A_to_B,
        reentry_prior_A_to_B_diagnosed_history = cyc$reentry_prior_A_to_B_diagnosed_history,
        reclass_A_to_prior_A = rc$reclass_undiagnosed + rc$reclass_diagnosed, reclass_A_to_prior_A_on_bpg = rc$reclass_on_bpg,
        rhd_deaths_A = cyc$rhd_deaths_A, rhd_deaths_B = cyc$rhd_deaths_B,
        rhd_deaths_C = cyc$rhd_deaths_C, rhd_deaths_D = cyc$rhd_deaths_D,
        rhd_deaths_D_undiagnosed = cyc$rhd_deaths_D_undiagnosed, rhd_deaths_D_diagnosed = cyc$rhd_deaths_D_diagnosed,
        other_deaths_A = cyc$other_deaths_A, other_deaths_B = cyc$other_deaths_B,
        other_deaths_C = cyc$other_deaths_C, other_deaths_D = cyc$other_deaths_D,
        rhd_deaths = cyc$rhd_deaths, other_deaths = cyc$other_deaths, rhd_exits_total = cyc$rhd_exits_total,
        # programme volumes (costed once in 08; confirmation is uncosted)
        n_screened = cyc$n_screened, n_screen_true_positive = cyc$n_screen_true_positive,
        n_confirm_completed = cyc$n_confirm_completed,
        n_diagnosed = cyc$n_diagnosed, n_undiagnosed = cyc$n_undiagnosed, n_new_diagnoses = cyc$n_new_diagnoses,
        n_screen_detected = cyc$n_screen_detected, n_dx_screen = cyc$n_dx_screen, n_dx_clinical = cyc$n_dx_clinical,
        n_dx_A = cyc$n_dx_A, n_dx_B = cyc$n_dx_B, n_dx_C = cyc$n_dx_C, n_dx_D = cyc$n_dx_D,
        n_bpg_eligible_diagnosed = cyc$n_bpg_eligible_diagnosed, n_bpg_initiations = cyc$n_bpg_initiations,
        n_bpg_dropouts = cyc$n_bpg_dropouts, n_on_bpg_A = cyc$n_on_bpg_A, n_on_bpg_B = cyc$n_on_bpg_B,
        n_on_optimal_treatment = cyc$n_on_optimal_treatment, n_effectively_protected = cyc$n_effectively_protected,
        C_requiring_surgery = cyc$C_requiring_surgery, D_requiring_surgery = cyc$D_requiring_surgery,
        surgeries_C = cyc$surgeries_C, surgeries_D = cyc$surgeries_D, total_surgeries = cyc$total_surgeries,
        # stocks + diagnostics
        D_start = cyc$D_start,
        well_deaths = cyc$well_deaths, never_onset_deaths = cyc$never_onset_deaths, recovered_deaths = cyc$recovered_deaths,
        prior_A_deaths = cyc$prior_A_deaths,
        never_onset_start = cyc$never_onset_start, recovered_start = cyc$recovered_start, prior_A_start = cyc$prior_A_start,
        never_onset_end = cyc$never_onset_end, normalised_end = cyc$normalised_end, prior_A_end = cyc$prior_A_end,
        undiagnosed_end = cyc$undiagnosed_end, diagnosed_end = cyc$diagnosed_end, on_bpg_end = cyc$on_bpg_end,
        never_screened_share = cyc$never_screened_share, screen_annual_prob = cyc$screen_annual_prob,
        diagnosed_share = cyc$diagnosed_share, on_bpg_share = cyc$on_bpg_share,
        sap_effective_protection = cyc$sap_c_eff,
        sap_prog_mult_A_to_B = cyc$sap_prog_mult_A_to_B, sap_prog_mult_B_to_C = cyc$sap_prog_mult_B_to_C,
        sap_prog_mult_C_to_D = cyc$sap_prog_mult_C_to_D + popm * 0,
        demog_adjust_never_onset = cyc$demog_adjust_never_onset, demog_adjust_recovered = cyc$demog_adjust_recovered))[
          , `:=`(screen_window_prob = covs$screen_window_prob[iy], clinical_detection = covs$clinical_detection[iy],
                 bpg_initiation = covs$bpg_initiation[iy], bpg_retention = covs$bpg_retention[iy],
                 adherence_good = covs$adherence_good[iy], adherence_partial = covs$adherence_partial[iy],
                 sap_adherence_factor = cyc$adherence_factor,
                 surgery_coverage = covs$surgery[iy],
                 effective_surgery_reach_C = cyc$reach_C, effective_surgery_reach_D = cyc$reach_D,
                 eff_cf = cyc$sap_mult,
                 eff_C_to_D_surgery = cyc$eff_C_to_D_surg, eff_D_to_rhd_death_surgery = cyc$eff_D_death_surg,
                 d_start_zero_flag = as.numeric(cyc$D_start <= 0))]
    }
    wsd <- rbindlist(wsd_l); stg <- rbindlist(stg_l)
    wsd[, `:=`(scenario = scenario, location = loc, cause = CAUSE, intervention = unname(ilab[[scenario]]))]
    stg[, `:=`(scenario = scenario, location = loc, cause = CAUSE, intervention = unname(ilab[[scenario]]))]
    list(wsd = wsd, stages = stg, pop_deficit = proj$pop_deficit_total, care_violation = proj$care_violation_total,
         onset_cap_hits = proj$onset_cap_hits, demographic = cbind(scenario = scenario, proj$demographic_residual))
  }

  runs <- lapply(SCEN, run_one)
  wsd_all <- rbindlist(lapply(runs, `[[`, "wsd"))
  stg_all <- rbindlist(lapply(runs, `[[`, "stages"))
  pop_deficit_total <- sum(vapply(runs, `[[`, numeric(1), "pop_deficit"))
  care_violation_total <- sum(vapply(runs, `[[`, numeric(1), "care_violation"))
  onset_cap_total <- sum(vapply(runs, `[[`, numeric(1), "onset_cap_hits"))
  demog <- rbindlist(lapply(runs, `[[`, "demographic"))

  # -------------------- per-location SANITY CHECKS --------------------------
  fail <- function(cond, msg) if (isTRUE(cond)) stop("[", loc, "] ", msg, call. = FALSE)
  exp_rows <- length(SCEN) * n_years * n_age * n_sex
  fail(nrow(wsd_all) != exp_rows, "WSD table has an incomplete scenario/age/year grid.")
  fail(nrow(stg_all) != exp_rows, "stage table has an incomplete scenario/age/year grid.")
  fail(!identical(sort(unique(as.integer(wsd_all$year))), as.integer(years)), "WSD years do not match the analysis horizon.")
  fail(!identical(sort(unique(as.integer(stg_all$year))), as.integer(years)), "stage years do not match the analysis horizon.")
  num_wsd <- wsd_all[, .SD, .SDcols = is.numeric]
  fail(anyNA(num_wsd) || any(!is.finite(as.matrix(num_wsd))), "WSD table contains NA / non-finite values.")
  fail(any(vapply(num_wsd, function(x) any(x < -1e-6), logical(1))), "WSD table contains negative values.")
  signed <- c("demog_adjust_never_onset", "demog_adjust_recovered")      # backbone adjustments may be negative
  num_stg <- stg_all[, .SD, .SDcols = setdiff(names(stg_all)[vapply(stg_all, is.numeric, logical(1))], signed)]
  fail(anyNA(stg_all[, .SD, .SDcols = is.numeric]) || any(!is.finite(as.matrix(stg_all[, .SD, .SDcols = is.numeric]))),
       "stage table contains NA / non-finite values.")
  fail(any(vapply(num_stg, function(x) any(x < -1e-6), logical(1))), "stage table contains negative stocks/flows/volumes.")

  # WSD identities + population backbone
  fail(wsd_all[, any(well + sick > pop + 1e-3)], "well + sick exceeds pop somewhere.")
  fail(wsd_all[, any(all.mx + 1e-6 < dead)],     "all-cause deaths < RHD deaths somewhere.")
  fail(wsd_all[, any(eff_ir < 0 | eff_ir > 1 | eff_cf < 0 | eff_cf > 1)], "eff_ir / eff_cf outside [0,1].")
  fail(pop_deficit_total > 1e-6, sprintf("living RHD stock exceeded population (total deficit %.3g).", pop_deficit_total))
  fail(care_violation_total > 1e-6, sprintf("invalid pooled care state (violation %.3g).", care_violation_total))
  fail(onset_cap_total > 0, "first-onset probability hit the engine cap.")
  fail(stg_all[, max(abs(never_onset_start + living_rhd_start + recovered_start + prior_A_start - pop))] > 1e-4,
       "never-onset + classifiable RHD + normalised + prior-A != population at the cycle start.")
  # within-cycle conservation: start = end alive + all deaths
  alive <- merge(stg_all[, .(scenario, sex, age, year, start = never_onset_start + living_rhd_start + recovered_start + prior_A_start,
                             end = never_onset_end + A + B + C + D + normalised_end + prior_A_end)],
                 wsd_all[, .(scenario, sex, age, year, all.mx)], by = c("scenario", "sex", "age", "year"))
  fail(alive[, max(abs(start - end - all.mx))] > 1e-4, "within-cycle population conservation violated.")
  chk_mx <- merge(stg_all[, .(scenario, sex, age, year, recon = rhd_deaths + other_deaths + well_deaths + prior_A_deaths)],
                  wsd_all[, .(scenario, sex, age, year, all.mx)], by = c("scenario", "sex", "age", "year"))
  fail(chk_mx[, max(abs(recon - all.mx))] > 1e-6, "all.mx != RHD deaths + other deaths (classifiable) + never-onset/normalised + prior-A deaths.")
  # classifiable stock-flow identity (before ageing; reclassification to prior-A happens at ageing)
  fail(stg_all[, max(abs((A + B + C + D) - (living_rhd_start + new_rhd_total + reentry_prior_A_to_B - A_to_no_rhd - B_to_no_rhd -
                                             rhd_deaths - other_deaths)))] > 1e-6,
       "classifiable stock-flow identity violated (end = start + new disease + prior-A re-entry - normalisation - RHD deaths - other deaths).")
  fail(stg_all[, max(abs(new_rhd_total - (new_rhd_A + new_rhd_B + new_rhd_C + new_rhd_D)))] > 1e-6, "entry-stage new disease does not sum to total.")
  nc <- merge(wsd_all[, .(scenario, sex, age, year, newcases)],
              stg_all[, .(scenario, sex, age, year, new_rhd_total)], by = c("scenario", "sex", "age", "year"))
  fail(nc[, max(abs(newcases - new_rhd_total))] > 1e-9, "$wsd newcases != total new modelled disease.")
  # attained-age structure
  fail(stg_all[age > y_last, any(A != 0 | A_to_no_rhd != 0 | B_to_A != 0 | C_to_B != 0 | n_on_bpg_A != 0 | new_rhd_A != 0)],
       "classifiable A, A on BPG, new A or childhood regression flows at adult ages.")
  fail(stg_all[age <= y_last, any(prior_A_end != 0 | reentry_prior_A_to_B != 0 | B_to_no_rhd != 0 | new_rhd_B != 0)],
       "prior-A, prior-A re-entry, B normalisation or new B entry at child ages.")
  fail(stg_all[, any(new_rhd_C != 0 | new_rhd_D != 0)], "new disease entered C/D (excluded by the entry mixes).")
  fail(stg_all[age != y_last, any(reclass_A_to_prior_A != 0)], "A -> prior-A reclassification outside the boundary age.")
  fail(stg_all[, any(B_to_no_rhd > B_start * p_Bn * (1 + 1e-9) + 1e-9)], "realised adult B normalisation exceeds its isolated input risk.")
  # D-only + entry timing: deaths in D never exceed the D START stock x the D input risk
  fail(stg_all[, any(rhd_deaths_A != 0 | rhd_deaths_B != 0 | rhd_deaths_C != 0)], "D_only: non-zero RHD deaths in A/B/C.")
  fail(stg_all[, max(abs(rhd_deaths - rhd_deaths_D_undiagnosed - rhd_deaths_D_diagnosed))] > 1e-9, "D deaths by diagnosis status do not sum to RHD deaths.")
  fail(stg_all[, any(rhd_deaths_D > D_start * nhe$pd * (1 + 1e-9) + 1e-9)],
       "D deaths exceed start-of-cycle D x input risk (a new D case may not die in its entry cycle).")
  fail(stg_all[, sum(other_deaths_A + other_deaths_B + other_deaths_C)] <= 0, "background deaths absent from A/B/C.")
  fail(stg_all[d_start_zero_flag == 1 & rhd_deaths_D > 0, .N] > 0, "RHD deaths from D with zero D start stock.")
  # SAP progression mechanism
  fail(stg_all[, any(abs(eff_cf - 1) > 1e-12)], "progression mode: direct SAP mortality multiplier must be 1.")
  fail(stg_all[, any(sap_prog_mult_A_to_B > 1 + 1e-12 | sap_prog_mult_B_to_C > 1 + 1e-12 | sap_prog_mult_C_to_D > 1 + 1e-12)],
       "progression multipliers exceed 1.")
  # care nesting: on BPG <= eligible diagnosed <= diagnosed <= classifiable; protection <= on-BPG;
  # screen diagnoses <= confirmations <= true positives; first diagnoses split by route and stage
  fail(stg_all[, any(n_on_optimal_treatment > n_bpg_eligible_diagnosed * (1 + 1e-9) + 1e-9 |
                     n_bpg_eligible_diagnosed > n_diagnosed * (1 + 1e-9) + 1e-9 | n_diagnosed > living_rhd_start * (1 + 1e-9) + 1e-9)],
       "care state violates on BPG <= eligible diagnosed <= diagnosed <= classifiable RHD.")
  fail(stg_all[, max(abs(n_diagnosed + n_undiagnosed - living_rhd_start))] > 1e-6, "diagnosed + undiagnosed != classifiable RHD.")
  fail(stg_all[, any(n_dx_screen > n_confirm_completed * (1 + 1e-9) + 1e-9 | n_confirm_completed > n_screen_true_positive * (1 + 1e-9) + 1e-9)],
       "screen diagnoses > confirmations > true positives ordering violated.")
  fail(stg_all[, max(abs(n_new_diagnoses - n_dx_screen - n_dx_clinical), abs(n_new_diagnoses - n_dx_A - n_dx_B - n_dx_C - n_dx_D))] > 1e-6,
       "first diagnoses do not split exactly by route and by stage.")
  fail(stg_all[, any(n_effectively_protected > n_on_optimal_treatment * (1 + 1e-12) + 1e-12)], "protected volume exceeds on-BPG volume.")
  fail(stg_all[, any(never_screened_share < -1e-12 | never_screened_share > 1 + 1e-12)], "never-screened share outside [0,1].")
  fail(stg_all[age < win_lo | age > win_hi, any(n_screened != 0)], "screening encounters outside the screening age window.")
  fail(stg_all[, any(n_screened > pop * screen_annual_prob * (1 + 1e-9) + 1e-9)],
       "screening encounters exceed population x annual screening probability.")
  fail(stg_all[, max(abs(n_on_optimal_treatment - living_rhd_start * on_bpg_share))] > 1e-6, "on-BPG volume != living RHD x on-BPG share.")
  fail(stg_all[, max(abs(n_on_optimal_treatment - n_on_bpg_A - n_on_bpg_B))] > 1e-6, "on-BPG volume != on BPG in A + B.")
  # surgery is a service, not a stock
  fail(stg_all[, any(surgeries_C > C_requiring_surgery + 1e-6)], "surgeries_C exceeds C_requiring_surgery.")
  fail(stg_all[, any(surgeries_D > D_requiring_surgery + 1e-6)], "surgeries_D exceeds D_requiring_surgery.")
  fail(stg_all[, any(total_surgeries > C_requiring_surgery + D_requiring_surgery + 1e-6)], "total_surgeries exceeds total requiring surgery.")
  # stages reconstruct WSD sick exactly
  chk <- merge(stg_all[, .(sick_stg = A + B + C + D, sick_care = undiagnosed_end + diagnosed_end), by = .(scenario, sex, age, year)],
               wsd_all[, .(scenario, sex, age, year, sick)], by = c("scenario", "sex", "age", "year"))
  fail(chk[, max(abs(sick_stg - sick), abs(sick_care - sick))] > 1e-6, "stage / care-state stocks do not reconstruct WSD sick.")

  # intervention effect direction: scale-up averts RHD deaths vs reference (cumulative)
  d_ref <- wsd_all[scenario == "ref", sum(dead)]; d_sap <- wsd_all[scenario == "sap", sum(dead)]
  fail(d_sap > d_ref + 1e-6, "SAP scenario has MORE cumulative RHD deaths than reference.")

  # order-of-magnitude anchors at the base year (reference arm)
  by1 <- min(years)
  base <- wsd_all[scenario == "ref" & year == by1]
  rhd_d <- base[, sum(dead)]; allc_d <- base[, sum(all.mx)]
  rhd_lo  <- if (is.null(meta$rhd_death_lo))  1e3 else meta$rhd_death_lo
  rhd_hi  <- if (is.null(meta$rhd_death_hi))  5e4 else meta$rhd_death_hi
  allc_lo <- if (is.null(meta$allc_death_lo)) 5e5 else meta$allc_death_lo
  allc_hi <- if (is.null(meta$allc_death_hi)) 4e6 else meta$allc_death_hi
  fail(rhd_d  < rhd_lo  || rhd_d  > rhd_hi,  sprintf("base-year RHD deaths %.0f outside band %g-%g.", rhd_d, rhd_lo, rhd_hi))
  fail(allc_d < allc_lo || allc_d > allc_hi, sprintf("base-year all-cause deaths %.0f outside band %g-%g.", allc_d, allc_lo, allc_hi))
  attr(wsd_all, "deaths_averted_cum") <- d_ref - d_sap

  # --- additive diagnostics ($diag) ---------------------------------------------
  d_fat <- rbind(
    stg_all[, .(ages = "all", D_start = sum(D_start), rhd_deaths_D = sum(rhd_deaths_D)), by = .(scenario, year)],
    stg_all[age <= y_last, .(ages = sprintf("0-%d", y_last), D_start = sum(D_start), rhd_deaths_D = sum(rhd_deaths_D)), by = .(scenario, year)],
    stg_all[age > y_last, .(ages = sprintf("%d+", y_last + 1L), D_start = sum(D_start), rhd_deaths_D = sum(rhd_deaths_D)), by = .(scenario, year)])
  d_fat[, `:=`(input_D_death_prob = nhe$pd, realized_D_death_risk = fifelse(D_start > 0, rhd_deaths_D / D_start, NA_real_), D_start_zero = D_start <= 0)]
  sap_tr <- stg_all[, .(living_rhd_start = sum(living_rhd_start), diagnosed = sum(n_diagnosed), new_diagnoses = sum(n_new_diagnoses),
                        new_diagnoses_screen = sum(n_dx_screen), new_diagnoses_clinical = sum(n_dx_clinical),
                        screen_true_positives = sum(n_screen_true_positive), confirmations_completed = sum(n_confirm_completed),
                        bpg_eligible_diagnosed = sum(n_bpg_eligible_diagnosed),
                        on_bpg_person_years = sum(n_on_optimal_treatment), on_bpg_A = sum(n_on_bpg_A), on_bpg_B = sum(n_on_bpg_B),
                        on_bpg_B_adult = sum(n_on_bpg_B[age > y_last]),
                        effectively_protected_person_years = sum(n_effectively_protected),
                        adherence_factor = sap_adherence_factor[1], bpg_initiation = bpg_initiation[1], bpg_retention = bpg_retention[1],
                        eff_cf = mean(eff_cf)), by = .(scenario, year)]
  sap_tr[, `:=`(diagnosed_share = diagnosed / living_rhd_start, on_bpg_share = on_bpg_person_years / living_rhd_start,
                effective_protection = effectively_protected_person_years / living_rhd_start)]
  scr <- stg_all[age >= win_lo & age <= win_hi, .(screen_window_prob = screen_window_prob[1], annual_prob = screen_annual_prob[1],
                                                  pop_in_window = sum(pop), encounters = sum(n_screened),
                                                  ever_screened_share_in_window = 1 - sum(pop * never_screened_share) / sum(pop)),
                 by = .(scenario, year)]
  state_tr <- stg_all[, .(classifiable = sum(A + B + C + D), undiagnosed = sum(undiagnosed_end), diagnosed = sum(diagnosed_end),
                          on_bpg = sum(on_bpg_end), D = sum(D), prior_A = sum(prior_A_end), normalised = sum(normalised_end),
                          new_disease = sum(new_rhd_total), reentry_prior_A_to_B = sum(reentry_prior_A_to_B),
                          reclass_A_to_prior_A = sum(reclass_A_to_prior_A), A_resolution = sum(A_to_no_rhd), B_normalisation = sum(B_to_no_rhd),
                          first_diagnoses = sum(n_new_diagnoses), rhd_deaths = sum(rhd_deaths), rhd_deaths_0_19 = sum(rhd_deaths[age < 20]),
                          rhd_deaths_D_undiagnosed = sum(rhd_deaths_D_undiagnosed)), by = .(scenario, year)]
  mx_recon <- chk_mx[, .(all_mx = sum(all.mx), reconstructed = sum(recon)), by = .(scenario, year)]

  list(wsd = wsd_all, stages = stg_all,
       diag = list(deaths_averted_cum = d_ref - d_sap,
                   rhd_deaths_base = rhd_d, allcause_deaths_base = allc_d,
                   engine_spec = meta$engine_spec,
                   mortality_structure = meta$mortality_structure, sap_mechanism = eff$sap_mechanism,
                   d_fatality_by_year = d_fat,
                   d_fatality_convention = "deaths_D / D_start (start-of-cycle stock), all ages and by child/adult rule ages; NA with D_start_zero = TRUE when D_start = 0",
                   sap_trace_by_year = sap_tr, screening_by_year = scr, care_state_by_year = state_tr,
                   uncosted = "confirmatory echo (n_confirm_completed) has no unit cost in 00 block J",
                   all_cause_reconstruction_by_year = mx_recon,
                   demographic_adjustment_by_year = demog,
                   pop_deficit_total = pop_deficit_total, care_violation_total = care_violation_total))
}

# ==============================================================================
# RUN — parallel by location (fall back to sequential for a single location)
# ==============================================================================
use_par <- RUN_PAR && length(LOCATIONS) > 1L
if (use_par) {
  suppressPackageStartupMessages({ library(doParallel); library(parallel) })
  n_cores <- max(1L, min(MAX_CORES, length(LOCATIONS), parallel::detectCores() - 1L))
  message(sprintf("  Running %d locations in parallel on %d cores...", length(LOCATIONS), n_cores))
  cl <- makeCluster(n_cores); registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  results <- foreach(loc = LOCATIONS, .packages = "data.table",
                     .export = c("run_location", "baseline_state", ABCD_ENGINE_FUNS)) %dopar% {
    run_location(loc, baseline_state)
  }
} else {
  message("  Running sequentially (single location or parallel disabled)...")
  results <- foreach(loc = LOCATIONS) %do% run_location(loc, baseline_state)
}
names(results) <- LOCATIONS

# ==============================================================================
# WRITE one RDS per location + report
# ==============================================================================
for (loc in LOCATIONS) {
  res <- results[[loc]]
  out <- list(wsd = res$wsd, stages = res$stages, diag = res$diag, meta = baseline_state$meta)
  fn  <- paste0(OUT_DIR, gsub("[^A-Za-z0-9]+", "_", loc), ".rds")
  saveRDS(out, file = fn)
  message(sprintf("  %-12s -> %s | RHD deaths averted (cum, %d-%d) = %s | base-yr RHD deaths = %s",
                  loc, basename(fn), min(baseline_state$meta$years), max(baseline_state$meta$years),
                  formatC(round(res$diag$deaths_averted_cum), format = "d", big.mark = ","),
                  formatC(round(res$diag$rhd_deaths_base),    format = "d", big.mark = ",")))
}

message("── 06_run_prevention_model.R complete ─────────────────")
message(sprintf("  Wrote %d file(s) to %s ($wsd aggregate + $stages A/B/C/D + care-state detail + care volumes)",
                length(LOCATIONS), OUT_DIR))
message("  Next: 07_make_outputs.R (tidy long tables), 08_economic_evaluation.R (economics)")
