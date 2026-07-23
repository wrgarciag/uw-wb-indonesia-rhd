# ==============================================================================
# RHD secondary-prevention investment case: MODEL RUNNER (A/B/C/D)
# scripts/06_run_prevention_model.R
#
# Disease structure: WHF RHD stages A/B/C/D.
#   No RHD -> A <-> B <-> C -> D -> RHD death, competing other-cause death from
#   every living stage. Incident RHD enters stage A (no separate ARF state).
# Intervention package: echo screening + diagnosis + secondary antibiotic
#   prophylaxis (SAP / optimal treatment), which cuts RHD-specific MORTALITY.
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT (the actual model runner)
# ------------------------------------------------------------------------------
# Runs the reference and SAP scale-up scenarios starting from the initial state
# assembled by 05. It does NOT recompute inputs and contains NO monetary values
# (economics are in 08 only).
#
#   INPUT : data/baseline_state.rds                      (from 05)
#   OUTPUT: output/out_model/<location>.rds              (one RDS per location)
#
# THE ENGINE  (matrix form: [age x sex] arrays; Markov cycles as matrix ops)
#   Four LIVING stocks A, B, C, D advance annually by elementwise matrix
#   arithmetic (zero_mat/age_shift/melt_year, extended from three to four stocks):
#     * new incident RHD enters stage A via  new_rhd_A = no_rhd x IR (CALIBRATED IR);
#     * SAP reduces the RHD-specific DEATH probability of EVERY stage by
#         (1 - sap_rrr_rhd_death x effective_treatment_coverage);
#     * SURGERY is a clinical SERVICE (not a state): a fraction of the C/D stock
#       is operated each cycle; its only epidemiological effect is a risk
#       reduction on C->D progression and D->RHD-death (population-average form).
#       Surgery NEVER removes people from C or D;
#     * competing (other-cause) background mortality is data-fed and held fixed,
#       the SAME age x sex risk for every stage;
#     * cohorts age one year per cycle (age_shift); the terminal age is open 100+.
#   Only the care-cascade coverage differs between the reference and SAP arms
#   (surgery coverage/effects are held equal in both -> surgery is a background
#   cost, not a driver of incremental cost).
#
#   Population-average surgery formulation (documented):
#     effective_surgery_reach_C = frac_C_requiring_surgery x surgery_coverage
#     effective_surgery_reach_D = frac_D_requiring_surgery x surgery_coverage
#     p_C_to_D_eff       = p_C_to_D      x (1 - eff_surgery_C_to_D        x reach_C)
#     p_D_to_rhd_death_eff = p_rhd_death_D x sap_mult x (1 - eff_surgery_D_to_rhd_death x reach_D)
#
# TWO OUTPUT TABLES PER LOCATION
#   (1) $wsd    — WELL-SICK-DEAD aggregate (ref + sap row-bound, `scenario` col).
#                 sick = A + B + C + D. Consumed by 07 for the standard long table.
#                 eff_ir = 1 (no incidence effect); eff_cf = SAP RHD-mortality
#                 multiplier = 1 - sap_rrr_rhd_death x effective_treatment_coverage.
#   (2) $stages — the A/B/C/D stock-and-flow table: four stocks, every transition
#                 flow, incident inflow, per-stage RHD & other-cause deaths, the
#                 stated + effective cascade coverages, screening/diagnosis/
#                 treatment volumes, and the surgery TRACE (C/D requiring surgery,
#                 surgeries delivered to C/D, total surgeries, effective reach).
#                 Collapsing A+B+C+D reproduces $wsd `sick` exactly.
#
# PARALLELISM: by LOCATION (foreach/doParallel). Workers call setDTthreads(1).
# ==============================================================================

library(data.table)
library(foreach)

if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")
if (!exists("wd_outp")) wd_outp <- paste0(here::here("output"), "/")

# Shared one-cycle A/B/C/D engine (single source of truth; also used by 04b).
if (!exists("abcd_one_cycle")) {
  eng <- if (exists("wd")) paste0(wd, "R/abcd_engine.R") else here::here("R", "abcd_engine.R")
  source(eng)
}

RUN_PAR   <- if (exists("run_model_par")) isTRUE(run_model_par) else TRUE
MAX_CORES <- if (exists("MAX_CORES")) MAX_CORES else 14L

IN_FILE  <- paste0(wd_data, "baseline_state.rds")
OUT_DIR  <- paste0(wd_outp, "out_model/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(IN_FILE))
  stop("Missing baseline state:\n  ", IN_FILE,
       "\n  Run 05_build_baseline.R first.", call. = FALSE)

baseline_state <- readRDS(IN_FILE)
LOCATIONS <- baseline_state$locations
SCENARIOS <- baseline_state$meta$scenarios

# analysis-period integrity: horizon inherited from 05 must equal ANALYSIS_YEARS
# when 00_run_all.R is driving (else a stale baseline_state under a different window).
if (exists("ANALYSIS_YEARS", inherits = TRUE)) {
  ay <- sort(as.integer(get("ANALYSIS_YEARS", inherits = TRUE)))
  my <- sort(as.integer(baseline_state$meta$years))
  if (!identical(my, ay))
    stop("baseline_state horizon (", min(my), "-", max(my), ") != ANALYSIS_YEARS (",
         min(ay), "-", max(ay), "). Re-run 03-05 under the current window.", call. = FALSE)
}

message(sprintf("── 06_run_prevention_model.R (A/B/C/D) : %d location(s), scenarios %s ──",
                length(LOCATIONS), paste(SCENARIOS, collapse = ", ")))

# ==============================================================================
# ENGINE — self-contained per-location runner (safe to ship to a parallel worker)
# ==============================================================================
run_location <- function(loc, baseline_state) {
  setDTthreads(1)                                   # safe inside parallel workers
  meta  <- baseline_state$meta
  st    <- baseline_state$states[[loc]]
  AGES  <- meta$AGES; SEXES <- meta$SEXES; years <- meta$years
  n_age <- length(AGES); n_sex <- length(SEXES); n_years <- length(years)
  SCEN  <- meta$scenarios
  ilab  <- meta$intervention_labels
  CAUSE <- "Rheumatic heart disease"                 # defined inside (worker-safe)

  tr  <- st$transitions; pd <- st$p_rhd_death; eff <- st$effects
  surg <- st$surgery;    oth <- st$oth_mort

  zero_mat  <- function() matrix(0, n_age, n_sex, dimnames = list(AGES, SEXES))
  age_shift <- function(M) {                         # advance age a -> a+1; 100 open
    N <- zero_mat()
    N[2:n_age, ] <- M[1:(n_age - 1), ]
    N[n_age,  ]  <- N[n_age, ] + M[n_age, ]          # accumulate terminal 100+ group
    N
  }
  # optional school-age screening mask (OFF by default; sensitivity only)
  restrict_screen <- isTRUE(st$coverage$screen_age_restrict)
  screen_mask <- zero_mat() + 1                      # default: whole population screened
  if (restrict_screen) {
    screen_mask <- zero_mat()
    screen_ages <- as.character(st$coverage$screen_age_lo:st$coverage$screen_age_hi)
    screen_mask[intersect(screen_ages, rownames(screen_mask)), ] <- 1
  }

  # melt a named list of [age x sex] matrices into a tidy DT (adds age/sex/year)
  melt_year <- function(iy, mats) {
    dt <- data.table(age = rep(AGES, times = n_sex),
                     sex = rep(SEXES, each = n_age),
                     year = years[iy])
    for (nm in names(mats)) set(dt, j = nm, value = as.vector(mats[[nm]]))
    dt
  }

  run_one <- function(scenario) {
    covs <- st$coverage[[scenario]]
    A <- B <- C <- D <- zero_mat()
    wsd_l <- vector("list", n_years); stg_l <- vector("list", n_years)

    for (iy in seq_len(n_years)) {
      popm <- st$pop[, , iy]
      irm  <- st$ir[, , iy]                          # CALIBRATED incidence, trended
      c_screen <- covs$screen[iy]
      c_dx     <- covs$eff_diagnosis[iy]             # effective diagnosis coverage
      c_tx     <- covs$eff_treatment[iy]             # effective (optimal) treatment coverage
      c_surg   <- covs$surgery[iy]

      if (iy == 1L) { A <- st$seed$A; B <- st$seed$B; C <- st$seed$C; D <- st$seed$D }

      # --- ONE annual A/B/C/D cycle via the shared hazard engine ---------------
      #  Single source of truth (R/abcd_engine.R), also used by the structural
      #  calibration (04b). Competing risks are combined on the HAZARD scale, so
      #  each stage "stay" is a proper residual, not a floored subtraction.
      cyc <- abcd_one_cycle(A, B, C, D, popm, irm, oth, tr, pd, eff, surg,
                            cov = list(treatment = c_tx, surgery = c_surg))

      rhd_start <- cyc$rhd_start
      eff_ir <- cyc$eff_ir; new_rhd_A <- cyc$new_rhd_A
      sap_mult <- cyc$sap_mult; reach_C <- cyc$reach_C; reach_D <- cyc$reach_D
      eff_C_to_D_surg <- cyc$eff_C_to_D_surg; eff_D_death_surg <- cyc$eff_D_death_surg
      A_next <- cyc$A; B_next <- cyc$B; C_next <- cyc$C; D_next <- cyc$D
      A_to_no_rhd <- cyc$A_to_no_rhd; A_to_B <- cyc$A_to_B
      B_to_A <- cyc$B_to_A; B_to_C <- cyc$B_to_C
      C_to_B <- cyc$C_to_B; C_to_D <- cyc$C_to_D; D_to_C <- cyc$D_to_C
      A_death_rhd <- cyc$rhd_deaths_A; B_death_rhd <- cyc$rhd_deaths_B
      C_death_rhd <- cyc$rhd_deaths_C; D_death_rhd <- cyc$rhd_deaths_D
      A_death_oth <- cyc$other_deaths_A; B_death_oth <- cyc$other_deaths_B
      C_death_oth <- cyc$other_deaths_C; D_death_oth <- cyc$other_deaths_D
      rhd_deaths <- cyc$rhd_deaths; other_deaths <- cyc$other_deaths; all_mx <- cyc$all_mx
      sick_end <- cyc$sick_end; well_end <- cyc$well_end
      C_req_surg <- cyc$C_requiring_surgery; D_req_surg <- cyc$D_requiring_surgery
      surgeries_C <- cyc$surgeries_C; surgeries_D <- cyc$surgeries_D
      total_surg <- cyc$total_surgeries

      # --- program volumes (counts; costed in 08) ------------------------------
      n_screened <- popm * screen_mask * c_screen                  # total pop (default)
      n_diagnosed <- rhd_start * c_dx                              # effective diagnosis
      n_on_optimal_treatment <- rhd_start * c_tx                   # effective treatment

      # --- (1) WSD aggregate row set -------------------------------------------
      wsd_l[[iy]] <- melt_year(iy, list(
        well = well_end, sick = sick_end, newcases = new_rhd_A,
        dead = rhd_deaths, pop = popm, all.mx = all_mx))[
          , `:=`(eff_ir = eff_ir, eff_cf = sap_mult)]

      # --- (2) A/B/C/D stock-and-flow row set -----------------------------------
      #  living_rhd_start = start-of-cycle prevalent RHD (the denominator for
      #  diagnosis/treatment volumes and for surgery need in THIS cycle).
      stg_l[[iy]] <- melt_year(iy, list(
        A = A_next, B = B_next, C = C_next, D = D_next, pop = popm,
        living_rhd_start = rhd_start,
        new_rhd_A = new_rhd_A,
        A_to_no_rhd = A_to_no_rhd, A_to_B = A_to_B,
        B_to_A = B_to_A, B_to_C = B_to_C,
        C_to_B = C_to_B, C_to_D = C_to_D,
        D_to_C = D_to_C,
        rhd_deaths_A = A_death_rhd, rhd_deaths_B = B_death_rhd,
        rhd_deaths_C = C_death_rhd, rhd_deaths_D = D_death_rhd,
        other_deaths_A = A_death_oth, other_deaths_B = B_death_oth,
        other_deaths_C = C_death_oth, other_deaths_D = D_death_oth,
        rhd_deaths = rhd_deaths, other_deaths = other_deaths,
        n_screened = n_screened, n_diagnosed = n_diagnosed,
        n_on_optimal_treatment = n_on_optimal_treatment,
        C_requiring_surgery = C_req_surg, D_requiring_surgery = D_req_surg,
        surgeries_C = surgeries_C, surgeries_D = surgeries_D,
        total_surgeries = total_surg))[
          , `:=`(screen_coverage = c_screen,
                 diagnosis_coverage = covs$diagnosis[iy],
                 effective_diagnosis_coverage = c_dx,
                 optimal_treatment_coverage = covs$treatment[iy],
                 effective_treatment_coverage = c_tx,
                 surgery_coverage = c_surg,
                 effective_surgery_reach_C = reach_C,
                 effective_surgery_reach_D = reach_D,
                 eff_cf = sap_mult,
                 eff_C_to_D_surgery = eff_C_to_D_surg,
                 eff_D_to_rhd_death_surgery = eff_D_death_surg)]

      # --- age the surviving stocks one year for the next cycle -----------------
      A <- age_shift(A_next); B <- age_shift(B_next)
      C <- age_shift(C_next); D <- age_shift(D_next)
    }

    wsd <- rbindlist(wsd_l); stg <- rbindlist(stg_l)
    wsd[, `:=`(scenario = scenario, location = loc, cause = CAUSE,
               intervention = unname(ilab[[scenario]]))]
    stg[, `:=`(scenario = scenario, location = loc, cause = CAUSE,
               intervention = unname(ilab[[scenario]]))]
    list(wsd = wsd, stages = stg)
  }

  runs <- lapply(SCEN, run_one)
  wsd_all <- rbindlist(lapply(runs, `[[`, "wsd"))
  stg_all <- rbindlist(lapply(runs, `[[`, "stages"))

  # -------------------- per-location SANITY CHECKS --------------------------
  fail <- function(cond, msg) if (isTRUE(cond)) stop("[", loc, "] ", msg, call. = FALSE)

  # completeness: scenarios x years x ages x sexes
  exp_rows <- length(SCEN) * n_years * n_age * n_sex
  fail(nrow(wsd_all) != exp_rows, "WSD table has an incomplete scenario/age/year grid.")
  fail(nrow(stg_all) != exp_rows, "stage table has an incomplete scenario/age/year grid.")

  # horizon integrity: the model output must cover EXACTLY the analysis years
  # (meta$years), no more, no fewer.
  fail(!identical(sort(unique(as.integer(wsd_all$year))), as.integer(years)),
       "WSD years do not match the analysis horizon (meta$years).")
  fail(!identical(sort(unique(as.integer(stg_all$year))), as.integer(years)),
       "stage years do not match the analysis horizon (meta$years).")

  num_wsd <- wsd_all[, .SD, .SDcols = is.numeric]
  fail(anyNA(num_wsd), "WSD table contains NA.")
  fail(any(vapply(num_wsd, function(x) any(x < -1e-6), logical(1))),
       "WSD table contains negative values.")
  num_stg <- stg_all[, .SD, .SDcols = is.numeric]
  fail(anyNA(num_stg), "stage table contains NA.")
  fail(any(vapply(num_stg, function(x) any(x < -1e-6), logical(1))),
       "stage table contains negative values.")

  # WSD identities: well + sick <= pop ; all.mx >= dead ; eff_* in [0,1]
  fail(wsd_all[, any(well + sick > pop + 1e-3)], "well + sick exceeds pop somewhere.")
  fail(wsd_all[, any(all.mx + 1e-6 < dead)],     "all-cause deaths < RHD deaths somewhere.")
  fail(wsd_all[, any(eff_ir < 0 | eff_ir > 1 | eff_cf < 0 | eff_cf > 1)],
       "eff_ir / eff_cf outside [0,1].")

  # stages must reconstruct WSD sick exactly (A+B+C+D == sick)
  chk <- merge(
    stg_all[, .(sick_stg = A + B + C + D), by = .(scenario, sex, age, year)],
    wsd_all[, .(scenario, sex, age, year, sick)],
    by = c("scenario", "sex", "age", "year"))
  fail(chk[, max(abs(sick_stg - sick))] > 1e-6,
       "stage A+B+C+D does not reconstruct WSD sick.")

  # surgery is a service, not a stock: volumes never exceed the number requiring
  fail(stg_all[, any(surgeries_C > C_requiring_surgery + 1e-6)],
       "surgeries_C exceeds C_requiring_surgery.")
  fail(stg_all[, any(surgeries_D > D_requiring_surgery + 1e-6)],
       "surgeries_D exceeds D_requiring_surgery.")
  fail(stg_all[, any(total_surgeries > C_requiring_surgery + D_requiring_surgery + 1e-6)],
       "total_surgeries exceeds total requiring surgery.")

  # intervention effect direction: SAP averts RHD deaths vs reference (cumulative)
  d_ref <- wsd_all[scenario == "ref", sum(dead)]
  d_sap <- wsd_all[scenario == "sap", sum(dead)]
  fail(d_sap > d_ref + 1e-6, "SAP scenario has MORE cumulative RHD deaths than reference.")

  # order-of-magnitude anchors at the base year (reference arm)
  by1 <- min(years)
  base <- wsd_all[scenario == "ref" & year == by1]
  rhd_d  <- base[, sum(dead)]
  allc_d <- base[, sum(all.mx)]
  # COUNTRY-specific base-year death bands (from 05's meta; Indonesia defaults if absent)
  rhd_lo  <- if (is.null(meta$rhd_death_lo))  1e3 else meta$rhd_death_lo
  rhd_hi  <- if (is.null(meta$rhd_death_hi))  5e4 else meta$rhd_death_hi
  allc_lo <- if (is.null(meta$allc_death_lo)) 5e5 else meta$allc_death_lo
  allc_hi <- if (is.null(meta$allc_death_hi)) 4e6 else meta$allc_death_hi
  fail(rhd_d  < rhd_lo  || rhd_d  > rhd_hi,  sprintf("base-year RHD deaths %.0f outside band %g-%g.", rhd_d, rhd_lo, rhd_hi))
  fail(allc_d < allc_lo || allc_d > allc_hi, sprintf("base-year all-cause deaths %.0f outside band %g-%g.", allc_d, allc_lo, allc_hi))

  attr(wsd_all, "deaths_averted_cum") <- d_ref - d_sap
  list(wsd = wsd_all, stages = stg_all,
       diag = list(deaths_averted_cum = d_ref - d_sap,
                   rhd_deaths_base = rhd_d, allcause_deaths_base = allc_d))
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
                     .export = c("run_location", "baseline_state",
                                 "abcd_one_cycle", "abcd_compete")) %dopar% {
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
  out <- list(wsd = res$wsd, stages = res$stages, diag = res$diag,
              meta = baseline_state$meta)
  fn  <- paste0(OUT_DIR, gsub("[^A-Za-z0-9]+", "_", loc), ".rds")
  saveRDS(out, file = fn)
  message(sprintf("  %-12s -> %s | RHD deaths averted (cum, %d-%d) = %s | base-yr RHD deaths = %s",
                  loc, basename(fn), min(baseline_state$meta$years), max(baseline_state$meta$years),
                  formatC(round(res$diag$deaths_averted_cum), format = "d", big.mark = ","),
                  formatC(round(res$diag$rhd_deaths_base),    format = "d", big.mark = ",")))
}

message("── 06_run_prevention_model.R complete ─────────────────")
message(sprintf("  Wrote %d file(s) to %s ($wsd aggregate + $stages A/B/C/D + surgery trace)",
                length(LOCATIONS), OUT_DIR))
message("  Next: 07_make_outputs.R (tidy long tables), 08_economic_evaluation.R (economics)")
