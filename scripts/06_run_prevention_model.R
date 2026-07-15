# ==============================================================================
# RHD secondary-prevention investment case: MODEL RUNNER
# scripts/06_run_prevention_model.R
#
# Structure after: Coates et al., Lancet Glob Health 2021 (PMC9087136).
# Focus: scale-up of SECONDARY PREVENTION = echo screening + secondary antibiotic
#   prophylaxis (SAP) for asymptomatic (mild) RHD.
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT (the actual model runner)
# ------------------------------------------------------------------------------
# This script RUNS the reference and SAP scale-up scenarios starting from the
# initial state assembled by 05. It does NOT recompute any inputs and it contains
# NO monetary / cost values (economics are in 08 only). The former hard-coded
# Australia scalars (pop_2021, au_pop, prev_rhd_2021, seed_frac, rhd_incident_2021,
# ...) are gone; everything comes from data/baseline_state.rds.
#
#   INPUT : data/baseline_state.rds                      (from 05)
#   OUTPUT: output/out_model/<location>.rds              (one RDS per location)
#
# THE ENGINE  (matrix form: age x sex arrays, Markov cycles as matrix ops)
#   A well-sick-dead Markov model whose "sick" compartment is resolved into three
#   TUNNEL states — mild (asymptomatic) -> severe (heart failure) -> post-surgery:
#     * new incident (mild) RHD enters via  newcases = well x IR   (CALIBRATED IR);
#     * SAP acts on the mild -> severe PROGRESSION (the fatal pathway);
#     * HF management reduces severe RHD mortality; surgery moves severe -> post;
#     * competing (non-RHD) background mortality is data-fed and held fixed;
#     * cohorts age one year per cycle (age_shift); the terminal age is an open 100+.
#   Only SAP coverage differs between the reference and SAP arms.
#
# TWO OUTPUT TABLES PER LOCATION  (this is a model WITH TUNNEL STATES)
#   (1) $wsd    — the WELL-SICK-DEAD aggregate table, reference and SAP scenarios
#                 row-bound and distinguished by a `scenario` column. Here
#                 sick = mild + severe + post (all RHD cases collapsed). This is the
#                 table 07_make_outputs.R consumes to emit the standard long table.
#   (2) $tunnel — the SECOND, TUNNEL-STATE table required by the tunnel structure:
#                 the SAME rows but with the "sick" compartment DISAGGREGATED into
#                 its tunnel sub-states (mild / severe / post) plus the intervention
#                 FLOW volumes (newcases, mild_to_severe, surgeries, op_deaths,
#                 n_on_sap, n_screened, rhd_deaths) and the coverage levels. This is
#                 the companion table the tunnel decomposition produces; 08 uses its
#                 volumes for costing. (Collapsing $tunnel's mild+severe+post exactly
#                 reproduces $wsd's `sick`, by construction.)
#
# eff_ir / eff_cf  (reported RRR multipliers on [0,1]; 1 = no reduction):
#   eff_ir = 1 everywhere — secondary prevention does NOT reduce RHD INCIDENCE.
#   eff_cf = 1 - eff_sap_asymp * cov_sap — SAP acts on the fatal mild->severe
#            progression, i.e. the case-fatality lever in the WSD reduction. This is
#            exactly the multiplier the engine applies, so the reported column is
#            faithful to the mechanics.
#
# PARALLELISM: by LOCATION (foreach/doParallel). Workers call setDTthreads(1). The
#   loop is location-general (Indonesia only for now).
# ==============================================================================

library(data.table)
library(foreach)

if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")
if (!exists("wd_outp")) wd_outp <- paste0(here::here("output"), "/")

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

message(sprintf("── 06_run_prevention_model.R : %d location(s), scenarios %s ──",
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

  clin <- st$clinical; eff <- st$effects; oth <- st$oth_mort

  zero_mat  <- function() matrix(0, n_age, n_sex, dimnames = list(AGES, SEXES))
  age_shift <- function(M) {                         # advance age a -> a+1; 100 open
    N <- zero_mat()
    N[2:n_age, ] <- M[1:(n_age - 1), ]
    N[n_age,  ]  <- N[n_age, ] + M[n_age, ]          # accumulate terminal 100+ group
    N
  }
  # screening age mask (structural: who is echo-screened)
  screen_mask <- zero_mat()
  screen_ages <- as.character(st$coverage$screen_age_lo:st$coverage$screen_age_hi)
  screen_mask[intersect(screen_ages, rownames(screen_mask)), ] <- 1

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
    mild <- severe <- post <- zero_mat()
    wsd_l <- vector("list", n_years); tun_l <- vector("list", n_years)

    for (iy in seq_len(n_years)) {
      popm <- st$pop[, , iy]
      irm  <- st$ir[, , iy]                          # CALIBRATED incidence, trended
      cs <- covs$sap[iy]; ch <- covs$hf[iy]; cg <- covs$surg[iy]

      if (iy == 1L) { mild <- st$seed$mild; severe <- st$seed$severe; post <- st$seed$post }

      # susceptible ("well") = population minus the prevalent RHD (sick) stock
      sick_start <- mild + severe + post
      well       <- pmax(popm - sick_start, 0)

      # --- new incident asymptomatic (mild) RHD: well x IR (eff_ir = 1) ---------
      eff_ir   <- 1                                  # secondary prevention: no incidence effect
      newcases <- well * irm * eff_ir
      mild_pool <- mild + newcases

      # --- mild -> severe progression, slowed by SAP (the case-fatality lever) --
      eff_cf_mult    <- 1 - eff$eff_sap_asymp * cs   # = eff_cf reported below
      prog           <- clin$p_mild_to_severe * eff_cf_mult
      mild_to_severe <- mild_pool * prog
      mild_death_oth <- mild_pool * oth
      mild_next      <- pmax(mild_pool - mild_to_severe - mild_death_oth, 0)

      # --- severe RHD (heart failure): surgery, operative death, HF mortality ---
      severe_pool   <- severe + mild_to_severe
      surg_cand     <- severe_pool * clin$frac_severe_surg_elig
      surgeries     <- surg_cand * cg
      op_deaths     <- surgeries * clin$p_surg_op_mortality
      to_post       <- surgeries - op_deaths
      remain_severe <- severe_pool - surgeries
      sev_death_rhd <- remain_severe * clin$p_severe_death * (1 - eff$eff_hf_mgmt * ch)
      sev_death_oth <- remain_severe * oth
      severe_next   <- pmax(remain_severe - sev_death_rhd - sev_death_oth, 0)

      # --- post-surgery ---------------------------------------------------------
      post_pool      <- post + to_post
      post_death_rhd <- post_pool * clin$p_post_death_rhd
      post_death_oth <- post_pool * oth
      post_next      <- pmax(post_pool - post_death_rhd - post_death_oth, 0)

      # --- year-end quantities --------------------------------------------------
      rhd_deaths <- op_deaths + sev_death_rhd + post_death_rhd     # RHD cause-specific
      sick_end   <- mild_next + severe_next + post_next
      well_end   <- pmax(popm - sick_end, 0)
      # all-cause deaths = RHD deaths + background mortality of the whole population
      all_mx     <- rhd_deaths + oth * popm
      # intervention volumes (counts; costed in 08)
      n_on_sap   <- mild_pool * cs                                # SAP person-years
      n_screened <- popm * screen_mask * cs                       # echo screening volume

      # --- (1) WSD aggregate row set -------------------------------------------
      wsd_l[[iy]] <- melt_year(iy, list(
        well = well_end, sick = sick_end, newcases = newcases,
        dead = rhd_deaths, pop = popm, all.mx = all_mx))[
          , `:=`(eff_ir = eff_ir, eff_cf = eff_cf_mult)]

      # --- (2) tunnel-state row set (disaggregated sick + flow volumes) ---------
      tun_l[[iy]] <- melt_year(iy, list(
        mild = mild_next, severe = severe_next, post = post_next,
        newcases = newcases, mild_to_severe = mild_to_severe,
        surgeries = surgeries, op_deaths = op_deaths,
        n_on_sap = n_on_sap, n_screened = n_screened,
        rhd_deaths = rhd_deaths))[
          , `:=`(cov_sap = cs, cov_hf = ch, cov_surg = cg)]

      # --- age the surviving stocks one year for the next cycle -----------------
      mild   <- age_shift(mild_next)
      severe <- age_shift(severe_next)
      post   <- age_shift(post_next)
    }

    wsd <- rbindlist(wsd_l); tun <- rbindlist(tun_l)
    wsd[, `:=`(scenario = scenario, location = loc, cause = CAUSE,
               intervention = unname(ilab[[scenario]]))]
    tun[, `:=`(scenario = scenario, location = loc, cause = CAUSE,
               intervention = unname(ilab[[scenario]]))]
    list(wsd = wsd, tunnel = tun)
  }

  runs <- lapply(SCEN, run_one)
  wsd_all <- rbindlist(lapply(runs, `[[`, "wsd"))
  tun_all <- rbindlist(lapply(runs, `[[`, "tunnel"))

  # -------------------- per-location SANITY CHECKS --------------------------
  fail <- function(cond, msg) if (isTRUE(cond)) stop("[", loc, "] ", msg, call. = FALSE)

  # completeness: scenarios x years x ages x sexes
  exp_rows <- length(SCEN) * n_years * n_age * n_sex
  fail(nrow(wsd_all) != exp_rows, "WSD table has an incomplete scenario/age/year grid.")
  fail(nrow(tun_all) != exp_rows, "tunnel table has an incomplete scenario/age/year grid.")

  num_wsd <- wsd_all[, .SD, .SDcols = is.numeric]
  fail(anyNA(num_wsd), "WSD table contains NA.")
  fail(any(vapply(num_wsd, function(x) any(x < -1e-6), logical(1))),
       "WSD table contains negative values.")
  fail(anyNA(tun_all[, .SD, .SDcols = is.numeric]), "tunnel table contains NA.")

  # WSD identities: well + sick <= pop ; all.mx >= dead ; eff_* in [0,1]
  fail(wsd_all[, any(well + sick > pop + 1e-3)], "well + sick exceeds pop somewhere.")
  fail(wsd_all[, any(all.mx + 1e-6 < dead)],     "all-cause deaths < RHD deaths somewhere.")
  fail(wsd_all[, any(eff_ir < 0 | eff_ir > 1 | eff_cf < 0 | eff_cf > 1)],
       "eff_ir / eff_cf outside [0,1].")

  # tunnel must reconstruct WSD sick exactly (mild+severe+post == sick)
  chk <- merge(
    tun_all[, .(sick_tun = mild + severe + post), by = .(scenario, sex, age, year)],
    wsd_all[, .(scenario, sex, age, year, sick)],
    by = c("scenario", "sex", "age", "year"))
  fail(chk[, max(abs(sick_tun - sick))] > 1e-6,
       "tunnel mild+severe+post does not reconstruct WSD sick.")

  # intervention effect direction: SAP averts RHD deaths vs reference (cumulative)
  d_ref <- wsd_all[scenario == "ref", sum(dead)]
  d_sap <- wsd_all[scenario == "sap", sum(dead)]
  fail(d_sap > d_ref + 1e-6, "SAP scenario has MORE cumulative RHD deaths than reference.")

  # order-of-magnitude anchors at the base year (reference arm)
  by1 <- min(years)
  base <- wsd_all[scenario == "ref" & year == by1]
  rhd_d  <- base[, sum(dead)]
  allc_d <- base[, sum(all.mx)]
  fail(rhd_d  < 1e3 || rhd_d  > 5e4, sprintf("base-year RHD deaths %.0f outside band 1e3-5e4.", rhd_d))
  fail(allc_d < 5e5 || allc_d > 4e6, sprintf("base-year all-cause deaths %.0f outside band 5e5-4e6.", allc_d))

  attr(wsd_all, "deaths_averted_cum") <- d_ref - d_sap
  attr(wsd_all, "anchor_rhd_deaths_base") <- rhd_d
  list(wsd = wsd_all, tunnel = tun_all,
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
                     .export = c("run_location", "baseline_state")) %dopar% {
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
  out <- list(wsd = res$wsd, tunnel = res$tunnel, diag = res$diag,
              meta = baseline_state$meta)
  fn  <- paste0(OUT_DIR, gsub("[^A-Za-z0-9]+", "_", loc), ".rds")
  saveRDS(out, file = fn)
  message(sprintf("  %-12s -> %s | RHD deaths averted (cum, %d-%d) = %s | base-yr RHD deaths = %s",
                  loc, basename(fn), min(baseline_state$meta$years), max(baseline_state$meta$years),
                  formatC(round(res$diag$deaths_averted_cum), format = "d", big.mark = ","),
                  formatC(round(res$diag$rhd_deaths_base),    format = "d", big.mark = ",")))
}

message("── 06_run_prevention_model.R complete ─────────────────")
message(sprintf("  Wrote %d file(s) to %s", length(LOCATIONS), OUT_DIR))
message("  Next: 07_make_outputs.R (tidy long table), 08_economic_evaluation.R (economics)")
