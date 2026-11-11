# ==============================================================================
# SAP VOI analysis — 02: production-equivalence and analysis checks
# scripts/voi/02_validate_voi.R
#
# Before any PSA draw, per country:
#   V1  production engine (R/abcd_engine.R, isolated env) re-run on the saved baseline_state
#       reproduces the saved 06 outputs (annual deaths, prevalence and service volumes)
#   V2  VOI engine, production SHARED form at the point values == production engine (exact)
#   V3  VOI engine, SEPARATE form with M_AB = M_BC = 0.45 == shared form (exact, every cycle,
#       every cell, final state): the extension does not change production results
#   V4  economics recomputed from the saved applied 2025-US$ unit costs, VSL schedule and
#       discounting reproduce the saved 08 results (totals, BCR, cost/death, cost/DALY, annual)
#   V5  the two multipliers act separately (only the targeted transition is scaled)
#   V6  no effect (M_AB = M_BC = 1): deaths averted vanish; INB = - incremental cost
#   V7  protection-age mask (20 / 100) follows attained age
#   V8  partial-adherence protection rebuilds the care schedule (adherence factor per year)
#   V9  engine input validation refuses invalid / ambiguous SAP specifications
#   V10 effect matrices are re-prepared for every draw (no stale cache carried over)
#   V11 isolation: no engine objects in the global environment; both engines in own envs
# Hard checks stop the analysis on failure. Writes validation/validation_checks.csv,
# validation/point_run_annual.csv, validation/point_run_economics.csv.
# ==============================================================================

if (!exists("VOI_CFG")) { VOI_CONFIG_ONLY <- TRUE; source(file.path("scripts", "voi", "00_run_voi.R")); rm(VOI_CONFIG_ONLY) }
library(data.table)

voi_validate_country <- function(cfg, country) {
  p <- voi_paths(cfg, country)
  inp <- voi_load_inputs(cfg, country)
  checks <- list()
  add <- function(id, description, value, tolerance, pass, severity = "hard", note = "")
    checks[[length(checks) + 1L]] <<- data.table(country = country, check_id = id, description = description,
                                                 value = as.numeric(value), tolerance = as.numeric(tolerance),
                                                 pass = isTRUE(pass), severity = severity, note = note)
  rel <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))

  # --- engines in isolated environments ---------------------------------------------
  had_global <- exists("abcd_one_cycle", envir = globalenv(), inherits = FALSE)
  eng  <- voi_load_engine(cfg$engine_file)
  prod <- voi_load_engine(cfg$prod_engine_file)
  add("V11a", "engines loaded into isolated environments: no abcd_* objects created in the global environment",
      sum(grepl("^abcd_|^ABCD_", ls(globalenv()))), 0, !had_global && !any(grepl("^abcd_|^ABCD_", ls(globalenv()))))
  add("V11b", "production engine environment unchanged by loading the VOI engine (abcd_sap_hr bodies differ, production body = file)",
      0, 0, !identical(body(eng$abcd_sap_hr), body(prod$abcd_sap_hr)) &&
        identical(deparse(body(prod$abcd_sap_hr)), deparse(body(voi_load_engine(cfg$prod_engine_file)$abcd_sap_hr))))
  add("V0a", "engine spec: saved baseline == VOI engine == production engine", 0, 0,
      identical(inp$meta$engine_spec, eng$ABCD_ENGINE_SPEC) && identical(eng$ABCD_ENGINE_SPEC, prod$ABCD_ENGINE_SPEC))
  th0 <- voi_point_theta(inp)
  add("V0b", "saved production point values equal VOI_POINT (M 0.45 shared, protection age 39, partial 0.5)", 0, 0,
      isTRUE(all.equal(unname(unlist(th0)), unname(unlist(cfg$point)))) &&
        isTRUE(all.equal(1 - inp$st$effects$sap_progression_rrr, inp$st$effects$sap_progression_hazard_mult)))
  add("V0c", "fixed SAP specification: C->D, D->death and prior-A multipliers 1; transitions A_to_B, B_to_C", 0, 0,
      inp$st$effects$sap_hr_C_to_D == 1 && inp$st$effects$sap_hr_D_to_rhd_death == 1 && inp$st$effects$sap_hr_prior_A_to_B_adult == 1 &&
        identical(as.character(inp$st$effects$sap_effect_transitions), c("A_to_B", "B_to_C")))

  nhe_v <- eng$abcd_nh_build(inp$st$natural_history, inp$AGES, inp$SEXES)
  nhe_p <- prod$abcd_nh_build(inp$st$natural_history, inp$AGES, inp$SEXES)
  rel0 <- inp$st$coverage$constants$partial_relative_protection
  run_both <- function(engine, nhe, eff, relp, keep = FALSE) {
    a <- lapply(c("ref", "sap"), function(sc) voi_run_arm(engine, inp, nhe, eff, relp, sc, keep_proj = keep))
    names(a) <- c("ref", "sap"); a
  }
  ann <- function(a) { x <- rbindlist(lapply(a, `[[`, "annual")); x[, draw_id := 0L]; setcolorder(x, c("draw_id", "scenario", "year")); x }
  t0 <- proc.time()[["elapsed"]]
  r_prod   <- run_both(prod, nhe_p, inp$st$effects, rel0, keep = TRUE)
  t_prod <- proc.time()[["elapsed"]] - t0
  r_shared <- run_both(eng, nhe_v, inp$st$effects, rel0, keep = TRUE)
  r_split  <- run_both(eng, nhe_v, voi_make_eff(inp$st$effects, th0, split = TRUE), rel0, keep = TRUE)
  a_prod <- ann(r_prod); a_shared <- ann(r_shared); a_split <- ann(r_split)
  num <- setdiff(names(a_prod), c("draw_id", "scenario", "year"))

  # --- V1: production engine vs saved 06 outputs --------------------------------------
  om <- readRDS(inp$files[["out_model"]])
  wsd <- as.data.table(om$wsd); stg <- as.data.table(om$stages)
  saved <- merge(wsd[, .(rhd_deaths = sum(dead), sick_end = sum(sick)), by = .(scenario, year)],
                 stg[, .(n_screened = sum(n_screened), on_bpg_py = sum(n_on_optimal_treatment), surgeries = sum(total_surgeries),
                         confirmations = sum(n_confirm_completed), new_diagnoses = sum(n_new_diagnoses),
                         rhd_deaths_0_19 = sum(rhd_deaths[age < 20]), eff_protected_py = sum(n_effectively_protected)),
                     by = .(scenario, year)], by = c("scenario", "year"))
  cmp <- merge(saved, a_prod, by = c("scenario", "year"), suffixes = c("_saved", "_run"))
  v1_cols <- c("rhd_deaths", "sick_end", "n_screened", "on_bpg_py", "surgeries", "confirmations", "new_diagnoses", "rhd_deaths_0_19", "eff_protected_py")
  v1 <- max(vapply(v1_cols, function(cn) rel(cmp[[paste0(cn, "_run")]], cmp[[paste0(cn, "_saved")]]), numeric(1)))
  add("V1", "production engine re-run reproduces saved 06 outputs (annual deaths, deaths 0-19, prevalence, screening, on-BPG PY, surgeries, confirmations, diagnoses, protected PY)",
      v1, cfg$tol$saved_outputs, nrow(cmp) == 2L * length(inp$years) && v1 <= cfg$tol$saved_outputs)
  add("V1b", "production run flags: no population deficit, care violation or onset-cap hit", 0, 0,
      all(vapply(r_prod, function(a) all(a$checks == 0), logical(1))))

  # --- V2 / V3: engine identity (annual + every cycle matrix + final state) -----------
  same_cells <- function(x, y) {
    d <- 0
    for (sc in c("ref", "sap")) {
      px <- x[[sc]]$proj; py <- y[[sc]]$proj
      for (t in seq_along(px$cycles)) for (nm in c("rhd_deaths", "sick_end", "n_on_optimal_treatment", "n_screened", "total_surgeries",
                                                   "sap_prog_mult_A_to_B", "sap_prog_mult_B_to_C"))
        d <- max(d, max(abs(px$cycles[[t]][[nm]] - py$cycles[[t]][[nm]])))
      for (k in names(px$final_state)) d <- max(d, max(abs(px$final_state[[k]] - py$final_state[[k]])))
    }
    d
  }
  v2a <- max(abs(as.matrix(a_shared[, ..num]) - as.matrix(a_prod[, ..num]))); v2c <- same_cells(r_shared, r_prod)
  add("V2", "VOI engine (production shared multiplier) == production engine: annual sums, every cycle matrix, final state (max abs diff)",
      max(v2a, v2c), cfg$tol$engine_identity, max(v2a, v2c) <= cfg$tol$engine_identity)
  v3a <- max(abs(as.matrix(a_split[, ..num]) - as.matrix(a_shared[, ..num]))); v3c <- same_cells(r_split, r_shared)
  add("V3", "VOI engine SEPARATE multipliers at M_AB = M_BC = 0.45 == shared form (max abs diff, every cycle/cell/final state)",
      max(v3a, v3c), cfg$tol$engine_identity, max(v3a, v3c) <= cfg$tol$engine_identity)
  hr_s <- eng$abcd_sap_hr(inp$st$effects, nhe_v); hr_x <- eng$abcd_sap_hr(voi_make_eff(inp$st$effects, th0), nhe_v)
  hr_p <- prod$abcd_sap_hr(inp$st$effects, nhe_p)
  add("V3b", "effect matrices identical: production abcd_sap_hr vs VOI shared vs VOI separate at the point values", 0, 0,
      identical(hr_s, hr_p) && identical(hr_x, hr_p))

  # --- V4: economics vs saved 08 results ---------------------------------------------
  voi_add_costs(a_split, inp$econ)
  tot <- voi_econ_totals(a_split, inp$econ)
  s8 <- inp$econ$summary
  pairs <- list(c("C_d", "total_cost_disc"), c("B_d", "total_benefit_disc"), c("INB_vsl_d", "net_benefit_disc"), c("bcr_d", "bcr_disc"),
                c("C_u", "total_cost_undisc"), c("B_u", "total_benefit_undisc"), c("INB_vsl_u", "net_benefit_undisc"), c("bcr_u", "bcr_undisc"),
                c("DA", "deaths_averted"), c("cost_per_death", "cost_per_death"), c("cost_per_daly", "cost_per_daly"))
  v4 <- max(vapply(pairs, function(pp) rel(tot[[pp[1]]], s8[[pp[2]]]), numeric(1)))
  add("V4a", "economic totals reproduce saved 08 summary (disc/undisc cost, benefit, net benefit, BCR, deaths averted, cost/death, cost/DALY)",
      v4, cfg$tol$economics, v4 <= cfg$tol$economics)
  by_y <- voi_econ_by_year(a_split, inp$econ)
  ci <- merge(by_y, inp$econ$incremental[, .(year, s_cost_ref = cost_ref, s_cost_sap = cost_sap, s_death_ref = death_ref, s_death_sap = death_sap)], by = "year")
  v4b <- max(rel(ci$cost_ref, ci$s_cost_ref), rel(ci$cost_sap, ci$s_cost_sap), rel(ci$rhd_deaths_ref, ci$s_death_ref), rel(ci$rhd_deaths_sap, ci$s_death_sap))
  add("V4b", "annual arm costs and RHD deaths reproduce saved 08 $incremental (every year)", v4b, cfg$tol$economics, v4b <= cfg$tol$economics)
  add("V4c", "unit costs taken as saved APPLIED 2025 US$ (no repeated price conversion); VSL schedule = saved 08 schedule", 0, 0,
      all(is.finite(inp$econ$uc)) && identical(inp$econ$annual$vsl, inp$econ$incremental$vsl), "info",
      sprintf("screen %.4f, SAP %.2f per PY, surgery %.0f (constant 2025 US$)", inp$econ$uc[["screen"]], inp$econ$uc[["sap"]], inp$econ$uc[["surgery"]]))

  # --- V5: separate multipliers act separately ---------------------------------------
  run_theta <- function(th, keep = TRUE) voi_run_draw(eng, inp, nhe_v, th, 0L, split = TRUE, keep_proj = keep)
  mk <- function(...) { th <- th0; u <- list(...); for (nm in names(u)) th[[nm]] <- u[[nm]]; th }
  mult_check <- function(r, which_one, M, amax) {   # returns max deviation of the on-BPG multiplier from its definition
    dev <- 0
    for (sc in c("ref", "sap")) {
      pj <- r$arms[[sc]]$proj; adh <- r$arms[[sc]]$annual$adherence_factor
      for (t in seq_along(pj$cycles)) {
        cy <- pj$cycles[[t]]
        prot <- matrix(as.numeric(inp$AGES <= amax), length(inp$AGES), length(inp$SEXES))
        exp_ab <- prot * (1 - adh[t] * (1 - M[["AB"]])) + (1 - prot)
        exp_bc <- prot * (1 - adh[t] * (1 - M[["BC"]])) + (1 - prot)
        dev <- max(dev, max(abs(cy$sap_prog_mult_A_to_B - exp_ab)), max(abs(cy$sap_prog_mult_B_to_C - exp_bc)))
      }
    }
    dev
  }
  r_ab <- run_theta(mk(M_BC = 1)); r_bc <- run_theta(mk(M_AB = 1))
  d5a <- mult_check(r_ab, "AB", c(AB = 0.45, BC = 1), 39L); d5b <- mult_check(r_bc, "BC", c(AB = 1, BC = 0.45), 39L)
  add("V5a", "M_AB = 0.45, M_BC = 1: A->B multiplier = 1 - a(1 - 0.45) at age <= 39, B->C multiplier = 1 everywhere", d5a, 1e-15, d5a <= 1e-15)
  add("V5b", "M_AB = 1, M_BC = 0.45: A->B multiplier = 1 everywhere, B->C multiplier = 1 - a(1 - 0.45) at age <= 39", d5b, 1e-15, d5b <= 1e-15)
  da <- function(r) { a <- copy(r$annual); voi_add_costs(a, inp$econ); voi_econ_totals(a, inp$econ)$DA }
  da_pt <- tot$DA; da_ab <- da(r_ab); da_bc <- da(r_bc)
  add("V5c", "removing either channel lowers deaths averted (informational; deaths averted: point / A->B only / B->C only)",
      da_pt, NA, da_ab < da_pt && da_bc < da_pt, "info", sprintf("point %.1f | A->B effect only %.1f | B->C effect only %.1f", da_pt, da_ab, da_bc))

  # --- V6: no effect ----------------------------------------------------------------
  r_null <- run_theta(mk(M_AB = 1, M_BC = 1), keep = FALSE)
  a_null <- copy(r_null$annual); voi_add_costs(a_null, inp$econ); t_null <- voi_econ_totals(a_null, inp$econ)
  tot_deaths <- a_null[scenario == "ref", sum(rhd_deaths)]
  v6 <- max(abs(t_null$DA) / tot_deaths, abs(t_null$B_d) / abs(tot$B_d))
  add("V6", "M_AB = M_BC = 1 (no SAP effect): |deaths averted| / total deaths and |benefit| / point benefit ~ 0, so INB = -incremental cost",
      v6, cfg$tol$no_effect_deaths, v6 <= cfg$tol$no_effect_deaths && t_null$C_d > 0,
      note = sprintf("deaths averted %.3g; INB %.4g = -cost %.4g", t_null$DA, t_null$INB_vsl_d, t_null$C_d))

  # --- V7: protection-age mask -------------------------------------------------------
  r20 <- run_theta(mk(sap_protection_age_max = 20L)); r100 <- run_theta(mk(sap_protection_age_max = 100L))
  d7a <- mult_check(r20, "both", c(AB = 0.45, BC = 0.45), 20L); d7b <- mult_check(r100, "both", c(AB = 0.45, BC = 0.45), 100L)
  above20 <- max(vapply(r20$arms$sap$proj$cycles, function(cy) max(abs(cy$sap_prog_mult_B_to_C[inp$AGES > 20, ] - 1)), numeric(1)))
  add("V7a", "protection age 20: multipliers follow attained age <= 20 and equal 1 above 20", max(d7a, above20), 1e-15, d7a <= 1e-15 && above20 == 0)
  add("V7b", "protection age 100 (no cap): multipliers < 1 at every attained age (adherence > 0)", d7b, 1e-15,
      d7b <= 1e-15 && all(vapply(r100$arms$sap$proj$cycles, function(cy) all(cy$sap_prog_mult_B_to_C < 1), logical(1))))

  # --- V8: partial-adherence protection rebuilds the care schedule -------------------
  r_rel <- run_theta(mk(partial_adherence_relative_protection = 0.2))
  exp_adh <- function(sc) inp$st$coverage[[sc]]$adherence_good + 0.2 * inp$st$coverage[[sc]]$adherence_partial
  d8 <- max(abs(r_rel$annual[scenario == "ref", adherence_factor] - exp_adh("ref")), abs(r_rel$annual[scenario == "sap", adherence_factor] - exp_adh("sap")))
  d8b <- mult_check(r_rel, "both", c(AB = 0.45, BC = 0.45), 39L)
  add("V8", "partial protection 0.2: adherence factor rebuilt every year in both arms (ref 0.64; scale-up 0.64 -> 0.87) and used by the multipliers",
      max(d8, d8b), 1e-15, max(d8, d8b) <= 1e-15,
      note = sprintf("ref %.3f; scale-up %.3f -> %.3f", r_rel$annual[scenario == "ref", adherence_factor][1],
                     r_rel$annual[scenario == "sap", adherence_factor][1], tail(r_rel$annual[scenario == "sap", adherence_factor], 1)))

  # --- V9: engine input validation ----------------------------------------------------
  errs <- function(e) inherits(tryCatch(eng$abcd_sap_hr(e, nhe_v), error = function(x) x), "error")
  base_split <- voi_make_eff(inp$st$effects, th0, split = TRUE)
  bad <- list(
    M_AB_zero = within(base_split, sap_hr_A_to_B <- 0),
    M_BC_above_1 = within(base_split, sap_hr_B_to_C <- 1.2),
    M_AB_NA = within(base_split, sap_hr_A_to_B <- NA_real_),
    only_one_split = within(base_split, rm(sap_hr_B_to_C)),
    shared_and_split = c(base_split, list(sap_progression_hazard_mult = 0.45, sap_progression_rrr = 0.55)),
    per_transition_rrr_mismatch = c(base_split, list(sap_progression_rrr_A_to_B = 0.30)),
    protection_age_non_integer = within(base_split, sap_protection_age_max <- 39.5),
    protection_age_above_max = within(base_split, sap_protection_age_max <- 101L),
    C_to_D_not_1 = within(base_split, sap_hr_C_to_D <- 0.9),
    D_death_not_1 = within(base_split, sap_hr_D_to_rhd_death <- 0.9),
    shared_rrr_mismatch = within(inp$st$effects, sap_progression_rrr <- 0.5))
  refused <- vapply(bad, errs, logical(1))
  ok_split <- !errs(base_split) && !errs(c(base_split, list(sap_progression_rrr_A_to_B = 0.55, sap_progression_rrr_B_to_C = 0.55))) &&
    !errs(within(base_split, { sap_hr_A_to_B <- 1; sap_hr_B_to_C <- 0.01 }))
  add("V9", "engine refuses invalid / ambiguous SAP specifications and accepts valid ones", sum(!refused), 0, all(refused) && ok_split,
      note = if (all(refused)) paste("refused:", paste(names(bad), collapse = ", ")) else paste("NOT refused:", paste(names(bad)[!refused], collapse = ", ")))

  # --- V10: effect matrices prepared per draw (stale cache never used) ----------------
  stale <- eng$abcd_eff_prepare(voi_make_eff(inp$st$effects, mk(M_AB = 0.9, M_BC = 0.9)), nhe_v)
  stale$sap_hr_A_to_B <- 0.45; stale$sap_hr_B_to_C <- 0.45            # new draw values, stale .sap_hr still attached
  care_s <- eng$abcd_care_schedule(c(inp$st$coverage$sap[VOI_CARE_VEC], inp$st$coverage$constants), length(inp$years), inp$AGES, inp$SEXES)
  pj <- eng$abcd_project(init = inp$st$seed_state, pop = inp$st$pop, ir = inp$st$ir, oth = inp$st$oth_mort, nhe = nhe_v, eff = stale,
                         surg = inp$st$surgery, care_year = care_s, years = inp$years, light = FALSE)
  d10 <- max(abs(pj$cycles[[1]]$sap_prog_mult_B_to_C - r_split$sap$proj$cycles[[1]]$sap_prog_mult_B_to_C),
             abs(sum(pj$cycles[[length(inp$years)]]$rhd_deaths) - sum(r_split$sap$proj$cycles[[length(inp$years)]]$rhd_deaths)))
  add("V10", "abcd_project re-prepares the effect matrices for each draw (a stale cached .sap_hr from another draw is overwritten)", d10, 0, d10 == 0)
  add("V11c", "global environment still free of engine objects after all validation runs", 0, 0, !any(grepl("^abcd_|^ABCD_", ls(globalenv()))))

  res <- rbindlist(checks)
  voi_fwrite(res, file.path(p$validation, "validation_checks.csv"), cfg)
  voi_fwrite(a_split, file.path(p$validation, "point_run_annual.csv"), cfg)
  voi_fwrite(cbind(country = country, tot, run = "production point values via VOI engine (separate multipliers 0.45 / 0.45)"),
             file.path(p$validation, "point_run_economics.csv"), cfg)
  voi_manifest_update(p, cfg, validation = list(n_checks = nrow(res), n_hard_failed = res[severity == "hard" & !pass, .N],
                                                n_info_flagged = res[severity == "info" & !pass, .N],
                                                seconds_per_paired_run = t_prod,
                                                tolerances = cfg$tol),
                      input_fingerprints = inp$md5, engine_files = voi_md5(c(cfg$engine_file, cfg$prod_engine_file)),
                      engine_version = eng$ABCD_VOI_ENGINE_VERSION, engine_spec = eng$ABCD_ENGINE_SPEC,
                      fixed_assumptions = list(
                        calibration = "04b calibrated natural history, incidence and bridged 2026 engine state held FIXED (from baseline_state.rds); produced under the production point SAP values (reference care)",
                        calibration_adequacy = inp$calibration_adequacy,
                        natural_history = inp$st$natural_history[c("tr", "adult_prog_hazard_mult", "p_rhd_death_D_child", "p_rhd_death_D_adult",
                                                                    "d_death_adult_age_min", "entry_mix_child", "entry_mix_adult", "p_B_to_no_rhd_adult")],
                        care_constants = inp$st$coverage$constants, surgery = inp$st$surgery,
                        unit_costs_usd2025 = inp$econ$uc, disc_rate = inp$econ$disc_rate, discount_base_year = inp$econ$base_year,
                        dalys_per_death = inp$econ$dalys_per_death, horizon = range(inp$years),
                        sap_fixed = list(sap_hr_C_to_D = 1, sap_hr_D_to_rhd_death = 1, sap_hr_prior_A_to_B_adult = 1)))
  for (i in seq_len(nrow(res)))
    voi_log(p, sprintf("02 %-5s %-4s %s | value %s | %s", res$check_id[i], if (res$pass[i]) "PASS" else if (res$severity[i] == "hard") "FAIL" else "FLAG",
                       res$description[i], format(signif(res$value[i], 4)), res$note[i]))
  if (res[severity == "hard" & !pass, .N] > 0)
    stop("[", country, "] VOI validation FAILED (", paste(res[severity == "hard" & !pass, check_id], collapse = ", "),
         "); see ", file.path(p$validation, "validation_checks.csv"), ". The PSA is not run.", call. = FALSE)
  res
}

voi_run_validation <- function(cfg) rbindlist(lapply(cfg$countries, function(ct) voi_validate_country(cfg, ct)))

if (identical(get0("VOI_CURRENT_STAGE", ifnotfound = ""), "validate") || !exists("VOI_CURRENT_STAGE")) VOI_VALIDATION <- voi_run_validation(VOI_CFG)
