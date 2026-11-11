# ==============================================================================
# SAP VOI analysis: shared helpers (scripts/voi/R/voi_functions.R)
#
# Paths, logging, write guards, isolated engine loading, saved-input loading, the per-draw
# runner (both arms, compact annual summaries) and the economics recomputation. Contains NO
# engine code: every engine call goes through an environment returned by voi_load_engine().
# Sourced by scripts/voi/00_run_voi.R and by every PSA worker process.
# ==============================================================================

VOI_CARE_VEC <- c("screen_window_prob", "clinical_detection", "bpg_initiation", "bpg_retention",
                  "adherence_good", "adherence_partial", "surgery")

# ------------------------------------------------------------------------------
# paths, logging, write guards (every VOI write must land under output/voi/<run_id>/)
# ------------------------------------------------------------------------------
voi_paths <- function(cfg, country, create = TRUE) {
  base <- file.path(cfg$out_root, cfg$run_id, country)
  p <- list(base = base, draws = file.path(base, "draws"), validation = file.path(base, "validation"),
            pilot = file.path(base, "pilot"), pilot_checkpoints = file.path(base, "pilot", "checkpoints"),
            psa = file.path(base, "psa"), checkpoints = file.path(base, "psa", "checkpoints"),
            tables = file.path(base, "tables"), figures = file.path(base, "figures"), logs = file.path(base, "logs"))
  if (create) for (d in p) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  p
}
voi_assert_out <- function(path, cfg) {
  a <- normalizePath(path, winslash = "/", mustWork = FALSE)
  r <- normalizePath(file.path(cfg$out_root, cfg$run_id), winslash = "/", mustWork = FALSE)
  if (!startsWith(a, paste0(r, "/")))
    stop("VOI write guard: refusing to write outside ", r, ": ", a, call. = FALSE)
  invisible(a)
}
voi_fwrite <- function(x, path, cfg, ...) { voi_assert_out(path, cfg); data.table::fwrite(x, path, ...) }
voi_save   <- function(x, path, cfg)      { voi_assert_out(path, cfg); saveRDS(x, path) }
voi_log <- function(paths, ...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  message(msg)
  cat(msg, "\n", file = file.path(paths$logs, "voi_log.txt"), append = TRUE, sep = "")
}
voi_md5 <- function(files) { files <- files[file.exists(files)]; setNames(unname(tools::md5sum(files)), normalizePath(files, winslash = "/")) }
voi_hash_object <- function(x) {          # text deparse -> md5 (stable across sessions; used for prior / draw identity)
  tf <- tempfile(fileext = ".txt"); on.exit(unlink(tf))
  writeLines(deparse(x, control = c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17")), tf)
  unname(tools::md5sum(tf))
}

# manifest: one list per country, updated by each stage (rds = full record; csv = flat view)
voi_manifest_update <- function(paths, cfg, ...) {
  f <- file.path(paths$base, "manifest.rds")
  m <- if (file.exists(f)) readRDS(f) else list()
  upd <- list(...)
  for (nm in names(upd)) m[[nm]] <- upd[[nm]]
  m$last_updated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  voi_save(m, f, cfg)
  flat <- function(x, pre = "") {
    if (is.list(x) && !is.data.frame(x)) {
      if (!length(x)) return(NULL)
      nms <- names(x); if (is.null(nms)) nms <- seq_along(x)
      rbindlist(lapply(seq_along(x), function(i) flat(x[[i]], paste0(pre, if (nzchar(pre)) "." else "", nms[i]))), fill = TRUE)
    } else if (is.data.frame(x)) data.table(field = pre, value = sprintf("<table %d x %d>", nrow(x), ncol(x)))
    else data.table(field = pre, value = paste(format(x, digits = 15), collapse = "; "))
  }
  voi_fwrite(flat(m), file.path(paths$base, "manifest.csv"), cfg)
  invisible(m)
}

# ------------------------------------------------------------------------------
# isolated engine loading: the engine lives in its own environment (never globalenv)
# ------------------------------------------------------------------------------
voi_load_engine <- function(file) {
  if (!file.exists(file)) stop("Engine file not found: ", file, call. = FALSE)
  env <- new.env(parent = globalenv())
  sys.source(file, envir = env, keep.source = FALSE)
  env
}

# ------------------------------------------------------------------------------
# production integrity: fingerprints of every production file the VOI run could touch
# ------------------------------------------------------------------------------
voi_production_fingerprint <- function(cfg) {
  ls_files <- function(d, pattern = NULL) { if (!dir.exists(d)) return(character(0))
    f <- list.files(d, pattern = pattern, full.names = TRUE); f[!dir.exists(f)] }
  files <- c(cfg$prod_engine_file, ls_files(file.path(cfg$root, "R")), ls_files(file.path(cfg$root, "scripts"), "\\.R$"),
             ls_files(file.path(cfg$root, "output"), "^rhd_country_comparison"))
  for (ct in cfg$countries)
    files <- c(files, ls_files(file.path(cfg$root, "data", ct)), ls_files(file.path(cfg$root, "output", ct, "out_model")),
               ls_files(file.path(cfg$root, "output", ct, "tables")))
  files <- unique(normalizePath(files[file.exists(files)], winslash = "/"))
  data.table(file = files, md5 = unname(tools::md5sum(files)), size = file.size(files))
}
voi_check_integrity <- function(cfg, before, after) {
  m <- merge(before, after, by = "file", all = TRUE, suffixes = c("_before", "_after"))
  m[, unchanged := !is.na(md5_before) & !is.na(md5_after) & md5_before == md5_after]
  for (ct in cfg$countries) {
    p <- voi_paths(cfg, ct)
    voi_fwrite(m, file.path(p$validation, "production_integrity.csv"), cfg)
  }
  if (!all(m$unchanged))
    stop("PRODUCTION INTEGRITY FAILURE: production files changed or disappeared during the VOI run:\n  ",
         paste(m[unchanged == FALSE, file], collapse = "\n  "), call. = FALSE)
  message(sprintf("  production integrity: %d production files unchanged (md5) before vs after the VOI run.", nrow(m)))
  invisible(m)
}

# ------------------------------------------------------------------------------
# saved production inputs (read-only). Missing inputs are a BLOCKER: never regenerated.
# ------------------------------------------------------------------------------
voi_input_files <- function(cfg, country)
  c(baseline_state = file.path(cfg$root, "data", country, "baseline_state.rds"),
    economic_results = file.path(cfg$root, "output", country, "tables", "rhd_economic_results.rds"),
    out_model = file.path(cfg$root, "output", country, "out_model", paste0(gsub("[^A-Za-z0-9]+", "_", country), ".rds")))

voi_load_inputs <- function(cfg, country) {
  f <- voi_input_files(cfg, country)
  miss <- f[!file.exists(f)]
  if (length(miss))
    stop("BLOCKER [", country, "]: required saved production input(s) missing:\n  ", paste(miss, collapse = "\n  "),
         "\n  The VOI analysis never regenerates production outputs; produce them with the production pipeline first.", call. = FALSE)
  bs <- readRDS(f[["baseline_state"]])
  if (!identical(bs$meta$engine_spec, cfg$engine_spec_expected))
    stop("BLOCKER [", country, "]: baseline_state engine '", bs$meta$engine_spec, "' is not the expected production engine '",
         cfg$engine_spec_expected, "'.", call. = FALSE)
  if (!identical(as.character(bs$locations), country))
    stop("BLOCKER [", country, "]: baseline_state holds locations ", paste(bs$locations, collapse = ", "), call. = FALSE)
  er <- readRDS(f[["economic_results"]])
  need_par <- c("cost_screen_per_person", "cost_sap_per_year", "cost_surgery", "disc_rate", "dalys_per_death")
  if (is.null(er$params) || !all(need_par %in% names(er$params)) || is.null(er$discount_base_year) ||
      is.null(er$incremental) || is.null(er$vsl_derivation$anchor$intY_to_usd2025) || is.null(er$summary))
    stop("BLOCKER [", country, "]: rhd_economic_results.rds lacks the applied unit costs / VSL schedule / discounting fields.", call. = FALSE)
  if (!identical(as.character(er$parameters$country), country))
    stop("BLOCKER [", country, "]: rhd_economic_results.rds is for country '", er$parameters$country, "'.", call. = FALSE)
  years <- as.integer(bs$meta$years)
  inc <- as.data.table(er$incremental)[location == country][order(year)]
  if (!identical(as.integer(inc$year), years))
    stop("BLOCKER [", country, "]: economic results horizon differs from the baseline horizon.", call. = FALSE)
  r <- er$params$disc_rate; b0 <- as.integer(er$discount_base_year)
  df <- 1 / (1 + r)^(years - b0)
  if (max(abs(df - inc$df)) > 1e-12) stop("[", country, "] saved discount factors do not match disc_rate / discount_base_year.", call. = FALSE)
  vsl_ann <- as.data.table(er$vsl_derivation$annual)[order(year)]
  if (max(abs(vsl_ann[year %in% years, vsl_usd2025] - inc$vsl)) > 0)
    stop("[", country, "] saved VSL schedule differs between $incremental and $vsl_derivation.", call. = FALSE)
  smry <- as.data.table(er$summary)[location == country]
  gni_ref <- inc[year == min(years), gdp_pc] * er$vsl_derivation$anchor$intY_to_usd2025
  st <- bs$states[[country]]
  list(country = country, files = f, md5 = voi_md5(f), bs = bs, st = st, meta = bs$meta,
       AGES = bs$meta$AGES, SEXES = bs$meta$SEXES, years = years,
       econ = list(uc = c(screen = er$params$cost_screen_per_person, sap = er$params$cost_sap_per_year,
                          surgery = er$params$cost_surgery),
                   disc_rate = r, base_year = b0, dalys_per_death = er$params$dalys_per_death,
                   annual = data.table(year = years, df = df, vsl = inc$vsl),
                   gni_pc_ref_usd2025 = gni_ref, N0 = smry$N0, summary = smry, incremental = inc,
                   currency = er$params$currency, price_year = er$params$price_year,
                   calibration_status = er$calibration_status),
       calibration_adequacy = bs$meta$calibration_adequacy)
}

# production point values of the four SAP inputs, as stored in the saved baseline state
voi_point_theta <- function(inp) {
  e <- inp$st$effects
  list(M_AB = e$sap_progression_hazard_mult, M_BC = e$sap_progression_hazard_mult,
       sap_protection_age_max = as.integer(e$sap_protection_age_max),
       partial_adherence_relative_protection = inp$st$coverage$constants$partial_relative_protection)
}

# effect list for one draw: SEPARATE multipliers (split = TRUE) or the production shared form
voi_make_eff <- function(effects, theta, split = TRUE) {
  eff <- effects
  eff$.sap_hr <- NULL                                   # never carry cached matrices across draws
  if (split) {
    eff$sap_progression_hazard_mult <- NULL; eff$sap_progression_rrr <- NULL
    eff$sap_hr_A_to_B <- as.numeric(theta$M_AB); eff$sap_hr_B_to_C <- as.numeric(theta$M_BC)
  } else {
    if (!isTRUE(all.equal(theta$M_AB, theta$M_BC))) stop("shared SAP form requires M_AB == M_BC.", call. = FALSE)
    eff$sap_progression_hazard_mult <- as.numeric(theta$M_AB); eff$sap_progression_rrr <- 1 - as.numeric(theta$M_AB)
  }
  eff$sap_protection_age_max <- as.integer(theta$sap_protection_age_max)
  eff
}

# ------------------------------------------------------------------------------
# one arm: same call path as production 06 (abcd_care_schedule -> abcd_project, full mode);
# the care schedule is REBUILT for the draw's partial-adherence protection
# ------------------------------------------------------------------------------
voi_run_arm <- function(eng, inp, nhe, eff, rel, scenario, keep_proj = FALSE) {
  st <- inp$st; years <- inp$years
  const <- st$coverage$constants
  const$partial_relative_protection <- as.numeric(rel)
  care <- eng$abcd_care_schedule(c(st$coverage[[scenario]][VOI_CARE_VEC], const), length(years), inp$AGES, inp$SEXES)
  proj <- eng$abcd_project(init = st$seed_state, pop = st$pop, ir = st$ir, oth = st$oth_mort, nhe = nhe, eff = eff,
                           surg = st$surgery, care_year = care, years = years, light = FALSE)
  cyc <- proj$cycles
  S <- function(nm) vapply(cyc, function(cy) sum(cy[[nm]]), numeric(1))
  ch <- inp$AGES < 20
  out <- data.table(scenario = scenario, year = years,
                    rhd_deaths = S("rhd_deaths"),
                    rhd_deaths_0_19 = vapply(cyc, function(cy) sum(cy$rhd_deaths[ch, , drop = FALSE]), numeric(1)),
                    n_screened = S("n_screened"), on_bpg_py = S("n_on_optimal_treatment"),
                    surgeries = S("total_surgeries"), confirmations = S("n_confirm_completed"),
                    new_diagnoses = S("n_new_diagnoses"), sick_end = S("sick_end"), stage_D_end = S("D"),
                    eff_protected_py = S("n_effectively_protected"),
                    flow_A_to_B = S("A_to_B"), flow_B_to_C = S("B_to_C"),
                    adherence_factor = vapply(cyc, function(cy) cy$adherence_factor, numeric(1)))
  chk <- c(pop_deficit = proj$pop_deficit_total, care_violation = proj$care_violation_total,
           onset_cap_hits = proj$onset_cap_hits)
  list(annual = out, checks = chk, proj = if (keep_proj) proj else NULL, care = if (keep_proj) care else NULL)
}

# one paired draw: the SAME sampled SAP inputs in both arms
voi_run_draw <- function(eng, inp, nhe, theta, draw_id = NA_integer_, split = TRUE, keep_proj = FALSE) {
  eff <- voi_make_eff(inp$st$effects, theta, split)
  arms <- lapply(c("ref", "sap"), function(sc)
    voi_run_arm(eng, inp, nhe, eff, theta$partial_adherence_relative_protection, sc, keep_proj))
  names(arms) <- c("ref", "sap")
  annual <- rbindlist(lapply(arms, `[[`, "annual"))
  annual[, draw_id := as.integer(draw_id)]
  setcolorder(annual, c("draw_id", "scenario", "year"))
  chk <- rbind(arms$ref$checks, arms$sap$checks)
  num <- annual[, .SD, .SDcols = setdiff(names(annual), c("draw_id", "scenario", "year"))]
  flags <- list(pop_deficit = sum(chk[, "pop_deficit"]), care_violation = sum(chk[, "care_violation"]),
                onset_cap_hits = sum(chk[, "onset_cap_hits"]),
                nonfinite = sum(!is.finite(as.matrix(num))), negative = sum(as.matrix(num) < -1e-6))
  list(annual = annual, flags = flags, arms = if (keep_proj) arms else NULL)
}

# ------------------------------------------------------------------------------
# economics: the saved APPLIED 2025-US$ unit costs, the saved VSL schedule and discounting
# (no price conversion is repeated). Same arithmetic as production 08.
# ------------------------------------------------------------------------------
voi_add_costs <- function(annual, econ) {
  annual[, cost := n_screened * econ$uc[["screen"]] + on_bpg_py * econ$uc[["sap"]] + surgeries * econ$uc[["surgery"]]]
  invisible(annual)
}
voi_econ_by_year <- function(annual, econ) {
  if (!"cost" %in% names(annual)) { annual <- copy(annual); voi_add_costs(annual, econ) }
  w <- dcast(annual, draw_id + year ~ scenario, value.var = c("cost", "rhd_deaths"))
  w <- merge(w, econ$annual, by = "year")
  setorder(w, draw_id, year)
  w[, `:=`(inc_cost = cost_sap - cost_ref, deaths_avert = rhd_deaths_ref - rhd_deaths_sap)]
  w
}
voi_econ_totals <- function(annual, econ) {
  w <- voi_econ_by_year(annual, econ)
  tot <- w[, .(C_d = sum(df * inc_cost), C_u = sum(inc_cost),
               B_d = sum(df * (deaths_avert * vsl)), B_u = sum(deaths_avert * vsl),
               DA = sum(deaths_avert), DA_d = sum(df * deaths_avert),
               cost_ref_d = sum(df * cost_ref), cost_sap_d = sum(df * cost_sap)), by = draw_id]
  tot[, `:=`(DALY_d = DA_d * econ$dalys_per_death,
             INB_vsl_d = B_d - C_d, INB_vsl_u = B_u - C_u,
             bcr_d = fifelse(C_d > 0, B_d / C_d, NA_real_), bcr_u = fifelse(C_u > 0, B_u / C_u, NA_real_))]
  tot[, `:=`(cost_per_death = fifelse(DA > 0, C_d / DA, NA_real_),
             cost_per_daly = fifelse(DA_d > 0, C_d / DALY_d, NA_real_))]
  tot[]
}
