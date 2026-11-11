# ==============================================================================
# SAP VOI analysis — 03: paired reference / scale-up projections and economics
# scripts/voi/03_run_psa.R
#
# For every draw (one row of draws/sap_draws.csv) BOTH arms (ref, sap) are projected
# 2026-2100 from the saved 2026 engine state with the SAME sampled SAP inputs; the care
# schedule is rebuilt for the draw's partial-adherence protection and the effect matrices
# are re-prepared by abcd_project. Service volumes (screening encounters, on-BPG person-
# years, surgeries; confirmations uncosted, as in production) are recomputed per draw and
# costed with the saved applied 2025-US$ unit costs; deaths are valued with the saved VSL
# schedule and discounting.
#
# Stage "pilot": draws 1..N_PILOT -> pilot/ (timing, numerical checks, projected runtime).
# Stage "psa"  : draws 1..N_DRAWS -> psa/   (resumable chunk checkpoints; pilot draws must be
#                reproduced exactly).
# Workers are separate R processes (PSOCK); each loads the VOI engine into its own
# environment. A checkpoint is reused only if its signature (draw ids + parameter values +
# input / engine / helper fingerprints) matches, so a later run with a larger N_DRAWS only
# simulates the new draws.
#
# Writes (pilot/ or psa/): annual_summary.rds + annual_summary.csv.gz (draw x arm x year:
#   deaths, deaths 0-19, volumes, cost, prevalence, flows, adherence factor),
#   psa_draw_results.csv/.rds (draw parameters + discounted / undiscounted totals),
#   run_flags.csv, failures.csv; validation/pilot_reproducibility.csv (psa stage).
# ==============================================================================

if (!exists("VOI_CFG")) { VOI_CONFIG_ONLY <- TRUE; source(file.path("scripts", "voi", "00_run_voi.R")); rm(VOI_CONFIG_ONLY) }
library(data.table)

# worker: one chunk of draws (runs inside a PSOCK worker that holds inp, cfg, VOI_ENG, VOI_NHE)
voi_psa_chunk <- function(job) {
  t0 <- proc.time()[["elapsed"]]
  ann <- list(); fl <- list(); fails <- list()
  for (i in seq_len(nrow(job$draws))) {
    d <- job$draws[i]
    th <- list(M_AB = d$M_AB, M_BC = d$M_BC, sap_protection_age_max = as.integer(d$sap_protection_age_max),
               partial_adherence_relative_protection = d$partial_adherence_relative_protection)
    t1 <- proc.time()[["elapsed"]]
    r <- tryCatch(voi_run_draw(VOI_ENG, inp, VOI_NHE, th, d$draw_id, split = TRUE), error = function(e) e)
    if (inherits(r, "error")) { fails[[length(fails) + 1L]] <- data.table(draw_id = d$draw_id, error = conditionMessage(r)); next }
    ann[[length(ann) + 1L]] <- r$annual
    fl[[length(fl) + 1L]] <- as.data.table(c(list(draw_id = d$draw_id, seconds = proc.time()[["elapsed"]] - t1), r$flags))
  }
  out <- list(signature = job$signature, ids = job$ids, annual = rbindlist(ann), flags = rbindlist(fl),
              failures = if (length(fails)) rbindlist(fails) else data.table(draw_id = integer(0), error = character(0)),
              engine_version = VOI_ENG$ABCD_VOI_ENGINE_VERSION, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              worker_pid = Sys.getpid())
  voi_save(out, job$file, cfg)
  list(ids = job$ids, n_fail = length(fails), seconds = proc.time()[["elapsed"]] - t0)
}

voi_run_psa_country <- function(cfg, country, stage = c("pilot", "psa")) {
  stage <- match.arg(stage)
  p <- voi_paths(cfg, country)
  inp <- voi_load_inputs(cfg, country)
  f_draws <- file.path(p$draws, "sap_draws.csv")
  if (!file.exists(f_draws)) stop("[", country, "] no sampled draws; run stage 'sample' first.", call. = FALSE)
  if (!file.exists(file.path(p$validation, "validation_checks.csv"))) stop("[", country, "] validation not run; run stage 'validate' first.", call. = FALSE)
  vchk <- fread(file.path(p$validation, "validation_checks.csv"))
  if (vchk[severity == "hard" & pass == FALSE, .N] > 0) stop("[", country, "] validation has hard failures; PSA not run.", call. = FALSE)
  draws <- fread(f_draws)
  ids <- if (stage == "pilot") seq_len(cfg$n_pilot) else seq_len(cfg$n_draws)
  if (!all(ids %in% draws$draw_id)) stop("[", country, "] draws file lacks ids up to ", max(ids), "; re-run stage 'sample'.", call. = FALSE)
  ck_dir <- if (stage == "pilot") p$pilot_checkpoints else p$checkpoints
  out_dir <- if (stage == "pilot") p$pilot else p$psa
  par_cols <- c("draw_id", "M_AB", "M_BC", "sap_protection_age_max", "partial_adherence_relative_protection")
  fp <- unname(c(inp$md5, voi_md5(c(cfg$engine_file, cfg$functions_file))))
  chunks <- split(ids, (ids - 1L) %/% cfg$chunk)
  ck_file <- function(ch) file.path(ck_dir, sprintf("chunk_%05d_%05d.rds", min(ch), max(ch)))
  sig <- function(ch) voi_hash_object(list(as.list(draws[draw_id %in% ch, ..par_cols]), fp))
  valid_ck <- function(ch) { f <- ck_file(ch); file.exists(f) && identical(tryCatch(readRDS(f)$signature, error = function(e) ""), sig(ch)) }
  todo <- Filter(function(ch) !valid_ck(ch), chunks)
  voi_log(p, sprintf("03 %s: %d draws in %d chunk(s); %d chunk(s) to simulate, %d reused from checkpoints",
                     stage, length(ids), length(chunks), length(todo), length(chunks) - length(todo)))
  t_run <- 0; n_workers <- 0L
  if (length(todo)) {
    n_workers <- as.integer(max(1L, min(cfg$max_cores, length(todo), parallel::detectCores() - 1L)))
    jobs <- lapply(todo, function(ch) list(ids = ch, draws = draws[draw_id %in% ch], file = ck_file(ch), signature = sig(ch)))
    t0 <- proc.time()[["elapsed"]]
    cl <- parallel::makeCluster(n_workers)
    res <- tryCatch({
      parallel::clusterExport(cl, c("inp", "cfg", "voi_psa_chunk"), envir = environment())
      parallel::clusterEvalQ(cl, {
        library(data.table); setDTthreads(1L)
        source(cfg$functions_file)
        VOI_ENG <- voi_load_engine(cfg$engine_file)                      # isolated env inside the worker process
        VOI_NHE <- VOI_ENG$abcd_nh_build(inp$st$natural_history, inp$AGES, inp$SEXES)
        NULL
      })
      parallel::parLapplyLB(cl, jobs, voi_psa_chunk)
    }, finally = parallel::stopCluster(cl))
    t_run <- proc.time()[["elapsed"]] - t0
    voi_log(p, sprintf("03 %s: simulated %d draws on %d worker(s) in %.1f s (%d failure(s))",
                       stage, sum(lengths(lapply(res, `[[`, "ids"))), n_workers, t_run, sum(vapply(res, `[[`, 0L, "n_fail"))))
  }

  # ---- combine checkpoints for exactly these ids ----
  cks <- lapply(chunks, function(ch) { if (!valid_ck(ch)) stop("[", country, "] checkpoint missing/invalid after run: ", ck_file(ch), call. = FALSE); readRDS(ck_file(ch)) })
  annual <- rbindlist(lapply(cks, `[[`, "annual"))
  flags <- rbindlist(lapply(cks, `[[`, "flags"))
  failures <- rbindlist(lapply(cks, `[[`, "failures"))
  flags[, numerical_failure := pop_deficit > 1e-6 | care_violation > 1e-6 | onset_cap_hits > 0 | nonfinite > 0 | negative > 0]
  bad_ids <- sort(unique(c(failures$draw_id, flags[numerical_failure == TRUE, draw_id])))
  # paired completeness: each successful draw has both arms x every year
  cnt <- annual[, .N, by = .(draw_id, scenario)]
  incomplete <- setdiff(unique(annual$draw_id), cnt[N == length(inp$years), .N, by = draw_id][N == 2L, draw_id])
  bad_ids <- sort(unique(c(bad_ids, incomplete)))
  if (length(bad_ids)) voi_log(p, sprintf("03 %s: WARNING %d draw(s) failed or flagged and are EXCLUDED: %s", stage, length(bad_ids), paste(bad_ids, collapse = ", ")))
  annual_ok <- annual[!draw_id %in% bad_ids]
  voi_add_costs(annual_ok, inp$econ)
  tot <- voi_econ_totals(annual_ok, inp$econ)
  res_dt <- merge(draws[draw_id %in% ids], tot, by = "draw_id", all.x = TRUE)
  res_dt[, `:=`(status = fifelse(draw_id %in% bad_ids, "failed_or_flagged", "ok"), country = country)]
  setcolorder(res_dt, c("country", "draw_id", "status"))
  if (res_dt[status == "ok", any(!is.finite(C_d) | !is.finite(B_d) | !is.finite(DALY_d))])
    stop("[", country, "] non-finite economic totals among successful draws.", call. = FALSE)

  # ---- write ----
  voi_save(annual_ok, file.path(out_dir, "annual_summary.rds"), cfg)
  voi_fwrite(annual_ok, file.path(out_dir, "annual_summary.csv.gz"), cfg, compress = "gzip")
  voi_save(res_dt, file.path(out_dir, "psa_draw_results.rds"), cfg)
  voi_fwrite(res_dt, file.path(out_dir, "psa_draw_results.csv"), cfg)
  voi_fwrite(flags, file.path(out_dir, "run_flags.csv"), cfg)
  voi_fwrite(if (nrow(failures)) failures else data.table(draw_id = integer(0), error = character(0)), file.path(out_dir, "failures.csv"), cfg)
  sec <- flags$seconds
  info <- list(stage = stage, n_requested = length(ids), n_ok = res_dt[status == "ok", .N], n_failed_or_flagged = length(bad_ids),
               seconds_simulated = t_run, workers = n_workers, mean_seconds_per_paired_draw = mean(sec),
               engine_version = unique(vapply(cks, `[[`, "", "engine_version")),
               input_fingerprints = inp$md5, engine_md5 = voi_md5(cfg$engine_file), functions_md5 = voi_md5(cfg$functions_file),
               script_md5 = voi_md5(file.path(cfg$voi_dir, c("00_run_voi.R", "01_sample_sap_parameters.R", "02_validate_voi.R", "03_run_psa.R", "04_analyse_and_plot.R"))))
  if (stage == "pilot") {
    info$projected_main_minutes <- mean(sec) * cfg$n_draws / max(1L, min(cfg$max_cores, parallel::detectCores() - 1L)) / 60
    voi_log(p, sprintf("03 pilot: %d/%d draws ok | %.2f s per paired draw | INB (VSL) range %.3g to %.3g | projected main run ~%.1f min",
                       info$n_ok, length(ids), mean(sec), min(res_dt$INB_vsl_d, na.rm = TRUE), max(res_dt$INB_vsl_d, na.rm = TRUE), info$projected_main_minutes))
  } else {
    # pilot draws must be reproduced exactly (deterministic engine, same streams)
    pf <- file.path(p$pilot, "psa_draw_results.rds")
    if (file.exists(pf)) {
      pr <- readRDS(pf)
      cmp <- merge(pr[status == "ok", .(draw_id, C_d, B_d, DALY_d, DA)], res_dt[, .(draw_id, C_d, B_d, DALY_d, DA)], by = "draw_id")
      dmax <- cmp[, max(abs(C_d.x - C_d.y), abs(B_d.x - B_d.y), abs(DALY_d.x - DALY_d.y), abs(DA.x - DA.y))]
      rep_tab <- data.table(country = country, n_pilot_compared = nrow(cmp), max_abs_diff = dmax, reproduced_exactly = dmax == 0)
      voi_fwrite(rep_tab, file.path(p$validation, "pilot_reproducibility.csv"), cfg)
      voi_log(p, sprintf("03 psa: pilot draws reproduced in the main run: %d compared, max abs diff %g", nrow(cmp), dmax))
      if (dmax != 0) stop("[", country, "] main PSA does not reproduce the pilot draws exactly.", call. = FALSE)
    }
    voi_log(p, sprintf("03 psa: %d/%d draws ok | E[INB, VSL] = %.4g | P(INB > 0) = %.3f | mean %.2f s per paired draw",
                       info$n_ok, length(ids), res_dt[status == "ok", mean(INB_vsl_d)], res_dt[status == "ok", mean(INB_vsl_d > 0)], mean(sec)))
  }
  m <- list(); m[[paste0("run_", stage)]] <- info
  do.call(voi_manifest_update, c(list(paths = p, cfg = cfg), m))
  invisible(res_dt)
}

voi_run_psa <- function(cfg, stage) for (ct in cfg$countries) voi_run_psa_country(cfg, ct, stage)

.voi_stage <- get0("VOI_CURRENT_STAGE", ifnotfound = "")
if (.voi_stage %in% c("pilot", "psa")) {
  voi_run_psa(VOI_CFG, .voi_stage)
} else if (!exists("VOI_CURRENT_STAGE")) {          # standalone: pilot, then the main PSA
  voi_run_psa(VOI_CFG, "pilot"); voi_run_psa(VOI_CFG, "psa")
}
rm(.voi_stage)
