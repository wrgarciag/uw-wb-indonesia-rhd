# ==============================================================================
# SAP VOI analysis — 04: information-value estimates, summaries and figures
# scripts/voi/04_analyse_and_plot.R
#
# Decision: ref vs sap. Incremental net benefit per draw (constant 2025 US$, discounted to
# 2026, national population, 2026-2100), from the saved per-draw totals (no projection is
# re-run here):
#   VSL family : INB(k)      = k x B - C   (B = sum_t df_t VSL_t x deaths averted_t)
#   DALY family: INB(lambda) = lambda x D - C   (D = sum_t df_t deaths averted_t x 30: PROXY)
#   C = sum_t df_t (cost_sap_t - cost_ref_t)
# Because the four SAP inputs are the ONLY uncertain inputs and the model is deterministic
# given them, joint SAP EVPPI = EVPI conditional on the fixed model assumptions, computed
# DIRECTLY from the simulated net benefits: EVPI = mean(max(0, INB)) - max(0, mean(INB)).
#
# Single-parameter EVPPI (learn one input; the others integrated over their CONDITIONAL
# joint prior, so the M_AB / M_BC dependence is preserved — never a one-way analysis):
#   continuous (M_AB, M_BC on the log scale; partial protection on its own scale):
#     single-input GAM (mgcv, thin-plate, REML) of each NB component (B, D, C) on the input;
#     g(theta) = a x E[B or D | theta] - E[C | theta] (exact linear combination, so one set
#     of fits serves every threshold); EVPPI = mean(max(0, g)) - max(0, mean(g)).
#   categorical (protection age 20 / 39 / 100): STRATIFIED estimate from category means and
#     the sample category frequencies (Monte Carlo estimate, not an exact calculation).
# Precision: bootstrap of COMPLETE paired draw rows (VOI_N_BOOT replicates; every estimator
# refitted). Checks: GAM basis adequacy (k-index; k doubled once if inadequate), R^2,
# binned-estimator cross-check, split-half agreement, running estimates (convergence),
# 0 <= EVPPI_i <= EVPI (apparent violations investigated with the bootstrap, never clipped).
# Individual EVPPIs are NOT additive: never summed or stacked.
# Exploratory prior-scenario sensitivity: the regime weight w_trial re-weighted to 0.25 /
# 0.75 by importance weights (joint quantities only).
#
# Writes (per country, tables/ and figures/): inb_summary.csv, psa_outcome_summary.csv,
#   decision_curves.csv, evppi_comparison.csv (reference thresholds; tidy, ranked),
#   evppi_curves.csv (all thresholds), evppi_bootstrap_ranks.csv, gam_diagnostics.csv,
#   evppi_binned_check.csv, evppi_split_half.csv, age_category_diagnostics.csv,
#   convergence.csv, prior_scenario_reweighting.csv, mc_limitations.csv;
#   figures evppi_comparison.png, evppi_threshold_curves.png, evppi_ranked_reference.png,
#   decision_uncertainty.png, gam_fits.png, convergence.png.
#
# CONFERENCE PATH (2026-2050; section at the end of this file): recomputes every economic total,
# BCR, decision probability and VOI estimate for 2026-2050 INCLUSIVE from the saved annual paired
# arm-level results (no projection, calibration or PSA is re-run), adds BCR summaries, the
# cost-reweighted BCR-hurdle objective G_r = B - r C, and normalises VOI by gross programme cost,
# GDP and domestic general government health expenditure (GGHE-D). Writes ONLY to
# output/voi/<run_id>/<country>/conference_2026_2050/ and output/voi/<run_id>/conference_2026_2050/;
# the 2026-2100 outputs above are not touched. reports/slideck-voi.Rmd reads those files only.
# Select the path with VOI_ANALYSIS_PATH before sourcing: "full" (default; the 2026-2100 analysis
# above), "conference" (2026-2050 path only) or "both". Conference-only command (repository root):
#   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" -e "VOI_ANALYSIS_PATH <- 'conference'; source('scripts/voi/04_analyse_and_plot.R')"
# ==============================================================================

if (!exists("VOI_CFG")) { VOI_CONFIG_ONLY <- TRUE; source(file.path("scripts", "voi", "00_run_voi.R")); rm(VOI_CONFIG_ONLY) }
library(data.table)
suppressPackageStartupMessages({ library(mgcv); library(ggplot2); library(patchwork) })

VOI_PARS <- c("M_AB", "M_BC", "sap_protection_age_max", "partial_adherence_relative_protection")
VOI_PAR_TYPE <- c(M_AB = "log", M_BC = "log", sap_protection_age_max = "category", partial_adherence_relative_protection = "identity")
VOI_PAR_LABEL <- c(M_AB = "M_AB (A->B multiplier)", M_BC = "M_BC (B->C multiplier)",
                   sap_protection_age_max = "Protection age cap", partial_adherence_relative_protection = "Partial-adherence protection",
                   joint_SAP = "Joint SAP (= EVPI)")
# reference categorical slots 1-4 (dataviz reference palette, fixed order); joint = primary ink
VOI_PAR_COL <- c(M_AB = "#2a78d6", M_BC = "#eb6834", sap_protection_age_max = "#1baf7a",
                 partial_adherence_relative_protection = "#eda100", joint_SAP = "#0b0b0b")
VOI_PAR_SHAPE <- c(M_AB = 16, M_BC = 17, sap_protection_age_max = 15, partial_adherence_relative_protection = 18, joint_SAP = 23)
VOI_NOT_ADDITIVE <- paste("Individual EVPPIs are NOT additive: they overlap through correlation (M_AB, M_BC) and interactions,",
                          "so they must not be summed or stacked. Each is the value of learning that input alone.")

# ------------------------------------------------------------------------------
# estimators
# ------------------------------------------------------------------------------
voi_thresholds <- function(cfg, gni) {
  near <- function(a, b) vapply(a, function(z) any(abs(z - b) < 1e-10), logical(1))
  v <- data.table(family = "VSL multiple", x = cfg$vsl_mult_grid, gni_mult = NA_real_, is_ref = near(cfg$vsl_mult_grid, cfg$vsl_mult_ref))
  v[, `:=`(unit = "multiple k of the 08 reference-case VSL schedule (k = 1: production valuation)",
           label = paste0("VSL x ", format(x, drop0trailing = TRUE, scientific = FALSE, trim = TRUE)))]
  v[abs(x - 1) < 1e-10, label := "VSL x 1 (reference case)"]
  d <- data.table(family = "US$ per DALY (proxy)", x = cfg$daly_gni_grid * gni, gni_mult = cfg$daly_gni_grid,
                  is_ref = near(cfg$daly_gni_grid, cfg$daly_gni_ref))
  d[, `:=`(unit = "constant 2025 US$ per DALY averted (deaths x 30 proxy)",
           label = sprintf("%s x GNI pc = US$%s/DALY", format(gni_mult, drop0trailing = TRUE, trim = TRUE), formatC(round(x), format = "d", big.mark = ",")))]
  th <- rbind(v, d)
  th[, tid := .I]
  th[]
}
voi_nb <- function(B, D, C, th) {                 # n x T matrix of INB
  # "VSL multiple": k B - C | "US$ per DALY (proxy)": lambda D - C | "BCR hurdle" (conference path only):
  # G_r = B - r C, a cost-reweighted objective (not a welfare measure). Element-wise identical to the
  # original two-family form for the first two families.
  ish <- th$family == "BCR hurdle"
  isv <- th$family == "VSL multiple"
  isd <- !isv & !ish
  m <- matrix(0, length(B), nrow(th))
  if (any(isv)) m[, isv] <- outer(B, th$x[isv]) - C
  if (any(isd)) m[, isd] <- outer(D, th$x[isd]) - C
  if (any(ish)) m[, ish] <- B - outer(C, th$x[ish])
  m
}
voi_evpi_cols <- function(m, w = NULL) {
  if (is.null(w)) return(colMeans(pmax(m, 0)) - pmax(colMeans(m), 0))
  w <- w / sum(w)
  colSums(w * pmax(m, 0)) - pmax(colSums(w * m), 0)
}
voi_fit_gam <- function(y, x, k) mgcv::gam(y ~ s(x, k = k, bs = "tp"), method = "REML")
voi_xform <- function(v, type) if (type == "log") log(v) else as.numeric(v)
# conditional means of the NB components given ONE input (others integrated over their conditional prior)
voi_cond <- function(dat, par, type, kf = NULL, diag = FALSE) {
  if (type == "category") {
    g <- dat[[par]]
    return(list(fit = list(B = stats::ave(dat$B, g), D = stats::ave(dat$D, g), C = stats::ave(dat$C, g)), k = NULL, diag = NULL))
  }
  x <- voi_xform(dat[[par]], type)
  fit <- list(); kk <- list(); dg <- list()
  for (cmp in c("B", "D", "C")) {
    k <- if (is.null(kf)) 10L else kf[[cmp]]
    m <- voi_fit_gam(dat[[cmp]], x, k); n_ref <- 0L
    if (diag) {
      set.seed(1L); kc <- mgcv::k.check(m)
      while (is.null(kf) && kc[1, "p-value"] < 0.05 && kc[1, "edf"] > 0.8 * (k - 1) && k < 40L) {
        k <- 2L * k; n_ref <- n_ref + 1L; m <- voi_fit_gam(dat[[cmp]], x, k); set.seed(1L); kc <- mgcv::k.check(m)
      }
      sm <- summary(m)
      dg[[cmp]] <- data.table(parameter = par, component = cmp, scale = type, k_final = k, k_refits = n_ref,
                              edf = kc[1, "edf"], k_index = kc[1, "k-index"], k_index_p = kc[1, "p-value"],
                              r_sq = sm$r.sq, dev_expl = sm$dev.expl, resid_sd = stats::sd(stats::residuals(m)),
                              fitted_mean_minus_y_mean = mean(fitted(m)) - mean(dat[[cmp]]))
    }
    fit[[cmp]] <- as.numeric(fitted(m)); kk[[cmp]] <- k
  }
  list(fit = fit, k = kk, diag = if (diag) rbindlist(dg) else NULL)
}
voi_evppi_from_fit <- function(fit, th, w = NULL) voi_evpi_cols(voi_nb(fit$B, fit$D, fit$C, th), w)
# all estimates (joint + 4 individual) for a data subset, with fixed basis sizes
voi_all_estimates <- function(d, th, kfs) {
  out <- matrix(NA_real_, length(VOI_PARS) + 1L, nrow(th), dimnames = list(c("joint_SAP", VOI_PARS), NULL))
  out["joint_SAP", ] <- voi_evpi_cols(voi_nb(d$B, d$D, d$C, th))
  for (pp in VOI_PARS) {
    f <- tryCatch(voi_cond(d, pp, VOI_PAR_TYPE[[pp]], kfs[[pp]])$fit, error = function(e) NULL)
    if (!is.null(f)) out[pp, ] <- voi_evppi_from_fit(f, th)
  }
  out
}
voi_boot_block <- function(bids, idx, dat, th, kfs) {
  arr <- array(NA_real_, c(length(bids), length(VOI_PARS) + 1L, nrow(th)))
  for (j in seq_along(bids)) arr[j, , ] <- voi_all_estimates(dat[idx[, bids[j]]], th, kfs)
  arr
}
# binned (equal-count) stratified cross-check for continuous inputs
voi_binned <- function(dat, par, th, nb = 10L) {
  br <- unique(stats::quantile(dat[[par]], seq(0, 1, length.out = nb + 1L), names = FALSE))
  g <- cut(dat[[par]], br, include.lowest = TRUE)
  voi_evppi_from_fit(list(B = stats::ave(dat$B, g), D = stats::ave(dat$D, g), C = stats::ave(dat$C, g)), th)
}

# ------------------------------------------------------------------------------
# figures
# ------------------------------------------------------------------------------
voi_theme <- function() theme_minimal(base_size = 10) + theme(
  plot.background = element_rect(fill = "#fcfcfb", colour = NA), panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3), axis.line = element_line(colour = "#c3c2b7", linewidth = 0.3),
  axis.text = element_text(colour = "#52514e"), axis.title = element_text(colour = "#52514e", size = 9),
  plot.title = element_text(face = "bold", size = 10.5, colour = "#0b0b0b"), plot.subtitle = element_text(colour = "#52514e", size = 8.5),
  plot.caption = element_text(colour = "#898781", size = 7.5, hjust = 0), legend.position = "bottom",
  legend.text = element_text(size = 8.5, colour = "#0b0b0b"), strip.text = element_text(face = "bold", colour = "#0b0b0b", size = 8.5))
voi_ggsave <- function(fig, file, cfg, w, h) { voi_assert_out(file, cfg); ggsave(file, fig, width = w, height = h, dpi = 150, bg = "#fcfcfb") }
voi_direct_labels <- function(dt, xvar, yvar) {           # label each series once, at its maximum
  dt <- copy(dt)
  short <- c(M_AB = "M_AB", M_BC = "M_BC", sap_protection_age_max = "Age cap", partial_adherence_relative_protection = "Partial adh.",
             joint_SAP = "Joint (EVPI)")
  dt[, is_peak := seq_len(.N) == which.max(get(yvar)) & max(get(yvar)) > 0, by = parameter]
  dt[, txt := fifelse(is_peak, short[as.character(parameter)], "")]  # empty labels: curve points the labels repel from
  if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(data = dt, aes(x = .data[[xvar]], y = .data[[yvar]], label = txt), colour = "#0b0b0b", size = 2.8,
                             min.segment.length = 0, segment.colour = "#898781", segment.size = 0.3, box.padding = 0.5, point.padding = 0.2,
                             force = 2, max.overlaps = Inf, seed = 1, inherit.aes = FALSE)
  else geom_text(data = dt[is_peak == TRUE], aes(x = .data[[xvar]], y = .data[[yvar]], label = txt), colour = "#0b0b0b", size = 2.8, vjust = -0.6, inherit.aes = FALSE)
}

# ------------------------------------------------------------------------------
# per-country analysis
# ------------------------------------------------------------------------------
voi_analyse_country <- function(cfg, country) {
  p <- voi_paths(cfg, country)
  inp <- voi_load_inputs(cfg, country)
  f_res <- file.path(p$psa, "psa_draw_results.rds")
  if (!file.exists(f_res)) stop("[", country, "] no PSA results; run stage 'psa' first.", call. = FALSE)
  res <- readRDS(f_res)
  n_fail <- res[status != "ok", .N]
  dat <- res[status == "ok"][order(draw_id)]
  dat[, `:=`(B = B_d, D = DALY_d, C = C_d)]
  N <- nrow(dat); NB <- cfg$n_boot
  gni <- inp$econ$gni_pc_ref_usd2025; N0 <- inp$econ$N0
  th <- voi_thresholds(cfg, gni)
  ref <- th[is_ref == TRUE]
  point <- fread(file.path(p$validation, "point_run_economics.csv"))
  voi_log(p, sprintf("04 analysis: %d paired draws (%d excluded) | %d thresholds | GNI pc reference US$%.0f | bootstrap %d", N, n_fail, nrow(th), gni, NB))

  # ---------- joint: INB summaries, EVPI (direct) ----------
  nb <- voi_nb(dat$B, dat$D, dat$C, th)
  nb_pt <- as.numeric(voi_nb(point$B_d, point$DALY_d, point$C_d, th))
  einb <- colMeans(nb); sdinb <- apply(nb, 2, stats::sd); ppos <- colMeans(nb > 0)
  eol_ref <- colMeans(pmax(nb, 0)); eol_sap <- colMeans(pmax(-nb, 0))
  evpi <- voi_evpi_cols(nb)

  # ---------- single-parameter fits on the full sample (diagnostics, basis sizes) ----------
  fits <- list(); kfs <- list(); gdiag <- list()
  for (pp in VOI_PARS) {
    cf <- voi_cond(dat, pp, VOI_PAR_TYPE[[pp]], diag = TRUE)
    fits[[pp]] <- cf$fit; kfs[[pp]] <- cf$k; gdiag[[pp]] <- cf$diag
  }
  est <- voi_all_estimates(dat, th, kfs)
  stopifnot(max(abs(est["joint_SAP", ] - evpi)) == 0)

  # ---------- bootstrap of complete paired rows ----------
  set.seed(cfg$boot_seed)
  idx <- matrix(sample.int(N, N * NB, replace = TRUE), N, NB)
  blocks <- split(seq_len(NB), cut(seq_len(NB), max(1L, min(cfg$max_cores, parallel::detectCores() - 1L)), labels = FALSE))
  t0 <- proc.time()[["elapsed"]]
  cl <- parallel::makeCluster(length(blocks))
  barr <- tryCatch({
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({ library(data.table); library(mgcv) }); setDTthreads(1L); NULL })
    parallel::clusterExport(cl, c("VOI_PARS", "VOI_PAR_TYPE", "voi_nb", "voi_evpi_cols", "voi_fit_gam", "voi_xform", "voi_cond",
                                  "voi_evppi_from_fit", "voi_all_estimates", "voi_boot_block"), envir = globalenv())
    dsmall <- dat[, c("draw_id", VOI_PARS, "B", "D", "C"), with = FALSE]
    parts <- parallel::parLapply(cl, blocks, voi_boot_block, idx = idx, dat = dsmall, th = th, kfs = kfs)
    a <- array(NA_real_, c(NB, length(VOI_PARS) + 1L, nrow(th)))
    for (i in seq_along(blocks)) a[blocks[[i]], , ] <- parts[[i]]
    a
  }, finally = parallel::stopCluster(cl))
  dimnames(barr) <- list(NULL, c("joint_SAP", VOI_PARS), NULL)
  n_boot_fail <- sum(is.na(barr[, VOI_PARS, 1]))
  voi_log(p, sprintf("04 bootstrap: %d replicates x %d estimators x %d thresholds in %.1f s (%d failed fits)",
                     NB, length(VOI_PARS) + 1L, nrow(th), proc.time()[["elapsed"]] - t0, n_boot_fail))

  # ---------- split-half, binned check, convergence ----------
  half <- dat$draw_id <= stats::median(dat$draw_id)
  est_h1 <- voi_all_estimates(dat[half], th, kfs); est_h2 <- voi_all_estimates(dat[!half], th, kfs)
  binned <- rbindlist(lapply(c("M_AB", "M_BC", "partial_adherence_relative_protection"), function(pp)
    data.table(parameter = pp, tid = th$tid, evppi_binned10 = voi_binned(dat, pp, th), evppi_gam = est[pp, ])))[tid %in% ref$tid]
  binned <- merge(th[, .(tid, family, label)], binned, by = "tid")
  conv_n <- sort(unique(c(seq(100L, N, by = 100L), N)))
  conv <- rbindlist(lapply(conv_n, function(n) {
    d <- dat[seq_len(n)]; e <- voi_all_estimates(d, th[is_ref == TRUE], kfs); m <- voi_nb(d$B, d$D, d$C, th[is_ref == TRUE])
    rbind(data.table(n = n, tid = ref$tid, parameter = rep(rownames(e), each = ncol(e)), value = as.numeric(t(e)), quantity = "EVPPI / EVPI"),
          data.table(n = n, tid = ref$tid, parameter = "INB", value = colMeans(m), quantity = "E[INB]"),
          data.table(n = n, tid = ref$tid, parameter = "INB", value = colMeans(m > 0), quantity = "P(INB > 0)"))
  }))
  conv <- merge(th[, .(tid, family, label)], conv, by = "tid")
  setorder(conv, tid, quantity, parameter, n)

  # ---------- tidy EVPPI table (all thresholds) ----------
  bq <- function(a, pr) apply(a, 2, stats::quantile, probs = pr, na.rm = TRUE, names = FALSE)
  curves <- rbindlist(lapply(c("joint_SAP", VOI_PARS), function(pp) {
    b <- barr[, pp, ]
    data.table(tid = th$tid, parameter = pp, estimate = est[pp, ], boot_mean = colMeans(b, na.rm = TRUE),
               boot_se = apply(b, 2, stats::sd, na.rm = TRUE), ci_lo = bq(b, 0.025), ci_hi = bq(b, 0.975),
               split_half_1 = est_h1[pp, ], split_half_2 = est_h2[pp, ],
               p_boot_exceeds_joint = if (pp == "joint_SAP") NA_real_ else colMeans(b > barr[, "joint_SAP", ] + 1e-9 * pmax(1, abs(barr[, "joint_SAP", ])), na.rm = TRUE),
               diff_joint_ci_lo = if (pp == "joint_SAP") NA_real_ else bq(b - barr[, "joint_SAP", ], 0.025))
  }))
  curves <- merge(th, curves, by = "tid")
  curves[, `:=`(evpi = evpi[tid], ratio_to_joint = fifelse(evpi[tid] > 0, estimate / evpi[tid], NA_real_),
                rel_se = fifelse(estimate > 0, boot_se / estimate, NA_real_),
                estimate_per_capita = estimate / N0, n_draws = N, n_boot = NB,
                method = fcase(parameter == "joint_SAP", "direct: EVPI = mean(max(0, INB)) - max(0, mean(INB)) (joint SAP EVPPI)",
                               parameter %in% c("M_AB", "M_BC"), "single-input GAM (thin-plate, REML) on log(input); NB components B, D, C fitted separately",
                               parameter == "partial_adherence_relative_protection", "single-input GAM (thin-plate, REML); NB components fitted separately",
                               default = "stratified: category means x sample category frequencies"))]
  curves[, parameter_label := VOI_PAR_LABEL[parameter]]
  curves[, exceeds_joint := parameter != "joint_SAP" & estimate > evpi * (1 + 1e-9) + 1e-6]
  curves[, bound_check := fcase(parameter == "joint_SAP", "reference (joint)",
                                estimate < 0, "NEGATIVE: investigate",
                                !exceeds_joint, "0 <= EVPPI <= joint: OK",
                                exceeds_joint & diff_joint_ci_lo <= 0, "exceeds joint, within bootstrap error",
                                default = "BOUND VIOLATION: exceeds joint beyond bootstrap error (investigate)")]
  curves[, status := fcase(evpi == 0, "EVPI = 0: no decision uncertainty in the sample (every draw favours one option)",
                           ci_hi < 0.05 * evpi, "negligible: bootstrap 97.5% bound < 5% of the joint value",
                           parameter != "joint_SAP" & estimate == 0,
                           "EXPLORATORY: point estimate 0 (learning this input alone never reverses the decision in the fitted conditional mean) but bootstrap replicates are positive; more draws needed",
                           rel_se <= 0.20, "stable at this sample size (bootstrap relative SE <= 0.20)",
                           default = "EXPLORATORY: unstable at this sample size (bootstrap relative SE > 0.20); more draws needed")]
  setorder(curves, family, x, parameter)

  # ---------- reference thresholds: ranks + bootstrap rank probabilities ----------
  cmpt <- curves[is_ref == TRUE]
  cmpt[parameter != "joint_SAP", rank_among_individual := frank(-estimate, ties.method = "min"), by = tid]
  cmpt[parameter != "joint_SAP", rank_note := fifelse(estimate == 0, "tied at 0 (rank not informative)",
                                                      fifelse(duplicated(rank_among_individual) | duplicated(rank_among_individual, fromLast = TRUE), "tied", "")), by = tid]
  rk <- rbindlist(lapply(ref$tid, function(t) {
    b <- barr[, VOI_PARS, t, drop = FALSE][, , 1]
    r <- t(apply(-b, 1, rank, ties.method = "min"))
    data.table(tid = t, parameter = VOI_PARS, p_rank_first = colMeans(r == 1, na.rm = TRUE), median_boot_rank = apply(r, 2, stats::median, na.rm = TRUE))
  }))
  cmpt <- merge(cmpt, rk, by = c("tid", "parameter"), all.x = TRUE)
  cmpt[, country := country]
  setcolorder(cmpt, c("country", "family", "label", "x", "gni_mult", "unit", "parameter", "parameter_label", "estimate", "ci_lo", "ci_hi", "boot_se", "rel_se",
                      "evpi", "ratio_to_joint", "rank_among_individual", "rank_note", "p_rank_first", "median_boot_rank", "bound_check", "status", "estimate_per_capita",
                      "split_half_1", "split_half_2", "method", "n_draws", "n_boot"))
  setorder(cmpt, family, x, rank_among_individual, na.last = FALSE)

  # ---------- INB / decision summaries ----------
  dec <- data.table(tid = th$tid, E_INB = einb, mc_se_E_INB = sdinb / sqrt(N), sd_INB = sdinb, P_INB_pos = ppos,
                    mc_se_P = sqrt(ppos * (1 - ppos) / N), EOL_ref = eol_ref, EOL_sap = eol_sap, EVPI = evpi,
                    INB_point_values = nb_pt, optimal_current_info = ifelse(einb > 0, "sap (scale-up)", "ref"))
  dec <- merge(th, dec, by = "tid")
  dec <- merge(dec, curves[parameter == "joint_SAP", .(tid, EVPI_ci_lo = ci_lo, EVPI_ci_hi = ci_hi, EVPI_boot_se = boot_se)], by = "tid")
  dec[, `:=`(E_INB_ci_lo = E_INB - 1.96 * mc_se_E_INB, E_INB_ci_hi = E_INB + 1.96 * mc_se_E_INB,
             P_INB_pos_note = fifelse(P_INB_pos %in% c(0, 1), sprintf("no draw on the other side; rule-of-three 95%% bound on the other side ~ %.4f", 3 / N), ""),
             EVPI_per_capita = EVPI / N0, country = country)]
  inb_ref <- dec[is_ref == TRUE]
  setcolorder(inb_ref, c("country", "family", "label", "x", "E_INB", "E_INB_ci_lo", "E_INB_ci_hi", "mc_se_E_INB", "P_INB_pos", "mc_se_P", "P_INB_pos_note",
                         "EOL_ref", "EOL_sap", "EVPI", "EVPI_ci_lo", "EVPI_ci_hi", "EVPI_per_capita", "optimal_current_info", "INB_point_values"))
  q <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
  outc <- rbindlist(lapply(list(c("C_d", "incremental cost, discounted (US$)"), c("B_d", "VSL benefit, discounted (US$)"),
                                c("INB_vsl_d", "INB at VSL x 1, discounted (US$)"), c("INB_vsl_u", "INB at VSL x 1, undiscounted (US$)"),
                                c("DA", "RHD deaths averted, 2026-2100 (undiscounted)"), c("DALY_d", "DALYs averted proxy (discounted deaths x 30)"),
                                c("bcr_d", "benefit-cost ratio, discounted"), c("cost_per_daly", "cost per DALY averted proxy (US$)")), function(v) {
    x <- dat[[v[1]]]; qq <- q(x)
    data.table(country = country, outcome = v[1], description = v[2], mean = mean(x), mc_se = stats::sd(x) / sqrt(N), sd = stats::sd(x),
               q025 = qq[1], median = qq[2], q975 = qq[3], production_point_value = point[[v[1]]])
  }))

  # ---------- age-category diagnostics ----------
  agecat <- rbindlist(lapply(ref$tid, function(t) {
    m <- nb[, t]
    out <- dat[, .(idx = .I), by = sap_protection_age_max][, .(n = .N, mean_INB = mean(m[idx]), se_mean_INB = stats::sd(m[idx]) / sqrt(.N)), by = sap_protection_age_max]
    out[, `:=`(share = n / N, tid = t, prior_prob = cfg$prior$protection_age$probs[match(sap_protection_age_max, cfg$prior$protection_age$values)])]
  }))
  agecat <- merge(th[, .(tid, family, label)], agecat, by = "tid")

  # ---------- exploratory prior-scenario re-weighting (regime weight) ----------
  w0 <- cfg$prior$w_trial; z <- dat$regime == "trial_efficacy"
  rew <- rbindlist(lapply(cfg$prior_reweight_w, function(w1) {
    w <- ifelse(z, w1 / w0, (1 - w1) / (1 - w0)); wn <- w / sum(w)
    m <- nb[, ref$tid, drop = FALSE]
    data.table(tid = ref$tid, w_trial = w1, ess = sum(w)^2 / sum(w^2), E_INB = colSums(wn * m), P_INB_pos = colSums(wn * (m > 0)),
               EVPI = voi_evpi_cols(m, w))
  }))
  rew <- merge(th[, .(tid, family, label)], rew, by = "tid")[, note := "importance re-weighting of the same draws; joint quantities only; exploratory"]

  # ---------- Monte Carlo limitations ----------
  mcl <- data.table(country = country, n_draws = N, n_excluded = n_fail, n_boot = NB, n_boot_failed_fits = n_boot_fail,
                    min_category_n = min(agecat$n), regime_trial_n = sum(z), regime_recurrence_n = sum(!z),
                    n_reference_estimates = cmpt[parameter != "joint_SAP", .N],
                    n_exploratory = cmpt[grepl("^EXPLORATORY", status), .N], n_stable = cmpt[grepl("^stable", status), .N],
                    n_negligible = cmpt[grepl("^negligible", status), .N], n_evpi_zero = cmpt[grepl("^EVPI = 0", status) & parameter != "joint_SAP", .N],
                    n_exceeds_joint = cmpt[exceeds_joint == TRUE, .N], n_bound_violation = cmpt[grepl("^BOUND", bound_check), .N],
                    max_rel_change_E_INB_400_to_500 = conv[quantity == "E[INB]" & n %in% c(max(conv_n[conv_n < N]), N),
                                                           abs(diff(value)) / abs(value[2]), by = tid][, max(V1)],
                    note = "N = 500 is the initial sample; no automatic increase. Estimates flagged EXPLORATORY need more draws (extend by raising N_DRAWS under the same run_id).")

  # ---------- write tables ----------
  tw <- function(x, f) voi_fwrite(x, file.path(p$tables, f), cfg)
  tw(inb_ref, "inb_summary.csv"); tw(outc, "psa_outcome_summary.csv")
  tw(dec[, c("country", setdiff(names(dec), "country")), with = FALSE], "decision_curves.csv")
  tw(cmpt, "evppi_comparison.csv"); tw(cbind(country = country, curves), "evppi_curves.csv")
  tw(cbind(country = country, rk), "evppi_bootstrap_ranks.csv"); tw(cbind(country = country, rbindlist(gdiag)), "gam_diagnostics.csv")
  tw(cbind(country = country, binned), "evppi_binned_check.csv")
  tw(cbind(country = country, curves[is_ref == TRUE, .(family, label, parameter, estimate, split_half_1, split_half_2)]), "evppi_split_half.csv")
  tw(cbind(country = country, agecat), "age_category_diagnostics.csv"); tw(cbind(country = country, conv), "convergence.csv")
  tw(cbind(country = country, rew), "prior_scenario_reweighting.csv"); tw(mcl, "mc_limitations.csv")
  voi_save(list(estimates = est, bootstrap = barr, thresholds = th, kfs = kfs, boot_index_seed = cfg$boot_seed),
           file.path(p$tables, "evppi_bootstrap_array.rds"), cfg)

  # ---------- figures ----------
  sc_m <- function(x) x / 1e6
  ylab_m <- "EVPPI (US$ million, constant 2025, discounted to 2026; national, 2026-2100)"
  cols <- VOI_PAR_COL; labs_p <- VOI_PAR_LABEL[names(cols)]
  cv <- copy(curves); cv[, parameter := factor(parameter, levels = c("joint_SAP", VOI_PARS))]
  refl <- function(fam) ref[family == fam]
  mk_curve <- function(fam, logx) {
    d <- cv[family == fam]
    g <- ggplot(d, aes(x, sc_m(estimate), colour = parameter, fill = parameter)) +
      geom_vline(data = refl(fam), aes(xintercept = x), colour = "#898781", linewidth = 0.3, linetype = "dashed") +
      geom_ribbon(aes(ymin = sc_m(ci_lo), ymax = sc_m(ci_hi)), alpha = 0.10, colour = NA) +
      geom_line(aes(linewidth = parameter == "joint_SAP")) +
      scale_linewidth_manual(values = c(`FALSE` = 0.7, `TRUE` = 1.0), guide = "none") +
      scale_colour_manual(values = cols, labels = labs_p, name = NULL) + scale_fill_manual(values = cols, labels = labs_p, name = NULL) +
      voi_direct_labels(d[, .(parameter, x, y = sc_m(estimate))], "x", "y") +
      labs(y = ylab_m) + voi_theme()
    if (logx) g + scale_x_log10(breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2), labels = function(b) format(b, drop0trailing = TRUE, scientific = FALSE)) +
      labs(x = "VSL multiple k (k = 1: 08 reference-case VSL schedule; log scale)", title = "A. Valuation by VSL (deaths averted x k x VSL_t)",
           subtitle = "Dashed lines: reference thresholds k = 0.05, 0.10, 1")
    else g + scale_x_continuous(labels = function(b) paste0("$", formatC(b, format = "d", big.mark = ","))) +
      labs(x = "Threshold, US$ per DALY averted (proxy: discounted deaths averted x 30)",
           title = "B. Valuation per DALY (US$/DALY threshold)",
           subtitle = sprintf("Dashed lines: 0.5 and 1 x GNI per capita (US$%s, 2026, from 08's income path)", formatC(round(gni), format = "d", big.mark = ",")))
  }
  gA <- mk_curve("VSL multiple", TRUE); gB <- mk_curve("US$ per DALY (proxy)", FALSE)
  rk_d <- copy(cmpt)
  rk_d[, panel := factor(label, levels = unique(cmpt[order(family == "US$ per DALY (proxy)", x), label]))]
  rk_d[, ordkey := fifelse(parameter == "joint_SAP", Inf, estimate)]
  rk_d[, yy := factor(paste(panel, parameter, sep = "___"), levels = unique(paste(panel, parameter, sep = "___")[order(panel, ordkey)]))]
  rk_d[, lab_txt := fcase(parameter == "joint_SAP" & estimate == 0, sprintf("0 (INB > 0 in all %d draws)", N),
                          parameter == "joint_SAP", sprintf("%.1f", sc_m(estimate)),
                          estimate == 0, "0",
                          default = sprintf("#%d  %.1f", rank_among_individual, sc_m(estimate)))]
  rk_d[, flag := fifelse(grepl("^EXPLORATORY", status), " (exploratory)", "")]
  gC <- ggplot(rk_d, aes(sc_m(estimate), yy, colour = parameter, shape = parameter)) +
    geom_linerange(aes(xmin = sc_m(ci_lo), xmax = sc_m(ci_hi)), linewidth = 0.6) +
    geom_point(aes(fill = parameter), size = 2.4) +
    geom_text(aes(x = sc_m(ci_hi), label = paste0(lab_txt, flag)), hjust = -0.12, size = 2.5, colour = "#52514e") +
    facet_wrap(~panel, scales = "free", ncol = 3) +
    scale_y_discrete(labels = function(s) VOI_PAR_LABEL[sub("^.*___", "", s)]) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.45))) +
    scale_colour_manual(values = cols, labels = labs_p, name = NULL) + scale_fill_manual(values = cols, labels = labs_p, name = NULL) +
    scale_shape_manual(values = VOI_PAR_SHAPE, labels = labs_p, name = NULL) +
    labs(x = "US$ million (point estimate and bootstrap 95% interval)", y = NULL,
         title = "C. Ranked comparison at the reference thresholds (#rank among the four inputs; joint value shown separately, not a sum)") +
    voi_theme() + theme(panel.grid.major.y = element_blank()) + guides(colour = "none", fill = "none", shape = "none")
  cap <- paste0(VOI_NOT_ADDITIVE, "\nJoint SAP EVPPI = EVPI conditional on the fixed calibration, initial state, care cascade, surgery and unit costs. ",
                sprintf("N = %d paired draws; %d bootstrap replicates. Provisional, not expert-agreed priors. Upper bound on the value of research; not EVSI.", N, NB),
                "\nModel baseline: calibration fit INADEQUATE (outside GBD diagnostic tolerances; production meta$calibration_adequacy), so all values are preliminary.")
  if (!grepl("^INADEQUATE", inp$calibration_adequacy)) cap <- sub("calibration fit INADEQUATE.*$", paste0("calibration status: ", substr(inp$calibration_adequacy, 1, 80)), cap)
  cap <- paste(unlist(lapply(strsplit(cap, "\n", fixed = TRUE)[[1]], strwrap, width = 175)), collapse = "\n")
  # caption on the bottom panel (left-aligned by voi_theme); patchwork annotation themes are not applied with this ggplot2 version
  fig_cmp <- (gA | gB) / (gC + labs(caption = cap)) + plot_layout(heights = c(1, 1.15), guides = "collect") +
    plot_annotation(title = sprintf("%s: value of information for the four SAP inputs (ref vs scale-up)", country)) &
    theme(legend.position = "bottom")
  voi_ggsave(fig_cmp, file.path(p$figures, "evppi_comparison.png"), cfg, 12, 11)
  voi_ggsave(((gA + labs(caption = paste(strwrap(VOI_NOT_ADDITIVE, 95), collapse = "\n"))) | gB) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
             file.path(p$figures, "evppi_threshold_curves.png"), cfg, 12, 5.2)
  voi_ggsave(gC + labs(caption = VOI_NOT_ADDITIVE), file.path(p$figures, "evppi_ranked_reference.png"), cfg, 12, 6)

  # decision uncertainty: E[INB] and P(INB > 0) in separate panels (no dual axis)
  dd <- copy(dec)
  mk_dec <- function(fam, logx, what) {
    d <- dd[family == fam]
    g <- ggplot(d, aes(x)) + geom_vline(data = refl(fam), aes(xintercept = x), colour = "#898781", linewidth = 0.3, linetype = "dashed")
    g <- if (what == "E") g + geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.4) +
      geom_ribbon(aes(ymin = sc_m(E_INB_ci_lo), ymax = sc_m(E_INB_ci_hi)), fill = "#2a78d6", alpha = 0.15) +
      geom_line(aes(y = sc_m(E_INB)), colour = "#2a78d6", linewidth = 0.8) +
      geom_line(aes(y = sc_m(INB_point_values)), colour = "#52514e", linewidth = 0.5, linetype = "dotted") +
      labs(y = "E[INB], US$ million (95% MC interval)")
    else g + geom_line(aes(y = P_INB_pos), colour = "#2a78d6", linewidth = 0.8) + scale_y_continuous(limits = c(0, 1)) +
      labs(y = "P(INB > 0): scale-up preferred")
    g <- g + voi_theme()
    if (logx) g + scale_x_log10(breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2), labels = function(b) format(b, drop0trailing = TRUE, scientific = FALSE)) +
      labs(x = "VSL multiple k (log scale)") else g + scale_x_continuous(labels = function(b) paste0("$", formatC(b, format = "d", big.mark = ","))) + labs(x = "US$ per DALY (proxy)")
  }
  fig_dec <- (mk_dec("VSL multiple", TRUE, "E") + labs(title = "VSL valuation") | mk_dec("US$ per DALY (proxy)", FALSE, "E") + labs(title = "US$/DALY valuation (proxy)")) /
    (mk_dec("VSL multiple", TRUE, "P") | mk_dec("US$ per DALY (proxy)", FALSE, "P")) +
    plot_annotation(title = sprintf("%s: expected incremental net benefit and probability scale-up is preferred", country),
                    caption = sprintf("Blue: PSA mean over %d draws (ribbon: 95%% Monte Carlo interval). Dotted grey: deterministic production point values (M = 0.45, age 39, partial 0.5).", N))
  voi_ggsave(fig_dec, file.path(p$figures, "decision_uncertainty.png"), cfg, 11, 7.5)

  # GAM fits at the reference threshold with the largest EVPI
  tbest <- ref[which.max(evpi[ref$tid]), tid]
  fit_rows <- rbindlist(lapply(VOI_PARS, function(pp) {
    a <- th$x[tbest]; isv <- th$family[tbest] == "VSL multiple"
    g <- a * (if (isv) fits[[pp]]$B else fits[[pp]]$D) - fits[[pp]]$C
    data.table(parameter = pp, value = as.numeric(dat[[pp]]), inb = nb[, tbest], g = g)
  }))
  fit_rows[, parameter := factor(parameter, levels = VOI_PARS, labels = VOI_PAR_LABEL[VOI_PARS])]
  fig_fit <- ggplot(fit_rows, aes(value, sc_m(inb))) + geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.4) +
    geom_point(colour = "#898781", size = 0.7, alpha = 0.6) +
    geom_line(data = fit_rows[order(parameter, value)], aes(y = sc_m(g)), colour = "#2a78d6", linewidth = 0.8) +
    facet_wrap(~parameter, scales = "free_x") + labs(x = "Sampled value (M on log scale in the fit; protection age categorical)", y = "INB, US$ million",
      title = sprintf("%s: INB at %s vs each input, with the estimated conditional mean E[INB | input] (blue)", country, th$label[tbest]),
      caption = "Scatter around the blue curve is the variation due to the other three inputs (integrated over their conditional prior), not noise.") + voi_theme()
  voi_ggsave(fig_fit, file.path(p$figures, "gam_fits.png"), cfg, 10, 7)

  # convergence (running estimates at the reference thresholds)
  cvg <- conv[quantity == "EVPPI / EVPI"]
  cvg[, parameter := factor(parameter, levels = c("joint_SAP", VOI_PARS))]
  cvg[, panel := factor(label, levels = levels(rk_d$panel))]
  fig_conv <- ggplot(cvg, aes(n, sc_m(value), colour = parameter, shape = parameter)) + geom_line(linewidth = 0.6) + geom_point(size = 1.8, aes(fill = parameter)) +
    facet_wrap(~panel, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = cols, labels = labs_p, name = NULL) + scale_fill_manual(values = cols, labels = labs_p, name = NULL) +
    scale_shape_manual(values = VOI_PAR_SHAPE, labels = labs_p, name = NULL) +
    labs(x = "Number of paired draws used (draw ids 1..n)", y = "US$ million",
         title = sprintf("%s: running EVPPI / EVPI estimates (fixed GAM basis sizes)", country),
         caption = "Stabilisation across n indicates Monte Carlo adequacy; drift between n = 400 and 500 signals that more draws are needed.") + voi_theme()
  voi_ggsave(fig_conv, file.path(p$figures, "convergence.png"), cfg, 11, 6.5)

  voi_manifest_update(p, cfg, analysis = list(n_draws = N, n_excluded = n_fail, n_boot = NB, boot_seed = cfg$boot_seed,
                                              gni_pc_ref_usd2025 = gni, N0 = N0, reference_thresholds = ref$label,
                                              gam_basis = kfs, completed = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  for (i in seq_len(nrow(inb_ref)))
    voi_log(p, sprintf("04 %-38s E[INB] %10.1f m (MC se %.1f) | P(INB>0) %.3f | EVPI %8.2f m (95%% %.2f-%.2f)", inb_ref$label[i], inb_ref$E_INB[i] / 1e6,
                       inb_ref$mc_se_E_INB[i] / 1e6, inb_ref$P_INB_pos[i], inb_ref$EVPI[i] / 1e6, inb_ref$EVPI_ci_lo[i] / 1e6, inb_ref$EVPI_ci_hi[i] / 1e6))
  for (i in which(cmpt$parameter != "joint_SAP"))
    voi_log(p, sprintf("04   %-38s #%d %-38s %8.2f m (95%% %.2f-%.2f) | %s | %s", cmpt$label[i], cmpt$rank_among_individual[i], cmpt$parameter_label[i],
                       cmpt$estimate[i] / 1e6, cmpt$ci_lo[i] / 1e6, cmpt$ci_hi[i] / 1e6, cmpt$bound_check[i], substr(cmpt$status[i], 1, 40)))
  invisible(list(inb = inb_ref, evppi = cmpt))
}

voi_run_analysis <- function(cfg) for (ct in cfg$countries) voi_analyse_country(cfg, ct)

# ==============================================================================
# CONFERENCE PATH: 2026-2050 recomputation for the 15-minute VOI talk (reports/slideck-voi.Rmd)
#
# Horizon 2026-2050 INCLUSIVE. Benefits and costs accrued through 2050 only; effects after 2050
# are excluded (nothing after 2050 enters any total, ratio, probability or VOI estimate).
# Source data (read-only): psa/annual_summary.rds (draw x arm x year), draws/sap_draws.csv (draw
# parameters, joined by draw id), validation/point_run_annual.csv (deterministic point run), and the
# saved production economics via voi_load_inputs() (APPLIED 2025-US$ unit costs, saved VSL schedule,
# 3% discounting to 2026). No price conversion is repeated on the costs: the saved annual cost is
# checked against volumes x the saved applied unit costs.
#
# Per draw (s = sap scale-up arm, r = reference care, df_t = 1.03^-(t - 2026)):
#   C = sum_t df_t (cost_s,t - cost_r,t)            discounted incremental programme cost
#   B = sum_t df_t VSL_t (deaths_r,t - deaths_s,t)   discounted incremental monetised mortality benefit
#   D = sum_t df_t deaths averted_t x 30             DALY PROXY (08 convention; not modelled DALYs)
#   G = sum_t df_t cost_s,t                          discounted GROSS scale-up programme cost
#   BCR = B / C (incremental, discounted), defined only when C > 0; draws with C <= 0 are counted and
#   handled through net benefit (never by reversing inequalities or dropping them).
# Decision objectives: INB(k) = k B - C (k = 1: reference-case VSL; INB > 0 <=> BCR > 1 when C > 0),
#   lambda D - C (US$/DALY proxy) and, descriptively, the BCR hurdle G_r = B - r C = r [(1/r) B - C]
#   (r > 1 re-weights costs: a DIFFERENT objective; its VOI is not a welfare gain).
# VOI estimators, bootstrap, diagnostics: as in the 2026-2100 analysis above (same functions, same
# bootstrap seed, complete paired draw rows), applied to the 2026-2050 totals.
# Normalisation (scale comparisons only; NOT fiscal savings, research budgets or research BCRs):
#   VOI / E[G] (gross programme cost, same draws), VOI / one year of GDP (2023) and VOI / one year of
#   GGHE-D (2023), all in constant 2025 US$; bootstrap bounds divided by the SAME fixed denominator.
#   VOI / E[C] (incremental cost) is written separately and only when E[C] > 0.
# ==============================================================================

# Macro observations for the GDP / GGHE-D denominators. Retrieved 2026-11-10 (no network call at run
# time; 04 stops on a missing value). Common observed reference year 2023: the latest year WHO GHED
# reports GGHE-D for both countries (2024-2025 not yet published at retrieval); no forecast is used.
VOI_CONF_WDI_SRC <- "World Bank World Development Indicators (api.worldbank.org/v2, source=2), last updated 2026-10-08"
VOI_CONF_GHED_SRC <- "WHO Global Health Expenditure Database via the WHO GHO OData API (ghoapi.azureedge.net), indicator GHED_GGHE-DGDP_SHA2011, data dated 2025-12-05; identical to WDI SH.XPD.GHED.GD.ZS"
VOI_CONF_MACRO <- data.table::rbindlist(list(
  list("NY.GDP.MKTP.CN", "IDN", 2023L, 20892312300000000, "current IDR (GDP, local currency units)", VOI_CONF_WDI_SRC, "denominator"),
  list("NY.GDP.MKTP.CN", "UGA", 2023L, 183004386312800, "current UGX (GDP, local currency units)", VOI_CONF_WDI_SRC, "denominator"),
  list("GHED_GGHE-DGDP_SHA2011", "IDN", 2023L, 1.44882917, "% of GDP (GGHE-D: domestic general government health expenditure)", VOI_CONF_GHED_SRC, "denominator"),
  list("GHED_GGHE-DGDP_SHA2011", "UGA", 2023L, 0.92339057, "% of GDP (GGHE-D: domestic general government health expenditure)", VOI_CONF_GHED_SRC, "denominator"),
  list("FP.CPI.TOTL", "IDN", 2023L, 169.053520242749, "index, 2010 = 100", VOI_CONF_WDI_SRC, "price conversion"),
  list("FP.CPI.TOTL", "IDN", 2025L, 176.046565851866, "index, 2010 = 100", VOI_CONF_WDI_SRC, "price conversion (must equal 08's value)"),
  list("FP.CPI.TOTL", "UGA", 2023L, 209.931348222213, "index, 2010 = 100", VOI_CONF_WDI_SRC, "price conversion"),
  list("FP.CPI.TOTL", "UGA", 2025L, 224.675273352364, "index, 2010 = 100", VOI_CONF_WDI_SRC, "price conversion (must equal 08's value)"),
  list("PA.NUS.FCRF", "IDN", 2025L, 16477.8651797835, "IDR per US$, 2025 average", VOI_CONF_WDI_SRC, "price conversion (must equal 08's value)"),
  list("PA.NUS.FCRF", "UGA", 2025L, 3602.98994206813, "UGX per US$, 2025 average", VOI_CONF_WDI_SRC, "price conversion (must equal 08's value)"),
  list("PA.NUS.FCRF", "IDN", 2023L, 15236.8846620506, "IDR per US$, 2023 average", VOI_CONF_WDI_SRC, "cross-check only"),
  list("PA.NUS.FCRF", "UGA", 2023L, 3726.14045997085, "UGX per US$, 2023 average", VOI_CONF_WDI_SRC, "cross-check only"),
  list("NY.GDP.MKTP.CD", "IDN", 2023L, 1371166925749.88, "current US$ (GDP)", VOI_CONF_WDI_SRC, "cross-check only"),
  list("NY.GDP.MKTP.CD", "UGA", 2023L, 48768955863.3175, "current US$ (GDP)", VOI_CONF_WDI_SRC, "cross-check only")),
  use.names = FALSE)
data.table::setnames(VOI_CONF_MACRO, c("series", "iso3", "year", "value", "units", "source", "role"))
VOI_CONF_MACRO[, retrieved := "2026-11-10"]

VOI_CONF <- list(
  year_start = 2026L, year_end = 2050L, subdir = "conference_2026_2050",
  hurdle_voi   = c(1, 1.5, 2, 3, 4, 5),   # r for the G_r = B - r C VOI family (r = 1 must equal VSL x 1)
  hurdle_ref   = c(2, 3, 5),              # reference hurdles (descriptive)
  hurdle_prob  = c(1, 2, 3, 5),           # P(BCR > r) in the BCR summary
  hurdle_curve = sort(unique(round(c(seq(0.1, 8, by = 0.05), 1, 2, 3, 5), 8))),
  zero_tol_usd = 1,                       # |VOI| below US$1 is treated as zero for ranking only
  denominator_year = 2023L, price_year = 2025L,
  iso3 = c(Indonesia = "IDN", Uganda = "UGA"),
  macro = VOI_CONF_MACRO)

voi_conf_log <- function(od, ...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  message(msg)
  cat(msg, "\n", file = file.path(od, "conference_log.txt"), append = TRUE, sep = "")
}

# GDP and GGHE-D (2023, one year) in constant 2025 US$: LCU_2023 x CPI_2025 / CPI_2023 / XR_2025
# (own-country CPI + market exchange rate, the same local-CPI / XR basis 08 uses for unit costs).
voi_conf_denominators <- function(cc, country, er_macro) {
  iso <- cc$iso3[[country]]; y0 <- cc$denominator_year; y1 <- cc$price_year
  g <- function(s, y) {
    v <- cc$macro[series == s & iso3 == iso & year == y, value]
    if (length(v) != 1L || !is.finite(v)) stop("[", country, "] conference denominators: missing macro observation ", s, " ", iso, " ", y,
                                                " in VOI_CONF_MACRO (scripts/voi/04_analyse_and_plot.R).", call. = FALSE)
    v
  }
  for (s in c("FP.CPI.TOTL", "PA.NUS.FCRF")) {           # same 2025 basis as production 08
    v08 <- er_macro[series == s & iso3 == iso & year == y1, value]
    if (length(v08) != 1L || abs(v08 - g(s, y1)) > 1e-9 * abs(v08))
      stop("[", country, "] conference denominators: ", s, " ", y1, " differs from the value used by production 08.", call. = FALSE)
  }
  cpi_ratio <- g("FP.CPI.TOTL", y1) / g("FP.CPI.TOTL", y0); xr <- g("PA.NUS.FCRF", y1)
  gdp_lcu <- g("NY.GDP.MKTP.CN", y0); share <- g("GHED_GGHE-DGDP_SHA2011", y0) / 100
  xchk <- gdp_lcu / g("PA.NUS.FCRF", y0) / g("NY.GDP.MKTP.CD", y0) - 1   # LCU / XR_2023 vs WDI current-US$ GDP
  if (abs(xchk) > 0.01) stop("[", country, "] conference denominators: GDP LCU / XR does not match WDI current-US$ GDP (", signif(xchk, 3), ").", call. = FALSE)
  conv <- sprintf("2023 %s x CPI_2025/CPI_2023 (%.6f) / XR_2025 (%.4f %s per US$)", if (iso == "IDN") "IDR" else "UGX", cpi_ratio, xr, if (iso == "IDN") "IDR" else "UGX")
  data.table(country = country, iso3 = iso,
             denominator = c("gdp_2023", "gghed_2023"),
             value_usd2025 = c(gdp_lcu * cpi_ratio / xr, share * gdp_lcu * cpi_ratio / xr),
             reference_year = y0, price_basis = "constant 2025 US$ (market exchange rate; not PPP)",
             label = c("one year of GDP (2023)", "one year of domestic general government health expenditure (GGHE-D, 2023)"),
             definition = c("Gross domestic product, 2023, current local currency (WDI NY.GDP.MKTP.CN), converted to constant 2025 US$",
                            paste("Domestic general government health expenditure (GGHE-D) = WHO GHED GGHE-D % of GDP (2023) x WDI GDP in current LCU (2023),",
                                  "converted to constant 2025 US$. Excludes external (donor) funds, private and out-of-pocket spending; not total health",
                                  "expenditure and not a ministry-of-health budget")),
             indicator = c("NY.GDP.MKTP.CN", "GHED_GGHE-DGDP_SHA2011 x NY.GDP.MKTP.CN"),
             source = c(VOI_CONF_WDI_SRC, paste(VOI_CONF_GHED_SRC, "|", VOI_CONF_WDI_SRC)),
             original_units = c(sprintf("%s, current LCU, 2023: %s", "GDP", format(gdp_lcu, big.mark = ",", scientific = FALSE)),
                                sprintf("GGHE-D %% of GDP, 2023: %.8f", share * 100)),
             price_conversion = paste(conv, "(local-CPI / exchange-rate method, as production 08 for unit costs)"),
             crosscheck_gdp_lcu_over_xr2023_vs_wdi_usd = xchk,
             note = "Contextual scale comparison: 2026-2050 VOI / ONE year of the aggregate (not a share of cumulative 25-year resources)")
}

voi_conf_thresholds <- function(cfg, gni, cc) {
  th <- voi_thresholds(cfg, gni)[, tid := NULL]
  near <- function(a, b) vapply(a, function(z) any(abs(z - b) < 1e-10), logical(1))
  h <- data.table(family = "BCR hurdle", x = cc$hurdle_voi, gni_mult = NA_real_, is_ref = near(cc$hurdle_voi, cc$hurdle_ref),
                  unit = "BCR hurdle r: objective G_r = B - r C = r [(1/r) B - C]; re-weights costs (a different objective, not a welfare gain)")
  h[, label := paste0("BCR hurdle r = ", format(x, drop0trailing = TRUE, trim = TRUE))]
  th <- rbind(th, h, use.names = TRUE)
  th[, tid := .I]
  th[]
}

voi_conf_bootstrap <- function(cfg, dat, th, kfs) {
  N <- nrow(dat); NB <- cfg$n_boot
  set.seed(cfg$boot_seed)
  idx <- matrix(sample.int(N, N * NB, replace = TRUE), N, NB)
  blocks <- split(seq_len(NB), cut(seq_len(NB), max(1L, min(cfg$max_cores, parallel::detectCores() - 1L)), labels = FALSE))
  cl <- parallel::makeCluster(length(blocks))
  barr <- tryCatch({
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({ library(data.table); library(mgcv) }); setDTthreads(1L); NULL })
    parallel::clusterExport(cl, c("VOI_PARS", "VOI_PAR_TYPE", "voi_nb", "voi_evpi_cols", "voi_fit_gam", "voi_xform", "voi_cond",
                                  "voi_evppi_from_fit", "voi_all_estimates", "voi_boot_block"), envir = globalenv())
    dsmall <- dat[, c("draw_id", VOI_PARS, "B", "D", "C"), with = FALSE]
    parts <- parallel::parLapply(cl, blocks, voi_boot_block, idx = idx, dat = dsmall, th = th, kfs = kfs)
    a <- array(NA_real_, c(NB, length(VOI_PARS) + 1L, nrow(th)))
    for (i in seq_along(blocks)) a[blocks[[i]], , ] <- parts[[i]]
    a
  }, finally = parallel::stopCluster(cl))
  dimnames(barr) <- list(NULL, c("joint_SAP", VOI_PARS), NULL)
  barr
}

voi_conf_country <- function(cfg, country, cc = VOI_CONF) {
  p <- voi_paths(cfg, country)
  od <- file.path(p$base, cc$subdir)
  voi_assert_out(file.path(od, "x"), cfg)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  yrs <- cc$year_start:cc$year_end
  hz <- sprintf("%d-%d", cc$year_start, cc$year_end)
  inp <- voi_load_inputs(cfg, country)
  if (!all(yrs %in% inp$years)) stop("[", country, "] conference horizon ", hz, " is not inside the saved horizon.", call. = FALSE)
  if (!identical(as.integer(inp$econ$price_year), as.integer(cc$price_year)) || !identical(inp$econ$currency, "USD"))
    stop("[", country, "] saved economics are not in constant ", cc$price_year, " US$.", call. = FALSE)

  # ---------- saved PSA outputs; requested / successful / excluded; paired completeness over 2026-2050 ----------
  f_res <- file.path(p$psa, "psa_draw_results.rds"); f_ann <- file.path(p$psa, "annual_summary.rds")
  f_draws <- file.path(p$draws, "sap_draws.csv"); f_pt <- file.path(p$validation, "point_run_annual.csv")
  for (f in c(f_res, f_ann, f_draws, f_pt))
    if (!file.exists(f)) stop("[", country, "] missing saved VOI output: ", f, "\n  Produce it with scripts/voi/00_run_voi.R (stages sample, validate, pilot, psa) for run '",
                              cfg$run_id, "'; the conference path never re-runs the PSA.", call. = FALSE)
  res <- readRDS(f_res); ann_all <- readRDS(f_ann); draws <- fread(f_draws)
  ids_req <- seq_len(cfg$n_draws)
  miss_res <- setdiff(ids_req, res$draw_id)
  bad <- res[draw_id %in% ids_req & status != "ok", draw_id]
  ann <- ann_all[draw_id %in% ids_req & year %in% yrs]
  cnt <- ann[, .(n_years = uniqueN(year), n_rows = .N), by = .(draw_id, scenario)]
  complete_ids <- cnt[n_years == length(yrs) & n_rows == length(yrs) & scenario %in% c("ref", "sap"), .N, by = draw_id][N == 2L, draw_id]
  excluded <- sort(unique(c(miss_res, bad, setdiff(ids_req, complete_ids))))
  if (length(excluded))
    stop("[", country, "] ", length(excluded), " of ", length(ids_req), " requested draws failed or lack paired ", hz, " annual results (draw ids ",
         paste(utils::head(excluded, 20), collapse = ", "), if (length(excluded) > 20) ", ..." else "", ").\n  A complete-case VOI would be biased and is not ",
         "computed. Re-run stage 'psa' of scripts/voi/00_run_voi.R for run '", cfg$run_id, "' (checkpoints are reused), then re-run this path.", call. = FALSE)
  if (!setequal(draws[draw_id %in% ids_req, draw_id], ids_req)) stop("[", country, "] sap_draws.csv lacks requested draw ids.", call. = FALSE)

  # costs: saved annual cost == volumes x saved APPLIED 2025-US$ unit costs (no repeated price conversion)
  ck <- copy(ann); saved_cost <- ck$cost; voi_add_costs(ck, inp$econ)
  dcost <- max(abs(ck$cost - saved_cost) / pmax(1, abs(saved_cost)))
  if (dcost > 1e-12) stop("[", country, "] saved annual costs do not equal volumes x applied unit costs (max rel diff ", signif(dcost, 3), ").", call. = FALSE)

  # per-draw 2026-2050 totals; the same function on the full saved horizon must reproduce the saved per-draw totals
  tot <- voi_econ_totals(ann, inp$econ)
  tot_full <- voi_econ_totals(ann_all[draw_id %in% ids_req], inp$econ)
  dfull <- merge(tot_full[, .(draw_id, C_d, B_d, DA)], res[draw_id %in% ids_req, .(draw_id, C_d, B_d, DA)], by = "draw_id")
  rfull <- dfull[, max(abs(C_d.x - C_d.y) / abs(C_d.y), abs(B_d.x - B_d.y) / pmax(1, abs(B_d.y)), abs(DA.x - DA.y) / pmax(1, abs(DA.y)))]
  if (nrow(dfull) != length(ids_req) || rfull > 1e-9)
    stop("[", country, "] recomputed full-horizon totals do not reproduce psa_draw_results.rds (max rel diff ", signif(rfull, 3), ").", call. = FALSE)

  # join the draw parameters by draw id (checked against the copy stored with the PSA results)
  par_cols <- c(VOI_PARS, "regime")
  dat <- merge(draws[draw_id %in% ids_req, c("draw_id", par_cols), with = FALSE], tot, by = "draw_id")
  pchk <- merge(dat[, c("draw_id", VOI_PARS), with = FALSE], res[, c("draw_id", VOI_PARS), with = FALSE], by = "draw_id")
  if (nrow(pchk) != length(ids_req) || any(vapply(VOI_PARS, function(v) max(abs(pchk[[paste0(v, ".x")]] - pchk[[paste0(v, ".y")]])), 0) > 1e-12))
    stop("[", country, "] draw parameters in sap_draws.csv differ from those stored with the PSA results.", call. = FALSE)
  setorder(dat, draw_id)
  dat[, `:=`(B = B_d, D = DALY_d, C = C_d, G = cost_sap_d)]
  N <- nrow(dat); NB <- cfg$n_boot
  if (any(!is.finite(as.matrix(dat[, .(B, D, C, G)])))) stop("[", country, "] non-finite 2026-2050 totals.", call. = FALSE)

  # deterministic point run (production point SAP values), same horizon; cross-checked against production 08
  pa <- fread(f_pt)[year %in% yrs]
  if (pa[, .(n = uniqueN(year)), by = scenario][, !(all(n == length(yrs)) && .N == 2L)]) stop("[", country, "] point run lacks paired ", hz, " years.", call. = FALSE)
  pt <- voi_econ_totals(pa, inp$econ)
  e08 <- inp$econ$incremental[year %in% yrs]
  r08 <- c(C = sum(e08$df * e08$inc_cost), B = sum(e08$df * e08$benefit_vsl), DA = sum(e08$deaths_avert), G = sum(e08$df * e08$cost_sap))
  d08 <- max(abs(c(pt$C_d, pt$B_d, pt$DA, pt$cost_sap_d) - r08) / pmax(1, abs(r08)))
  if (d08 > 1e-9) stop("[", country, "] point-run ", hz, " totals differ from production 08 annual results (max rel diff ", signif(d08, 3), ").", call. = FALSE)

  gni <- inp$econ$gni_pc_ref_usd2025
  th <- voi_conf_thresholds(cfg, gni, cc)
  ref <- th[is_ref == TRUE]
  voi_conf_log(od, sprintf("conference %s: %d requested, %d successful, %d excluded paired draws | %d thresholds | bootstrap %d | point BCR %.3f",
                           hz, length(ids_req), N, length(excluded), nrow(th), NB, pt$B_d / pt$C_d))

  # ---------- decision quantities and EVPI ----------
  nb <- voi_nb(dat$B, dat$D, dat$C, th)
  nb_pt <- as.numeric(voi_nb(pt$B_d, pt$DALY_d, pt$C_d, th))
  einb <- colMeans(nb); sdinb <- apply(nb, 2, stats::sd); ppos <- colMeans(nb > 0)
  evpi <- voi_evpi_cols(nb)
  t_v1 <- th[family == "VSL multiple" & abs(x - 1) < 1e-10, tid]; t_h1 <- th[family == "BCR hurdle" & abs(x - 1) < 1e-10, tid]
  if (max(abs(nb[, t_v1] - nb[, t_h1])) > 1e-6 * max(1, max(abs(nb[, t_v1])))) stop("G_1 differs from INB at VSL x 1.", call. = FALSE)

  # ---------- single-parameter fits, estimates, bootstrap ----------
  fits <- list(); kfs <- list(); gdiag <- list()
  for (pp in VOI_PARS) { cf <- voi_cond(dat, pp, VOI_PAR_TYPE[[pp]], diag = TRUE); fits[[pp]] <- cf$fit; kfs[[pp]] <- cf$k; gdiag[[pp]] <- cf$diag }
  est <- voi_all_estimates(dat, th, kfs)
  if (max(abs(est["joint_SAP", ] - evpi)) != 0) stop("joint estimate differs from the direct EVPI.", call. = FALSE)
  t0 <- proc.time()[["elapsed"]]
  barr <- voi_conf_bootstrap(cfg, dat, th, kfs)
  n_boot_fail <- sum(is.na(barr[, VOI_PARS, 1]))
  voi_conf_log(od, sprintf("bootstrap: %d replicates x %d estimators x %d thresholds in %.1f s (%d failed fits)", NB, length(VOI_PARS) + 1L, nrow(th),
                           proc.time()[["elapsed"]] - t0, n_boot_fail))
  half <- dat$draw_id <= stats::median(dat$draw_id)
  est_h1 <- voi_all_estimates(dat[half], th, kfs); est_h2 <- voi_all_estimates(dat[!half], th, kfs)

  # ---------- tidy VOI table (all thresholds) ----------
  bq <- function(a, pr) apply(a, 2, stats::quantile, probs = pr, na.rm = TRUE, names = FALSE)
  curves <- rbindlist(lapply(c("joint_SAP", VOI_PARS), function(pp) {
    b <- barr[, pp, ]
    data.table(tid = th$tid, parameter = pp, estimate = est[pp, ], boot_mean = colMeans(b, na.rm = TRUE),
               boot_se = apply(b, 2, stats::sd, na.rm = TRUE), ci_lo = bq(b, 0.025), ci_hi = bq(b, 0.975),
               split_half_1 = est_h1[pp, ], split_half_2 = est_h2[pp, ],
               diff_joint_ci_lo = if (pp == "joint_SAP") NA_real_ else bq(b - barr[, "joint_SAP", ], 0.025))
  }))
  curves <- merge(th, curves, by = "tid")
  curves[, `:=`(evpi = evpi[tid], ratio_to_joint = fifelse(evpi[tid] > 0, estimate / evpi[tid], NA_real_),
                rel_se = fifelse(estimate > 0, boot_se / estimate, NA_real_), n_draws = N, n_boot = NB, horizon = hz,
                interval_type = "bootstrap 95% interval of the VOI ESTIMATOR (Monte Carlo precision; complete paired draw rows resampled), not a PSA outcome interval",
                method = fcase(parameter == "joint_SAP", "direct: EVPI = mean(max(0, NB)) - max(0, mean(NB)) (joint SAP EVPPI)",
                               parameter %in% c("M_AB", "M_BC"), "single-input GAM (thin-plate, REML) on log(input); NB components B, D, C fitted separately",
                               parameter == "partial_adherence_relative_protection", "single-input GAM (thin-plate, REML); NB components fitted separately",
                               default = "stratified: category means x sample category frequencies"))]
  curves[, parameter_label := VOI_PAR_LABEL[parameter]]
  curves[, exceeds_joint := parameter != "joint_SAP" & estimate > evpi * (1 + 1e-9) + 1e-6]
  curves[, bound_check := fcase(parameter == "joint_SAP", "reference (joint)",
                                estimate < 0, "NEGATIVE: investigate",
                                !exceeds_joint, "0 <= EVPPI <= joint: OK",
                                exceeds_joint & diff_joint_ci_lo <= 0, "exceeds joint, within bootstrap error",
                                default = "BOUND VIOLATION: exceeds joint beyond bootstrap error (investigate)")]
  curves[, status := fcase(evpi == 0, "EVPI = 0: no decision uncertainty in the sample (every draw favours the same option)",
                           ci_hi < 0.05 * evpi, "negligible: bootstrap 97.5% bound < 5% of the joint value",
                           parameter != "joint_SAP" & estimate == 0,
                           "EXPLORATORY: point estimate 0 but bootstrap replicates are positive; more draws needed",
                           rel_se <= 0.20, "stable at this sample size (bootstrap relative SE <= 0.20)",
                           default = "EXPLORATORY: unstable at this sample size (bootstrap relative SE > 0.20); more draws needed")]
  curves[, value_type := fifelse(family == "BCR hurdle", "VOI for the cost-reweighted objective G_r = B - r C (not a welfare gain)",
                                 "expected value of perfect information for the ref vs scale-up decision (constant 2025 US$, discounted to 2026)")]
  setorder(curves, family, x, parameter)
  if (curves[grepl("^NEGATIVE|^BOUND", bound_check), .N] > 0)
    stop("[", country, "] VOI bound check failed: ", paste(curves[grepl("^NEGATIVE|^BOUND", bound_check), unique(paste(label, parameter))], collapse = "; "), call. = FALSE)

  # ---------- reference thresholds: tie-aware ranks + bootstrap rank probabilities ----------
  tol <- cc$zero_tol_usd
  cmpt <- curves[is_ref == TRUE]
  cmpt[parameter != "joint_SAP", rank_among_individual := {
    r <- rep(NA_integer_, .N); pos <- estimate > tol
    if (any(pos)) r[pos] <- as.integer(frank(-estimate[pos], ties.method = "min"))
    r
  }, by = tid]
  cmpt[parameter != "joint_SAP", rank_note := fcase(evpi <= tol, "not ranked: EVPI = 0 (no decision uncertainty)",
                                                    estimate <= tol, "not ranked: estimate 0",
                                                    default = "")]
  cmpt[parameter != "joint_SAP" & !is.na(rank_among_individual),
       rank_note := fifelse(duplicated(rank_among_individual) | duplicated(rank_among_individual, fromLast = TRUE), "tied", rank_note), by = tid]
  rk <- rbindlist(lapply(ref$tid, function(t) {
    b <- barr[, VOI_PARS, t, drop = FALSE][, , 1]
    mx <- apply(b, 1, max, na.rm = TRUE)
    top <- (b >= mx - 1e-9 * pmax(1, abs(mx))) & (mx > tol)
    ntop <- rowSums(top)
    data.table(tid = t, parameter = VOI_PARS, p_rank_first_unique = colMeans(top & ntop == 1L, na.rm = TRUE),
               p_tied_first = colMeans(top & ntop > 1L, na.rm = TRUE), p_all_zero = mean(mx <= tol))
  }))
  cmpt <- merge(cmpt, rk, by = c("tid", "parameter"), all.x = TRUE)
  cmpt[, country := country]
  setcolorder(cmpt, c("country", "horizon", "family", "label", "x", "gni_mult", "unit", "parameter", "parameter_label", "estimate", "ci_lo", "ci_hi",
                      "boot_se", "rel_se", "evpi", "ratio_to_joint", "rank_among_individual", "rank_note", "p_rank_first_unique", "p_tied_first",
                      "p_all_zero", "bound_check", "status", "split_half_1", "split_half_2", "method", "interval_type", "value_type", "n_draws", "n_boot"))
  setorder(cmpt, family, x, rank_among_individual, na.last = TRUE)

  # ---------- decision summaries (all thresholds) ----------
  dec <- data.table(tid = th$tid, E_NB = einb, mc_se_E_NB = sdinb / sqrt(N), sd_NB = sdinb, P_NB_pos = ppos,
                    mc_se_P = sqrt(ppos * (1 - ppos) / N), EVPI = evpi, NB_point_values = nb_pt,
                    optimal_current_info = ifelse(einb > 0, "sap (scale-up)", "ref (reference care)"),
                    point_value_decision = ifelse(nb_pt > 0, "sap (scale-up)", "ref (reference care)"))
  dec <- merge(th, dec, by = "tid")
  dec <- merge(dec, curves[parameter == "joint_SAP", .(tid, EVPI_ci_lo = ci_lo, EVPI_ci_hi = ci_hi, EVPI_boot_se = boot_se)], by = "tid")
  dec[, `:=`(E_NB_ci_lo = E_NB - 1.96 * mc_se_E_NB, E_NB_ci_hi = E_NB + 1.96 * mc_se_E_NB,
             P_draws_favour_optimal = fifelse(E_NB > 0, P_NB_pos, 1 - P_NB_pos),
             P_NB_pos_note = fifelse(P_NB_pos %in% c(0, 1), sprintf("no draw on the other side; rule-of-three 95%% bound ~ %.4f", 3 / N), ""),
             country = country, horizon = hz, n_draws = N)]
  setcolorder(dec, c("country", "horizon", "family", "label", "x", "is_ref", "E_NB", "E_NB_ci_lo", "E_NB_ci_hi", "mc_se_E_NB", "P_NB_pos", "mc_se_P",
                     "P_draws_favour_optimal", "P_NB_pos_note", "optimal_current_info", "point_value_decision", "EVPI", "EVPI_ci_lo", "EVPI_ci_hi"))

  # ---------- BCR summary (incremental, discounted) ----------
  cpos <- dat$C > 0
  bcr <- dat[cpos, B / C]
  q3 <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
  qb <- if (length(bcr)) q3(bcr) else rep(NA_real_, 3)
  pr <- vapply(cc$hurdle_prob, function(r) mean(dat$B - r * dat$C > 0), 0)
  bsum <- data.table(country = country, horizon = hz,
                     definition = "BCR = discounted incremental monetised mortality benefit (deaths averted x reference-case VSL_t) / discounted incremental programme cost, screening + SAP scale-up vs reference care; 3% to 2026; constant 2025 US$; benefits and costs accrued 2026-2050 only",
                     n_requested = length(ids_req), n_successful = N, n_excluded = length(excluded),
                     n_C_positive = sum(cpos), n_C_nonpositive = sum(!cpos), n_C_nonpositive_B_positive = sum(!cpos & dat$B > 0),
                     bcr_median = qb[2], bcr_q025 = qb[1], bcr_q975 = qb[3],
                     bcr_mean_of_draw_ratios = if (length(bcr)) mean(bcr) else NA_real_,
                     E_B = mean(dat$B), E_C = mean(dat$C), ratio_of_expectations = if (mean(dat$C) > 0) mean(dat$B) / mean(dat$C) else NA_real_,
                     E_INB = mean(dat$B - dat$C), mc_se_E_INB = stats::sd(dat$B - dat$C) / sqrt(N),
                     point_B = pt$B_d, point_C = pt$C_d, point_bcr = if (pt$C_d > 0) pt$B_d / pt$C_d else NA_real_, point_DA = pt$DA,
                     decision_current_info = if (mean(dat$B - dat$C) > 0) "sap (scale-up)" else "ref (reference care)",
                     note = paste("E[B]/E[C] (ratio of expectations) is the decision-relevant ratio; the mean of draw-level ratios E[B/C] is descriptive only.",
                                  "With C > 0, BCR > r <=> G_r = B - r C > 0; probabilities are computed as P(G_r > 0) over ALL draws",
                                  if (all(cpos)) "(every draw has C > 0)." else "(draws with C <= 0 handled through net benefit, counted above)."))
  for (i in seq_along(cc$hurdle_prob)) {
    r <- cc$hurdle_prob[i]; nm <- paste0("p_bcr_gt_", format(r, drop0trailing = TRUE))
    bsum[, (nm) := pr[i]]; bsum[, (paste0("mc_se_", nm)) := sqrt(pr[i] * (1 - pr[i]) / N)]
  }
  hcurve <- data.table(country = country, horizon = hz, r = cc$hurdle_curve,
                       p_G_r_pos = vapply(cc$hurdle_curve, function(r) mean(dat$B - r * dat$C > 0), 0))
  hcurve[, mc_se := sqrt(p_G_r_pos * (1 - p_G_r_pos) / N)]

  outc <- rbindlist(lapply(list(c("C_d", "incremental programme cost, discounted (US$)"), c("B_d", "incremental monetised mortality benefit (VSL x 1), discounted (US$)"),
                                c("INB_vsl_d", "incremental net benefit at VSL x 1, discounted (US$)"), c("cost_sap_d", "gross scale-up programme cost, discounted (US$)"),
                                c("cost_ref_d", "reference-care programme cost, discounted (US$)"),
                                c("DA", "RHD deaths averted, 2026-2050 (undiscounted)"), c("DALY_d", "DALYs averted PROXY (discounted deaths averted x 30; not modelled DALYs)"),
                                c("bcr_d", "benefit-cost ratio, discounted (draws with C > 0)"), c("cost_per_daly", "cost per DALY averted proxy (US$)")), function(v) {
    x <- dat[[v[1]]]; x <- x[is.finite(x)]; qq <- q3(x)
    data.table(country = country, horizon = hz, outcome = v[1], description = v[2], n = length(x), mean = mean(x), mc_se = stats::sd(x) / sqrt(length(x)),
               sd = stats::sd(x), q025 = qq[1], median = qq[2], q975 = qq[3], point_value = pt[[v[1]]],
               interval_type = "PSA outcome distribution (2.5th-97.5th percentiles over draws)")
  }))

  # ---------- normalisation: gross programme cost, GDP, GGHE-D (and incremental cost, separately) ----------
  er_macro <- as.data.table(readRDS(inp$files[["economic_results"]])$macro_observations)
  den_mac <- voi_conf_denominators(cc, country, er_macro)
  den <- rbind(data.table(country = country, iso3 = cc$iso3[[country]], denominator = c("programme_cost_gross", "programme_cost_incremental"),
                          value_usd2025 = c(mean(dat$G), mean(dat$C)), reference_year = NA_integer_,
                          price_basis = "constant 2025 US$, discounted at 3% to 2026 (same basis as the VOI numerator)",
                          label = c("expected gross scale-up programme cost, 2026-2050", "expected incremental programme cost vs reference care, 2026-2050"),
                          definition = c(sprintf("E[sum_t df_t cost_sap,t], t = %s, over the %d valid draws: screening encounters, on-BPG person-years and surgeries in the scale-up arm x saved applied 2025-US$ unit costs (GROSS programme cost; confirmatory echoes uncosted)", hz, N),
                                         sprintf("E[sum_t df_t (cost_sap,t - cost_ref,t)], t = %s, over the %d valid draws (INCREMENTAL cost; backup comparison only)", hz, N)),
                          indicator = c("cost_sap_d", "C_d"), source = "this analysis (saved annual arm-level PSA results)",
                          original_units = "constant 2025 US$", price_conversion = "none (unit costs already in constant 2025 US$ from production 08)",
                          crosscheck_gdp_lcu_over_xr2023_vs_wdi_usd = NA_real_, note = c("Denominator for VOI / programme cost", "Shown only in backup; requires E[C] > 0")),
               den_mac, use.names = TRUE)
  if (den[denominator == "programme_cost_gross", value_usd2025] <= 0) stop("[", country, "] non-positive gross programme cost.", call. = FALSE)
  nrm_src <- cmpt[family != "BCR hurdle"]
  nrm <- rbindlist(lapply(seq_len(nrow(den)), function(i) {
    dd <- den[i]
    if (dd$denominator == "programme_cost_incremental" && !(dd$value_usd2025 > 0)) return(NULL)
    nrm_src[, .(country, horizon, family, label, x, is_ref, parameter, parameter_label, estimate_usd = estimate, ci_lo_usd = ci_lo, ci_hi_usd = ci_hi,
                denominator = dd$denominator, denominator_label = dd$label, denominator_usd2025 = dd$value_usd2025,
                multiple = estimate / dd$value_usd2025, multiple_ci_lo = ci_lo / dd$value_usd2025, multiple_ci_hi = ci_hi / dd$value_usd2025)]
  }))
  nrm[, `:=`(percent = 100 * multiple, percent_ci_lo = 100 * multiple_ci_lo, percent_ci_hi = 100 * multiple_ci_hi,
             interpretation = "scale comparison only: not fiscal savings, an available research budget or a research benefit-cost ratio")]
  nrm[denominator %in% c("gdp_2023", "gghed_2023"), interpretation := paste("2026-2050 VOI / ONE year of the aggregate;", interpretation)]
  nrm[denominator == "programme_cost_incremental", interpretation := paste("VOI / expected INCREMENTAL cost (distinct from gross programme cost);", interpretation)]

  # ---------- convergence, binned cross-check, age categories, prior-scenario re-weighting ----------
  conv_n <- sort(unique(c(seq(100L, N, by = 100L), N)))
  conv <- rbindlist(lapply(conv_n, function(n) {
    d <- dat[seq_len(n)]; e <- voi_all_estimates(d, ref, kfs); m <- voi_nb(d$B, d$D, d$C, ref)
    rbind(data.table(n = n, tid = ref$tid, parameter = rep(rownames(e), each = ncol(e)), value = as.numeric(t(e)), quantity = "EVPPI / EVPI"),
          data.table(n = n, tid = ref$tid, parameter = "NB", value = colMeans(m), quantity = "E[NB]"),
          data.table(n = n, tid = ref$tid, parameter = "NB", value = colMeans(m > 0), quantity = "P(NB > 0)"))
  }))
  conv <- merge(th[, .(tid, family, label)], conv, by = "tid")
  setorder(conv, tid, quantity, parameter, n)
  conv[, horizon := hz]
  binned <- rbindlist(lapply(c("M_AB", "M_BC", "partial_adherence_relative_protection"), function(pp)
    data.table(parameter = pp, tid = th$tid, evppi_binned10 = voi_binned(dat, pp, th), evppi_gam = est[pp, ])))[tid %in% ref$tid]
  binned <- merge(th[, .(tid, family, label)], binned, by = "tid")
  agecat <- rbindlist(lapply(ref$tid, function(t) {
    m <- nb[, t]
    out <- dat[, .(idx = .I), by = sap_protection_age_max][, .(n = .N, mean_NB = mean(m[idx]), se_mean_NB = stats::sd(m[idx]) / sqrt(.N)), by = sap_protection_age_max]
    out[, `:=`(share = n / N, tid = t)]
  }))
  agecat <- merge(th[, .(tid, family, label)], agecat, by = "tid")
  w0 <- cfg$prior$w_trial; z <- dat$regime == "trial_efficacy"
  rew <- rbindlist(lapply(cfg$prior_reweight_w, function(w1) {
    w <- ifelse(z, w1 / w0, (1 - w1) / (1 - w0)); wn <- w / sum(w); m <- nb[, ref$tid, drop = FALSE]
    data.table(tid = ref$tid, w_trial = w1, ess = sum(w)^2 / sum(w^2), E_NB = colSums(wn * m), P_NB_pos = colSums(wn * (m > 0)), EVPI = voi_evpi_cols(m, w))
  }))
  rew <- merge(th[, .(tid, family, label)], rew, by = "tid")[, note := "importance re-weighting of the same draws; joint quantities only; exploratory"]
  mcl <- data.table(country = country, horizon = hz, n_requested = length(ids_req), n_successful = N, n_excluded = length(excluded), n_boot = NB,
                    n_boot_failed_fits = n_boot_fail, min_category_n = min(agecat$n), regime_trial_n = sum(z), regime_recurrence_n = sum(!z),
                    n_reference_individual = cmpt[parameter != "joint_SAP", .N], n_exploratory = cmpt[grepl("^EXPLORATORY", status), .N],
                    n_stable = cmpt[grepl("^stable", status), .N], n_negligible = cmpt[grepl("^negligible", status), .N],
                    n_evpi_zero = cmpt[grepl("^EVPI = 0", status) & parameter != "joint_SAP", .N], n_exceeds_joint = cmpt[exceeds_joint == TRUE, .N],
                    max_rel_change_E_NB_400_to_500 = conv[quantity == "E[NB]" & n %in% c(max(conv_n[conv_n < N]), N), abs(diff(value)) / abs(value[2]), by = tid][, max(V1)],
                    max_abs_cost_check = dcost, full_horizon_reproduction_max_rel_diff = rfull, point_vs_08_max_rel_diff = d08,
                    note = "N = 500 initial sample (not extended). EXPLORATORY estimates need more draws (raise N_DRAWS under the same run id).")

  # ---------- write (all under <country>/conference_2026_2050/) ----------
  tw <- function(x, f) voi_fwrite(x, file.path(od, f), cfg)
  tw(bsum, "bcr_summary.csv"); tw(hcurve, "bcr_hurdle_curve.csv"); tw(outc, "psa_outcome_summary.csv")
  tw(dec, "decision_curves.csv"); tw(dec[is_ref == TRUE], "decision_summary_reference.csv")
  tw(cmpt, "evppi_reference.csv"); tw(cbind(country = country, curves), "evppi_curves.csv")
  tw(den, "denominators.csv"); tw(nrm, "voi_normalised.csv")
  tw(cbind(country = country, horizon = hz, rbindlist(gdiag)), "gam_diagnostics.csv"); tw(cbind(country = country, horizon = hz, binned), "evppi_binned_check.csv")
  tw(cbind(country = country, conv), "convergence.csv"); tw(cbind(country = country, horizon = hz, agecat), "age_category_diagnostics.csv")
  tw(cbind(country = country, horizon = hz, rew), "prior_scenario_reweighting.csv"); tw(mcl, "mc_limitations.csv")
  tw(dat[, c("draw_id", par_cols, "C_d", "B_d", "INB_vsl_d", "DA", "DA_d", "DALY_d", "cost_sap_d", "cost_ref_d", "bcr_d"), with = FALSE], "draw_totals_2026_2050.csv")
  voi_save(list(estimates = est, bootstrap = barr, thresholds = th, kfs = kfs, boot_seed = cfg$boot_seed, horizon = hz),
           file.path(od, "evppi_bootstrap_array.rds"), cfg)
  man <- readRDS(file.path(p$base, "manifest.rds"))
  srcs <- c(f_ann, f_res, f_draws, f_pt, inp$files)
  prov <- list(run_id = cfg$run_id, country = country, horizon = c(cc$year_start, cc$year_end), created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
               r_version = R.version.string, n_requested = length(ids_req), n_successful = N, n_excluded = length(excluded), excluded_ids = excluded,
               prior_hash = man$prior_hash, prior_status = man$prior_status, engine_version = man$engine_version, engine_spec = man$engine_spec,
               calibration_adequacy = inp$calibration_adequacy, currency = inp$econ$currency, price_year = inp$econ$price_year,
               disc_rate = inp$econ$disc_rate, discount_base_year = inp$econ$base_year, unit_costs_usd2025 = inp$econ$uc,
               vsl_source = "saved annual VSL schedule from rhd_economic_results.rds ($incremental$vsl; reference case, Robinson et al. 2019)",
               dalys_per_death = inp$econ$dalys_per_death, gni_pc_ref_usd2025 = gni, n_boot = NB, boot_seed = cfg$boot_seed, gam_basis = kfs,
               source_md5 = voi_md5(srcs), script_md5 = voi_md5(file.path(cfg$voi_dir, "04_analyse_and_plot.R")),
               denominators = den, macro_observations = cc$macro[iso3 == cc$iso3[[country]]],
               checks = list(cost_check_max_rel_diff = dcost, full_horizon_reproduction_max_rel_diff = rfull, point_vs_08_max_rel_diff = d08))
  voi_save(prov, file.path(od, "provenance.rds"), cfg)
  for (i in seq_len(nrow(dec[is_ref == TRUE]))) with(dec[is_ref == TRUE][i],
    voi_conf_log(od, sprintf("%-32s E[NB] %9.2f m | P(NB>0) %.3f | EVPI %7.3f m (95%% %.3f-%.3f)", label, E_NB / 1e6, P_NB_pos, EVPI / 1e6, EVPI_ci_lo / 1e6, EVPI_ci_hi / 1e6)))
  voi_conf_log(od, sprintf("BCR %s: median %.2f (95%% %.2f-%.2f) | E[B]/E[C] %.2f | P(BCR>1) %.3f | point %.2f | gross programme cost %.1f m",
                           hz, bsum$bcr_median, bsum$bcr_q025, bsum$bcr_q975, bsum$ratio_of_expectations, bsum$p_bcr_gt_1, bsum$point_bcr, mean(dat$G) / 1e6))
  list(bcr = bsum, hurdle_curve = hcurve, outcomes = outc, decision = dec, evppi_reference = cmpt, curves = cbind(country = country, curves),
       denominators = den, normalised = nrm, convergence = cbind(country = country, conv), gam = cbind(country = country, rbindlist(gdiag)),
       binned = cbind(country = country, binned), agecat = cbind(country = country, agecat), reweighting = cbind(country = country, rew),
       mc_limitations = mcl, draw_totals = cbind(country = country, dat[, c("draw_id", par_cols, "C_d", "B_d", "DA", "DALY_d", "cost_sap_d"), with = FALSE]),
       provenance = prov)
}

voi_run_conference <- function(cfg, cc = VOI_CONF) {
  res <- lapply(cfg$countries, function(ct) voi_conf_country(cfg, ct, cc)); names(res) <- cfg$countries
  od <- file.path(cfg$out_root, cfg$run_id, cc$subdir)
  voi_assert_out(file.path(od, "x"), cfg)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  pick <- function(nm) rbindlist(lapply(res, `[[`, nm), use.names = TRUE, fill = TRUE)
  tabs <- c("bcr", "hurdle_curve", "outcomes", "decision", "evppi_reference", "curves", "denominators", "normalised", "convergence", "gam",
            "binned", "agecat", "reweighting", "mc_limitations", "draw_totals")
  comb <- setNames(lapply(tabs, pick), tabs)
  comb$macro_observations <- cc$macro
  comb$meta <- list(run_id = cfg$run_id, horizon = c(cc$year_start, cc$year_end), countries = cfg$countries, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                    n_draws = cfg$n_draws, n_boot = cfg$n_boot, reference_vsl_mult = 1, vsl_mult_ref = cfg$vsl_mult_ref, daly_gni_ref = cfg$daly_gni_ref,
                    hurdle_ref = cc$hurdle_ref, prior = cfg$prior, point = cfg$point,
                    provenance = lapply(res, `[[`, "provenance"),
                    script = "scripts/voi/04_analyse_and_plot.R (VOI_ANALYSIS_PATH = 'conference')")
  voi_save(comb, file.path(od, "voi_conference_2026_2050.rds"), cfg)
  for (nm in c("bcr", "decision", "evppi_reference", "denominators", "normalised", "mc_limitations"))
    voi_fwrite(comb[[nm]][, lapply(.SD, function(v) if (is.list(v)) vapply(v, paste, "", collapse = "; ") else v)],
               file.path(od, paste0(nm, ".csv")), cfg)
  voi_fwrite(cc$macro, file.path(od, "macro_observations.csv"), cfg)
  message(sprintf("== conference path %d-%d written: %s ==", cc$year_start, cc$year_end, od))
  invisible(comb)
}

# Which analysis to run when this file is sourced: "full" (default: the 2026-2100 analysis), "conference"
# (2026-2050 conference path only, reading the saved PSA outputs) or "both".
VOI_ANALYSIS_PATH <- get0("VOI_ANALYSIS_PATH", ifnotfound = "full")
if (!VOI_ANALYSIS_PATH %in% c("full", "conference", "both"))
  stop("VOI_ANALYSIS_PATH must be 'full', 'conference' or 'both'.", call. = FALSE)
if (identical(get0("VOI_CURRENT_STAGE", ifnotfound = ""), "analyse") || !exists("VOI_CURRENT_STAGE")) {
  if (VOI_ANALYSIS_PATH %in% c("full", "both")) voi_run_analysis(VOI_CFG)
  if (VOI_ANALYSIS_PATH %in% c("conference", "both")) voi_run_conference(VOI_CFG)
}
