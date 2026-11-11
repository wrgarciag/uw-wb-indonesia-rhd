# ==============================================================================
# SAP VOI analysis — 01: joint sampling of the four SAP inputs + prior diagnostics
# scripts/voi/01_sample_sap_parameters.R
#
# Implements the prior declared in 00_run_voi.R (VOI_PRIOR; PROVISIONAL, not agreed):
#   shared regime indicator  Z ~ Bernoulli(w_trial)          (fixed weight in this scenario)
#   Z = 1 (trial efficacy):  (log M_AB, log M_BC) bivariate normal, rho = 0.5
#   Z = 0 (recurrence only): M_AB = M_BC = 1 - s (1 - RR_rec), RR_rec ~ LN, s ~ Beta
#   M > 1 redrawn (truncation at 1; rejections counted per draw)
#   sap_protection_age_max   ~ categorical {20, 39, 100}
#   partial_adherence_relative_protection ~ Beta
# Draw i is generated from its OWN L'Ecuyer-CMRG stream (i-th nextRNGStream of the base
# seed), so a later run with a larger N_DRAWS reproduces draws 1..N exactly and only adds
# new ids. Same draws for both countries (written to each country folder).
#
# Writes (per country): draws/sap_draws.csv (+ .rds), draws/prior_specification.rds,
#   tables/prior_diagnostics.csv, tables/prior_category_counts.csv, figures/prior_draws.png
# ==============================================================================

if (!exists("VOI_CFG")) { VOI_CONFIG_ONLY <- TRUE; source(file.path("scripts", "voi", "00_run_voi.R")); rm(VOI_CONFIG_ONLY) }
library(data.table)

# ------------------------------------------------------------------------------
# sampling functions
# ------------------------------------------------------------------------------
# one draw from the CURRENT RNG state (stream already selected by the caller)
voi_sample_one <- function(prior) {
  z_trial <- stats::runif(1) < prior$w_trial
  n_rej <- 0L
  repeat {
    if (z_trial) {
      e <- stats::rnorm(2)
      z1 <- e[1]; z2 <- prior$rho_log_trial * e[1] + sqrt(1 - prior$rho_log_trial^2) * e[2]
      M_AB <- exp(prior$M_AB_trial[["meanlog"]] + prior$M_AB_trial[["sdlog"]] * z1)
      M_BC <- exp(prior$M_BC_trial[["meanlog"]] + prior$M_BC_trial[["sdlog"]] * z2)
      RR <- NA_real_; s <- NA_real_
    } else {
      RR <- exp(prior$RR_recurrence[["meanlog"]] + prior$RR_recurrence[["sdlog"]] * stats::rnorm(1))
      s  <- stats::rbeta(1, prior$share_recurrence[["shape1"]], prior$share_recurrence[["shape2"]])
      M_AB <- M_BC <- 1 - s * (1 - RR)
    }
    if (M_AB > 0 && M_BC > 0 && M_AB <= prior$M_upper && M_BC <= prior$M_upper) break
    n_rej <- n_rej + 1L
    if (n_rej > prior$max_rejections) stop("voi_sample_one: truncation rejected more than ", prior$max_rejections, " times.", call. = FALSE)
  }
  age <- prior$protection_age$values[sample.int(length(prior$protection_age$values), 1L, prob = prior$protection_age$probs)]
  rel <- stats::rbeta(1, prior$partial_protection[["shape1"]], prior$partial_protection[["shape2"]])
  list(M_AB = M_AB, M_BC = M_BC, sap_protection_age_max = as.integer(age), partial_adherence_relative_protection = rel,
       regime = if (z_trial) "trial_efficacy" else "recurrence_only", RR_recurrence = RR, share_recurrence = s,
       n_truncation_rejections = n_rej)
}

# draws for the requested ids, each from its own stream; the caller's RNG state is restored
voi_sample_draws <- function(ids, seed, prior) {
  ids <- as.integer(ids)
  stopifnot(length(ids) > 0L, all(ids >= 1L), !anyDuplicated(ids))
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  s <- get(".Random.seed", envir = globalenv())
  streams <- vector("list", max(ids))
  for (i in seq_len(max(ids))) { s <- parallel::nextRNGStream(s); streams[[i]] <- s }
  rbindlist(lapply(ids, function(i) {
    assign(".Random.seed", streams[[i]], envir = globalenv())
    d <- voi_sample_one(prior)
    as.data.table(c(list(draw_id = i), d))
  }))
}

# large-sample vectorised prior Monte Carlo (diagnostics only; separate seed, never used in the PSA)
voi_prior_mc <- function(n, prior, seed) {
  old_kind <- RNGkind(); on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
  z <- stats::runif(n) < prior$w_trial
  M_AB <- M_BC <- numeric(n); todo <- rep(TRUE, n); n_rej <- 0L
  while (any(todo)) {
    i <- which(todo & z); j <- which(todo & !z)
    if (length(i)) {
      e1 <- stats::rnorm(length(i)); e2 <- stats::rnorm(length(i))
      M_AB[i] <- exp(prior$M_AB_trial[["meanlog"]] + prior$M_AB_trial[["sdlog"]] * e1)
      M_BC[i] <- exp(prior$M_BC_trial[["meanlog"]] + prior$M_BC_trial[["sdlog"]] * (prior$rho_log_trial * e1 + sqrt(1 - prior$rho_log_trial^2) * e2))
    }
    if (length(j)) {
      RR <- exp(prior$RR_recurrence[["meanlog"]] + prior$RR_recurrence[["sdlog"]] * stats::rnorm(length(j)))
      s  <- stats::rbeta(length(j), prior$share_recurrence[["shape1"]], prior$share_recurrence[["shape2"]])
      M_AB[j] <- M_BC[j] <- 1 - s * (1 - RR)
    }
    bad <- todo & (M_AB > prior$M_upper | M_BC > prior$M_upper | M_AB <= 0 | M_BC <= 0)
    n_rej <- n_rej + sum(bad); todo <- bad
  }
  data.table(M_AB = M_AB, M_BC = M_BC, regime = ifelse(z, "trial_efficacy", "recurrence_only"),
             sap_protection_age_max = prior$protection_age$values[sample.int(3L, n, replace = TRUE, prob = prior$protection_age$probs)],
             partial_adherence_relative_protection = stats::rbeta(n, prior$partial_protection[["shape1"]], prior$partial_protection[["shape2"]]),
             n_rej_total = n_rej)
}

voi_prior_summary <- function(d, label) {
  q <- function(x, p) as.numeric(stats::quantile(x, p, names = FALSE))
  rbindlist(lapply(c("M_AB", "M_BC", "partial_adherence_relative_protection", "sap_protection_age_max"), function(v) {
    x <- d[[v]]
    data.table(sample = label, parameter = v, n = length(x), mean = mean(x), sd = stats::sd(x),
               q025 = q(x, 0.025), q25 = q(x, 0.25), median = q(x, 0.5), q75 = q(x, 0.75), q975 = q(x, 0.975),
               min = min(x), max = max(x),
               p_ge_production_point = if (v %in% c("M_AB", "M_BC")) mean(x >= 0.45) else NA_real_)
  }))
}

# ------------------------------------------------------------------------------
# run: sample (or extend) the draws, then diagnostics, for every country
# ------------------------------------------------------------------------------
voi_run_sampling <- function(cfg) {
  prior_hash <- voi_hash_object(cfg$prior)
  ids_all <- seq_len(max(cfg$n_draws, cfg$n_pilot))
  draws <- NULL
  for (ct in cfg$countries) {
    p <- voi_paths(cfg, ct)
    f_csv <- file.path(p$draws, "sap_draws.csv"); f_spec <- file.path(p$draws, "prior_specification.rds")
    if (file.exists(f_spec)) {
      old <- readRDS(f_spec)
      if (!identical(old$prior_hash, prior_hash) || !identical(old$seed, cfg$seed))
        stop("Run '", cfg$run_id, "' [", ct, "] already holds draws from a DIFFERENT prior or seed. ",
             "Draws can only be extended under the same prior and seed; use a new VOI_RUN_ID.", call. = FALSE)
    }
    if (is.null(draws)) draws <- voi_sample_draws(ids_all, cfg$seed, cfg$prior)
    if (file.exists(f_csv)) {                                 # extension: existing ids must be reproduced exactly
      old_d <- fread(f_csv)
      chk <- merge(old_d[, .(draw_id, M_AB, M_BC, sap_protection_age_max, partial_adherence_relative_protection)],
                   draws[, .(draw_id, M_AB, M_BC, sap_protection_age_max, partial_adherence_relative_protection)], by = "draw_id")
      if (nrow(chk) != min(nrow(old_d), nrow(draws)) ||
          chk[, max(abs(M_AB.x - M_AB.y), abs(M_BC.x - M_BC.y), abs(partial_adherence_relative_protection.x - partial_adherence_relative_protection.y),
                    abs(sap_protection_age_max.x - sap_protection_age_max.y))] > 1e-12)
        stop("[", ct, "] re-sampled draws do not reproduce the saved draws (stream mismatch).", call. = FALSE)
      voi_log(p, sprintf("draws: %d saved draws reproduced exactly; total now %d", nrow(old_d), nrow(draws)))
    }
    voi_fwrite(draws, f_csv, cfg)
    voi_save(draws, file.path(p$draws, "sap_draws.rds"), cfg)
    voi_save(list(prior = cfg$prior, prior_hash = prior_hash, seed = cfg$seed, rng = "L'Ecuyer-CMRG; draw i = i-th nextRNGStream of set.seed(seed)",
                  n_draws = nrow(draws), point_values = cfg$point,
                  status = cfg$prior$status), f_spec, cfg)

    # ---- diagnostics: realised sample vs large-sample prior Monte Carlo ----
    mc <- voi_prior_mc(cfg$prior_diag_n, cfg$prior, cfg$seed + 1L)
    smp <- draws[draw_id <= cfg$n_draws]
    diag <- rbind(voi_prior_summary(smp, sprintf("PSA sample (n = %d)", nrow(smp))),
                  voi_prior_summary(mc, sprintf("prior Monte Carlo (n = %d)", nrow(mc))))
    dep <- function(d, label) data.table(
      sample = label,
      share_trial_regime = mean(d$regime == "trial_efficacy"),
      pearson_log_MAB_MBC = stats::cor(log(d$M_AB), log(d$M_BC)),
      spearman_MAB_MBC = stats::cor(d$M_AB, d$M_BC, method = "spearman"),
      pearson_log_within_trial = d[regime == "trial_efficacy", stats::cor(log(M_AB), log(M_BC))],
      share_MAB_eq_MBC = mean(abs(d$M_AB - d$M_BC) < 1e-15),
      spearman_MAB_partial = stats::cor(d$M_AB, d$partial_adherence_relative_protection, method = "spearman"),
      spearman_MAB_age = stats::cor(d$M_AB, d$sap_protection_age_max, method = "spearman"),
      min_M = min(c(d$M_AB, d$M_BC)), max_M = max(c(d$M_AB, d$M_BC)),
      truncation_rejections = if ("n_truncation_rejections" %in% names(d)) sum(d$n_truncation_rejections) else d$n_rej_total[1])
    dep_tab <- rbind(dep(smp, sprintf("PSA sample (n = %d)", nrow(smp))), dep(mc, sprintf("prior Monte Carlo (n = %d)", nrow(mc))))
    cats <- merge(smp[, .(n = .N), by = sap_protection_age_max][, share := n / sum(n)],
                  data.table(sap_protection_age_max = cfg$prior$protection_age$values, prior_prob = cfg$prior$protection_age$probs),
                  by = "sap_protection_age_max", all = TRUE)
    # bounds / validity of the PSA sample (stop on any violation)
    viol <- smp[!(M_AB > 0 & M_AB <= 1 & M_BC > 0 & M_BC <= 1 & partial_adherence_relative_protection >= 0 &
                    partial_adherence_relative_protection <= 1 & sap_protection_age_max %in% cfg$prior$protection_age$values)]
    if (nrow(viol)) stop("[", ct, "] ", nrow(viol), " sampled draws violate parameter bounds.", call. = FALSE)
    if (smp[regime == "recurrence_only", any(abs(M_AB - M_BC) > 0)]) stop("[", ct, "] recurrence regime must give M_AB == M_BC.", call. = FALSE)
    voi_fwrite(diag, file.path(p$tables, "prior_diagnostics.csv"), cfg)
    voi_fwrite(dep_tab, file.path(p$tables, "prior_dependence.csv"), cfg)
    voi_fwrite(cats, file.path(p$tables, "prior_category_counts.csv"), cfg)

    # ---- figure ----
    if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)) {
      library(ggplot2)
      col_reg <- c(trial_efficacy = "#2a78d6", recurrence_only = "#eb6834")   # reference categorical slots 1-2
      th <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3), plot.background = element_rect(fill = "#fcfcfb", colour = NA),
                                                  plot.title = element_text(size = 10, face = "bold"))
      long <- melt(smp[, .(draw_id, regime, M_AB, M_BC, partial = partial_adherence_relative_protection)],
                   id.vars = c("draw_id", "regime"))
      g1 <- ggplot(long, aes(value, fill = regime)) + geom_histogram(bins = 40, position = "stack", colour = NA) +
        geom_vline(data = data.table(variable = c("M_AB", "M_BC", "partial"), v = c(0.45, 0.45, 0.5)), aes(xintercept = v),
                   linetype = "dashed", linewidth = 0.4) +
        facet_wrap(~variable, scales = "free", labeller = as_labeller(c(M_AB = "M_AB (A->B, full adherence)", M_BC = "M_BC (B->C, full adherence)",
                                                                         partial = "Partial-adherence relative protection"))) +
        scale_fill_manual(values = col_reg, name = "Regime") + labs(x = NULL, y = "Draws",
          title = "Sampled SAP inputs (dashed = production point value)") + th
      g2 <- ggplot(smp, aes(M_AB, M_BC, colour = regime)) + geom_point(size = 1, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.3) +
        annotate("point", x = 0.45, y = 0.45, shape = 4, size = 3, stroke = 1) +
        scale_x_log10() + scale_y_log10() + scale_colour_manual(values = col_reg, name = "Regime") +
        labs(title = "Joint draws of M_AB and M_BC (log scales; x = production point)") + th
      g3 <- ggplot(cats, aes(factor(sap_protection_age_max), share)) + geom_col(fill = "#898781", width = 0.6) +
        geom_point(aes(y = prior_prob), shape = 23, size = 2.5, fill = "white") +
        labs(x = "sap_protection_age_max (100 = no cap)", y = "Share of draws",
             title = "Protection age (bars = sample, diamonds = prior)") + th
      fig <- (g1 / (g2 | g3)) + patchwork::plot_annotation(
        title = sprintf("Provisional SAP prior (not expert-agreed): %d draws, seed %d", nrow(smp), cfg$seed),
        caption = "Trial-efficacy regime: GOAL-based; recurrence-only regime: M = 1 - s (1 - RR_rec). Shared regime indicator for M_AB and M_BC.")
      f_png <- file.path(p$figures, "prior_draws.png"); voi_assert_out(f_png, cfg)
      ggsave(f_png, fig, width = 9, height = 7.5, dpi = 150, bg = "white")
    }
    voi_manifest_update(p, cfg, run_id = cfg$run_id, country = ct, seed = cfg$seed, n_draws = cfg$n_draws, n_pilot = cfg$n_pilot,
                        prior = cfg$prior, prior_hash = prior_hash, prior_status = cfg$prior$status,
                        rng = "L'Ecuyer-CMRG streams (draw i = i-th nextRNGStream)",
                        r_version = R.version.string)
    voi_log(p, sprintf("01 sampling: %d draws | trial-regime share %.3f | M_AB median %.3f (95%% %.3f-%.3f) | M_BC median %.3f (95%% %.3f-%.3f) | truncation rejections %d",
                       nrow(smp), mean(smp$regime == "trial_efficacy"), median(smp$M_AB), quantile(smp$M_AB, 0.025), quantile(smp$M_AB, 0.975),
                       median(smp$M_BC), quantile(smp$M_BC, 0.025), quantile(smp$M_BC, 0.975), sum(smp$n_truncation_rejections)))
  }
  invisible(draws)
}

if (identical(get0("VOI_CURRENT_STAGE", ifnotfound = ""), "sample") || !exists("VOI_CURRENT_STAGE")) voi_run_sampling(VOI_CFG)
