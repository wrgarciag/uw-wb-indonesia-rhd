# ==============================================================================
# SHARED A/B/C/D ONE-CYCLE ENGINE  (single source of truth)
# R/abcd_engine.R
#
# One annual A/B/C/D Markov cycle used by BOTH the production runner
# (06_run_prevention_model.R) and the structural calibration
# (04b_calibrate_structural.R), so there is exactly one implementation of the
# stage dynamics. The small projection wrapper (abcd_project) and the ageing /
# burn-in helpers live here too, so calibration and production share the SAME
# initialisation, ageing and observation logic.
#
# Disease structure: No RHD -> A <-> B <-> C -> D -> RHD death, with competing
# other-cause death from every living stage (and from the No-RHD pool); incident
# RHD enters stage A.
#
# COMPETING RISKS ON THE HAZARD SCALE. Within a cycle each living stage faces
# several competing annual events (progression, regression, RHD death, other-cause
# death). Given their individual annual probabilities p_i (each p_i = 1 - exp(-h_i)):
#     H      = sum_i h_i,  h_i = -log(1 - p_i)
#     p_any  = 1 - exp(-H)                 (leave the stage via ANY event)
#     p_i^out = p_any * h_i / H            (leave via event i; 0 if H = 0)
#     stay    = exp(-H)                    (a PROPER residual, never floored)
# so sum_i p_i^out + stay == 1 exactly and every outflow is non-negative without
# clipping. The No-RHD pool faces {incident RHD, other-cause death} the same way,
# so all-cause deaths are RECONSTRUCTED from flows (well_deaths + other deaths
# among RHD + RHD deaths), not as rhd_deaths + oth x population.
#
# Surgery is a SERVICE: it scales the C->D and D->RHD-death event probabilities
# (population-average form, 1 - RRR x reach) and is traced, but never removes stock.
#
# SAP (secondary antibiotic prophylaxis) — SELECTABLE MECHANISM (eff$sap_mechanism):
#   "progression"      (primary) : pathway-specific HAZARD-RATIO effects on the
#                                  supported progression transitions. For pathway
#                                  i->j with hazard ratio HR_ij and effective
#                                  protection c (= treatment coverage x adherence):
#                                     h_ij(c) = h_ij(0) * [1 - c * (1 - HR_ij)]
#                                  i.e. an annual population-average hazard
#                                  multiplier (an explicit approximation, not an
#                                  exact mixture of annual risks), then converted to
#                                  an engine probability p = 1 - (1 - p0)^mult.
#                                  NO direct mortality effect (sap_mult == 1).
#   "direct_mortality" (legacy)  : the former comparator — every stage's RHD-death
#                                  probability is multiplied by
#                                  1 - sap_rrr_rhd_death x c (probability scale).
# Mortality STRUCTURE is carried entirely by `pd`: with pd[A]=pd[B]=pd[C]=0 (the
# D-only specification) the hazard split gives EXACTLY zero RHD deaths in A-C
# (h = -log1p(-0) = 0 -> outflow 0), while background deaths remain in every state.
#
# Inputs (per cycle):
#   A,B,C,D : [age x sex] living-stage stocks entering the cycle (START stocks)
#   popm    : [age x sex] population this year (exogenous backbone)
#   irm     : [age x sex] CALIBRATED incidence probability per SUSCEPTIBLE (No RHD -> A)
#   oth     : [age x sex] background (non-RHD) other-cause mortality probability
#   tr      : list of transition probabilities (p_A_to_no_rhd, p_A_to_B, p_B_to_A,
#             p_B_to_C, p_C_to_B, p_C_to_D, p_D_to_C) — scalars OR [age x sex] matrices
#   pd      : per-stage RHD-death probabilities c(A,B,C,D) — scalars (exact zeros allowed)
#   eff     : list(sap_rrr_rhd_death, eff_surgery_C_to_D, eff_surgery_D_to_rhd_death,
#                  sap_mechanism ["direct_mortality" if absent], sap_hr_A_to_B,
#                  sap_hr_B_to_C, sap_hr_C_to_D [1 if absent])
#   surg    : list(frac_C_requiring_surgery, frac_D_requiring_surgery)
#   cov     : list(treatment = effective optimal-treatment coverage [scalar],
#                  surgery   = surgery coverage [scalar],
#                  adherence = annual population-average adherence modifier [1 if absent])
#
# Returns a named list of [age x sex] matrices: next-cycle stocks (A,B,C,D,
# pre-ageing), incident inflow, every flow, per-stage RHD & other deaths, deaths in
# the No-RHD pool, the aggregate deaths, the surgery trace, well/sick end stocks,
# D_start (start-of-cycle D stock, the denominator of the realised D fatality) and
# the scalar effect multipliers. Ageing is abcd_age_shift(); melting stays in 06.
# ==============================================================================

# multi-way competing-risk split on the hazard scale.
# `probs` : named list of annual event probabilities (scalar or [age x sex] matrix).
# returns list(out = named list of per-event OUTFLOW probabilities, stay = matrix).
abcd_compete <- function(probs) {
  n <- length(probs)
  hz <- vector("list", n); H <- 0
  for (i in seq_len(n)) { h <- -log1p(-pmin(pmax(probs[[i]], 0), 1 - 1e-12)); hz[[i]] <- h; H <- H + h }
  stay <- exp(-H)
  pany <- 1 - stay
  scale <- pany / ifelse(H > 0, H, 1)         # where H == 0, pany == 0 so outflow == 0
  out <- vector("list", n); names(out) <- names(probs)
  for (i in seq_len(n)) out[[i]] <- scale * hz[[i]]
  list(out = out, stay = stay)
}

# hazard-scale multiplier on an annual probability: h -> h * m  <=>  p -> 1 - (1-p)^m
abcd_hazard_scale <- function(p, m) 1 - (1 - pmin(pmax(p, 0), 1 - 1e-12))^m

# SAP effect multipliers for one cycle (mechanism-aware). Returns the mortality
# multiplier (probability scale, legacy mechanism only), the three progression
# HAZARD multipliers, the effective protection c and the mechanism label.
abcd_sap_effects <- function(eff, cov) {
  mech <- if (is.null(eff$sap_mechanism)) "direct_mortality" else as.character(eff$sap_mechanism)
  adh  <- if (is.null(cov$adherence)) 1 else cov$adherence
  c_eff <- cov$treatment * adh                        # effective protection in the eligible pool
  if (mech == "direct_mortality") {
    list(mechanism = mech, c_eff = c_eff,
         sap_mult = 1 - eff$sap_rrr_rhd_death * c_eff,
         m_A_to_B = 1, m_B_to_C = 1, m_C_to_D = 1)
  } else if (mech == "progression") {
    hr <- function(nm) if (is.null(eff[[nm]])) 1 else eff[[nm]]
    list(mechanism = mech, c_eff = c_eff,
         sap_mult = 1,                                  # NO direct mortality effect in progression mode
         m_A_to_B = 1 - c_eff * (1 - hr("sap_hr_A_to_B")),
         m_B_to_C = 1 - c_eff * (1 - hr("sap_hr_B_to_C")),
         m_C_to_D = 1 - c_eff * (1 - hr("sap_hr_C_to_D")))
  } else stop("Unknown eff$sap_mechanism '", mech, "' (use 'progression' or 'direct_mortality').",
              call. = FALSE)
}

abcd_one_cycle <- function(A, B, C, D, popm, irm, oth, tr, pd, eff, surg, cov, light = FALSE) {
  # start-of-cycle prevalent RHD and susceptible ("No RHD") pools. A stock larger
  # than the population is a VIOLATION: it is measured (pop_deficit) and surfaced by
  # the callers, not hidden (the pmax only keeps the arithmetic finite).
  rhd_start   <- A + B + C + D
  pop_deficit <- pmax(rhd_start - popm, 0)
  no_rhd      <- pmax(popm - rhd_start, 0)

  # No-RHD pool: {incident RHD, other-cause death} compete on the hazard scale
  # (eff_ir = 1: secondary prevention does not cut incidence)
  eff_ir <- 1
  sW <- abcd_compete(list(inc = irm * eff_ir, oth = oth))
  new_rhd_A   <- no_rhd * sW$out$inc
  well_deaths <- no_rhd * sW$out$oth

  # intervention multipliers ---------------------------------------------------
  sap      <- abcd_sap_effects(eff, cov)
  sap_mult <- sap$sap_mult                                           # RHD-mortality multiplier (legacy mechanism)
  reach_C  <- surg$frac_C_requiring_surgery * cov$surgery            # population-average surgery reach
  reach_D  <- surg$frac_D_requiring_surgery * cov$surgery
  eff_C_to_D_surg  <- 1 - eff$eff_surgery_C_to_D         * reach_C
  eff_D_death_surg <- 1 - eff$eff_surgery_D_to_rhd_death * reach_D

  # effective per-stage RHD-death probabilities (exact zeros stay exact)
  pA_death <- pd[["A"]] * sap_mult
  pB_death <- pd[["B"]] * sap_mult
  pC_death <- pd[["C"]] * sap_mult
  pD_death <- pd[["D"]] * sap_mult * eff_D_death_surg
  # effective progression probabilities: surgery (probability scale, legacy service
  # form) then SAP progression effect on the HAZARD scale
  pA_to_B <- abcd_hazard_scale(tr$p_A_to_B, sap$m_A_to_B)
  pB_to_C <- abcd_hazard_scale(tr$p_B_to_C, sap$m_B_to_C)
  pC_to_D <- abcd_hazard_scale(tr$p_C_to_D * eff_C_to_D_surg, sap$m_C_to_D)

  # --- stage A: {-> No RHD, -> B, RHD death, other death} ----------------------
  sA <- abcd_compete(list(no_rhd = tr$p_A_to_no_rhd, toB = pA_to_B,
                          rhd = pA_death, oth = oth))
  A_to_no_rhd <- A * sA$out$no_rhd; A_to_B <- A * sA$out$toB
  A_death_rhd <- A * sA$out$rhd;    A_death_oth <- A * sA$out$oth
  A_stay      <- A * sA$stay

  # --- stage B: {-> A, -> C, RHD death, other death} --------------------------
  sB <- abcd_compete(list(toA = tr$p_B_to_A, toC = pB_to_C,
                          rhd = pB_death, oth = oth))
  B_to_A <- B * sB$out$toA; B_to_C <- B * sB$out$toC
  B_death_rhd <- B * sB$out$rhd; B_death_oth <- B * sB$out$oth
  B_stay <- B * sB$stay

  # --- stage C: {-> B, -> D (surgery + SAP lower), RHD death, other death} ----
  sC <- abcd_compete(list(toB = tr$p_C_to_B, toD = pC_to_D,
                          rhd = pC_death, oth = oth))
  C_to_B <- C * sC$out$toB; C_to_D <- C * sC$out$toD
  C_death_rhd <- C * sC$out$rhd; C_death_oth <- C * sC$out$oth
  C_stay <- C * sC$stay

  # --- stage D: {-> C, RHD death (x surgery [x SAP, legacy]), other death} -----
  sD <- abcd_compete(list(toC = tr$p_D_to_C, rhd = pD_death, oth = oth))
  D_to_C <- D * sD$out$toC
  D_death_rhd <- D * sD$out$rhd; D_death_oth <- D * sD$out$oth
  D_stay <- D * sD$stay

  # --- next-cycle stocks (year-end, pre-ageing) -------------------------------
  A_next <- A_stay + B_to_A + new_rhd_A
  B_next <- B_stay + A_to_B + C_to_B
  C_next <- C_stay + B_to_C + D_to_C
  D_next <- D_stay + C_to_D

  rhd_deaths <- A_death_rhd + B_death_rhd + C_death_rhd + D_death_rhd

  # FAST PATH for the structural calibration hot loop: same stage arithmetic, but
  # skip the surgery trace / program volumes / full output list. Returns only the
  # quantities the calibration scores (mathematically identical to the full path).
  if (isTRUE(light))
    return(list(A = A_next, B = B_next, C = C_next, D = D_next,
                rhd_deaths = rhd_deaths, rhd_deaths_D = D_death_rhd,
                new_rhd_A = new_rhd_A, A_to_no_rhd = A_to_no_rhd,
                other_deaths = A_death_oth + B_death_oth + C_death_oth + D_death_oth,
                sick_end = A_next + B_next + C_next + D_next,
                pop_deficit = pop_deficit))

  # --- surgery TRACE (service; on START stocks; NEVER subtracted from a stock) -
  C_req_surg  <- C * surg$frac_C_requiring_surgery
  D_req_surg  <- D * surg$frac_D_requiring_surgery
  surgeries_C <- C_req_surg * cov$surgery
  surgeries_D <- D_req_surg * cov$surgery

  # --- deaths (all-cause RECONSTRUCTED from competing-risk flows) -------------
  other_deaths <- A_death_oth + B_death_oth + C_death_oth + D_death_oth   # among RHD
  all_mx       <- rhd_deaths + other_deaths + well_deaths                # whole population

  sick_end <- A_next + B_next + C_next + D_next
  well_end <- pmax(popm - sick_end, 0)

  list(
    A = A_next, B = B_next, C = C_next, D = D_next,
    well_end = well_end, sick_end = sick_end,
    rhd_start = rhd_start, no_rhd = no_rhd, new_rhd_A = new_rhd_A,
    D_start = D, pop_deficit = pop_deficit,
    A_to_no_rhd = A_to_no_rhd, A_to_B = A_to_B, B_to_A = B_to_A, B_to_C = B_to_C,
    C_to_B = C_to_B, C_to_D = C_to_D, D_to_C = D_to_C,
    rhd_deaths_A = A_death_rhd, rhd_deaths_B = B_death_rhd,
    rhd_deaths_C = C_death_rhd, rhd_deaths_D = D_death_rhd,
    other_deaths_A = A_death_oth, other_deaths_B = B_death_oth,
    other_deaths_C = C_death_oth, other_deaths_D = D_death_oth,
    rhd_deaths = rhd_deaths, other_deaths = other_deaths,
    well_deaths = well_deaths, all_mx = all_mx,
    C_requiring_surgery = C_req_surg, D_requiring_surgery = D_req_surg,
    surgeries_C = surgeries_C, surgeries_D = surgeries_D,
    total_surgeries = surgeries_C + surgeries_D,
    sap_mult = sap_mult, sap_mechanism = sap$mechanism, sap_c_eff = sap$c_eff,
    sap_prog_mult_A_to_B = sap$m_A_to_B, sap_prog_mult_B_to_C = sap$m_B_to_C,
    sap_prog_mult_C_to_D = sap$m_C_to_D,
    reach_C = reach_C, reach_D = reach_D,
    eff_ir = eff_ir, eff_C_to_D_surg = eff_C_to_D_surg, eff_D_death_surg = eff_D_death_surg)
}

# ------------------------------------------------------------------------------
# AGEING: advance every cohort one year; the terminal row is the OPEN group
# (100+ in production; the last calibration age otherwise) and accumulates.
# Works on an [age x sex] matrix (or a one-column matrix).
# ------------------------------------------------------------------------------
abcd_age_shift <- function(M) {
  n <- nrow(M)
  N <- M * 0
  N[2:n, ] <- M[1:(n - 1), , drop = FALSE]
  N[n, ]   <- N[n, ] + M[n, ]
  N
}

# ------------------------------------------------------------------------------
# BURN-IN to the transition-implied stage equilibrium under CONSTANT inputs
# (first-window-year population, incidence, background mortality and care).
# Returns the stocks entering the first scored year and a convergence table
# (relative change of each stage total over the last cycle and vs the half-way
# state) so "equilibrium" is CHECKED, not asserted.
# ------------------------------------------------------------------------------
abcd_burn_in <- function(pop1, ir1, oth1, tr, pd, eff, surg, cov1, n_cycles, init = NULL) {
  z <- pop1 * 0
  A <- if (is.null(init)) z else init$A; B <- if (is.null(init)) z else init$B
  C <- if (is.null(init)) z else init$C; D <- if (is.null(init)) z else init$D
  half <- NULL; prev <- NULL
  for (b in seq_len(n_cycles)) {
    cyc <- abcd_one_cycle(A, B, C, D, pop1, ir1, oth1, tr, pd, eff, surg, cov1, light = TRUE)
    if (b == n_cycles) prev <- list(A = A, B = B, C = C, D = D)
    A <- abcd_age_shift(cyc$A); B <- abcd_age_shift(cyc$B)
    C <- abcd_age_shift(cyc$C); D <- abcd_age_shift(cyc$D)
    if (b == max(1L, n_cycles %/% 2L)) half <- list(A = A, B = B, C = C, D = D)
  }
  fin <- list(A = A, B = B, C = C, D = D)
  rel <- function(x, y) { sx <- sum(x); if (sx <= 0) 0 else abs(sx - sum(y)) / sx }
  conv <- data.frame(stage = c("A", "B", "C", "D"),
                     rel_change_last_cycle = c(rel(fin$A, prev$A), rel(fin$B, prev$B),
                                               rel(fin$C, prev$C), rel(fin$D, prev$D)),
                     rel_change_vs_half   = c(rel(fin$A, half$A), rel(fin$B, half$B),
                                               rel(fin$C, half$C), rel(fin$D, half$D)))
  list(state = fin, convergence = conv, n_cycles = n_cycles)
}

# ------------------------------------------------------------------------------
# PROJECTION WRAPPER shared by calibration (light = TRUE) and production
# (light = FALSE). ONE time convention everywhere:
#   start stock (entering year t) -> within-year flows -> end stock (pre-ageing)
#   -> age_shift -> start stock of year t+1.
# Observation of PREVALENCE for year t (prev_obs): "start" = start stock,
# "end" = end stock (pre-ageing), "mid" = (start + end)/2 (a mid-year point-
# prevalence approximation to GBD's annual prevalence). Deaths and incident cases
# are within-year FLOWS of year t.
#
#   init      : list(A,B,C,D) [age x sex] stocks ENTERING years[1]
#   pop, ir   : [age x sex x year] arrays (dim 3 == length(years))
#   oth       : [age x sex] matrix (held) or [age x sex x year] array
#   cov_year  : list(treatment, surgery, adherence) — numeric vectors, one value per year
#   tr        : scalar list, or a list whose entries are [age x sex] matrices
# Returns list(start = per-year start stocks, cycles = per-year engine output,
#              prev_obs = per-year [age x sex] prevalence observation,
#              final_state = stocks entering years[n]+1 (after ageing),
#              pop_deficit_total = sum of stock-over-population violations).
# ------------------------------------------------------------------------------
abcd_project <- function(init, pop, ir, oth, tr, pd, eff, surg, cov_year, years,
                         light = TRUE, prev_obs = c("mid", "start", "end")) {
  prev_obs <- match.arg(prev_obs)
  n <- length(years)
  stopifnot(dim(pop)[3] == n, dim(ir)[3] == n)
  oth_is_arr <- length(dim(oth)) == 3L
  if (is.null(cov_year$adherence)) cov_year$adherence <- rep(1, n)
  A <- init$A; B <- init$B; C <- init$C; D <- init$D
  starts <- vector("list", n); cycles <- vector("list", n); pobs <- vector("list", n)
  deficit <- 0
  for (t in seq_len(n)) {
    popm <- pop[, , t]; irm <- ir[, , t]
    othm <- if (oth_is_arr) oth[, , t] else oth
    if (is.null(dim(popm))) { dm <- dim(pop)[1:2]; dim(popm) <- dm; dim(irm) <- dm; dim(othm) <- dm }
    cv <- list(treatment = cov_year$treatment[t], surgery = cov_year$surgery[t],
               adherence = cov_year$adherence[t])
    cyc <- abcd_one_cycle(A, B, C, D, popm, irm, othm, tr, pd, eff, surg, cv, light = light)
    starts[[t]] <- list(A = A, B = B, C = C, D = D)
    cycles[[t]] <- cyc
    start_stock <- A + B + C + D
    pobs[[t]] <- switch(prev_obs, start = start_stock, end = cyc$sick_end,
                        mid = 0.5 * (start_stock + cyc$sick_end))
    deficit <- deficit + sum(cyc$pop_deficit)
    A <- abcd_age_shift(cyc$A); B <- abcd_age_shift(cyc$B)
    C <- abcd_age_shift(cyc$C); D <- abcd_age_shift(cyc$D)
  }
  list(start = starts, cycles = cycles, prev_obs = pobs,
       final_state = list(A = A, B = B, C = C, D = D),
       pop_deficit_total = deficit, years = years, prev_obs_convention = prev_obs)
}
