# ==============================================================================
# SHARED A/B/C/D ONE-CYCLE ENGINE  (single source of truth)
# R/abcd_engine.R
#
# One annual A/B/C/D Markov cycle used by BOTH the production runner
# (06_run_prevention_model.R) and the structural calibration
# (04b_calibrate_structural.R), so there is exactly one implementation of the
# stage dynamics.
#
# Disease structure: No RHD -> A <-> B <-> C -> D -> RHD death, with competing
# other-cause death from every living stage; incident RHD enters stage A.
#
# COMPETING RISKS ON THE HAZARD SCALE (replaces the former additive-probability +
# pmax(...,0) flooring). Within a cycle each living stage faces several competing
# annual events (progression, regression, RHD death, other-cause death). Given
# their individual annual probabilities p_i (each p_i = 1 - exp(-h_i)):
#     H      = sum_i h_i,  h_i = -log(1 - p_i)
#     p_any  = 1 - exp(-H)                 (leave the stage via ANY event)
#     p_i^out = p_any * h_i / H            (leave via event i; 0 if H = 0)
#     stay    = exp(-H)                    (a PROPER residual, never floored)
# so sum_i p_i^out + stay == 1 exactly and every outflow is non-negative without
# clipping. Surgery is a SERVICE: it scales the C->D and D->RHD-death event
# probabilities and is traced, but never removes stock.
#
# SAP reduces the RHD-death event probability of EVERY stage by the mortality
# multiplier sap_mult = 1 - sap_rrr_rhd_death * effective_treatment_coverage.
#
# Inputs (per cycle):
#   A,B,C,D : [age x sex] living-stage stocks entering the cycle
#   popm    : [age x sex] population this year
#   irm     : [age x sex] CALIBRATED incidence probability (No RHD -> A)
#   oth     : [age x sex] background (non-RHD) other-cause mortality probability
#   tr      : list of transition probabilities (p_A_to_no_rhd, p_A_to_B, p_B_to_A,
#             p_B_to_C, p_C_to_B, p_C_to_D, p_D_to_C) — scalars
#   pd      : per-stage RHD-death probabilities c(A,B,C,D) — scalars
#   eff     : list(sap_rrr_rhd_death, eff_surgery_C_to_D, eff_surgery_D_to_rhd_death)
#   surg    : list(frac_C_requiring_surgery, frac_D_requiring_surgery)
#   cov     : list(treatment = effective optimal-treatment coverage [scalar],
#                  surgery   = surgery coverage [scalar])
#
# Returns a named list of [age x sex] matrices: next-cycle stocks (A,B,C,D,
# pre-ageing), incident inflow, every flow, per-stage RHD & other deaths, the
# aggregate deaths, the surgery trace, well/sick end stocks, and the scalar
# effect multipliers (sap_mult, reach_C, reach_D, eff_C_to_D_surg,
# eff_D_death_surg). Ageing (age_shift) and melting to tidy tables stay in 06.
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

abcd_one_cycle <- function(A, B, C, D, popm, irm, oth, tr, pd, eff, surg, cov, light = FALSE) {
  # start-of-cycle prevalent RHD and susceptible ("No RHD") pools
  rhd_start <- A + B + C + D
  no_rhd    <- pmax(popm - rhd_start, 0)

  # incident RHD -> stage A (eff_ir = 1: secondary prevention does not cut incidence)
  eff_ir    <- 1
  new_rhd_A <- no_rhd * irm * eff_ir

  # intervention multipliers ---------------------------------------------------
  sap_mult <- 1 - eff$sap_rrr_rhd_death * cov$treatment              # SAP RHD-mortality multiplier
  reach_C  <- surg$frac_C_requiring_surgery * cov$surgery            # population-average surgery reach
  reach_D  <- surg$frac_D_requiring_surgery * cov$surgery
  eff_C_to_D_surg  <- 1 - eff$eff_surgery_C_to_D         * reach_C
  eff_D_death_surg <- 1 - eff$eff_surgery_D_to_rhd_death * reach_D

  # effective per-stage RHD-death / C->D event probabilities
  pA_death <- pd[["A"]] * sap_mult
  pB_death <- pd[["B"]] * sap_mult
  pC_death <- pd[["C"]] * sap_mult
  pD_death <- pd[["D"]] * sap_mult * eff_D_death_surg
  pC_to_D  <- tr$p_C_to_D * eff_C_to_D_surg

  # --- stage A: {-> No RHD, -> B, RHD death, other death} ----------------------
  sA <- abcd_compete(list(no_rhd = tr$p_A_to_no_rhd, toB = tr$p_A_to_B,
                          rhd = pA_death, oth = oth))
  A_to_no_rhd <- A * sA$out$no_rhd; A_to_B <- A * sA$out$toB
  A_death_rhd <- A * sA$out$rhd;    A_death_oth <- A * sA$out$oth
  A_stay      <- A * sA$stay

  # --- stage B: {-> A, -> C, RHD death, other death} --------------------------
  sB <- abcd_compete(list(toA = tr$p_B_to_A, toC = tr$p_B_to_C,
                          rhd = pB_death, oth = oth))
  B_to_A <- B * sB$out$toA; B_to_C <- B * sB$out$toC
  B_death_rhd <- B * sB$out$rhd; B_death_oth <- B * sB$out$oth
  B_stay <- B * sB$stay

  # --- stage C: {-> B, -> D (surgery lowers), RHD death, other death} ---------
  sC <- abcd_compete(list(toB = tr$p_C_to_B, toD = pC_to_D,
                          rhd = pC_death, oth = oth))
  C_to_B <- C * sC$out$toB; C_to_D <- C * sC$out$toD
  C_death_rhd <- C * sC$out$rhd; C_death_oth <- C * sC$out$oth
  C_stay <- C * sC$stay

  # --- stage D: {-> C, RHD death (SAP x surgery), other death} ----------------
  sD <- abcd_compete(list(toC = tr$p_D_to_C, rhd = pD_death, oth = oth))
  D_to_C <- D * sD$out$toC
  D_death_rhd <- D * sD$out$rhd; D_death_oth <- D * sD$out$oth
  D_stay <- D * sD$stay

  # --- next-cycle stocks (year-end, pre-ageing) -------------------------------
  A_next <- A_stay + B_to_A + new_rhd_A
  B_next <- B_stay + A_to_B + C_to_B
  C_next <- C_stay + B_to_C + D_to_C
  D_next <- D_stay + C_to_D

  # FAST PATH for the structural calibration hot loop: same stage arithmetic, but
  # skip the surgery trace / program volumes / 40-element output list. Returns
  # only the next stocks + RHD deaths (mathematically identical to the full path).
  if (isTRUE(light))
    return(list(A = A_next, B = B_next, C = C_next, D = D_next,
                rhd_deaths = A_death_rhd + B_death_rhd + C_death_rhd + D_death_rhd,
                sick_end = A_next + B_next + C_next + D_next))

  # --- surgery TRACE (service; on START stocks; NEVER subtracted from a stock) -
  C_req_surg  <- C * surg$frac_C_requiring_surgery
  D_req_surg  <- D * surg$frac_D_requiring_surgery
  surgeries_C <- C_req_surg * cov$surgery
  surgeries_D <- D_req_surg * cov$surgery

  # --- deaths -----------------------------------------------------------------
  rhd_deaths   <- A_death_rhd + B_death_rhd + C_death_rhd + D_death_rhd
  other_deaths <- A_death_oth + B_death_oth + C_death_oth + D_death_oth   # among RHD
  all_mx       <- rhd_deaths + oth * popm            # whole population all-cause deaths

  sick_end <- A_next + B_next + C_next + D_next
  well_end <- pmax(popm - sick_end, 0)

  list(
    A = A_next, B = B_next, C = C_next, D = D_next,
    well_end = well_end, sick_end = sick_end,
    rhd_start = rhd_start, no_rhd = no_rhd, new_rhd_A = new_rhd_A,
    A_to_no_rhd = A_to_no_rhd, A_to_B = A_to_B, B_to_A = B_to_A, B_to_C = B_to_C,
    C_to_B = C_to_B, C_to_D = C_to_D, D_to_C = D_to_C,
    rhd_deaths_A = A_death_rhd, rhd_deaths_B = B_death_rhd,
    rhd_deaths_C = C_death_rhd, rhd_deaths_D = D_death_rhd,
    other_deaths_A = A_death_oth, other_deaths_B = B_death_oth,
    other_deaths_C = C_death_oth, other_deaths_D = D_death_oth,
    rhd_deaths = rhd_deaths, other_deaths = other_deaths, all_mx = all_mx,
    C_requiring_surgery = C_req_surg, D_requiring_surgery = D_req_surg,
    surgeries_C = surgeries_C, surgeries_D = surgeries_D,
    total_surgeries = surgeries_C + surgeries_D,
    sap_mult = sap_mult, reach_C = reach_C, reach_D = reach_D,
    eff_ir = eff_ir, eff_C_to_D_surg = eff_C_to_D_surg, eff_D_death_surg = eff_D_death_surg)
}
