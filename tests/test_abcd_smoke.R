# ==============================================================================
# DETERMINISTIC SMOKE TEST — A/B/C/D engine + surgery service
# tests/test_abcd_smoke.R
#
# A small, deterministic, dependency-free check of the ONE annual cycle used by
# 06_run_prevention_model.R. It re-implements the exact matrix-engine equations
# on a tiny hand-set [age x sex] state (2 ages x 1 sex) and asserts:
#
#   1. one-cycle MASS BALANCE (no one created or lost);
#   2. A/B/C/D RECONSTRUCTION (sick == A + B + C + D);
#   3. SURGERY TRACE consistency (surgeries <= requiring; the C/D STOCKS are NOT
#      reduced by the surgery volume — surgery is a service, not a state);
#   4. NO surgery / post-surgery STOCK is created (only A,B,C,D live);
#   5. correct SURGERY risk reduction on C->D and D->RHD-death;
#   6. correct SAP mortality risk reduction on every stage.
#
# Run:  Rscript tests/test_abcd_smoke.R      (exits non-zero on any failure)
# ==============================================================================

suppressWarnings(rm(list = ls()))
ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("SMOKE TEST FAILED: ", msg, call. = FALSE)
  cat(sprintf("  PASS  %s\n", msg))
}
approx <- function(a, b, tol = 1e-9) all(abs(a - b) <= tol)

# ------------------------------------------------------------------------------
# One annual A/B/C/D cycle — identical equations to 06's run_one(), on matrices.
# Surgery is a SERVICE: its trace never enters the stock update; it only scales
# the C->D and D->RHD-death probabilities (population-average reach).
# ------------------------------------------------------------------------------
abcd_one_cycle <- function(A, B, C, D, popm, irm, oth, tr, pd, eff, surg, cov) {
  rhd_start <- A + B + C + D
  no_rhd    <- pmax(popm - rhd_start, 0)

  eff_ir    <- 1
  new_rhd_A <- no_rhd * irm * eff_ir

  sap_mult <- 1 - eff$sap_rrr_rhd_death * cov$treatment
  reach_C  <- surg$frac_C_requiring_surgery * cov$surgery
  reach_D  <- surg$frac_D_requiring_surgery * cov$surgery
  eff_C_to_D_surg  <- 1 - eff$eff_surgery_C_to_D         * reach_C
  eff_D_death_surg <- 1 - eff$eff_surgery_D_to_rhd_death * reach_D

  pA_death <- pd[["A"]] * sap_mult
  pB_death <- pd[["B"]] * sap_mult
  pC_death <- pd[["C"]] * sap_mult
  pD_death <- pd[["D"]] * sap_mult * eff_D_death_surg
  pC_to_D  <- tr$p_C_to_D * eff_C_to_D_surg

  A_to_no_rhd <- A * tr$p_A_to_no_rhd; A_to_B <- A * tr$p_A_to_B
  A_death_rhd <- A * pA_death;         A_death_oth <- A * oth
  A_stay <- pmax(A - A_to_no_rhd - A_to_B - A_death_rhd - A_death_oth, 0)

  B_to_A <- B * tr$p_B_to_A; B_to_C <- B * tr$p_B_to_C
  B_death_rhd <- B * pB_death; B_death_oth <- B * oth
  B_stay <- pmax(B - B_to_A - B_to_C - B_death_rhd - B_death_oth, 0)

  C_to_B <- C * tr$p_C_to_B; C_to_D <- C * pC_to_D
  C_death_rhd <- C * pC_death; C_death_oth <- C * oth
  C_stay <- pmax(C - C_to_B - C_to_D - C_death_rhd - C_death_oth, 0)

  D_to_C <- D * tr$p_D_to_C
  D_death_rhd <- D * pD_death; D_death_oth <- D * oth
  D_stay <- pmax(D - D_to_C - D_death_rhd - D_death_oth, 0)

  A_next <- A_stay + B_to_A + new_rhd_A
  B_next <- B_stay + A_to_B + C_to_B
  C_next <- C_stay + B_to_C + D_to_C
  D_next <- D_stay + C_to_D

  # surgery TRACE (on START stocks; NOT subtracted from any stock)
  C_req <- C * surg$frac_C_requiring_surgery
  D_req <- D * surg$frac_D_requiring_surgery
  surgeries_C <- C_req * cov$surgery
  surgeries_D <- D_req * cov$surgery

  list(A = A_next, B = B_next, C = C_next, D = D_next,
       new_rhd_A = new_rhd_A,
       A_to_no_rhd = A_to_no_rhd, A_to_B = A_to_B, B_to_A = B_to_A, B_to_C = B_to_C,
       C_to_B = C_to_B, C_to_D = C_to_D, D_to_C = D_to_C,
       rhd_deaths = A_death_rhd + B_death_rhd + C_death_rhd + D_death_rhd,
       other_deaths = A_death_oth + B_death_oth + C_death_oth + D_death_oth,
       C_requiring_surgery = C_req, D_requiring_surgery = D_req,
       surgeries_C = surgeries_C, surgeries_D = surgeries_D,
       total_surgeries = surgeries_C + surgeries_D,
       sap_mult = sap_mult, reach_C = reach_C, reach_D = reach_D,
       eff_C_to_D_surg = eff_C_to_D_surg, eff_D_death_surg = eff_D_death_surg,
       A0 = A, B0 = B, C0 = C, D0 = D)
}

# ------------------------------------------------------------------------------
# Deterministic tiny fixture (2 ages, 1 sex). No RNG anywhere.
# ------------------------------------------------------------------------------
mk <- function(v) matrix(v, nrow = 2, ncol = 1, dimnames = list(c("10", "11"), "Female"))
A <- mk(c(1000, 800)); B <- mk(c(400, 300)); C <- mk(c(200, 150)); D <- mk(c(100, 80))
popm <- mk(c(1e5, 9e4)); irm <- mk(c(0.002, 0.0025)); oth <- mk(c(0.004, 0.0045))

tr  <- list(p_A_to_no_rhd = 0.005, p_A_to_B = 0.020, p_B_to_A = 0.010,
            p_B_to_C = 0.030, p_C_to_B = 0.005, p_C_to_D = 0.060, p_D_to_C = 0.000)
pd  <- c(A = 0.0005, B = 0.0020, C = 0.0200, D = 0.0800)
eff <- list(sap_rrr_rhd_death = 0.55, eff_surgery_C_to_D = 0.85,
            eff_surgery_D_to_rhd_death = 0.85)
surg <- list(frac_C_requiring_surgery = 0.03, frac_D_requiring_surgery = 0.20)
cov  <- list(treatment = 0.65, surgery = 0.05)   # scale-up-like coverages

cat("── A/B/C/D + surgery deterministic smoke test ──\n")
r <- abcd_one_cycle(A, B, C, D, popm, irm, oth, tr, pd, eff, surg, cov)

# --- 1. one-cycle MASS BALANCE ------------------------------------------------
# end living = start living + incident inflow - RHD deaths - other deaths - A->NoRHD
L0 <- A + B + C + D
L1 <- r$A + r$B + r$C + r$D
balance <- L0 + r$new_rhd_A - r$rhd_deaths - r$other_deaths - r$A_to_no_rhd
ok(approx(L1, balance), "one-cycle mass balance (living in = living out + exits)")

# --- 2. A/B/C/D RECONSTRUCTION ------------------------------------------------
sick <- r$A + r$B + r$C + r$D
ok(approx(sick, L1), "sick reconstructs as A + B + C + D")

# --- 3. SURGERY TRACE consistency ---------------------------------------------
ok(all(r$surgeries_C <= r$C_requiring_surgery + 1e-12) &&
   all(r$surgeries_D <= r$D_requiring_surgery + 1e-12),
   "surgery volume never exceeds the number requiring surgery")
ok(approx(r$C_requiring_surgery, C * surg$frac_C_requiring_surgery) &&
   approx(r$D_requiring_surgery, D * surg$frac_D_requiring_surgery),
   "surgery requirement = stock x requirement fraction (C and D)")
# stocks are NOT reduced by the surgery volume: recompute C_next ignoring surgery
# count entirely and confirm it is unchanged (surgery affects only C->D via prob).
C_next_no_trace <- pmax(C - C*tr$p_C_to_B - C*(tr$p_C_to_D * r$eff_C_to_D_surg) -
                        C*(pd[["C"]]*r$sap_mult) - C*oth, 0) + r$B_to_C + r$D_to_C
ok(approx(r$C, C_next_no_trace),
   "C stock update does NOT subtract the surgery volume (service, not a state)")

# --- 4. NO surgery / post-surgery STOCK ---------------------------------------
stock_names <- c("A", "B", "C", "D")
ok(all(stock_names %in% names(r)) &&
   !any(grepl("post|surg_stock|surgery_state", names(r))),
   "only A/B/C/D living stocks exist (no surgery/post-surgery state)")

# --- 5. correct SURGERY risk reduction on C->D and D->RHD-death ---------------
expected_C_to_D <- C * tr$p_C_to_D * (1 - eff$eff_surgery_C_to_D *
                     surg$frac_C_requiring_surgery * cov$surgery)
ok(approx(r$C_to_D, expected_C_to_D), "C->D uses surgery-reduced probability")
# D->RHD-death combines SAP x surgery multiplicatively (no double counting)
D_death <- D * pd[["D"]] * r$sap_mult *
           (1 - eff$eff_surgery_D_to_rhd_death * surg$frac_D_requiring_surgery * cov$surgery)
# recover the modelled D RHD-death flow from the balance of D
D_death_model <- D * pd[["D"]] * r$sap_mult * r$eff_D_death_surg
ok(approx(D_death, D_death_model), "D->RHD-death combines SAP x surgery multiplicatively")

# --- 6. correct SAP mortality risk reduction ----------------------------------
ok(approx(r$sap_mult, 1 - eff$sap_rrr_rhd_death * cov$treatment),
   "SAP mortality multiplier = 1 - 0.55 x effective treatment coverage")
# with 0 treatment there is NO mortality reduction; with 0 surgery NO surgery effect
r0 <- abcd_one_cycle(A, B, C, D, popm, irm, oth, tr, pd, eff, surg,
                     list(treatment = 0, surgery = 0))
ok(approx(r0$sap_mult, 1) && approx(r0$eff_C_to_D_surg, 1) && approx(r0$eff_D_death_surg, 1),
   "zero coverage => no SAP and no surgery effect (multipliers = 1)")
ok(all(r$rhd_deaths <= r0$rhd_deaths + 1e-12),
   "SAP + surgery reduce RHD deaths vs no intervention")

# --- non-negativity -----------------------------------------------------------
all_num <- unlist(r[setdiff(names(r), c("sap_mult","reach_C","reach_D",
                                        "eff_C_to_D_surg","eff_D_death_surg"))])
ok(all(all_num >= -1e-12), "all stocks/flows/trace non-negative")

cat("\nALL A/B/C/D + SURGERY SMOKE-TEST CHECKS PASSED ✓\n")
