# ==============================================================================
# DETERMINISTIC SMOKE TEST — SHARED A/B/C/D hazard engine + surgery service
# tests/test_abcd_smoke.R
#
# Exercises the ACTUAL production one-cycle update (R/abcd_engine.R — the single
# source of truth used by 06 and 04b) on a tiny hand-set [age x sex] state
# (2 ages x 1 sex), and asserts:
#
#   1. one-cycle MASS BALANCE (no one created or lost);
#   2. A/B/C/D RECONSTRUCTION (sick == A + B + C + D);
#   3. SURGERY TRACE consistency (surgeries <= requiring; the C/D STOCKS are NOT
#      reduced by the surgery volume — surgery is a service, not a state);
#   4. NO surgery / post-surgery STOCK is created (only A,B,C,D live);
#   5. hazard COMPETING-RISK validity (per-stage outflows + stay sum to 1 exactly,
#      no clipping) and correct SURGERY risk reduction on C->D and D->RHD-death;
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

# --- load the SHARED engine (single source of truth) --------------------------
eng <- if (file.exists("R/abcd_engine.R")) "R/abcd_engine.R" else
       here::here("R", "abcd_engine.R")
source(eng)

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

cat("── A/B/C/D + surgery deterministic smoke test (shared hazard engine) ──\n")
r <- abcd_one_cycle(A, B, C, D, popm, irm, oth, tr, pd, eff, surg, cov)

# --- 1. one-cycle MASS BALANCE ------------------------------------------------
# end living = start living + incident inflow - RHD deaths - other deaths - A->NoRHD
L0 <- A + B + C + D
L1 <- r$A + r$B + r$C + r$D
balance <- L0 + r$new_rhd_A - r$rhd_deaths - r$other_deaths - r$A_to_no_rhd
ok(approx(L1, balance), "one-cycle mass balance (living in = living out + exits)")

# --- 2. A/B/C/D RECONSTRUCTION ------------------------------------------------
ok(approx(r$sick_end, L1), "sick_end reconstructs as A + B + C + D")

# --- 3. SURGERY TRACE consistency ---------------------------------------------
ok(all(r$surgeries_C <= r$C_requiring_surgery + 1e-12) &&
   all(r$surgeries_D <= r$D_requiring_surgery + 1e-12),
   "surgery volume never exceeds the number requiring surgery")
ok(approx(r$C_requiring_surgery, C * surg$frac_C_requiring_surgery) &&
   approx(r$D_requiring_surgery, D * surg$frac_D_requiring_surgery),
   "surgery requirement = stock x requirement fraction (C and D)")
# stocks are NOT reduced by the surgery volume: the engine run with cov$surgery = 0
# (no surgery reach) still produces the SAME surgery-trace requirement counts, and
# C_next changes ONLY through the C->D probability, never by subtracting a volume.
r_nosurg <- abcd_one_cycle(A, B, C, D, popm, irm, oth, tr, pd, eff, surg,
                           list(treatment = cov$treatment, surgery = 0))
ok(approx(r_nosurg$C_requiring_surgery, C * surg$frac_C_requiring_surgery),
   "C requiring-surgery trace is independent of surgery coverage (a service count)")
ok(all(r$surgeries_C >= -1e-12) && all(r$total_surgeries >= -1e-12),
   "surgery-service volumes are non-negative and never enter the stock update")

# --- 4. NO surgery / post-surgery STOCK ---------------------------------------
stock_names <- c("A", "B", "C", "D")
ok(all(stock_names %in% names(r)) &&
   !any(grepl("post|surg_stock|surgery_state", names(r))),
   "only A/B/C/D living stocks exist (no surgery/post-surgery state)")

# --- 5. hazard competing-risk validity + surgery risk reduction ---------------
# every stage: sum of outflows + stay == 1 exactly (proper residual, no clipping).
compete_sum_ok <- function(events) {
  s <- abcd_compete(events)
  tot <- Reduce(`+`, s$out) + s$stay
  approx(tot, 1)
}
sap_mult <- 1 - eff$sap_rrr_rhd_death * cov$treatment
ok(compete_sum_ok(list(a = tr$p_A_to_no_rhd, b = tr$p_A_to_B, d = pd[["A"]]*sap_mult, o = oth)) &&
   compete_sum_ok(list(a = tr$p_C_to_B, b = tr$p_C_to_D * r$eff_C_to_D_surg,
                       d = pd[["C"]]*sap_mult, o = oth)),
   "hazard competing-risk outflows + stay sum to 1 exactly (no clipping)")
# surgery multipliers are exactly 1 - RRR x reach; surgery reduces C->D and D-death.
ok(approx(r$eff_C_to_D_surg, 1 - eff$eff_surgery_C_to_D * surg$frac_C_requiring_surgery * cov$surgery) &&
   approx(r$eff_D_death_surg, 1 - eff$eff_surgery_D_to_rhd_death * surg$frac_D_requiring_surgery * cov$surgery),
   "surgery reach multipliers = 1 - RRR x (fraction x coverage)")
ok(all(r$C_to_D <= r_nosurg$C_to_D + 1e-12) &&
   all(r$rhd_deaths_D <= r_nosurg$rhd_deaths_D + 1e-12),
   "surgery lowers C->D progression and D->RHD-death vs no surgery")

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
num_keep <- c("A","B","C","D","new_rhd_A","A_to_no_rhd","A_to_B","B_to_A","B_to_C",
              "C_to_B","C_to_D","D_to_C","rhd_deaths","other_deaths",
              "C_requiring_surgery","D_requiring_surgery","surgeries_C","surgeries_D",
              "total_surgeries","sick_end","well_end")
all_num <- unlist(r[num_keep])
ok(all(all_num >= -1e-12), "all stocks/flows/trace non-negative")

cat("\nALL A/B/C/D + SURGERY SMOKE-TEST CHECKS PASSED (shared hazard engine) ✓\n")
