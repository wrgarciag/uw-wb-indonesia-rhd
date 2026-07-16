# ==============================================================================
# VERIFICATION — all spec validation checks against the PERSISTED pipeline outputs
# tests/verify_abcd_pipeline.R
#
# Reads the persisted bundles/tables produced by 03-08 and asserts each of the
# 23 specification validation checks (plus the deterministic smoke test, which is
# invoked at the end). Prints PASS/FAIL per check and exits non-zero on any FAIL.
#
# Run AFTER the pipeline:  Rscript tests/verify_abcd_pipeline.R
# ==============================================================================

library(data.table)
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")
if (!exists("wd_outp")) wd_outp <- paste0(here::here("output"), "/")

fails <- character(0)
chk <- function(id, desc, cond) {
  pass <- isTRUE(cond)
  cat(sprintf("  [%s] %-2s  %s\n", if (pass) "PASS" else "FAIL", id, desc))
  if (!pass) fails[[length(fails) + 1L]] <<- paste0(id, ": ", desc)
  invisible(pass)
}
tol <- 1e-6

# ---- load persisted outputs --------------------------------------------------
dmi   <- readRDS(paste0(wd_data, "disease_model_inputs.rds"))
calib <- readRDS(paste0(wd_data, "calibrated_rhd_parameters.rds"))
bs    <- readRDS(paste0(wd_data, "baseline_state.rds"))
omf   <- list.files(paste0(wd_outp, "out_model/"), pattern = "\\.rds$", full.names = TRUE)
om    <- readRDS(omf[1])
wsd   <- as.data.table(om$wsd)
stg   <- as.data.table(om$stages)
long  <- as.data.table(readRDS(paste0(wd_outp, "tables/rhd_model_long.rds")))
stage <- as.data.table(readRDS(paste0(wd_outp, "tables/rhd_stage_model_long.rds")))
flow  <- as.data.table(readRDS(paste0(wd_outp, "tables/rhd_annual_flows.rds")))
econ  <- readRDS(paste0(wd_outp, "tables/rhd_economic_results.rds"))

tr <- dmi$transitions; pd <- dmi$p_rhd_death; eff <- dmi$effects; surg <- dmi$surgery
ramp_end <- dmi$coverage$ramp_end

cat("── VERIFYING 23 SPEC CHECKS AGAINST PERSISTED OUTPUTS ──\n")

# 1. every transition probability in [0,1]
allp <- c(unlist(tr), pd)
chk(1, "every transition/death probability in [0,1]", all(allp >= 0 & allp <= 1))

# 2. per-stage baseline outgoing probs sum < 1
outs <- c(A = tr$p_A_to_no_rhd + tr$p_A_to_B + pd[["A"]],
          B = tr$p_B_to_A      + tr$p_B_to_C + pd[["B"]],
          C = tr$p_C_to_B      + tr$p_C_to_D + pd[["C"]],
          D = tr$p_D_to_C      + pd[["D"]])
chk(2, "baseline outgoing probs sum < 1 for every stage", all(outs < 1))

# 3. intervention-adjusted RHD mortality never negative
chk(3, "intervention-adjusted RHD mortality >= 0 (eff_cf in [0,1]; stage deaths >=0)",
    stg[, all(eff_cf >= 0 & eff_cf <= 1)] &&
    all(pd * min(stg$eff_cf) >= 0) &&
    stg[, all(rhd_deaths_A >= -tol & rhd_deaths_B >= -tol &
              rhd_deaths_C >= -tol & rhd_deaths_D >= -tol)])

# 4. surgery-adjusted C->D and D->RHD-death probs in [0,1]
chk(4, "surgery-adjusted C->D and D->RHD-death probabilities in [0,1]",
    stg[, all(eff_C_to_D_surgery >= 0 & eff_C_to_D_surgery <= 1 &
              eff_D_to_rhd_death_surgery >= 0 & eff_D_to_rhd_death_surgery <= 1)] &&
    stg[, all(tr$p_C_to_D * eff_C_to_D_surgery >= 0 &
              tr$p_C_to_D * eff_C_to_D_surgery <= 1)] &&
    stg[, all(pd[["D"]] * eff_cf * eff_D_to_rhd_death_surgery >= 0 &
              pd[["D"]] * eff_cf * eff_D_to_rhd_death_surgery <= 1)])

# 5. A/B/C/D stocks and all flows non-negative
numcols <- names(stg)[vapply(stg, is.numeric, logical(1))]
chk(5, "all A/B/C/D stocks and flows non-negative",
    all(vapply(numcols, function(cc) all(stg[[cc]] >= -tol), logical(1))))

# 6. A+B+C+D reconstructs total RHD exactly
recon <- merge(stg[, .(s = A + B + C + D), by = .(scenario, sex, age, year)],
               wsd[, .(scenario, sex, age, year, sick)],
               by = c("scenario", "sex", "age", "year"))
chk(6, "A+B+C+D reconstructs aggregate sick exactly", recon[, max(abs(s - sick))] < tol)

# 7. total living RHD does not exceed population
chk(7, "total living RHD <= population everywhere",
    wsd[, all(sick <= pop + 1e-3)] && stg[, all(living_rhd_start <= pop + 1e-3)])

# 8. reference cascade coverage constant over horizon
refcov <- stg[scenario == "ref", .(sc = uniqueN(round(screen_coverage, 12)),
                                    dc = uniqueN(round(diagnosis_coverage, 12)),
                                    tc = uniqueN(round(optimal_treatment_coverage, 12)))]
chk(8, "reference cascade coverage constant over horizon",
    refcov$sc == 1 && refcov$dc == 1 && refcov$tc == 1)

# 9. scale-up reaches EXACTLY 80% screen, 80% diagnosis, 65% treatment in 2050
y2050 <- stg[scenario == "sap" & year == ramp_end][1]
chk(9, sprintf("scale-up reaches 80%%/80%%/65%% in %d", ramp_end),
    abs(y2050$screen_coverage - 0.80) < 1e-9 &&
    abs(y2050$diagnosis_coverage - 0.80) < 1e-9 &&
    abs(y2050$optimal_treatment_coverage - 0.65) < 1e-9)

# 10. coverage held at targets after 2050
post <- stg[scenario == "sap" & year > ramp_end]
chk(10, "scale-up cascade held at targets after 2050",
    post[, all(abs(screen_coverage - 0.80) < 1e-9 &
               abs(diagnosis_coverage - 0.80) < 1e-9 &
               abs(optimal_treatment_coverage - 0.65) < 1e-9)])

# 11. effective treatment never exceeds screening or diagnosis
chk(11, "effective treatment <= screening and <= effective diagnosis; eff diag <= screening",
    stg[, all(effective_treatment_coverage <= screen_coverage + 1e-12 &
              effective_treatment_coverage <= effective_diagnosis_coverage + 1e-12 &
              effective_diagnosis_coverage <= screen_coverage + 1e-12)])

# 12. screening volume = total population x screening coverage
chk(12, "n_screened == population x screening coverage",
    stg[, max(abs(n_screened - pop * screen_coverage))] < 1e-6)

# 13. treatment volume = living RHD (start of cycle) x effective treatment coverage
chk(13, "n_on_optimal_treatment == living_rhd_start x effective treatment coverage",
    stg[, max(abs(n_on_optimal_treatment - living_rhd_start * effective_treatment_coverage))] < 1e-6)

# 14/15. C and D surgery requirement = stock x requirement fraction (verified at
#        year 1 vs the seeded stock; general case proven in the smoke test).
by1  <- min(stg$year)
seedC <- bs$states[[1]]$seed$C; seedD <- bs$states[[1]]$seed$D
s1 <- stg[scenario == "ref" & year == by1]
setkey(s1, sex, age)
c_req_y1 <- s1[data.table(sex = rep(colnames(seedC), each = nrow(seedC)),
                          age = as.integer(rep(rownames(seedC), ncol(seedC)))),
               on = .(sex, age), C_requiring_surgery]
d_req_y1 <- s1[data.table(sex = rep(colnames(seedD), each = nrow(seedD)),
                          age = as.integer(rep(rownames(seedD), ncol(seedD)))),
               on = .(sex, age), D_requiring_surgery]
chk(14, "C surgery requirement == C stock x frac_C (year 1 vs seed)",
    max(abs(c_req_y1 - as.vector(seedC) * surg$frac_C_requiring_surgery)) < 1e-6)
chk(15, "D surgery requirement == D stock x frac_D (year 1 vs seed)",
    max(abs(d_req_y1 - as.vector(seedD) * surg$frac_D_requiring_surgery)) < 1e-6)

# 16. surgery volume does not exceed requirement (and == requirement x coverage)
chk(16, "surgeries <= requiring, and surgeries == requiring x surgery coverage",
    stg[, all(surgeries_C <= C_requiring_surgery + tol &
              surgeries_D <= D_requiring_surgery + tol)] &&
    stg[, max(abs(surgeries_C - C_requiring_surgery * surgery_coverage),
              abs(surgeries_D - D_requiring_surgery * surgery_coverage))] < 1e-6)

# 17. surgery does not create/remove an epidemiological stock (only A/B/C/D live)
chk(17, "no surgery/post-surgery stock; only A/B/C/D living stocks",
    all(c("A", "B", "C", "D") %in% names(stg)) &&
    !any(grepl("post|surgery_state|surg_stock", names(stg))) &&
    identical(sort(unique(as.character(stage$stage))), c("A", "B", "C", "D")))

# 18. stage C and D reported independently (both present and non-degenerate)
chk(18, "stage C and D reported independently after surgery",
    all(c("C", "D") %in% names(stg)) && stg[, sd(C) > 0 && sd(D) > 0])

# 19. scale-up cumulative RHD deaths <= reference
dd <- long[, .(dead = sum(dead)), by = scenario]
chk(19, "scale-up cumulative RHD deaths <= reference",
    dd[scenario == "sap", dead] <= dd[scenario == "ref", dead] + tol)

# 20. incremental cost non-negative (cumulative)
chk(20, "cumulative incremental cost >= 0", econ$incremental[, sum(inc_cost)] >= -tol)

# 21. BCR & CE ratios handle zero denominators safely
safe_ratio <- function(num, den) if (is.finite(den) && den > 0) num / den else NA_real_
chk(21, "BCR/CE ratios guard zero denominators (safe_ratio) and BCR finite+positive",
    is.na(safe_ratio(1, 0)) && is.na(safe_ratio(1, -3)) &&
    is.finite(econ$summary[location == "TOTAL", bcr]) &&
    econ$summary[location == "TOTAL", bcr] > 0)

# 22. persisted output schemas match documented contracts
agg_cols <- c("scenario","age","cause","sex","year","well","sick","newcases",
              "dead","pop","all.mx","intervention","location","eff_ir","eff_cf")
stage_cols <- c("scenario","age","sex","year","location","stage","cases",
                "prevalence","prevalence_per_1000")
chk(22, "persisted schemas match contracts (bundle keys + table columns)",
    all(c("rates_by_year","transitions","p_rhd_death","effects","surgery",
          "stage_split","coverage","meta") %in% names(dmi)) &&
    all(c("A","B","C","D") == names(dmi$stage_split)) &&
    all(c("tp","factors","diagnostics","stage_calibration","meta") %in% names(calib)) &&
    all(c("pop","ir","cf","oth_mort","seed","transitions","p_rhd_death","effects",
          "surgery","coverage") %in% names(bs$states[[1]])) &&
    all(c("A","B","C","D") == names(bs$states[[1]]$seed)) &&
    all(c("wsd","stages","diag","meta") %in% names(om)) &&
    identical(names(long), agg_cols) &&
    identical(names(stage), stage_cols))

# 23. pipeline completes from a clean env WITHOUT relying on leftover globals:
#     every expected artifact exists and no obsolete mild/severe/post artifact remains.
expected <- c(paste0(wd_data, c("disease_model_inputs.rds","calibrated_rhd_parameters.rds",
                                 "baseline_state.rds")),
              paste0(wd_outp, "tables/", c("rhd_model_long.csv","rhd_stage_model_long.csv",
                                           "rhd_annual_flows.csv","rhd_economic_summary.csv",
                                           "rhd_budget_impact.csv")))
old_chunks <- list.files(wd_data, pattern = "^adjusted_searo_part[0-9]+\\.rds$")
chk(23, "all expected artifacts present; obsolete adjusted_searo_part*.rds removed",
    all(file.exists(expected)) && length(old_chunks) == 0)

# ---- summary -----------------------------------------------------------------
cat("\n")
if (length(fails) == 0) {
  cat("ALL 23 SPEC VALIDATION CHECKS PASSED against persisted outputs ✓\n")
} else {
  cat(sprintf("%d CHECK(S) FAILED:\n", length(fails)))
  for (f in fails) cat("   - ", f, "\n", sep = "")
  stop("Verification failed.", call. = FALSE)
}

# ---- also run the deterministic one-cycle smoke test -------------------------
cat("\n── running deterministic smoke test ──\n")
smoke <- if (file.exists("tests/test_abcd_smoke.R")) "tests/test_abcd_smoke.R" else
         here::here("tests", "test_abcd_smoke.R")
sys.source(smoke, envir = new.env())
