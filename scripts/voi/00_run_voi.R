# ==============================================================================
# RHD SAP VALUE-OF-INFORMATION ANALYSIS (separate from production)
# scripts/voi/00_run_voi.R  —  independent configuration and entry point
#
# Smallest SAP-focused VOI analysis for Indonesia and Uganda. Two options only:
#   ref  = reference care (production 06 'ref')
#   sap  = current screening + SAP scale-up (production 06 'sap')
# Four uncertain SAP inputs (everything else fixed at the production values):
#   M_AB   full-adherence progression-hazard multiplier on treated A->B
#   M_BC   full-adherence progression-hazard multiplier on treated B->C
#   sap_protection_age_max                (attained-age cap on the SAP benefit)
#   partial_adherence_relative_protection (protection under partial adherence)
# C->D and direct D-death SAP multipliers stay 1. Calibrated natural history, the 2026
# initial state, care-cascade inputs, surgery and unit costs are held at the saved
# production values; service volumes and costs are recomputed for every draw.
#
# ISOLATION: this script NEVER sources scripts/00_run_all.R, never re-runs production and
# never writes outside output/voi/<run_id>/<country>/. It READS (read-only):
#   data/<COUNTRY>/baseline_state.rds                    (production 05)
#   output/<COUNTRY>/tables/rhd_economic_results.rds     (production 08)
#   output/<COUNTRY>/out_model/<COUNTRY>.rds             (production 06; validation only)
#   R/abcd_engine.R                                      (production engine; validation only,
#                                                         loaded into an isolated environment)
# The VOI engine (scripts/voi/R/abcd_engine_voi.R) is always loaded into its own
# environment (or its own worker process), so production functions/constants are never
# overwritten. Production file fingerprints are checked before and after the run.
#
# RUN (from the repository root, R 4.5.1):
#   source("scripts/voi/00_run_voi.R")
# or  "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" scripts/voi/00_run_voi.R
# Stages (VOI_STAGES): sample -> validate -> pilot -> psa -> analyse. Set VOI_STAGES before
# sourcing to run a subset; set VOI_CONFIG_ONLY <- TRUE to define the configuration only
# (the stage scripts do this when sourced standalone).
#
# EXTENDING THE DRAWS LATER: keep VOI_RUN_ID, VOI_SEED and the prior unchanged and raise
# N_DRAWS. Draw i always comes from its own L'Ecuyer-CMRG stream, so draws 1..500 are
# reproduced exactly, only the new draws are simulated (chunk checkpoints are reused) and
# the analysis is re-run on the enlarged sample. A changed prior or seed is refused for an
# existing run_id (start a new run_id instead).
#
# PENDING (later extensions; documented here and in README.md only, NOT implemented):
#   * a third option "scale-up with BPG continued after B->C" plus uncertainty in the C->D
#     SAP effect: needs an on-BPG C stock (T_C), continuation / dropout rules for C, the
#     corresponding SAP service costs and a C->D multiplier < 1 in the engine;
#   * calibration sensitivity (this analysis conditions on the saved 04b calibration and the
#     bridged 2026 initial state, which were produced under the point SAP values);
#   * trial-specific expected value of SAMPLE information (EVSI).
# ==============================================================================

library(data.table)

# ------------------------------------------------------------------------------
# 0. PATHS
# ------------------------------------------------------------------------------
VOI_REPO_ROOT <- Sys.getenv("RHD_REPO_ROOT", unset = "")
if (!nzchar(VOI_REPO_ROOT))
  VOI_REPO_ROOT <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
VOI_REPO_ROOT <- normalizePath(VOI_REPO_ROOT, winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(VOI_REPO_ROOT, "scripts", "voi", "00_run_voi.R")))
  stop("VOI_REPO_ROOT (", VOI_REPO_ROOT, ") is not the repository root; set RHD_REPO_ROOT.", call. = FALSE)

VOI_DIR        <- file.path(VOI_REPO_ROOT, "scripts", "voi")
VOI_ENGINE     <- file.path(VOI_DIR, "R", "abcd_engine_voi.R")
VOI_FUNCTIONS  <- file.path(VOI_DIR, "R", "voi_functions.R")
PROD_ENGINE    <- file.path(VOI_REPO_ROOT, "R", "abcd_engine.R")      # read-only, validation only
VOI_OUT_ROOT   <- file.path(VOI_REPO_ROOT, "output", "voi")           # every VOI file is written below here

# ------------------------------------------------------------------------------
# 1. RUN SETTINGS
# ------------------------------------------------------------------------------
VOI_COUNTRIES <- c("Indonesia", "Uganda")
VOI_RUN_ID    <- "sap4_n500_2026-11-10"   # output/voi/<VOI_RUN_ID>/<country>/
N_DRAWS       <- 500L                     # paired draws per country (each run under ref AND sap).
                                          # Do NOT raise automatically; raise here later to extend.
N_PILOT       <- 20L                      # pilot draws (draw ids 1..N_PILOT), run first, separately
VOI_SEED      <- 20261110L                # base seed: draw i uses the i-th L'Ecuyer-CMRG stream
VOI_CHUNK     <- 25L                      # draws per checkpoint file
VOI_MAX_CORES <- 12L                      # PSOCK workers (each loads the VOI engine in its own process)
VOI_N_BOOT    <- 1000L                    # bootstrap replicates (complete paired draw rows)
VOI_BOOT_SEED <- 20261111L
VOI_PRIOR_DIAG_N <- 200000L               # large-sample prior Monte Carlo for prior diagnostics only

# Validation tolerances (relative, |a - b| / max(1, |b|)), documented in README
VOI_TOL <- list(
  engine_identity   = 0,       # VOI engine (shared or split at the point values) vs production engine: exact
  saved_outputs     = 1e-9,    # production-engine re-run vs saved out_model annual sums (summation order only)
  economics         = 1e-9,    # recomputed economics vs saved rhd_economic_results.rds
  no_effect_deaths  = 1e-9)    # M_AB = M_BC = 1: deaths averted must vanish (relative to total deaths)

# ------------------------------------------------------------------------------
# 2. PRIOR SPECIFICATION (FINAL FOR THIS SEPARATE ANALYSIS; PROVISIONAL, NOT AGREED)
#    Starting point: dev/RHD_EVPPI_parameter_register.xlsx (SAP_effect_EVPPI S1-S3, S8, S9),
#    consulted as a reference only (never read at run time). Mapping resolved against the
#    current engine (see README "Evidence-to-parameter mapping"). Production point values
#    (00 block H/I): M_AB = M_BC = 0.45 (shared M), protection age 39, partial 0.5.
#
#    M = FULL-ADHERENCE multiplier on the treated progression HAZARD; the engine applies
#    h_BPG = h x [1 - a (1 - M)] with programme adherence a = good + rel x partial modelled
#    SEPARATELY (reference a = 0.6 + 0.2 rel; scale-up target 0.85 + 0.1 rel). Efficacy
#    extrapolation (M) and adherence response (rel, a) are therefore distinct inputs.
#
#    Two-regime mixture with ONE shared regime indicator for M_AB and M_BC (fixed weight):
#    Regime T ("trial efficacy extrapolates", weight w_trial):
#      M_AB ~ LN, 95% 0.03-0.29 (median 0.093): Bray 2024 Cochrane rapid review RR for early/
#             mild RHD, i.e. GOAL pooled latent-RHD progression 3/399 vs 33/400 (constant-
#             hazard HR 0.088). POOLED endpoint (borderline ~A and mild definite ~B; any
#             echocardiographic progression), NOT an A->B estimate as such: it is used for
#             A->B because the B-subgroup effect is weaker (HR 0.20, below) and the pooled
#             interval contains the borderline-subgroup value implied by the GOAL risks
#             recorded in 00 block H (0.075 vs 0.003 -> HR 0.039; unverified).
#      M_BC ~ LN, 95% 0.068-0.59 (median 0.20): GOAL mild-definite (~B) subgroup 2-year risks
#             0.119 vs 0.025 recorded in 00 block H (HR 0.200; subgroup CI not recorded:
#             log-sd chosen for small numbers; verify against the NEJM supplement).
#      log-scale correlation rho_trial = 0.5 within the regime (same trial, same children).
#    Regime R ("benefit only through prevented RF recurrences", weight 1 - w_trial):
#      An RF-recurrence RR is NOT a progression-hazard multiplier. Bridge: if a share s of
#      the untreated progression hazard of early disease is driven by recurrent ARF and BPG
#      multiplies the recurrence hazard by RR_rec, then h_BPG = h x [(1 - s) + s RR_rec],
#      i.e. M = 1 - s (1 - RR_rec), the SAME for A->B and B->C (one mechanism) ->
#      M_AB = M_BC in this regime.
#      RR_rec ~ LN, 95% 0.22-0.69 (Bray 2024: RF recurrence, 6 RCTs; RR ~ HR at low risk).
#      s ~ Beta(8, 2) (mean 0.80, 95% 0.52-0.97): NEW judgement (no data identified);
#      s = 1 reproduces the workbook's recurrence component (M = RR_rec).
#    Truncation: M in (0, 1]; a draw with any M > 1 is redrawn (rejections counted).
#    sap_protection_age_max: categorical 20 / 39 / 100 (100 = top model age = no cap) with
#      P = 0.25 / 0.50 / 0.25 (workbook S8); independent of M.
#    partial_adherence_relative_protection ~ Beta(4.537, 4.537) (mean 0.5, 95% 0.20-0.80;
#      workbook S9, de Dassel 2018); independent of M.
#    Mixture weight w_trial = 0.5 FIXED within the base prior scenario; 0.25 / 0.75 are
#    explored only by importance re-weighting of the same draws (joint quantities only).
#    The same draws are used for both countries (the SAP efficacy prior is not country-specific).
# ------------------------------------------------------------------------------
ln_from_ci <- function(lo, hi, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  c(meanlog = (log(lo) + log(hi)) / 2, sdlog = (log(hi) - log(lo)) / (2 * z))
}
VOI_PRIOR <- list(
  scenario  = "base",
  status    = "PROVISIONAL: proposed for this exploratory VOI analysis; NOT expert-agreed (no SAP-effect distribution has been agreed)",
  w_trial   = 0.5,
  M_AB_trial = c(ln_from_ci(0.03, 0.29), ci_lo = 0.03, ci_hi = 0.29),
  M_BC_trial = c(ln_from_ci(0.068, 0.59), ci_lo = 0.068, ci_hi = 0.59),
  rho_log_trial = 0.5,
  RR_recurrence = c(ln_from_ci(0.22, 0.69), ci_lo = 0.22, ci_hi = 0.69),
  share_recurrence = c(shape1 = 8, shape2 = 2),
  M_upper = 1,
  max_rejections = 1000L,
  protection_age = list(values = c(20L, 39L, 100L), probs = c(0.25, 0.50, 0.25)),
  partial_protection = c(shape1 = 4.537, shape2 = 4.537),
  sources = c(
    M_AB_trial = "Bray et al. 2024 (Cochrane rapid review) early/mild RHD RR 0.09 (0.03-0.29) = GOAL pooled latent-RHD progression (Beaton 2022); pooled endpoint used for A->B (see README)",
    M_BC_trial = "GOAL mild-definite subgroup 2-y risks 0.119 vs 0.025 as recorded in 00_run_all.R block H (HR 0.200); interval from workbook S2 (log-sd 0.551)",
    rho_log_trial = "workbook S2 (within-regime log correlation 0.5)",
    RR_recurrence = "Bray et al. 2024 RF recurrence RR 0.39 (0.22-0.69), 6 RCTs",
    share_recurrence = "NEW provisional judgement for this analysis (not in the workbook); requires elicitation",
    w_trial = "workbook S3 (0.5; scenarios 0.25 / 0.75)",
    protection_age = "workbook S8",
    partial_protection = "workbook S9; de Dassel et al. 2018"))
VOI_PRIOR_REWEIGHT_W <- c(0.25, 0.50, 0.75)      # exploratory prior-scenario re-weighting (joint quantities)

# Production point values (validated against baseline_state in 02; listed for the record)
VOI_POINT <- list(M_AB = 0.45, M_BC = 0.45, sap_protection_age_max = 39L, partial_adherence_relative_protection = 0.5)

# ------------------------------------------------------------------------------
# 3. VALUATION THRESHOLDS
#    VSL family : INB(k) = k x sum_t df_t VSL_t x deaths averted_t - sum_t df_t dC_t,
#                 VSL_t = 08 reference-case schedule (constant 2025 US$); k = 1 is the
#                 production valuation.
#    DALY family: INB(lambda) = lambda x sum_t df_t deaths averted_t x dalys_per_death -
#                 sum_t df_t dC_t; dalys_per_death = 30 is the existing deaths-to-DALYs
#                 PROXY (08 cost per DALY), not a DALY model. lambda in constant 2025 US$
#                 per DALY, expressed as multiples of a GNI-per-capita reference derived
#                 from 08's own income path (2026, constant 2025 US$). Illustrative, not an
#                 endorsed willingness-to-pay.
# ------------------------------------------------------------------------------
VOI_VSL_MULT_REF  <- c(0.05, 0.10, 1.00)
VOI_VSL_MULT_GRID <- sort(unique(round(c(exp(seq(log(0.005), log(2), length.out = 71)), VOI_VSL_MULT_REF), 8)))
VOI_DALY_GNI_REF  <- c(0.5, 1.0)
VOI_DALY_GNI_GRID <- sort(unique(round(c(seq(0, 3, by = 0.05), VOI_DALY_GNI_REF), 8)))

# ------------------------------------------------------------------------------
# 4. STAGES + configuration object
# ------------------------------------------------------------------------------
if (!exists("VOI_STAGES")) VOI_STAGES <- c("sample", "validate", "pilot", "psa", "analyse")

VOI_CFG <- list(
  root = VOI_REPO_ROOT, voi_dir = VOI_DIR, engine_file = VOI_ENGINE, functions_file = VOI_FUNCTIONS,
  prod_engine_file = PROD_ENGINE, out_root = VOI_OUT_ROOT,
  countries = VOI_COUNTRIES, run_id = VOI_RUN_ID, n_draws = N_DRAWS, n_pilot = N_PILOT,
  seed = VOI_SEED, chunk = VOI_CHUNK, max_cores = VOI_MAX_CORES, n_boot = VOI_N_BOOT, boot_seed = VOI_BOOT_SEED,
  prior_diag_n = VOI_PRIOR_DIAG_N, tol = VOI_TOL, prior = VOI_PRIOR, prior_reweight_w = VOI_PRIOR_REWEIGHT_W,
  point = VOI_POINT,
  vsl_mult_ref = VOI_VSL_MULT_REF, vsl_mult_grid = VOI_VSL_MULT_GRID,
  daly_gni_ref = VOI_DALY_GNI_REF, daly_gni_grid = VOI_DALY_GNI_GRID,
  scenarios = c("ref", "sap"),
  params = c("M_AB", "M_BC", "sap_protection_age_max", "partial_adherence_relative_protection"),
  engine_spec_expected = "abcd_first_diagnosis_v1_adult_prog_mult_sap55_dmort_age")

source(VOI_FUNCTIONS)          # shared helpers (paths, logging, runner, economics); defines no engine code

if (isTRUE(get0("VOI_CONFIG_ONLY", ifnotfound = FALSE))) {
  message("VOI configuration loaded (VOI_CONFIG_ONLY = TRUE); no stage run.")
} else {
  t_start <- Sys.time()
  integrity0 <- voi_production_fingerprint(VOI_CFG)        # production files before the run
  message(sprintf("== SAP VOI run '%s': %s | N_DRAWS = %d | stages: %s ==", VOI_RUN_ID,
                  paste(VOI_COUNTRIES, collapse = ", "), N_DRAWS, paste(VOI_STAGES, collapse = " -> ")))
  stage_files <- c(sample = "01_sample_sap_parameters.R", validate = "02_validate_voi.R",
                   pilot = "03_run_psa.R", psa = "03_run_psa.R", analyse = "04_analyse_and_plot.R")
  for (stg in VOI_STAGES) {
    VOI_CURRENT_STAGE <- stg
    message("\n---- stage: ", stg, " ----")
    source(file.path(VOI_DIR, stage_files[[stg]]))
  }
  rm(VOI_CURRENT_STAGE)
  integrity1 <- voi_production_fingerprint(VOI_CFG)        # production files after the run
  voi_check_integrity(VOI_CFG, integrity0, integrity1)
  message(sprintf("== SAP VOI run '%s' finished in %.1f min ==", VOI_RUN_ID, as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
