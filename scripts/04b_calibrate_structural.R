#===============================================================================
# 04b_calibrate_structural.R
#-------------------------------------------------------------------------------
# LAYER-2 STRUCTURAL A/B/C/D CALIBRATION — JOINT INCIDENCE + NATURAL HISTORY,
# PREVALENCE-LED, ON THE PRODUCTION ENGINE
#
# Runs the ACTUAL production engine (R/abcd_engine.R: abcd_one_cycle + the shared
# abcd_project / abcd_burn_in / abcd_age_shift wrappers, the same ones 06 uses)
# inside the objective, and JOINTLY estimates
#   * low-dimensional INCIDENCE corrections (one log-multiplier per broad age band
#     x sex on the GBD incidence pattern; Layer-1 alpha are STARTING VALUES only),
#   * the stage TRANSITIONS, and
#   * the RHD MORTALITY level under the selected MORTALITY STRUCTURE:
#       "D_only"   (PRIMARY) p_rhd_death_A = _B = _C = 0 EXACTLY; one stage-D level;
#                  residual multiplier FIXED at 1 (documented fixed field);
#       "all_stage" (LEGACY) severity-ordered m_A <= m_B < m_C < m_D (nested logistic
#                  ratios) + bounded residual multiplier (bound-hit = FAILURE),
# against GBD RHD PREVALENCE (the principal target, incl. its age pattern) and
# RHD-SPECIFIC DEATHS, with incidence REPORTED as a validation output (direct
# incidence weight 0 in the primary objective; a weak non-zero weight is a labelled
# sensitivity). Reference care (treatment coverage x adherence, surgery coverage)
# and the selected SAP mechanism are applied INSIDE the objective exactly as in
# production, so the fitted parameters are untreated natural history and nothing
# is applied twice.
#
# OBJECTIVE (transparent; every term exported):
#   L = wP * MSLE(prevalence) + wD * MSLE(RHD deaths) + wI * MSLE(incidence)
#       + LAMBDA * [ prior(mortality) + prior(transitions) + prior(alpha) (+ prior(resid)) ]
#       + LAMBDA_SMOOTH * smooth(alpha)  + wS * stage-target error (if usable targets)
#   MSLE = weighted mean over MASKED cells (both sexes pooled) of
#          (log(model + eps) - log(target + eps))^2, weights = age x year x optional IV.
#   Masked cells: missing observations (NA in GBD) and counts < TARGET_MIN_COUNT.
#   wP/wD = STRUCT_PREV_DEATH_RATIO (2.5 primary), wP + wD fixed (controlled scale).
#
# TIME CONVENTION (shared with 06 via abcd_project): start stock -> within-year
# flows -> end stock -> ageing. Prevalence observation = STRUCT_PREV_OBS ("mid" =
# mean of start and end stock, a mid-year point-prevalence approximation to GBD's
# annual prevalence); deaths/incident cases = within-year flows.
#
# INITIALISATION: zero-start burn-in under first-window-year inputs to the
# transition-implied equilibrium; convergence is CHECKED (and a 2x burn-in run) and
# reported, not asserted. After the fit, the end-of-window state is FORWARD-
# PROJECTED through the bridge years to the ANALYSIS start (reference care,
# calibrated incidence x trend anchored at CAL_YEAR_END, held background
# mortality) and persisted as the 06 seed, with a reconciliation table vs GBD.
#
# STAGE TARGETS: a supplied, VALIDATED stage-target table (04) is READ, mapped
# (stage / ages / sex / year / denominator), scored and INCLUDED in the loss when
# usable rows exist; metadata distinguish detected / usable / USED, and the status
# says "with_echo_targets" ONLY if rows actually entered the objective.
#
# VALIDATION + FINALISATION: separate structural / numerical / epidemiological-fit
# status fields. A candidate failing a STRUCTURAL check (NA/negative flows, exact
# zeros violated, residual at bound in the legacy branch, stock > population) is
# NOT written as best_parameters: the bundle status becomes FAILED_*, the failed
# candidate and reason are saved separately, and the script stops. The previous
# valid bundle (if any) is kept as calibrated_rhd_parameters_last_valid.rds.
#
# OUTPUT: extends data/<COUNTRY>/calibrated_rhd_parameters.rds $stage_calibration
#   (best_parameters incl. the FINAL incidence, initial_state for 06, accepted sets,
#    objective decomposition, priors, bounds, identifiability, statuses, attribution)
#   + calibration_structural_*.csv diagnostics. 05 consumes $stage_calibration.
#
# Source AFTER 03 (disease inputs) and 04 (Layer-1 bundle incl. $layer2_targets).
#===============================================================================

library(data.table)
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

if (!exists("wd_raw"))  wd_raw  <- paste0(here::here("data-raw"), "/")
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

# shared engine (single source of truth) — always (re)source so the wrappers exist
eng <- if (exists("wd")) paste0(wd, "R/abcd_engine.R") else here::here("R", "abcd_engine.R")
source(eng)

#===============================================================================
# 0. SETTINGS (honour 00_run_all.R; else documented defaults)
#===============================================================================
LOCATION  <- getp("LOCATION", "Indonesia")
CAL_YEAR_START <- as.integer(getp("CAL_YEAR_START", 2000L))
CAL_YEAR_END   <- as.integer(getp("CAL_YEAR_END",   2019L))
MATCH_FROM <- as.integer(getp("STRUCT_MATCH_FROM", 2010L))
SEED <- as.integer(getp("SEED", 42L))
RUN_STRUCT <- isTRUE(getp("run_structural_calibration", TRUE))
CALIB_PDC  <- isTRUE(getp("calibrate_p_D_to_C", FALSE))
MORT_STRUCT <- match.arg(getp("mortality_structure", "D_only"), c("D_only", "all_stage"))

## priors (sd) + weights
SIGMA_MORT  <- getp("SIGMA_MORT",  0.35)
SIGMA_TRANS <- getp("SIGMA_TRANS", 1.50)
SIGMA_RESID <- getp("SIGMA_RESID", 0.25)
LAMBDA_PRIOR_STRUCT <- getp("LAMBDA_PRIOR_STRUCT", 0.05)   # structural priors (transitions, alpha, age multipliers)
LAMBDA_PRIOR_MORT   <- getp("LAMBDA_PRIOR_MORT",   1.0)    # mortality (+ legacy residual) prior keeps full strength
PD_RATIO <- getp("STRUCT_PREV_DEATH_RATIO", 2.5)
W_TOTAL  <- getp("STRUCT_DATA_WEIGHT_TOTAL", 3.5)
W_PREV_STRUCT  <- getp("W_PREV_STRUCT",  W_TOTAL * PD_RATIO / (1 + PD_RATIO))
W_DEATH_STRUCT <- getp("W_DEATH_STRUCT", W_TOTAL / (1 + PD_RATIO))
W_INC_STRUCT   <- getp("W_INC_STRUCT", 0)            # 0 = incidence is a VALIDATION output
W_STAGE_STRUCT <- getp("W_STAGE_STRUCT", 1)          # weight on usable stage targets
SIGMA_ALPHA_STRUCT   <- getp("SIGMA_ALPHA_STRUCT", 1.0)
LAMBDA_SMOOTH_STRUCT <- getp("LAMBDA_SMOOTH_STRUCT", 0.5)
ALPHA_LO <- getp("STRUCT_ALPHA_LO", 0.2); ALPHA_HI <- getp("STRUCT_ALPHA_HI", 5.0)
USE_IV_W <- isTRUE(getp("USE_IV_WEIGHTS_STRUCT", FALSE)); IV_CAP <- getp("STRUCT_IV_WEIGHT_CAP", 5)
AGE_W    <- getp("STRUCT_AGE_WEIGHTS", NULL)
PREV_OBS <- match.arg(getp("STRUCT_PREV_OBS", "mid"), c("mid", "start", "end"))
PREV_FIRST <- isTRUE(getp("STRUCT_PREV_FIRST", TRUE))
FIT_TOL_PREV  <- getp("STRUCT_FIT_TOL_PREV", 0.15)
FIT_TOL_DEATH <- getp("STRUCT_FIT_TOL_DEATH", 0.25)
TARGET_MIN_COUNT <- getp("TARGET_MIN_COUNT", 1)
EPS_LOG <- getp("STRUCT_EPS_COUNT", 1)                 # count-scale epsilon: one case (log(x + 1))

## bounds
MD_LO <- getp("MD_LO", 0.03); MD_HI <- getp("MD_HI", 0.25)
RESID_LO <- getp("RESID_MULT_LO", 0.5); RESID_HI <- getp("RESID_MULT_HI", 2.0)
TRANS_BOUNDS <- getp("TRANS_BOUNDS", list(
  p_A_to_no_rhd = c(0.005, 0.60), p_A_to_B = c(0.001, 0.25), p_B_to_A = c(0.002, 0.40),
  p_B_to_C = c(0.001, 0.20), p_C_to_B = c(0.001, 0.15), p_C_to_D = c(0.001, 0.30),
  p_D_to_C = c(0.000, 0.10)))
N_STARTS <- as.integer(getp("STRUCT_N_STARTS", 2L))
NEAROPT_DELTA <- getp("STRUCT_NEAROPT_DELTA", 0.05)
BURN_YEARS <- as.integer(getp("STRUCT_BURN_YEARS", 80L))
BRIDGE <- isTRUE(getp("STRUCT_BRIDGE_TO_ANALYSIS_START", TRUE))
ANALYSIS_YEAR_START <- as.integer(getp("ANALYSIS_YEAR_START", 2026L))

## files
BUNDLE_FILE  <- paste0(wd_data, "calibrated_rhd_parameters.rds")
LASTVALID_FILE <- paste0(wd_data, "calibrated_rhd_parameters_last_valid.rds")
DISEASE_FILE <- paste0(wd_data, "disease_model_inputs.rds")
POP_OBS_FILE <- paste0(wd_data, "pop_observed.rds")
POP_PROJ_FILE <- paste0(wd_data, "pop_projection.rds")
for (f in c(BUNDLE_FILE, DISEASE_FILE))
  if (!file.exists(f)) stop("04b: missing input ", f, " (run 03 and 04 first).", call. = FALSE)

calib <- readRDS(BUNDLE_FILE)
dmi   <- readRDS(DISEASE_FILE)
if (is.null(calib$layer2_targets) || is.null(calib$layer2_inputs))
  stop("04b: the bundle lacks $layer2_targets / $layer2_inputs -- re-run 04_calibration_random_tp.R.", call. = FALSE)

# spec consistency: 03's inputs were built under the SAME mortality structure / SAP mechanism
if (!identical(dmi$meta$mortality_structure, MORT_STRUCT))
  stop(sprintf("04b: mortality_structure '%s' but disease_model_inputs were built under '%s' -- re-run 03.",
               MORT_STRUCT, dmi$meta$mortality_structure), call. = FALSE)
SAP_MECH <- dmi$effects$sap_mechanism
if (!identical(SAP_MECH, getp("sap_mechanism", SAP_MECH)))
  stop("04b: sap_mechanism in 00 differs from disease_model_inputs -- re-run 03.", call. = FALSE)

tr_prior <- dmi$transitions                              # transition prior centers (scalars)
pd_prior <- dmi$p_rhd_death                              # mortality prior centers c(A,B,C,D)
dfrac_prior <- dmi$meta$rhd_d_fraction
eff_cal <- dmi$effects; surg_cal <- dmi$surgery
REM_MULT <- if (is.null(dmi$meta$remission_adult_mult)) 1 else dmi$meta$remission_adult_mult
REM_CUT  <- if (is.null(dmi$meta$remission_age_cutoff)) 20L else dmi$meta$remission_age_cutoff
# parsimonious AGE STRUCTURE (00 block G): calibrate an adult (age >= REM_CUT)
# multiplier on the progression hazards and one on remission; else fixed.
CALIB_AGE <- isTRUE(getp("STRUCT_CALIB_AGE_STRUCTURE", TRUE))
SIGMA_AGE_MULT <- getp("SIGMA_AGE_MULT", 1.0)
AGE_MULT_LO <- getp("AGE_MULT_LO", 0.05); AGE_MULT_HI <- getp("AGE_MULT_HI", 3.0)
# REFERENCE care applied in the objective (untreated parameters + reference care in
# the engine, exactly as 06's reference arm): treatment coverage x adherence; surgery.
C0_TREAT <- dmi$coverage$treatment_ref_baseline
C0_ADH   <- if (is.null(dmi$coverage$adherence_ref_baseline)) 1 else dmi$coverage$adherence_ref_baseline
C0_SURG  <- dmi$surgery$coverage_ref_baseline
# cascade capping as in 05 (effective treatment <= diagnosis <= screening)
C0_TREAT_EFF <- min(dmi$coverage$screen_ref_baseline, dmi$coverage$diagnosis_ref_baseline, C0_TREAT)

if (!RUN_STRUCT) {
  message("04b: run_structural_calibration = FALSE -> leaving stage params at 03 priors.")
  calib$stage_calibration$status <- "structural_disabled_using_priors"
  saveRDS(calib, BUNDLE_FILE)
  message("04b: wrote passthrough status; 05 will use 03 stage params.")
} else {

message(sprintf("── 04b structural A/B/C/D calibration : %s | match %d-%d | mortality=%s | SAP=%s | engine=R/abcd_engine.R ──",
                LOCATION, MATCH_FROM, CAL_YEAR_END, MORT_STRUCT, SAP_MECH))

#===============================================================================
# 1. GRIDS, INPUT ARRAYS AND TARGET ARRAYS (from 04's $layer2_inputs / $layer2_targets)
#===============================================================================
age_to_gbd_group <- function(a) fcase(
  a < 1,"<1 year", a < 2,"12-23 months", a < 5,"2-4 years", a < 10,"5-9 years",
  a < 15,"10-14 years", a < 20,"15-19 years", a < 25,"20-24 years", a < 30,"25-29 years",
  a < 35,"30-34 years", a < 40,"35-39 years", a < 45,"40-44 years", a < 50,"45-49 years",
  a < 55,"50-54 years", a < 60,"55-59 years", a < 65,"60-64 years", a < 70,"65-69 years",
  a < 75,"70-74 years", a < 80,"75-79 years", a < 85,"80-84 years", a < 90,"85-89 years",
  a < 95,"90-94 years", default = "95+ years")

L2I <- as.data.table(calib$layer2_inputs)[location == LOCATION & year >= CAL_YEAR_START & year <= CAL_YEAR_END]
L2T <- as.data.table(calib$layer2_targets)[location == LOCATION & year >= CAL_YEAR_START & year <= CAL_YEAR_END]
ages <- sort(unique(L2I$age)); n_age <- length(ages)
if (!identical(as.integer(ages), as.integer(dmi$meta$AGES)))
  stop("04b: calibration age grid differs from the production grid (dmi$meta$AGES). Set AGE_HI = 100 in 04.", call. = FALSE)
SEXES <- intersect(c("Female", "Male"), unique(L2I$sex)); n_sex <- length(SEXES)
yrs  <- CAL_YEAR_START:CAL_YEAR_END; n_yr <- length(yrs)
grp_lab_age <- age_to_gbd_group(ages)
grp_present <- unique(grp_lab_age)
grp_int     <- match(grp_lab_age, grp_present); G <- length(grp_present)
match_cols  <- yrs >= MATCH_FROM
# the 95+ GBD group is compared with the model AGGREGATED over ages 95..100
stopifnot(sum(grp_lab_age == "95+ years") == sum(ages >= 95))
band_edges <- calib$meta$band_edges; band_labels <- calib$meta$bands
n_band <- length(band_edges); band_of <- function(a) findInterval(a, band_edges)
band_age <- band_of(ages)

# input arrays [age x sex x year]
arr0 <- function() array(0, dim = c(n_age, n_sex, n_yr), dimnames = list(ages, SEXES, yrs))
to_arr <- function(col) { a <- arr0()
  a[cbind(match(L2I$age, ages), match(L2I$sex, SEXES), match(L2I$year, yrs))] <- L2I[[col]]; a }
POP  <- to_arr("pop_wpp"); IRG <- to_arr("ir_gbd_total_pop"); PREVG <- to_arr("prev_gbd"); OTH <- to_arr("bg_mx")
if (any(!is.finite(POP)) || any(POP < 0) || any(!is.finite(IRG)) || any(!is.finite(OTH)))
  stop("04b: non-finite/negative Layer-2 inputs.", call. = FALSE)

# target arrays [G x sex x year]: value used (NA = missing), mask, weights
tarr <- function(meas, col) { a <- array(NA_real_, dim = c(G, n_sex, n_yr), dimnames = list(grp_present, SEXES, yrs))
  d <- L2T[measure == meas]
  a[cbind(match(d$age_group, grp_present), match(d$sex, SEXES), match(d$year, yrs))] <- d[[col]]; a }
TGT <- list(prev = tarr("Prevalence", "value_used"), death = tarr("Deaths", "value_used"),
            inc = tarr("Incidence", "value_used"))
TGT_NUM <- list(prev = tarr("Prevalence", "number"), death = tarr("Deaths", "number"), inc = tarr("Incidence", "number"))
POP_G <- tarr("Prevalence", "pop_wpp")                   # WPP population per target cell (rates)
year_w <- as.numeric(match_cols)                        # years before MATCH_FROM carry weight 0 (held-out)
age_w  <- rep(1, G)
if (!is.null(AGE_W)) {                                   # documented age/evidence weights by band label
  grp_lo_age <- vapply(grp_present, function(g) min(ages[grp_lab_age == g]), numeric(1))
  bl <- band_labels[band_of(grp_lo_age)]
  age_w <- ifelse(bl %in% names(AGE_W), as.numeric(AGE_W[bl]), 1)
}
mk_w <- function(meas) {
  w <- array(1, dim = c(G, n_sex, n_yr)) * age_w            # recycles over the first dim
  w <- sweep(w, 3, year_w, `*`)
  if (USE_IV_W) {                                           # inverse-variance (log scale), floor + cap
    lv <- tarr(meas, "logvar"); iv <- 1 / lv; iv[!is.finite(iv)] <- NA
    if (any(is.finite(iv))) { iv <- iv / mean(iv, na.rm = TRUE); iv[is.na(iv)] <- 1; iv <- pmin(iv, IV_CAP); w <- w * iv }
  }
  w
}
WT  <- list(prev = mk_w("Prevalence"), death = mk_w("Deaths"), inc = mk_w("Incidence"))
MSK <- lapply(TGT, function(a) !is.na(a) & a >= TARGET_MIN_COUNT)
MSK <- lapply(MSK, function(m) { m[is.na(m)] <- FALSE; m })
# explicit death-target age floor (00 block G2; 0 = none). Masked cells stay in the
# fit tables as validation (in_loss = FALSE).
DEATH_AGE_MIN <- as.integer(getp("STRUCT_DEATH_TARGET_AGE_MIN", 0L))
grp_lo_age_all <- vapply(grp_present, function(g) min(ages[grp_lab_age == g]), numeric(1))
n_death_age_masked <- 0L
if (DEATH_AGE_MIN > 0) { drop_g <- grp_lo_age_all < DEATH_AGE_MIN
  n_death_age_masked <- sum(MSK$death[drop_g, , ] & WT$death[drop_g, , ] > 0); MSK$death[drop_g, , ] <- FALSE }
message(sprintf("  targets (%s basis): prevalence cells in loss %d | deaths %d (age floor %d: %d cells masked) | incidence (validation) %d | missing masked: %d",
                calib$meta$layer2_target_basis, sum(MSK$prev & WT$prev > 0), sum(MSK$death & WT$death > 0),
                DEATH_AGE_MIN, n_death_age_masked, sum(MSK$inc & WT$inc > 0), sum(is.na(TGT$prev)) + sum(is.na(TGT$death))))

# stage targets (validated by 04; USED only here)
STG <- calib$stage_calibration$stage_targets
STG <- if (!is.null(STG)) as.data.table(STG)[reason == "usable"] else NULL
n_stage_used <- if (is.null(STG)) 0L else nrow(STG)

#===============================================================================
# 2. PARAMETER PACK / UNPACK (branch by mortality structure)
#===============================================================================
abc_raw <- c(A = 0.565, B = 0.272, C = 0.162); abc <- abc_raw / sum(abc_raw)  # Cannon split (legacy stage_split_of)
trans_names <- c("p_A_to_no_rhd","p_A_to_B","p_B_to_A","p_B_to_C","p_C_to_B","p_C_to_D")
if (CALIB_PDC) trans_names <- c(trans_names, "p_D_to_C")
n_tr <- length(trans_names); n_alpha <- n_band * n_sex
D_ONLY <- MORT_STRUCT == "D_only"
# layout: [logit(trans) x n_tr, log(mD), {thetaC, thetaB, thetaA, log(resid)} (legacy only),
#          {log(adult_prog_mult), log(adult_remission_mult)} (if CALIB_AGE), log(alpha) x n_alpha]
idx <- list(tr = seq_len(n_tr), mD = n_tr + 1L)
k <- n_tr + 1L
if (!D_ONLY) { idx$rC <- k + 1L; idx$rB <- k + 2L; idx$rA <- k + 3L; idx$res <- k + 4L; k <- k + 4L }
if (CALIB_AGE) { idx$agp <- k + 1L; idx$agr <- k + 2L; k <- k + 2L }
idx$alpha <- k + seq_len(n_alpha)
npar <- k + n_alpha
alpha_names <- as.vector(outer(band_labels, SEXES, function(b, s) paste0("alpha[", b, "|", s, "]")))
param_labels <- c(trans_names, "log_mD", if (!D_ONLY) c("ratio_C", "ratio_B", "ratio_A", "log_resid"),
                  if (CALIB_AGE) c("log_adult_prog_mult", "log_adult_remission_mult"), alpha_names)

tb <- function(nm) { b <- TRANS_BOUNDS[[nm]]; if (is.null(b)) c(1e-4, 0.6) else b }
clamp_logit <- function(p, nm) qlogis(pmin(pmax(p, max(tb(nm)[1], 1e-6)), min(tb(nm)[2], 1 - 1e-6)))

unpack <- function(theta) {
  tr <- lapply(seq_len(n_tr), function(i) plogis(theta[idx$tr[i]])); names(tr) <- trans_names
  if (!CALIB_PDC) tr$p_D_to_C <- tr_prior$p_D_to_C
  mD <- exp(theta[idx$mD])
  if (D_ONLY) {
    pd <- c(A = 0, B = 0, C = 0, D = mD); ratios <- c(A = 0, B = 0, C = 0); resid <- 1
  } else {
    rC <- plogis(theta[idx$rC]); rB <- rC * plogis(theta[idx$rB]); rA <- rB * plogis(theta[idx$rA])
    pd <- c(A = rA * mD, B = rB * mD, C = rC * mD, D = mD); ratios <- c(A = rA, B = rB, C = rC)
    resid <- exp(theta[idx$res])
  }
  alpha <- matrix(exp(theta[idx$alpha]), n_band, n_sex, dimnames = list(band_labels, SEXES))
  age <- if (CALIB_AGE) c(prog = exp(theta[idx$agp]), rem = exp(theta[idx$agr])) else c(prog = 1, rem = REM_MULT)
  list(tr = tr, pd = pd, mD = mD, ratios = ratios, resid = resid, alpha = alpha, age = age)
}
# engine-side transition list: the scalar (child / reference-age) probabilities,
# with the explicit ADULT (age >= REM_CUT) multipliers applied on the HAZARD scale
# to the progression transitions (A->B, B->C, C->D) and to A->No-RHD remission.
age_mult_hazard <- function(p, m) 1 - (1 - p)^m
tr_engine <- function(tr, age = c(prog = 1, rem = REM_MULT)) {
  adult <- ages >= REM_CUT
  mk <- function(p, m) { M <- matrix(p, n_age, n_sex, dimnames = list(ages, SEXES))
    M[adult, ] <- age_mult_hazard(p, m); M }
  if (age[["rem"]]  != 1) tr$p_A_to_no_rhd <- mk(tr$p_A_to_no_rhd, age[["rem"]])
  if (age[["prog"]] != 1) { tr$p_A_to_B <- mk(tr$p_A_to_B, age[["prog"]]); tr$p_B_to_C <- mk(tr$p_B_to_C, age[["prog"]])
                            tr$p_C_to_D <- mk(tr$p_C_to_D, age[["prog"]]) }
  tr
}
# incidence per SUSCEPTIBLE person: GBD per-total-population rate x alpha(band, sex)
# / (1 - GBD prevalence), applied ONCE here (05 consumes the final values)
ir_engine_arr <- function(alpha, irg = IRG, prevg = PREVG) {
  a <- irg
  for (s in seq_len(n_sex)) a[, s, ] <- irg[, s, ] * alpha[band_age, s] / pmax(1 - prevg[, s, ], 0.5)
  pmin(a, 0.99)
}

# bounds (unconstrained scale)
lo <- numeric(npar); hi <- numeric(npar)
for (i in seq_len(n_tr)) { b <- tb(trans_names[i]); lo[idx$tr[i]] <- qlogis(max(b[1], 1e-6)); hi[idx$tr[i]] <- qlogis(min(b[2], 1 - 1e-6)) }
lo[idx$mD] <- log(MD_LO); hi[idx$mD] <- log(MD_HI)
if (!D_ONLY) { lo[c(idx$rC, idx$rB, idx$rA)] <- -8; hi[c(idx$rC, idx$rB, idx$rA)] <- 8
               lo[idx$res] <- log(RESID_LO); hi[idx$res] <- log(RESID_HI) }
lo[idx$alpha] <- log(ALPHA_LO); hi[idx$alpha] <- log(ALPHA_HI)
if (CALIB_AGE) { lo[c(idx$agp, idx$agr)] <- log(AGE_MULT_LO); hi[c(idx$agp, idx$agr)] <- log(AGE_MULT_HI) }

# prior-centre start (theta0): transitions at 03 centres; mortality at the prior
# centre (D-only: the provisional col-K level; legacy: recover ratios); alpha at 1.
theta0 <- numeric(npar)
for (i in seq_len(n_tr)) theta0[idx$tr[i]] <- clamp_logit(tr_prior[[trans_names[i]]], trans_names[i])
mD0 <- pd_prior[["D"]]
if (mD0 < MD_LO || mD0 > MD_HI)
  stop(sprintf("Stage-D mortality prior centre %.4f lies outside the bounds [%.3f, %.3f]; adjust MD_LO/MD_HI deliberately (documented) rather than truncating.", mD0, MD_LO, MD_HI), call. = FALSE)
theta0[idx$mD] <- log(mD0)
if (!D_ONLY) {
  rC0 <- pmin(pmax(pd_prior[["C"]]/mD0, 1e-3), 1-1e-3)
  rB0 <- pmin(pmax(pd_prior[["B"]]/pd_prior[["C"]], 1e-3), 1-1e-3)
  rA0 <- pmin(pmax(pd_prior[["A"]]/pd_prior[["B"]], 1e-3), 1-1e-3)
  theta0[idx$rC] <- qlogis(rC0); theta0[idx$rB] <- qlogis(rB0); theta0[idx$rA] <- qlogis(rA0); theta0[idx$res] <- 0
}
theta0[idx$alpha] <- 0
if (CALIB_AGE) { theta0[idx$agp] <- 0; theta0[idx$agr] <- log(max(REM_MULT, AGE_MULT_LO)) }
theta0 <- pmin(pmax(theta0, lo), hi)
# Layer-1 alpha as an alternative START (never a prior centre)
alpha_L1 <- matrix(1, n_band, n_sex, dimnames = list(band_labels, SEXES))
fac <- as.data.table(calib$factors)[location == LOCATION]
for (s in SEXES) { f <- fac[sex == s]; if (nrow(f)) alpha_L1[match(f$band, band_labels), s] <- f$alpha }
theta_L1 <- theta0; theta_L1[idx$alpha] <- pmin(pmax(log(as.vector(alpha_L1)), lo[idx$alpha]), hi[idx$alpha])

#===============================================================================
# 3. ENGINE PROJECTION (shared wrapper) + OBSERVATIONS + OBJECTIVE
#===============================================================================
cov_ref <- list(treatment = rep(C0_TREAT_EFF, n_yr), surgery = rep(C0_SURG, n_yr), adherence = rep(C0_ADH, n_yr))
cov_ref1 <- list(treatment = C0_TREAT_EFF, surgery = C0_SURG, adherence = C0_ADH)
agg_grp <- function(M) rowsum(M, grp_int, reorder = FALSE)   # [age x sex] -> [G x sex]

project_window <- function(p, burn_years = BURN_YEARS, keep_states = FALSE) {
  tr <- tr_engine(p$tr, p$age); pd <- p$pd * p$resid; IR <- ir_engine_arr(p$alpha)
  bi <- abcd_burn_in(POP[, , 1], IR[, , 1], OTH[, , 1], tr, pd, eff_cal, surg_cal, cov_ref1, burn_years)
  pr <- abcd_project(bi$state, POP, IR, OTH, tr, pd, eff_cal, surg_cal, cov_ref, yrs,
                     light = TRUE, prev_obs = PREV_OBS)
  prev <- array(0, dim = c(G, n_sex, n_yr)); death <- prev; inc <- prev; deathD <- prev
  stg <- lapply(c("A", "B", "C", "D"), function(s) vector("list", n_yr)); names(stg) <- c("A", "B", "C", "D")
  for (t in seq_len(n_yr)) {
    cyc <- pr$cycles[[t]]; st <- pr$start[[t]]
    prev[, , t]  <- agg_grp(pr$prev_obs[[t]])
    death[, , t] <- agg_grp(cyc$rhd_deaths); deathD[, , t] <- agg_grp(cyc$rhd_deaths_D)
    inc[, , t]   <- agg_grp(cyc$new_rhd_A)
    for (s in c("A", "B", "C", "D"))                 # stage observation on the same convention
      stg[[s]][[t]] <- switch(PREV_OBS, start = st[[s]], end = cyc[[s]], mid = 0.5 * (st[[s]] + cyc[[s]]))
  }
  out <- list(prev = prev, death = death, deathD = deathD, inc = inc, stages = stg,
              final_state = pr$final_state, burn = bi, pop_deficit = pr$pop_deficit_total)
  if (keep_states) out$start <- pr$start
  out
}

loss_measure <- function(model, tgt, msk, wt) {
  m <- msk & wt > 0
  if (!any(m)) return(0)
  r <- (log(model[m] + EPS_LOG) - log(tgt[m] + EPS_LOG))^2 * wt[m]
  sum(r) / sum(wt[m])
}
# stage-target loss on the USABLE rows (model observation built to the target definition)
stage_loss <- function(pw) {
  if (n_stage_used == 0L) return(0)
  num <- 0; den <- 0
  for (i in seq_len(nrow(STG))) {
    r <- STG[i]; t <- match(r$year, yrs); ai <- which(ages >= r$age_lo & ages <= r$age_hi)
    si <- if (r$sex == "Both") seq_len(n_sex) else match(r$sex, SEXES)
    cases <- sum(pw$stages[[r$stage]][[t]][ai, si])
    den_m <- if (r$target_type == "prevalence_fraction") sum(POP[ai, si, t]) else
      sum(vapply(c("A", "B", "C", "D"), function(s) sum(pw$stages[[s]][[t]][ai, si]), numeric(1)))
    obs <- if (den_m > 0) cases / den_m else 0
    num <- num + r$weight * (log(obs + 1e-6) - log(r$target_prevalence + 1e-6))^2; den <- den + r$weight
  }
  num / den
}
prior_terms <- function(p, theta) {
  pri_mort <- if (D_ONLY) (log(p$pd[["D"]]) - log(pd_prior[["D"]]))^2 / SIGMA_MORT^2 else
    sum((log(p$pd) - log(pd_prior[c("A","B","C","D")]))^2) / SIGMA_MORT^2
  pri_tr <- sum(vapply(trans_names, function(nm) (qlogis(p$tr[[nm]]) - clamp_logit(tr_prior[[nm]], nm))^2, numeric(1))) / SIGMA_TRANS^2
  la <- log(p$alpha)
  pri_alpha <- sum((la / SIGMA_ALPHA_STRUCT)^2)
  smo_alpha <- if (n_band >= 3) sum(apply(la, 2, function(v) sum(diff(v, differences = 2)^2))) else 0
  pri_res <- if (D_ONLY) 0 else (log(p$resid)^2) / SIGMA_RESID^2
  pri_age <- if (CALIB_AGE) ((log(p$age[["prog"]]))^2 + (log(p$age[["rem"]]) - log(max(REM_MULT, AGE_MULT_LO)))^2) / SIGMA_AGE_MULT^2 else 0
  c(prior_mortality = LAMBDA_PRIOR_MORT * pri_mort, prior_transitions = LAMBDA_PRIOR_STRUCT * pri_tr,
    prior_incidence_alpha = LAMBDA_PRIOR_STRUCT * pri_alpha, smooth_incidence_alpha = LAMBDA_SMOOTH_STRUCT * smo_alpha,
    prior_residual = LAMBDA_PRIOR_MORT * pri_res, prior_age_structure = LAMBDA_PRIOR_STRUCT * pri_age)
}
decompose <- function(theta, wP = W_PREV_STRUCT, wD = W_DEATH_STRUCT, wI = W_INC_STRUCT, wS = W_STAGE_STRUCT) {
  p <- unpack(theta); pw <- project_window(p)
  LP <- loss_measure(pw$prev,  TGT$prev,  MSK$prev,  WT$prev)
  LD <- loss_measure(pw$death, TGT$death, MSK$death, WT$death)
  LI <- loss_measure(pw$inc,   TGT$inc,   MSK$inc,   WT$inc)
  LS <- stage_loss(pw)
  pri <- prior_terms(p, theta)
  terms <- c(prevalence_unweighted = LP, deaths_unweighted = LD, incidence_unweighted = LI, stage_unweighted = LS,
             prevalence_weighted = wP * LP, deaths_weighted = wD * LD, incidence_weighted = wI * LI, stage_weighted = wS * LS,
             pri)
  total <- sum(terms[c("prevalence_weighted", "deaths_weighted", "incidence_weighted", "stage_weighted")]) + sum(pri)
  list(total = total, terms = terms, pw = pw, p = p)
}
objective_factory <- function(wP, wD, wI, wS) function(theta) {
  v <- tryCatch(decompose(theta, wP, wD, wI, wS)$total, error = function(e) NA_real_)
  if (!is.finite(v)) 1e12 else v
}
objective      <- objective_factory(W_PREV_STRUCT, W_DEATH_STRUCT, W_INC_STRUCT, W_STAGE_STRUCT)
objective_prev <- objective_factory(W_PREV_STRUCT, 0, 0, 0)

#===============================================================================
# 4. OPTIMISE — multi-start L-BFGS-B (prior-centre, Layer-1 and prevalence-first
#    starts always evaluated); convergence codes kept
#===============================================================================
run_lbfgs <- function(p0, fn = objective, maxit = 150L) tryCatch(
  optim(p0, fn, method = "L-BFGS-B", lower = lo, upper = hi, control = list(factr = 1e9, maxit = maxit)),
  error = function(e) NULL)

t_start <- Sys.time()
starts <- list(prior_center = theta0, layer1_alpha = theta_L1)
if (PREV_FIRST) {                                    # prevalence-only preliminary fit -> a START only
  op <- run_lbfgs(theta0, objective_prev, maxit = 80L)
  if (!is.null(op)) starts$prevalence_first <- op$par
}
set.seed(SEED)
for (i in seq_len(N_STARTS)) starts[[paste0("random_", i)]] <- runif(npar, lo, hi)
fits <- lapply(starts, run_lbfgs)
ok_fit <- !vapply(fits, is.null, logical(1))
fits <- fits[ok_fit]; starts <- starts[ok_fit]
vals <- vapply(fits, function(f) f$value, numeric(1))
conv <- data.table(start = names(fits), objective = vals,
                   convergence_code = vapply(fits, function(f) f$convergence, numeric(1)),
                   message = vapply(fits, function(f) if (is.null(f$message)) "" else f$message, character(1)),
                   n_fn_evals = vapply(fits, function(f) f$counts[["function"]], numeric(1)))
best <- fits[[which.min(vals)]]
base_val <- objective(theta0)                        # guarantee not worse than the prior-centre baseline
if (base_val < best$value) best <- list(par = theta0, value = base_val, convergence = 0L)
theta_hat <- best$par; L_min <- best$value
dec <- decompose(theta_hat); bp <- dec$p; pw <- dec$pw
elapsed_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("  optimum L = %.4f (prior-centre baseline L = %.4f) | starts %d | spread %.4f/%.4f/%.4f | %.1f min\n",
            L_min, base_val, length(fits), min(vals), stats::median(vals), max(vals), elapsed_min))
cat("  loss decomposition at the optimum:\n"); print(round(dec$terms, 4))
prior_over_data <- sum(dec$terms[c("prior_mortality","prior_transitions","prior_incidence_alpha","smooth_incidence_alpha","prior_residual","prior_age_structure")]) /
  max(sum(dec$terms[c("prevalence_weighted","deaths_weighted","incidence_weighted","stage_weighted")]), 1e-12)
cat(sprintf("  summed prior+smoothness penalties / weighted data terms at optimum = %.3f (>1 => priors dominate)\n", prior_over_data))

#===============================================================================
# 5. FIT TABLES, ATTRIBUTION, STAGE STRUCTURE, D-FATALITY DIAGNOSTICS
#===============================================================================
ly <- n_yr
tot <- function(a) sum(a[, , ly]); tot_t <- function(a) sum(a[, , ly], na.rm = TRUE)
model_prev_tot <- tot(pw$prev); model_death_tot <- tot(pw$death); model_inc_tot <- tot(pw$inc)
gbd_prev_tot <- tot_t(TGT$prev); gbd_death_tot <- tot_t(TGT$death); gbd_inc_tot <- tot_t(TGT$inc)
prev_ratio <- model_prev_tot / gbd_prev_tot; death_ratio <- model_death_tot / gbd_death_tot
inc_ratio  <- model_inc_tot / gbd_inc_tot

# observed vs modelled by age group x sex x year, counts AND rates per 100k (WPP pop)
fit_long <- rbindlist(lapply(c("prev", "death", "inc"), function(m) {
  meas <- c(prev = "Prevalence", death = "Deaths", inc = "Incidence")[[m]]
  idxg <- expand.grid(g = seq_len(G), s = seq_len(n_sex), t = seq_len(n_yr))
  data.table(location = LOCATION, measure = meas, age_group = grp_present[idxg$g], sex = SEXES[idxg$s], year = yrs[idxg$t],
             model = pw[[m]][cbind(idxg$g, idxg$s, idxg$t)],
             target_used = TGT[[m]][cbind(idxg$g, idxg$s, idxg$t)],
             target_number = TGT_NUM[[m]][cbind(idxg$g, idxg$s, idxg$t)],
             pop_wpp = POP_G[cbind(idxg$g, idxg$s, idxg$t)],
             in_loss = (MSK[[m]] & WT[[m]] > 0)[cbind(idxg$g, idxg$s, idxg$t)] & (if (m == "inc") W_INC_STRUCT > 0 else TRUE),
             weight = WT[[m]][cbind(idxg$g, idxg$s, idxg$t)])
}))
fit_long[, `:=`(model_rate_per100k = fifelse(pop_wpp > 0, model / pop_wpp * 1e5, NA_real_),
                target_rate_per100k = fifelse(pop_wpp > 0, target_used / pop_wpp * 1e5, NA_real_),
                log_ratio = fifelse(!is.na(target_used) & target_used > 0, log((model + EPS_LOG) / (target_used + EPS_LOG)), NA_real_))]
fit_long[, held_out := year < MATCH_FROM]
# national summaries by year (both sexes) + by age group at the last year
fit_year <- fit_long[, .(model = sum(model), target = sum(target_used, na.rm = TRUE)), by = .(measure, year, held_out)]
fit_year[, ratio := model / target]
fit_age_last <- fit_long[year == CAL_YEAR_END, .(model = sum(model), target = sum(target_used, na.rm = TRUE)), by = .(measure, age_group)]
fit_age_last[, `:=`(ratio = model / target, abs_log_ratio = abs(log(pmax(model, EPS_LOG) / pmax(target, EPS_LOG))))]
# age-tail sensitivity: loss excluding the 85+ groups
tail_grp <- grp_present %in% c("85-89 years", "90-94 years", "95+ years")
msk_notail <- lapply(MSK, function(m) { m2 <- m; m2[tail_grp, , ] <- FALSE; m2 })
loss_notail <- c(prevalence = loss_measure(pw$prev, TGT$prev, msk_notail$prev, WT$prev),
                 deaths = loss_measure(pw$death, TGT$death, msk_notail$death, WT$death))

# stage stocks / shares / flows at the last year, per age (seed fallback for 05) ----
pw_s <- project_window(bp, keep_states = TRUE)
fs <- pw_s$start[[n_yr]]                                   # stocks ENTERING the final window year
stage_shares_by_age <- rbindlist(lapply(seq_len(n_sex), function(si) {
  tot_s <- fs$A[, si] + fs$B[, si] + fs$C[, si] + fs$D[, si]
  data.table(sex = SEXES[si], age = ages,
             fracA = ifelse(tot_s > 0, fs$A[, si] / tot_s, 1), fracB = ifelse(tot_s > 0, fs$B[, si] / tot_s, 0),
             fracC = ifelse(tot_s > 0, fs$C[, si] / tot_s, 0), fracD = ifelse(tot_s > 0, fs$D[, si] / tot_s, 0))
}))
tot_stock <- c(A = sum(fs$A), B = sum(fs$B), C = sum(fs$C), D = sum(fs$D))
stage_shares <- tot_stock / sum(tot_stock)
bp$dfrac <- unname(stage_shares[["D"]])
stage_stocks_year <- rbindlist(lapply(seq_len(n_yr), function(t) data.table(
  year = yrs[t], A = sum(pw$stages$A[[t]]), B = sum(pw$stages$B[[t]]), C = sum(pw$stages$C[[t]]), D = sum(pw$stages$D[[t]]),
  rhd_deaths = sum(pw$death[, , t]), rhd_deaths_D = sum(pw$deathD[, , t]), incident = sum(pw$inc[, , t]))))
stage_stocks_year[, `:=`(share_D = D / (A + B + C + D), deaths_outside_D = rhd_deaths - rhd_deaths_D)]
# realised D fatality on START-of-cycle D stock (undefined -> NA with flag; zero deaths required)
d_fat <- rbindlist(lapply(seq_len(n_yr), function(t) {
  Ds <- sum(pw_s$start[[t]]$D); dd <- sum(pw$deathD[, , t])
  data.table(year = yrs[t], D_start = Ds, rhd_deaths_D = dd,
             realized_D_death_risk = if (Ds > 0) dd / Ds else NA_real_,
             D_start_zero = Ds <= 0, zero_deaths_when_undefined_ok = if (Ds > 0) NA else dd == 0)
}))
cf_implied <- model_death_tot / model_prev_tot
gbd_cf     <- gbd_death_tot / gbd_prev_tot
cf_gap     <- cf_implied / gbd_cf
resid_at_bound <- if (D_ONLY) FALSE else (abs(bp$resid - RESID_LO) < 1e-3) || (abs(bp$resid - RESID_HI) < 1e-3)
# CF floor: everyone in stage A at the (effective) A mortality. In D-only mode this
# is exactly 0 (any aggregate CF is attainable with a small enough D share).
min_plausible_cf <- unname(bp$pd[["A"]] * bp$resid)
mD_at_bound <- (abs(log(bp$mD) - log(MD_LO)) < 1e-3) || (abs(log(bp$mD) - log(MD_HI)) < 1e-3)

attribution <- if (resid_at_bound)
    "residual mortality multiplier AT BOUND (VALIDATION FAILURE — revisit priors/bounds/stock)" else
  if (abs(log(death_ratio)) < FIT_TOL_DEATH && abs(log(prev_ratio)) < FIT_TOL_PREV)
    sprintf("GBD prevalence AND RHD deaths reproduced within tolerance (|log ratio| < %.2f / %.2f) under %s mortality; D share %.1f%% vs GBD HF-due-to-RHD proxy %s",
            FIT_TOL_PREV, FIT_TOL_DEATH, MORT_STRUCT, 100 * stage_shares[["D"]],
            if (is.na(dmi$meta$d_share_proxy_gbd_hf)) "n/a" else sprintf("%.1f%%", 100 * dmi$meta$d_share_proxy_gbd_hf)) else
  if (!D_ONLY && gbd_cf < min_plausible_cf)
    sprintf(paste0("GBD aggregate case fatality (%.2g/yr) is BELOW the lowest clinically-plausible ",
      "stage mortality (m_A=%.2g/yr): even a mostly-A equilibrium (A=%.0f%%) implies CF=%.2g/yr (%.1fx GBD)."),
      gbd_cf, min_plausible_cf, 100*stage_shares[["A"]], cf_implied, cf_gap) else
  if (abs(log(prev_ratio)) >= FIT_TOL_PREV)
    sprintf("TOTAL PREVALENCE outside tolerance (model %.2fx GBD): incidence (alpha %.2f-%.2f) x dwell (p_A_to_no_rhd=%.3f) vs GBD prevalence.",
            prev_ratio, min(bp$alpha), max(bp$alpha), bp$tr$p_A_to_no_rhd) else
  if (abs(log(death_ratio)) >= FIT_TOL_DEATH)
    sprintf("RHD DEATHS outside tolerance (model %.2fx GBD) with m_D=%.3f (prior %.3f%s) and D share %.2f%%.",
            death_ratio, bp$mD, pd_prior[["D"]], if (mD_at_bound) ", AT BOUND" else "", 100 * stage_shares[["D"]]) else
  "RHD deaths reconciled by stock/progression at prior-consistent mortality"

#===============================================================================
# 6. PARTIAL IDENTIFICATION — near-optimal set + one-at-a-time flatness sweep
#===============================================================================
acc <- which(vals <= L_min * (1 + NEAROPT_DELTA))
accepted <- rbindlist(lapply(acc, function(k) {
  pk <- unpack(fits[[k]]$par)
  data.table(L = vals[k], m_A = pk$pd[["A"]], m_B = pk$pd[["B"]], m_C = pk$pd[["C"]],
             m_D = pk$pd[["D"]], resid = pk$resid,
             p_A_to_no_rhd = pk$tr$p_A_to_no_rhd, p_B_to_C = pk$tr$p_B_to_C, p_C_to_D = pk$tr$p_C_to_D)
}), fill = TRUE)
accepted_full <- rbindlist(lapply(acc, function(k) {   # COMPLETE parameter sets (incl. incidence)
  pk <- unpack(fits[[k]]$par)
  cbind(data.table(start = names(fits)[k], L = vals[k], convergence_code = fits[[k]]$convergence),
        as.data.table(as.list(unlist(pk$tr))), data.table(m_A = pk$pd[["A"]], m_B = pk$pd[["B"]], m_C = pk$pd[["C"]], m_D = pk$pd[["D"]], resid = pk$resid,
                                                           adult_prog_mult = pk$age[["prog"]], adult_remission_mult = pk$age[["rem"]]),
        as.data.table(as.list(setNames(as.vector(pk$alpha), alpha_names))))
}), fill = TRUE)
sweep_frac <- function(j, npts = 11) {
  grid <- seq(lo[j], hi[j], length.out = npts)
  Ls <- vapply(grid, function(v) { th <- theta_hat; th[j] <- v; objective(th) }, numeric(1))
  mean(Ls <= L_min * (1 + NEAROPT_DELTA))
}
identifiability <- data.table(parameter = param_labels, flatness = vapply(seq_len(npar), sweep_frac, numeric(1)))
identifiability[, status := fifelse(flatness < 0.5, "data-identified", "prior-dominated")]
identifiability[, note := "one-at-a-time profile flatness of the PENALISED objective (others at optimum); not joint identification and not proof of data identification"]

# parameter table: value, prior centre, bounds, boundary hits
par_tab <- data.table(parameter = param_labels, theta = theta_hat, lower = lo, upper = hi)
par_tab[, `:=`(at_lower = abs(theta - lower) < 1e-3, at_upper = abs(theta - upper) < 1e-3)]
par_tab[, value := c(unlist(bp$tr[trans_names]), bp$mD, if (!D_ONLY) c(bp$ratios[c("C","B","A")], bp$resid),
                     if (CALIB_AGE) c(bp$age[["prog"]], bp$age[["rem"]]), as.vector(bp$alpha))]
par_tab[, prior_center := c(unlist(tr_prior[trans_names]), pd_prior[["D"]],
                            if (!D_ONLY) c(pd_prior[["C"]]/pd_prior[["D"]], pd_prior[["B"]]/pd_prior[["C"]], pd_prior[["A"]]/pd_prior[["B"]], 1),
                            if (CALIB_AGE) c(1, REM_MULT),
                            rep(1, n_alpha))]
par_tab[, source := c(rep("00 block G prior centre; bounds = literature workbook annualised ranges (widened)", n_tr),
                      dmi$meta$p_rhd_death_D_prior_source,
                      if (!D_ONLY) rep("legacy nested-logistic severity ratio / residual", 4),
                      if (CALIB_AGE) rep(sprintf("adult (age >= %d) hazard multiplier; weak prior (log-sd SIGMA_AGE_MULT); child-cohort evidence base", REM_CUT), 2),
                      rep("weak prior on GBD incidence pattern (log-sd SIGMA_ALPHA_STRUCT); Layer-1 alpha = start only", n_alpha))]

#===============================================================================
# 7. INITIALISATION CHECKS + BRIDGE TO THE ANALYSIS START (06 seed)
#===============================================================================
burn_conv <- as.data.table(pw$burn$convergence)
pw_long <- project_window(bp, burn_years = 2L * BURN_YEARS)
burnin_check <- data.table(stage = c("A", "B", "C", "D"),
  stock_burn = c(sum(pw$stages$A[[1]]), sum(pw$stages$B[[1]]), sum(pw$stages$C[[1]]), sum(pw$stages$D[[1]])),
  stock_burn_2x = c(sum(pw_long$stages$A[[1]]), sum(pw_long$stages$B[[1]]), sum(pw_long$stages$C[[1]]), sum(pw_long$stages$D[[1]])))
burnin_check[, rel_diff_2x := abs(stock_burn_2x - stock_burn) / pmax(stock_burn, 1e-9)]
burnin_check <- merge(burnin_check, burn_conv, by = "stage")
L_burn2x <- loss_measure(pw_long$prev, TGT$prev, MSK$prev, WT$prev) * W_PREV_STRUCT +
            loss_measure(pw_long$death, TGT$death, MSK$death, WT$death) * W_DEATH_STRUCT
burn_converged <- all(burnin_check$rel_diff_2x < 0.01)

initial_state <- NULL; bridge_recon <- NULL
if (BRIDGE) {
  bridge_years <- (CAL_YEAR_END + 1L):(ANALYSIS_YEAR_START - 1L)
  if (length(bridge_years) && bridge_years[1] <= bridge_years[length(bridge_years)]) {
    pop_obs  <- as.data.table(readRDS(POP_OBS_FILE))[location == LOCATION]
    pop_proj <- as.data.table(readRDS(POP_PROJ_FILE))[location == LOCATION]
    pop_all  <- rbind(pop_obs[year %in% bridge_years], pop_proj[year %in% bridge_years & !year %in% pop_obs$year])
    miss_y <- setdiff(bridge_years, unique(pop_all$year))
    if (length(miss_y)) stop("04b bridge: population missing for year(s) ", paste(miss_y, collapse = ","), call. = FALSE)
    nb <- length(bridge_years)
    POPb <- array(0, dim = c(n_age, n_sex, nb), dimnames = list(ages, SEXES, bridge_years))
    POPb[cbind(match(pop_all$age, ages), match(pop_all$sex, SEXES), match(pop_all$year, bridge_years))] <- pop_all$Nx
    if (any(POPb <= 0)) stop("04b bridge: zero/negative population cells.", call. = FALSE)
    trend <- dmi$meta$incidence_trend
    ir_anchor <- ir_engine_arr(bp$alpha)[, , n_yr]                 # per-susceptible, calibrated, at CAL_YEAR_END
    IRb <- array(0, dim = dim(POPb), dimnames = dimnames(POPb))
    for (t in seq_len(nb)) IRb[, , t] <- ir_anchor * trend^(bridge_years[t] - CAL_YEAR_END)
    OTHb <- dmi$rates_by_year$oth_mort[, , 1]                       # production background pattern (held)
    covb <- list(treatment = rep(C0_TREAT_EFF, nb), surgery = rep(C0_SURG, nb), adherence = rep(C0_ADH, nb))
    prb <- abcd_project(pw$final_state, POPb, IRb, OTHb, tr_engine(bp$tr, bp$age), bp$pd * bp$resid, eff_cal, surg_cal,
                        covb, bridge_years, light = TRUE, prev_obs = PREV_OBS)
    initial_state <- list(year = ANALYSIS_YEAR_START, A = prb$final_state$A, B = prb$final_state$B,
                          C = prb$final_state$C, D = prb$final_state$D,
                          method = sprintf("calibrated end-of-%d state forward-projected %d-%d with reference care, calibrated incidence x trend^(y-%d), held 03 background mortality, WPP population",
                                           CAL_YEAR_END, min(bridge_years), max(bridge_years), CAL_YEAR_END),
                          pop_deficit_total = prb$pop_deficit_total)
    for (s in c("A", "B", "C", "D")) dimnames(initial_state[[s]]) <- list(ages, SEXES)
    # reconciliation vs GBD in the bridge years where GBD exists (Number and Rate x WPP pop)
    gbd_all <- as.data.table(readRDS(paste0(wd_data, "temp_baseline_rates_gbd.rds")))
    gb <- gbd_all[location_name == LOCATION & cause_name == "Rheumatic heart disease" & year %in% bridge_years &
                    measure_name %in% c("Prevalence", "Deaths"), .(measure = measure_name, metric = metric_name, sex = sex_name, age_group = age_name, year, val)]
    popg <- pop_all[, .(pop = sum(Nx)), by = .(sex, age_group = age_to_gbd_group(age), year)]
    gbn <- gb[metric == "Number", .(gbd_number = sum(val)), by = .(measure, year)]
    gbr <- merge(gb[metric == "Rate"], popg, by = c("sex", "age_group", "year"))[, .(gbd_rate_x_wpp = sum(val / 1e5 * pop)), by = .(measure, year)]
    modb <- rbindlist(lapply(seq_len(nb), function(t) data.table(
      year = bridge_years[t], measure = c("Prevalence", "Deaths"),
      model = c(sum(prb$prev_obs[[t]]), sum(prb$cycles[[t]]$rhd_deaths)))))
    bridge_recon <- merge(merge(modb, gbn, by = c("measure", "year"), all.x = TRUE), gbr, by = c("measure", "year"), all.x = TRUE)
    bridge_recon[, `:=`(ratio_model_over_number = model / gbd_number, ratio_model_over_rate_x_wpp = model / gbd_rate_x_wpp)]
    seed_tot <- sum(initial_state$A + initial_state$B + initial_state$C + initial_state$D)
    legacy_seed_tot <- sum(dmi$rates_by_year$prev_seed * {   # legacy: GBD 2023 prevalence x analysis-start WPP pop
      p1 <- pop_proj[year == ANALYSIS_YEAR_START]; m <- matrix(0, n_age, n_sex, dimnames = list(ages, SEXES))
      m[cbind(match(p1$age, ages), match(p1$sex, SEXES))] <- p1$Nx; m })
    bridge_recon <- rbind(bridge_recon, data.table(measure = "Prevalence (seed entering analysis start)", year = ANALYSIS_YEAR_START,
      model = seed_tot, gbd_number = NA_real_, gbd_rate_x_wpp = legacy_seed_tot,
      ratio_model_over_number = NA_real_, ratio_model_over_rate_x_wpp = seed_tot / legacy_seed_tot), fill = TRUE)
    cat(sprintf("  bridge %d-%d: seed entering %d = %.0f living RHD (%.2fx the legacy GBD-%d-rate x pop seed); D share %.2f%%\n",
                min(bridge_years), max(bridge_years), ANALYSIS_YEAR_START, seed_tot, seed_tot / legacy_seed_tot,
                dmi$meta$RATE_BASE_YEAR, 100 * sum(initial_state$D) / seed_tot))
  }
}

#===============================================================================
# 8. STATUSES + VALIDATION (structural = hard; numerical/fit = labels)
#===============================================================================
structural_violations <- character(0)
pdv <- bp$pd * bp$resid
if (D_ONLY) {
  if (any(pdv[c("A","B","C")] != 0)) structural_violations <- c(structural_violations, "D-only: A/B/C RHD mortality not exactly zero")
  if (!(pdv[["D"]] > 0 && pdv[["D"]] < 1)) structural_violations <- c(structural_violations, "D-only: invalid D mortality")
  if (any(stage_stocks_year$deaths_outside_D > 1e-9)) structural_violations <- c(structural_violations, "D-only: RHD deaths occurred outside stage D")
  if (bp$resid != 1) structural_violations <- c(structural_violations, "D-only: residual multiplier not fixed at 1")
} else {
  mm <- bp$pd
  if (!((mm[["A"]] <= mm[["B"]] + 1e-12) && (mm[["B"]] < mm[["C"]]) && (mm[["C"]] < mm[["D"]])))
    structural_violations <- c(structural_violations, "all_stage: mortality ordering 0<=m_A<=m_B<m_C<m_D violated")
  if (resid_at_bound) structural_violations <- c(structural_violations, sprintf("all_stage: residual multiplier at bound (%.3f in [%.2f,%.2f])", bp$resid, RESID_LO, RESID_HI))
}
if (any(!is.finite(unlist(stage_stocks_year))) || any(unlist(stage_stocks_year[, .(A, B, C, D, rhd_deaths, incident)]) < -1e-9))
  structural_violations <- c(structural_violations, "non-finite or negative stocks/flows")
if (pw$pop_deficit > 1e-6) structural_violations <- c(structural_violations, sprintf("RHD stock exceeded population (total deficit %.3g)", pw$pop_deficit))
if (d_fat[D_start_zero == TRUE & zero_deaths_when_undefined_ok == FALSE, .N] > 0)
  structural_violations <- c(structural_violations, "RHD deaths from D with zero D start stock")
if (!is.null(initial_state) && initial_state$pop_deficit_total > 1e-6)
  structural_violations <- c(structural_violations, "bridge: RHD stock exceeded population")
structural_status <- if (length(structural_violations)) paste("FAILED:", paste(structural_violations, collapse = "; ")) else "ok"
numerical_status  <- sprintf("best start '%s' convergence code %d; %d/%d starts converged (code 0); spread %.4f-%.4f",
                             names(fits)[which.min(vals)], as.integer(best$convergence), sum(conv$convergence_code == 0), nrow(conv), min(vals), max(vals))
fit_status <- sprintf("prevalence %s tolerance (|log ratio| %.3f vs %.2f); deaths %s tolerance (%.3f vs %.2f); incidence ratio %.2f (validation only)",
                      if (abs(log(prev_ratio)) < FIT_TOL_PREV) "within" else "OUTSIDE", abs(log(prev_ratio)), FIT_TOL_PREV,
                      if (abs(log(death_ratio)) < FIT_TOL_DEATH) "within" else "OUTSIDE", abs(log(death_ratio)), FIT_TOL_DEATH, inc_ratio)
boundary_hits <- par_tab[at_lower | at_upper, parameter]

status <- if (length(structural_violations)) "FAILED_structural_validation" else
  if (n_stage_used > 0L) "calibrated_with_echo_targets" else "calibrated_no_echo_targets_partial_identification"

# human-readable diagnostics (written in every case, incl. failure)
fwrite(fit_long,        paste0(wd_data, "calibration_structural_fit_by_group_year.csv"))
fwrite(fit_year,        paste0(wd_data, "calibration_structural_fit_by_year.csv"))
fwrite(fit_age_last,    paste0(wd_data, "calibration_structural_fit_by_age_last_year.csv"))
fwrite(stage_stocks_year, paste0(wd_data, "calibration_structural_stage_stocks_by_year.csv"))
fwrite(d_fat,           paste0(wd_data, "calibration_structural_d_fatality.csv"))
fwrite(par_tab,         paste0(wd_data, "calibration_structural_parameters.csv"))
fwrite(conv,            paste0(wd_data, "calibration_structural_convergence.csv"))
fwrite(burnin_check,    paste0(wd_data, "calibration_structural_burnin_check.csv"))
fwrite(data.table(term = names(dec$terms), value = dec$terms), paste0(wd_data, "calibration_structural_loss_decomposition.csv"))
if (!is.null(bridge_recon)) fwrite(bridge_recon, paste0(wd_data, "calibration_structural_bridge_reconciliation.csv"))

if (length(structural_violations)) {
  failed <- list(status = status, reason = structural_violations, candidate = bp, theta = theta_hat, L = L_min,
                 decomposition = dec$terms, spec = list(mortality_structure = MORT_STRUCT, sap_mechanism = SAP_MECH))
  saveRDS(failed, paste0(wd_data, "calibration_structural_failed_candidate.rds"))
  fwrite(data.table(parameter = par_tab$parameter, value = par_tab$value), paste0(wd_data, "calibration_structural_failed_candidate.csv"))
  calib$stage_calibration$status <- status
  calib$stage_calibration$failure_reason <- structural_violations
  calib$stage_calibration$best_parameters <- NULL          # a failed candidate is NEVER written as calibrated
  calib$meta$stage_calibrated <- FALSE
  saveRDS(calib, BUNDLE_FILE)
  stop("04b STRUCTURAL VALIDATION FAILED: ", paste(structural_violations, collapse = "; "),
       "\n  Failed candidate saved to calibration_structural_failed_candidate.rds; bundle status = ", status,
       if (file.exists(LASTVALID_FILE)) paste0("\n  Last valid bundle kept at ", basename(LASTVALID_FILE)) else "", call. = FALSE)
}

#===============================================================================
# 9. ASSEMBLE best_parameters + write into the bundle $stage_calibration
#===============================================================================
ir_anchor_final <- ir_engine_arr(bp$alpha)[, , n_yr]; dimnames(ir_anchor_final) <- list(ages, SEXES)
alpha_tab <- data.table(location = LOCATION, band = rep(band_labels, times = n_sex), band_lo = rep(band_edges, times = n_sex),
                        sex = rep(SEXES, each = n_band), alpha = as.vector(bp$alpha),
                        alpha_layer1_start = as.vector(alpha_L1))
best_parameters <- list(
  transitions = bp$tr,                      # scalar annual probabilities (legacy contract)
  # p_rhd_death consumed by 06 = untreated estimate x residual multiplier. D-only:
  # A/B/C are EXACT structural zeros and the residual is FIXED at 1 (one D scale).
  p_rhd_death = bp$pd * bp$resid,
  p_rhd_death_untreated = bp$pd,
  rhd_d_fraction = bp$dfrac,                # = model D-share entering the last window year (derived)
  residual_mortality_mult = bp$resid,       # D-only: fixed at 1 (documented fixed field, not estimated)
  residual_mortality_mult_role = if (D_ONLY) "fixed at 1 (not estimated; one D-mortality scale)" else "estimated, bounded",
  stage_split = stage_shares,               # aggregate model A/B/C/D shares entering the last window year
  stage_shares_by_age = stage_shares_by_age,# per-age shares (legacy seed fallback for 05)
  aggregate_stage_shares = stage_shares,
  # FINAL incidence (joint Layer-2 estimate; applied ONCE; 05 consumes these)
  incidence_alpha = alpha_tab,
  ir_calibrated = ir_anchor_final,          # [age x sex] per-SUSCEPTIBLE annual probability at ir_anchor_year
  ir_anchor_year = CAL_YEAR_END,
  ir_denominator = "per susceptible (No-RHD) person: IR_GBD(total pop) x alpha / (1 - prev_GBD)",
  initial_state = initial_state,            # stocks ENTERING ANALYSIS_YEAR_START (NULL if bridge disabled)
  mortality_structure = MORT_STRUCT,
  sap_mechanism = SAP_MECH,
  remission_adult_mult = REM_MULT, remission_age_cutoff = REM_CUT,
  # explicit AGE STRUCTURE applied by 05/06 exactly as in calibration (hazard-scale
  # multipliers at ages >= adult_age_cutoff); transitions above = child/reference values
  age_structure = list(adult_prog_mult = unname(bp$age[["prog"]]), adult_remission_mult = unname(bp$age[["rem"]]),
                       adult_age_cutoff = REM_CUT, calibrated = CALIB_AGE,
                       progression_transitions = c("p_A_to_B", "p_B_to_C", "p_C_to_D"), remission_transition = "p_A_to_no_rhd",
                       scale = "hazard: p_adult = 1 - (1 - p)^mult"),
  reference_care_in_calibration = list(treatment_effective = C0_TREAT_EFF, adherence = C0_ADH, surgery = C0_SURG,
                                       interpretation = dmi$meta$baseline_care_interpretation))

calib$stage_calibration$status              <- status
calib$stage_calibration$failure_reason      <- NULL
calib$stage_calibration$best_parameters     <- best_parameters
calib$stage_calibration$accepted_parameter_sets <- accepted
calib$stage_calibration$accepted_parameter_sets_full <- accepted_full
calib$stage_calibration$identifiability     <- identifiability
calib$stage_calibration$parameters          <- par_tab
calib$stage_calibration$convergence         <- conv
calib$stage_calibration$objective_components <- list(
  L_min = L_min, prior_center_L = base_val, n_starts = N_STARTS,
  multistart_spread = c(min = min(vals), median = stats::median(vals), max = max(vals)),
  weights = c(W_PREV = W_PREV_STRUCT, W_DEATH = W_DEATH_STRUCT, W_INC = W_INC_STRUCT, W_STAGE = W_STAGE_STRUCT),
  prev_death_ratio = W_PREV_STRUCT / W_DEATH_STRUCT, data_weight_total = W_PREV_STRUCT + W_DEATH_STRUCT,
  match_from = MATCH_FROM, prev_obs_convention = PREV_OBS, target_basis = calib$meta$layer2_target_basis,
  loss_decomposition = dec$terms, prior_over_data_ratio = prior_over_data,
  loss_excluding_85plus = loss_notail, L_with_2x_burnin = L_burn2x,
  formula = "L = wP*MSLE(prev) + wD*MSLE(deaths) + wI*MSLE(inc) + LAMBDA*[priors] + LAMBDA_SMOOTH*smooth(alpha) + wS*stage; MSLE = weighted mean over masked cells (both sexes pooled; age x year [x capped IV] weights) of (log(m+eps)-log(t+eps))^2",
  masking = sprintf("missing GBD cells masked (NA, never zero-filled); counts < %g masked; years < %d weight 0 (held-out); eps = %g; death-target age floor = %d (%d cells masked, reported as validation)",
                    TARGET_MIN_COUNT, MATCH_FROM, EPS_LOG, DEATH_AGE_MIN, n_death_age_masked),
  death_target_age_min = DEATH_AGE_MIN, lambda_prior_structural = LAMBDA_PRIOR_STRUCT, lambda_prior_mortality = LAMBDA_PRIOR_MORT,
  iv_weights = list(used = USE_IV_W, cap = IV_CAP), age_weights = if (is.null(AGE_W)) "equal" else AGE_W,
  elapsed_min = elapsed_min)
calib$stage_calibration$priors <- list(
  mortality_centers = pd_prior, sigma_mort = SIGMA_MORT,
  transition_centers = tr_prior, sigma_trans = SIGMA_TRANS,
  sigma_resid = SIGMA_RESID,
  sigma_alpha = SIGMA_ALPHA_STRUCT, lambda_smooth_alpha = LAMBDA_SMOOTH_STRUCT, alpha_centers = "1 (GBD incidence pattern)",
  rhd_d_fraction_note = "derived as the model D-share entering the last window year (not a free parameter)",
  mortality_note = if (D_ONLY) "D-only: single log-normal prior on m_D; A/B/C fixed at exact 0 (no prior terms, no log(0))" else "all_stage: log-normal priors on all four stage mortalities",
  burn_years = BURN_YEARS)
calib$stage_calibration$bounds <- list(
  m_D = c(MD_LO, MD_HI), residual_mult = if (D_ONLY) c(1, 1) else c(RESID_LO, RESID_HI),
  transitions = TRANS_BOUNDS[trans_names], alpha = c(ALPHA_LO, ALPHA_HI), boundary_hits = boundary_hits)
calib$stage_calibration$excess_death_attribution <- list(
  model_prev = model_prev_tot, gbd_prev = gbd_prev_tot, prev_ratio = prev_ratio,
  model_deaths = model_death_tot, gbd_deaths = gbd_death_tot, death_ratio = death_ratio,
  model_inc = model_inc_tot, gbd_inc = gbd_inc_tot, inc_ratio = inc_ratio,
  cf_implied = cf_implied, gbd_cf = gbd_cf, cf_gap = cf_gap,
  min_plausible_cf = min_plausible_cf, stage_shares_last = stage_shares,
  d_share_proxy_gbd_hf = dmi$meta$d_share_proxy_gbd_hf,
  residual_mult = bp$resid, residual_at_bound = resid_at_bound, mD_at_bound = mD_at_bound,
  classification = attribution)
calib$stage_calibration$fit_by_group_year <- fit_long
calib$stage_calibration$fit_by_year <- fit_year
calib$stage_calibration$stage_stocks_by_year <- stage_stocks_year
calib$stage_calibration$d_fatality <- d_fat
calib$stage_calibration$burnin_check <- burnin_check
calib$stage_calibration$bridge_reconciliation <- bridge_recon
calib$stage_calibration$stage_targets_n_used <- n_stage_used
calib$stage_calibration$statuses <- list(structural = structural_status, numerical = numerical_status,
                                         epidemiological_fit = fit_status, burn_in_converged = burn_converged)
calib$stage_calibration$spec <- list(mortality_structure = MORT_STRUCT, sap_mechanism = SAP_MECH,
                                     sap_hr = unlist(eff_cal[c("sap_hr_A_to_B", "sap_hr_B_to_C", "sap_hr_C_to_D")]),
                                     reference_care = best_parameters$reference_care_in_calibration,
                                     prev_obs = PREV_OBS, target_basis = calib$meta$layer2_target_basis,
                                     incidence_free = W_INC_STRUCT == 0,
                                     layer1_role = "aggregate proxy; alpha used as starting values only; CF diagnostic/legacy")
calib$meta$stage_calibrated <- TRUE
calib$meta$layer2_spec <- calib$stage_calibration$spec

saveRDS(calib, BUNDLE_FILE)
saveRDS(calib, LASTVALID_FILE)                       # last VALIDATED bundle (never overwritten by a failed run)

fwrite(accepted,        paste0(wd_data, "calibration_structural_accepted_sets.csv"))
fwrite(accepted_full,   paste0(wd_data, "calibration_structural_accepted_sets_full.csv"))
fwrite(identifiability, paste0(wd_data, "calibration_structural_identifiability.csv"))
fwrite(alpha_tab,       paste0(wd_data, "calibration_structural_incidence_alpha.csv"))
fwrite(data.table(parameter = c(names(unlist(bp$tr)), names(bp$pd), "rhd_d_fraction", "residual_mult", alpha_names),
                  value = c(unlist(bp$tr), bp$pd * bp$resid, bp$dfrac, bp$resid, as.vector(bp$alpha))),
       paste0(wd_data, "calibration_structural_best_parameters.csv"))

#===============================================================================
# 10. REPORT
#===============================================================================
cat("\n", strrep("=", 70), "\nSTRUCTURAL (LAYER-2) VALIDATION — ", MORT_STRUCT, " mortality | SAP ", SAP_MECH, "\n", strrep("=", 70), "\n", sep = "")
mm <- bp$pd * bp$resid
cat(sprintf("Stage RHD mortality (/yr, effective): A=%.4f B=%.4f C=%.4f D=%.4f | %s\n", mm[["A"]], mm[["B"]], mm[["C"]], mm[["D"]],
            if (D_ONLY) "exact A=B=C=0 verified" else "ordering A<=B<C<D verified"))
cat(sprintf("m_D in [%.2f,%.2f] (prior %.3f%s) | residual mult = %.3f (%s)\n", MD_LO, MD_HI, pd_prior[["D"]],
            if (mD_at_bound) ", AT BOUND" else "", bp$resid, best_parameters$residual_mortality_mult_role))
cat(sprintf("Transitions (age < %d): %s\n", REM_CUT, paste(sprintf("%s=%.4f", names(unlist(bp$tr)), unlist(bp$tr)), collapse=", ")))
cat(sprintf("Adult (age >= %d) hazard multipliers: progression x%.3f | remission x%.3f (%s)\n", REM_CUT, bp$age[["prog"]], bp$age[["rem"]],
            if (CALIB_AGE) "calibrated" else "fixed"))
cat(sprintf("Incidence alpha (band x sex): %s\n", paste(sprintf("%s=%.2f", alpha_names, as.vector(bp$alpha)), collapse = ", ")))
cat(sprintf("Model stage shares entering %d: A=%.1f%% B=%.1f%% C=%.1f%% D=%.2f%% (GBD HF-due-to-RHD proxy D share: %s)\n", CAL_YEAR_END,
            100*stage_shares[["A"]],100*stage_shares[["B"]],100*stage_shares[["C"]],100*stage_shares[["D"]],
            if (is.na(dmi$meta$d_share_proxy_gbd_hf)) "n/a" else sprintf("%.2f%%", 100*dmi$meta$d_share_proxy_gbd_hf)))
cat(sprintf("\nFit at %d (totals): prev model=%.0f vs target %.0f (%.2fx) | RHD deaths %.0f vs %.0f (%.2fx) | incidence %.0f vs %.0f (%.2fx, validation)\n",
            CAL_YEAR_END, model_prev_tot, gbd_prev_tot, prev_ratio, model_death_tot, gbd_death_tot, death_ratio, model_inc_tot, gbd_inc_tot, inc_ratio))
cat("Age-specific fit at the last year (model/target ratio):\n")
print(dcast(fit_age_last[, .(measure, age_group, ratio = round(ratio, 2))], age_group ~ measure, value.var = "ratio")[match(grp_present, age_group)])
cat(sprintf("Implied total CF = %.5f vs GBD CF %.5f | realised D death risk (last yr) = %s\n", cf_implied, gbd_cf,
            if (is.na(d_fat[.N, realized_D_death_risk])) "undefined (D_start = 0)" else sprintf("%.4f", d_fat[.N, realized_D_death_risk])))
cat(sprintf("EXCESS-DEATH ATTRIBUTION: %s\n", attribution))
cat(sprintf("STATUS structural: %s\nSTATUS numerical : %s\nSTATUS fit       : %s\nBurn-in converged (2x check, <1%% change): %s | boundary hits: %s\n",
            structural_status, numerical_status, fit_status, burn_converged,
            if (length(boundary_hits)) paste(boundary_hits, collapse = ", ") else "none"))
cat("\nParameter identifiability (one-at-a-time flatness; not joint identification):\n"); print(identifiability[, .(parameter, flatness, status)])
cat(sprintf("\nStage targets: detected %d | usable %d | USED in objective %d -> status: %s | near-optimal sets kept: %d\n",
            calib$stage_calibration$stage_targets_n_detected, calib$stage_calibration$stage_targets_n_usable, n_stage_used, status, nrow(accepted)))
if (n_stage_used == 0L)
  cat("PARTIAL IDENTIFICATION: no usable local echo stage targets -> stage split is prior-informed;\n",
      "  see accepted_parameter_sets(_full) + identifiability. Populate calibration_targets_stage.csv to pin it.\n", sep="")
cat("04b complete. 05 will bake best_parameters (incl. the FINAL incidence and the analysis-start initial state) into the engine.\n")

}  # end if RUN_STRUCT
