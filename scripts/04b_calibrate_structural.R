#===============================================================================
# 04b_calibrate_structural.R
#-------------------------------------------------------------------------------
# LAYER-2 STRUCTURAL CALIBRATION — POPULATION BURDEN + FIRST DIAGNOSIS,
# PREVALENCE-LED, ON THE PRODUCTION ENGINE
#
# Runs the ACTUAL production engine (R/abcd_engine.R: abcd_one_cycle + the shared
# abcd_burn_in / abcd_project / abcd_age_state wrappers and the reference care, the
# same ones 06 uses) inside the objective, and JOINTLY estimates
#   * low-dimensional INCIDENCE corrections (one log-multiplier per broad age band
#     x sex on the GBD incidence-per-total-population pattern; Layer-1 alpha are
#     STARTING VALUES only), and
#   * the six non-zero annual transition probabilities (shared across sexes; the
#     CHILDHOOD values: adult (> 20) B->C and C->D = 1 - (1 - p)^adult_prog_hazard_mult
#     (fixed m, 00 block G); A->B also drives adult prior-A -> B unscaled; prior
#     centres = 00 block G; bounds TRANS_BOUNDS on the childhood values),
# against GBD RHD PREVALENCE (principal target: TOTAL classifiable disease, diagnosed +
# undiagnosed; adult prior-A and normalised history excluded) and RHD-SPECIFIC DEATHS
# (diagnosed + undiagnosed D), with incidence REPORTED as a validation output.
# FIXED (never fitted): D mortality (one input for every history), entry mixes, adult
# B normalisation, the adult progression hazard multiplier, D->C = 0 and the residual
# multiplier 1. Reference care is applied
# INSIDE the objective exactly as in production; fitted parameters are untreated.
#
# OBJECTIVE (every term exported):
#   L = wP * L_prev + wD * L_death + wI * MSLE(inc) + LAMBDA*[prior(trans) +
#       prior(alpha)] + LAMBDA_SMOOTH * smooth(alpha) + wS * stage-target error
#   L_prev  = p * MSLE(prevalence, ages 0-19) + (1 - p) * MSLE(prevalence, ages 20+)
#   L_death = s * MSLE(deaths, ages 0-19)     + (1 - s) * MSLE(deaths, ages 20+)
#   MSLE = weighted mean over USABLE cells (both sexes pooled) of
#          (log(model + eps) - log(target + eps))^2. Usable = observed (NA masked),
#          count >= TARGET_MIN_COUNT, scored year (>= STRUCT_MATCH_FROM).
#   Each component is normalised over its OWN usable cells; the script stops if any
#   component has none (never replaced by a zero loss).
#
# TIME CONVENTION (shared with 06 via abcd_project): start state -> care update ->
# within-year flows -> end state -> ageing (+ A -> adult prior-A reclassification).
#
# INITIALISATION: zero-start burn-in under first-window-year inputs with the
# provisional initial care shares held as static labels; at CAL_YEAR_START the
# dynamic care state is initialised ONCE (diagnosed share of every stage stock, on-BPG
# share of the SAP-eligible stage stocks, nobody ever screened) and carried forward
# through the window and the bridge (2020-2025) into the analysis start — never reset.
# Burn-in convergence (stage / prior-A / normalised totals, age profiles, deaths) is
# CHECKED against a 2x run at the best solution. The bridged state (all stage,
# diagnosis, BPG, prior-A, normalised and screening histories) is the 06 seed.
# Burn-in reconciliation follows STRUCT_BURNIN_RECONCILIATION
# ("per_capita" default: all count stocks rescaled to the constant burn-in population
# each cycle; the adjustments are saved as initialisation balancing and the hand-off
# accounting is a structural check). Observed years and the bridge never rescale
# disease stocks.
#
# VALIDATION + FINALISATION: separate structural / numerical / epidemiological-fit
# statuses. A candidate failing a STRUCTURAL check is NOT written as best_parameters
# (status FAILED_*, failed candidate saved, script stops; the last valid bundle is
# kept as calibrated_rhd_parameters_last_valid.rds). Fit outside the diagnostic
# tolerances is labelled in the status (never forced to look successful).
#
# OUTPUT: data/<COUNTRY>/calibrated_rhd_parameters.rds $stage_calibration
#   (+ calibration_structural_*.csv diagnostics). 05 consumes $stage_calibration.
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

SIGMA_TRANS <- getp("SIGMA_TRANS", 2.0)
LAMBDA_PRIOR_STRUCT <- getp("LAMBDA_PRIOR_STRUCT", 0.05)
PD_RATIO <- getp("STRUCT_PREV_DEATH_RATIO", 2.5)
W_TOTAL  <- getp("STRUCT_DATA_WEIGHT_TOTAL", 3.5)
W_PREV_STRUCT  <- getp("W_PREV_STRUCT",  W_TOTAL * PD_RATIO / (1 + PD_RATIO))
W_DEATH_STRUCT <- getp("W_DEATH_STRUCT", W_TOTAL / (1 + PD_RATIO))
W_INC_STRUCT   <- getp("W_INC_STRUCT", 0)            # 0 = incidence is a VALIDATION output
W_STAGE_STRUCT <- getp("W_STAGE_STRUCT", 1)          # weight on usable stage targets
CHILD_SHARE    <- getp("STRUCT_DEATH_CHILD_SHARE", 0.5)
PREV_CHILD_SHARE <- getp("STRUCT_PREV_CHILD_SHARE", 0.5)
CHILD_AGE_MAX  <- as.integer(getp("STRUCT_DEATH_CHILD_AGE_MAX", 19L))
CONV_ADULT_MIN <- as.integer(getp("STRUCT_DIAG_ADULT_CONVENTION_AGE_MIN", 21L))   # 21+ diagnostic (reported only)
for (v in c(CHILD_SHARE, PREV_CHILD_SHARE)) if (!is.finite(v) || v < 0 || v > 1)
  stop("04b: STRUCT_PREV_CHILD_SHARE / STRUCT_DEATH_CHILD_SHARE must be in [0, 1].", call. = FALSE)
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
TRANS_BOUNDS <- getp("TRANS_BOUNDS", list(
  p_A_to_no_rhd = c(0.005, 0.60), p_A_to_B = c(0.001, 0.25), p_B_to_A = c(0.002, 0.40),
  p_B_to_C = c(0.001, 0.20), p_C_to_B = c(0.001, 0.15), p_C_to_D = c(0.001, 0.30), p_D_to_C = c(0, 0)))
N_STARTS <- as.integer(getp("STRUCT_N_STARTS", 2L))
NEAROPT_DELTA <- getp("STRUCT_NEAROPT_DELTA", 0.05)
BURN_YEARS <- as.integer(getp("STRUCT_BURN_YEARS", 120L))
BURNIN_RECON <- getp("STRUCT_BURNIN_RECONCILIATION", "per_capita")
if (length(BURNIN_RECON) != 1L || !BURNIN_RECON %in% ABCD_BURNIN_RECONCILIATION_MODES)
  stop("04b: STRUCT_BURNIN_RECONCILIATION must be one of ", paste(ABCD_BURNIN_RECONCILIATION_MODES, collapse = ", "), ".", call. = FALSE)
BRIDGE <- isTRUE(getp("STRUCT_BRIDGE_TO_ANALYSIS_START", TRUE))
ANALYSIS_YEAR_START <- as.integer(getp("ANALYSIS_YEAR_START", 2026L))
DEATH_AGE_MIN <- as.integer(getp("STRUCT_DEATH_TARGET_AGE_MIN", 0L))

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

# spec consistency: 03's inputs were built for THIS engine and default specification
if (!identical(dmi$meta$engine_spec, ABCD_ENGINE_SPEC) || is.null(dmi$natural_history))
  stop("04b: disease_model_inputs.rds was not built for engine '", ABCD_ENGINE_SPEC, "' -- re-run 03.", call. = FALSE)
MORT_STRUCT <- dmi$meta$mortality_structure; SAP_MECH <- dmi$effects$sap_mechanism
if (!identical(MORT_STRUCT, "D_only") || !identical(SAP_MECH, "progression"))
  stop("04b: the default calibration requires D_only mortality and progression SAP.", call. = FALSE)
if (!identical(MORT_STRUCT, getp("mortality_structure", MORT_STRUCT)) || !identical(SAP_MECH, getp("sap_mechanism", SAP_MECH)))
  stop("04b: 00 specification differs from disease_model_inputs -- re-run 03.", call. = FALSE)
if (!identical(TRANS_BOUNDS$p_D_to_C, c(0, 0)) || dmi$transitions$p_D_to_C != 0)
  stop("04b: D->C is a fixed structural zero (TRANS_BOUNDS$p_D_to_C must be c(0, 0)).", call. = FALSE)

NH0      <- dmi$natural_history                         # fixed inputs + transition centres
tr_prior <- dmi$transitions                             # transition prior centres
eff_cal  <- dmi$effects; surg_cal <- dmi$surgery
cv <- dmi$coverage
# REFERENCE care applied in burn-in (static initial shares), window and bridge
# (dynamic): the reference-arm baseline values, i.e. exactly 06's reference arm.
REF_CARE <- c(lapply(cv$care_inputs, function(x) x$ref[["baseline"]]),
              list(surgery = surg_cal$coverage_ref_baseline, sensitivity = cv$screen_sensitivity,
                   confirm_completion = cv$confirm_completion, confirm_sensitivity = cv$confirm_sensitivity,
                   partial_relative_protection = cv$partial_adherence_relative_protection,
                   screen_age_lo = cv$screen_age_lo, screen_age_hi = cv$screen_age_hi))
REF_ADH  <- REF_CARE$adherence_good + REF_CARE$partial_relative_protection * REF_CARE$adherence_partial
INIT <- cv$init
CARE_BURN <- abcd_care_static(INIT$diagnosed_share, INIT$on_bpg_share, REF_ADH, REF_CARE$surgery)

if (!RUN_STRUCT) {
  message("04b: run_structural_calibration = FALSE -> no calibrated natural history is written.")
  calib$stage_calibration$status <- "structural_disabled_using_priors"
  calib$stage_calibration$best_parameters <- NULL
  saveRDS(calib, BUNDLE_FILE)
} else {

message(sprintf("── 04b structural calibration : %s | scored %d-%d (window %d-%d) | engine %s | D_only + progression SAP ──",
                LOCATION, MATCH_FROM, CAL_YEAR_END, CAL_YEAR_START, CAL_YEAR_END, ABCD_ENGINE_SPEC))

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
grp_lo_age  <- vapply(grp_present, function(g) min(ages[grp_lab_age == g]), numeric(1))
grp_hi_age  <- vapply(grp_present, function(g) max(ages[grp_lab_age == g]), numeric(1))
child_grp   <- grp_hi_age <= CHILD_AGE_MAX                  # GBD groups wholly within ages 0..CHILD_AGE_MAX
match_cols  <- yrs >= MATCH_FROM
stopifnot(sum(grp_lab_age == "95+ years") == sum(ages >= 95))   # 95+ compared with the model aggregated over 95..100
band_edges <- calib$meta$band_edges; band_labels <- calib$meta$bands
n_band <- length(band_edges); band_of <- function(a) findInterval(a, band_edges)
band_age <- band_of(ages)
CARE_WIN <- abcd_care_schedule(REF_CARE, n_yr, ages, SEXES)

# input arrays [age x sex x year]
arr0 <- function() array(0, dim = c(n_age, n_sex, n_yr), dimnames = list(ages, SEXES, yrs))
to_arr <- function(col) { a <- arr0()
  a[cbind(match(L2I$age, ages), match(L2I$sex, SEXES), match(L2I$year, yrs))] <- L2I[[col]]; a }
POP  <- to_arr("pop_wpp"); IRG <- to_arr("ir_gbd_total_pop"); OTH <- to_arr("bg_mx")
if (any(!is.finite(POP)) || any(POP < 0) || any(!is.finite(IRG)) || any(!is.finite(OTH)))
  stop("04b: non-finite/negative Layer-2 inputs.", call. = FALSE)

# target arrays [G x sex x year]: value used (NA = missing), mask, weights
tarr <- function(meas, col) { a <- array(NA_real_, dim = c(G, n_sex, n_yr), dimnames = list(grp_present, SEXES, yrs))
  d <- L2T[measure == meas]
  a[cbind(match(d$age_group, grp_present), match(d$sex, SEXES), match(d$year, yrs))] <- d[[col]]; a }
TGT <- list(prev = tarr("Prevalence", "value_used"), death = tarr("Deaths", "value_used"),
            inc = tarr("Incidence", "value_used"))
TGT_NUM <- list(prev = tarr("Prevalence", "number"), death = tarr("Deaths", "number"), inc = tarr("Incidence", "number"))
POP_G <- tarr("Prevalence", "pop_wpp")
year_w <- as.numeric(match_cols)                        # years before MATCH_FROM carry weight 0 (not scored)
age_w  <- rep(1, G)
if (!is.null(AGE_W)) {
  bl <- band_labels[band_of(grp_lo_age)]
  age_w <- ifelse(bl %in% names(AGE_W), as.numeric(AGE_W[bl]), 1)
}
mk_w <- function(meas) {
  w <- array(1, dim = c(G, n_sex, n_yr)) * age_w
  w <- sweep(w, 3, year_w, `*`)
  if (USE_IV_W) {
    lv <- tarr(meas, "logvar"); iv <- 1 / lv; iv[!is.finite(iv)] <- NA
    if (any(is.finite(iv))) { iv <- iv / mean(iv, na.rm = TRUE); iv[is.na(iv)] <- 1; iv <- pmin(iv, IV_CAP); w <- w * iv }
  }
  w
}
WT  <- list(prev = mk_w("Prevalence"), death = mk_w("Deaths"), inc = mk_w("Incidence"))
MSK <- lapply(TGT, function(a) !is.na(a) & a >= TARGET_MIN_COUNT)
MSK <- lapply(MSK, function(m) { m[is.na(m)] <- FALSE; m })
n_death_age_masked <- 0L
if (DEATH_AGE_MIN > 0) { drop_g <- grp_lo_age < DEATH_AGE_MIN
  n_death_age_masked <- sum(MSK$death[drop_g, , ] & WT$death[drop_g, , ] > 0); MSK$death[drop_g, , ] <- FALSE }
MSK$death_child <- MSK$death & array(child_grp, dim = dim(MSK$death))
MSK$death_adult <- MSK$death & array(!child_grp, dim = dim(MSK$death))
MSK$prev_child  <- MSK$prev & array(child_grp, dim = dim(MSK$prev))
MSK$prev_adult  <- MSK$prev & array(!child_grp, dim = dim(MSK$prev))
n_cells <- c(prev = sum(MSK$prev & WT$prev > 0), prev_child = sum(MSK$prev_child & WT$prev > 0),
             prev_adult = sum(MSK$prev_adult & WT$prev > 0), death_child = sum(MSK$death_child & WT$death > 0),
             death_adult = sum(MSK$death_adult & WT$death > 0), inc = sum(MSK$inc & WT$inc > 0))
if (n_cells[["prev_child"]] == 0L || n_cells[["prev_adult"]] == 0L)
  stop("04b: the prevalence loss needs usable childhood (0-19) AND adult cells (not replaced by a zero loss).", call. = FALSE)
# childhood cells below the count floor (reported; masked from the CELL loss only)
child_obs <- !is.na(TGT$death) & array(child_grp, dim = dim(TGT$death)) & array(rep(match_cols, each = G * n_sex), dim = dim(TGT$death))
n_child_floor_masked <- sum(child_obs & TGT$death < TARGET_MIN_COUNT)
if (n_cells[["death_child"]] == 0L || n_cells[["death_adult"]] == 0L)
  stop(sprintf("04b: the death loss needs usable childhood (0-%d) AND adult cells; got %d childhood / %d adult (age floor %d). Not replaced by a zero loss.",
               CHILD_AGE_MAX, n_cells[["death_child"]], n_cells[["death_adult"]], DEATH_AGE_MIN), call. = FALSE)
message(sprintf("  targets (%s basis): prevalence cells 0-%d %d / adult %d | deaths 0-%d: %d cells (%d observed childhood cells below count floor %g, reported) | deaths %d+: %d | incidence (validation) %d",
                calib$meta$layer2_target_basis, CHILD_AGE_MAX, n_cells[["prev_child"]], n_cells[["prev_adult"]], CHILD_AGE_MAX, n_cells[["death_child"]], n_child_floor_masked,
                TARGET_MIN_COUNT, CHILD_AGE_MAX + 1L, n_cells[["death_adult"]], n_cells[["inc"]]))

# stage targets (validated by 04; USED only here)
STG <- calib$stage_calibration$stage_targets
STG <- if (!is.null(STG)) as.data.table(STG)[reason == "usable"] else NULL
n_stage_used <- if (is.null(STG)) 0L else nrow(STG)

#===============================================================================
# 2. PARAMETER PACK / UNPACK: logit(6 transitions) + log(alpha band x sex)
#===============================================================================
trans_names <- c("p_A_to_no_rhd", "p_A_to_B", "p_B_to_A", "p_B_to_C", "p_C_to_B", "p_C_to_D")
n_tr <- length(trans_names); n_alpha <- n_band * n_sex
idx <- list(tr = seq_len(n_tr), alpha = n_tr + seq_len(n_alpha))
npar <- n_tr + n_alpha
alpha_names <- as.vector(outer(band_labels, SEXES, function(b, s) paste0("alpha[", b, "|", s, "]")))
param_labels <- c(trans_names, alpha_names)
tb <- function(nm) TRANS_BOUNDS[[nm]]
for (nm in trans_names) { b <- tb(nm)
  if (is.null(b) || length(b) != 2L || !(b[1] > 0 && b[2] < 1 && b[1] < b[2]))
    stop("04b: TRANS_BOUNDS$", nm, " must be 0 < lo < hi < 1 for a calibrated transition.", call. = FALSE)
  if (tr_prior[[nm]] < b[1] || tr_prior[[nm]] > b[2])
    stop(sprintf("04b: prior centre %s = %.4f lies outside its bounds [%.3f, %.3f].", nm, tr_prior[[nm]], b[1], b[2]), call. = FALSE) }
clamp_logit <- function(p, nm) qlogis(pmin(pmax(p, tb(nm)[1]), tb(nm)[2]))

unpack <- function(theta) {
  tr <- lapply(seq_len(n_tr), function(i) plogis(theta[idx$tr[i]])); names(tr) <- trans_names
  tr$p_D_to_C <- 0                                       # fixed structural zero (never transformed)
  alpha <- matrix(exp(theta[idx$alpha]), n_band, n_sex, dimnames = list(band_labels, SEXES))
  list(tr = tr, alpha = alpha)
}
nh_of <- function(tr) { nh <- NH0; nh$tr <- tr[c(trans_names, "p_D_to_C")]; nh$derived <- NULL; nh }
# incidence per TOTAL population: GBD rate x alpha(band, sex) (applied ONCE; the
# engine converts it to a first-onset probability per never-onset susceptible)
ir_total_arr <- function(alpha, irg = IRG) {
  a <- irg
  for (s in seq_len(n_sex)) a[, s, ] <- irg[, s, ] * alpha[band_age, s]
  a
}

lo <- numeric(npar); hi <- numeric(npar)
for (i in seq_len(n_tr)) { b <- tb(trans_names[i]); lo[idx$tr[i]] <- qlogis(b[1]); hi[idx$tr[i]] <- qlogis(b[2]) }
lo[idx$alpha] <- log(ALPHA_LO); hi[idx$alpha] <- log(ALPHA_HI)
theta0 <- numeric(npar)
for (i in seq_len(n_tr)) theta0[idx$tr[i]] <- clamp_logit(tr_prior[[trans_names[i]]], trans_names[i])
theta0[idx$alpha] <- 0
theta0 <- pmin(pmax(theta0, lo), hi)
alpha_L1 <- matrix(1, n_band, n_sex, dimnames = list(band_labels, SEXES))
fac <- as.data.table(calib$factors)[location == LOCATION]
for (s in SEXES) { f <- fac[sex == s]; if (nrow(f)) alpha_L1[match(f$band, band_labels), s] <- f$alpha }
theta_L1 <- theta0; theta_L1[idx$alpha] <- pmin(pmax(log(as.vector(alpha_L1)), lo[idx$alpha]), hi[idx$alpha])

#===============================================================================
# 3. ENGINE PROJECTION (shared wrappers) + OBSERVATIONS + OBJECTIVE
#===============================================================================
agg_grp <- function(M) rowsum(M, grp_int, reorder = FALSE)   # [age x sex] -> [G x sex]
obs_conv <- function(st, en) switch(PREV_OBS, start = st, end = en, mid = 0.5 * (st + en))

conv_age <- ages >= CONV_ADULT_MIN                          # 21+ diagnostic (single-year model ages)
project_window <- function(p, burn_years = BURN_YEARS, keep_states = FALSE) {
  nhe <- abcd_nh_build(nh_of(p$tr), ages, SEXES)
  IR  <- ir_total_arr(p$alpha)
  bi  <- abcd_burn_in(POP[, , 1], IR[, , 1], OTH[, , 1], nhe, eff_cal, surg_cal, CARE_BURN, burn_years, reconcile = BURNIN_RECON)
  st0 <- abcd_care_initialize(bi$state, INIT$diagnosed_share, INIT$on_bpg_share, INIT$ever_screened_share, nhe)
  # light engine inside the objective; the FULL engine (same code path) for the traced optimum
  pr  <- abcd_project(st0, POP, IR, OTH, nhe, eff_cal, surg_cal, CARE_WIN, yrs, light = !keep_states, prev_obs = PREV_OBS)
  prev <- array(0, dim = c(G, n_sex, n_yr)); death <- prev; inc <- prev
  conv21 <- array(0, dim = c(2L, n_sex, n_yr), dimnames = list(c("prev", "death"), SEXES, yrs))
  need_stg <- keep_states || n_stage_used > 0L
  stg <- if (need_stg) lapply(c(A = "A", B = "B", C = "C", D = "D"), function(s) vector("list", n_yr)) else NULL
  for (t in seq_len(n_yr)) {
    cyc <- pr$cycles[[t]]
    prev[, , t]  <- agg_grp(pr$prev_obs[[t]])
    death[, , t] <- agg_grp(cyc$rhd_deaths)
    inc[, , t]   <- agg_grp(cyc$new_rhd_total)
    conv21["prev", , t]  <- colSums(pr$prev_obs[[t]][conv_age, , drop = FALSE])
    conv21["death", , t] <- colSums(cyc$rhd_deaths[conv_age, , drop = FALSE])
    if (need_stg) { sa <- abcd_stage_aggregate(pr$start[[t]]); ea <- abcd_stage_aggregate(cyc$state)
      for (s in c("A", "B", "C", "D")) stg[[s]][[t]] <- obs_conv(sa[[s]], ea[[s]]) }
  }
  out <- list(prev = prev, death = death, inc = inc, conv21 = conv21, stages = stg,
              final_state = pr$final_state, burn = bi, nhe = nhe,
              pop_deficit = pr$pop_deficit_total, pop_deficit_cells = pr$pop_deficit_cells, care_violation = pr$care_violation_total,
              onset_cap_hits = pr$onset_cap_hits, demographic_residual = pr$demographic_residual)
  if (keep_states) { out$start <- pr$start; out$cycles <- pr$cycles; out$boundary_transfer <- pr$boundary_transfer
    out$reclass <- pr$reclass; out$prev_obs_cells <- pr$prev_obs }
  out
}

loss_measure <- function(model, tgt, msk, wt) {
  m <- msk & wt > 0
  if (!any(m)) return(0)
  r <- (log(model[m] + EPS_LOG) - log(tgt[m] + EPS_LOG))^2 * wt[m]
  sum(r) / sum(wt[m])
}
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
prior_terms <- function(p) {
  pri_tr <- sum(vapply(trans_names, function(nm) (qlogis(p$tr[[nm]]) - clamp_logit(tr_prior[[nm]], nm))^2, numeric(1))) / SIGMA_TRANS^2
  la <- log(p$alpha)
  pri_alpha <- sum((la / SIGMA_ALPHA_STRUCT)^2)
  smo_alpha <- if (n_band >= 3) sum(apply(la, 2, function(v) sum(diff(v, differences = 2)^2))) else 0
  c(prior_transitions = LAMBDA_PRIOR_STRUCT * pri_tr, prior_incidence_alpha = LAMBDA_PRIOR_STRUCT * pri_alpha,
    smooth_incidence_alpha = LAMBDA_SMOOTH_STRUCT * smo_alpha)
}
decompose <- function(theta, wP = W_PREV_STRUCT, wD = W_DEATH_STRUCT, wI = W_INC_STRUCT, wS = W_STAGE_STRUCT,
                      burn_years = BURN_YEARS) {
  p <- unpack(theta); pw <- project_window(p, burn_years = burn_years)
  LPc <- loss_measure(pw$prev,  TGT$prev,  MSK$prev_child, WT$prev)
  LPa <- loss_measure(pw$prev,  TGT$prev,  MSK$prev_adult, WT$prev)
  LP  <- PREV_CHILD_SHARE * LPc + (1 - PREV_CHILD_SHARE) * LPa
  LDc <- loss_measure(pw$death, TGT$death, MSK$death_child, WT$death)
  LDa <- loss_measure(pw$death, TGT$death, MSK$death_adult, WT$death)
  LD  <- CHILD_SHARE * LDc + (1 - CHILD_SHARE) * LDa
  LI  <- loss_measure(pw$inc,   TGT$inc,   MSK$inc,   WT$inc)
  LS  <- stage_loss(pw)
  pri <- prior_terms(p)
  terms <- c(prevalence_0_19_unweighted = LPc, prevalence_20plus_unweighted = LPa, prevalence_unweighted = LP,
             deaths_0_19_unweighted = LDc, deaths_20plus_unweighted = LDa,
             deaths_unweighted = LD, incidence_unweighted = LI, stage_unweighted = LS,
             prevalence_0_19_weighted = wP * PREV_CHILD_SHARE * LPc, prevalence_20plus_weighted = wP * (1 - PREV_CHILD_SHARE) * LPa,
             deaths_0_19_weighted = wD * CHILD_SHARE * LDc, deaths_20plus_weighted = wD * (1 - CHILD_SHARE) * LDa,
             prevalence_weighted = wP * LP, deaths_weighted = wD * LD, incidence_weighted = wI * LI,
             stage_weighted = wS * LS, pri)
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
base_dec <- decompose(theta0)                        # prior-centre baseline (fails loudly if the engine errors)
base_val <- base_dec$total
starts <- list(prior_center = theta0, layer1_alpha = theta_L1)
if (PREV_FIRST) starts$prevalence_first <- theta0    # prevalence-only pre-fit from the prior centre, then the joint fit
set.seed(SEED)
for (i in seq_len(N_STARTS)) starts[[paste0("random_", i)]] <- runif(npar, lo, hi)
run_start <- function(nm) {
  if (nm == "prevalence_first") {
    op <- run_lbfgs(theta0, objective_prev, maxit = 80L)
    if (is.null(op)) return(NULL)
    r <- run_lbfgs(op$par); if (!is.null(r)) r$start_par <- op$par
    return(r)
  }
  run_lbfgs(starts[[nm]])
}
# independent deterministic starts -> optionally run on parallel workers (identical results)
PAR_STARTS <- isTRUE(getp("STRUCT_PARALLEL_STARTS", TRUE)) && length(starts) > 1L
n_workers <- min(length(starts), as.integer(getp("MAX_CORES", 14L)), max(1L, parallel::detectCores() - 1L))
if (PAR_STARTS && n_workers > 1L) {
  message(sprintf("  running %d optimiser starts on %d parallel workers", length(starts), n_workers))
  cl <- parallel::makeCluster(n_workers)
  fits <- tryCatch({
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)); NULL })
    parallel::clusterExport(cl, ls(globalenv(), all.names = TRUE), envir = globalenv())
    parallel::parLapply(cl, names(starts), run_start)
  }, finally = parallel::stopCluster(cl))
} else fits <- lapply(names(starts), run_start)
names(fits) <- names(starts)
for (nm in names(fits)) if (!is.null(fits[[nm]]$start_par)) starts[[nm]] <- fits[[nm]]$start_par
ok_fit <- !vapply(fits, is.null, logical(1))
if (!any(ok_fit)) stop("04b: every optimiser start failed.", call. = FALSE)
fits <- fits[ok_fit]; starts <- starts[ok_fit]
vals <- vapply(fits, function(f) f$value, numeric(1))
conv <- data.table(start = names(fits), objective = vals,
                   convergence_code = vapply(fits, function(f) f$convergence, numeric(1)),
                   message = vapply(fits, function(f) if (is.null(f$message)) "" else f$message, character(1)),
                   n_fn_evals = vapply(fits, function(f) f$counts[["function"]], numeric(1)))
best <- fits[[which.min(vals)]]; best_start <- names(fits)[which.min(vals)]
if (base_val < best$value) { best <- list(par = theta0, value = base_val, convergence = 0L); best_start <- "prior_center_unoptimised" }
theta_hat <- best$par; L_min <- best$value
dec <- decompose(theta_hat); bp <- dec$p
pw <- project_window(bp, keep_states = TRUE)           # full trace at the optimum (light engine)
elapsed_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("  optimum L = %.4f (prior-centre baseline L = %.4f) | starts %d | spread %.4f/%.4f/%.4f | %.1f min\n",
            L_min, base_val, length(fits), min(vals), stats::median(vals), max(vals), elapsed_min))
cat("  loss decomposition at the optimum:\n"); print(round(dec$terms, 4))
pri_names <- c("prior_transitions", "prior_incidence_alpha", "smooth_incidence_alpha")
prior_over_data <- sum(dec$terms[pri_names]) /
  max(sum(dec$terms[c("prevalence_weighted", "deaths_weighted", "incidence_weighted", "stage_weighted")]), 1e-12)
cat(sprintf("  summed prior+smoothness penalties / weighted data terms at optimum = %.3f (>1 => priors dominate)\n", prior_over_data))

#===============================================================================
# 5. FIT TABLES (all ages / childhood / adult; Rate x WPP target and GBD Number)
#===============================================================================
ly <- n_yr
fit_long <- rbindlist(lapply(c("prev", "death", "inc"), function(m) {
  meas <- c(prev = "Prevalence", death = "Deaths", inc = "Incidence")[[m]]
  idxg <- expand.grid(g = seq_len(G), s = seq_len(n_sex), t = seq_len(n_yr))
  in_l <- if (m == "death") (MSK$death_child | MSK$death_adult) & WT$death > 0 else
          if (m == "inc") MSK$inc & WT$inc > 0 & W_INC_STRUCT > 0 else MSK$prev & WT$prev > 0
  data.table(location = LOCATION, measure = meas, age_group = grp_present[idxg$g], sex = SEXES[idxg$s], year = yrs[idxg$t],
             age_component = fifelse(child_grp[idxg$g], sprintf("0-%d", CHILD_AGE_MAX), sprintf("%d+", CHILD_AGE_MAX + 1L)),
             model = pw[[m]][cbind(idxg$g, idxg$s, idxg$t)],
             target_used = TGT[[m]][cbind(idxg$g, idxg$s, idxg$t)],
             target_number = TGT_NUM[[m]][cbind(idxg$g, idxg$s, idxg$t)],
             pop_wpp = POP_G[cbind(idxg$g, idxg$s, idxg$t)],
             in_loss = in_l[cbind(idxg$g, idxg$s, idxg$t)],
             weight = WT[[m]][cbind(idxg$g, idxg$s, idxg$t)])
}))
fit_long[, `:=`(model_rate_per100k = fifelse(pop_wpp > 0, model / pop_wpp * 1e5, NA_real_),
                target_rate_per100k = fifelse(pop_wpp > 0, target_used / pop_wpp * 1e5, NA_real_),
                log_ratio = fifelse(!is.na(target_used) & target_used > 0, log((model + EPS_LOG) / (target_used + EPS_LOG)), NA_real_),
                abs_error = fifelse(!is.na(target_used), model - target_used, NA_real_),
                target_minus_number = target_used - target_number)]
fit_long[, scored := year >= MATCH_FROM]
fit_long[, held_out := !scored]

# aggregated fit: counts aggregated FIRST, then ratios (zero target -> ratio NA, abs error kept)
agg_fit <- function(d, by) {
  a <- d[, .(model = sum(model), target = sum(target_used, na.rm = TRUE), gbd_number = sum(target_number, na.rm = TRUE),
             n_cells = .N, n_cells_observed = sum(!is.na(target_used)), n_cells_in_loss = sum(in_loss),
             n_cells_below_floor = sum(!is.na(target_used) & target_used < TARGET_MIN_COUNT)), by = by]
  a[, `:=`(ratio = fifelse(target > 0, model / target, NA_real_), abs_error = model - target,
           log_ratio = fifelse(target > 0 & model > 0, log(model / target), NA_real_),
           ratio_vs_gbd_number = fifelse(gbd_number > 0, model / gbd_number, NA_real_))]
  a[, abs_log_ratio := abs(log_ratio)]
  a[]
}
fl_d <- fit_long[measure == "Deaths"]
death_fit <- rbindlist(list(
  agg_fit(fl_d[age_component == sprintf("0-%d", CHILD_AGE_MAX)], c("year", "sex", "age_group"))[, level := "childhood GBD group"],
  agg_fit(fl_d, c("year", "sex", "age_component"))[, level := "age component"],
  agg_fit(copy(fl_d)[, age_component := "All ages"], c("year", "sex", "age_component"))[, level := "all ages"],
  agg_fit(copy(fl_d)[, sex := "Both"], c("year", "sex", "age_component"))[, level := "age component"],
  agg_fit(copy(fl_d)[, `:=`(sex = "Both", age_component = "All ages")], c("year", "sex", "age_component"))[, level := "all ages"],
  agg_fit(copy(fl_d)[age_component == sprintf("0-%d", CHILD_AGE_MAX)][, sex := "Both"], c("year", "sex", "age_group"))[, level := "childhood GBD group"]),
  fill = TRUE, use.names = TRUE)
death_fit[is.na(age_component), age_component := age_group]
death_fit[, scored := year >= MATCH_FROM]
death_fit[, within_tol := fifelse(is.na(abs_log_ratio), NA, abs_log_ratio < FIT_TOL_DEATH)]
child_share <- merge(death_fit[sex == "Both" & age_component == sprintf("0-%d", CHILD_AGE_MAX), .(year, child_model = model, child_target = target, child_number = gbd_number)],
                     death_fit[sex == "Both" & age_component == "All ages", .(year, all_model = model, all_target = target, all_number = gbd_number)], by = "year")
child_share[, `:=`(child_share_model = child_model / all_model, child_share_target = child_target / all_target,
                   child_share_gbd_number = child_number / all_number)]
death_fit <- merge(death_fit, child_share[, .(year, child_share_model, child_share_target, child_share_gbd_number)], by = "year", all.x = TRUE)
setorder(death_fit, level, sex, age_component, age_group, year)

fl_p <- fit_long[measure == "Prevalence"]
# 21+ DIAGNOSTIC (isolates the <= 20 / > 20 modelling convention from the GBD 0-19 / 20+
# split): model = single ages >= CONV_ADULT_MIN; target = 20+ target minus the age-20 share
# of the 20-24 group (GBD 20-24 rate x WPP age-20 population; GBD Number split by the same
# population share -- an approximation, labelled).
g2024 <- which(grp_present == "20-24 years"); a_lo <- which(ages == 20L):which(ages == CONV_ADULT_MIN - 1L)
conv21_tab <- rbindlist(lapply(c(prev = "Prevalence", death = "Deaths"), function(meas) {
  m <- c(Prevalence = "prev", Deaths = "death")[[meas]]
  rbindlist(lapply(seq_len(n_yr), function(t) {
    share <- colSums(matrix(POP[a_lo, , t], nrow = length(a_lo))) / POP_G[g2024, , t]
    tg <- sum(TGT[[m]][!child_grp, , t], na.rm = TRUE) - sum(TGT[[m]][g2024, , t] * share, na.rm = TRUE)
    gn <- sum(TGT_NUM[[m]][!child_grp, , t], na.rm = TRUE) - sum(TGT_NUM[[m]][g2024, , t] * share, na.rm = TRUE)
    data.table(measure = sprintf("%s %d+ (diagnostic)", meas, CONV_ADULT_MIN), year = yrs[t],
               model = sum(pw$conv21[m, , t]), target = tg, gbd_number = gn)
  }))
}))
conv21_tab[, `:=`(ratio = model / target, abs_error = model - target, log_ratio = log(model / target),
                  ratio_vs_gbd_number = model / gbd_number)]
conv21_tab[, abs_log_ratio := abs(log_ratio)]
fit_year <- rbindlist(list(
  agg_fit(fit_long, c("measure", "year")),
  agg_fit(fl_d, c("measure", "year", "age_component"))[, measure := paste0("Deaths ", age_component)][, age_component := NULL],
  agg_fit(fl_p, c("measure", "year", "age_component"))[, measure := paste0("Prevalence ", age_component)][, age_component := NULL],
  conv21_tab),
  fill = TRUE, use.names = TRUE)
fit_year[, `:=`(scored = year >= MATCH_FROM, held_out = year < MATCH_FROM)]
fit_age_last <- agg_fit(fit_long[year == CAL_YEAR_END], c("measure", "age_group"))
prev_fit_age_sex <- agg_fit(fit_long[measure == "Prevalence"], c("year", "sex", "age_group"))
fit_sex_component <- rbindlist(list(
  agg_fit(fl_p, c("year", "sex", "age_component"))[, measure := "Prevalence"],
  agg_fit(copy(fl_p)[, age_component := "All ages"], c("year", "sex", "age_component"))[, measure := "Prevalence"],
  agg_fit(fl_d, c("year", "sex", "age_component"))[, measure := "Deaths"],
  agg_fit(copy(fl_d)[, age_component := "All ages"], c("year", "sex", "age_component"))[, measure := "Deaths"]), use.names = TRUE)

# tolerance status across ALL scored years (not one favourable endpoint). role:
# "checkpoint" quantities decide fit adequacy; "diagnostic" rows are reported only.
cA <- sprintf("0-%d", CHILD_AGE_MAX); aA <- sprintf("%d+", CHILD_AGE_MAX + 1L)
tq <- function(meas, q, tol, role) fit_year[measure == meas, .(quantity = q, year, ratio, abs_log_ratio, tol = tol, role = role)]
tol_tab <- rbindlist(list(
  tq("Prevalence", "prevalence, all ages", FIT_TOL_PREV, "checkpoint"),
  tq(paste("Prevalence", aA), sprintf("prevalence, ages %s", aA), FIT_TOL_PREV, "checkpoint"),
  tq("Deaths", "RHD deaths, all ages", FIT_TOL_DEATH, "checkpoint"),
  tq(paste("Deaths", cA), sprintf("RHD deaths, ages %s", cA), FIT_TOL_DEATH, "checkpoint"),
  tq(paste("Deaths", aA), sprintf("RHD deaths, ages %s", aA), FIT_TOL_DEATH, "checkpoint"),
  tq(paste("Prevalence", cA), sprintf("prevalence, ages %s", cA), FIT_TOL_PREV, "diagnostic"),
  tq(sprintf("Prevalence %d+ (diagnostic)", CONV_ADULT_MIN), sprintf("prevalence, ages %d+ (convention diagnostic)", CONV_ADULT_MIN), FIT_TOL_PREV, "diagnostic"),
  tq(sprintf("Deaths %d+ (diagnostic)", CONV_ADULT_MIN), sprintf("RHD deaths, ages %d+ (convention diagnostic)", CONV_ADULT_MIN), FIT_TOL_DEATH, "diagnostic"),
  tq("Incidence", "incidence (validation)", NA_real_, "validation")))
tol_tab[, within := fifelse(is.na(tol) | is.na(abs_log_ratio), NA, abs_log_ratio < tol)]
tol_summary <- tol_tab[year >= MATCH_FROM, .(role = role[1], ratio_min = min(ratio, na.rm = TRUE), ratio_max = max(ratio, na.rm = TRUE),
                                            ratio_last_year = ratio[year == CAL_YEAR_END],
                                            max_abs_log_ratio = max(abs_log_ratio, na.rm = TRUE), tol = tol[1],
                                            years_within = sum(within, na.rm = TRUE), years_scored = .N), by = quantity]
tol_summary[, status := fifelse(is.na(tol), "validation only",
                                fifelse(years_within == years_scored, "within tolerance in every scored year",
                                        sprintf("OUTSIDE tolerance in %d of %d scored years", years_scored - years_within, years_scored)))]
gtot <- function(q) tol_tab[quantity == q & year == CAL_YEAR_END]
prev_ratio <- gtot("prevalence, all ages")$ratio; death_ratio <- gtot("RHD deaths, all ages")$ratio
inc_ratio  <- gtot("incidence (validation)")$ratio
child_ratio <- gtot(sprintf("RHD deaths, ages %s", cA))$ratio
adult_ratio <- gtot(sprintf("RHD deaths, ages %s", aA))$ratio
fit_ok <- all(tol_summary[role == "checkpoint", years_within == years_scored])
# age-tail sensitivity: loss excluding the 85+ groups
tail_grp <- grp_present %in% c("85-89 years", "90-94 years", "95+ years")
msk_notail <- lapply(MSK, function(m) { m2 <- m; m2[tail_grp, , ] <- FALSE; m2 })
loss_notail <- c(prevalence_0_19 = loss_measure(pw$prev, TGT$prev, msk_notail$prev_child, WT$prev),
                 prevalence_20plus = loss_measure(pw$prev, TGT$prev, msk_notail$prev_adult, WT$prev),
                 deaths_0_19 = loss_measure(pw$death, TGT$death, msk_notail$death_child, WT$death),
                 deaths_20plus = loss_measure(pw$death, TGT$death, msk_notail$death_adult, WT$death))

#===============================================================================
# 6. STOCKS, FLOWS, ACCOUNTING, D FATALITY, STAGE AT FIRST DIAGNOSIS
#    (traced optimum: FULL engine output, same code path as the light objective)
#===============================================================================
nhe_best <- pw$nhe
b_age <- nhe_best$boundary_age; ch_a <- ages <= b_age; ad_a <- !ch_a
w_lo <- cv$screen_age_lo; w_hi <- cv$screen_age_hi; in_w <- ages >= w_lo & ages <= w_hi
stage_stocks_year <- rbindlist(lapply(seq_len(n_yr), function(t) {
  st <- pw$start[[t]]; cyc <- pw$cycles[[t]]; rc <- pw$reclass[[t]]; ag <- abcd_stage_aggregate(st)
  living_start <- sum(cyc$rhd_start + cyc$prior_A_start + cyc$never_onset_start + cyc$recovered_start)
  data.table(year = yrs[t],
             A = sum(pw$stages$A[[t]]), B = sum(pw$stages$B[[t]]), C = sum(pw$stages$C[[t]]), D = sum(pw$stages$D[[t]]),   # observed (mid)
             A_start = sum(ag$A), B_start = sum(ag$B), C_start = sum(ag$C), D_start = sum(ag$D),
             D_start_adult = sum(ag$D[ad_a, ]), D_start_child = sum(ag$D[ch_a, ]),
             undiagnosed_start = sum(abcd_undiagnosed_total(st)), diagnosed_start = sum(abcd_diagnosed_total(st)),
             on_bpg_start = sum(st$T_A + st$T_B), prior_A_start = sum(abcd_prior_A_total(st)),
             normalised_start = sum(abcd_normalised_total(st)), never_onset_start = sum(abcd_never_onset_total(st)),
             A_adult_start = sum(ag$A[ad_a, ]), T_A_adult_start = sum(st$T_A[ad_a, ]),
             P_child_start = sum(abcd_prior_A_total(st)[ch_a, ]),
             rhd_deaths = sum(cyc$rhd_deaths), rhd_deaths_D = sum(cyc$rhd_deaths_D),
             rhd_deaths_D_undiagnosed = sum(cyc$rhd_deaths_D_undiagnosed), rhd_deaths_D_diagnosed = sum(cyc$rhd_deaths_D_diagnosed),
             rhd_deaths_adult = sum(cyc$rhd_deaths[ad_a, ]),
             new_disease = sum(cyc$new_rhd_total), new_disease_child = sum(cyc$new_rhd_total[ch_a, ]),
             new_disease_C_or_D = sum(cyc$new_rhd_C + cyc$new_rhd_D),
             reentry_prior_A_to_B = sum(cyc$reentry_prior_A_to_B),
             reclass_A_to_prior_A = sum(rc$reclass_undiagnosed + rc$reclass_diagnosed), reclass_on_bpg_A = sum(rc$reclass_on_bpg),
             B_to_C = sum(cyc$B_to_C), B_to_C_adult = sum(cyc$B_to_C[ad_a, ]),
             C_to_D = sum(cyc$C_to_D), C_to_D_adult = sum(cyc$C_to_D[ad_a, ]),
             A_to_normalised = sum(cyc$A_to_no_rhd), B_to_normalised = sum(cyc$B_to_no_rhd),
             other_deaths_rhd = sum(cyc$other_deaths),
             first_dx = sum(cyc$n_new_diagnoses), first_dx_screen = sum(cyc$n_dx_screen), first_dx_clinical = sum(cyc$n_dx_clinical),
             first_dx_A = sum(cyc$n_dx_A), first_dx_B = sum(cyc$n_dx_B), first_dx_C = sum(cyc$n_dx_C), first_dx_D = sum(cyc$n_dx_D),
             screen_encounters = sum(cyc$n_screened), screen_true_positive = sum(cyc$n_screen_true_positive),
             confirm_completed = sum(cyc$n_confirm_completed), on_bpg_person_years = sum(cyc$n_on_optimal_treatment),
             screen_outside_window = sum(cyc$n_screened[!in_w, ]),
             ever_screened_5_39 = 1 - sum((cyc$never_screened_share * POP[, , t])[in_w, ]) / sum(POP[in_w, , t]),
             # accounting: start (after backbone reconciliation) - all deaths = end (pre-ageing); classifiable identity
             mass_balance_gap = living_start - sum(cyc$all_mx) - sum(abcd_living_total(cyc$state)),
             classifiable_identity_gap = sum(cyc$sick_end) - (sum(cyc$rhd_start) + sum(cyc$new_rhd_total) + sum(cyc$reentry_prior_A_to_B) -
                                           sum(cyc$A_to_no_rhd) - sum(cyc$B_to_no_rhd) - sum(cyc$rhd_deaths) - sum(cyc$other_deaths)),
             # nesting: diagnosed <= classifiable; on BPG <= eligible diagnosed; no adult A on BPG;
             # first diagnoses <= undiagnosed disease entering the cycle (each person diagnosed once)
             dx_nesting_violation = max(0, max(cyc$n_diagnosed - cyc$rhd_start), max(cyc$n_on_optimal_treatment - cyc$n_bpg_eligible_diagnosed),
                                       max(cyc$n_on_bpg_A[ad_a, ])),
             first_dx_exceeds_undiagnosed = max(0, max(cyc$n_new_diagnoses - abcd_undiagnosed_total(st))))
}))
stage_stocks_year[, `:=`(living_rhd_start = A_start + B_start + C_start + D_start,
                         deaths_outside_D = rhd_deaths - rhd_deaths_D)]
stage_stocks_year[, `:=`(share_D = D / (A + B + C + D), diagnosed_share = diagnosed_start / living_rhd_start,
                         on_bpg_share = on_bpg_start / living_rhd_start)]
# stage stocks (observed convention) + flows by AGE COMPONENT, every window year: GBD split
# (0-19 / 20+) and the engine convention (0-20 / 21+) -- child vs adult D share, deaths per case
stage_by_age_component <- rbindlist(lapply(seq_len(n_yr), function(t) {
  cyc <- pw$cycles[[t]]; st <- pw$start[[t]]
  rbindlist(lapply(setNames(list(ages <= CHILD_AGE_MAX, ages > CHILD_AGE_MAX, ch_a, ad_a, rep(TRUE, n_age)),
                            c(sprintf("0-%d", CHILD_AGE_MAX), sprintf("%d+", CHILD_AGE_MAX + 1L),
                              sprintf("0-%d", b_age), sprintf("%d+", b_age + 1L), "All ages")),
    function(sel) data.table(year = yrs[t],
      A = sum(pw$stages$A[[t]][sel, ]), B = sum(pw$stages$B[[t]][sel, ]), C = sum(pw$stages$C[[t]][sel, ]), D = sum(pw$stages$D[[t]][sel, ]),
      prev_obs = sum(pw$prev_obs_cells[[t]][sel, ]), D_start = sum(cyc$D_start[sel, ]),
      rhd_deaths = sum(cyc$rhd_deaths[sel, ]), prior_A_start = sum(abcd_prior_A_total(st)[sel, ]),
      reentry_prior_A_to_B = sum(cyc$reentry_prior_A_to_B[sel, ]), new_disease = sum(cyc$new_rhd_total[sel, ]),
      B_to_C = sum(cyc$B_to_C[sel, ]), C_to_D = sum(cyc$C_to_D[sel, ]), B_to_normalised = sum(cyc$B_to_no_rhd[sel, ]))),
    idcol = "age_component")
}))
stage_by_age_component[, `:=`(share_D = D / (A + B + C + D), deaths_per_prev_case = rhd_deaths / prev_obs)]
d_fat <- rbindlist(lapply(seq_len(n_yr), function(t) {
  cyc <- pw$cycles[[t]]
  rbindlist(lapply(list(all = rep(TRUE, n_age), "0-20" = ch_a, "21+" = ad_a), function(sel) {
    Ds <- sum(cyc$D_start[sel, ]); dd <- sum(cyc$rhd_deaths_D[sel, ])
    data.table(year = yrs[t], D_start = Ds, rhd_deaths_D = dd, input_D_death_prob = nhe_best$pd,
               realized_D_death_risk = if (Ds > 0) dd / Ds else NA_real_,
               D_start_zero = Ds <= 0, zero_deaths_when_undefined_ok = if (Ds > 0) NA else dd == 0)
  }), idcol = "ages")
}))
# stage at FIRST CONFIRMED DIAGNOSIS by age at diagnosis (scored years) vs the PROVISIONAL benchmark
dx_comp <- rbindlist(lapply(list("<=20" = ch_a, ">20" = ad_a), function(sel) {
  v <- vapply(c("A", "B", "C", "D"), function(s) sum(vapply(which(match_cols), function(t) sum(pw$cycles[[t]][[paste0("n_dx_", s)]][sel, ]), 0)), 0)
  data.table(stage = c("A", "B", "C", "D"), first_dx = v, model_share = v / sum(v))
}), idcol = "age_at_diagnosis")
dx_comp[, benchmark_share := c(dmi$dx_mix_benchmark$child, dmi$dx_mix_benchmark$adult)[.I]]
dx_comp[, note := "benchmark = PROVISIONAL diagnosis-composition reference (not a fit term or national stage target)"]
fs <- pw$start[[n_yr]]                                     # stocks ENTERING the final window year
fsa <- abcd_stage_aggregate(fs)
tot_stock <- vapply(fsa, sum, numeric(1)); stage_shares <- tot_stock / sum(tot_stock)
care_stage_shares <- vapply(ABCD_DISEASE_STOCKS, function(k) sum(fs[[k]]), numeric(1)) / sum(tot_stock)
stage_shares_by_age <- rbindlist(lapply(seq_len(n_sex), function(si) {
  tot_s <- fsa$A[, si] + fsa$B[, si] + fsa$C[, si] + fsa$D[, si]
  data.table(sex = SEXES[si], age = ages,
             fracA = ifelse(tot_s > 0, fsa$A[, si] / tot_s, 0), fracB = ifelse(tot_s > 0, fsa$B[, si] / tot_s, 0),
             fracC = ifelse(tot_s > 0, fsa$C[, si] / tot_s, 0), fracD = ifelse(tot_s > 0, fsa$D[, si] / tot_s, 0),
             diagnosed_share = ifelse(tot_s > 0, abcd_diagnosed_total(fs)[, si] / tot_s, NA_real_),
             prior_A = abcd_prior_A_total(fs)[, si], normalised = abcd_normalised_total(fs)[, si])
}))
model_prev_tot <- sum(pw$prev[, , ly]); model_death_tot <- sum(pw$death[, , ly])
gbd_prev_tot <- sum(TGT$prev[, , ly], na.rm = TRUE); gbd_death_tot <- sum(TGT$death[, , ly], na.rm = TRUE)
cf_implied <- model_death_tot / model_prev_tot; gbd_cf <- gbd_death_tot / gbd_prev_tot

attribution <- if (fit_ok)
  sprintf("prevalence (all ages, %s) and RHD deaths (all ages, %s and %s) within the checkpoint tolerances in every scored year", aA, cA, aA) else
  paste0("OUTSIDE checkpoint tolerance: ", paste(tol_summary[role == "checkpoint" & years_within < years_scored,
          sprintf("%s (ratio %.2f-%.2f over %d-%d)", quantity, ratio_min, ratio_max, MATCH_FROM, CAL_YEAR_END)], collapse = "; "),
         sprintf(". Fixed inputs that bound the fit: D death %.4f (all histories), entry mix child %s / adult %s, adult B normalisation %.3f; D share entering %d = %.1f%% (GBD HF proxy %s).",
                 nhe_best$pd, paste(NH0$entry_mix_child, collapse = "/"), paste(NH0$entry_mix_adult, collapse = "/"),
                 NH0$p_B_to_no_rhd_adult, CAL_YEAR_END, 100 * stage_shares[["D"]],
                 if (is.na(dmi$meta$d_share_proxy_gbd_hf)) "n/a" else sprintf("%.1f%%", 100 * dmi$meta$d_share_proxy_gbd_hf)))

#===============================================================================
# 7. PARTIAL IDENTIFICATION — near-optimal set + one-at-a-time flatness sweep;
#    fitted vs fixed parameter table with boundary hits
#===============================================================================
acc <- which(vals <= L_min * (1 + NEAROPT_DELTA))
accepted_full <- rbindlist(lapply(acc, function(k) {
  pk <- unpack(fits[[k]]$par)
  cbind(data.table(start = names(fits)[k], L = vals[k], convergence_code = fits[[k]]$convergence),
        as.data.table(as.list(unlist(pk$tr[trans_names]))),
        as.data.table(as.list(setNames(as.vector(pk$alpha), alpha_names))))
}), fill = TRUE)
accepted <- accepted_full[, c("start", "L", trans_names), with = FALSE]
sweep_frac <- function(j, npts = 11) {
  grid <- seq(lo[j], hi[j], length.out = npts)
  Ls <- vapply(grid, function(v) { th <- theta_hat; th[j] <- v; objective(th) }, numeric(1))
  mean(Ls <= L_min * (1 + NEAROPT_DELTA))
}
identifiability <- data.table(parameter = param_labels, flatness = vapply(seq_len(npar), sweep_frac, numeric(1)))
identifiability[, status := fifelse(flatness < 0.5, "data-identified", "prior-dominated")]
identifiability[, note := "one-at-a-time profile flatness of the PENALISED objective (others at optimum); not joint identification"]

par_tab <- data.table(parameter = param_labels, role = "calibrated", theta = theta_hat, lower = lo, upper = hi)
par_tab[, `:=`(at_lower = abs(theta - lower) < 1e-3, at_upper = abs(theta - upper) < 1e-3)]
par_tab[, value := c(unlist(bp$tr[trans_names]), as.vector(bp$alpha))]
par_tab[, prior_center := c(unlist(tr_prior[trans_names]), rep(1, n_alpha))]
par_tab[, `:=`(bound_lo = c(vapply(trans_names, function(nm) tb(nm)[1], numeric(1)), rep(ALPHA_LO, n_alpha)),
               bound_hi = c(vapply(trans_names, function(nm) tb(nm)[2], numeric(1)), rep(ALPHA_HI, n_alpha)))]
par_tab[, source := c("age <= 20 annual probability [PROVISIONAL centre 00 block G]; bounds TRANS_BOUNDS",
                      "age <= 20 A->B and adult prior-A->B (shared; adult use = provisional extrapolation) [PROVISIONAL centre]",
                      "age <= 20 annual probability (adult structural zero) [PROVISIONAL centre]",
                      "B->C at age <= 20; adult (> 20) = 1 - (1 - p)^adult_prog_hazard_mult (fixed multiplier) [PROVISIONAL centre]",
                      "age <= 20 annual probability (adult structural zero) [PROVISIONAL centre]",
                      "C->D at age <= 20; adult (> 20) = 1 - (1 - p)^adult_prog_hazard_mult (fixed multiplier) [PROVISIONAL centre]",
                      rep("weak prior on GBD incidence-per-total-population pattern (log-sd SIGMA_ALPHA_STRUCT); Layer-1 alpha = start only", n_alpha))]
eff_hr <- dmi$effects
fixed_tab <- data.table(
  parameter = c("child_last_age (modelling convention)", "p_D_to_C", "p_rhd_death_A/B/C", "p_rhd_death_D (all histories)",
                "late-entry D-death multiplier (removed)", "adult_prog_hazard_mult (adult > 20 B->C and C->D hazard multiplier)",
                "residual_rhd_mortality_mult", "p_B_to_no_rhd_adult",
                paste0("entry_mix_child[", names(NH0$entry_mix_child), "]"), paste0("entry_mix_adult[", names(NH0$entry_mix_adult), "]"),
                "sap_progression_rrr", "sap_progression_hazard_mult (M)", "sap_protection_age_max (attained age)",
                "sap_hr_prior_A_to_B_adult", "sap_hr_C_to_D", "sap_hr_D_to_rhd_death"),
  role = c("fixed", "fixed (structural zero)", "fixed (structural zeros)", "fixed", "fixed at 1", "fixed", "fixed", "fixed",
           rep("fixed (new simplifying assumption)", 8), rep("fixed (applied in reference care)", 6)),
  value = c(NH0$child_last_age, 0, 0, nhe_best$pd, 1, NH0$adult_prog_hazard_mult, NH0$residual_rhd_mortality_mult, NH0$p_B_to_no_rhd_adult,
            NH0$entry_mix_child, NH0$entry_mix_adult,
            eff_hr$sap_progression_rrr, eff_hr$sap_progression_hazard_mult, eff_hr$sap_protection_age_max,
            eff_hr$sap_hr_prior_A_to_B_adult, eff_hr$sap_hr_C_to_D, eff_hr$sap_hr_D_to_rhd_death),
  source = c(rep("00 block G/H [PROVISIONAL reference input / country input / modelling assumption]; not estimated", 16),
             rep(paste0("00 block H: assumed 55% reduction in progression hazards (expert-informed scenario), informed by evidence on prevention of recurrent rheumatic fever ",
                        "[ASSUMPTION; not estimated]: M on treated ", paste(eff_hr$sap_effect_transitions, collapse = " and "),
                        " hazards (same M on both = assumption) at attained age <= sap_protection_age_max; Manyemba et al. RR 0.45 is for RF recurrence (not a HR, not valve progression)"), 3),
             rep("00 block H: multiplier 1 (not SAP-eligible / no direct mortality effect)", 3)))
par_tab <- rbind(par_tab, fixed_tab, fill = TRUE)
# effective ADULT (> 20) annual probabilities the engine actually runs (derived, not fitted)
par_tab <- rbind(par_tab, data.table(
  parameter = c("p_B_to_C_adult (effective, age > 20)", "p_C_to_D_adult (effective, age > 20)", "p_prior_A_to_B_adult (effective, age > 20)"),
  role = "derived",
  value = c(nhe_best$derived$p_B_to_C_adult, nhe_best$derived$p_C_to_D_adult, nhe_best$derived$p_prior_A_to_B_adult),
  source = sprintf("1 - (1 - p_child)^%g for B->C and C->D (adult_prog_hazard_mult); prior-A -> B = p_A_to_B (not scaled)", NH0$adult_prog_hazard_mult)),
  fill = TRUE)
boundary_hits <- par_tab[role == "calibrated" & (at_lower | at_upper), parameter]

#===============================================================================
# 8. INITIALISATION CHECK (2x burn-in at the best solution) + BRIDGE (06 seed)
#===============================================================================
pw_long <- project_window(bp, burn_years = 2L * BURN_YEARS, keep_states = TRUE)
s1 <- pw$start[[1]]; s2 <- pw_long$start[[1]]
grp_tot <- function(st) c(vapply(abcd_stage_aggregate(st), sum, numeric(1)), undiagnosed = sum(abcd_undiagnosed_total(st)),
                          diagnosed = sum(abcd_diagnosed_total(st)), prior_A = sum(abcd_prior_A_total(st)),
                          normalised = sum(abcd_normalised_total(st)))
g1 <- grp_tot(s1); g2 <- grp_tot(s2)
burnin_check <- data.table(check = names(g1), value_burn = unname(g1), value_burn_2x = unname(g2))
prof1 <- agg_grp(abcd_rhd_total(s1)); prof2 <- agg_grp(abcd_rhd_total(s2))
big <- prof1 > 1e-6 * sum(prof1)
burnin_check <- rbind(burnin_check, data.table(
  check = c("age_profile_living_rhd (max over GBD group x sex)", sprintf("rhd_deaths_0_%d_first_window_year", CHILD_AGE_MAX),
            sprintf("rhd_deaths_%dplus_first_window_year", CHILD_AGE_MAX + 1L), sprintf("rhd_deaths_0_%d_last_window_year", CHILD_AGE_MAX),
            "prevalence_last_window_year"),
  value_burn = c(NA, sum(pw$death[child_grp, , 1]), sum(pw$death[!child_grp, , 1]), sum(pw$death[child_grp, , ly]), sum(pw$prev[, , ly])),
  value_burn_2x = c(NA, sum(pw_long$death[child_grp, , 1]), sum(pw_long$death[!child_grp, , 1]), sum(pw_long$death[child_grp, , ly]), sum(pw_long$prev[, , ly]))))
burnin_check[, rel_diff_2x := abs(value_burn_2x - value_burn) / pmax(abs(value_burn), 1e-9)]
burnin_check[check == "age_profile_living_rhd (max over GBD group x sex)", rel_diff_2x := max(abs(prof2[big] - prof1[big]) / prof1[big])]
burn_conv <- as.data.table(pw$burn$convergence)
burnin_check <- merge(burnin_check, burn_conv[, .(check = stock, rel_change_last_cycle, rel_change_vs_half)], by = "check", all = TRUE)
L_burn2x <- decompose(theta_hat, burn_years = 2L * BURN_YEARS)$total
burn_converged <- all(burnin_check[!is.na(rel_diff_2x), rel_diff_2x] < 0.01)
cat(sprintf("  burn-in %d vs %d cycles (%s reconciliation): max rel. difference %.2e (stocks, age profile, childhood/adult deaths) | L %.5f vs %.5f -> %s\n",
            BURN_YEARS, 2L * BURN_YEARS, BURNIN_RECON, max(burnin_check$rel_diff_2x, na.rm = TRUE), L_min, L_burn2x,
            if (burn_converged) "converged" else "NOT converged (initialisation-dependent)"))
# initialisation balancing (per_capita only): per-cycle adjustments by stock group + hand-off accounting
init_balancing <- if (!is.null(pw$burn$init_balancing)) as.data.table(pw$burn$init_balancing) else NULL
handoff_gap <- NA_real_
if (identical(BURNIN_RECON, "per_capita")) {
  hs <- pw$start[[1]]                                          # state entering CAL_YEAR_START (after care initialisation)
  handoff_gap <- max(abs(abcd_living_total(hs)[ages >= 1, ] - POP[ages >= 1, , 1]))
  cat(sprintf("  initialisation balancing: equilibrium adjustment per burn-in cycle disease %.0f (%.2f%% of stock), prior-A %.0f, normalised %.0f, never-onset %.0f | hand-off max |carried - WPP %d| = %.2g persons\n",
              tail(init_balancing$disease, 1), 100 * tail(init_balancing$disease_rel_adj, 1), tail(init_balancing$prior_A, 1),
              tail(init_balancing$normalised, 1), tail(init_balancing$never_onset, 1), CAL_YEAR_START, handoff_gap))
}

initial_state <- NULL; bridge_recon <- NULL; bridge_years <- integer(0); bridge_handoff_gap <- NA_real_; bridge_trace <- NULL
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
    ir_anchor <- ir_total_arr(bp$alpha)[, , n_yr]                 # per TOTAL population, calibrated, at CAL_YEAR_END
    IRb <- array(0, dim = dim(POPb), dimnames = dimnames(POPb))
    for (t in seq_len(nb)) IRb[, , t] <- ir_anchor * trend^(bridge_years[t] - CAL_YEAR_END)
    OTHb <- dmi$rates_by_year$oth_mort[, , 1]                       # production background pattern (held)
    care_b <- abcd_care_schedule(REF_CARE, nb, ages, SEXES)         # reference care continues (no reset)
    prb <- abcd_project(pw$final_state, POPb, IRb, OTHb, nhe_best, eff_cal, surg_cal, care_b, bridge_years,
                        light = FALSE, prev_obs = PREV_OBS)
    fsb <- prb$final_state
    for (k in ABCD_STATE_FIELDS) dimnames(fsb[[k]]) <- list(ages, SEXES)
    agg_b <- abcd_stage_aggregate(fsb)
    initial_state <- list(year = ANALYSIS_YEAR_START, A = agg_b$A, B = agg_b$B, C = agg_b$C, D = agg_b$D,
                          engine_state = fsb,                    # full engine state (stage x care stocks, prior-A, normalised, never-onset, screening history)
                          method = sprintf("calibrated end-of-%d state (care initialised once in %d) forward-projected %d-%d with reference care, calibrated incidence x trend^(y-%d), held 03 background mortality, WPP population",
                                           CAL_YEAR_END, CAL_YEAR_START, min(bridge_years), max(bridge_years), CAL_YEAR_END),
                          pop_deficit_total = prb$pop_deficit_total, care_violation_total = prb$care_violation_total,
                          onset_cap_hits = prb$onset_cap_hits)
    # handoff continuity: the window's final state is the bridge's first start state (no reseeding)
    bridge_handoff_gap <- max(vapply(ABCD_STATE_FIELDS, function(k) max(abs(prb$start[[1]][[k]] - pw$final_state[[k]])), 0))
    bridge_trace <- rbindlist(lapply(seq_len(nb), function(t) { cyc <- prb$cycles[[t]]
      data.table(year = bridge_years[t], classifiable = sum(abcd_rhd_total(prb$start[[t]])), D_start = sum(cyc$D_start),
                 rhd_deaths = sum(cyc$rhd_deaths), first_dx = sum(cyc$n_new_diagnoses), diagnosed = sum(cyc$n_diagnosed),
                 on_bpg = sum(cyc$n_on_optimal_treatment), screen_encounters = sum(cyc$n_screened),
                 prior_A = sum(abcd_prior_A_total(prb$start[[t]])), normalised = sum(abcd_normalised_total(prb$start[[t]])),
                 mass_balance_gap = sum(cyc$rhd_start + cyc$prior_A_start + cyc$never_onset_start + cyc$recovered_start) -
                   sum(cyc$all_mx) - sum(abcd_living_total(cyc$state))) }))
    # reconciliation vs GBD in the bridge years where GBD exists (Number and Rate x WPP pop)
    gbd_all <- as.data.table(readRDS(paste0(wd_data, "temp_baseline_rates_gbd.rds")))
    gb <- gbd_all[location_name == LOCATION & cause_name == "Rheumatic heart disease" & year %in% bridge_years &
                    measure_name %in% c("Prevalence", "Deaths"), .(measure = measure_name, metric = metric_name, sex = sex_name, age_group = age_name, year, val)]
    gb[, child := age_group %in% grp_present[child_grp]]
    popg <- pop_all[, .(pop = sum(Nx)), by = .(sex, age_group = age_to_gbd_group(age), year)]
    gbn <- rbind(gb[metric == "Number", .(gbd_number = sum(val)), by = .(measure, year)],
                 gb[metric == "Number" & measure == "Deaths" & child == TRUE, .(gbd_number = sum(val)), by = .(year)][, measure := sprintf("Deaths 0-%d", CHILD_AGE_MAX)])
    gbr0 <- merge(gb[metric == "Rate"], popg, by = c("sex", "age_group", "year"))
    gbr <- rbind(gbr0[, .(gbd_rate_x_wpp = sum(val / 1e5 * pop)), by = .(measure, year)],
                 gbr0[measure == "Deaths" & child == TRUE, .(gbd_rate_x_wpp = sum(val / 1e5 * pop)), by = .(year)][, measure := sprintf("Deaths 0-%d", CHILD_AGE_MAX)])
    modb <- rbindlist(lapply(seq_len(nb), function(t) data.table(
      year = bridge_years[t], measure = c("Prevalence", "Deaths", sprintf("Deaths 0-%d", CHILD_AGE_MAX)),
      model = c(sum(prb$prev_obs[[t]]), sum(prb$cycles[[t]]$rhd_deaths), sum(prb$cycles[[t]]$rhd_deaths[ages <= CHILD_AGE_MAX, ])))))
    bridge_recon <- merge(merge(modb, gbn, by = c("measure", "year"), all.x = TRUE), gbr, by = c("measure", "year"), all.x = TRUE)
    bridge_recon[, `:=`(ratio_model_over_number = model / gbd_number, ratio_model_over_rate_x_wpp = model / gbd_rate_x_wpp)]
    seed_tot <- sum(abcd_rhd_total(fsb))
    legacy_seed_tot <- sum(dmi$rates_by_year$prev_seed * {
      p1 <- pop_proj[year == ANALYSIS_YEAR_START]; m <- matrix(0, n_age, n_sex, dimnames = list(ages, SEXES))
      m[cbind(match(p1$age, ages), match(p1$sex, SEXES))] <- p1$Nx; m })
    bridge_recon <- rbind(bridge_recon, data.table(measure = "Prevalence (seed entering analysis start)", year = ANALYSIS_YEAR_START,
      model = seed_tot, gbd_number = NA_real_, gbd_rate_x_wpp = legacy_seed_tot,
      ratio_model_over_number = NA_real_, ratio_model_over_rate_x_wpp = seed_tot / legacy_seed_tot), fill = TRUE)
    cat(sprintf("  bridge %d-%d: seed entering %d = %.0f classifiable RHD (%.2fx GBD-%d-rate x pop); D share %.2f%% | diagnosed %.1f%% | on BPG %.1f%% | prior-A %.0f | handoff gap %.2g\n",
                min(bridge_years), max(bridge_years), ANALYSIS_YEAR_START, seed_tot, seed_tot / legacy_seed_tot,
                dmi$meta$RATE_BASE_YEAR, 100 * sum(agg_b$D) / seed_tot,
                100 * sum(abcd_diagnosed_total(fsb)) / seed_tot, 100 * sum(fsb$T_A + fsb$T_B) / seed_tot,
                sum(abcd_prior_A_total(fsb)), bridge_handoff_gap))
  }
}

#===============================================================================
# 9. STATUSES + VALIDATION (structural = hard; numerical / fit = labels)
#===============================================================================
structural_violations <- character(0)
viol <- function(cond, msg) if (isTRUE(cond)) structural_violations <<- c(structural_violations, msg)
viol(any(abs(stage_stocks_year$deaths_outside_D) > 1e-9 * (1 + stage_stocks_year$rhd_deaths)),
     "D-only: RHD deaths occurred outside stage D")   # engine has no A/B/C death pathway; tolerance = summation round-off
viol(!isTRUE(all.equal(nhe_best$pd, NH0$p_rhd_death_D)), "fixed D-mortality input not preserved")
viol(!identical(abcd_nh_fixed(nh_of(bp$tr)), dmi$meta$natural_history_fixed), "fixed natural-history inputs changed during calibration")
viol(any(!is.finite(unlist(stage_stocks_year))) ||
       any(unlist(stage_stocks_year[, .(A, B, C, D, undiagnosed_start, diagnosed_start, prior_A_start, normalised_start, rhd_deaths, new_disease)]) < -1e-9),
     "non-finite or negative stocks/flows")
pop_scale <- sum(POP[, , n_yr])
viol(max(abs(stage_stocks_year$mass_balance_gap)) > 1e-8 * pop_scale, sprintf("population mass balance violated (max gap %.3g persons)", max(abs(stage_stocks_year$mass_balance_gap))))
viol(max(abs(stage_stocks_year$classifiable_identity_gap)) > 1e-6 * (1 + max(stage_stocks_year$living_rhd_start)),
     "classifiable stock-flow identity violated (end = start + new + prior-A re-entry - normalisation - deaths)")
viol(any(stage_stocks_year$dx_nesting_violation > 1e-6), "diagnosis / BPG counts not nested within their eligible stocks")
viol(any(stage_stocks_year$first_dx_exceeds_undiagnosed > 1e-6), "first diagnoses exceed the undiagnosed disease available")
viol(any(stage_stocks_year$P_child_start != 0), "prior-A history at ages <= the child boundary")
viol(any(stage_stocks_year$new_disease_C_or_D > 0) && all(c(NH0$entry_mix_child[c("C", "D")], NH0$entry_mix_adult[c("C", "D")]) == 0),
     "new disease entered C/D although the entry mixes exclude it")
viol(any(stage_stocks_year$screen_outside_window > 0), "screening encounters outside the screening age window")
viol(!is.null(bridge_trace) && max(abs(bridge_trace$mass_balance_gap)) > 1e-8 * pop_scale, "bridge: population mass balance violated")
viol(is.finite(bridge_handoff_gap) && bridge_handoff_gap > 1e-9, "calibration -> bridge hand-off is not continuous (state reseeded or altered)")
deficit_where <- function(cells) if (is.null(cells)) "" else {
  cc <- as.data.table(cells)[order(-deficit)][seq_len(min(.N, 3L))]
  paste0(" at ", paste(sprintf("%s age %d %d (living RHD exceeds a backbone population of %.2f by %.2f)", cc$sex, cc$age, cc$year, cc$population, cc$deficit), collapse = "; "))
}
viol(pw$pop_deficit > 1e-6, sprintf("living RHD stock exceeded the WPP population (total deficit %.3g)%s", pw$pop_deficit, deficit_where(pw$pop_deficit_cells)))
viol(pw$care_violation > 1e-6, sprintf("invalid pooled care state (violation %.3g)", pw$care_violation))
viol(pw$onset_cap_hits > 0, sprintf("first-onset probability hit the engine cap in %d cells", pw$onset_cap_hits))
viol(any(stage_stocks_year$A_adult_start != 0) || any(stage_stocks_year$T_A_adult_start != 0),
     "classifiable A (or A on BPG) not exactly zero at adult ages (must be prior-A history)")
viol(d_fat[D_start_zero == TRUE & zero_deaths_when_undefined_ok == FALSE, .N] > 0, "RHD deaths from D with zero D start stock")
viol(identical(BURNIN_RECON, "per_capita") && !(is.finite(handoff_gap) && handoff_gap < 1e-4),
     sprintf("per-capita burn-in hand-off does not reproduce the %d population (max gap %.3g)", CAL_YEAR_START, handoff_gap))
if (!is.null(initial_state)) {
  viol(initial_state$pop_deficit_total > 1e-6, "bridge: living RHD stock exceeded the WPP population")
  viol(initial_state$care_violation_total > 1e-6, "bridge: invalid pooled care state")
  viol(initial_state$onset_cap_hits > 0, "bridge: onset-probability cap hit")
}
structural_status <- if (length(structural_violations)) paste("FAILED:", paste(structural_violations, collapse = "; ")) else "ok"
numerical_status  <- sprintf("best start '%s' convergence code %d; %d/%d starts converged (code 0); spread %.4f-%.4f; burn-in 2x %s",
                             best_start, as.integer(best$convergence), sum(conv$convergence_code == 0), nrow(conv), min(vals), max(vals),
                             if (burn_converged) "converged" else "NOT converged")
fit_status <- paste(tol_summary[, sprintf("%s: %s (ratio %.2f-%.2f; %d: %.2f)", quantity, status, ratio_min, ratio_max, CAL_YEAR_END, ratio_last_year)],
                    collapse = " | ")
base_status <- if (length(structural_violations)) "FAILED_structural_validation" else
  if (n_stage_used > 0L) "calibrated_with_echo_targets" else "calibrated_no_echo_targets_partial_identification"
status <- if (length(structural_violations)) base_status else
  paste0(base_status, if (fit_ok) "" else "__fit_outside_diagnostic_tolerance")

# human-readable diagnostics (written in every case, incl. failure)
fwrite(fit_long,          paste0(wd_data, "calibration_structural_fit_by_group_year.csv"))
fwrite(fit_year,          paste0(wd_data, "calibration_structural_fit_by_year.csv"))
fwrite(fit_age_last,      paste0(wd_data, "calibration_structural_fit_by_age_last_year.csv"))
fwrite(death_fit,         paste0(wd_data, "calibration_structural_death_age_fit.csv"))
fwrite(prev_fit_age_sex,  paste0(wd_data, "calibration_structural_prevalence_age_sex_fit.csv"))
fwrite(fit_sex_component, paste0(wd_data, "calibration_structural_fit_by_sex_component.csv"))
fwrite(dx_comp,           paste0(wd_data, "calibration_structural_first_diagnosis_stage.csv"))
if (!is.null(bridge_trace)) fwrite(bridge_trace, paste0(wd_data, "calibration_structural_bridge_trace.csv"))
fwrite(tol_tab,           paste0(wd_data, "calibration_structural_fit_tolerance_by_year.csv"))
fwrite(stage_stocks_year, paste0(wd_data, "calibration_structural_stage_stocks_by_year.csv"))
fwrite(stage_by_age_component, paste0(wd_data, "calibration_structural_stage_by_age_component.csv"))
fwrite(d_fat,             paste0(wd_data, "calibration_structural_d_fatality.csv"))
fwrite(par_tab,           paste0(wd_data, "calibration_structural_parameters.csv"))
fwrite(conv,              paste0(wd_data, "calibration_structural_convergence.csv"))
fwrite(burnin_check,      paste0(wd_data, "calibration_structural_burnin_check.csv"))
if (!is.null(init_balancing)) fwrite(init_balancing, paste0(wd_data, "calibration_structural_init_balancing.csv"))
fwrite(data.table(term = c(names(dec$terms), paste0("n_cells_", names(n_cells)), "n_child_cells_below_count_floor"),
                  value = c(dec$terms, n_cells, n_child_floor_masked)),
       paste0(wd_data, "calibration_structural_loss_decomposition.csv"))
if (!is.null(bridge_recon)) fwrite(bridge_recon, paste0(wd_data, "calibration_structural_bridge_reconciliation.csv"))

if (length(structural_violations)) {
  failed <- list(status = status, reason = structural_violations, candidate = bp, theta = theta_hat, L = L_min,
                 decomposition = dec$terms, stage_stocks_by_year = stage_stocks_year,
                 spec = list(engine_spec = ABCD_ENGINE_SPEC, mortality_structure = MORT_STRUCT, sap_mechanism = SAP_MECH))
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
# 10. ASSEMBLE best_parameters + write into the bundle $stage_calibration
#===============================================================================
ir_anchor_final <- ir_total_arr(bp$alpha)[, , n_yr]; dimnames(ir_anchor_final) <- list(ages, SEXES)
alpha_tab <- data.table(location = LOCATION, band = rep(band_labels, times = n_sex), band_lo = rep(band_edges, times = n_sex),
                        sex = rep(SEXES, each = n_band), alpha = as.vector(bp$alpha),
                        alpha_layer1_start = as.vector(alpha_L1))
nh_best <- nh_of(bp$tr); nh_best$derived <- nhe_best$derived
ref_care_record <- c(REF_CARE, list(adherence_factor = REF_ADH, init = INIT,
                                    interpretation = dmi$meta$baseline_care_interpretation))
best_parameters <- list(
  transitions = bp$tr,                      # scalar CHILDHOOD annual probabilities (incl. p_D_to_C = 0); adult B->C / C->D in natural_history$derived
  natural_history = nh_best,                # full natural history the engine ran (fixed + calibrated + derived)
  p_rhd_death = c(A = 0, B = 0, C = 0, D = nhe_best$pd),     # legacy c(A,B,C,D): the single D input (all histories)
  p_rhd_death_note = "D-only; one fixed D input for every disease / diagnosis history (no late-entry multiplier)",
  residual_mortality_mult = 1, residual_mortality_mult_role = "fixed at 1 (not estimated)",
  aggregate_stage_shares = stage_shares,    # derived model shares entering CAL_YEAR_END (diagnostic)
  care_stage_shares = care_stage_shares,    # undiagnosed / diagnosed / on-BPG x stage shares (diagnostic)
  stage_shares_by_age = stage_shares_by_age,
  first_diagnosis_stage = dx_comp,          # stage at first confirmed diagnosis vs PROVISIONAL benchmark (diagnostic)
  # FINAL incidence (joint Layer-2 estimate; applied ONCE; 05 consumes these)
  incidence_alpha = alpha_tab,
  ir_calibrated = ir_anchor_final,          # [age x sex] incidence per TOTAL population at ir_anchor_year
  ir_anchor_year = CAL_YEAR_END,
  ir_denominator = "per TOTAL population: IR_GBD x alpha; the engine converts IR x pop to a first-onset probability per never-onset susceptible each cycle",
  initial_state = initial_state,            # state ENTERING ANALYSIS_YEAR_START (aggregate A-D + $engine_state full state)
  engine_spec = ABCD_ENGINE_SPEC,
  mortality_structure = MORT_STRUCT,
  sap_mechanism = SAP_MECH,
  natural_history_fixed = abcd_nh_fixed(nh_best),
  reference_care_in_calibration = ref_care_record,
  # burn-in initialisation actually used (05 refuses a bundle built under another setting)
  initialisation = list(burnin_reconciliation = BURNIN_RECON, burn_years = BURN_YEARS,
                        handoff_max_gap_persons = handoff_gap,
                        note = if (identical(BURNIN_RECON, "per_capita"))
                          "zero-start burn-in on the first-window-year cross-section; all count stocks rescaled to it after each ageing step (per-capita equilibrium); adjustments = initialisation balancing ($stage_calibration$init_balancing), not flows; observed years / bridge never rescale disease stocks"
                          else "zero-start burn-in; only No-RHD pools reconciled (engine-v2 rule); disease stocks not rescaled"))

calib$stage_calibration$status              <- status
calib$stage_calibration$failure_reason      <- NULL
calib$stage_calibration$best_parameters     <- best_parameters
calib$stage_calibration$accepted_parameter_sets <- accepted
calib$stage_calibration$accepted_parameter_sets_full <- accepted_full
calib$stage_calibration$identifiability     <- identifiability
calib$stage_calibration$parameters          <- par_tab
calib$stage_calibration$convergence         <- conv
calib$stage_calibration$objective_components <- list(
  L_min = L_min, prior_center_L = base_val, prior_center_terms = base_dec$terms, n_starts = N_STARTS,
  multistart_spread = c(min = min(vals), median = stats::median(vals), max = max(vals)),
  weights = c(W_PREV = W_PREV_STRUCT, W_DEATH = W_DEATH_STRUCT, W_INC = W_INC_STRUCT, W_STAGE = W_STAGE_STRUCT,
              prev_child_share = PREV_CHILD_SHARE, death_child_share = CHILD_SHARE),
  prev_death_ratio = W_PREV_STRUCT / W_DEATH_STRUCT, data_weight_total = W_PREV_STRUCT + W_DEATH_STRUCT,
  match_from = MATCH_FROM, window = c(CAL_YEAR_START, CAL_YEAR_END), prev_obs_convention = PREV_OBS,
  target_basis = calib$meta$layer2_target_basis, loss_decomposition = dec$terms, prior_over_data_ratio = prior_over_data,
  n_cells = n_cells, n_child_cells_below_count_floor = n_child_floor_masked,
  loss_excluding_85plus = loss_notail, L_with_2x_burnin = L_burn2x,
  formula = "L = wP*[p*MSLE(prev 0-19) + (1-p)*MSLE(prev 20+)] + wD*[s*MSLE(deaths 0-19) + (1-s)*MSLE(deaths 20+)] + wI*MSLE(inc) + LAMBDA*[priors] + LAMBDA_SMOOTH*smooth(alpha) + wS*stage; each MSLE = weighted mean over its own usable cells (both sexes pooled) of (log(m+eps)-log(t+eps))^2",
  masking = sprintf("missing GBD cells masked (NA, never zero-filled); counts < %g masked from the cell loss (%d observed childhood cells, reported); years < %d weight 0 (not scored); eps = %g; death-target age floor = %d",
                    TARGET_MIN_COUNT, n_child_floor_masked, MATCH_FROM, EPS_LOG, DEATH_AGE_MIN),
  death_target_age_min = DEATH_AGE_MIN, lambda_prior_structural = LAMBDA_PRIOR_STRUCT,
  iv_weights = list(used = USE_IV_W, cap = IV_CAP), age_weights = if (is.null(AGE_W)) "equal" else AGE_W,
  elapsed_min = elapsed_min)
calib$stage_calibration$priors <- list(
  transition_centers = tr_prior, sigma_trans = SIGMA_TRANS,
  sigma_alpha = SIGMA_ALPHA_STRUCT, lambda_smooth_alpha = LAMBDA_SMOOTH_STRUCT, alpha_centers = "1 (GBD incidence pattern)",
  mortality_note = "D-only: one FIXED D mortality input for every history (no prior, not estimated, no late-entry multiplier); A/B/C exact 0",
  burn_years = BURN_YEARS)
calib$stage_calibration$bounds <- list(transitions = TRANS_BOUNDS, alpha = c(ALPHA_LO, ALPHA_HI), boundary_hits = boundary_hits)
calib$stage_calibration$excess_death_attribution <- list(
  model_prev = model_prev_tot, gbd_prev = gbd_prev_tot, prev_ratio = prev_ratio,
  model_deaths = model_death_tot, gbd_deaths = gbd_death_tot, death_ratio = death_ratio,
  child_death_ratio = child_ratio, adult_death_ratio = adult_ratio, inc_ratio = inc_ratio,
  cf_implied = cf_implied, gbd_cf = gbd_cf, stage_shares_last = stage_shares, care_stage_shares = care_stage_shares,
  d_share_proxy_gbd_hf = dmi$meta$d_share_proxy_gbd_hf, classification = attribution)
calib$stage_calibration$fit_by_group_year <- fit_long
calib$stage_calibration$fit_by_year <- fit_year
calib$stage_calibration$death_age_fit <- death_fit
calib$stage_calibration$fit_tolerance_by_year <- tol_tab
calib$stage_calibration$fit_tolerance_summary <- tol_summary
calib$stage_calibration$stage_stocks_by_year <- stage_stocks_year
calib$stage_calibration$stage_by_age_component <- stage_by_age_component
calib$stage_calibration$d_fatality <- d_fat
calib$stage_calibration$burnin_check <- burnin_check
calib$stage_calibration$init_balancing <- init_balancing      # per burn-in cycle (per_capita); NULL for no_rhd_only
calib$stage_calibration$bridge_reconciliation <- bridge_recon
calib$stage_calibration$bridge_trace <- bridge_trace
calib$stage_calibration$bridge_handoff_gap <- bridge_handoff_gap
calib$stage_calibration$fit_by_sex_component <- fit_sex_component
calib$stage_calibration$first_diagnosis_stage <- dx_comp
calib$stage_calibration$demographic_residual_window <- pw$demographic_residual
calib$stage_calibration$stage_targets_n_used <- n_stage_used
calib$stage_calibration$statuses <- list(structural = structural_status, numerical = numerical_status,
                                         epidemiological_fit = fit_status, fit_within_tolerance = fit_ok,
                                         burn_in_converged = burn_converged)
calib$stage_calibration$spec <- list(engine_spec = ABCD_ENGINE_SPEC, mortality_structure = MORT_STRUCT, sap_mechanism = SAP_MECH,
                                     child_last_age = NH0$child_last_age, entry_mix_child = NH0$entry_mix_child,
                                     entry_mix_adult = NH0$entry_mix_adult,
                                     adult_prog_hazard_mult = NH0$adult_prog_hazard_mult,
                                     adult_progression = "adult (> child_last_age) B->C and C->D: p_adult = 1 - (1 - p_child)^adult_prog_hazard_mult",
                                     sap_hr = eff_cal[c("sap_progression_rrr", "sap_progression_hazard_mult", "sap_protection_age_max",
                                                        "sap_effect_transitions", "sap_hr_prior_A_to_B_adult",
                                                        "sap_hr_C_to_D", "sap_hr_D_to_rhd_death")],   # typed named list (not unlist())
                                     sap_eligibility = eff_cal$sap_eligibility,
                                     prevalence_mapping = "classifiable RHD = undiagnosed + diagnosed (off / on BPG), stages A-D; adult prior-A and normalised history excluded",
                                     death_mapping = "RHD deaths from diagnosed + undiagnosed D only",
                                     reference_care = ref_care_record, prev_obs = PREV_OBS,
                                     calendar = list(window = c(CAL_YEAR_START, CAL_YEAR_END), match_from = MATCH_FROM,
                                                     bridge = range(bridge_years), analysis_start = ANALYSIS_YEAR_START),
                                     fit_tolerances = c(prevalence = FIT_TOL_PREV, deaths = FIT_TOL_DEATH),
                                     target_basis = calib$meta$layer2_target_basis, incidence_free = W_INC_STRUCT == 0,
                                     burnin_reconciliation = BURNIN_RECON,
                                     layer1_role = "aggregate proxy; alpha used as starting values only; CF diagnostic/legacy")
calib$meta$stage_calibrated <- TRUE
calib$meta$layer2_spec <- calib$stage_calibration$spec

saveRDS(calib, BUNDLE_FILE)
saveRDS(calib, LASTVALID_FILE)                       # last VALIDATED bundle (never overwritten by a failed run)

fwrite(accepted,        paste0(wd_data, "calibration_structural_accepted_sets.csv"))
fwrite(accepted_full,   paste0(wd_data, "calibration_structural_accepted_sets_full.csv"))
fwrite(identifiability, paste0(wd_data, "calibration_structural_identifiability.csv"))
fwrite(alpha_tab,       paste0(wd_data, "calibration_structural_incidence_alpha.csv"))
fwrite(par_tab[, .(parameter, role, value)], paste0(wd_data, "calibration_structural_best_parameters.csv"))

#===============================================================================
# 11. REPORT
#===============================================================================
cat("\n", strrep("=", 74), "\nSTRUCTURAL (LAYER-2) CALIBRATION — population burden + first diagnosis (", ABCD_ENGINE_SPEC, ")\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("Fixed: D death %.4f (all histories) | adult B->normalised %.3f | adult B->C / C->D hazard multiplier %.2f | entry mix child %s, adult %s | child rules at age <= %d\n",
            nhe_best$pd, NH0$p_B_to_no_rhd_adult, NH0$adult_prog_hazard_mult, paste(NH0$entry_mix_child, collapse = "/"),
            paste(NH0$entry_mix_adult, collapse = "/"), NH0$child_last_age))
cat(sprintf("Calibrated transitions: %s\n",
            paste(sprintf("%s=%.4f (centre %.3f)", trans_names, unlist(bp$tr[trans_names]), unlist(tr_prior[trans_names])), collapse = ", ")))
cat(sprintf("  (childhood values; adult > %d: B->C %.4f, C->D %.4f via hazard x %.2f; A->B also drives adult prior-A -> B unscaled)\n",
            NH0$child_last_age, nhe_best$derived$p_B_to_C_adult, nhe_best$derived$p_C_to_D_adult, NH0$adult_prog_hazard_mult))
cat(sprintf("Incidence alpha (band x sex): %s\n", paste(sprintf("%s=%.2f", alpha_names, as.vector(bp$alpha)), collapse = ", ")))
cat(sprintf("Model stage shares entering %d: A %.1f%% B %.1f%% C %.1f%% D %.2f%% (GBD HF proxy D share %s) | diagnosed %.1f%% | prior-A %.0f | normalised %.0f\n",
            CAL_YEAR_END, 100 * stage_shares[["A"]], 100 * stage_shares[["B"]], 100 * stage_shares[["C"]], 100 * stage_shares[["D"]],
            if (is.na(dmi$meta$d_share_proxy_gbd_hf)) "n/a" else sprintf("%.2f%%", 100 * dmi$meta$d_share_proxy_gbd_hf),
            100 * sum(abcd_diagnosed_total(fs)) / sum(tot_stock), sum(abcd_prior_A_total(fs)), sum(abcd_normalised_total(fs))))
cat("Stage at first confirmed diagnosis (scored years) vs PROVISIONAL benchmark:\n"); print(dx_comp[, .(age_at_diagnosis, stage, first_dx, model_share, benchmark_share)], digits = 3)
cat("\nFit across scored years (aggregate counts, target = GBD Rate x WPP population):\n"); print(tol_summary)
cat("\nDeaths by age component, both sexes (model / target / GBD Number):\n")
print(dcast(death_fit[sex == "Both" & level != "childhood GBD group" & scored == TRUE],
            year ~ age_component, value.var = c("model", "target", "ratio")), digits = 3)
cat("\nChildhood GBD groups, both sexes, last year:\n")
print(death_fit[sex == "Both" & level == "childhood GBD group" & year == CAL_YEAR_END, .(age_group, model, target, gbd_number, ratio, n_cells_below_floor)], digits = 3)
cat(sprintf("Childhood share of RHD deaths %d: model %.3f vs target %.3f\n", CAL_YEAR_END,
            child_share[year == CAL_YEAR_END, child_share_model], child_share[year == CAL_YEAR_END, child_share_target]))
cat("Age-specific fit at the last year (model/target ratio):\n")
print(dcast(fit_age_last[, .(measure, age_group, ratio = round(ratio, 2))], age_group ~ measure, value.var = "ratio")[match(grp_present, age_group)])
cat(sprintf("Realised D death risk %d: all %s | ages 0-20 %s | 21+ %s (input %.4f)\n", CAL_YEAR_END,
            format(round(d_fat[year == CAL_YEAR_END & ages == "all", realized_D_death_risk], 4)),
            format(round(d_fat[year == CAL_YEAR_END & ages == "0-20", realized_D_death_risk], 4)),
            format(round(d_fat[year == CAL_YEAR_END & ages == "21+", realized_D_death_risk], 4)), nhe_best$pd))
cat(sprintf("ATTRIBUTION: %s\n", attribution))
cat(sprintf("STATUS: %s\n  structural: %s\n  numerical : %s\n  boundary hits: %s\n", status, structural_status, numerical_status,
            if (length(boundary_hits)) paste(boundary_hits, collapse = ", ") else "none"))
cat("\nParameter identifiability (one-at-a-time flatness; not joint identification):\n"); print(identifiability[, .(parameter, flatness, status)])
cat(sprintf("\nStage targets: detected %d | usable %d | USED %d | near-optimal sets kept: %d\n",
            calib$stage_calibration$stage_targets_n_detected, calib$stage_calibration$stage_targets_n_usable, n_stage_used, nrow(accepted)))
cat("04b complete. 05 bakes best_parameters (FINAL incidence, natural history and the analysis-start engine state) into the engine.\n")

}  # end if RUN_STRUCT
