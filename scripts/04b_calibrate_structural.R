#===============================================================================
# 04b_calibrate_structural.R
#-------------------------------------------------------------------------------
# LAYER-2 STRUCTURAL A/B/C/D CALIBRATION  (the Stage-2 fix)
#
# Runs the ACTUAL production engine (R/abcd_engine.R — the same one-cycle update
# used by 06) inside the calibration objective, and calibrates the stage
# transitions + a parsimonious severity-ratio stage mortality + rhd_d_fraction to
# reproduce GBD RHD prevalence AND GBD RHD-SPECIFIC deaths, using the Layer-1
# calibrated incidence (04) as the No-RHD -> A inflow.
#
# THIS REPLACES the hand-tuned rhd_mortality_calibration_mult (1/2.5 Indonesia,
# 1/8 Uganda) that formerly reconciled deaths. Now:
#   * stage mortalities carry TIGHT clinical priors (they stay near clinical
#     values, severity-ordered 0 <= m_A <= m_B < m_C < m_D by construction);
#   * transitions + rhd_d_fraction carry LOOSER priors, so the penalized objective
#     reconciles GBD RHD deaths by moving the STOCK / progression (esp. the
#     regression p_A_to_no_rhd that reconciles high GBD incidence with observed
#     prevalence — the mass balance), NOT by suppressing mortality (Part 10);
#   * an OPTIONAL residual mortality multiplier is bounded to [RESID_MULT_LO,
#     RESID_MULT_HI]; a solution sitting AT that bound is a VALIDATION FAILURE
#     (guards against 1/8 reappearing under a new name).
#
# PARTIAL IDENTIFICATION: GBD gives only TOTAL prevalence + TOTAL RHD deaths (no
# stage-specific echo targets), so the stage split is partially identified. This
# script does NOT fabricate stage prevalences: it uses clinical priors/bounds,
# saves ALL near-optimal parameter sets (L <= L_min x (1+delta)), and reports
# which parameters are data-identified vs prior-dominated (a flatness sweep).
# If a local echo target file IS present, its stage-prevalence loss is added.
#
# OUTPUT: extends data/<COUNTRY>/calibrated_rhd_parameters.rds $stage_calibration
#   with best_parameters, accepted_parameter_sets, objective_components, priors,
#   bounds, identifiability, excess_death_attribution, status. Writes CSV diags.
#   05 consumes $stage_calibration$best_parameters (baked ONCE into the engine).
#
# Source AFTER 03 (disease inputs) and 04 (Layer-1 bundle).
#===============================================================================

library(data.table)
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

if (!exists("wd_raw"))  wd_raw  <- paste0(here::here("data-raw"), "/")
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

# shared engine (single source of truth)
if (!exists("abcd_one_cycle")) {
  eng <- if (exists("wd")) paste0(wd, "R/abcd_engine.R") else here::here("R", "abcd_engine.R")
  source(eng)
}

#===============================================================================
# 0. SETTINGS (honour 00_run_all.R; else documented defaults)
#    NOTE: this Layer-2 calibration operates entirely on the CALIBRATION window
#    (CAL_YEAR_START..CAL_YEAR_END, matched from STRUCT_MATCH_FROM) — a period
#    conceptually SEPARATE from the model ANALYSIS period. Its calibrated stage
#    parameters are horizon-agnostic (per-cycle probabilities), so 04b needs no
#    ANALYSIS_YEARS input; 05 applies them over whatever analysis horizon 03 set.
#===============================================================================
LOCATION  <- getp("LOCATION", "Indonesia")
AGE_LO <- as.integer(getp("AGE_LO", 0L)); AGE_HI <- as.integer(getp("AGE_HI", 95L))
CAL_YEAR_START <- as.integer(getp("CAL_YEAR_START", 2000L))
CAL_YEAR_END   <- as.integer(getp("CAL_YEAR_END",   2019L))
MATCH_FROM <- as.integer(getp("STRUCT_MATCH_FROM", 2010L))
SEED <- as.integer(getp("SEED", 42L))
RUN_STRUCT <- isTRUE(getp("run_structural_calibration", TRUE))
CALIB_PDC  <- isTRUE(getp("calibrate_p_D_to_C", FALSE))

## priors (sd) + weights
SIGMA_MORT  <- getp("SIGMA_MORT",  0.35)
SIGMA_TRANS <- getp("SIGMA_TRANS", 1.50)
SIGMA_RESID <- getp("SIGMA_RESID", 0.25)
LAMBDA_PRIOR_STRUCT <- getp("LAMBDA_PRIOR_STRUCT", 1.0)
W_PREV_STRUCT  <- getp("W_PREV_STRUCT",  1)
W_DEATH_STRUCT <- getp("W_DEATH_STRUCT", 2)
EPS_LOG <- 1e-3

## bounds
MD_LO <- getp("MD_LO", 0.03); MD_HI <- getp("MD_HI", 0.25)
RESID_LO <- getp("RESID_MULT_LO", 0.5); RESID_HI <- getp("RESID_MULT_HI", 2.0)
N_STARTS <- as.integer(getp("STRUCT_N_STARTS", 10L))
NEAROPT_DELTA <- getp("STRUCT_NEAROPT_DELTA", 0.05)

## prior CENTERS (from 00 / 03 disease inputs)
BUNDLE_FILE  <- paste0(wd_data, "calibrated_rhd_parameters.rds")
DISEASE_FILE <- paste0(wd_data, "disease_model_inputs.rds")
GBD_FILE     <- paste0(wd_data, "temp_baseline_rates_gbd.rds")
STAGE_TARGET_CANDIDATES <- c(paste0(wd_data, "calibration_targets_stage.csv"),
                             paste0(wd_raw,  "calibration_targets_stage.csv"))

for (f in c(BUNDLE_FILE, DISEASE_FILE, GBD_FILE))
  if (!file.exists(f)) stop("04b: missing input ", f, " (run 03 and 04 first).", call. = FALSE)

calib <- readRDS(BUNDLE_FILE)
dmi   <- readRDS(DISEASE_FILE)
gbd   <- as.data.table(readRDS(GBD_FILE))

tr_prior <- dmi$transitions                              # transition prior centers
pd_prior <- dmi$p_rhd_death                              # clinical mortality prior centers c(A,B,C,D)
dfrac_prior <- dmi$meta$rhd_d_fraction
eff03  <- dmi$effects; surg03 <- dmi$surgery             # (inert here: cov = 0)

if (!RUN_STRUCT) {
  message("04b: run_structural_calibration = FALSE -> leaving stage params at 03 priors.")
  calib$stage_calibration$status <- "structural_disabled_using_priors"
  saveRDS(calib, BUNDLE_FILE)
  message("04b: wrote passthrough status; 05 will use 03 stage params.")
} else {

message(sprintf("── 04b structural A/B/C/D calibration : %s | match %d-%d | engine=R/abcd_engine.R ──",
                LOCATION, MATCH_FROM, CAL_YEAR_END))

#===============================================================================
# 1. AGE-GROUP MAP + GBD TARGETS (Number: RHD prevalence + RHD deaths)
#===============================================================================
age_to_gbd_group <- function(a) fcase(
  a < 1,"<1 year", a < 2,"12-23 months", a < 5,"2-4 years", a < 10,"5-9 years",
  a < 15,"10-14 years", a < 20,"15-19 years", a < 25,"20-24 years", a < 30,"25-29 years",
  a < 35,"30-34 years", a < 40,"35-39 years", a < 45,"40-44 years", a < 50,"45-49 years",
  a < 55,"50-54 years", a < 60,"55-59 years", a < 65,"60-64 years", a < 70,"65-69 years",
  a < 75,"70-74 years", a < 80,"75-79 years", a < 85,"80-84 years", a < 90,"85-89 years",
  a < 95,"90-94 years", default = "95+ years")

ages <- AGE_LO:AGE_HI; n_age <- length(ages)
yrs  <- CAL_YEAR_START:CAL_YEAR_END; n_yr <- length(yrs)
grp_lab_age <- age_to_gbd_group(ages)
grp_present <- unique(grp_lab_age)
grp_int     <- match(grp_lab_age, grp_present)
G <- length(grp_present)
match_cols <- which(yrs >= MATCH_FROM)                   # year columns entering the loss

gbd_t <- gbd[metric_name == "Number" & location_name == LOCATION &
               cause_name == "Rheumatic heart disease" &
               measure_name %in% c("Prevalence", "Deaths") &
               year %in% yrs & age_name %in% grp_present,
             .(sex = sex_name, age = age_name, year, measure = measure_name, val)]
SEXES <- intersect(c("Female", "Male"), unique(gbd_t$sex))

# GBD target [G x n_yr] per (sex, measure)
gbd_mat <- function(sx, meas) {
  d <- gbd_t[sex == sx & measure == meas]
  m <- matrix(0, G, n_yr, dimnames = list(grp_present, as.character(yrs)))
  ri <- match(d$age, grp_present); ci <- match(d$year, yrs); ok <- !is.na(ri) & !is.na(ci)
  m[cbind(ri[ok], ci[ok])] <- d$val[ok]; m
}

#===============================================================================
# 2. PER-SEX MODEL INPUTS from the Layer-1 bundle ($tp: IR, PREVt0, BG.mx, Nx)
#===============================================================================
tp <- as.data.table(calib$tp)[location == LOCATION & age %in% ages & year %in% yrs]
sex_data <- lapply(SEXES, function(sx) {
  d <- tp[sex == sx]
  to_mat <- function(col) {
    m <- matrix(0, n_age, n_yr, dimnames = list(ages, as.character(yrs)))
    m[cbind(match(d$age, ages), match(d$year, yrs))] <- d[[col]]; m
  }
  list(IR = to_mat("IR"), PREV0 = to_mat("PREVt0"), OTH = to_mat("BG.mx"),
       POP = to_mat("Nx"),
       gbd_prev = gbd_mat(sx, "Prevalence"), gbd_death = gbd_mat(sx, "Deaths"))
})
names(sex_data) <- SEXES

#===============================================================================
# 3. PARAMETER PACK / UNPACK  (unconstrained theta -> named param set)
#    Stage mortality: m_D = level; ratios via nested logistic so 0<=m_A<=m_B<m_C<m_D.
#===============================================================================
abc_raw <- c(A = 0.565, B = 0.272, C = 0.162); abc <- abc_raw / sum(abc_raw)  # Cannon split
stage_split_of <- function(dfrac) c(abc * (1 - dfrac), D = dfrac)

trans_names <- c("p_A_to_no_rhd","p_A_to_B","p_B_to_A","p_B_to_C","p_C_to_B","p_C_to_D")
if (CALIB_PDC) trans_names <- c(trans_names, "p_D_to_C")
n_tr <- length(trans_names)

# layout: [logit(trans) x n_tr, log(mD), thetaC, thetaB, thetaA, log(resid)]
# NOTE: rhd_d_fraction is NOT a free parameter — the zero-start burn-in makes the
# initial split irrelevant; the seed D-share is the transition-implied EQUILIBRIUM
# D-share, reported as the calibrated rhd_d_fraction after optimisation.
idx <- list(tr = 1:n_tr, mD = n_tr + 1, rC = n_tr + 2, rB = n_tr + 3, rA = n_tr + 4,
            res = n_tr + 5)
npar <- n_tr + 5

unpack <- function(theta) {
  tr <- as.list(plogis(theta[idx$tr])); names(tr) <- trans_names
  if (!CALIB_PDC) tr$p_D_to_C <- getp("p_D_to_C", 0)
  mD <- exp(theta[idx$mD])
  rC <- plogis(theta[idx$rC]); rB <- rC * plogis(theta[idx$rB]); rA <- rB * plogis(theta[idx$rA])
  pd <- c(A = rA * mD, B = rB * mD, C = rC * mD, D = mD)
  resid <- exp(theta[idx$res])
  list(tr = tr, pd = pd, mD = mD, ratios = c(A = rA, B = rB, C = rC), resid = resid)
}

# bounds (unconstrained scale)
lo <- numeric(npar); hi <- numeric(npar)
lo[idx$tr] <- qlogis(1e-4);  hi[idx$tr] <- qlogis(0.6)
lo[idx$mD] <- log(MD_LO);    hi[idx$mD] <- log(MD_HI)
lo[c(idx$rC, idx$rB, idx$rA)] <- -8; hi[c(idx$rC, idx$rB, idx$rA)] <- 8
lo[idx$res] <- log(RESID_LO); hi[idx$res] <- log(RESID_HI)

# prior-center start (theta0): transitions at 03 centers; mortality at clinical
# priors (recover ratios from clinical pd); resid = 1.
theta0 <- numeric(npar)
theta0[idx$tr] <- qlogis(pmin(pmax(unlist(tr_prior[trans_names]), 1e-4), 0.6))
mD0 <- pmin(pmax(pd_prior[["D"]], MD_LO), MD_HI)
theta0[idx$mD] <- log(mD0)
rC0 <- pmin(pmax(pd_prior[["C"]]/mD0, 1e-3), 1-1e-3)
rB0 <- pmin(pmax(pd_prior[["B"]]/pd_prior[["C"]], 1e-3), 1-1e-3)
rA0 <- pmin(pmax(pd_prior[["A"]]/pd_prior[["B"]], 1e-3), 1-1e-3)
theta0[idx$rC] <- qlogis(rC0); theta0[idx$rB] <- qlogis(rB0); theta0[idx$rA] <- qlogis(rA0)
theta0[idx$res] <- 0
theta0 <- pmin(pmax(theta0, lo), hi)

BURN_YEARS <- as.integer(getp("STRUCT_BURN_YEARS", 80L))

#===============================================================================
# 4. ENGINE PROJECTION (per sex) + OBJECTIVE
#===============================================================================
eff0  <- list(sap_rrr_rhd_death = 0, eff_surgery_C_to_D = 0, eff_surgery_D_to_rhd_death = 0)
surg0 <- list(frac_C_requiring_surgery = 0, frac_D_requiring_surgery = 0)
cov0  <- list(treatment = 0, surgery = 0)
age_shift1 <- function(M) { n <- nrow(M); N <- matrix(0, n, 1); N[2:n, 1] <- M[1:(n-1), 1]
                            N[n, 1] <- N[n, 1] + M[n, 1]; N }

# zero-start burn-in with year-1 rates to reach the transition-implied stage
# EQUILIBRIUM (removes the Cannon-seed artifact; prevalence is GENERATED, not
# seeded). Returns the equilibrated stocks entering the window.
burn_in <- function(sd, tr, pd) {
  A <- B <- C <- D <- matrix(0, n_age, 1)
  pop1 <- sd$POP[, 1, drop = FALSE]; ir1 <- sd$IR[, 1, drop = FALSE]; oth1 <- sd$OTH[, 1, drop = FALSE]
  for (b in seq_len(BURN_YEARS)) {
    cyc <- abcd_one_cycle(A, B, C, D, pop1, ir1, oth1, tr, pd, eff0, surg0, cov0, light = TRUE)
    A <- age_shift1(cyc$A); B <- age_shift1(cyc$B); C <- age_shift1(cyc$C); D <- age_shift1(cyc$D)
  }
  list(A = A, B = B, C = C, D = D)
}

project_structural <- function(sd, tr, pd) {
  s <- burn_in(sd, tr, pd); A <- s$A; B <- s$B; C <- s$C; D <- s$D
  prev_rec <- matrix(0, n_age, n_yr); death_rec <- matrix(0, n_age, n_yr)
  for (t in seq_len(n_yr)) {
    popm <- sd$POP[, t, drop = FALSE]; irm <- sd$IR[, t, drop = FALSE]
    oth  <- sd$OTH[, t, drop = FALSE]
    cyc <- abcd_one_cycle(A, B, C, D, popm, irm, oth, tr, pd, eff0, surg0, cov0, light = TRUE)
    prev_rec[, t]  <- (A + B + C + D)[, 1]
    death_rec[, t] <- cyc$rhd_deaths[, 1]
    A <- age_shift1(cyc$A); B <- age_shift1(cyc$B); C <- age_shift1(cyc$C); D <- age_shift1(cyc$D)
  }
  list(prev = rowsum(prev_rec, grp_int), death = rowsum(death_rec, grp_int))
}

# equilibrated per-age stocks entering the FINAL window year — the calibrated,
# dynamically self-consistent stage distribution used to SEED 06 (its SHARES,
# applied to GBD prevalence in 05).
engine_final_state <- function(sd, tr, pd) {
  s <- burn_in(sd, tr, pd); A <- s$A; B <- s$B; C <- s$C; D <- s$D
  As <- Bs <- Cs <- Ds <- NULL
  for (t in seq_len(n_yr)) {
    if (t == n_yr) { As <- A[, 1]; Bs <- B[, 1]; Cs <- C[, 1]; Ds <- D[, 1] }
    cyc <- abcd_one_cycle(A, B, C, D, sd$POP[, t, drop = FALSE], sd$IR[, t, drop = FALSE],
                          sd$OTH[, t, drop = FALSE], tr, pd, eff0, surg0, cov0, light = TRUE)
    A <- age_shift1(cyc$A); B <- age_shift1(cyc$B); C <- age_shift1(cyc$C); D <- age_shift1(cyc$D)
  }
  list(A = As, B = Bs, C = Cs, D = Ds)
}

# optional local echo stage targets (partial-identification loss add-on)
stage_target_file <- STAGE_TARGET_CANDIDATES[file.exists(STAGE_TARGET_CANDIDATES)][1]
have_stage_targets <- !is.na(stage_target_file)

loss_measure <- function(model, gbd) {
  msk <- gbd >= 1
  msk[, -match_cols] <- FALSE                # only score years >= MATCH_FROM
  if (!any(msk)) return(0)
  mean((log(model[msk] + EPS_LOG) - log(gbd[msk] + EPS_LOG))^2)
}

objective <- function(theta) {
  p <- unpack(theta)
  LP <- 0; LD <- 0
  for (sx in SEXES) {
    pr <- project_structural(sex_data[[sx]], p$tr, p$pd * p$resid)
    LP <- LP + loss_measure(pr$prev,  sex_data[[sx]]$gbd_prev)
    LD <- LD + loss_measure(pr$death, sex_data[[sx]]$gbd_death)
  }
  # priors: mortality TIGHT (log), transitions LOOSE (logit), resid modest (log)
  pri_mort <- sum((log(p$pd) - log(pd_prior[c("A","B","C","D")]))^2) / SIGMA_MORT^2
  pri_tr   <- sum((qlogis(unlist(p$tr[trans_names])) -
                   qlogis(pmin(pmax(unlist(tr_prior[trans_names]),1e-4),0.6)))^2) / SIGMA_TRANS^2
  pri_res  <- (log(p$resid)^2) / SIGMA_RESID^2
  data_term  <- W_PREV_STRUCT * LP + W_DEATH_STRUCT * LD
  prior_term <- LAMBDA_PRIOR_STRUCT * (pri_mort + pri_tr + pri_res)
  attr_val <- data_term + prior_term
  if (!is.finite(attr_val)) return(1e12)
  attr_val
}

#===============================================================================
# 5. OPTIMISE — multi-start L-BFGS-B; prior-center start always evaluated
#===============================================================================
run_lbfgs <- function(p0) tryCatch(
  optim(p0, objective, method = "L-BFGS-B", lower = lo, upper = hi,
        control = list(factr = 1e9, maxit = 120)),
  error = function(e) NULL)

set.seed(SEED)
starts <- c(list(theta0), lapply(seq_len(N_STARTS), function(i) runif(npar, lo, hi)))
fits <- lapply(starts, run_lbfgs)
fits <- Filter(Negate(is.null), fits)
vals <- vapply(fits, function(f) f$value, numeric(1))
best <- fits[[which.min(vals)]]
# guarantee not worse than the prior-center baseline
base_val <- objective(theta0)
if (base_val < best$value) best <- list(par = theta0, value = base_val)
theta_hat <- best$par; L_min <- best$value
bp <- unpack(theta_hat)

cat(sprintf("  optimum L = %.4f (prior-center baseline L = %.4f) | multi-start spread %.4f/%.4f/%.4f\n",
            L_min, base_val, min(vals), stats::median(vals), max(vals)))

#===============================================================================
# 6. FIT DECOMPOSITION + EXCESS-DEATH ATTRIBUTION (at the last matched year)
#===============================================================================
model_prev_tot <- 0; model_death_tot <- 0; gbd_prev_tot <- 0; gbd_death_tot <- 0
for (sx in SEXES) {
  pr <- project_structural(sex_data[[sx]], bp$tr, bp$pd * bp$resid)
  ly <- n_yr
  model_prev_tot  <- model_prev_tot  + sum(pr$prev[, ly])
  model_death_tot <- model_death_tot + sum(pr$death[, ly])
  gbd_prev_tot    <- gbd_prev_tot    + sum(sex_data[[sx]]$gbd_prev[, ly])
  gbd_death_tot   <- gbd_death_tot   + sum(sex_data[[sx]]$gbd_death[, ly])
}
# calibrated, dynamically self-consistent (equilibrium) per-age stage fractions,
# per sex -> used to SEED 06 (so base-year deaths reflect the CALIBRATED
# distribution). Aggregate prevalence-weighted shares feed the attribution.
final_states <- lapply(SEXES, function(sx)
  engine_final_state(sex_data[[sx]], bp$tr, bp$pd * bp$resid))
names(final_states) <- SEXES
stage_shares_by_age <- rbindlist(lapply(SEXES, function(sx) {
  fs <- final_states[[sx]]; tot <- fs$A + fs$B + fs$C + fs$D
  data.table(sex = sx, age = ages,
             fracA = ifelse(tot > 0, fs$A / tot, 1), fracB = ifelse(tot > 0, fs$B / tot, 0),
             fracC = ifelse(tot > 0, fs$C / tot, 0), fracD = ifelse(tot > 0, fs$D / tot, 0))
}))
tot_stock <- Reduce(`+`, lapply(final_states, function(fs) c(A=sum(fs$A),B=sum(fs$B),C=sum(fs$C),D=sum(fs$D))))
stage_shares <- tot_stock / sum(tot_stock)
bp$dfrac <- unname(stage_shares[["D"]])   # calibrated rhd_d_fraction = equilibrium D-share
cf_implied <- model_death_tot / model_prev_tot
gbd_cf     <- gbd_death_tot / gbd_prev_tot            # GBD implied aggregate case fatality
cf_gap     <- cf_implied / gbd_cf
prev_ratio  <- model_prev_tot  / gbd_prev_tot
death_ratio <- model_death_tot / gbd_death_tot
resid_at_bound <- (abs(bp$resid - RESID_LO) < 1e-3) || (abs(bp$resid - RESID_HI) < 1e-3)
min_plausible_cf <- unname(bp$pd[["A"]] * bp$resid)  # CF floor: everyone in stage A at clinical m_A

# attribution: name the dominant driver of any RHD-deaths mismatch (Part 10/13).
attribution <- if (resid_at_bound)
    "residual mortality multiplier AT BOUND (VALIDATION FAILURE — revisit priors/bounds/stock)" else
  if (abs(log(death_ratio)) < 0.25 && abs(log(prev_ratio)) < 0.25)
    "GBD prevalence AND deaths both reproduced with clinically-plausible mortality" else
  if (gbd_cf < min_plausible_cf)
    sprintf(paste0("GBD aggregate case fatality (%.2g/yr) is BELOW the lowest clinically-plausible ",
      "stage mortality (m_A=%.2g/yr): even the calibrated mostly-A equilibrium (A=%.0f%%) implies ",
      "CF=%.2g/yr (%.1fx GBD). The RHD-deaths gap traces to GBD's very high echo-prevalence:death ",
      "ratio (mild-disease-dominated prevalence), NOT to model mortality (kept clinical); closing it ",
      "would require implausibly low mortality."),
      gbd_cf, min_plausible_cf, 100*stage_shares[["A"]], cf_implied, cf_gap) else
  if (abs(log(prev_ratio)) > 0.25)
    sprintf("TOTAL PREVALENCE mismatch (model %.2fx GBD): incidence x dwell vs GBD prevalence; the regression rate p_A_to_no_rhd=%.3f sets the dwell.", prev_ratio, bp$tr$p_A_to_no_rhd) else
  if ((stage_shares[["C"]] + stage_shares[["D"]]) > 0.30)
    sprintf("excess C/D stage SHARE (%.0f%%): progression p_C_to_D=%.3f drives severe-stage stock.",
            100*(stage_shares[["C"]]+stage_shares[["D"]]), bp$tr$p_C_to_D) else
  "RHD deaths reconciled by stock/progression at clinical mortality"

#===============================================================================
# 7. PARTIAL IDENTIFICATION — near-optimal set + per-parameter flatness sweep
#===============================================================================
# near-optimal set from the multi-starts (L <= L_min x (1+delta))
acc <- which(vals <= L_min * (1 + NEAROPT_DELTA))
accepted <- rbindlist(lapply(acc, function(k) {
  pk <- unpack(fits[[k]]$par)
  data.table(L = vals[k], m_A = pk$pd[["A"]], m_B = pk$pd[["B"]], m_C = pk$pd[["C"]],
             m_D = pk$pd[["D"]], resid = pk$resid,
             p_A_to_no_rhd = pk$tr$p_A_to_no_rhd, p_B_to_C = pk$tr$p_B_to_C,
             p_C_to_D = pk$tr$p_C_to_D)
}), fill = TRUE)

# flatness: sweep each param across its bounds (others at optimum); fraction of the
# range within delta of L_min -> wide = prior-dominated, narrow = data-identified.
sweep_frac <- function(j, npts = 21) {
  grid <- seq(lo[j], hi[j], length.out = npts)
  Ls <- vapply(grid, function(v) { th <- theta_hat; th[j] <- v; objective(th) }, numeric(1))
  mean(Ls <= L_min * (1 + NEAROPT_DELTA))
}
param_labels <- c(trans_names, "log_mD", "ratio_C", "ratio_B", "ratio_A", "log_resid")
identifiability <- data.table(
  parameter = param_labels,
  flatness = vapply(seq_len(npar), sweep_frac, numeric(1)))
identifiability[, status := fifelse(flatness < 0.5, "data-identified", "prior-dominated")]

#===============================================================================
# 8. ASSEMBLE best_parameters + write into the bundle $stage_calibration
#===============================================================================
best_parameters <- list(
  transitions = bp$tr,
  # p_rhd_death consumed by 06 is the EFFECTIVE per-stage mortality = untreated
  # clinical estimate x the (off-bound) residual multiplier, baked in ONCE here.
  p_rhd_death = bp$pd * bp$resid,
  p_rhd_death_untreated = bp$pd,           # clinical untreated estimate (pre-residual), for reporting
  rhd_d_fraction = bp$dfrac,               # = equilibrium D-share (derived, not a free param)
  residual_mortality_mult = bp$resid,
  stage_split = stage_shares,                       # aggregate calibrated equilibrium A/B/C/D shares
  stage_shares_by_age = stage_shares_by_age,        # CALIBRATED per-age seed for 06 (05 uses this)
  aggregate_stage_shares = stage_shares)

status <- if (have_stage_targets) "calibrated_with_echo_targets" else
          "calibrated_no_echo_targets_partial_identification"

calib$stage_calibration$status              <- status
calib$stage_calibration$best_parameters     <- best_parameters
calib$stage_calibration$accepted_parameter_sets <- accepted
calib$stage_calibration$identifiability     <- identifiability
calib$stage_calibration$objective_components <- list(
  L_min = L_min, prior_center_L = base_val, n_starts = N_STARTS,
  multistart_spread = c(min = min(vals), median = stats::median(vals), max = max(vals)),
  weights = c(W_PREV = W_PREV_STRUCT, W_DEATH = W_DEATH_STRUCT),
  match_from = MATCH_FROM)
calib$stage_calibration$priors <- list(
  mortality_centers = pd_prior, sigma_mort = SIGMA_MORT,
  transition_centers = tr_prior, sigma_trans = SIGMA_TRANS,
  sigma_resid = SIGMA_RESID,
  rhd_d_fraction_note = "derived as the equilibrium D-share (not a free parameter)",
  burn_years = BURN_YEARS)
calib$stage_calibration$bounds <- list(
  m_D = c(MD_LO, MD_HI), residual_mult = c(RESID_LO, RESID_HI))
calib$stage_calibration$excess_death_attribution <- list(
  model_prev = model_prev_tot, gbd_prev = gbd_prev_tot, prev_ratio = prev_ratio,
  model_deaths = model_death_tot, gbd_deaths = gbd_death_tot, death_ratio = death_ratio,
  cf_implied = cf_implied, gbd_cf = gbd_cf, cf_gap = cf_gap,
  min_plausible_cf = min_plausible_cf, stage_shares_last = stage_shares,
  residual_mult = bp$resid, residual_at_bound = resid_at_bound,
  classification = attribution)
calib$meta$stage_calibrated <- TRUE

saveRDS(calib, BUNDLE_FILE)

# human-readable diagnostics
fwrite(accepted,        paste0(wd_data, "calibration_structural_accepted_sets.csv"))
fwrite(identifiability, paste0(wd_data, "calibration_structural_identifiability.csv"))
fwrite(data.table(parameter = c(names(unlist(bp$tr)), names(bp$pd), "rhd_d_fraction", "residual_mult"),
                  value = c(unlist(bp$tr), bp$pd, bp$dfrac, bp$resid)),
       paste0(wd_data, "calibration_structural_best_parameters.csv"))

#===============================================================================
# 9. VALIDATION + REPORT
#===============================================================================
cat("\n", strrep("=", 70), "\nSTRUCTURAL (LAYER-2) VALIDATION\n", strrep("=", 70), "\n", sep = "")
mm <- bp$pd
ordering_ok <- (mm[["A"]] <= mm[["B"]] + 1e-12) && (mm[["B"]] < mm[["C"]]) && (mm[["C"]] < mm[["D"]])
if (!ordering_ok) stop("Stage mortality ordering 0<=m_A<=m_B<m_C<m_D violated.", call. = FALSE)
if (resid_at_bound)
  stop(sprintf(paste0("Residual mortality multiplier sat at its bound (%.3f in [%.2f,%.2f]) -- ",
       "this is a VALIDATION FAILURE (1/8-style scaling reappearing). Revisit priors/bounds/stock."),
       bp$resid, RESID_LO, RESID_HI), call. = FALSE)

cat(sprintf("Stage mortality (untreated, /yr): A=%.4f B=%.4f C=%.4f D=%.4f | ordering OK\n",
            mm[["A"]], mm[["B"]], mm[["C"]], mm[["D"]]))
cat(sprintf("m_D in [%.2f,%.2f]? %s | residual mult = %.3f (bound [%.2f,%.2f], at bound: %s)\n",
            MD_LO, MD_HI, mm[["D"]] >= MD_LO && mm[["D"]] <= MD_HI,
            bp$resid, RESID_LO, RESID_HI, resid_at_bound))
cat(sprintf("Transitions: %s\n", paste(sprintf("%s=%.4f", names(unlist(bp$tr)), unlist(bp$tr)), collapse=", ")))
cat(sprintf("rhd_d_fraction (equilibrium D-share) = %.3f | calibrated equilibrium split A/B/C/D = %s\n",
            bp$dfrac, paste(sprintf("%.0f%%", 100*stage_shares[c("A","B","C","D")]), collapse="/")))
cat(sprintf("\nFit at %d: model prev=%.0f vs GBD %.0f (ratio %.2f) | model RHD deaths=%.0f vs GBD %.0f (ratio %.2f)\n",
            CAL_YEAR_END, model_prev_tot, gbd_prev_tot, prev_ratio,
            model_death_tot, gbd_death_tot, death_ratio))
cat(sprintf("Implied total CF = %.5f | Layer-1 anchor CF (mean) = %.5f\n",
            cf_implied, mean(as.data.table(calib$tp)[year==CAL_YEAR_END, CF], na.rm=TRUE)))
cat(sprintf("Model stage shares @%d: A=%.0f%% B=%.0f%% C=%.0f%% D=%.0f%%\n", CAL_YEAR_END,
            100*stage_shares[["A"]],100*stage_shares[["B"]],100*stage_shares[["C"]],100*stage_shares[["D"]]))
cat(sprintf("EXCESS-DEATH ATTRIBUTION: %s\n", attribution))
cat("\nParameter identifiability (flatness = fraction of range within near-optimal delta):\n")
print(identifiability)
cat(sprintf("\nStage-calibration status: %s | near-optimal sets kept: %d\n", status, nrow(accepted)))
if (!have_stage_targets)
  cat("PARTIAL IDENTIFICATION: no local echo stage targets -> stage split is prior-informed;\n",
      "  see accepted_parameter_sets + identifiability. Populate calibration_targets_stage.csv to pin it.\n", sep="")
cat("04b complete. 05 will bake best_parameters into the engine (no 1/8-style scalar).\n")

}  # end if RUN_STRUCT
