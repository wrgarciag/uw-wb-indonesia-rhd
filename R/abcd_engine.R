# ==============================================================================
# SHARED A/B/C/D ONE-CYCLE ENGINE  (single source of truth)
# R/abcd_engine.R
#
# One annual Markov cycle used by the production runner (06) and the structural
# calibration (04b: burn-in, scored window, bridge), so there is exactly one
# implementation of the stage dynamics, care update, ageing and observation.
#
# DEFAULT SPECIFICATION (ABCD_ENGINE_SPEC = "abcd_first_diagnosis_v1")
#   National population model that separates DISEASE HISTORY from FIRST CONFIRMED
#   DIAGNOSIS. Per [age x sex] cell the state holds (all counts of persons):
#     N_ns, N_sc       never-onset (no prior modelled disease), never / ever screened
#     U_s_ns, U_s_sc   UNDIAGNOSED classifiable RHD in stage s = A/B/C/D, never screened /
#                      screened without a confirmed diagnosis (missed or unconfirmed)
#     X_s              CONFIRMED DIAGNOSIS (diagnosed / care cohort), stage s, not on BPG
#     T_A, T_B         confirmed diagnosis ON BPG (SAP-eligible stages: A at age <= 20, B)
#     P_ns, P_sc, P_dx ADULT PRIOR-A lesion history (surviving A reclassified at the
#                      20 -> 21 ageing step): NOT classifiable RHD prevalence, not SAP-
#                      eligible; split by screening / prior-diagnosis history
#     R_ns, R_sc, R_dx NORMALISED history (no recurrence; background mortality only)
#   Classifiable RHD prevalence (GBD mapping) = sum of U, X and T stocks (A-D); prior-A
#   and normalised history are excluded. First confirmed diagnosis moves U -> X once;
#   nobody re-enters U after diagnosis, so a person is first-diagnosed at most once.
#
#   NATURAL HISTORY (attained age; boundary = child_last_age, default 20):
#     age <= 20 : A->normalised, A->B, B->A, B->C, C->B, C->D, D->RHD death
#     age  > 20 : prior-A->B (re-entry, shared A->B parameter), B->C, B->normalised
#                 (fixed), C->D, D->RHD death; no adult A, B->A = C->B = 0
#     D->C = 0 exactly; RHD death ONLY from D (A/B/C exact zeros), one base D input for
#     every disease / diagnosis history (no late-entry multiplier, no adult progression
#     multiplier). Hazards depend on attained age only, never on when diagnosis occurred.
#   NEW MODELLED DISEASE (first classifiable entry) from the never-onset pool enters the
#     stage given by the age-specific ENTRY MIX (00: child A = 1, adult B = 1) at year
#     end (undiagnosed, carrying its screening history); it first faces stage hazards in
#     the next cycle. Screening never creates disease.
#   AGE BOUNDARY (ageing 20 -> 21): surviving A (undiagnosed, diagnosed, on BPG) becomes
#     adult prior-A (people conserved; BPG ends; recorded as reclassification, not
#     remission or death). Prior-A -> B is recorded as RE-ENTRY into classifiable
#     disease, separate from first entry; previously diagnosed prior-A re-enters the
#     diagnosed cohort (no repeat first diagnosis).
#
# CARE UPDATE (start of each cycle, dynamic mode):
#   screening among NEVER-SCREENED people without a prior confirmed diagnosis at window
#   ages: annual probability q; encounters = q x (N_ns + U_ns + P_ns + R_ns);
#   true positives = q x sens x U_ns (perfect specificity); confirmation completed =
#   true positives x confirm_completion; first confirmed (screen) = completed x
#   confirm_sensitivity. Screened people without a confirmed diagnosis stay undiagnosed
#   (U_sc): disease is never removed. Clinical detection k (confirmed clinical
#   diagnosis) applies to all remaining undiagnosed disease (no overlap with screening).
#   BPG: on-BPG retained with probability r; dropouts rejoin the diagnosed-off-BPG pool;
#   initiation i among ELIGIBLE diagnosed people not on BPG (A at age <= 20, B).
#   Protection for on-BPG people: h x [1 - adh x (1 - HR)] with adh = good + rel x
#   partial (coverage, retention and adherence have distinct denominators/roles).
#   Burn-in (static) mode: each stage's living total is re-split into undiagnosed /
#   diagnosed / on-BPG labels with the provisional initial shares (labels only; stage
#   totals unchanged; no screening history accrues).
#
# COMPETING RISKS ON THE HAZARD SCALE: annual probabilities -> hazards once,
#   h = -log1p(-p); H = sum h; p_i = (1 - exp(-H)) h_i / H; stay = exp(-H). Exact.
# POPULATION BACKBONE: population exogenous (WPP). Each cycle the No-RHD pools (never-
#   onset + normalised) are reconciled with population - classifiable RHD - prior-A
#   (births at age 0; residual shared proportionally). Disease / prior-A stocks are
#   never rescaled; living disease above population is a measured violation. Burn-in
#   only, explicit: "per_capita" rescales every count stock to the constant burn-in
#   cross-section after each ageing step (initialisation balancing, recorded).
# SURGERY is a SERVICE (population-average reach = requirement x coverage) scaling
#   C->D and D->death on the probability scale; traced, never removes stock.
# ==============================================================================

ABCD_ENGINE_SPEC <- "abcd_first_diagnosis_v1"
ABCD_BURNIN_RECONCILIATION_MODES <- c("per_capita", "no_rhd_only")
ABCD_STAGES             <- c("A", "B", "C", "D")
ABCD_UNDX_NS            <- paste0("U_", ABCD_STAGES, "_ns")
ABCD_UNDX_SC            <- paste0("U_", ABCD_STAGES, "_sc")
ABCD_DX_OFF             <- paste0("X_", ABCD_STAGES)
ABCD_DX_ON              <- c("T_A", "T_B")
ABCD_DISEASE_STOCKS     <- c(ABCD_UNDX_NS, ABCD_UNDX_SC, ABCD_DX_OFF, ABCD_DX_ON)   # classifiable RHD
ABCD_PRIOR_A_STOCKS     <- c("P_ns", "P_sc", "P_dx")
ABCD_NORMALISED_STOCKS  <- c("R_ns", "R_sc", "R_dx")
ABCD_NEVER_ONSET_STOCKS <- c("N_ns", "N_sc")
ABCD_BACKBONE_STOCKS    <- c(ABCD_NEVER_ONSET_STOCKS, ABCD_NORMALISED_STOCKS)        # reconciled with WPP
ABCD_COUNT_STOCKS       <- c(ABCD_DISEASE_STOCKS, ABCD_PRIOR_A_STOCKS, ABCD_NORMALISED_STOCKS, ABCD_NEVER_ONSET_STOCKS)
ABCD_STATE_FIELDS       <- ABCD_COUNT_STOCKS
ABCD_SAP_ELIGIBLE       <- c("A", "B")      # A only exists at age <= child_last_age
ABCD_ONSET_CAP          <- 0.99             # cap on the per-susceptible onset probability (hits counted, never silent)

# ------------------------------------------------------------------------------
# competing-risk helpers
# ------------------------------------------------------------------------------
abcd_p2h   <- function(p) -log1p(-pmin(pmax(p, 0), 1 - 1e-12))
abcd_p2h_v <- function(p) -log1p(-p)          # hot path: inputs already validated in [0, 1)
abcd_split_h <- function(hz) {
  H <- 0
  for (h in hz) H <- H + h
  stay <- exp(-H)
  scale <- (1 - stay) / (H + (H == 0))
  out <- hz
  for (i in seq_along(hz)) out[[i]] <- scale * hz[[i]]
  list(out = out, stay = stay)
}
abcd_compete <- function(probs) abcd_split_h(lapply(probs, abcd_p2h))
abcd_hazard_scale <- function(p, m) 1 - (1 - pmin(pmax(p, 0), 1 - 1e-12))^m
abcd_screen_annual_prob <- function(c, window_years) 1 - (1 - c)^(1 / window_years)

# ------------------------------------------------------------------------------
# NATURAL HISTORY -> engine-ready [age x sex] matrices (built once per parameter set)
#   nh: list(child_last_age, tr = list(p_A_to_no_rhd, p_A_to_B, p_B_to_A, p_B_to_C,
#            p_C_to_B, p_C_to_D, p_D_to_C), p_B_to_no_rhd_adult, p_rhd_death_D,
#            residual_rhd_mortality_mult, entry_mix_child, entry_mix_adult)
# ------------------------------------------------------------------------------
abcd_nh_build <- function(nh, ages, sexes) {
  n_age <- length(ages); n_sex <- length(sexes)
  b <- as.integer(nh$child_last_age)
  b_row <- match(b, ages)
  if (length(b) != 1L || is.na(b_row) || b_row >= n_age)
    stop("abcd_nh_build: child_last_age must be an interior model age.", call. = FALSE)
  tr <- nh$tr
  need <- c("p_A_to_no_rhd", "p_A_to_B", "p_B_to_A", "p_B_to_C", "p_C_to_B", "p_C_to_D", "p_D_to_C")
  if (!all(need %in% names(tr))) stop("abcd_nh_build: tr lacks ", paste(setdiff(need, names(tr)), collapse = ", "), call. = FALSE)
  pr <- unlist(tr[need])
  if (any(!is.finite(pr)) || any(pr < 0 | pr >= 1)) stop("abcd_nh_build: invalid transition probability.", call. = FALSE)
  if (tr$p_D_to_C != 0) stop("abcd_nh_build: D->C must be an exact structural zero.", call. = FALSE)
  mC <- nh$entry_mix_child; mA <- nh$entry_mix_adult
  for (m in list(mC, mA))
    if (!identical(names(m), ABCD_STAGES) || any(m < 0) || abs(sum(m) - 1) > 1e-12)
      stop("abcd_nh_build: entry mixes must be named A/B/C/D, non-negative, summing to 1.", call. = FALSE)
  if (mA[["A"]] != 0) stop("abcd_nh_build: no adult A state (entry_mix_adult['A'] must be 0; adult A is prior-A history).", call. = FALSE)
  if (!isTRUE(all.equal(nh$residual_rhd_mortality_mult, 1)))
    stop("abcd_nh_build: residual RHD mortality multiplier is fixed at 1 in the default.", call. = FALSE)
  for (v in c(nh$p_B_to_no_rhd_adult, nh$p_rhd_death_D))
    if (length(v) != 1L || !is.finite(v) || v < 0 || v >= 1) stop("abcd_nh_build: invalid fixed probability input.", call. = FALSE)
  if (!(nh$p_rhd_death_D > 0)) stop("abcd_nh_build: D_only requires p_rhd_death_D > 0.", call. = FALSE)

  child <- ages <= b; adult <- !child
  M  <- function(v) matrix(v, n_age, n_sex, dimnames = list(ages, sexes))
  ch <- M(as.numeric(child)); ad <- M(as.numeric(adult))
  p <- list(A_to_R = ch * tr$p_A_to_no_rhd,                # childhood resolution to normalised history
            A_to_B = ch * tr$p_A_to_B,
            P_to_B = ad * tr$p_A_to_B,                     # adult prior-A -> B (shared parameter; extrapolation)
            B_to_A = ch * tr$p_B_to_A,
            B_to_C = M(tr$p_B_to_C),                       # one shared parameter at every attained age
            B_to_R = ad * nh$p_B_to_no_rhd_adult,           # fixed adult normalisation
            C_to_B = ch * tr$p_C_to_B,
            C_to_D = M(tr$p_C_to_D))                       # one shared parameter at every attained age
  pd <- nh$p_rhd_death_D * nh$residual_rhd_mortality_mult
  entry <- lapply(setNames(ABCD_STAGES, ABCD_STAGES), function(s) ch * mC[[s]] + ad * mA[[s]])
  list(nh = nh, ages = ages, sexes = sexes, boundary_age = b, boundary_row = b_row,
       child = child, child_m = ch, adult_m = ad,
       p = p, h = lapply(p, abcd_p2h_v), pd = pd, entry = entry,
       mix_child = mC, mix_adult = mA,
       derived = list(p_prior_A_to_B_adult = tr$p_A_to_B, p_B_to_C_adult = tr$p_B_to_C,
                      p_C_to_D_adult = tr$p_C_to_D, p_rhd_death_D = pd))
}

# record of the FIXED natural-history inputs (stale-bundle / spec-consistency check
# 03 -> 04b -> 05: a calibration is consumed only under identical fixed inputs)
abcd_nh_fixed <- function(nh)
  c(nh[c("child_last_age", "p_B_to_no_rhd_adult", "p_rhd_death_D", "residual_rhd_mortality_mult",
         "entry_mix_child", "entry_mix_adult")],
    list(p_D_to_C = nh$tr$p_D_to_C, engine_spec = ABCD_ENGINE_SPEC))

# SAP hazard ratios as [age x sex] matrices (validated): A->B child; B->C child / adult.
# C->D and D->death HRs must be 1 in this default (C/D are not SAP-eligible).
# abcd_eff_prepare() validates once and caches the matrices in eff$.sap_hr (hot path).
abcd_eff_prepare <- function(eff, nhe) { eff$.sap_hr <- abcd_sap_hr(eff, nhe); eff }
abcd_sap_hr <- function(eff, nhe) {
  mech <- if (is.null(eff$sap_mechanism)) "progression" else as.character(eff$sap_mechanism)
  if (mech != "progression") stop("abcd_engine: only sap_mechanism = 'progression' is implemented.", call. = FALSE)
  g <- function(nm) { v <- eff[[nm]]; if (is.null(v) || length(v) != 1L || !is.finite(v) || v < 0 || v > 1)
    stop("abcd_engine: SAP hazard ratio '", nm, "' missing or outside [0, 1].", call. = FALSE); v }
  hAB <- g("sap_hr_A_to_B_child"); hBCc <- g("sap_hr_B_to_C_child"); hBCa <- g("sap_hr_B_to_C_adult")
  for (nm in c("sap_hr_C_to_D", "sap_hr_D_to_rhd_death", "sap_hr_prior_A_to_B_adult"))
    if (!isTRUE(all.equal(g(nm), 1)))
      stop("abcd_engine: ", nm, " must be 1 (C, D and adult prior-A are not SAP-eligible in this default).", call. = FALSE)
  list(A_to_B = nhe$child_m * hAB + nhe$adult_m,
       B_to_C = nhe$child_m * hBCc + nhe$adult_m * hBCa)
}

# ------------------------------------------------------------------------------
# CARE schedules
#   care_in: per-year vectors screen_window_prob, clinical_detection, bpg_initiation,
#            bpg_retention, adherence_good, adherence_partial, surgery (length n or 1)
#            + scalars sensitivity, confirm_completion, confirm_sensitivity,
#            partial_relative_protection, screen_age_lo, screen_age_hi.
# ------------------------------------------------------------------------------
abcd_care_schedule <- function(care_in, n, ages, sexes) {
  n_age <- length(ages); n_sex <- length(sexes)
  vec_nm <- c("screen_window_prob", "clinical_detection", "bpg_initiation", "bpg_retention",
              "adherence_good", "adherence_partial", "surgery")
  for (nm in vec_nm) {
    v <- care_in[[nm]]
    if (is.null(v)) stop("abcd_care_schedule: missing care input '", nm, "'.", call. = FALSE)
    if (length(v) == 1L) v <- rep(v, n)
    if (length(v) != n || any(!is.finite(v)) || any(v < 0 | v > 1))
      stop("abcd_care_schedule: care input '", nm, "' must be ", n, " probabilities in [0,1].", call. = FALSE)
    care_in[[nm]] <- v
  }
  if (any(care_in$adherence_good + care_in$adherence_partial > 1 + 1e-12))
    stop("abcd_care_schedule: good + partial adherence shares exceed 1.", call. = FALSE)
  for (nm in c("sensitivity", "confirm_completion", "confirm_sensitivity", "partial_relative_protection")) {
    v <- care_in[[nm]]
    if (is.null(v) || length(v) != 1L || !is.finite(v) || v < 0 || v > 1)
      stop("abcd_care_schedule: '", nm, "' must be one value in [0,1].", call. = FALSE)
  }
  lo <- as.integer(care_in$screen_age_lo); hi <- as.integer(care_in$screen_age_hi)
  if (lo > hi || lo < min(ages) || hi > max(ages)) stop("abcd_care_schedule: invalid screening window.", call. = FALSE)
  W <- hi - lo + 1L
  mask <- matrix(as.numeric(ages >= lo & ages <= hi), n_age, n_sex, dimnames = list(ages, sexes))
  lapply(seq_len(n), function(t) list(
    dynamic = TRUE,
    window_prob = care_in$screen_window_prob[t],
    q = mask * abcd_screen_annual_prob(care_in$screen_window_prob[t], W),
    sens = care_in$sensitivity, cc = care_in$confirm_completion, cs = care_in$confirm_sensitivity,
    k = care_in$clinical_detection[t],
    init = care_in$bpg_initiation[t], retain = care_in$bpg_retention[t],
    adh_factor = care_in$adherence_good[t] + care_in$partial_relative_protection * care_in$adherence_partial[t],
    surgery = care_in$surgery[t]))
}

# burn-in care: provisional initial shares held as static pre-history care
abcd_care_static <- function(d0, b0, adh_factor, surgery) {
  if (!(d0 >= 0 && d0 <= 1 && b0 >= 0 && b0 <= d0)) stop("abcd_care_static: need 0 <= b0 <= d0 <= 1.", call. = FALSE)
  list(dynamic = FALSE, d0 = d0, b0 = b0, adh_factor = adh_factor, surgery = surgery,
       q = 0, sens = 0, cc = 0, cs = 0, k = 0, init = 0, retain = 0, window_prob = 0)
}

# allocate care LABELS within each stage (stage totals unchanged): undiagnosed (1 - d0),
# diagnosed d0 of which on BPG b0 (of the stage total) in SAP-eligible stages (A at
# age <= child_last_age, B). Undiagnosed are placed in the never-screened stock.
abcd_care_allocate <- function(st, d0, b0, nhe) {
  ch <- nhe$child_m
  tA <- st$U_A_ns + st$U_A_sc + st$X_A + st$T_A; tB <- st$U_B_ns + st$U_B_sc + st$X_B + st$T_B
  tC <- st$U_C_ns + st$U_C_sc + st$X_C;           tD <- st$U_D_ns + st$U_D_sc + st$X_D
  z <- tA * 0
  st$U_A_ns <- tA * (1 - d0); st$U_A_sc <- z; st$X_A <- tA * (d0 - b0 * ch); st$T_A <- tA * b0 * ch
  st$U_B_ns <- tB * (1 - d0); st$U_B_sc <- z; st$X_B <- tB * (d0 - b0);      st$T_B <- tB * b0
  st$U_C_ns <- tC * (1 - d0); st$U_C_sc <- z; st$X_C <- tC * d0
  st$U_D_ns <- tD * (1 - d0); st$U_D_sc <- z; st$X_D <- tD * d0
  st
}

# initialise the dynamic care state ONCE (first observed-history year): allocate the
# provisional diagnosed / on-BPG shares to eligible stage stocks and the initial
# ever-screened share among people without a prior diagnosis
abcd_care_initialize <- function(st, d0, b0, ever_screened0, nhe) {
  if (!(d0 >= 0 && d0 <= 1 && b0 >= 0 && b0 <= d0 && ever_screened0 >= 0 && ever_screened0 <= 1))
    stop("abcd_care_initialize: need 0 <= b0 <= d0 <= 1 and ever_screened0 in [0, 1].", call. = FALSE)
  st <- abcd_care_allocate(st, d0, b0, nhe)
  mv <- function(st, ns, sc) { x <- st[[ns]] * ever_screened0; st[[ns]] <- st[[ns]] - x; st[[sc]] <- st[[sc]] + x; st }
  for (p in list(c("N_ns", "N_sc"), c("P_ns", "P_sc"), c("R_ns", "R_sc")))
    st <- mv(st, p[1], p[2])
  for (s in ABCD_STAGES) st <- mv(st, paste0("U_", s, "_ns"), paste0("U_", s, "_sc"))
  st
}

# ------------------------------------------------------------------------------
# state helpers
# ------------------------------------------------------------------------------
abcd_state_zero <- function(template) {
  z <- template * 0
  setNames(rep(list(z), length(ABCD_COUNT_STOCKS)), ABCD_COUNT_STOCKS)
}
abcd_stage_total <- function(st, s) switch(s,
  A = st$U_A_ns + st$U_A_sc + st$X_A + st$T_A, B = st$U_B_ns + st$U_B_sc + st$X_B + st$T_B,
  C = st$U_C_ns + st$U_C_sc + st$X_C,          D = st$U_D_ns + st$U_D_sc + st$X_D)
abcd_rhd_total <- function(st)              # classifiable RHD (GBD prevalence mapping)
  st$U_A_ns + st$U_B_ns + st$U_C_ns + st$U_D_ns + st$U_A_sc + st$U_B_sc + st$U_C_sc + st$U_D_sc +
  st$X_A + st$X_B + st$X_C + st$X_D + st$T_A + st$T_B
abcd_undiagnosed_total <- function(st)
  st$U_A_ns + st$U_B_ns + st$U_C_ns + st$U_D_ns + st$U_A_sc + st$U_B_sc + st$U_C_sc + st$U_D_sc
abcd_diagnosed_total <- function(st) st$X_A + st$X_B + st$X_C + st$X_D + st$T_A + st$T_B
abcd_prior_A_total <- function(st) st$P_ns + st$P_sc + st$P_dx
abcd_normalised_total <- function(st) st$R_ns + st$R_sc + st$R_dx
abcd_never_onset_total <- function(st) st$N_ns + st$N_sc
abcd_stage_aggregate <- function(st) lapply(setNames(ABCD_STAGES, ABCD_STAGES), function(s) abcd_stage_total(st, s))
abcd_living_total <- function(st) { x <- 0; for (k in ABCD_COUNT_STOCKS) x <- x + st[[k]]; x }

# one undiagnosed screening-history group (never screened "ns" / screened "sc"): the
# within-cycle stage flows (top-level so it is byte-compiled once; hot path)
abcd_undx_next <- function(UA, UB, UC, UD, P, R, on, sA, sB, sC, sD, sP, stay_bg, E)
  list(A = UA * sA$stay + UB * sB$out$toA + on * E$A,
       B = UB * sB$stay + UA * sA$out$toB + UC * sC$out$toB + P * sP$out$toB + on * E$B,
       C = UC * sC$stay + UB * sB$out$toC + on * E$C,
       D = UD * sD$stay + UC * sC$out$toD + on * E$D,
       P = P * sP$stay,
       R = R * stay_bg + UA * sA$out$toR + UB * sB$out$toR)

# per-capita burn-in balancing (initialisation only): rescale every count stock to pop1
abcd_burnin_balance <- function(st, pop1) {
  carried <- abcd_living_total(st)
  nc <- carried == 0
  f1 <- (!nc) * pop1 / (carried + nc) + nc - 1
  dis <- abcd_rhd_total(st); pa <- abcd_prior_A_total(st); nr <- abcd_normalised_total(st); nv <- abcd_never_onset_total(st)
  rec <- c(disease = sum(dis * f1), prior_A = sum(pa * f1), normalised = sum(nr * f1), never_onset = sum(nv * f1),
           carried_total = sum(carried[!nc]), backbone_total = sum(pop1[!nc]),
           disease_rel_adj = if (sum(dis) > 0) sum(dis * f1) / sum(dis) else 0)
  f <- f1 + 1
  for (k in ABCD_COUNT_STOCKS) st[[k]] <- st[[k]] * f
  list(st = st, rec = rec)
}

# ------------------------------------------------------------------------------
# ONE ANNUAL CYCLE
#   state : list of [age x sex] matrices (ABCD_STATE_FIELDS) ENTERING the cycle
#   popm  : population this year (exogenous backbone)
#   irm   : incidence per TOTAL population (calibrated); converted per never-onset person
#   oth   : background (non-RHD) mortality probability
#   nhe   : abcd_nh_build() output ; eff : SAP HRs + surgery RRRs ; surg : requirement fractions
#   care  : one element of abcd_care_schedule() (or abcd_care_static())
# Returns the next-cycle state (PRE-ageing) + flows. light = TRUE returns only what the
# calibration scores; the arithmetic is the same code path (identical numbers).
# ------------------------------------------------------------------------------
abcd_one_cycle <- function(state, popm, irm, oth, nhe, eff, surg, care, light = FALSE) {
  s <- state
  dynamic <- isTRUE(care$dynamic)
  if (!dynamic) s <- abcd_care_allocate(s, care$d0, care$b0, nhe)      # static pre-history labels (burn-in)
  rhd_start <- abcd_rhd_total(s)
  P_start   <- s$P_ns + s$P_sc + s$P_dx

  # --- 0. reconcile the No-RHD pools (never-onset + normalised) with the backbone ---
  nr_backbone <- popm - rhd_start - P_start
  neg_nr <- nr_backbone < 0
  pop_deficit <- -nr_backbone * neg_nr
  nr_backbone <- nr_backbone * !neg_nr
  N0 <- s$N_ns + s$N_sc; R0 <- s$R_ns + s$R_sc + s$R_dx
  nr_model <- N0 + R0
  no_carry <- nr_model == 0
  f_dem <- (!no_carry) * nr_backbone / (nr_model + no_carry)
  for (k in ABCD_BACKBONE_STOCKS) s[[k]] <- s[[k]] * f_dem
  s$N_ns <- s$N_ns + no_carry * nr_backbone                           # births (age 0): never onset, never screened
  demog_adj_N <- (s$N_ns + s$N_sc) - N0
  demog_adj_R <- (s$R_ns + s$R_sc + s$R_dx) - R0
  rec_start <- s$R_ns + s$R_sc + s$R_dx

  # --- 1. care update (start of cycle) -----------------------------------------
  z <- popm * 0
  dx_scr <- z; dx_clin <- z; tp <- z; conf_done <- z; encounters <- z
  dx_stage <- list(A = z, B = z, C = z, D = z)
  init_flow <- z; dropouts <- z; q <- 0
  if (dynamic) {
    q <- care$q; k <- care$k
    p_tp <- q * care$sens; p_cd <- p_tp * care$cc; p_dx <- p_cd * care$cs
    U_ns_tot <- z
    for (st_ in ABCD_STAGES) {
      ns <- paste0("U_", st_, "_ns"); sc <- paste0("U_", st_, "_sc"); xx <- paste0("X_", st_)
      Uns <- s[[ns]]; Usc <- s[[sc]]
      U_ns_tot <- U_ns_tot + Uns
      d_scr  <- Uns * p_dx                                  # first confirmed diagnosis via screening
      missed <- Uns * q - d_scr                             # screened, not confirmed: stays undiagnosed
      unscr  <- Uns - Uns * q
      d_clin <- k * (missed + unscr + Usc)                  # confirmed clinical diagnosis (no overlap)
      s[[ns]] <- unscr * (1 - k)
      s[[sc]] <- (missed + Usc) * (1 - k)
      s[[xx]] <- s[[xx]] + d_scr + d_clin
      dx_scr <- dx_scr + d_scr; dx_clin <- dx_clin + d_clin
      dx_stage[[st_]] <- d_scr + d_clin
      tp <- tp + Uns * p_tp; conf_done <- conf_done + Uns * p_cd
    }
    encounters <- q * (s$N_ns + U_ns_tot + s$P_ns + s$R_ns)    # people without a prior confirmed diagnosis
    for (p in list(c("N_ns", "N_sc"), c("P_ns", "P_sc"), c("R_ns", "R_sc"))) {
      x <- s[[p[1]]] * q; s[[p[1]]] <- s[[p[1]]] - x; s[[p[2]]] <- s[[p[2]]] + x }
    # BPG: retention, dropouts back to diagnosed-off-BPG, initiation among eligible diagnosed
    for (st_ in ABCD_SAP_ELIGIBLE) {
      tt <- paste0("T_", st_); xx <- paste0("X_", st_)
      elig <- if (st_ == "A") nhe$child_m else 1
      drop <- s[[tt]] * (1 - care$retain)
      s[[tt]] <- s[[tt]] - drop; s[[xx]] <- s[[xx]] + drop
      ini <- s[[xx]] * care$init * elig
      s[[xx]] <- s[[xx]] - ini; s[[tt]] <- s[[tt]] + ini
      init_flow <- init_flow + ini; dropouts <- dropouts + drop
    }
  }
  care_violation <- max(0, -min(unlist(s, use.names = FALSE)) - 1e-9)   # any negative stock after the care update
  adh <- care$adh_factor
  hr <- if (is.null(eff$.sap_hr)) abcd_sap_hr(eff, nhe) else eff$.sap_hr
  mAB <- 1 - adh * (1 - hr$A_to_B); mBC <- 1 - adh * (1 - hr$B_to_C)

  # --- 2. surgery service ------------------------------------------------------
  reach_C <- surg$frac_C_requiring_surgery * care$surgery
  reach_D <- surg$frac_D_requiring_surgery * care$surgery
  eff_C_to_D_surg  <- 1 - eff$eff_surgery_C_to_D         * reach_C
  eff_D_death_surg <- 1 - eff$eff_surgery_D_to_rhd_death * reach_D

  # --- 3. never-onset pool: first modelled disease vs background death ---------
  h_oth <- abcd_p2h_v(oth)
  N_tot <- s$N_ns + s$N_sc
  p_on  <- (N_tot > 0) * irm * popm / (N_tot + (N_tot == 0))
  onset_cap_hits <- sum(p_on > ABCD_ONSET_CAP)
  if (onset_cap_hits > 0) p_on <- pmin(p_on, ABCD_ONSET_CAP)
  sN <- abcd_split_h(list(on = abcd_p2h_v(p_on), oth = h_oth))
  on_ns <- s$N_ns * sN$out$on; on_sc <- s$N_sc * sN$out$on
  onset <- on_ns + on_sc
  N_deaths <- N_tot * sN$out$oth
  stay_bg <- exp(-h_oth)

  # --- 4. stage dynamics (untreated / on-BPG splits; attained-age rules) ---------
  h <- nhe$h
  hCD <- abcd_p2h_v(nhe$p$C_to_D * eff_C_to_D_surg)
  hD  <- abcd_p2h_v(nhe$pd * eff_D_death_surg)
  sA  <- abcd_split_h(list(toR = h$A_to_R, toB = h$A_to_B, oth = h_oth))
  sAt <- abcd_split_h(list(toR = h$A_to_R, toB = h$A_to_B * mAB, oth = h_oth))
  sB  <- abcd_split_h(list(toA = h$B_to_A, toC = h$B_to_C, toR = h$B_to_R, oth = h_oth))
  sBt <- abcd_split_h(list(toA = h$B_to_A, toC = h$B_to_C * mBC, toR = h$B_to_R, oth = h_oth))
  sC  <- abcd_split_h(list(toB = h$C_to_B, toD = hCD, oth = h_oth))
  sD  <- abcd_split_h(list(rhd = hD, oth = h_oth))
  sP  <- abcd_split_h(list(toB = h$P_to_B, oth = h_oth))

  E <- nhe$entry
  nxt <- s
  # undiagnosed (never screened / screened) -- identical rules, own screening history
  u <- abcd_undx_next(s$U_A_ns, s$U_B_ns, s$U_C_ns, s$U_D_ns, s$P_ns, s$R_ns, on_ns, sA, sB, sC, sD, sP, stay_bg, E)
  nxt$U_A_ns <- u$A; nxt$U_B_ns <- u$B; nxt$U_C_ns <- u$C; nxt$U_D_ns <- u$D; nxt$P_ns <- u$P; nxt$R_ns <- u$R
  u <- abcd_undx_next(s$U_A_sc, s$U_B_sc, s$U_C_sc, s$U_D_sc, s$P_sc, s$R_sc, on_sc, sA, sB, sC, sD, sP, stay_bg, E)
  nxt$U_A_sc <- u$A; nxt$U_B_sc <- u$B; nxt$U_C_sc <- u$C; nxt$U_D_sc <- u$D; nxt$P_sc <- u$P; nxt$R_sc <- u$R
  # diagnosed (off BPG / on BPG)
  nxt$X_A <- s$X_A * sA$stay + s$X_B * sB$out$toA
  nxt$X_B <- s$X_B * sB$stay + s$X_A * sA$out$toB + s$X_C * sC$out$toB + s$P_dx * sP$out$toB
  nxt$X_C <- s$X_C * sC$stay + s$X_B * sB$out$toC + s$T_B * sBt$out$toC     # C is not SAP-eligible: BPG ends
  nxt$X_D <- s$X_D * sD$stay + s$X_C * sC$out$toD
  nxt$T_A <- s$T_A * sAt$stay + s$T_B * sBt$out$toA
  nxt$T_B <- s$T_B * sBt$stay + s$T_A * sAt$out$toB
  nxt$P_dx <- s$P_dx * sP$stay
  nxt$R_dx <- s$R_dx * stay_bg + s$X_A * sA$out$toR + s$T_A * sAt$out$toR + s$X_B * sB$out$toR + s$T_B * sBt$out$toR
  nxt$N_ns <- s$N_ns * sN$stay; nxt$N_sc <- s$N_sc * sN$stay

  UD_all <- s$U_D_ns + s$U_D_sc
  rhd_deaths <- (UD_all + s$X_D) * sD$out$rhd
  sick_end <- abcd_rhd_total(nxt)

  if (isTRUE(light))
    return(list(state = nxt, sick_end = sick_end, rhd_deaths = rhd_deaths, new_rhd_total = onset,
                pop_deficit = pop_deficit, care_violation = care_violation, onset_cap_hits = onset_cap_hits,
                demog_adj_births = sum(demog_adj_N[1, ]), demog_adj_N = sum(demog_adj_N[-1, ]),
                demog_adj_R = sum(demog_adj_R)))

  # ------------------------------ full output ----------------------------------
  A0 <- s$U_A_ns + s$U_A_sc + s$X_A; B0 <- s$U_B_ns + s$U_B_sc + s$X_B
  C0 <- s$U_C_ns + s$U_C_sc + s$X_C; D0 <- UD_all + s$X_D
  P0 <- s$P_ns + s$P_sc + s$P_dx
  R_deaths <- rec_start * (1 - stay_bg)
  other_A <- A0 * sA$out$oth + s$T_A * sAt$out$oth
  other_B <- B0 * sB$out$oth + s$T_B * sBt$out$oth
  other_C <- C0 * sC$out$oth; other_D <- D0 * sD$out$oth
  P_deaths <- P0 * sP$out$oth
  other_rhd <- other_A + other_B + other_C + other_D
  A_to_R <- A0 * sA$out$toR + s$T_A * sAt$out$toR
  B_to_R <- B0 * sB$out$toR + s$T_B * sBt$out$toR
  all_mx <- rhd_deaths + other_rhd + P_deaths + N_deaths + R_deaths
  C_req_surg <- C0 * surg$frac_C_requiring_surgery
  D_req_surg <- D0 * surg$frac_D_requiring_surgery
  dx_after <- abcd_diagnosed_total(s); tx_after <- s$T_A + s$T_B
  rhd0 <- rhd_start == 0
  no_prior_dx <- s$N_ns + s$N_sc + abcd_undiagnosed_total(s) + s$P_ns + s$P_sc + s$R_ns + s$R_sc
  ns_after <- s$N_ns + s$U_A_ns + s$U_B_ns + s$U_C_ns + s$U_D_ns + s$P_ns + s$R_ns
  agg_n <- abcd_stage_aggregate(nxt)
  list(
    state = nxt,
    A = agg_n$A, B = agg_n$B, C = agg_n$C, D = agg_n$D,
    sick_end = sick_end, well_end = pmax(popm - sick_end, 0),
    undiagnosed_end = abcd_undiagnosed_total(nxt), diagnosed_end = abcd_diagnosed_total(nxt),
    on_bpg_end = nxt$T_A + nxt$T_B, prior_A_end = abcd_prior_A_total(nxt), normalised_end = abcd_normalised_total(nxt),
    rhd_start = rhd_start, prior_A_start = P_start, never_onset_start = N_tot,
    never_onset_end = nxt$N_ns + nxt$N_sc, recovered_start = rec_start,
    pop_deficit = pop_deficit,
    demog_adjust_never_onset = demog_adj_N, demog_adjust_recovered = demog_adj_R, demog_scale_no_rhd = f_dem,
    demog_adj_births = sum(demog_adj_N[1, ]), demog_adj_N = sum(demog_adj_N[-1, ]), demog_adj_R = sum(demog_adj_R),
    # NEW MODELLED DISEASE (first classifiable entry; never diagnoses) by entry stage
    new_rhd_A = onset * E$A, new_rhd_B = onset * E$B, new_rhd_C = onset * E$C, new_rhd_D = onset * E$D,
    new_rhd_total = onset, onset_prob_susceptible = p_on, onset_cap_hits = onset_cap_hits,
    # re-entry of adult prior-A into classifiable B (separate from first entry / first diagnosis)
    reentry_prior_A_to_B = P0 * sP$out$toB, reentry_prior_A_to_B_diagnosed_history = s$P_dx * sP$out$toB,
    # transitions (all care states)
    A_to_no_rhd = A_to_R, B_to_no_rhd = B_to_R,
    A_to_B = A0 * sA$out$toB + s$T_A * sAt$out$toB, B_to_A = B0 * sB$out$toA + s$T_B * sBt$out$toA,
    B_to_C = B0 * sB$out$toC + s$T_B * sBt$out$toC, C_to_B = C0 * sC$out$toB, C_to_D = C0 * sC$out$toD, D_to_C = z,
    rhd_deaths_A = z, rhd_deaths_B = z, rhd_deaths_C = z, rhd_deaths_D = rhd_deaths,
    rhd_deaths_D_undiagnosed = UD_all * sD$out$rhd, rhd_deaths_D_diagnosed = s$X_D * sD$out$rhd,
    other_deaths_A = other_A, other_deaths_B = other_B, other_deaths_C = other_C, other_deaths_D = other_D,
    rhd_deaths = rhd_deaths, other_deaths = other_rhd, rhd_exits_total = A_to_R + B_to_R + rhd_deaths + other_rhd,
    prior_A_deaths = P_deaths, never_onset_deaths = N_deaths, recovered_deaths = R_deaths,
    well_deaths = N_deaths + R_deaths, all_mx = all_mx,
    D_start = D0,
    # surgery TRACE (service; stocks after the care update = start-of-cycle disease; never subtracted)
    C_requiring_surgery = C_req_surg, D_requiring_surgery = D_req_surg,
    surgeries_C = C_req_surg * care$surgery, surgeries_D = D_req_surg * care$surgery,
    total_surgeries = (C_req_surg + D_req_surg) * care$surgery,
    # care volumes (start of cycle, after the care update)
    screen_annual_prob = q + z, n_screened = encounters,
    never_screened_share = ns_after / (no_prior_dx + (no_prior_dx == 0)) + (no_prior_dx == 0),
    n_screen_true_positive = tp, n_confirm_completed = conf_done,
    n_dx_screen = dx_scr, n_dx_clinical = dx_clin, n_new_diagnoses = dx_scr + dx_clin,
    n_screen_detected = dx_scr,
    n_dx_A = dx_stage$A, n_dx_B = dx_stage$B, n_dx_C = dx_stage$C, n_dx_D = dx_stage$D,
    n_diagnosed = dx_after, n_undiagnosed = abcd_undiagnosed_total(s),
    n_diagnosed_A = s$X_A + s$T_A, n_diagnosed_B = s$X_B + s$T_B, n_diagnosed_C = s$X_C, n_diagnosed_D = s$X_D,
    n_on_optimal_treatment = tx_after, n_on_bpg_A = s$T_A, n_on_bpg_B = s$T_B,
    n_bpg_eligible_diagnosed = s$X_A * nhe$child_m + s$X_B + tx_after,
    n_bpg_initiations = init_flow, n_bpg_dropouts = dropouts,
    n_effectively_protected = tx_after * adh,
    diagnosed_share = dx_after / (rhd_start + rhd0), on_bpg_share = tx_after / (rhd_start + rhd0),
    care_violation = care_violation,
    sap_mult = 1, sap_mechanism = "progression", sap_c_eff = tx_after * adh / (rhd_start + rhd0),
    sap_prog_mult_A_to_B = mAB, sap_prog_mult_B_to_C = mBC, sap_prog_mult_C_to_D = 1 + z,
    adherence_factor = adh,
    reach_C = reach_C, reach_D = reach_D, eff_ir = 1,
    eff_C_to_D_surg = eff_C_to_D_surg, eff_D_death_surg = eff_D_death_surg)
}

# ------------------------------------------------------------------------------
# AGEING: advance every cohort one year; the terminal row is the OPEN group.
# ------------------------------------------------------------------------------
abcd_age_shift <- function(M) {
  n <- nrow(M)
  N <- M * 0
  N[2:n, ] <- M[1:(n - 1), , drop = FALSE]
  N[n, ]   <- N[n, ] + M[n, ]
  N
}

# age the whole state; at the child_last_age -> +1 step surviving A becomes adult
# prior-A (undiagnosed -> P_ns / P_sc by screening history; diagnosed incl. on BPG ->
# P_dx; BPG ends). Returns the reclassified counts at the boundary age (people conserved).
abcd_age_state <- function(st, boundary_row) {
  out <- st
  for (nm in ABCD_COUNT_STOCKS) out[[nm]] <- abcd_age_shift(st[[nm]])
  b2 <- boundary_row + 1L
  rc_ns <- out$U_A_ns[b2, ]; rc_sc <- out$U_A_sc[b2, ]; rc_x <- out$X_A[b2, ]; rc_t <- out$T_A[b2, ]
  out$P_ns[b2, ] <- out$P_ns[b2, ] + rc_ns; out$U_A_ns[b2, ] <- 0
  out$P_sc[b2, ] <- out$P_sc[b2, ] + rc_sc; out$U_A_sc[b2, ] <- 0
  out$P_dx[b2, ] <- out$P_dx[b2, ] + rc_x + rc_t; out$X_A[b2, ] <- 0; out$T_A[b2, ] <- 0
  z <- st$U_A_ns * 0
  m_u <- z; m_u[boundary_row, ] <- rc_ns + rc_sc
  m_d <- z; m_d[boundary_row, ] <- rc_x + rc_t
  m_t <- z; m_t[boundary_row, ] <- rc_t
  list(state = out, reclass_undiagnosed = m_u, reclass_diagnosed = m_d, reclass_on_bpg = m_t,
       boundary_transfer = m_u + m_d)
}

# ------------------------------------------------------------------------------
# BURN-IN under CONSTANT inputs (first-window-year population, incidence, background
# mortality, static pre-history care). Convergence is CHECKED (stage / prior-A /
# normalised totals over the last cycle and vs the half-way state, age profile of
# classifiable RHD, childhood and adult RHD deaths).
#   reconcile = "per_capita" : after every burn-in ageing step EVERY count stock is
#     rescaled by f = pop1 / carried living total (per-capita equilibrium on the one
#     constant cross-section); adjustments = INITIALISATION BALANCING ($init_balancing),
#     never deaths or flows. "no_rhd_only": only the No-RHD pools absorb the residual.
# Observed years, bridge and production never rescale disease stocks.
# ------------------------------------------------------------------------------
abcd_burn_in <- function(pop1, ir1, oth1, nhe, eff, surg, care1, n_cycles, reconcile, init = NULL) {
  if (missing(reconcile) || length(reconcile) != 1L || !reconcile %in% ABCD_BURNIN_RECONCILIATION_MODES)
    stop("abcd_burn_in: `reconcile` must be one of ", paste(ABCD_BURNIN_RECONCILIATION_MODES, collapse = ", "),
         " (set STRUCT_BURNIN_RECONCILIATION in 00 block G2).", call. = FALSE)
  st <- if (is.null(init)) abcd_state_zero(pop1) else init
  eff <- abcd_eff_prepare(eff, nhe)
  half <- NULL; prev <- NULL; d_prev <- NULL; d_last <- NULL
  per_capita <- reconcile == "per_capita"
  bal <- if (per_capita) matrix(NA_real_, n_cycles, 7L) else NULL
  for (b in seq_len(n_cycles)) {
    cyc <- abcd_one_cycle(st, pop1, ir1, oth1, nhe, eff, surg, care1, light = TRUE)
    if (b == n_cycles) { prev <- st; d_last <- cyc$rhd_deaths }
    if (b == n_cycles - 1L) d_prev <- cyc$rhd_deaths
    st <- abcd_age_state(cyc$state, nhe$boundary_row)$state
    if (per_capita) { bb <- abcd_burnin_balance(st, pop1); st <- bb$st; bal[b, ] <- bb$rec; cn <- names(bb$rec) }
    if (b == max(1L, n_cycles %/% 2L)) half <- st
  }
  if (per_capita) bal <- data.frame(cycle = seq_len(n_cycles), setNames(as.data.frame(bal), cn))
  rel <- function(x, y) { sx <- sum(x); if (sx <= 0) 0 else abs(sx - sum(y)) / sx }
  tot <- function(st) c(lapply(abcd_stage_aggregate(st), identity),
                        list(prior_A = abcd_prior_A_total(st), normalised = abcd_normalised_total(st)))
  t_end <- tot(st); t_prev <- tot(prev); t_half <- tot(half)
  conv <- data.frame(stock = names(t_end),
                     rel_change_last_cycle = vapply(names(t_end), function(k) rel(t_end[[k]], t_prev[[k]]), numeric(1)),
                     rel_change_vs_half    = vapply(names(t_end), function(k) rel(t_end[[k]], t_half[[k]]), numeric(1)),
                     row.names = NULL)
  ch <- nhe$ages < 20
  prof <- abcd_rhd_total(st); prof_prev <- abcd_rhd_total(prev)
  big <- prof > 1e-6 * max(sum(prof), 1)
  extra <- data.frame(
    stock = c("age_profile_max_rel_change_last_cycle", "child_deaths_rel_change_last_cycle",
              "adult_deaths_rel_change_last_cycle"),
    rel_change_last_cycle = c(if (any(big)) max(abs(prof[big] - prof_prev[big]) / prof[big]) else 0,
                              rel(d_last[ch, , drop = FALSE], if (is.null(d_prev)) d_last[ch, , drop = FALSE] else d_prev[ch, , drop = FALSE]),
                              rel(d_last[!ch, , drop = FALSE], if (is.null(d_prev)) d_last[!ch, , drop = FALSE] else d_prev[!ch, , drop = FALSE])),
    rel_change_vs_half = NA_real_)
  list(state = st, convergence = rbind(conv, extra), n_cycles = n_cycles, reconcile = reconcile,
       init_balancing = if (per_capita) bal else NULL)
}

# ------------------------------------------------------------------------------
# PROJECTION WRAPPER shared by calibration (light = TRUE) and production (light =
# FALSE). ONE time convention everywhere:
#   start state (entering year t) -> No-RHD reconciliation -> care update -> within-
#   year competing events -> new disease entries at year end -> end state (pre-ageing)
#   -> ageing (+ A -> prior-A reclassification at the boundary) -> state entering t+1.
# Prevalence observation (prev_obs): "start", "end" (pre-ageing) or "mid" = mean, of
# CLASSIFIABLE RHD (prior-A and normalised history excluded).
# ------------------------------------------------------------------------------
abcd_project <- function(init, pop, ir, oth, nhe, eff, surg, care_year, years,
                         light = TRUE, prev_obs = c("mid", "start", "end")) {
  prev_obs <- match.arg(prev_obs)
  n <- length(years)
  stopifnot(dim(pop)[3] == n, dim(ir)[3] == n, length(care_year) == n)
  oth_is_arr <- length(dim(oth)) == 3L
  eff <- abcd_eff_prepare(eff, nhe)
  st <- init
  starts <- vector("list", n); cycles <- vector("list", n); pobs <- vector("list", n)
  reclass <- vector("list", n)
  deficit <- 0; care_viol <- 0; cap_hits <- 0; deficit_cells <- list()
  demog <- data.frame(year = years, births = NA_real_, adjust_never_onset_age1plus = NA_real_, adjust_recovered = NA_real_)
  for (t in seq_len(n)) {
    popm <- pop[, , t]; irm <- ir[, , t]
    othm <- if (oth_is_arr) oth[, , t] else oth
    if (is.null(dim(popm))) { dm <- dim(pop)[1:2]; dim(popm) <- dm; dim(irm) <- dm; dim(othm) <- dm }
    cyc <- abcd_one_cycle(st, popm, irm, othm, nhe, eff, surg, care_year[[t]], light = light)
    starts[[t]] <- st
    cycles[[t]] <- cyc
    start_stock <- abcd_rhd_total(st)
    pobs[[t]] <- switch(prev_obs, start = start_stock, end = cyc$sick_end,
                        mid = 0.5 * (start_stock + cyc$sick_end))
    deficit <- deficit + sum(cyc$pop_deficit)
    if (any(cyc$pop_deficit > 0)) {
      w <- which(cyc$pop_deficit > 0, arr.ind = TRUE)
      deficit_cells[[length(deficit_cells) + 1L]] <- data.frame(year = years[t], age = nhe$ages[w[, 1]], sex = nhe$sexes[w[, 2]],
                                                                deficit = cyc$pop_deficit[w], population = popm[w])
    }
    care_viol <- care_viol + cyc$care_violation
    cap_hits <- cap_hits + cyc$onset_cap_hits
    demog$births[t] <- cyc$demog_adj_births; demog$adjust_never_onset_age1plus[t] <- cyc$demog_adj_N
    demog$adjust_recovered[t] <- cyc$demog_adj_R
    aged <- abcd_age_state(cyc$state, nhe$boundary_row)
    reclass[[t]] <- aged[c("reclass_undiagnosed", "reclass_diagnosed", "reclass_on_bpg")]
    st <- aged$state
  }
  list(start = starts, cycles = cycles, prev_obs = pobs, reclass = reclass,
       boundary_transfer = lapply(reclass, function(r) r$reclass_undiagnosed + r$reclass_diagnosed),
       final_state = st, pop_deficit_total = deficit,
       pop_deficit_cells = if (length(deficit_cells)) do.call(rbind, deficit_cells) else NULL,
       care_violation_total = care_viol,
       onset_cap_hits = cap_hits, demographic_residual = demog,
       years = years, prev_obs_convention = prev_obs)
}

ABCD_ENGINE_FUNS <- c("abcd_compete", "abcd_split_h", "abcd_p2h", "abcd_p2h_v", "abcd_hazard_scale",
                      "abcd_screen_annual_prob", "abcd_nh_build", "abcd_nh_fixed", "abcd_sap_hr", "abcd_eff_prepare",
                      "abcd_undx_next", "abcd_burnin_balance",
                      "abcd_care_schedule", "abcd_care_static", "abcd_care_allocate", "abcd_care_initialize",
                      "abcd_state_zero", "abcd_stage_total", "abcd_rhd_total", "abcd_undiagnosed_total",
                      "abcd_diagnosed_total", "abcd_prior_A_total", "abcd_normalised_total", "abcd_never_onset_total",
                      "abcd_stage_aggregate", "abcd_living_total", "abcd_one_cycle", "abcd_age_shift",
                      "abcd_age_state", "abcd_burn_in", "abcd_project",
                      "ABCD_ENGINE_SPEC", "ABCD_BURNIN_RECONCILIATION_MODES", "ABCD_STAGES", "ABCD_UNDX_NS",
                      "ABCD_UNDX_SC", "ABCD_DX_OFF", "ABCD_DX_ON", "ABCD_DISEASE_STOCKS", "ABCD_PRIOR_A_STOCKS",
                      "ABCD_NORMALISED_STOCKS", "ABCD_NEVER_ONSET_STOCKS", "ABCD_BACKBONE_STOCKS",
                      "ABCD_COUNT_STOCKS", "ABCD_STATE_FIELDS", "ABCD_SAP_ELIGIBLE", "ABCD_ONSET_CAP")
