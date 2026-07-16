# ==============================================================================
# RHD secondary-prevention investment case: DISEASE-MODEL INPUT BUILDER
# scripts/03_build_disease_model.R
#
# Disease structure: World Heart Federation RHD stages A/B/C/D.
#   No RHD -> A <-> B <-> C -> D -> RHD death, with competing other-cause death
#   from every living stage. Incident RHD enters stage A (no separate ARF state).
# Focus: scale-up of SECONDARY PREVENTION =
#   (a) echocardiographic screening, (b) diagnosis, and (c) secondary antibiotic
#   prophylaxis / optimal treatment (SAP) which reduces RHD-specific MORTALITY.
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT
# ------------------------------------------------------------------------------
# This script ONLY BUILDS AND VALIDATES the disease-model INPUTS. It does NOT run
# the baseline projection or the intervention scenarios — that is done by the
# initial-state assembler (05_build_baseline.R) and the model runner
# (06_run_prevention_model.R). Economics live in 08 only.
#
# It assembles data-fed age x sex ( x year) rate arrays PLUS the tagged A/B/C/D
# natural-history, intervention-effect, surgery-service and care-cascade
# parameter blocks, validates them, and PERSISTS a single self-contained input
# bundle to disk so 04/05/06 consume it without recomputation.
#
#   INPUTS
#     data-raw/temp_baseline_rates_gbd.rds   (from 01_prepare_inputs.R)
#         GBD 2023: RHD + All-causes; Deaths / Prevalence / Incidence; Number+Rate;
#         1990-2023; 22 GBD age groups. Used (metric = Rate, per capita):
#           RHD Incidence  -> incident inflow of new RHD (enters stage A)
#           RHD Prevalence -> prevalent RHD seed at the first model year
#           RHD Deaths     -> RHD cause-specific mortality anchor
#           All-cause Deaths - RHD Deaths -> background (other-cause) mortality
#     data/pop_projection_2025_2100.rds      (from 02_build_demography.R)
#         single-year population by age x sex x year (persons), Indonesia. Read here
#         ONLY to establish the projection horizon (year grid) and a base-year
#         order-of-magnitude anchor. 05 re-reads the full population itself.
#
#   OUTPUT (written to wd_data = data/):
#     disease_model_inputs.rds  — a named list (see section 8) holding:
#       rates_by_year : ir/mort_rhd/mort_all/oth_mort as [age x sex x year] arrays,
#                       plus prev_seed [age x sex] (single-year prevalent seed);
#       transitions   : A/B/C/D annual transition probabilities (tagged);
#       p_rhd_death   : per-stage RHD-specific death probabilities c(A,B,C,D);
#       effects       : SAP RRR on RHD mortality + surgery RRRs (tagged);
#       surgery       : surgery-need fractions + coverage (ref/scale-up) [service];
#       stage_split   : prevalent-pool split across A/B/C/D (Cannon + rhd_d_fraction);
#       coverage      : 3-metric care-cascade STRUCTURE (baselines, targets, ramp
#                       window) — realised into trajectories by 05;
#       meta          : AGES, SEXES, years, LOCATION, RATE_BASE_YEAR, incidence_trend.
#
# GBD gives TOTAL RHD prevalence (not split by stage) and rates only to 2023:
#   * the prevalent pool is split into A/B/C/D by a tagged vector (Cannon et al
#     multi-state severity split for A/B/C, normalised to the non-D share; D from
#     the Indonesia RHD-with-heart-failure fraction, rhd_d_fraction [CALIBRATE]);
#   * for the projection horizon the GBD age-sex RATE pattern is held at its 2023
#     level, with an explicit incidence secular trend ([CALIBRATE]) baked into the
#     incidence array (anchored at the first projection year). Both are flagged.
#
# PARAMETER TAGS: [PAPER] article main text; [LIT] literature value (source in
#   comment); [CALIBRATE] tune to setting/appendix. Edit the parameters in
#   00_run_all.R; NO monetary values live here (unit costs / VSL / DALY in 08 only).
#
# SURGERY IS NOT A HEALTH STATE. There is no surgery / post-surgery stock. Surgery
#   is a clinical SERVICE (a fraction of the C/D stock treated each cycle) whose
#   only epidemiological effect is a risk reduction on C->D and D->RHD-death.
# ==============================================================================

library(data.table)

# ------------------------------------------------------------------------------
# 0. PATHS + DIMENSIONS
# ------------------------------------------------------------------------------
if (!exists("wd_raw"))  wd_raw  <- paste0(here::here("data-raw"), "/")
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

# getp(): use the 00_run_all.R global when set, else the documented standalone
# default below. Every user-controllable parameter is funnelled through this so
# 00_run_all.R can act as the single control panel while 03 still runs standalone.
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

LOCATION <- getp("LOCATION", "Indonesia")   # engine runs per location
AGES     <- 0:100                # single-year ages (0..100, 100 = 100+ open group)
SEXES    <- c("Female", "Male")
n_age    <- length(AGES)

RATE_BASE_YEAR <- getp("RATE_BASE_YEAR", 2023L)  # latest GBD year
POP_PROJ_FILE  <- paste0(wd_data, "pop_projection_2025_2100.rds")

# ------------------------------------------------------------------------------
# 1. RAMP WINDOW (structure only; no discounting/economics here)
# ------------------------------------------------------------------------------
ramp_start <- as.integer(getp("ramp_start", 2026L))   # first year cascade begins rising
ramp_end   <- as.integer(getp("ramp_end",   2050L))   # year target cascade is reached
incidence_trend <- getp("incidence_trend", 0.985)     # [CALIBRATE] ~15%/decade decline

# ------------------------------------------------------------------------------
# 2. A/B/C/D NATURAL HISTORY  (annual transition probabilities; tagged)
#    Competing OTHER-cause mortality is data-fed (oth_mort, below) and NOT set
#    here. Per-stage RHD-specific death probabilities are the mortality lever
#    that SAP (all stages) and surgery (stage D) act on.
# ------------------------------------------------------------------------------
transitions <- list(
  p_A_to_no_rhd = getp("p_A_to_no_rhd", 0.005),  # [CALIBRATE] A regression to No RHD
  p_A_to_B      = getp("p_A_to_B",      0.020),  # [CALIBRATE] A -> B
  p_B_to_A      = getp("p_B_to_A",      0.010),  # [CALIBRATE] B -> A regression
  p_B_to_C      = getp("p_B_to_C",      0.030),  # [CALIBRATE] B -> C
  p_C_to_B      = getp("p_C_to_B",      0.005),  # [CALIBRATE] C -> B regression
  p_C_to_D      = getp("p_C_to_D",      0.060),  # [CALIBRATE] C -> D (surgery lowers this)
  p_D_to_C      = getp("p_D_to_C",      0.000)   # [CALIBRATE] D -> C regression (default 0)
)

p_rhd_death <- c(
  A = getp("p_rhd_death_A", 0.0005),  # [CALIBRATE] untreated RHD death /yr, stage A
  B = getp("p_rhd_death_B", 0.0020),  # [CALIBRATE]                         stage B
  C = getp("p_rhd_death_C", 0.0200),  # [CALIBRATE]                         stage C
  D = getp("p_rhd_death_D", 0.0800)   # [CALIBRATE]                         stage D
)

# ------------------------------------------------------------------------------
# 3. INTERVENTION EFFECTS  (relative risk reductions on [0,1])
# ------------------------------------------------------------------------------
effects <- list(
  sap_rrr_rhd_death          = getp("sap_rrr_rhd_death", 0.55),  # [PAPER] SAP: 55% RRR on RHD mortality (all stages)
  eff_surgery_C_to_D         = getp("eff_surgery_C_to_D", 0.85), # [LIT] surgery RRR on C -> D progression
  eff_surgery_D_to_rhd_death = getp("eff_surgery_D_to_rhd_death", 0.85) # [LIT] surgery RRR on D -> RHD death
)

# ------------------------------------------------------------------------------
# 4. SURGERY SERVICE  (fractions requiring surgery + coverage; NOT a state)
#    Held equal in both arms by default (tertiary care fixed) so surgery is a
#    background cost, not a driver of incremental program cost.
# ------------------------------------------------------------------------------
surgery <- list(
  frac_C_requiring_surgery = getp("frac_C_requiring_surgery", 0.03), # [CALIBRATE]
  frac_D_requiring_surgery = getp("frac_D_requiring_surgery", 0.20), # [CALIBRATE]
  coverage_ref_baseline = getp("cov_surgery_ref",      0.05),        # reference arm (held)
  coverage_ref_target   = getp("cov_surgery_ref",      0.05),
  coverage_up_baseline  = getp("cov_surgery_scale_up", 0.05),        # scale-up arm (held)
  coverage_up_target    = getp("cov_surgery_scale_up", 0.05)
)

# ------------------------------------------------------------------------------
# 5. STAGE SPLIT of the prevalent RHD pool at the seed year  (A/B/C/D)
#    A/B/C from Cannon et al (multi-state model, n=591): 56.5 / 27.2 / 16.2%,
#    normalised and applied to the non-D share; D = rhd_d_fraction [CALIBRATE]
#    (Indonesia RHD-with-heart-failure / complicated-RHD fraction; GBD-informed).
# ------------------------------------------------------------------------------
rhd_d_fraction <- getp("rhd_d_fraction", 0.10)   # [CALIBRATE]
if (rhd_d_fraction < 0 || rhd_d_fraction >= 1)
  stop("rhd_d_fraction must be in [0, 1).", call. = FALSE)
abc_raw    <- c(A = 0.565, B = 0.272, C = 0.162)  # [LIT] Cannon et al
abc        <- abc_raw / sum(abc_raw)              # normalise (they sum to ~0.999)
stage_split <- c(abc * (1 - rhd_d_fraction), D = rhd_d_fraction)
stopifnot(abs(sum(stage_split) - 1) < 1e-9)

# ------------------------------------------------------------------------------
# 6. CARE-CASCADE STRUCTURE  (3 metrics; realised into trajectories by 05)
#    Cumulative population coverages. Reference holds baselines flat over the
#    horizon; scale-up ramps each linearly from the 2026 baseline to the 2050
#    target (held to 2100). Effective diagnosis/treatment are capped by earlier
#    cascade stages in 05 (effective_treatment <= diagnosis <= screening).
# ------------------------------------------------------------------------------
coverage <- list(
  # screening
  screen_ref_baseline = getp("screen_base",   0.05), screen_ref_target = getp("screen_base",   0.05),
  screen_up_baseline  = getp("screen_base",   0.05), screen_up_target  = getp("screen_target", 0.80),
  # diagnosis
  diagnosis_ref_baseline = getp("diagnosis_base",   0.05), diagnosis_ref_target = getp("diagnosis_base",   0.05),
  diagnosis_up_baseline  = getp("diagnosis_base",   0.05), diagnosis_up_target  = getp("diagnosis_target", 0.80),
  # optimal treatment (SAP)
  treatment_ref_baseline = getp("treatment_base",   0.04), treatment_ref_target = getp("treatment_base",   0.04),
  treatment_up_baseline  = getp("treatment_base",   0.04), treatment_up_target  = getp("treatment_target", 0.65),
  ramp_start = ramp_start,
  ramp_end   = ramp_end,
  # optional school-age-only screening mask (OFF by default; sensitivity only)
  screen_age_restrict = isTRUE(getp("screen_age_restrict", FALSE)),
  screen_age_lo = as.integer(getp("screen_age_lo",  5L)),
  screen_age_hi = as.integer(getp("screen_age_hi", 15L))
)
stopifnot(coverage$ramp_end > coverage$ramp_start)

# ==============================================================================
# 7. BUILD DATA-FED RATE INPUTS  (single age x sex arrays for LOCATION)
# ==============================================================================

## 7a. GBD 2023 base-year rates (per capita) at single-year age -----------------
gbd <- as.data.table(readRDS(paste0(wd_raw, "temp_baseline_rates_gbd.rds")))
gbd <- gbd[location_name == LOCATION & metric_name == "Rate" & year == RATE_BASE_YEAR,
           .(sex = sex_name, age_group = age_name,
             cause = cause_name, measure = measure_name, rate = val / 1e5)]

# single-year age -> GBD age-group label (matches 01's 22 groups incl. <5 split).
# This SAME mapping is mirrored in 04_calibration_random_tp.R.
age_to_gbd_group <- function(a) {
  fcase(
    a < 1,  "<1 year",
    a < 2,  "12-23 months",
    a < 5,  "2-4 years",
    a < 10, "5-9 years",
    a < 15, "10-14 years", a < 20, "15-19 years", a < 25, "20-24 years",
    a < 30, "25-29 years", a < 35, "30-34 years", a < 40, "35-39 years",
    a < 45, "40-44 years", a < 50, "45-49 years", a < 55, "50-54 years",
    a < 60, "55-59 years", a < 65, "60-64 years", a < 70, "65-69 years",
    a < 75, "70-74 years", a < 80, "75-79 years", a < 85, "80-84 years",
    a < 90, "85-89 years", a < 95, "90-94 years",
    default = "95+ years"
  )
}
age_map <- data.table(age = AGES, age_group = age_to_gbd_group(AGES))

# helper: pull one (cause, measure) rate -> [age x sex] matrix, missing -> 0
rate_matrix <- function(cause_name, measure_name) {
  d <- gbd[cause == cause_name & measure == measure_name,
           .(sex, age_group, rate)]
  d <- age_map[d, on = "age_group", allow.cartesian = TRUE]   # expand to single age
  m <- matrix(0, n_age, length(SEXES), dimnames = list(AGES, SEXES))
  for (s in SEXES) {
    ds <- d[sex == s]
    if (nrow(ds)) m[as.character(ds$age), s] <- ds$rate
  }
  m
}

ir_rhd0   <- rate_matrix("Rheumatic heart disease", "Incidence")   # incidence rate
prev_rhd0 <- rate_matrix("Rheumatic heart disease", "Prevalence")  # prevalent fraction
mort_rhd0 <- rate_matrix("Rheumatic heart disease", "Deaths")      # RHD death rate
mort_all0 <- rate_matrix("All causes",              "Deaths")      # all-cause death rate
oth_mort0 <- pmax(mort_all0 - mort_rhd0, 0)   # background (non-RHD) competing mortality

## 7b. establish the projection horizon (year grid) from 02's projection table --
if (!file.exists(POP_PROJ_FILE))
  stop("Missing population projection input:\n  ", POP_PROJ_FILE,
       "\n  Run 02_build_demography.R first.", call. = FALSE)
pop_proj <- as.data.table(readRDS(POP_PROJ_FILE))
pop_proj <- pop_proj[location == LOCATION & age %in% AGES & sex %in% SEXES]
years    <- sort(unique(pop_proj$year))     # 2025..2100  (HORIZON driven by input)
n_years  <- length(years)

## 7c. broadcast base-year [age x sex] matrices to [age x sex x year] arrays -----
##     incidence carries the secular trend (anchored at the first projection year);
##     the mortality envelope is held at the base-year pattern across the horizon.
broadcast_years <- function(mat, mult = rep(1, n_years)) {
  arr <- array(0, dim = c(n_age, length(SEXES), n_years),
               dimnames = list(AGES, SEXES, years))
  for (iy in seq_len(n_years)) arr[, , iy] <- mat * mult[iy]
  arr
}
trend_mult <- incidence_trend^(seq_len(n_years) - 1L)   # 1.0 at first projection year

ir_rhd_yr   <- broadcast_years(ir_rhd0, trend_mult)   # incidence, trend baked in
mort_rhd_yr <- broadcast_years(mort_rhd0)             # RHD death rate (held)
mort_all_yr <- broadcast_years(mort_all0)             # all-cause death rate (held)
oth_mort_yr <- broadcast_years(oth_mort0)             # background mortality (held)

# base-year population [age x sex] for the order-of-magnitude anchor below
py1  <- pop_proj[year == years[1]]
pop1 <- matrix(0, n_age, length(SEXES), dimnames = list(AGES, SEXES))
for (s in SEXES) {
  ps <- py1[sex == s]
  pop1[as.character(ps$age), s] <- ps$Nx
}

# ------------------------------------------------------------------------------
# 8. INPUT VALIDATION  (fail loudly BEFORE persisting)
# ------------------------------------------------------------------------------
chk_rate <- function(x, nm, prob = TRUE) {
  if (any(is.na(x)))  stop("Input '", nm, "' contains NA.", call. = FALSE)
  if (any(x < 0))     stop("Input '", nm, "' contains negative values.", call. = FALSE)
  if (prob && any(x > 1))
    stop("Input '", nm, "' has values > 1 (should be a per-capita probability/rate).",
         call. = FALSE)
}
chk_rate(ir_rhd_yr,   "RHD incidence rate")
chk_rate(prev_rhd0,   "RHD prevalence rate")
chk_rate(mort_rhd_yr, "RHD death rate")
chk_rate(mort_all_yr, "all-cause death rate")
chk_rate(oth_mort_yr, "background (other-cause) mortality")

# arrays complete & correctly shaped
for (nm in c("ir_rhd_yr", "mort_rhd_yr", "mort_all_yr", "oth_mort_yr")) {
  a <- get(nm)
  if (!identical(dim(a), c(n_age, length(SEXES), n_years)))
    stop("Array '", nm, "' has wrong dimensions.", call. = FALSE)
}
if (any(is.na(pop1)) || any(pop1 < 0) || sum(pop1) <= 0)
  stop("Base-year population slice is empty/invalid.", call. = FALSE)
if (!all(SEXES %in% pop_proj$sex))
  stop("Population projection is missing a sex.", call. = FALSE)

# GBD rates must actually carry RHD signal (not all zero)
if (sum(prev_rhd0) == 0 || sum(ir_rhd0) == 0)
  stop("GBD RHD prevalence/incidence rates are all zero -- check inputs.", call. = FALSE)

# competing-risk headroom at the base year: IR + background <= 1
if (any(ir_rhd0 + oth_mort0 > 1 + 1e-9))
  stop("IR + background mortality exceeds 1 at some age-sex cell.", call. = FALSE)

# --- A/B/C/D parameter-block sanity -------------------------------------------
# every transition / death / effect probability in [0,1]
# cascade coverage keys are exactly the *_baseline / *_target entries (this
# excludes screen_age_lo/hi/restrict, which are not probabilities).
cascade_keys <- grep("_(baseline|target)$", names(coverage), value = TRUE)
prob_params <- c(unlist(transitions), p_rhd_death, unlist(effects),
                 surgery$frac_C_requiring_surgery, surgery$frac_D_requiring_surgery,
                 surgery$coverage_ref_baseline, surgery$coverage_ref_target,
                 surgery$coverage_up_baseline,  surgery$coverage_up_target,
                 unlist(coverage[cascade_keys]))
if (any(prob_params < 0 | prob_params > 1))
  stop("A transition/effect/surgery/coverage parameter is outside [0,1].", call. = FALSE)

# baseline OUTGOING disease-side probabilities from each stage must sum to < 1
# (competing other-cause mortality is data-fed age x sex and enforced in the
# engine via a non-negative 'stay' floor; excluded from this scalar check).
out_A <- transitions$p_A_to_no_rhd + transitions$p_A_to_B + p_rhd_death[["A"]]
out_B <- transitions$p_B_to_A      + transitions$p_B_to_C + p_rhd_death[["B"]]
out_C <- transitions$p_C_to_B      + transitions$p_C_to_D + p_rhd_death[["C"]]
out_D <- transitions$p_D_to_C      + p_rhd_death[["D"]]
bad_stage <- c(A = out_A, B = out_B, C = out_C, D = out_D)
if (any(bad_stage >= 1))
  stop("Baseline outgoing probabilities sum to >= 1 for stage(s): ",
       paste(names(bad_stage)[bad_stage >= 1], collapse = ", "), call. = FALSE)

# stage split is a proper distribution
if (abs(sum(stage_split) - 1) > 1e-9 || any(stage_split < 0))
  stop("stage_split is not a valid distribution over A/B/C/D.", call. = FALSE)

# screening age window valid (used only if screen_age_restrict = TRUE)
if (coverage$screen_age_lo < min(AGES) || coverage$screen_age_hi > max(AGES) ||
    coverage$screen_age_lo > coverage$screen_age_hi)
  stop("Screening age window is invalid.", call. = FALSE)

# order-of-magnitude anchors vs GBD (base-year rate x base-year population)
rhd_prev_count  <- sum(prev_rhd0 * pop1)
rhd_death_count <- sum(mort_rhd0 * pop1)
all_death_count <- sum(mort_all0 * pop1)
message(sprintf(
  "Data-fed inputs OK | horizon %d-%d | ages %d-%d\n  base-year(%d): RHD prevalence ~ %s | RHD deaths ~ %s | all-cause deaths ~ %s",
  as.integer(min(years)), as.integer(max(years)),
  as.integer(min(AGES)), as.integer(max(AGES)), as.integer(years[1]),
  formatC(round(rhd_prev_count),  format = "d", big.mark = ","),
  formatC(round(rhd_death_count), format = "d", big.mark = ","),
  formatC(round(all_death_count), format = "d", big.mark = ",")))
# sane bands for Indonesia RHD (order of magnitude): prevalence 1e5-1e7, deaths 1e3-1e5
if (rhd_prev_count < 1e5 || rhd_prev_count > 1e7)
  stop("Base-year RHD prevalence count ", round(rhd_prev_count),
       " is outside the sane band 1e5-1e7.", call. = FALSE)
if (rhd_death_count < 1e3 || rhd_death_count > 1e5)
  stop("Base-year RHD death count ", round(rhd_death_count),
       " is outside the sane band 1e3-1e5.", call. = FALSE)

# ==============================================================================
# 9. ASSEMBLE + PERSIST THE DISEASE-MODEL INPUT BUNDLE
# ==============================================================================
disease_model_inputs <- list(
  rates_by_year = list(
    ir_rhd   = ir_rhd_yr,     # [age x sex x year] incidence (secular trend baked in)
    mort_rhd = mort_rhd_yr,   # [age x sex x year] RHD cause-specific death rate (held)
    mort_all = mort_all_yr,   # [age x sex x year] all-cause death rate (held)
    oth_mort = oth_mort_yr,   # [age x sex x year] background (non-RHD) mortality (held)
    prev_seed = prev_rhd0     # [age x sex] prevalent RHD fraction — single-year seed
  ),
  transitions = transitions,  # A/B/C/D annual transition probabilities [CALIBRATE]
  p_rhd_death = p_rhd_death,   # per-stage RHD-specific death probability c(A,B,C,D)
  effects     = effects,       # SAP RRR on RHD mortality + surgery RRRs [PAPER]/[LIT]
  surgery     = surgery,       # surgery-service fractions + coverage (ref/scale-up)
  stage_split = stage_split,   # prevalent-pool split over A/B/C/D (Cannon + rhd_d_fraction)
  coverage    = coverage,      # 3-metric care-cascade STRUCTURE (baselines/targets/window)
  meta = list(
    LOCATION       = LOCATION,
    AGES           = AGES,
    SEXES          = SEXES,
    years          = years,
    RATE_BASE_YEAR = as.integer(RATE_BASE_YEAR),
    incidence_trend = incidence_trend,
    rhd_d_fraction  = rhd_d_fraction,
    built_from = c("data-raw/temp_baseline_rates_gbd.rds",
                   basename(POP_PROJ_FILE))
  )
)

OUT_FILE <- paste0(wd_data, "disease_model_inputs.rds")
saveRDS(disease_model_inputs, file = OUT_FILE)

message(sprintf(
  "03_build_disease_model.R complete | wrote %s\n  arrays [age=%d x sex=%d x year=%d]; A/B/C/D stage_split = %s",
  basename(OUT_FILE), n_age, length(SEXES), n_years,
  paste(sprintf("%s %.1f%%", names(stage_split), 100 * stage_split), collapse = " | ")))
message(sprintf("  stage outgoing-prob sums (excl. background): A %.3f | B %.3f | C %.3f | D %.3f",
                out_A, out_B, out_C, out_D))
message("  Next: 04_calibration_random_tp.R (calibrate IR/CF), then 05_build_baseline.R.")
