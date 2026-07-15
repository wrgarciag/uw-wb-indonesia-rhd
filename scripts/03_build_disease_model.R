# ==============================================================================
# RHD secondary-prevention investment case: cohort state-transition model
# scripts/03_build_disease_model.R
#
# Structure after: Coates et al., Lancet Glob Health 2021 (PMC9087136).
# Focus: scale-up of SECONDARY PREVENTION =
#   (a) echocardiographic screening to detect asymptomatic (mild) RHD, and
#   (b) secondary antibiotic prophylaxis (SAP) for screen-detected mild RHD.
#
# ------------------------------------------------------------------------------
# WHAT CHANGED IN THIS REFACTOR (data-fed, age-sex structured)
# ------------------------------------------------------------------------------
# The former hard-coded Australia-flavoured scalars (pop_2021 <- 288e6, au_pop,
# gdp_pc_2019 <- 1800, prev_rhd_2021 <- 1.4e6, rhd_incident_2021 <- 66000,
# seed_frac, ...) are removed. Demographic and epidemiological MAGNITUDES are now
# joined in from upstream, by single-year age x sex x year x location:
#
#   INPUTS
#     data-raw/temp_baseline_rates_gbd.rds   (from 01_prepare_inputs.R)
#         GBD 2023: RHD + All-causes; Deaths / Prevalence / Incidence; Number+Rate;
#         1990-2023; 22 GBD age groups. Used (metric = Rate, per 100 000):
#           RHD Incidence  -> incident inflow of new asymptomatic RHD
#           RHD Prevalence -> prevalent RHD seed at the first model year
#           RHD Deaths     -> validation anchor for model-produced RHD deaths
#           All-cause & RHD Deaths -> background (other-cause) mortality = all - RHD
#     data/pop_projection_2026_2100.rds      (from 02_build_demography.R)
#         single-year population by age x sex x year (persons), Indonesia.
#     data/pop_observed_1990_2023.rds        (from 02; used only for a 2023 sanity print)
#
# The seed prevalent pool, the incident inflow, the age-sex structure, and the
# competing (other-cause) mortality are ALL data-fed. Only the clinical
# progression/case-fatality parameters, the intervention effect sizes, the
# coverage ramps and the unit costs remain as explicit tagged parameters
# ([PAPER]/[LIT]/[CALIBRATE]) -- exactly the levers the interventions act on.
#
# GBD gives TOTAL RHD prevalence (not split by severity) and rates only to 2023:
#   * the prevalent pool is split into mild/severe/post by a tagged [LIT] vector
#     (asymptomatic RHD dominates prevalence);
#   * for 2026-2100 the GBD age-sex RATE pattern is held at its 2023 level, with
#     an explicit incidence secular trend ([CALIBRATE]).  Both are flagged.
#
#   HORIZON: 2026-2100 (driven by the demography projection table). Reference vs
#     SAP scale-up; SAP ramps over ramp_start..ramp_end. Cohorts age one year per
#     cycle; new incident cases enter at every age each year.
#
#   OUTPUTS (in-memory, consumed by 07_make_outputs.R -- CONTRACT PRESERVED):
#     ref, sap  : data.table national totals by year with columns
#                 year, pop, new_rhd, mild, severe, post, rhd_deaths, surgeries,
#                 n_on_sap, n_screened, c_screen, c_sap, c_hf, c_surg, cost
#     plus scalar globals years, disc_rate, gdp_pc_2019, gdp_growth, vsl_mult,
#     dalys_per_death (07 reads these).
#
# PARAMETER TAGS: [PAPER] article main text; [LIT] literature value (source in
#   comment); [CALIBRATE] tune to setting/appendix. Edit only PARAMETERS/COVERAGE/
#   COST blocks; the demographic & epidemiological inputs come from the tables.
# ==============================================================================

library(data.table)
set.seed(1)

# ------------------------------------------------------------------------------
# 0. PATHS + DIMENSIONS
# ------------------------------------------------------------------------------
if (!exists("wd_raw"))  wd_raw  <- paste0(here::here("data-raw"), "/")
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")

LOCATION <- "Indonesia"          # parameterised: the engine runs per location
AGES     <- 0:100                # single-year ages (0..100, 100 = 100+ open group)
SEXES    <- c("Female", "Male")
n_age    <- length(AGES)

# ------------------------------------------------------------------------------
# 1. TIME, RAMP WINDOW, DISCOUNTING  (horizon driven by the demography input)
# ------------------------------------------------------------------------------
disc_rate <- 0.03                # [PAPER] 3% discount on costs and benefits

# --- flexible scale-up window (change these two lines only) --------------------
ramp_start <- 2026               # first year coverage begins rising
ramp_end   <- 2030               # year target coverage is reached
# ------------------------------------------------------------------------------

RATE_BASE_YEAR <- 2023           # latest GBD year; forward horizon holds this pattern

# ------------------------------------------------------------------------------
# 2. DISEASE / CLINICAL PARAMETERS  (annual transition probabilities; tagged)
#    RHD-specific mortality & progression stay parameters because the SAP / HF /
#    surgery levers act on them; competing (other-cause) mortality is data-fed.
# ------------------------------------------------------------------------------
p <- list(
  incidence_trend      = 0.985,  # [CALIBRATE] ~15%/decade incidence decline (Table 1 footnote)

  p_mild_to_severe     = 0.010,  # [CALIBRATE] asymptomatic -> heart failure /yr
                                 #   (tuned so uncalibrated baseline RHD deaths sit
                                 #    near the GBD 2023 rate-implied level; 04 refines)
  p_severe_death       = 0.09,   # [LIT] REMEDY ~17%/2yr -> ~9%/yr untreated severe
  p_surg_op_mortality  = 0.03,   # [PAPER] 3% operative mortality (Table 1 footnote)
  p_post_death_rhd     = 0.020,  # [LIT] residual RHD mortality after valve surgery

  eff_sap_asymp        = 0.55,   # [PAPER] Table 1 #2b: SAP in asymptomatic RHD, 55% (7-78)
  eff_hf_mgmt          = 0.60,   # [PAPER] Table 1 #3: HF management, 60% (30-80)
  eff_surgery          = 0.85,   # [PAPER] Table 1 #4: surgery, 85% (70-92)

  frac_severe_surg_elig = 0.50   # [CALIBRATE] share of severe RHD surgery-eligible
)

# prevalent RHD severity split (GBD gives only TOTAL prevalence) -- [LIT]/[CALIBRATE]
# Echo-detected RHD is overwhelmingly subclinical/asymptomatic; symptomatic HF
# (severe) and post-surgical stocks are small shares. Tuned so the uncalibrated
# baseline RHD deaths sit near the GBD 2023 rate-implied level (04 refines).
seed_split <- c(mild = 0.96, severe = 0.03, post = 0.01)   # [LIT]/[CALIBRATE] must sum to 1
stopifnot(abs(sum(seed_split) - 1) < 1e-9)

# ------------------------------------------------------------------------------
# 3. ECONOMIC PARAMETERS  (Indonesia; consumed by 07_make_outputs.R)
# ------------------------------------------------------------------------------
gdp_pc_2019     <- 4150          # [LIT] Indonesia GDP per capita, 2019 US$ (~US$4,150)
gdp_growth      <- 0.03          # [CALIBRATE] real per-capita growth for VSL/GDP benefit
vsl_mult        <- 30            # [PAPER] "VSL is roughly 30 times GDP per capita"
dalys_per_death <- 30            # [CALIBRATE] undiscounted DALYs per RHD death

# ------------------------------------------------------------------------------
# 4. UNIT COSTS (2019 US$)  --  [LIT]/[CALIBRATE]; replace with appendix values
# ------------------------------------------------------------------------------
costs <- list(
  cost_screen_per_person = 12,   # [LIT] echo screening cost/person (handheld/decentralised)
  cost_sap_per_year      = 45,   # [LIT] annual SAP cost/person (monthly BPG + visits)
  cost_hf_per_year       = 120,  # [LIT] annual HF medical management cost/person
  cost_surgery           = 9000, # [LIT] valve surgery + first-year post-op cost/person
  screen_age_lo          = 5,    # [LIT] school-based echo screening age window
  screen_age_hi          = 15
)

# ==============================================================================
# 5. BUILD DATA-FED INPUTS  (single age x sex x year arrays for LOCATION)
# ==============================================================================

## 5a. GBD 2023 base-year rates (per capita) at single-year age -----------------
gbd <- as.data.table(readRDS(paste0(wd_raw, "temp_baseline_rates_gbd.rds")))
gbd <- gbd[location_name == LOCATION & metric_name == "Rate" & year == RATE_BASE_YEAR,
           .(sex = sex_name, age_group = age_name,
             cause = cause_name, measure = measure_name, rate = val / 1e5)]

# single-year age -> GBD age-group label (matches 01's 22 groups incl. <5 split)
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

## 5b. population projection 2026-2100 (single age x sex) ----------------------
pop_proj <- as.data.table(readRDS(paste0(wd_data, "pop_projection_2026_2100.rds")))
pop_proj <- pop_proj[location == LOCATION & age %in% AGES & sex %in% SEXES]

years   <- sort(unique(pop_proj$year))     # 2026..2100  (HORIZON driven by input)
n_years <- length(years)

# population array [age x sex x year]
pop_arr <- array(0, dim = c(n_age, length(SEXES), n_years),
                 dimnames = list(AGES, SEXES, years))
for (iy in seq_along(years)) {
  py <- pop_proj[year == years[iy]]
  m  <- matrix(0, n_age, length(SEXES), dimnames = list(AGES, SEXES))
  for (s in SEXES) {
    ps <- py[sex == s]
    m[as.character(ps$age), s] <- ps$Nx
  }
  pop_arr[, , iy] <- m
}

# ------------------------------------------------------------------------------
# 6. INPUT VALIDATION  (fail loudly BEFORE running the engine)
# ------------------------------------------------------------------------------
chk <- function(x, nm) {
  if (any(is.na(x)))  stop("Input '", nm, "' contains NA.", call. = FALSE)
  if (any(x < 0))     stop("Input '", nm, "' contains negative values.", call. = FALSE)
}
chk(ir_rhd0, "RHD incidence rate");  chk(prev_rhd0, "RHD prevalence rate")
chk(mort_rhd0, "RHD death rate");    chk(mort_all0, "all-cause death rate")
chk(oth_mort0, "other-cause mortality"); chk(pop_arr, "population projection")

# every (age, sex, year) population cell present and positive somewhere
if (any(colSums(apply(pop_arr, 3, colSums)) <= 0))
  stop("Population projection has an empty sex-year slice.", call. = FALSE)
if (!all(SEXES %in% pop_proj$sex))
  stop("Population projection is missing a sex.", call. = FALSE)
# GBD rates must actually carry RHD signal (not all zero)
if (sum(prev_rhd0) == 0 || sum(ir_rhd0) == 0)
  stop("GBD RHD prevalence/incidence rates are all zero -- check inputs.", call. = FALSE)

message(sprintf("Data-fed inputs OK | horizon %d-%d | ages %d-%d | RHD prev(all-age, %d) ~ %s",
                as.integer(min(years)), as.integer(max(years)),
                as.integer(min(AGES)), as.integer(max(AGES)), as.integer(RATE_BASE_YEAR),
                formatC(round(sum(prev_rhd0 * pop_arr[, , 1])), format = "d", big.mark = ",")))

# ------------------------------------------------------------------------------
# 7. COVERAGE TRAJECTORIES  (flexible ramp window; unchanged logic)
# ------------------------------------------------------------------------------
ramp <- function(baseline, target, start = ramp_start, end = ramp_end) {
  frac <- (years - start) / (end - start)
  frac <- pmin(pmax(frac, 0), 1)
  baseline + (target - baseline) * frac
}
cov <- list(
  sap_asymp_ref = rep(0.05, n_years),      # Table 1 #2b: 5.0% baseline
  sap_asymp_up  = ramp(0.05, 0.40),        #              -> 40% scaled up
  hf_ref   = rep(0.08, n_years), hf_up   = rep(0.08, n_years),   # Table 1 #3
  surg_ref = rep(0.05, n_years), surg_up = rep(0.05, n_years)    # Table 1 #4
)

# ==============================================================================
# 8. MODEL ENGINE -- one scenario, age-sex-structured cohorts, returns by year
#    Same well-known mild -> severe -> post logic as before, now applied per
#    [age x sex] cell with data-fed incidence, seed and competing mortality, and
#    with cohorts ageing one year per cycle.
# ==============================================================================
zero_mat  <- function() matrix(0, n_age, length(SEXES), dimnames = list(AGES, SEXES))
age_shift <- function(M) {                    # advance age a -> a+1; 100 is open
  N <- zero_mat()
  N[2:n_age, ] <- M[1:(n_age - 1), ]
  N[n_age,  ]  <- N[n_age, ] + M[n_age, ]     # accumulate terminal 100+ group
  N
}

run_scenario <- function(cov_sap_asymp, cov_hf, cov_surg) {
  mild <- severe <- post <- zero_mat()
  out  <- vector("list", n_years)

  for (i in seq_len(n_years)) {
    popm <- pop_arr[, , i]
    othm <- oth_mort0                       # background (non-RHD) mortality, per age-sex
    tr   <- p$incidence_trend^(i - 1)       # secular incidence trend from base year

    # seed the prevalent pool at the first model year (data-fed from GBD prevalence)
    if (i == 1) {
      prevm  <- prev_rhd0 * popm
      mild   <- prevm * seed_split[["mild"]]
      severe <- prevm * seed_split[["severe"]]
      post   <- prevm * seed_split[["post"]]
    }

    # --- new incident asymptomatic (mild) RHD (data-fed incidence x population) --
    new_mild  <- ir_rhd0 * popm * tr
    mild_pool <- mild + new_mild
    prog_mild <- p$p_mild_to_severe * (1 - p$eff_sap_asymp * cov_sap_asymp[i])  # SAP slows progression
    mild_to_severe <- mild_pool * prog_mild
    mild_death_oth <- mild_pool * othm
    mild_next      <- pmax(mild_pool - mild_to_severe - mild_death_oth, 0)

    # --- severe RHD (heart failure) ---------------------------------------------
    severe_pool     <- severe + mild_to_severe
    surg_candidates <- severe_pool * p$frac_severe_surg_elig
    surgeries       <- surg_candidates * cov_surg[i]
    op_deaths       <- surgeries * p$p_surg_op_mortality
    to_post         <- surgeries - op_deaths
    remain_severe   <- severe_pool - surgeries
    sev_death_rhd   <- remain_severe * p$p_severe_death * (1 - p$eff_hf_mgmt * cov_hf[i])
    sev_death_oth   <- remain_severe * othm
    severe_next     <- pmax(remain_severe - sev_death_rhd - sev_death_oth, 0)

    # --- post-surgery -----------------------------------------------------------
    post_pool    <- post + to_post
    post_death_r <- post_pool * p$p_post_death_rhd
    post_death_o <- post_pool * othm
    post_next    <- pmax(post_pool - post_death_r - post_death_o, 0)

    rhd_deaths_m <- op_deaths + sev_death_rhd + post_death_r

    # --- costs this year (national) ---------------------------------------------
    school_pop <- sum(popm[as.character(costs$screen_age_lo:costs$screen_age_hi), ])
    n_screened <- school_pop * cov_sap_asymp[i]            # school-based echo screening
    n_on_sap   <- sum(mild_pool) * cov_sap_asymp[i]
    c_screen   <- n_screened * costs$cost_screen_per_person
    c_sap      <- n_on_sap   * costs$cost_sap_per_year
    c_hf       <- sum(remain_severe) * cov_hf[i] * costs$cost_hf_per_year
    c_surg     <- sum(surgeries) * costs$cost_surgery
    cost_total <- c_screen + c_sap + c_hf + c_surg

    out[[i]] <- data.table(
      year = years[i], pop = sum(popm),
      new_rhd = sum(new_mild), mild = sum(mild_next),
      severe = sum(severe_next), post = sum(post_next),
      rhd_deaths = sum(rhd_deaths_m), surgeries = sum(surgeries),
      n_on_sap = n_on_sap, n_screened = n_screened,
      c_screen = c_screen, c_sap = c_sap, c_hf = c_hf, c_surg = c_surg,
      cost = cost_total
    )

    # --- age the surviving stocks one year for the next cycle -------------------
    mild   <- age_shift(mild_next)
    severe <- age_shift(severe_next)
    post   <- age_shift(post_next)
  }
  rbindlist(out)
}

# ==============================================================================
# 9. RUN reference vs secondary-prevention scale-up
# ==============================================================================
ref <- run_scenario(cov$sap_asymp_ref, cov$hf_ref, cov$surg_ref)
sap <- run_scenario(cov$sap_asymp_up,  cov$hf_ref, cov$surg_up)

# ------------------------------------------------------------------------------
# 10. OUTPUT VALIDATION + a light GBD anchor check
# ------------------------------------------------------------------------------
for (nm in c("ref", "sap")) {
  d <- get(nm)
  if (nrow(d) != n_years) stop(nm, ": wrong number of years.", call. = FALSE)
  if (anyNA(d))           stop(nm, ": contains NA.", call. = FALSE)
  num <- d[, .SD, .SDcols = is.numeric]
  if (any(vapply(num, function(x) any(x < 0), logical(1))))
    stop(nm, ": contains negative values.", call. = FALSE)
}

# model RHD deaths in the first year vs GBD observed RHD deaths (order-of-magnitude)
gbd_rhd_deaths_base <- sum(mort_rhd0 * pop_arr[, , 1])
message(sprintf(
  "Model check | ref RHD deaths %d (%d) vs GBD-rate implied %d | deaths averted by 2100 (cum) = %s",
  as.integer(round(ref$rhd_deaths[1])), as.integer(years[1]),
  as.integer(round(gbd_rhd_deaths_base)),
  formatC(round(sum(ref$rhd_deaths - sap$rhd_deaths)), format = "d", big.mark = ",")))

message(sprintf("03_build_disease_model.R complete | scenarios ref & sap built for %d-%d.",
                as.integer(min(years)), as.integer(max(years))))
