# ==============================================================================
# RHD secondary-prevention investment case: cohort state-transition model (data.table)
#
# Structure after: Coates et al., Lancet Glob Health 2021 (PMC9087136).
#
# Focus: scale-up of SECONDARY PREVENTION only =
#   (a) RHD screening (echocardiography) to detect asymptomatic (mild) RHD, and
#   (b) secondary antibiotic prophylaxis (SAP) for screen-detected asymptomatic RHD.
#
# CHANGES IN THIS VERSION
#   1) ARF state and all ARF parameters removed. Incident RHD now enters the model
#      directly as new asymptomatic (mild) cases. SAP acts on mild -> severe only.
#   2) The model is SEEDED with an existing prevalent RHD population at the start
#      of 2021 (split across mild/severe/post-surgery); new incident cases are then
#      added each cycle on top of that stock.
#   3) The scale-up (ramp) window is now flexible: run 2021-2030 but ramp coverage
#      over any sub-period, e.g. 2026-2030. Coverage = baseline before the window,
#      ramps linearly across it, then holds at target afterwards.
#
# PARAMETER TAGS: [PAPER] from the article's main text (Table 1 etc.);
#   [LIT] plausible literature value (source in comment); [CALIBRATE] tune to
#   your setting / the appendix. Edit only the PARAMETERS/COVERAGE/COST blocks.
# ==============================================================================

library(data.table)
set.seed(1)

# ------------------------------------------------------------------------------
# 1. TIME, RAMP WINDOW, DISCOUNTING, DEMOGRAPHY
# ------------------------------------------------------------------------------
years      <- 2017:2050         # [PAPER] full model horizon
n_years    <- length(years)
disc_rate  <- 0.03               # [PAPER] 3% discount on costs and benefits

# --- flexible scale-up window (change these two lines only) -------------------
ramp_start <- 2026              # first year coverage begins rising (>= min(years))
ramp_end   <- 2030               # year target coverage is reached (<= max(years))
# e.g. to scale up only in the second half: ramp_start <- 2026; ramp_end <- 2030
# ------------------------------------------------------------------------------

pop_2021   <- 288e6             # [LIT] UN WPP 2019, AU ~1.36bn in 2021
pop_growth <- 0.025              # [LIT] ~2.5%/yr AU population growth
au_pop     <- pop_2021 * (1 + pop_growth)^(0:(n_years - 1))

gdp_pc_2019 <- 1800              # [LIT] AU GDP per capita, 2019 US$
gdp_growth  <- 0.03              # [CALIBRATE] real per-capita growth for VSL/GDP benefit

# ------------------------------------------------------------------------------
# 2. INITIAL PREVALENT RHD POPULATION (seed at start of 2021)
# ------------------------------------------------------------------------------
prev_rhd_2021 <- 1.4e6          # [CALIBRATE] total prevalent RHD in AU at 2021 start
seed_frac <- c(mild = 0.80, severe = 0.18, post = 0.02)  # [CALIBRATE] stock split
seed <- list(mild   = prev_rhd_2021 * seed_frac[["mild"]],
             severe = prev_rhd_2021 * seed_frac[["severe"]],
             post   = prev_rhd_2021 * seed_frac[["post"]])

# ------------------------------------------------------------------------------
# 3. DISEASE PARAMETERS (annual transition probabilities; ARF removed)
# ------------------------------------------------------------------------------
p <- list(
  # incidence inflow of NEW asymptomatic RHD (annual, AU-wide)
  rhd_incident_2021    = 66000,   # [CALIBRATE] new asymptomatic RHD/yr in 2021
  incidence_trend      = 0.985,    # [PAPER-ADJ] ~15% incidence decline/decade (Table 1 footnote)
  
  # progression / mortality (baseline = untreated)
  p_mild_to_severe     = 0.045,    # [CALIBRATE] asymptomatic -> heart failure /yr
  p_mild_death_other   = 0.006,    # [LIT] background mortality of mild-RHD population
  
  p_severe_death       = 0.02,    # [LIT] REMEDY ~17%/2yr -> ~9%/yr untreated severe
  p_severe_death_other = 0.010,    # [LIT] competing background mortality
  p_surg_op_mortality  = 0.03,     # [PAPER] 3% operative mortality (Table 1 footnote)
  
  p_post_death_rhd     = 0.020,    # [LIT] residual RHD mortality after valve surgery
  p_post_death_other   = 0.010,    # [LIT] background mortality post-surgery
  
  # intervention effect sizes (relative risk reductions)
  eff_sap_asymp        = 0.55,     # [PAPER] Table 1 #2b: SAP in asymptomatic RHD, 55% (7-78)
  eff_hf_mgmt          = 0.60,     # [PAPER] Table 1 #3: HF management, 60% (30-80)
  eff_surgery          = 0.85,     # [PAPER] Table 1 #4: surgery, 85% (70-92)
  
  frac_severe_surg_elig = 0.55     # [CALIBRATE] share of severe RHD surgery-eligible (aged 10-40)
)

# ------------------------------------------------------------------------------
# 4. COVERAGE TRAJECTORIES  (flexible ramp window)
# ------------------------------------------------------------------------------
# Coverage in year y:
#   baseline                        for y <= ramp_start
#   linear baseline -> target       for ramp_start < y < ramp_end
#   target                          for y >= ramp_end
ramp <- function(baseline, target, start = ramp_start, end = ramp_end) {
  frac <- (years - start) / (end - start)
  frac <- pmin(pmax(frac, 0), 1)
  baseline + (target - baseline) * frac
}

cov <- list(
  # Screening-enabled SAP for asymptomatic RHD -- Table 1 #2b: 5.0% -> 40%
  sap_asymp_ref = rep(0.05, n_years),
  sap_asymp_up  = ramp(0.05, 0.40),
  
  # Tertiary care held at baseline here (isolating secondary prevention).
  # Flip *_up to ramp() versions to add integrated care.
  hf_ref   = rep(0.08, n_years), hf_up   = rep(0.08, n_years),  # Table 1 #3: 8% -> (55%)
  surg_ref = rep(0.05, n_years), surg_up = rep(0.05, n_years)   # Table 1 #4: 5% -> (25%)
)

# ------------------------------------------------------------------------------
# 5. UNIT COSTS (2019 US$)  --  [CALIBRATE]/[LIT]; replace with appendix values
# ------------------------------------------------------------------------------
costs <- list(
  cost_screen_per_person = 12,     # [LIT] echo screening cost/person (decentralised/handheld)
  n_screened_per_year    = 8e6,    # [CALIBRATE] population echo-screened/yr at full coverage
  cost_sap_per_year      = 45,     # [LIT] annual SAP cost/person (monthly BPG + visits)
  cost_hf_per_year       = 120,    # [LIT] annual HF medical management cost/person
  cost_surgery           = 9000    # [LIT] valve surgery + first-year post-op cost/person
)

# ------------------------------------------------------------------------------
# 6. VALUE OF A STATISTICAL LIFE (full-income benefit)
# ------------------------------------------------------------------------------
vsl_mult        <- 30            # [PAPER] "VSL is roughly 30 times GDP per capita"
dalys_per_death <- 30            # [CALIBRATE] undiscounted DALYs per RHD death

# ==============================================================================
# 7. MODEL ENGINE -- one scenario, returns a data.table by year
# ==============================================================================
run_scenario <- function(cov_sap_asymp, cov_hf, cov_surg) {
  
  # seed the living stocks with the existing prevalent population
  st <- list(mild = seed$mild, severe = seed$severe, post = seed$post)
  
  out <- vector("list", n_years)
  
  for (i in seq_len(n_years)) {
    tr       <- p$incidence_trend^(i - 1)
    new_mild <- p$rhd_incident_2021 * tr           # NEW incident asymptomatic RHD added to stock
    
    # --- Mild (asymptomatic) RHD: carried stock + new incidence ---
    mild_pool <- st$mild + new_mild
    prog_mild <- p$p_mild_to_severe * (1 - p$eff_sap_asymp * cov_sap_asymp[i])  # SAP slows progression
    mild_to_severe <- mild_pool * prog_mild
    mild_death_oth <- mild_pool * p$p_mild_death_other
    st$mild <- mild_pool - mild_to_severe - mild_death_oth
    
    # --- Severe RHD (heart failure) ---
    severe_pool     <- st$severe + mild_to_severe
    surg_candidates <- severe_pool * p$frac_severe_surg_elig
    surgeries       <- surg_candidates * cov_surg[i]
    op_deaths       <- surgeries * p$p_surg_op_mortality
    to_post         <- surgeries - op_deaths
    remain_severe   <- severe_pool - surgeries
    sev_death_rhd   <- remain_severe * p$p_severe_death * (1 - p$eff_hf_mgmt * cov_hf[i])
    sev_death_oth   <- remain_severe * p$p_severe_death_other
    st$severe <- remain_severe - sev_death_rhd - sev_death_oth
    
    # --- Post-surgery ---
    post_pool    <- st$post + to_post
    post_death_r <- post_pool * p$p_post_death_rhd
    post_death_o <- post_pool * p$p_post_death_other
    st$post <- post_pool - post_death_r - post_death_o
    
    rhd_deaths <- op_deaths + sev_death_rhd + post_death_r
    
    # --- costs this year ---
    n_screened <- costs$n_screened_per_year * cov_sap_asymp[i]  # screening scales with 2b coverage
    n_on_sap   <- mild_pool * cov_sap_asymp[i]
    c_screen   <- n_screened * costs$cost_screen_per_person
    c_sap      <- n_on_sap   * costs$cost_sap_per_year
    c_hf       <- remain_severe * cov_hf[i] * costs$cost_hf_per_year
    c_surg     <- surgeries * costs$cost_surgery
    cost_total <- c_screen + c_sap + c_hf + c_surg
    
    out[[i]] <- data.table(
      year = years[i], pop = au_pop[i],
      new_rhd = new_mild, mild = st$mild, severe = st$severe, post = st$post,
      rhd_deaths = rhd_deaths, surgeries = surgeries,
      n_on_sap = n_on_sap, n_screened = n_screened,
      c_screen = c_screen, c_sap = c_sap, c_hf = c_hf, c_surg = c_surg, cost = cost_total
    )
  }
  rbindlist(out)
}

# ==============================================================================
# 8. RUN reference vs secondary-prevention scale-up
# ==============================================================================
ref <- run_scenario(cov$sap_asymp_ref, cov$hf_ref, cov$surg_ref)
sap <- run_scenario(cov$sap_asymp_up,  cov$hf_ref, cov$surg_up)

