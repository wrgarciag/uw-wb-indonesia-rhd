# ==============================================================================
# RHD secondary-prevention investment case: cohort state-transition model (data.table)
#
# Replicates the structure of:
#   Coates et al. "An investment case for the prevention and management of
#   rheumatic heart disease in the African Union 2021-30." Lancet Glob Health 2021.
#   https://pmc.ncbi.nlm.nih.gov/articles/PMC9087136/
#
# Focus of THIS script: scale-up of SECONDARY PREVENTION only, defined as
#   (a) RHD screening (echocardiography) to detect asymptomatic (mild) RHD, and
#   (b) secondary antibiotic prophylaxis (SAP, monthly benzathine penicillin G)
#       given to (i) people after ARF and (ii) screen-detected asymptomatic RHD.
#
# Outputs: incident RHD, severe RHD, RHD deaths; costs; monetised benefits
#          (full-income); incremental cost-effectiveness (cost per death and per
#          DALY averted); benefit-cost ratio; and a year-by-year budget impact.
#
# IMPORTANT ON PARAMETERS
#   The paper's exact per-unit costs and transition probabilities live in an
#   appendix that is not reproduced in the main text. Every number the main
#   article gives (Table 1 effect sizes and coverages, headline results) is used
#   directly and tagged  [PAPER]. Everything else is a plausible literature value
#   tagged  [LIT]  (with a source in the comment) or  [CALIBRATE]  where you
#   should tune it to your setting / the appendix. Change values only in the
#   PARAMETERS block; the engine below is generic.
# ==============================================================================

library(data.table)
set.seed(1)

# ------------------------------------------------------------------------------
# 1. TIME, DISCOUNTING, DEMOGRAPHY----
# ------------------------------------------------------------------------------
years      <- 2021:2030          # [PAPER] 10-year scale-up window
n_years    <- length(years)
disc_rate  <- 0.03               # [PAPER] 3% discount on costs and benefits

# Simple AU demographic backdrop (denominators for rates only).
pop_2021   <- 1.36e9             # [LIT] UN WPP 2019, AU ~1.36bn in 2021
pop_growth <- 0.025              # [LIT] ~2.5%/yr AU population growth
au_pop     <- pop_2021 * (1 + pop_growth)^(0:(n_years - 1))

gdp_pc_2019   <- 1800            # [LIT] AU GDP per capita, 2019 US$ (~World Bank)
gdp_growth    <- 0.03            # [CALIBRATE] real per-capita GDP growth for VSL/GDP benefit

# ------------------------------------------------------------------------------
# 2. DISEASE PARAMETERS (annual transition probabilities)----
#    Sub-model of the RHD-affected population. New disease enters via ARF and via
#    incident asymptomatic RHD; secondary prevention acts on the two progression
#    steps 2a (post-ARF) and 2b (asymptomatic RHD) from Table 1.
# ------------------------------------------------------------------------------
p <- list(

  # --- incidence inflows (annual, AU-wide) ---
  arf_incident_2021      = 340000,   # [CALIBRATE] new ARF cases/yr; tune to GBD ARF burden
  rhd_incident_2021      = 250000,   # [CALIBRATE] new asymptomatic RHD/yr; tune so that primary
                                     #   prevention's 7.6% incidence drop ~ paper's 187k averted
  incidence_trend        = 0.985,    # [PAPER-ADJ] ~15% RHD incidence decline over decade from
                                     #   improving living conditions -> ~1.5%/yr (Table 1 footnote)

  # --- progression / mortality (baseline, i.e. untreated) ---
  p_arf_to_rhd           = 0.30,     # [LIT] fraction of ARF episodes progressing to RHD
  p_arf_death            = 0.010,    # [LIT] acute case fatality of ARF
  p_arf_resolve          = 0.69,     # remainder -> remission/healthy (1 - the two above)

  p_mild_to_severe       = 0.045,    # [CALIBRATE] asymptomatic -> heart-failure/yr (subclinical
                                     #   RHD progression cohorts; tune to appendix)
  p_mild_death_other     = 0.006,    # [LIT] background mortality of mild-RHD population

  p_severe_death         = 0.090,    # [LIT] REMEDY ~17% mortality/2yr -> ~9%/yr untreated severe
  p_severe_death_other   = 0.010,    # [LIT] competing background mortality
  p_surg_op_mortality    = 0.03,     # [PAPER] 3% operative mortality (Table 1 footnote)

  p_post_death_rhd       = 0.020,    # [LIT] residual RHD mortality after valve surgery
  p_post_death_other     = 0.010,    # [LIT] background mortality post-surgery

  # --- intervention EFFECT SIZES (relative risk reductions) ---
  eff_sap_postarf        = 0.55,     # [PAPER] Table 1 #2a: SAP after ARF, 55% (8-78)
  eff_sap_asymp          = 0.55,     # [PAPER] Table 1 #2b: SAP in asymptomatic RHD, 55% (7-78)
  eff_hf_mgmt            = 0.60,     # [PAPER] Table 1 #3: HF management, 60% (30-80)
  eff_surgery            = 0.85,     # [PAPER] Table 1 #4: surgery, 85% (70-92)

  # --- eligibility fractions ---
  frac_severe_surg_elig  = 0.55      # [CALIBRATE] share of severe RHD aged 10-40 (surgery-eligible)
)

# ------------------------------------------------------------------------------
# 3. COVERAGE TRAJECTORIES (linear scale-up from baseline 2020 to target 2030)
#    Baselines and targets from Table 1. Secondary prevention = SAP (2a & 2b);
#    "screening" governs how much asymptomatic RHD is *found* so SAP can be given.
# ------------------------------------------------------------------------------
# Reference (no scale-up): coverage held at baseline for all years.
# Scale-up: linear ramp baseline -> target across 2021..2030.
ramp <- function(baseline, target) {
  # value in each modelled year (target reached in the final year)
  baseline + (target - baseline) * (seq_len(n_years) - 1) / (n_years - 1)
}

cov <- list(
  # Secondary prevention (SAP) after ARF -- Table 1 #2a: 5.0% -> 40%
  sap_postarf_ref  = rep(0.05, n_years),
  sap_postarf_up   = ramp(0.05, 0.40),

  # Screening-enabled SAP for asymptomatic RHD -- Table 1 #2b: 5.0% -> 40%
  # Interpreted as (screening detection) x (SAP uptake among detected).
  sap_asymp_ref    = rep(0.05, n_years),
  sap_asymp_up     = ramp(0.05, 0.40),

  # The following are held at BASELINE here because this script isolates
  # secondary prevention. Flip *_up to the ramp() versions to add integrated care.
  hf_ref   = rep(0.08, n_years),  hf_up   = rep(0.08, n_years),  # Table 1 #3: 8% -> (55%)
  surg_ref = rep(0.05, n_years),  surg_up = rep(0.05, n_years)   # Table 1 #4: 5% -> (25%)
)

# ------------------------------------------------------------------------------
# 4. UNIT COSTS (2019 US$)  --  all [CALIBRATE]/[LIT]; replace with appendix values
# ------------------------------------------------------------------------------
costs <- list(
  cost_screen_per_person = 12,     # [LIT] echo screening cost/person (handheld/decentralised)
  n_screened_per_year    = 8e6,    # [CALIBRATE] target population echo-screened per year at
                                   #   full coverage (e.g. school programmes); scaled by coverage
  cost_sap_per_year      = 45,     # [LIT] annual SAP cost/person (monthly BPG + visits, WHO-CHOICE)
  cost_hf_per_year       = 120,    # [LIT] annual HF medical management cost/person
  cost_surgery           = 9000,   # [LIT] valve surgery + first-year post-op cost/person
  cost_arf_hospitalisation = 200   # [LIT] ARF admission cost (for averted-hospitalisation benefit)
)

# ------------------------------------------------------------------------------
# 5. VALUE OF A STATISTICAL LIFE (full-income benefit)
# ------------------------------------------------------------------------------
# Paper: adjust a US VSL to the AU using GNI ratio and income elasticity, and
# notes VSL ends up ~30x GDP per capita. We follow the shortcut VSL = mult x GDPpc.
vsl_mult      <- 30              # [PAPER] "VSL is roughly 30 times GDP per capita"
dalys_per_death <- 30            # [CALIBRATE] undiscounted DALYs per RHD/ARF death (young deaths)

# ==============================================================================
# 6. MODEL ENGINE  -- runs one scenario, returns a data.table by year
# ==============================================================================
run_scenario <- function(cov_sap_postarf, cov_sap_asymp, cov_hf, cov_surg) {

  # living state stocks carried across cycles
  st <- list(ARF = 0, mild = 0, severe = 0, post = 0)

  out <- vector("list", n_years)

  for (i in seq_len(n_years)) {
    tr <- p$incidence_trend^(i - 1)                    # incidence trend multiplier
    new_arf  <- p$arf_incident_2021 * tr
    new_mild <- p$rhd_incident_2021 * tr               # incident asymptomatic RHD inflow

    # --- ARF cohort this year (new + carried) ---
    arf_pool <- st$ARF + new_arf
    # SAP after ARF reduces progression ARF -> RHD
    prog_arf <- p$p_arf_to_rhd * (1 - p$eff_sap_postarf * cov_sap_postarf[i])
    arf_to_rhd  <- arf_pool * prog_arf
    arf_deaths  <- arf_pool * p$p_arf_death
    # remainder resolves (no long-term ARF stock in this parsimonious version)
    st$ARF <- 0

    # --- Mild (asymptomatic) RHD ---
    mild_pool <- st$mild + new_mild + arf_to_rhd
    # SAP in screen-detected asymptomatic RHD reduces mild -> severe
    prog_mild <- p$p_mild_to_severe * (1 - p$eff_sap_asymp * cov_sap_asymp[i])
    mild_to_severe <- mild_pool * prog_mild
    mild_death_oth <- mild_pool * p$p_mild_death_other
    st$mild <- mild_pool - mild_to_severe - mild_death_oth

    # --- Severe RHD (heart failure) ---
    severe_pool <- st$severe + mild_to_severe
    # surgery among eligible severe cases
    surg_candidates <- severe_pool * p$frac_severe_surg_elig
    surgeries       <- surg_candidates * cov_surg[i]
    op_deaths       <- surgeries * p$p_surg_op_mortality
    to_post         <- surgeries - op_deaths
    # HF management reduces death among remaining severe cases
    remain_severe   <- severe_pool - surgeries
    sev_death_rhd   <- remain_severe * p$p_severe_death * (1 - p$eff_hf_mgmt * cov_hf[i])
    sev_death_oth   <- remain_severe * p$p_severe_death_other
    st$severe <- remain_severe - sev_death_rhd - sev_death_oth

    # --- Post-surgery ---
    post_pool     <- st$post + to_post
    post_death_r  <- post_pool * p$p_post_death_rhd
    post_death_o  <- post_pool * p$p_post_death_other
    st$post <- post_pool - post_death_r - post_death_o

    rhd_deaths <- op_deaths + sev_death_rhd + post_death_r

    # --- costs this year ---
    n_screened <- costs$n_screened_per_year * (cov_sap_asymp[i])   # screening scaled w/ 2b coverage
    n_on_sap   <- arf_pool * cov_sap_postarf[i] + mild_pool * cov_sap_asymp[i]
    c_screen   <- n_screened * costs$cost_screen_per_person
    c_sap      <- n_on_sap   * costs$cost_sap_per_year
    c_hf       <- remain_severe * cov_hf[i] * costs$cost_hf_per_year
    c_surg     <- surgeries * costs$cost_surgery
    cost_total <- c_screen + c_sap + c_hf + c_surg

    out[[i]] <- data.table(
      year = years[i],
      pop  = au_pop[i],
      new_rhd = new_mild + arf_to_rhd,
      mild = st$mild, severe = st$severe, post = st$post,
      arf_deaths = arf_deaths, rhd_deaths = rhd_deaths,
      surgeries = surgeries,
      n_on_sap = n_on_sap, n_screened = n_screened,
      c_screen = c_screen, c_sap = c_sap, c_hf = c_hf, c_surg = c_surg,
      cost = cost_total
    )
  }
  rbindlist(out)
}

# ==============================================================================
# 7. RUN reference vs secondary-prevention scale-up
# ==============================================================================
ref <- run_scenario(cov$sap_postarf_ref, cov$sap_asymp_ref, cov$hf_ref, cov$surg_ref)
sap <- run_scenario(cov$sap_postarf_up,  cov$sap_asymp_up,  cov$hf_ref, cov$surg_up)

# incremental (scale-up minus reference)
inc <- data.table(
  year            = years,
  rhd_deaths_avert = ref$rhd_deaths - sap$rhd_deaths,
  arf_deaths_avert = ref$arf_deaths - sap$arf_deaths,
  rhd_cases_avert  = ref$new_rhd   - sap$new_rhd,
  inc_cost         = sap$cost - ref$cost
)
inc[, deaths_avert := rhd_deaths_avert + arf_deaths_avert]

# discount factors (base year 2021)
df <- 1 / (1 + disc_rate)^(years - years[1])

# ------------------------------------------------------------------------------
# 8. MONETISED BENEFITS (full income) & COST-EFFECTIVENESS
# ------------------------------------------------------------------------------
gdp_pc  <- gdp_pc_2019 * (1 + gdp_growth)^(years - 2019)
vsl     <- vsl_mult * gdp_pc

benefit_vsl  <- inc$deaths_avert * vsl        # value of health (dominant term, ~93% in paper)
benefit_gdp  <- inc$deaths_avert * gdp_pc     # productivity proxy (population/GDP difference)

# Averted ARF hospitalisations. NOTE: because ARF incidence is exogenous here,
# SAP prevents PROGRESSION but not recurrent ARF, so this term is ~0 in this
# parsimonious version. Set arf_hosp_frac > 0 if you model ARF recurrence.
arf_hosp_frac  <- 0.0                          # [CALIBRATE] averted ARF admissions per averted case
arf_hosp_avert <- inc$rhd_cases_avert * arf_hosp_frac
benefit_hosp   <- arf_hosp_avert * costs$cost_arf_hospitalisation

benefit_total <- benefit_vsl + benefit_gdp + benefit_hosp

# discounted aggregates
tot_cost_disc    <- sum(inc$inc_cost   * df)
tot_benefit_disc <- sum(benefit_total  * df)
tot_deaths_avert <- sum(inc$deaths_avert)                   # undiscounted (health outcomes)
tot_dalys_avert  <- tot_deaths_avert * dalys_per_death

bcr        <- tot_benefit_disc / tot_cost_disc
net_ben    <- tot_benefit_disc - tot_cost_disc
cost_per_death <- tot_cost_disc / tot_deaths_avert
cost_per_daly  <- tot_cost_disc / (sum(inc$deaths_avert * df) * dalys_per_death)

# ==============================================================================
# 9. REPORT
# ==============================================================================
cat("\n============================================================\n")
cat(" RHD SECONDARY PREVENTION (screening + SAP) : 2021-2030 scale-up\n")
cat("============================================================\n\n")

cat("--- BUDGET IMPACT (incremental cost of scale-up, US$ millions) ---\n")
bi <- sap[, .(year,
              screening = c_screen/1e6, SAP = c_sap/1e6,
              HF = c_hf/1e6, surgery = c_surg/1e6, total = cost/1e6)]
bi_ref <- ref[, .(year, total_ref = cost/1e6)]
bi <- merge(bi, bi_ref, by = "year")
bi[, incremental := total - total_ref]
print(bi[, lapply(.SD, function(x) if (is.numeric(x)) round(x,1) else x)])

cat("\n--- HEALTH IMPACT (cumulative 2021-2030) ---\n")
cat(sprintf("  RHD deaths averted        : %s\n", format(round(sum(inc$rhd_deaths_avert)), big.mark=",")))
cat(sprintf("  ARF deaths averted        : %s\n", format(round(sum(inc$arf_deaths_avert)), big.mark=",")))
cat(sprintf("  Total deaths averted      : %s\n", format(round(tot_deaths_avert), big.mark=",")))
cat(sprintf("  RHD incident cases averted: %s\n", format(round(sum(inc$rhd_cases_avert)), big.mark=",")))

cat("\n--- ECONOMIC RESULTS (3% discounting; 2019 US$) ---\n")
cat(sprintf("  Total incremental cost    : $%.2f billion\n", tot_cost_disc/1e9))
cat(sprintf("  Total monetised benefit   : $%.2f billion\n", tot_benefit_disc/1e9))
cat(sprintf("  Net benefit               : $%.2f billion\n", net_ben/1e9))
cat(sprintf("  Benefit-cost ratio        : %.2f\n", bcr))
cat(sprintf("  Cost per death averted    : $%s\n", format(round(cost_per_death), big.mark=",")))
cat(sprintf("  Cost per DALY averted      : $%s\n", format(round(cost_per_daly), big.mark=",")))
cat("\n(Compare with paper's integrated secondary+tertiary care: BCR 4.7 to 2030,\n")
cat(" ~59,500 RHD deaths averted, cost per death averted ~$14,800.)\n\n")

# tidy machine-readable output
results <- list(budget_impact = bi, incremental = inc,
                summary = data.table(bcr = bcr, net_benefit = net_ben,
                                     total_cost = tot_cost_disc,
                                     total_benefit = tot_benefit_disc,
                                     deaths_averted = tot_deaths_avert,
                                     cost_per_death = cost_per_death,
                                     cost_per_daly = cost_per_daly))
saveRDS(results, "rhd_results.rds")
fwrite(bi,  "rhd_budget_impact.csv")
fwrite(inc, "rhd_incremental.csv")
cat("Saved: rhd_budget_impact.csv, rhd_incremental.csv, rhd_results.rds\n")
