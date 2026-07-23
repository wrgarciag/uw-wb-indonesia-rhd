# ==============================================================================
# RHD secondary-prevention investment case: ECONOMIC EVALUATION (A/B/C/D)
# scripts/08_economic_evaluation.R
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT (the ONLY place economics live)
# ------------------------------------------------------------------------------
# Benefit-cost + cost-effectiveness of the SAP scale-up vs the reference scenario.
# ALL cost / monetary values live only here; parameters are read from the
# 00_run_all.R globals (with standalone fallbacks). It consumes 06's per-location
# outputs and applies unit costs to the STAGE-AND-FLOW trace.
#
#   INPUTS
#     output/out_model/<location>.rds        (from 06)
#       $stages — A/B/C/D stock-and-flow trace for costing: n_screened,
#                 n_on_optimal_treatment, total_surgeries, by scenario x age x sex x year;
#       $wsd    — well-sick-dead table for RHD deaths (`dead`), by scenario x ... .
#     Economic PARAMETERS (globals from 00_run_all.R; fallbacks if absent):
#       disc_rate, gdp_pc_base, gdp_pc_base_year, gdp_growth, vsl_mult,
#       dalys_per_death, and unit costs cost_screen_per_person, cost_sap_per_year,
#       cost_surgery.
#
#   OUTPUTS (written to output/tables/):
#     rhd_budget_impact.csv     per location x year: cost components (screen/SAP/
#                               surgery), ref/sap totals, incremental cost.
#     rhd_economic_summary.csv  per location + TOTAL: incremental cost, monetised
#                               benefit, net benefit, BCR, cost per death averted,
#                               cost per DALY averted (all discounted).
#     rhd_economic_results.rds  list bundle of the above + parameters used.
#
# COSTS (screening applies to the TOTAL population screened; treatment to the
#        living-RHD stock on optimal treatment; surgery to the surgery trace):
#   c_screen  = sum(n_screened)             x cost_screen_per_person
#   c_sap     = sum(n_on_optimal_treatment) x cost_sap_per_year
#   c_surgery = sum(total_surgeries)        x cost_surgery
#   cost      = c_screen + c_sap + c_surgery
#   (HF-management cost is retired — HF management is not a modelled intervention.
#    Surgery coverage/effects are held equal in both arms, so surgery is a
#    background cost; any incremental surgery cost is only the second-order effect
#    of scale-up changing the C/D stocks.)
#
# METHOD (discounted; matches the established investment-case calculation):
#   incremental cost   = cost(sap) - cost(ref)                        [per year]
#   deaths averted     = RHD deaths(ref) - RHD deaths(sap)            [per year]
#   monetised benefit  = deaths averted x (VSL + GDP-pc productivity)  [per year]
#     VSL_y   = vsl_mult x GDP-pc_y ;  GDP-pc_y grows at gdp_growth from base year
#   discount factor    = 1 / (1 + disc_rate)^(year - discount_base_year)
#   BCR = sum(disc benefit) / sum(disc cost) ; net = benefit - cost
#   cost per death averted = sum(disc cost) / sum(undiscounted deaths averted)
#   cost per DALY averted  = sum(disc cost) / (sum(disc deaths averted) x DALYs/death)
# ==============================================================================

library(data.table)

if (!exists("wd_outp")) wd_outp <- paste0(here::here("output"), "/")
IN_DIR  <- paste0(wd_outp, "out_model/")
OUT_DIR <- paste0(wd_outp, "tables/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 0. ECONOMIC PARAMETERS  (read from 00_run_all.R globals; else standalone default)
#    [PAPER]/[LIT]/[CALIBRATE] tags document provenance; edit these in 00.
# ------------------------------------------------------------------------------
getp <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

econ <- list(
  disc_rate        = getp("disc_rate",        0.03),   # [PAPER] 3% discount, costs & benefits
  gdp_pc_base      = getp("gdp_pc_base",       4150),   # [LIT] GDP pc base-yr US$ (4150 = Indonesia standalone fallback; 00 sets per COUNTRY)
  gdp_pc_base_year = getp("gdp_pc_base_year",  2019),   # [LIT] year gdp_pc_base refers to
  gdp_growth       = getp("gdp_growth",        0.03),   # [CALIBRATE] real per-capita growth
  vsl_mult         = getp("vsl_mult",          30),     # [PAPER] VSL ~ 30 x GDP per capita
  dalys_per_death  = getp("dalys_per_death",   30),     # [CALIBRATE] undiscounted DALYs per RHD death
  cost_screen_per_person = getp("cost_screen_per_person", 1.10),  # [PAPER] echo screening cost/person
  cost_sap_per_year      = getp("cost_sap_per_year",      110),   # [PAPER] annual optimal-treatment cost/person
  cost_surgery           = getp("cost_surgery",           9000)   # [LIT] valve surgery + first-yr post-op
)

# validate parameters
with(econ, {
  if (disc_rate < 0 || disc_rate > 0.2) stop("disc_rate outside plausible [0,0.2].", call. = FALSE)
  if (any(c(gdp_pc_base, vsl_mult, dalys_per_death,
            cost_screen_per_person, cost_sap_per_year, cost_surgery) < 0))
    stop("A negative economic parameter was supplied.", call. = FALSE)
})

message("── 08_economic_evaluation.R : economic parameters ─────")
message(sprintf("  disc=%.1f%% | GDPpc(%d)=$%s g=%.1f%% | VSL=%dx | DALYs/death=%d",
                100 * econ$disc_rate, econ$gdp_pc_base_year,
                formatC(econ$gdp_pc_base, format = "d", big.mark = ","),
                100 * econ$gdp_growth, econ$vsl_mult, econ$dalys_per_death))
message(sprintf("  unit costs (US$): screen=%.2f/person  SAP=%d/yr  surgery=%d",
                econ$cost_screen_per_person, econ$cost_sap_per_year, econ$cost_surgery))

# ------------------------------------------------------------------------------
# 1. LOAD model outputs (stage-flow volumes + WSD deaths)
# ------------------------------------------------------------------------------
files <- list.files(IN_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0)
  stop("No model outputs in ", IN_DIR, ".\n  Run 06_run_prevention_model.R first.",
       call. = FALSE)

stg <- rbindlist(lapply(files, function(f) as.data.table(readRDS(f)$stages)),
                 use.names = TRUE, fill = TRUE)
wsd <- rbindlist(lapply(files, function(f) as.data.table(readRDS(f)$wsd)),
                 use.names = TRUE, fill = TRUE)
if (is.null(stg) || !nrow(stg))
  stop("Model output has no $stages table — re-run 06_run_prevention_model.R.", call. = FALSE)

# expected model horizon (ANALYSIS years) from 06's persisted meta (asserted below).
analysis_years <- sort(as.integer(readRDS(files[1])$meta$years))

# ------------------------------------------------------------------------------
# 2. ANNUAL COST COMPONENTS per location x scenario x year (from stage-flow trace)
#    (No cost columns exist in the model output — costs are applied HERE only.)
# ------------------------------------------------------------------------------
cost_ty <- stg[, .(
  c_screen  = sum(n_screened)             * econ$cost_screen_per_person,
  c_sap     = sum(n_on_optimal_treatment) * econ$cost_sap_per_year,
  c_surgery = sum(total_surgeries)        * econ$cost_surgery
), by = .(location, scenario, year)]
cost_ty[, cost := c_screen + c_sap + c_surgery]

# RHD deaths per location x scenario x year (from WSD `dead`)
death_ty <- wsd[, .(rhd_deaths = sum(dead)), by = .(location, scenario, year)]

# ------------------------------------------------------------------------------
# 3. INCREMENTAL (sap - ref) per location x year, then discounting + benefits
# ------------------------------------------------------------------------------
wide_cost  <- dcast(cost_ty,  location + year ~ scenario, value.var = "cost")
wide_death <- dcast(death_ty, location + year ~ scenario, value.var = "rhd_deaths")
if (!all(c("ref", "sap") %in% names(wide_cost)))
  stop("Both 'ref' and 'sap' scenarios must be present in the model output.", call. = FALSE)

inc <- merge(wide_cost[,  .(location, year, cost_ref = ref, cost_sap = sap)],
             wide_death[, .(location, year, death_ref = ref, death_sap = sap)],
             by = c("location", "year"))
inc[, `:=`(inc_cost     = cost_sap - cost_ref,
           deaths_avert = death_ref - death_sap)]

# Economic discount reference year: set in 00_run_all.R (default ANALYSIS_YEAR_START);
# standalone fallback = the first analysis year in the data. Separate from the GDP-pc
# reference year (econ$gdp_pc_base_year) and the GBD rate reference year.
discount_base_year <- getp("discount_base_year", min(inc$year))
inc[, `:=`(
  df      = 1 / (1 + econ$disc_rate)^(year - discount_base_year),
  gdp_pc  = econ$gdp_pc_base * (1 + econ$gdp_growth)^(year - econ$gdp_pc_base_year)
)]
# NOTE: create `vsl` in its OWN := call before referencing it (columns created
# within a single := are not visible to sibling RHS expressions in that call).
inc[, vsl := econ$vsl_mult * gdp_pc]
inc[, `:=`(benefit_vsl = deaths_avert * vsl,
           benefit_gdp = deaths_avert * gdp_pc)]
inc[, benefit_total := benefit_vsl + benefit_gdp]

# ------------------------------------------------------------------------------
# 4. SUMMARISE per location, plus a pooled TOTAL  (safe zero-denominator guards)
# ------------------------------------------------------------------------------
safe_ratio <- function(num, den) if (is.finite(den) && den > 0) num / den else NA_real_

summarise_econ <- function(d, label) {
  tot_cost_disc    <- sum(d$inc_cost      * d$df)
  tot_benefit_disc <- sum(d$benefit_total * d$df)
  tot_deaths_avert <- sum(d$deaths_avert)
  disc_deaths      <- sum(d$deaths_avert  * d$df)
  data.table(
    location        = label,
    total_cost      = tot_cost_disc,
    total_benefit   = tot_benefit_disc,
    net_benefit     = tot_benefit_disc - tot_cost_disc,
    bcr             = safe_ratio(tot_benefit_disc, tot_cost_disc),
    deaths_averted  = tot_deaths_avert,
    cost_per_death  = safe_ratio(tot_cost_disc, tot_deaths_avert),
    cost_per_daly   = safe_ratio(tot_cost_disc, disc_deaths * econ$dalys_per_death)
  )
}

summary_by_loc <- rbindlist(lapply(split(inc, by = "location", keep.by = TRUE),
                                   function(d) summarise_econ(d, d$location[1])))
summary_total  <- summarise_econ(inc, "TOTAL")
econ_summary   <- rbindlist(list(summary_by_loc, summary_total), use.names = TRUE)

# budget-impact table (per location x year)
budget_impact <- inc[, .(location, year,
                         cost_ref, cost_sap, incremental = inc_cost,
                         deaths_avert)]
budget_impact <- merge(
  budget_impact,
  cost_ty[scenario == "sap", .(location, year, c_screen, c_sap, c_surgery)],
  by = c("location", "year"), all.x = TRUE)

# ------------------------------------------------------------------------------
# 5. VALIDATION  (fail loudly; deaths-averted / incremental-cost checks CUMULATIVE)
# ------------------------------------------------------------------------------
if (any(cost_ty$cost < -1e-6)) stop("Negative annual cost computed.", call. = FALSE)

# analysis-horizon integrity: the incremental (and hence budget/summary) table must
# span EXACTLY the analysis years — no economic row outside the configured window.
inc_years <- sort(unique(as.integer(inc$year)))
if (!identical(inc_years, analysis_years))
  stop("Economic table years (", min(inc_years), "-", max(inc_years),
       ") do not match the analysis horizon (", min(analysis_years), "-",
       max(analysis_years), ").", call. = FALSE)
if (!is.finite(discount_base_year))
  stop("discount_base_year is not finite.", call. = FALSE)

# per-year deaths-averted / incremental cost MAY dip negative in late years (a
# mortality-only intervention raises surviving prevalence, so later cohorts can
# briefly cost more or shift deaths); the investment-case guarantee is CUMULATIVE.
neg_da_years <- inc[deaths_avert < -1e-6, .N]
neg_ic_years <- inc[inc_cost     < -1e-6, .N]
if (neg_da_years > 0)
  message(sprintf("  NOTE: %d location-year(s) with negative per-year deaths averted (cumulative still checked).",
                  neg_da_years))
if (neg_ic_years > 0)
  message(sprintf("  NOTE: %d location-year(s) with negative per-year incremental cost.", neg_ic_years))

cum_da <- inc[, sum(deaths_avert)]
cum_ic <- inc[, sum(inc_cost)]
if (cum_da < -1e-6)
  stop("Cumulative deaths averted is negative (SAP worse than reference overall).", call. = FALSE)
if (cum_ic < -1e-6)
  stop("Cumulative incremental cost is negative (SAP cheaper than reference overall) — unexpected.",
       call. = FALSE)

if (anyNA(econ_summary[, .(total_cost, total_benefit, net_benefit, deaths_averted)]))
  stop("Economic summary contains NA in a required field.", call. = FALSE)
if (summary_total$total_cost <= 0) stop("Total incremental cost is non-positive.", call. = FALSE)
if (!is.finite(summary_total$bcr) || summary_total$bcr <= 0)
  stop("Benefit-cost ratio is non-finite or non-positive.", call. = FALSE)

# ------------------------------------------------------------------------------
# 6. WRITE + report
# ------------------------------------------------------------------------------
fwrite(budget_impact, paste0(OUT_DIR, "rhd_budget_impact.csv"))
fwrite(econ_summary,  paste0(OUT_DIR, "rhd_economic_summary.csv"))
saveRDS(list(summary = econ_summary, budget_impact = budget_impact,
             incremental = inc, params = econ,
             discount_base_year = discount_base_year),
        paste0(OUT_DIR, "rhd_economic_results.rds"))

t <- summary_total
message("── ECONOMIC RESULTS (discounted; base-year US$) ───────")
message(sprintf("  Total incremental cost : $%.2f billion", t$total_cost   / 1e9))
message(sprintf("  Total monetised benefit: $%.2f billion", t$total_benefit / 1e9))
message(sprintf("  Net benefit            : $%.2f billion", t$net_benefit  / 1e9))
message(sprintf("  Benefit-cost ratio     : %.2f", t$bcr))
message(sprintf("  RHD deaths averted     : %s", formatC(round(t$deaths_averted), format = "d", big.mark = ",")))
message(sprintf("  Cost per death averted : $%s", formatC(round(t$cost_per_death), format = "d", big.mark = ",")))
message(sprintf("  Cost per DALY averted  : $%s", formatC(round(t$cost_per_daly),  format = "d", big.mark = ",")))
message("── 08_economic_evaluation.R complete ──────────────────")
message("  Wrote: ", OUT_DIR, "rhd_budget_impact.csv, rhd_economic_summary.csv, rhd_economic_results.rds")
