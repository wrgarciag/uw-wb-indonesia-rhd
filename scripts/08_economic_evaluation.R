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
#       disc_rate, gni_pc_base + gni_pc_base_year (REQUIRED: no fallback, never
#       replaced by GDP), gdp_pc_base, gdp_pc_base_year, gdp_growth, vsl_mult,
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
# COSTS (each service volume from 06 is costed ONCE, never re-derived here):
#   screening : n_screened = screening ENCOUNTERS (one screen per person across ages
#               5-39 among people without a prior confirmed diagnosis: annual probability
#               x never-screened people; people without RHD included);
#   SAP       : n_on_optimal_treatment = on-BPG PERSON-YEARS of SAP-eligible confirmed
#               disease (A at age <= 20, B; start of cycle, after the care update) — NOT
#               the adherence-weighted protected volume (adherence enters the health
#               effect only). Adult B on BPG is costed although its progression HR is 1;
#   surgery   : the surgery trace (service; equal coverage in both arms).
#   UNCOSTED  : confirmatory echoes (n_confirm_completed) — no unit cost exists in 00
#               block J; their volume is reported and flagged, never priced here, so
#               total programme cost is understated by the confirmation cost.
#   c_screen  = sum(n_screened)             x cost_screen_per_person
#   c_sap     = sum(n_on_optimal_treatment) x cost_sap_per_year
#   c_surgery = sum(total_surgeries)        x cost_surgery
#   cost      = c_screen + c_sap + c_surgery
#   (HF-management cost is retired — HF management is not a modelled intervention.
#    Any incremental surgery cost is only the second-order effect of scale-up
#    changing the C/D stocks.)
#
# METHOD (discounted; incremental, not gross, programme costing):
#   incremental cost   = cost(sap) - cost(ref)                        [per year]
#   deaths averted     = RHD deaths(ref) - RHD deaths(sap)            [per year; signed, never clamped]
#   monetised benefit  = deaths averted x VSL_t                       [per year; VSL ONLY]
#     VSL_b = vsl_mult x GNI-pc_b  (b = gni_pc_base_year, the income anchor)
#     VSL_t = VSL_b x GDP-pc_t / GDP-pc_b, GDP-pc_t = gdp_pc_base x (1 + gdp_growth)^(t - gdp_pc_base_year)
#           = VSL_b x (1 + gdp_growth)^(t - b)   (GDP pc is a growth index only)
#     benefit_gdp is a RETIRED compatibility column, identically 0 (the former
#     additive one-year GDP-per-capita benefit is removed); benefit_total = benefit_vsl.
#   discount factor    = 1 / (1 + disc_rate)^(year - discount_base_year)
#   BCR = sum(disc benefit) / sum(disc incremental cost), defined ONLY when the
#         denominator is > 0 (NA otherwise: cost-neutral / cost-saving is reported
#         in incremental_cost_status, never given a fabricated ratio); net = benefit - cost
#   cost per death averted = sum(disc cost) / sum(undiscounted deaths averted)
#   cost per DALY averted  = sum(disc cost) / (sum(disc deaths averted) x DALYs/death)
#   TOTAL row pools benefits and costs across locations (not a mean of ratios).
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

# GNI per capita anchors VSL. REQUIRED, country-specific, no standalone fallback
# and never substituted by GDP per capita (00_run_all.R block J sets it per COUNTRY).
for (nm in c("gni_pc_base", "gni_pc_base_year"))
  if (!exists(nm, inherits = TRUE) || is.null(get(nm, inherits = TRUE)) ||
      length(get(nm, inherits = TRUE)) != 1L || !is.finite(get(nm, inherits = TRUE)))
    stop("Required economic parameter '", nm, "' is missing or not a finite scalar. Set the ",
         "country's GNI per capita (World Bank WDI NY.GNP.PCAP.CD) and its year in 00_run_all.R ",
         "block J; VSL is anchored on GNI per capita and never falls back to GDP per capita.", call. = FALSE)

econ <- list(
  disc_rate        = getp("disc_rate",        0.03),   # [PAPER] 3% discount, costs & benefits
  gni_pc_base      = getp("gni_pc_base",       NA_real_),     # [WDI] GNI pc (Atlas, US$) at the VSL anchor year (REQUIRED)
  gni_pc_base_year = as.integer(getp("gni_pc_base_year", NA_integer_)),  # GNI / VSL income anchor year b
  gni_pc_source    = getp("gni_pc_source",     "not recorded"),
  gdp_pc_base      = getp("gdp_pc_base",       4150),   # GDP pc path level (growth INDEX only; level cancels in VSL_t)
  gdp_pc_base_year = getp("gdp_pc_base_year",  2019),   # year gdp_pc_base refers to
  gdp_growth       = getp("gdp_growth",        0.03),   # [ASSUMPTION] constant real per-capita GDP growth
  vsl_mult         = getp("vsl_mult",          30),     # [PROVISIONAL] VSL = 30 x GNI per capita
  dalys_per_death  = getp("dalys_per_death",   30),     # [CALIBRATE] undiscounted DALYs per RHD death
  cost_screen_per_person = getp("cost_screen_per_person", 1.10),  # [PAPER] echo screening cost per encounter
  cost_sap_per_year      = getp("cost_sap_per_year",      110),   # [PAPER] cost per on-BPG person-year
  cost_surgery           = getp("cost_surgery",           9000),  # [LIT] valve surgery + first-yr post-op
  valuation_method = "VSL only: VSL_t = vsl_mult x GNI-pc_b x GDP-pc_t / GDP-pc_b; no additive GDP-per-capita benefit",
  uncosted_items = "confirmatory echo after a true screen positive (n_confirm_completed): no unit cost in 00 block J; omitted from total costs"
)

# validate parameters
with(econ, {
  if (disc_rate < 0 || disc_rate > 0.2) stop("disc_rate outside plausible [0,0.2].", call. = FALSE)
  if (any(c(gni_pc_base, gdp_pc_base, vsl_mult, dalys_per_death,
            cost_screen_per_person, cost_sap_per_year, cost_surgery) < 0))
    stop("A negative economic parameter was supplied.", call. = FALSE)
  if (gni_pc_base <= 0 || gdp_pc_base <= 0) stop("GNI and GDP per capita must be positive.", call. = FALSE)
  if (!is.finite(gdp_growth) || gdp_growth <= -1) stop("gdp_growth must be a finite rate > -1.", call. = FALSE)
  if (gni_pc_base_year < 1960L || gni_pc_base_year > 2100L)
    stop("gni_pc_base_year outside a plausible calendar range.", call. = FALSE)
})

message("── 08_economic_evaluation.R : economic parameters ─────")
message(sprintf("  disc=%.1f%% | VSL = %gx GNIpc(%d)=$%s -> VSL_%d=$%s, grown with GDP pc at g=%.1f%% | DALYs/death=%g",
                100 * econ$disc_rate, econ$vsl_mult, econ$gni_pc_base_year,
                formatC(econ$gni_pc_base, format = "d", big.mark = ","), econ$gni_pc_base_year,
                formatC(econ$vsl_mult * econ$gni_pc_base, format = "d", big.mark = ","),
                100 * econ$gdp_growth, econ$dalys_per_death))
message(sprintf("  unit costs (US$): screen=%.2f/encounter  SAP=%g/person-yr  surgery=%g",
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
# UNCOSTED volume (reported, flagged; never priced): completed confirmatory echoes
if (!"n_confirm_completed" %in% names(stg))
  stop("Model $stages lacks n_confirm_completed (re-run 06 under the first-diagnosis engine).", call. = FALSE)
conf_ty <- stg[, .(confirmations = sum(n_confirm_completed)), by = .(location, scenario, year)]

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
# reference year (econ$gdp_pc_base_year), the GNI/VSL anchor (econ$gni_pc_base_year)
# and the GBD rate reference year.
discount_base_year <- getp("discount_base_year", min(inc$year))
# GDP-per-capita path (growth index), ALIGNED to the GNI/VSL anchor year b before the
# ratio is taken: vsl_growth_index = GDPpc_t / GDPpc_b = (1 + g)^(t - b).
gdp_pc_at_anchor <- econ$gdp_pc_base * (1 + econ$gdp_growth)^(econ$gni_pc_base_year - econ$gdp_pc_base_year)
inc[, `:=`(
  df      = 1 / (1 + econ$disc_rate)^(year - discount_base_year),
  gdp_pc  = econ$gdp_pc_base * (1 + econ$gdp_growth)^(year - econ$gdp_pc_base_year)
)]
# NOTE: create each column in its OWN := call before referencing it (columns created
# within a single := are not visible to sibling RHS expressions in that call).
inc[, vsl := econ$vsl_mult * econ$gni_pc_base * (gdp_pc / gdp_pc_at_anchor)]
inc[, `:=`(benefit_vsl = deaths_avert * vsl,
           benefit_gdp = 0)]          # RETIRED compatibility field: the additive GDP-pc benefit is removed
inc[, benefit_total := benefit_vsl]
inc[, vsl_growth_index := gdp_pc / gdp_pc_at_anchor]   # additive column (appended; schema order kept)

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
    bcr             = safe_ratio(tot_benefit_disc, tot_cost_disc),   # NA (undefined) unless cost > 0
    deaths_averted  = tot_deaths_avert,
    cost_per_death  = safe_ratio(tot_cost_disc, tot_deaths_avert),
    cost_per_daly   = safe_ratio(tot_cost_disc, disc_deaths * econ$dalys_per_death),
    # additive: why bcr may be undefined (zero / negative discounted incremental cost)
    incremental_cost_status = if (!is.finite(tot_cost_disc)) "invalid" else if (tot_cost_disc > 0) "positive" else
      if (tot_cost_disc == 0) "cost_neutral" else "cost_saving"
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
# additive (appended): incremental confirmatory echoes, UNCOSTED (volume only; flagged)
conf_w <- dcast(conf_ty, location + year ~ scenario, value.var = "confirmations")
budget_impact <- merge(budget_impact, conf_w[, .(location, year, confirmations_incremental_uncosted = sap - ref)],
                       by = c("location", "year"), all.x = TRUE)

# ------------------------------------------------------------------------------
# 5. VALIDATION  (fail loudly; deaths-averted / incremental-cost checks CUMULATIVE)
# ------------------------------------------------------------------------------
if (any(cost_ty$cost < -1e-6)) stop("Negative annual cost computed.", call. = FALSE)
# service volumes costed here reconcile with 06's own diagnostics (each costed once;
# SAP cost uses on-BPG person-years, never the protection-weighted volume)
for (f in files) {
  o <- readRDS(f); d6 <- o$diag
  if (!is.null(d6$sap_trace_by_year) && !is.null(d6$screening_by_year)) {
    v <- merge(as.data.table(o$stages)[, .(sap_py = sum(n_on_optimal_treatment), enc = sum(n_screened)), by = .(scenario, year)],
               merge(as.data.table(d6$sap_trace_by_year)[, .(scenario, year, on_bpg_person_years)],
                     as.data.table(d6$screening_by_year)[, .(scenario, year, encounters)], by = c("scenario", "year")),
               by = c("scenario", "year"))
    if (v[, max(abs(sap_py - on_bpg_person_years) / pmax(on_bpg_person_years, 1), abs(enc - encounters) / pmax(encounters, 1))] > 1e-9)
      stop("Costed SAP person-years / screening encounters do not reconcile with 06 diagnostics (", basename(f), ").", call. = FALSE)
  }
}

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
# a cost-neutral or cost-saving scale-up is a legitimate result: reported, not a failure
if (cum_ic <= 0)
  message("  NOTE: cumulative undiscounted incremental cost is ", if (cum_ic == 0) "zero" else "negative",
          " (scale-up is cost-neutral or cost-saving overall).")

if (anyNA(econ_summary[, .(total_cost, total_benefit, net_benefit, deaths_averted)]) ||
    any(!is.finite(as.matrix(econ_summary[, .(total_cost, total_benefit, net_benefit, deaths_averted)]))))
  stop("Economic summary contains NA / non-finite values in a required field.", call. = FALSE)
for (i in seq_len(nrow(econ_summary))) {
  r <- econ_summary[i]
  if (r$total_cost > 0 && (!is.finite(r$bcr) || r$bcr <= 0))
    stop("Benefit-cost ratio for ", r$location, " is non-finite or non-positive despite positive incremental cost.",
         call. = FALSE)
  if (r$total_cost <= 0 && !is.na(r$bcr))
    stop("A BCR was reported for ", r$location, " although its discounted incremental cost is not positive.",
         call. = FALSE)
  if (r$total_cost <= 0)
    message(sprintf("  NOTE: %s discounted incremental cost is %s -> BCR undefined (%s).",
                    r$location, format(r$total_cost, big.mark = ","), r$incremental_cost_status))
}

# benefit / cost reconciliation (independent recomputation; VSL-only valuation)
if (any(inc$benefit_gdp != 0)) stop("Retired benefit_gdp must be exactly zero.", call. = FALSE)
if (inc[, max(abs(benefit_total - benefit_vsl))] > 0)
  stop("benefit_total must equal benefit_vsl (VSL-only valuation).", call. = FALSE)
vsl_chk <- econ$vsl_mult * econ$gni_pc_base * (1 + econ$gdp_growth)^(inc$year - econ$gni_pc_base_year)
if (max(abs(inc$vsl - vsl_chk) / vsl_chk) > 1e-10)
  stop("VSL path differs from vsl_mult x GNIpc_b x (1 + g)^(t - b).", call. = FALSE)
if (max(abs(inc$benefit_vsl - inc$deaths_avert * vsl_chk)) > 1e-6 * max(1, max(abs(inc$benefit_vsl))))
  stop("benefit_vsl != deaths averted x VSL.", call. = FALSE)
for (loc in c(unique(inc$location), "TOTAL")) {
  d <- if (loc == "TOTAL") inc else inc[location == loc]
  r <- econ_summary[location == loc]
  pvb <- sum((1 + econ$disc_rate)^-(d$year - discount_base_year) * d$deaths_avert * d$vsl)
  pvc <- sum((1 + econ$disc_rate)^-(d$year - discount_base_year) * (d$cost_sap - d$cost_ref))
  if (abs(r$total_benefit - pvb) > 1e-6 * max(1, abs(pvb)) || abs(r$total_cost - pvc) > 1e-6 * max(1, abs(pvc)) ||
      abs(r$net_benefit - (pvb - pvc)) > 1e-6 * max(1, abs(pvb), abs(pvc)) ||
      (pvc > 0 && abs(r$bcr - pvb / pvc) > 1e-10 * max(1, abs(pvb / pvc))))
    stop("Discounted benefits / costs / net benefit / BCR do not reconcile for ", loc, ".", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 6. WRITE + report
# ------------------------------------------------------------------------------
fwrite(budget_impact, paste0(OUT_DIR, "rhd_budget_impact.csv"))
fwrite(econ_summary,  paste0(OUT_DIR, "rhd_economic_summary.csv"))
saveRDS(list(summary = econ_summary, budget_impact = budget_impact,
             incremental = inc, params = econ,
             discount_base_year = discount_base_year,
             uncosted_confirmations = conf_ty),
        paste0(OUT_DIR, "rhd_economic_results.rds"))

t <- summary_total
message("── ECONOMIC RESULTS (discounted; base-year US$) ───────")
message(sprintf("  Total incremental cost : $%.2f billion", t$total_cost   / 1e9))
message(sprintf("  Total monetised benefit: $%.2f billion", t$total_benefit / 1e9))
message(sprintf("  Net benefit            : $%.2f billion", t$net_benefit  / 1e9))
message(sprintf("  Benefit-cost ratio     : %s", if (is.finite(t$bcr)) sprintf("%.2f", t$bcr) else
                paste0("undefined (", t$incremental_cost_status, ")")))
message(sprintf("  RHD deaths averted     : %s", formatC(round(t$deaths_averted), format = "d", big.mark = ",")))
message(sprintf("  Cost per death averted : $%s", formatC(round(t$cost_per_death), format = "d", big.mark = ",")))
message(sprintf("  Cost per DALY averted  : $%s", formatC(round(t$cost_per_daly),  format = "d", big.mark = ",")))
message(sprintf("  NOTE (UNCOSTED): %s incremental confirmatory echoes over %d-%d are NOT priced (no unit cost in 00 block J); total cost is understated by their cost.",
                formatC(round(sum(budget_impact$confirmations_incremental_uncosted)), format = "d", big.mark = ","),
                min(analysis_years), max(analysis_years)))
message("── 08_economic_evaluation.R complete ──────────────────")
message("  Wrote: ", OUT_DIR, "rhd_budget_impact.csv, rhd_economic_summary.csv, rhd_economic_results.rds")
