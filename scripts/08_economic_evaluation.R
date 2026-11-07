# ==============================================================================
# RHD secondary-prevention investment case: ECONOMIC EVALUATION (A/B/C/D)
# scripts/08_economic_evaluation.R
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT (the ONLY place economics live)
# ------------------------------------------------------------------------------
# Benefit-cost + cost-effectiveness of the SAP scale-up vs the reference scenario.
# It consumes 06's per-location outputs (never re-runs the model) and reads every
# economic input from the 00_run_all.R block-J globals (NO hidden defaults, NO
# network calls: missing inputs stop the run with a message naming them).
#
#   INPUTS
#     output/<COUNTRY>/out_model/<location>.rds        (from 06)
#       $stages — service volumes: n_screened, n_on_optimal_treatment, total_surgeries,
#                 n_confirm_completed, by scenario x age x sex x year;
#       $wsd    — RHD deaths (`dead`) and WPP population (`pop`) by scenario x age x sex x year;
#       $meta   — analysis years, country, calibration_adequacy (status carried verbatim).
#     00 block J globals: reporting unit (ECON_CURRENCY, PRICE_YEAR), ECON_MACRO_OBS
#       (WDI/BLS observations), cost_bundles + COST_PRICE_ADJUSTMENT_METHOD, VSL inputs
#       (vsl_*), gdp_growth (real income growth proxy), disc_rate, discount_base_year,
#       dalys_per_death, per-capita denominator, INT$ comparison settings.
#
#   UNITS: every monetary result is CONSTANT 2025 US$ (2025 market exchange rate);
#     discounted results are discounted to discount_base_year. INT$ columns
#     (suffix _int2025) are a final-results conversion for cross-country comparison.
#
# 1. UNIT-COST BUNDLES (local service bundles reported in US$; method "local_cpi_xr")
#      cost_LCU_y0   = cost_USD_y0 x XR_y0          (source-study rate, else y0 official
#                                                    annual average: flagged proxy)
#      cost_LCU_2025 = cost_LCU_y0 x CPI_c[2025] / CPI_c[y0]   (SOURCE country's CPI)
#      cost_USD_2025 = cost_LCU_2025 / XR_2025
#    LCU-originated bundles skip the first step; 2025-US$ bundles are kept unchanged.
#    The US CPI is NEVER used to inflate a local service bundle.
#    COSTS (each 06 volume costed once): c_screen = sum(n_screened) x screen_2025;
#      c_sap = sum(n_on_optimal_treatment) x sap_2025 (on-BPG person-years, not the
#      adherence-weighted volume); c_surgery = sum(total_surgeries) x surgery_2025.
#      UNCOSTED: confirmatory echoes (n_confirm_completed) — no unit cost; flagged.
#
# 2. VSL (reference-case preferred default; Y = vsl_income_anchor_year)
#      pf = CPI_US[Y] / CPI_US[2015];  VSL_US_Y = 9.4e6 pf;  inc_US_Y = 57,900 pf
#      raw_anchor (int$ of Y) = VSL_US_Y x (GNIpc_PPP_Y / inc_US_Y)^e_inc
#      G_t = (1 + g)^(t - Y);  raw_t = raw_anchor x G_t^e_time;  floor_t = 20 x GNIpc_PPP_Y x G_t
#      VSL_t (int$ Y) = max(raw_t, floor_t)  (floor applied every year; unfloored path kept)
#      VSL_t (US$ 2025) = VSL_t x PPP_Y x CPI_c[2025] / CPI_c[Y] / XR_2025
#    benefit_t = deaths averted_t x VSL_t ; benefit_gdp = 0 (retired).
#
# 3. RESULTS. C_t = cost_sap_t - cost_ref_t, B_t = benefit_t, df_t = (1 + r)^-(t - t0)
#      undiscounted: sum(C), sum(B), sum(B - C), bcr = sum(B) / sum(C)
#      discounted  : sum(df C), sum(df B), sum(df (B - C)), bcr = sum(df B) / sum(df C)
#      BCR = NA when its cost denominator <= 0 (status: cost_neutral / cost_saving).
#      total_cost / total_benefit / net_benefit / bcr = aliases of the discounted set.
#      cost per death = sum(df C) / sum(deaths averted); cost per DALY = sum(df C) /
#      (sum(df x deaths averted) x dalys_per_death)  (formulas unchanged).
#    Per capita: cumulative totals / N0 (WPP population, first analysis year, one arm).
#    TOTAL pools location totals and N0 (never a mean of ratios).
#
# OUTPUTS (output/<COUNTRY>/tables/): rhd_economic_summary.csv, rhd_budget_impact.csv,
#   rhd_economic_results.rds (existing names/columns kept; new columns appended),
#   rhd_economic_annual_per_capita.csv, rhd_vsl_derivation.csv,
#   rhd_price_adjustment_audit.csv, rhd_economic_parameter_register.csv;
#   output/rhd_country_comparison_2025.csv (+ .rds) once every country summary exists.
# SENSITIVITY ANALYSIS: PENDING; not run.
# ==============================================================================

library(data.table)

if (!exists("wd_outp")) wd_outp <- paste0(here::here("output"), "/")
IN_DIR  <- paste0(wd_outp, "out_model/")
OUT_DIR <- paste0(wd_outp, "tables/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SENSITIVITY_STATUS <- "PENDING; not run"

# ------------------------------------------------------------------------------
# 0. PARAMETERS — read once from the block-J globals into one named object
# ------------------------------------------------------------------------------
req <- function(nm, allow_null = FALSE) {
  if (!exists(nm, inherits = TRUE))
    stop("Required economic input '", nm, "' is not defined. Evaluate 00_run_all.R block J first ",
         "(08 has no hidden defaults).", call. = FALSE)
  v <- get(nm, inherits = TRUE)
  if (is.null(v) && !allow_null) stop("Required economic input '", nm, "' is NULL (00 block J).", call. = FALSE)
  v
}

par <- list(
  country = req("COUNTRY"), iso3 = req("ISO3"),
  currency = req("ECON_CURRENCY"), price_year = as.integer(req("PRICE_YEAR")),
  cmp_currency = req("ECON_COMPARISON_CURRENCY"), cmp_price_year = as.integer(req("ECON_COMPARISON_PRICE_YEAR")),
  cmp_ppp_series = req("ECON_COMPARISON_PPP_SERIES"),
  cmp_ppp_alt_source = req("ECON_COMPARISON_PPP_ALT_SOURCE", allow_null = TRUE),
  cmp_countries = req("ECON_COMPARISON_COUNTRIES"),
  disc_rate = req("disc_rate"), discount_base_year = as.integer(req("discount_base_year")),
  dalys_per_death = req("dalys_per_death"),
  pc_denominator = req("PER_CAPITA_DENOMINATOR"), pc_denominator_year = as.integer(req("PER_CAPITA_DENOMINATOR_YEAR")),
  macro_obs = req("ECON_MACRO_OBS"), macro_retrieved = req("ECON_MACRO_RETRIEVED"),
  cpi_series = req("ECON_LOCAL_CPI_SERIES"), xr_series = req("ECON_XR_SERIES"), ppp_series = req("ECON_PPP_SERIES"),
  price_method = req("COST_PRICE_ADJUSTMENT_METHOD"), bundles = req("cost_bundles"),
  vsl = list(
    method = req("VSL_METHOD"),
    ref_value = req("vsl_us_ref_value"), ref_price_year = as.integer(req("vsl_us_ref_price_year")),
    ref_income = req("vsl_us_ref_income"), ref_income_year = as.integer(req("vsl_us_ref_income_year")),
    ref_source = req("vsl_us_ref_source"), us_cpi_series = req("vsl_us_cpi_series"),
    anchor_year = as.integer(req("vsl_income_anchor_year")), income_series = req("vsl_income_series"),
    income_elasticity = req("vsl_income_elasticity"), time_elasticity = req("vsl_time_elasticity"),
    floor_mult = req("vsl_floor_mult"), floor_enabled = isTRUE(req("vsl_floor_enabled")),
    income_growth = req("gdp_growth"), income_growth_interpretation = req("income_growth_interpretation")),
  # DEPRECATED globals: carried into $params for report compatibility only (never used here)
  deprecated = list(gni_pc_base = req("gni_pc_base"), gni_pc_base_year = as.integer(req("gni_pc_base_year")),
                    gni_pc_source = req("gni_pc_source"), gdp_pc_base = req("gdp_pc_base"),
                    gdp_pc_base_year = as.integer(req("gdp_pc_base_year")), vsl_mult = req("vsl_mult"))
)

with(par, {
  if (!identical(currency, "USD")) stop("ECON_CURRENCY must be 'USD' (reporting currency).", call. = FALSE)
  if (!identical(price_year, 2025L) || !identical(cmp_price_year, 2025L))
    stop("PRICE_YEAR / ECON_COMPARISON_PRICE_YEAR must be 2025 (columns and labels are constant-2025 units); ",
         "08 never substitutes another reporting year.", call. = FALSE)
  if (!identical(cmp_currency, "INT$")) stop("ECON_COMPARISON_CURRENCY must be 'INT$'.", call. = FALSE)
  if (!identical(price_method, "local_cpi_xr"))
    stop("COST_PRICE_ADJUSTMENT_METHOD '", price_method, "' is not implemented (only 'local_cpi_xr').", call. = FALSE)
  if (!is.finite(disc_rate) || disc_rate < 0 || disc_rate > 0.2) stop("disc_rate outside plausible [0, 0.2].", call. = FALSE)
  if (!is.finite(dalys_per_death) || dalys_per_death <= 0) stop("dalys_per_death must be positive.", call. = FALSE)
  if (!is.finite(vsl$income_growth) || vsl$income_growth <= -1) stop("gdp_growth must be a finite rate > -1.", call. = FALSE)
  if (any(!is.finite(c(vsl$ref_value, vsl$ref_income, vsl$income_elasticity, vsl$time_elasticity, vsl$floor_mult))) ||
      vsl$ref_value <= 0 || vsl$ref_income <= 0 || vsl$floor_mult < 0)
    stop("VSL reference inputs must be finite and positive.", call. = FALSE)
  if (vsl$ref_price_year != vsl$ref_income_year)
    stop("The U.S. reference VSL and reference income must share one reference year.", call. = FALSE)
})

# ---- macro observations: one row per series x country x year (no network) -----
obs_dt <- rbindlist(lapply(names(par$macro_obs), function(k) {
  o <- par$macro_obs[[k]]
  data.table(obs_key = k, series = o$series, iso3 = o$iso3, year = as.integer(o$year),
             value = as.numeric(if (is.null(o$value)) NA_real_ else o$value),
             unit = o$unit, source = o$source, last_updated = o$last_updated, retrieved = o$retrieved,
             note = if (is.null(o$note)) "" else o$note)
}))
macro_value <- function(series_, iso3_, year_, required = TRUE, tbl = obs_dt) {
  r <- tbl[series == series_ & iso3 == iso3_ & year == as.integer(year_)]
  id <- sprintf("%s | %s | %d", series_, iso3_, as.integer(year_))
  if (nrow(r) > 1L) stop("Duplicate macro observation ", id, " in 00 block J ECON_MACRO_OBS.", call. = FALSE)
  if (nrow(r) == 0L || !is.finite(r$value)) {
    msg <- paste0("Macro observation ", id, " is missing or NA in 00 block J ECON_MACRO_OBS (TODO: supply).")
    if (required) stop(msg, call. = FALSE)
    message("  ", msg)
    return(NA_real_)
  }
  if (r$value <= 0) stop("Macro observation ", id, " is not positive.", call. = FALSE)
  r$value
}
if (anyDuplicated(obs_dt[, .(series, iso3, year)]))
  stop("Duplicate series x country x year entries in ECON_MACRO_OBS.", call. = FALSE)

# ------------------------------------------------------------------------------
# 1. UNIT-COST BUNDLES -> constant 2025 US$  (local CPI + exchange rates)
# ------------------------------------------------------------------------------
adjust_bundle <- function(nm, b, T = par$price_year) {
  y0 <- as.integer(b$original_price_year); iso <- b$source_iso3; v0 <- b$original_value
  if (!is.finite(v0) || v0 < 0) stop("Bundle '", nm, "' has an invalid original value.", call. = FALSE)
  if (!b$original_currency %in% c("USD", "LCU")) stop("Bundle '", nm, "': original_currency must be 'USD' or 'LCU'.", call. = FALSE)
  out <- data.table(bundle = nm, param = b$param, unit = b$unit, source_setting = b$source_setting, source_iso3 = iso,
                    source_setting_status = b$source_setting_status, original_value = v0,
                    original_currency = b$original_currency, original_price_year = y0,
                    price_year_status = b$price_year_status, method = par$price_method,
                    xr_y0_lcu_per_usd = NA_real_, xr_y0_source = NA_character_, xr_y0_proxy = NA,
                    cost_lcu_y0 = NA_real_, cpi_series = par$cpi_series, cpi_definition = NA_character_,
                    cpi_y0 = NA_real_, cpi_2025 = NA_real_, cost_lcu_2025 = NA_real_,
                    xr_2025_lcu_per_usd = NA_real_, adjustment_factor = NA_real_, adjusted_value_usd2025 = NA_real_,
                    reporting_currency = par$currency, reporting_price_year = T, path = NA_character_,
                    sources = b$source)
  if (b$original_currency == "USD" && y0 == T) {           # already 2025 US$: never inflated again
    out[, `:=`(adjustment_factor = 1, adjusted_value_usd2025 = v0,
               path = "already in reporting-year US$: retained unchanged")]
    return(out)
  }
  cpi0 <- macro_value(par$cpi_series, iso, y0); cpiT <- macro_value(par$cpi_series, iso, T)
  xrT  <- macro_value(par$xr_series, iso, T)
  if (b$original_currency == "USD") {
    has_src <- is.finite(b$orig_xr)
    xr0 <- if (has_src) b$orig_xr else macro_value(par$xr_series, iso, y0)
    lcu0 <- v0 * xr0
    out[, `:=`(xr_y0_lcu_per_usd = xr0, xr_y0_proxy = !has_src,
               xr_y0_source = if (has_src) "source-study conversion rate" else
                 paste0("PROXY: ", par$xr_series, " ", iso, " ", y0, " official annual average (source-study rate unavailable)"),
               path = "USD -> LCU (y0) -> local CPI to 2025 -> USD at 2025 rate")]
  } else {
    lcu0 <- v0
    out[, path := "LCU (y0) -> local CPI to 2025 -> USD at 2025 rate"]
  }
  lcuT <- lcu0 * cpiT / cpi0
  out[, `:=`(cost_lcu_y0 = lcu0, cpi_y0 = cpi0, cpi_2025 = cpiT, cost_lcu_2025 = lcuT, xr_2025_lcu_per_usd = xrT,
             adjusted_value_usd2025 = lcuT / xrT, adjustment_factor = (lcuT / xrT) / v0,
             cpi_definition = paste0(par$cpi_series, " (", iso, "): ",
                                     obs_dt[series == par$cpi_series & iso3 == iso][1]$unit, ", annual"))]
  out
}
price_audit <- rbindlist(lapply(names(par$bundles), function(nm) adjust_bundle(nm, par$bundles[[nm]])))
uc <- setNames(price_audit$adjusted_value_usd2025, price_audit$param)   # applied unit costs, 2025 US$
for (p in c("cost_screen_per_person", "cost_sap_per_year", "cost_surgery"))
  if (!p %in% names(uc) || !is.finite(uc[[p]])) stop("No adjusted 2025 unit cost for '", p, "'.", call. = FALSE)

# ------------------------------------------------------------------------------
# 2. LOAD model outputs (volumes, deaths, population, calibration status)
# ------------------------------------------------------------------------------
files <- list.files(IN_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0)
  stop("No model outputs in ", IN_DIR, ". 08 never re-runs the model; run 06 first.", call. = FALSE)
outs <- lapply(files, readRDS)
stg <- rbindlist(lapply(outs, function(o) as.data.table(o$stages)), use.names = TRUE, fill = TRUE)
wsd <- rbindlist(lapply(outs, function(o) as.data.table(o$wsd)), use.names = TRUE, fill = TRUE)
if (is.null(stg) || !nrow(stg))
  stop("Model output has no $stages table (incompatible output; re-run 06).", call. = FALSE)
for (cl in c("n_screened", "n_on_optimal_treatment", "total_surgeries", "n_confirm_completed"))
  if (!cl %in% names(stg)) stop("Model $stages lacks '", cl, "' (incompatible output; re-run 06).", call. = FALSE)
for (cl in c("dead", "pop", "age", "sex")) if (!cl %in% names(wsd)) stop("Model $wsd lacks '", cl, "'.", call. = FALSE)
for (o in outs) if (!is.null(o$meta$country) && !identical(as.character(o$meta$country), par$country))
  stop("Model output country '", o$meta$country, "' differs from COUNTRY '", par$country, "'.", call. = FALSE)

analysis_years <- sort(as.integer(outs[[1]]$meta$years))
cal_status <- vapply(outs, function(o) {
  s <- o$meta$calibration_adequacy
  if (is.null(s) || !length(s)) "unknown (not recorded in the model outputs)" else paste(as.character(s), collapse = " ")
}, "")
calibration_status <- paste(unique(cal_status), collapse = " | ")

# ------------------------------------------------------------------------------
# 3. ANNUAL COSTS per location x scenario x year (2025 US$)
# ------------------------------------------------------------------------------
cost_ty <- stg[, .(
  c_screen  = sum(n_screened)             * uc[["cost_screen_per_person"]],
  c_sap     = sum(n_on_optimal_treatment) * uc[["cost_sap_per_year"]],
  c_surgery = sum(total_surgeries)        * uc[["cost_surgery"]]
), by = .(location, scenario, year)]
cost_ty[, cost := c_screen + c_sap + c_surgery]
conf_ty  <- stg[, .(confirmations = sum(n_confirm_completed)), by = .(location, scenario, year)]
death_ty <- wsd[, .(rhd_deaths = sum(dead)), by = .(location, scenario, year)]

# population: one observation per location/age/sex/year/arm; arms must agree (never summed)
if (wsd[, .N, by = .(location, scenario, year, age, sex)][N > 1L, .N] > 0L)
  stop("Duplicate population rows per location/age/sex/year/arm in $wsd.", call. = FALSE)
pop_w <- dcast(wsd[, .(pop = sum(pop)), by = .(location, scenario, year)], location + year ~ scenario, value.var = "pop")
if (!all(c("ref", "sap") %in% names(pop_w)) || pop_w[, max(abs(ref - sap) / pmax(ref, 1))] > 1e-9)
  stop("Reference and scale-up arms do not share one population.", call. = FALSE)
pop_ty <- pop_w[, .(location, year, pop = ref)]

# ------------------------------------------------------------------------------
# 4. VSL path (reference-case default) in 2019 int$ -> constant 2025 US$
# ------------------------------------------------------------------------------
derive_vsl <- function(years, v = par$vsl, iso = par$iso3, T = par$price_year) {
  Y <- v$anchor_year
  cpi_us_ref <- macro_value(v$us_cpi_series, "USA", v$ref_price_year)
  cpi_us_Y   <- macro_value(v$us_cpi_series, "USA", Y)
  pf <- cpi_us_Y / cpi_us_ref                               # US$ of ref year -> US$ / int$ of Y
  vsl_us_Y <- v$ref_value * pf; inc_us_Y <- v$ref_income * pf
  inc_Y <- macro_value(v$income_series, iso, Y)             # int$ of Y per person
  income_ratio <- inc_Y / inc_us_Y
  raw_anchor <- vsl_us_Y * income_ratio^v$income_elasticity # int$ of Y (UNFLOORED; grown below)
  ppp_Y <- macro_value(par$ppp_series, iso, Y); cpi_Y <- macro_value(par$cpi_series, iso, Y)
  cpi_T <- macro_value(par$cpi_series, iso, T); xr_T  <- macro_value(par$xr_series, iso, T)
  to_usd <- ppp_Y * (cpi_T / cpi_Y) / xr_T                  # int$ of Y -> LCU of Y -> LCU of T -> US$ of T
  G <- (1 + v$income_growth)^(years - Y)
  tab <- data.table(year = years, income_growth_index = G, income_target_intY = inc_Y * G,
                    raw_vsl_intY = raw_anchor * G^v$time_elasticity)
  tab[, floor_vsl_intY := if (v$floor_enabled) v$floor_mult * income_target_intY else NA_real_]
  tab[, vsl_floor_binding := v$floor_enabled & floor_vsl_intY > raw_vsl_intY]
  tab[, vsl_intY := if (v$floor_enabled) pmax(raw_vsl_intY, floor_vsl_intY) else raw_vsl_intY]
  tab[, `:=`(vsl_lcu_Y = vsl_intY * ppp_Y)]
  tab[, vsl_lcu_2025 := vsl_lcu_Y * cpi_T / cpi_Y]
  tab[, vsl_usd2025 := vsl_lcu_2025 / xr_T]
  tab[, vsl_unfloored_usd2025 := raw_vsl_intY * to_usd]
  tab[, vsl_to_income_ratio := vsl_intY / income_target_intY]
  anchor <- list(income_anchor_year = Y, us_ref_vsl = v$ref_value, us_ref_income = v$ref_income,
                 us_ref_year = v$ref_price_year, us_cpi_series = v$us_cpi_series,
                 us_cpi_ref_year = cpi_us_ref, us_cpi_anchor_year = cpi_us_Y, reference_price_factor = pf,
                 vsl_us_ref_Y = vsl_us_Y, income_us_ref_Y = inc_us_Y, income_target_Y = inc_Y,
                 income_ratio = income_ratio, raw_vsl_anchor_intY = raw_anchor,
                 floor_vsl_anchor_intY = if (v$floor_enabled) v$floor_mult * inc_Y else NA_real_,
                 ppp_lcu_per_int_Y = ppp_Y, cpi_country_Y = cpi_Y, cpi_country_2025 = cpi_T,
                 xr_lcu_per_usd_2025 = xr_T, intY_to_usd2025 = to_usd,
                 income_elasticity = v$income_elasticity, time_elasticity = v$time_elasticity,
                 floor_mult = v$floor_mult, floor_enabled = v$floor_enabled, income_growth = v$income_growth,
                 income_growth_interpretation = v$income_growth_interpretation,
                 units = c(vsl_us_ref_Y = paste0("US$ of ", Y, " (= int$ of ", Y, ", US is the PPP numeraire)"),
                           income_us_ref_Y = paste0("int$ of ", Y, " per person (published reference observation re-expressed; not a new US income observation)"),
                           income_target_Y = paste0("int$ of ", Y, " per person (", v$income_series, ")"),
                           raw_vsl_anchor_intY = paste0("int$ of ", Y), vsl_lcu_Y = paste0("LCU of ", Y),
                           vsl_lcu_2025 = "LCU of 2025", vsl_usd2025 = "constant 2025 US$ (2025 market exchange rate)"))
  list(annual = tab, anchor = anchor)
}
vsl_d   <- derive_vsl(analysis_years)
vsl_tab <- vsl_d$annual

# ------------------------------------------------------------------------------
# 5. INCREMENTAL (sap - ref) per location x year
# ------------------------------------------------------------------------------
wide_cost  <- dcast(cost_ty,  location + year ~ scenario, value.var = "cost")
wide_death <- dcast(death_ty, location + year ~ scenario, value.var = "rhd_deaths")
if (!all(c("ref", "sap") %in% names(wide_cost)))
  stop("Both 'ref' and 'sap' scenarios must be present in the model output.", call. = FALSE)
inc <- merge(wide_cost[,  .(location, year, cost_ref = ref, cost_sap = sap)],
             wide_death[, .(location, year, death_ref = ref, death_sap = sap)], by = c("location", "year"))
inc[, `:=`(inc_cost = cost_sap - cost_ref, deaths_avert = death_ref - death_sap)]
discount_base_year <- par$discount_base_year
inc[, df := 1 / (1 + par$disc_rate)^(year - discount_base_year)]
inc <- merge(inc, vsl_tab[, .(year, gdp_pc = income_target_intY, vsl = vsl_usd2025,
                              vsl_growth_index = vsl_intY / vsl_d$anchor$raw_vsl_anchor_intY,
                              income_growth_index, vsl_unfloored_usd2025, vsl_floor_binding)], by = "year")
inc[, `:=`(benefit_vsl = deaths_avert * vsl, benefit_gdp = 0)]   # benefit_gdp: RETIRED, identically 0
inc[, benefit_total := benefit_vsl]
inc <- merge(inc, pop_ty, by = c("location", "year"))
inc[, `:=`(currency = par$currency, price_year = par$price_year)]
setcolorder(inc, c("location", "year", "cost_ref", "cost_sap", "death_ref", "death_sap", "inc_cost", "deaths_avert",
                   "df", "gdp_pc", "vsl", "benefit_vsl", "benefit_gdp", "benefit_total", "vsl_growth_index"))
setorder(inc, location, year)

# ------------------------------------------------------------------------------
# 6. SUMMARY per location + pooled TOTAL (USD), per capita, then the INT$ add-on
# ------------------------------------------------------------------------------
cost_status <- function(x) if (!is.finite(x)) "invalid" else if (x > 0) "positive" else if (x == 0) "cost_neutral" else "cost_saving"
safe_ratio  <- function(num, den) if (is.finite(den) && den > 0) num / den else NA_real_
econ_totals <- function(d, r = par$disc_rate, base = discount_base_year) {
  w <- 1 / (1 + r)^(d$year - base); C <- d$cost_sap - d$cost_ref; B <- d$deaths_avert * d$vsl
  list(cost_u = sum(C), ben_u = sum(B), net_u = sum(B - C),
       cost_d = sum(w * C), ben_d = sum(w * B), net_d = sum(w * (B - C)),
       da = sum(d$deaths_avert), dda = sum(w * d$deaths_avert))
}
N0_of <- function(d) d[year == par$pc_denominator_year, sum(pop)]   # inc: one row per location x year (single arm)

summarise_econ <- function(d, label) {
  t <- econ_totals(d); N0 <- N0_of(d)
  data.table(
    location = label,
    total_cost = t$cost_d, total_benefit = t$ben_d, net_benefit = t$net_d,
    bcr = safe_ratio(t$ben_d, t$cost_d),
    deaths_averted = t$da,
    cost_per_death = safe_ratio(t$cost_d, t$da),
    cost_per_daly  = safe_ratio(t$cost_d, t$dda * par$dalys_per_death),
    incremental_cost_status = cost_status(t$cost_d),
    # ---- appended (7 Nov 2026) ----
    currency = par$currency, price_year = par$price_year,
    analysis_start = min(analysis_years), analysis_end = max(analysis_years),
    disc_rate = par$disc_rate, discount_base_year = discount_base_year,
    total_cost_undisc = t$cost_u, total_benefit_undisc = t$ben_u, net_benefit_undisc = t$net_u,
    bcr_undisc = safe_ratio(t$ben_u, t$cost_u), incremental_cost_status_undisc = cost_status(t$cost_u),
    total_cost_disc = t$cost_d, total_benefit_disc = t$ben_d, net_benefit_disc = t$net_d,
    bcr_disc = safe_ratio(t$ben_d, t$cost_d),
    N0 = N0, N0_year = par$pc_denominator_year,
    cost_pc_undisc = t$cost_u / N0, benefit_pc_undisc = t$ben_u / N0, net_benefit_pc_undisc = t$net_u / N0,
    cost_pc_disc = t$cost_d / N0, benefit_pc_disc = t$ben_d / N0, net_benefit_pc_disc = t$net_d / N0,
    vsl_first_year_usd2025 = d[year == min(analysis_years), vsl][1],
    vsl_last_year_usd2025  = d[year == max(analysis_years), vsl][1]
  )
}
if (!par$pc_denominator_year %in% analysis_years)
  stop("PER_CAPITA_DENOMINATOR_YEAR is not an analysis year.", call. = FALSE)
summary_by_loc <- rbindlist(lapply(split(inc, by = "location", keep.by = TRUE), function(d) summarise_econ(d, d$location[1])))
summary_total  <- summarise_econ(inc[location != "TOTAL"], "TOTAL")   # pools totals and N0
econ_summary   <- rbindlist(list(summary_by_loc, summary_total), use.names = TRUE)
usd_snapshot   <- copy(econ_summary)

# ---- INT$ add-on (final results only): INT$_2025 = USD_2025 x XR_2025 / PPP_2025 ----
INT_COLS <- c("total_cost_undisc", "total_benefit_undisc", "net_benefit_undisc",
              "total_cost_disc", "total_benefit_disc", "net_benefit_disc",
              "cost_pc_undisc", "benefit_pc_undisc", "net_benefit_pc_undisc",
              "cost_pc_disc", "benefit_pc_disc", "net_benefit_pc_disc", "cost_per_death", "cost_per_daly")
resolve_int_factor <- function(tbl = obs_dt, iso = par$iso3, T = par$cmp_price_year) {
  xr  <- macro_value(par$xr_series, iso, T, tbl = tbl)                       # same XR as costs and VSL
  ppp <- macro_value(par$cmp_ppp_series, iso, T, required = FALSE, tbl = tbl)
  src <- if (is.finite(ppp)) {
    r <- tbl[series == par$cmp_ppp_series & iso3 == iso & year == T]
    paste0(par$cmp_ppp_series, " ", iso, " ", T, " (", r$source, ", updated ", r$last_updated, "; ", r$note, ") and ",
           par$xr_series, " ", iso, " ", T)
  } else {
    message("  INT$ columns set to NA: 2025 PPP factor ", par$cmp_ppp_series, " | ", iso, " | ", T,
            " missing", if (is.null(par$cmp_ppp_alt_source)) " and no alternative source is declared in block J." else ".")
    NA_character_
  }
  list(xr = xr, ppp = ppp, factor = xr / ppp, source = src)
}
add_int <- function(tab, f) {
  tab <- copy(tab)
  tab[, `:=`(comparison_currency = par$cmp_currency, comparison_price_year = par$cmp_price_year,
             xr_lcu_per_usd_2025 = f$xr, ppp_lcu_per_int_2025 = f$ppp, usd_to_int_2025 = f$factor,
             int_conversion_source = f$source)]
  for (nm in INT_COLS) set(tab, j = paste0(nm, "_int2025"), value = tab[[nm]] * f$factor)
  tab
}
intf <- resolve_int_factor()
econ_summary <- add_int(econ_summary, intf)
econ_summary[, `:=`(calibration_status = calibration_status, sensitivity_status = SENSITIVITY_STATUS)]

# ---- annual per-capita table (present-value contributions per resident of year t) ----
apc_of <- function(d) d[, .(N_t = sum(pop), C_t = sum(inc_cost), B_t = sum(benefit_total), df_t = df[1]), by = year]
annual_pc <- rbindlist(list(
  inc[, apc_of(.SD), by = location],
  apc_of(inc)[, location := "TOTAL"]), use.names = TRUE)
annual_pc[, `:=`(cost_pc_t = C_t / N_t, benefit_pc_t = B_t / N_t,
                 cost_pc_pv_t = df_t * C_t / N_t, benefit_pc_pv_t = df_t * B_t / N_t,
                 currency = par$currency, price_year = par$price_year)]

# ---- budget impact (2025 US$; undiscounted) ----
budget_impact <- inc[, .(location, year, cost_ref, cost_sap, incremental = inc_cost, deaths_avert)]
budget_impact <- merge(budget_impact, cost_ty[scenario == "sap", .(location, year, c_screen, c_sap, c_surgery)],
                       by = c("location", "year"), all.x = TRUE)
conf_w <- dcast(conf_ty, location + year ~ scenario, value.var = "confirmations")
budget_impact <- merge(budget_impact, conf_w[, .(location, year, confirmations_incremental_uncosted = sap - ref)],
                       by = c("location", "year"), all.x = TRUE)
budget_impact[, `:=`(currency = par$currency, price_year = par$price_year)]

# ---- INT$ values in the VSL derivation and price-adjustment audit ----
vsl_tab[, vsl_int2025 := vsl_usd2025 * intf$factor]
price_audit[, adjusted_value_int2025 := adjusted_value_usd2025 * intf$factor]

# ------------------------------------------------------------------------------
# 7. VALIDATION  (arithmetic and data validity; favourable results NOT required)
# ------------------------------------------------------------------------------
chk <- function(ok, msg) if (!isTRUE(ok)) stop("Validation failed: ", msg, call. = FALSE)
rel <- function(a, b) abs(a - b) / pmax(1, abs(b))
validation <- character(0); passed <- function(x) validation <<- c(validation, x)

# (a) service volumes reconcile with 06's diagnostics (each costed once)
for (o in outs) {
  d6 <- o$diag
  if (!is.null(d6$sap_trace_by_year) && !is.null(d6$screening_by_year)) {
    v <- merge(as.data.table(o$stages)[, .(sap_py = sum(n_on_optimal_treatment), enc = sum(n_screened)), by = .(scenario, year)],
               merge(as.data.table(d6$sap_trace_by_year)[, .(scenario, year, on_bpg_person_years)],
                     as.data.table(d6$screening_by_year)[, .(scenario, year, encounters)], by = c("scenario", "year")),
               by = c("scenario", "year"))
    chk(v[, max(abs(sap_py - on_bpg_person_years) / pmax(on_bpg_person_years, 1),
                abs(enc - encounters) / pmax(encounters, 1))] <= 1e-9,
        "costed SAP person-years / screening encounters do not reconcile with 06 diagnostics")
  }
}
passed("service volumes reconcile with 06 diagnostics")
chk(identical(sort(unique(as.integer(inc$year))), analysis_years), "economic years differ from the analysis horizon")
chk(all(cost_ty$cost >= -1e-6), "negative annual cost")
passed("horizon = analysis years; annual costs non-negative")

# (b) price adjustment: independent recomputation from the raw observations; local CPI only
for (i in seq_len(nrow(price_audit))) {
  a <- price_audit[i]; b <- par$bundles[[a$bundle]]
  if (a$original_currency == "USD" && a$original_price_year == par$price_year) {
    chk(a$adjusted_value_usd2025 == b$original_value, paste(a$bundle, "2025-US$ bundle was altered"))
  } else {
    o <- function(s, y) obs_dt[series == s & iso3 == b$source_iso3 & year == y, value]
    xr0 <- if (is.finite(b$orig_xr)) b$orig_xr else o("PA.NUS.FCRF", b$original_price_year)
    f_ind <- (if (b$original_currency == "USD") xr0 else 1) * (o("FP.CPI.TOTL", par$price_year) / o("FP.CPI.TOTL", b$original_price_year)) /
      o("PA.NUS.FCRF", par$price_year)
    chk(rel(a$adjusted_value_usd2025, b$original_value * f_ind) < 1e-12, paste(a$bundle, "price adjustment does not reconcile"))
    chk(a$cpi_series != par$vsl$us_cpi_series && !grepl("USA", a$cpi_definition, fixed = TRUE),
        paste(a$bundle, "local bundle inflated with the US CPI"))
  }
  chk(uc[[a$param]] == a$adjusted_value_usd2025, paste(a$bundle, "applied unit cost differs from the audit"))
}
# synthetic path checks: a 2025-US$ bundle is unchanged; an LCU bundle skips the USD->LCU step
tb <- par$bundles[[1]]; tb$original_price_year <- par$price_year
chk(adjust_bundle("selftest_usd2025", tb)$adjusted_value_usd2025 == tb$original_value, "2025-US$ self-test")
tb <- par$bundles[[1]]; tb$original_currency <- "LCU"; tb$original_value <- 1000
chk(rel(adjust_bundle("selftest_lcu", tb)$adjusted_value_usd2025,
        1000 * macro_value(par$cpi_series, tb$source_iso3, par$price_year) /
          macro_value(par$cpi_series, tb$source_iso3, tb$original_price_year) /
          macro_value(par$xr_series, tb$source_iso3, par$price_year)) < 1e-12, "LCU self-test")
passed("bundle conversions reconcile; 2025-US$ unchanged; LCU path; never the US CPI for local bundles")

# (c) VSL: reference-price alignment, unfloored growth, yearly floor, PPP->LCU->USD, positivity
A <- vsl_d$anchor
chk(rel(A$vsl_us_ref_Y / A$income_us_ref_Y, par$vsl$ref_value / par$vsl$ref_income) < 1e-12, "reference price alignment changed the VSL/income ratio")
chk(rel(A$reference_price_factor, macro_value(par$vsl$us_cpi_series, "USA", A$income_anchor_year) /
          macro_value(par$vsl$us_cpi_series, "USA", par$vsl$ref_price_year)) < 1e-12, "reference price factor")
chk(max(rel(vsl_tab$raw_vsl_intY, A$raw_vsl_anchor_intY * vsl_tab$income_growth_index^par$vsl$time_elasticity)) < 1e-12,
    "raw VSL path is not the unfloored anchor grown at the time elasticity")
if (par$vsl$floor_enabled)
  chk(max(rel(vsl_tab$vsl_intY, pmax(vsl_tab$raw_vsl_intY, par$vsl$floor_mult * vsl_tab$income_target_intY))) < 1e-12 &&
        all(vsl_tab$vsl_intY >= vsl_tab$floor_vsl_intY * (1 - 1e-12)), "floor not applied in every year")
chk(max(rel(vsl_tab$vsl_usd2025, vsl_tab$vsl_intY * A$ppp_lcu_per_int_Y * A$cpi_country_2025 / A$cpi_country_Y /
              A$xr_lcu_per_usd_2025)) < 1e-12, "PPP -> LCU -> 2025 US$ conversion")
chk(all(is.finite(vsl_tab$vsl_usd2025) & vsl_tab$vsl_usd2025 > 0), "VSL not positive and finite")
chk(max(rel(inc$benefit_total, inc$deaths_avert * inc$vsl)) < 1e-12 && all(inc$benefit_gdp == 0),
    "benefit != deaths averted x VSL, or benefit_gdp != 0")
passed("VSL reference alignment, unfloored growth, yearly floor, PPP/XR conversion, benefits = deaths averted x VSL")

# (d) discounted / undiscounted reconciliation (independent loop), and rate 0 => equal
for (loc in econ_summary$location) {
  d <- if (loc == "TOTAL") inc else inc[location == loc]; r <- econ_summary[location == loc]
  C <- 0; B <- 0; Cd <- 0; Bd <- 0
  for (k in seq_len(nrow(d))) {
    w <- (1 + par$disc_rate)^-(d$year[k] - discount_base_year)
    c1 <- d$cost_sap[k] - d$cost_ref[k]; b1 <- (d$death_ref[k] - d$death_sap[k]) * d$vsl[k]
    C <- C + c1; B <- B + b1; Cd <- Cd + w * c1; Bd <- Bd + w * b1
  }
  near <- function(a, b, scale) abs(a - b) <= 1e-9 * max(1, scale)
  chk(near(r$total_cost_undisc, C, abs(C)) && near(r$total_benefit_undisc, B, abs(B)) &&
        near(r$net_benefit_undisc, B - C, max(abs(B), abs(C))) &&
        near(r$total_cost_disc, Cd, abs(Cd)) && near(r$total_benefit_disc, Bd, abs(Bd)) &&
        near(r$net_benefit_disc, Bd - Cd, max(abs(Bd), abs(Cd))),
      paste(loc, "totals do not reconcile"))
  chk(identical(r$bcr_disc, if (Cd > 0) r$total_benefit_disc / r$total_cost_disc else NA_real_) &&
        identical(r$bcr_undisc, if (C > 0) r$total_benefit_undisc / r$total_cost_undisc else NA_real_),
      paste(loc, "BCR definition (NA when cost <= 0)"))
  chk(identical(c(r$total_cost, r$total_benefit, r$net_benefit, r$bcr),
                c(r$total_cost_disc, r$total_benefit_disc, r$net_benefit_disc, r$bcr_disc)),
      paste(loc, "headline columns are not aliases of the discounted results"))
  t0 <- econ_totals(d, r = 0)
  chk(t0$cost_d == t0$cost_u && t0$ben_d == t0$ben_u && t0$net_d == t0$net_u, paste(loc, "rate-0 discounting differs from undiscounted"))
  # per-capita x N0 reproduces totals
  chk(max(rel(unlist(r[, .(cost_pc_undisc, benefit_pc_undisc, net_benefit_pc_undisc, cost_pc_disc, benefit_pc_disc, net_benefit_pc_disc)]) * r$N0,
              unlist(r[, .(total_cost_undisc, total_benefit_undisc, net_benefit_undisc, total_cost_disc, total_benefit_disc, net_benefit_disc)]))) < 1e-9,
      paste(loc, "per-capita x N0 does not reproduce totals"))
}
chk(econ_summary[location == "TOTAL", N0] == sum(econ_summary[location != "TOTAL", N0]) &&
      econ_summary[location == "TOTAL", N0] == inc[year == par$pc_denominator_year, sum(pop)],
    "TOTAL N0 is not the pooled single-arm population")
chk(rel(annual_pc[location == "TOTAL", sum(C_t)], inc[, sum(inc_cost)]) < 1e-12 &&
      rel(annual_pc[location == "TOTAL", sum(df_t * B_t)], econ_summary[location == "TOTAL", total_benefit_disc]) < 1e-9 &&
      annual_pc[location == "TOTAL" & year == par$pc_denominator_year, N_t] == econ_summary[location == "TOTAL", N0],
    "annual per-capita table does not reconcile")
passed("discounted/undiscounted totals, BCR, rate-0 equality, per-capita x N0, pooled TOTAL")

# (e) INT$ add-on: exact scaling, BCR invariance per country, USD columns untouched, missing-PPP path
same_usd <- function(tab) identical(as.list(tab[, names(usd_snapshot), with = FALSE]), as.list(usd_snapshot))
chk(same_usd(econ_summary), "a USD column changed when the INT$ columns were added")
if (is.finite(intf$factor)) {
  for (nm in INT_COLS) {
    a <- econ_summary[[paste0(nm, "_int2025")]]; b <- econ_summary[[nm]] * intf$xr / intf$ppp
    chk(all((is.na(a) & is.na(b)) | rel(a, b) < 1e-12), paste("INT$ column", nm))
  }
  for (i in seq_len(nrow(econ_summary))) {
    r <- econ_summary[i]
    for (s in c("disc", "undisc")) {
      bu <- r[[paste0("bcr_", s)]]; ci <- r[[paste0("total_cost_", s, "_int2025")]]; bi <- r[[paste0("total_benefit_", s, "_int2025")]]
      chk((is.na(bu) && !(ci > 0)) || abs(bi / ci - bu) <= 1e-12 * max(1, abs(bu)), paste(r$location, s, "BCR differs between USD and INT$"))
    }
  }
}
# missing 2025 PPP factor: NA INT$ columns with a named message; USD results untouched
obs_wo <- obs_dt[!(series == par$cmp_ppp_series & iso3 == par$iso3 & year == par$cmp_price_year)]
msgs <- character(0)
f_na <- withCallingHandlers(resolve_int_factor(tbl = obs_wo),
                            message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
t_na <- add_int(usd_snapshot, f_na)
chk(is.na(f_na$factor) && all(vapply(paste0(INT_COLS, "_int2025"), function(j) all(is.na(t_na[[j]])), NA)) &&
      any(grepl(paste(par$cmp_ppp_series, "|", par$iso3, "|", par$cmp_price_year), msgs, fixed = TRUE)) &&
      same_usd(t_na),
    "missing-PPP path does not yield NA INT$ columns with a named message")
passed("INT$ = USD x XR/PPP; per-country BCR identical in USD and INT$; USD unchanged; missing PPP -> NA + message")

# (f) informational notes (no favourable-result requirement)
neg_da <- inc[deaths_avert < -1e-6, .N]; neg_ic <- inc[inc_cost < -1e-6, .N]
if (neg_da > 0) message(sprintf("  NOTE: %d location-year(s) with negative deaths averted.", neg_da))
if (neg_ic > 0) message(sprintf("  NOTE: %d location-year(s) with negative incremental cost.", neg_ic))
if (sum(inc$deaths_avert) < 0) message("  NOTE: cumulative deaths averted is negative (scale-up worse than reference overall).")
for (i in seq_len(nrow(econ_summary))) if (econ_summary$total_cost[i] <= 0)
  message(sprintf("  NOTE: %s discounted incremental cost is not positive -> BCR undefined (%s).",
                  econ_summary$location[i], econ_summary$incremental_cost_status[i]))
vsl_growth_rate <- (1 + par$vsl$income_growth)^par$vsl$time_elasticity - 1
vsl_vs_disc_note <- sprintf("Unfloored VSL grows at %.2f%% per year (= (1 + %.3f)^%.1f - 1) against a %.1f%% discount rate: discounted VSL %s over time.",
                            100 * vsl_growth_rate, par$vsl$income_growth, par$vsl$time_elasticity, 100 * par$disc_rate,
                            if (vsl_growth_rate > par$disc_rate) "rises" else "does not rise")

# ------------------------------------------------------------------------------
# 8. PARAMETER REGISTER, LABELS, DEFINITIONS
# ------------------------------------------------------------------------------
reg <- function(name, value, unit, role, source, year = NA_integer_, sens = FALSE)
  data.table(name = name, value = paste(format(value, digits = 15), collapse = "; "),
             value_numeric = if (is.numeric(value) && length(value) == 1L) as.numeric(value) else NA_real_,
             unit = unit, role = role, source = source, year = as.integer(year), future_sensitivity = sens)
param_register <- rbindlist(c(list(
  reg("ECON_CURRENCY", par$currency, "-", "reporting currency", "00 block J"),
  reg("PRICE_YEAR", par$price_year, "year", "reporting price year (constant prices)", "00 block J", par$price_year),
  reg("disc_rate", par$disc_rate, "per year", "discount rate, costs and benefits", "00 block J [PAPER]", NA, TRUE),
  reg("discount_base_year", discount_base_year, "year", "discount reference year", "00 block J", discount_base_year, TRUE),
  reg("dalys_per_death", par$dalys_per_death, "DALYs per death (undiscounted)", "cost per DALY (unchanged assumption)", "00 block J [CALIBRATE]"),
  reg("PER_CAPITA_DENOMINATOR", par$pc_denominator, "persons", "per-capita denominator convention", "00 block J", par$pc_denominator_year),
  reg("N0", econ_summary[location == "TOTAL", N0], "persons", "per-capita denominator value", "06 $wsd pop (WPP2024), ref arm", par$pc_denominator_year),
  reg("COST_PRICE_ADJUSTMENT_METHOD", par$price_method, "-", "aggregate local-service price adjustment", "Turner et al. 2019"),
  reg("VSL_METHOD", par$vsl$method, "-", "valuation method", par$vsl$ref_source),
  reg("vsl_us_ref_value", par$vsl$ref_value, paste0("US$ of ", par$vsl$ref_price_year), "U.S. reference VSL", par$vsl$ref_source, par$vsl$ref_price_year, TRUE),
  reg("vsl_us_ref_income", par$vsl$ref_income, paste0("int$ of ", par$vsl$ref_income_year, " per person"), "U.S. reference income", par$vsl$ref_source, par$vsl$ref_income_year, TRUE),
  reg("vsl_us_ref_price_year", par$vsl$ref_price_year, "year", "reference price year", par$vsl$ref_source, par$vsl$ref_price_year, TRUE),
  reg("vsl_income_anchor_year", par$vsl$anchor_year, "year", "target-country income anchor year Y", "00 block J", par$vsl$anchor_year, TRUE),
  reg("vsl_income_elasticity", par$vsl$income_elasticity, "-", "cross-country income elasticity", par$vsl$ref_source, NA, TRUE),
  reg("vsl_time_elasticity", par$vsl$time_elasticity, "-", "elasticity on real income growth over time", par$vsl$ref_source, NA, TRUE),
  reg("vsl_floor_mult", par$vsl$floor_mult, "x GNI per capita", "VSL floor (every year)", par$vsl$ref_source, NA, TRUE),
  reg("vsl_floor_enabled", par$vsl$floor_enabled, "-", "floor switch", "00 block J", NA, TRUE),
  reg("gdp_growth", par$vsl$income_growth, "per year (real)", par$vsl$income_growth_interpretation, "00 block J [ASSUMPTION]", NA, TRUE),
  reg("reference_price_factor", A$reference_price_factor, "CPI-U ratio", "US CPI-U Y / ref year (US reference values only)", "BLS CPI-U", A$income_anchor_year),
  reg("raw_vsl_anchor_intY", A$raw_vsl_anchor_intY, paste0("int$ of ", A$income_anchor_year), "unfloored VSL at the anchor (derived)", "08"),
  reg("vsl_first_year_usd2025", vsl_tab[1, vsl_usd2025], "constant 2025 US$", "VSL in the first analysis year (derived)", "08", min(analysis_years)),
  reg("vsl_last_year_usd2025", vsl_tab[.N, vsl_usd2025], "constant 2025 US$", "VSL in the last analysis year (derived)", "08", max(analysis_years)),
  reg("ECON_COMPARISON_CURRENCY", par$cmp_currency, "-", "secondary comparison unit (final results only)", "00 block J"),
  reg("ECON_COMPARISON_PPP_SERIES", par$cmp_ppp_series, "-", "2025 PPP conversion factor series", "00 block J", par$cmp_price_year),
  reg("ECON_COMPARISON_PPP_ALT_SOURCE", if (is.null(par$cmp_ppp_alt_source)) "none" else par$cmp_ppp_alt_source, "-", "fallback if 2025 PPP missing", "00 block J"),
  reg("usd_to_int_2025", intf$factor, "INT$ per US$", "XR_2025 / PPP_2025 (derived)", if (is.na(intf$source)) "missing" else intf$source, par$cmp_price_year),
  reg("sensitivity_analysis", SENSITIVITY_STATUS, "-", "status", "08", NA, FALSE)),
  lapply(seq_len(nrow(price_audit)), function(i) { a <- price_audit[i]; rbind(
    reg(paste0(a$param, "_original"), a$original_value, paste(a$original_currency, a$unit), paste0("original bundle value (", a$source_setting, ")"), a$sources, a$original_price_year, TRUE),
    reg(paste0(a$param, "_original_price_year"), a$original_price_year, "year", a$price_year_status, a$sources, a$original_price_year, TRUE),
    reg(paste0(a$param, "_usd2025"), a$adjusted_value_usd2025, paste("constant 2025 US$", a$unit), "applied unit cost (derived)", a$method, par$price_year)) }),
  lapply(seq_len(nrow(obs_dt)), function(i) { o <- obs_dt[i]
    reg(o$obs_key, o$value, o$unit, "macro observation", paste0(o$source, "; updated ", o$last_updated, "; retrieved ", o$retrieved,
                                                              if (nzchar(o$note)) paste0("; ", o$note) else ""), o$year) }),
  lapply(names(par$deprecated), function(nm)
    reg(nm, par$deprecated[[nm]], "-", "DEPRECATED: retained for report compatibility; not used by 08", "00 block J")))
)

report_labels <- list(
  currency = par$currency, price_year = par$price_year, currency_label = "constant 2025 US$",
  comparison_currency = par$cmp_currency, comparison_currency_label = "constant 2025 international $",
  discounting_label = sprintf("discounted to %d at %g%% per year", discount_base_year, 100 * par$disc_rate),
  vsl_method_label = sprintf("Reference-case VSL: PPP GNI, elasticity %g, floor %g x GNI per capita",
                             par$vsl$income_elasticity, par$vsl$floor_mult),
  analysis_period_label = sprintf("%d-%d", min(analysis_years), max(analysis_years)),
  per_capita_label = sprintf("cumulative, per resident at programme start (WPP population, all ages, both sexes, %d)", par$pc_denominator_year),
  income_growth_label = par$vsl$income_growth_interpretation,
  price_adjustment_label = "Unit-cost bundles (assumed 2019 US$) converted to 2025 US$ with local CPI and official exchange rates",
  sensitivity_status = SENSITIVITY_STATUS, calibration_status = calibration_status,
  vsl_vs_discount_note = vsl_vs_disc_note)

D <- function(formula, unit) list(formula = formula, unit = unit)
definitions <- list(
  summary = list(
    total_cost = D("alias of total_cost_disc", "constant 2025 US$, discounted"),
    total_benefit = D("alias of total_benefit_disc", "constant 2025 US$, discounted"),
    net_benefit = D("alias of net_benefit_disc", "constant 2025 US$, discounted"),
    bcr = D("alias of bcr_disc", "ratio"),
    cost_per_death = D("sum(df_t C_t) / sum(deaths averted_t) (undiscounted deaths)", "constant 2025 US$ per death averted"),
    cost_per_daly = D("sum(df_t C_t) / (sum(df_t deaths averted_t) x dalys_per_death)", "constant 2025 US$ per DALY averted"),
    currency = D("reporting currency", "-"), price_year = D("constant-price year", "year"),
    analysis_start = D("first analysis year", "year"), analysis_end = D("last analysis year", "year"),
    disc_rate = D("annual discount rate r", "per year"), discount_base_year = D("t0 in df_t = (1 + r)^-(t - t0)", "year"),
    total_cost_undisc = D("sum(C_t), C_t = cost_sap_t - cost_ref_t", "constant 2025 US$"),
    total_benefit_undisc = D("sum(B_t), B_t = deaths averted_t x VSL_t", "constant 2025 US$"),
    net_benefit_undisc = D("sum(B_t - C_t)", "constant 2025 US$"),
    bcr_undisc = D("sum(B_t) / sum(C_t); NA if sum(C_t) <= 0", "ratio"),
    incremental_cost_status_undisc = D("sign of sum(C_t): positive / cost_neutral / cost_saving", "-"),
    incremental_cost_status = D("sign of sum(df_t C_t)", "-"),
    total_cost_disc = D("sum(df_t C_t)", "constant 2025 US$, discounted to t0"),
    total_benefit_disc = D("sum(df_t B_t)", "constant 2025 US$, discounted to t0"),
    net_benefit_disc = D("sum(df_t (B_t - C_t))", "constant 2025 US$, discounted to t0"),
    bcr_disc = D("sum(df_t B_t) / sum(df_t C_t); NA if denominator <= 0", "ratio"),
    N0 = D("WPP population, all ages, both sexes, one arm, in N0_year (TOTAL: sum over locations)", "persons"),
    N0_year = D("first analysis year", "year"),
    cost_pc_undisc = D("total_cost_undisc / N0 (incremental, not gross spending)", "constant 2025 US$ per resident at programme start"),
    benefit_pc_undisc = D("total_benefit_undisc / N0", "constant 2025 US$ per resident at programme start"),
    net_benefit_pc_undisc = D("net_benefit_undisc / N0", "constant 2025 US$ per resident at programme start"),
    cost_pc_disc = D("total_cost_disc / N0 (same undiscounted N0)", "constant 2025 US$ (PV) per resident at programme start"),
    benefit_pc_disc = D("total_benefit_disc / N0", "constant 2025 US$ (PV) per resident at programme start"),
    net_benefit_pc_disc = D("net_benefit_disc / N0", "constant 2025 US$ (PV) per resident at programme start"),
    vsl_first_year_usd2025 = D("VSL_t at the first analysis year", "constant 2025 US$"),
    vsl_last_year_usd2025 = D("VSL_t at the last analysis year", "constant 2025 US$"),
    comparison_currency = D("secondary unit", "-"), comparison_price_year = D("INT$ price year", "year"),
    xr_lcu_per_usd_2025 = D("PA.NUS.FCRF 2025 (same as costs and VSL)", "LCU per US$"),
    ppp_lcu_per_int_2025 = D("PA.NUS.PPP 2025", "LCU per international $"),
    usd_to_int_2025 = D("xr_lcu_per_usd_2025 / ppp_lcu_per_int_2025", "INT$ per US$"),
    int_conversion_source = D("sources of the XR and PPP observations", "-"),
    `*_int2025` = D("matching USD column x usd_to_int_2025 (final-results conversion; per-country BCR unchanged)",
                    "constant 2025 international $ (purchasing-power comparison; not budget figures)"),
    calibration_status = D("model meta$calibration_adequacy, verbatim", "-"),
    sensitivity_status = D("sensitivity analysis status", "-")),
  incremental = list(
    df = D("1 / (1 + r)^(year - t0)", "-"),
    gdp_pc = D("retained column name: projected real GNI per capita, PPP = GNIpc_PPP_Y x (1 + g)^(t - Y) (g: GDP-pc growth proxy)", "int$ of the anchor year Y per person"),
    vsl = D("VSL_t = max(raw_t, floor_t) x PPP_Y x CPI_c[2025] / CPI_c[Y] / XR_2025", "constant 2025 US$"),
    vsl_growth_index = D("VSL_t (int$ Y) / unfloored anchor VSL", "index (anchor year = 1 if unfloored)"),
    income_growth_index = D("(1 + g)^(t - Y)", "index"),
    vsl_unfloored_usd2025 = D("raw_t converted as vsl", "constant 2025 US$"),
    vsl_floor_binding = D("floor_t > raw_t", "logical"),
    benefit_vsl = D("deaths_avert x vsl", "constant 2025 US$, undiscounted"),
    benefit_gdp = D("retired; identically 0", "-"),
    pop = D("WPP population, all ages, both sexes (ref arm = sap arm)", "persons")),
  annual_per_capita = list(
    N_t = D("WPP population in year t (one arm)", "persons"), C_t = D("incremental cost in t", "constant 2025 US$"),
    B_t = D("monetised benefit in t", "constant 2025 US$"), df_t = D("1 / (1 + r)^(t - t0)", "-"),
    cost_pc_t = D("C_t / N_t", "constant 2025 US$ per resident of year t"),
    benefit_pc_t = D("B_t / N_t", "constant 2025 US$ per resident of year t"),
    cost_pc_pv_t = D("df_t C_t / N_t: present-value contribution per resident of year t", "constant 2025 US$ (PV)"),
    benefit_pc_pv_t = D("df_t B_t / N_t: present-value contribution per resident of year t", "constant 2025 US$ (PV)")),
  vsl_derivation = c(lapply(vsl_d$anchor$units, function(u) D("see 08 section 4", u)),
    list(raw_vsl_intY = D("raw_anchor x income_growth_index^time_elasticity", "int$ of Y"),
         floor_vsl_intY = D("floor_mult x income_target_intY", "int$ of Y"),
         vsl_intY = D("max(raw_vsl_intY, floor_vsl_intY)", "int$ of Y"),
         vsl_int2025 = D("vsl_usd2025 x usd_to_int_2025", "constant 2025 international $"))),
  price_adjustment = list(
    adjusted_value_usd2025 = D("original_value x XR_y0 x CPI_c[2025] / CPI_c[y0] / XR_2025 (local_cpi_xr)", "constant 2025 US$"),
    adjustment_factor = D("adjusted_value_usd2025 / original_value", "-"),
    adjusted_value_int2025 = D("adjusted_value_usd2025 x usd_to_int_2025", "constant 2025 international $")),
  country_comparison = D("one row per country read from each saved rhd_economic_summary.csv; POOLED row computed separately in each unit from that unit's pooled totals (pooled BCR / per-capita differ between USD and INT$). INT$ values are purchasing-power-equivalent magnitudes for comparison, not budget figures; budget impact remains in 2025 US$.", "mixed"))

# ------------------------------------------------------------------------------
# 9. WRITE
# ------------------------------------------------------------------------------
# $params keeps every historical element name (report compatibility). cost_* are the
# APPLIED 2025-US$ unit costs (volumes x params reproduce cost_ref / cost_sap);
# gni_pc_base, gni_pc_base_year, gni_pc_source, gdp_pc_base, gdp_pc_base_year and
# vsl_mult are DEPRECATED pass-throughs, not used in any calculation.
params_compat <- list(
  disc_rate = par$disc_rate,
  gni_pc_base = par$deprecated$gni_pc_base, gni_pc_base_year = par$deprecated$gni_pc_base_year,
  gni_pc_source = par$deprecated$gni_pc_source,
  gdp_pc_base = par$deprecated$gdp_pc_base, gdp_pc_base_year = par$deprecated$gdp_pc_base_year,
  gdp_growth = par$vsl$income_growth, vsl_mult = par$deprecated$vsl_mult, dalys_per_death = par$dalys_per_death,
  cost_screen_per_person = uc[["cost_screen_per_person"]], cost_sap_per_year = uc[["cost_sap_per_year"]],
  cost_surgery = uc[["cost_surgery"]],
  valuation_method = paste0(par$vsl$method, "; VSL_t in constant 2025 US$ via PPP_Y, local CPI and XR_2025; no additive GDP-per-capita benefit"),
  uncosted_items = "confirmatory echo after a true screen positive (n_confirm_completed): no unit cost in 00 block J; omitted from total costs",
  # appended
  currency = par$currency, price_year = par$price_year,
  cost_screen_per_person_original = par$bundles$screen$original_value,
  cost_sap_per_year_original = par$bundles$sap$original_value,
  cost_surgery_original = par$bundles$surgery$original_value,
  deprecated_fields = c("gni_pc_base", "gni_pc_base_year", "gni_pc_source", "gdp_pc_base", "gdp_pc_base_year", "vsl_mult"))

fwrite(budget_impact, paste0(OUT_DIR, "rhd_budget_impact.csv"))
fwrite(econ_summary,  paste0(OUT_DIR, "rhd_economic_summary.csv"))
fwrite(annual_pc,     paste0(OUT_DIR, "rhd_economic_annual_per_capita.csv"))
vsl_out <- copy(vsl_tab)   # annual path + the anchor constants repeated on every row
for (nm in c("income_anchor_year", "us_ref_vsl", "us_ref_income", "us_ref_year", "us_cpi_ref_year", "us_cpi_anchor_year",
             "reference_price_factor", "vsl_us_ref_Y", "income_us_ref_Y", "income_target_Y", "income_ratio",
             "raw_vsl_anchor_intY", "floor_vsl_anchor_intY", "ppp_lcu_per_int_Y", "cpi_country_Y", "cpi_country_2025",
             "xr_lcu_per_usd_2025", "intY_to_usd2025", "income_elasticity", "time_elasticity", "floor_mult", "income_growth"))
  set(vsl_out, j = nm, value = vsl_d$anchor[[nm]])
vsl_out[, `:=`(usd_to_int_2025 = intf$factor, iso3 = par$iso3, currency = par$currency, price_year = par$price_year)]
fwrite(vsl_out, paste0(OUT_DIR, "rhd_vsl_derivation.csv"))
fwrite(price_audit,    paste0(OUT_DIR, "rhd_price_adjustment_audit.csv"))
fwrite(param_register, paste0(OUT_DIR, "rhd_economic_parameter_register.csv"))

results <- list(summary = econ_summary, budget_impact = budget_impact, incremental = inc, params = params_compat,
                discount_base_year = discount_base_year, uncosted_confirmations = conf_ty,
                # appended (7 Nov 2026)
                report_labels = report_labels, definitions = definitions,
                parameters = par[setdiff(names(par), "macro_obs")], macro_observations = obs_dt,
                parameter_register = param_register, price_adjustment = price_audit,
                vsl_derivation = list(anchor = vsl_d$anchor, annual = vsl_tab),
                annual_per_capita = annual_pc,
                int_conversion = list(currency = par$cmp_currency, price_year = par$cmp_price_year, xr_lcu_per_usd_2025 = intf$xr,
                                      ppp_lcu_per_int_2025 = intf$ppp, usd_to_int_2025 = intf$factor, source = intf$source),
                calibration_status = calibration_status, sensitivity_status = SENSITIVITY_STATUS,
                validation = validation)
saveRDS(results, paste0(OUT_DIR, "rhd_economic_results.rds"))

# ------------------------------------------------------------------------------
# 10. COUNTRY COMPARISON (reads the saved summaries; needs every country)
# ------------------------------------------------------------------------------
root_out <- normalizePath(file.path(wd_outp, ".."), winslash = "/")
cmp_files <- setNames(file.path(root_out, par$cmp_countries, "tables", "rhd_economic_summary.csv"), par$cmp_countries)
need_cols <- c("total_cost_undisc", "total_cost_disc", "usd_to_int_2025", "N0", "calibration_status")
cmp_rows <- lapply(names(cmp_files), function(ct) {
  f <- cmp_files[[ct]]
  if (!file.exists(f)) return(paste0(ct, ": summary missing (", f, ")"))
  s <- fread(f)
  if (!all(need_cols %in% names(s))) return(paste0(ct, ": summary predates the 2025-US$ update (re-run 08 for ", ct, ")"))
  s[location == ct]
})
miss <- vapply(cmp_rows, is.character, NA)
if (any(miss)) {
  message("  Country comparison NOT written: ", paste(unlist(cmp_rows[miss]), collapse = "; "))
} else {
  cmp <- rbindlist(cmp_rows, use.names = TRUE, fill = TRUE)
  tot <- c("total_cost_undisc", "total_benefit_undisc", "net_benefit_undisc", "total_cost_disc", "total_benefit_disc", "net_benefit_disc")
  pc  <- c("cost_pc_undisc", "benefit_pc_undisc", "net_benefit_pc_undisc", "cost_pc_disc", "benefit_pc_disc", "net_benefit_pc_disc")
  keep <- c("location", "N0", "N0_year", tot, paste0(tot, "_int2025"), pc, paste0(pc, "_int2025"),
            "bcr_undisc", "bcr_disc", "cost_per_death", "cost_per_daly", "cost_per_death_int2025", "cost_per_daly_int2025",
            "currency", "price_year", "disc_rate", "discount_base_year", "analysis_start", "analysis_end",
            "comparison_currency", "xr_lcu_per_usd_2025", "ppp_lcu_per_int_2025", "usd_to_int_2025", "int_conversion_source",
            "calibration_status")
  cmp <- cmp[, intersect(keep, names(cmp)), with = FALSE]
  setnames(cmp, "location", "country")
  cmp[, `:=`(row_type = "country", bcr_undisc_int2025 = NA_real_, bcr_disc_int2025 = NA_real_)]
  # POOLED row: each unit pooled separately from that unit's own totals (never mixed)
  pool <- data.table(country = paste(cmp$country, collapse = " + "), row_type = "pooled (computed separately in each unit)",
                     N0 = sum(cmp$N0), N0_year = cmp$N0_year[1])
  for (nm in c(tot, paste0(tot, "_int2025"))) pool[, (nm) := sum(cmp[[nm]])]
  for (i in seq_along(pc)) {
    pool[, (pc[i]) := pool[[tot[i]]] / N0]
    pool[, (paste0(pc[i], "_int2025")) := pool[[paste0(tot[i], "_int2025")]] / N0]
  }
  pool[, `:=`(bcr_undisc = safe_ratio(total_benefit_undisc, total_cost_undisc), bcr_disc = safe_ratio(total_benefit_disc, total_cost_disc),
              bcr_undisc_int2025 = safe_ratio(total_benefit_undisc_int2025, total_cost_undisc_int2025),
              bcr_disc_int2025 = safe_ratio(total_benefit_disc_int2025, total_cost_disc_int2025),
              currency = cmp$currency[1], price_year = cmp$price_year[1], comparison_currency = cmp$comparison_currency[1])]
  cmp <- rbindlist(list(cmp, pool), use.names = TRUE, fill = TRUE)
  setcolorder(cmp, c("country", "row_type"))
  fwrite(cmp, file.path(root_out, "rhd_country_comparison_2025.csv"))
  saveRDS(list(table = cmp, definitions = definitions$country_comparison,
               note = paste("Per-country BCRs are identical in USD and INT$ (bcr_*_int2025 is filled only on the pooled row,",
                            "where the two units differ). INT$ values are purchasing-power-equivalent magnitudes for comparison,",
                            "not budget figures; budget impact remains in constant 2025 US$."),
               sources = cmp_files),
          file.path(root_out, "rhd_country_comparison_2025.rds"))
  message("  Wrote country comparison: ", file.path(root_out, "rhd_country_comparison_2025.csv"))
}

# ------------------------------------------------------------------------------
# 11. REPORT
# ------------------------------------------------------------------------------
t <- econ_summary[location == "TOTAL"]
fm <- function(x, d = 2) if (is.finite(x)) formatC(x, format = "f", digits = d, big.mark = ",") else "NA"
message("── 08_economic_evaluation.R : ", par$country, " (constant 2025 US$) ─────")
message(sprintf("  Unit costs 2025 US$: screen=%s  SAP=%s/py  surgery=%s  (from %s / %s / %s assumed-2019 US$; factor %s)",
                fm(uc[[1]], 3), fm(uc[[2]]), fm(uc[[3]], 0), par$bundles$screen$original_value,
                par$bundles$sap$original_value, par$bundles$surgery$original_value, fm(price_audit$adjustment_factor[1], 4)))
message(sprintf("  VSL: raw anchor %s int$%d (floor %s) -> %d: US$%s, %d: US$%s; floor binding in %d year(s)",
                fm(A$raw_vsl_anchor_intY, 0), A$income_anchor_year, fm(A$floor_vsl_anchor_intY, 0),
                min(analysis_years), fm(vsl_tab[1, vsl_usd2025], 0), max(analysis_years), fm(vsl_tab[.N, vsl_usd2025], 0),
                sum(vsl_tab$vsl_floor_binding)))
message("  ", vsl_vs_disc_note)
message(sprintf("  %-26s %16s %16s", "", "undiscounted", report_labels$discounting_label))
for (k in list(c("Incremental cost ($m)", "total_cost"), c("Monetised benefit ($m)", "total_benefit"), c("Net benefit ($m)", "net_benefit")))
  message(sprintf("  %-26s %16s %16s", k[1], fm(t[[paste0(k[2], "_undisc")]] / 1e6, 1), fm(t[[paste0(k[2], "_disc")]] / 1e6, 1)))
message(sprintf("  %-26s %16s %16s", "Benefit-cost ratio", fm(t$bcr_undisc), fm(t$bcr_disc)))
message(sprintf("  %-26s %16s %16s", sprintf("Per resident %d (US$)", t$N0_year), fm(t$net_benefit_pc_undisc), fm(t$net_benefit_pc_disc)))
message(sprintf("  N0 = %s | deaths averted %s | cost/death $%s | cost/DALY $%s | INT$ per US$ %s",
                formatC(round(t$N0), format = "d", big.mark = ","), formatC(round(t$deaths_averted), format = "d", big.mark = ","),
                fm(t$cost_per_death, 0), fm(t$cost_per_daly, 0), fm(intf$factor, 4)))
message(sprintf("  NOTE (UNCOSTED): %s incremental confirmatory echoes over %d-%d are NOT priced; total cost is understated.",
                formatC(round(sum(budget_impact$confirmations_incremental_uncosted)), format = "d", big.mark = ","),
                min(analysis_years), max(analysis_years)))
message("  Calibration status (model outputs): ", substr(calibration_status, 1, 160), if (nchar(calibration_status) > 160) " ..." else "")
message("  Sensitivity analysis: ", SENSITIVITY_STATUS, " (VSL reference/elasticity/floor, income growth, discounting, bundle prices/price years)")
message("── 08_economic_evaluation.R complete ── wrote ", OUT_DIR, ": rhd_economic_summary.csv, rhd_budget_impact.csv, ",
        "rhd_economic_results.rds, rhd_economic_annual_per_capita.csv, rhd_vsl_derivation.csv, ",
        "rhd_price_adjustment_audit.csv, rhd_economic_parameter_register.csv")
