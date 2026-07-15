################################################################################
# INDONESIA RHD MODEL — POPULATION BACKBONE FROM WPP2024
# scripts/02_build_demography.R
# ------------------------------------------------------------------------------
# Builds the single-year-of-age population backbone the RHD disease model
# (03_build_disease_model.R) and calibration (04) consume. Two tidy tables plus
# a life-table helper and a country lookup.
#
# INPUTS  (wpp2024 R package datasets — data only; devtools::install_github('PPgp/wpp2024'))
#   popAge1dt      observed single-year population, 1949-2023, ages 0-100 (100 = 100+)
#   popprojAge1dt  projected single-year population (medium variant), 2024-2100, ages 0-100
#   mx1dt          single-year mortality rates (kept only for the get.lt helper)
#   UNlocations    country lookup (-> `locations`)
#
# OUTPUTS (written to wd_data = data/):
#   pop_observed_1990_2023.rds / .csv
#       tidy long {location, year, sex, age, Nx}; 34 yr x 2 sex x 101 age = 6,868 rows.
#       Single-year age distribution SPLINE-SMOOTHED (removes single-year sawtooth /
#       age-heaping), then rescaled so each (sex, year) sums EXACTLY to the observed
#       WPP all-age total for that sex-year (totals preserved).
#   pop_projection_2026_2100.rds / .csv
#       tidy long {location, year, sex, age, Nx}; 75 yr x 2 sex x 101 age = 15,150 rows.
#       WPP2024 medium-variant projection, taken directly (UN single-year projections
#       are already graduated, so no extra smoothing is applied).
#   wpp/indonesia_rhd_demography.Rda
#       get.lt (life-table fn), locations (lookup), and copies of both tables —
#       for any downstream script that still expects the .Rda bundle.
#
# YEAR RANGES : observed 1990-2023 ; projection 2026-2100. WPP2024 is internally
#   continuous across the 2024-2025 gap (same vintage), so the observed and
#   projected segments join without a level discontinuity (validated below).
# AGE         : single year 0..100; 100 is the WPP terminal OPEN age group (100+).
# Nx UNITS    : persons (WPP reports thousands; multiplied by 1e3 here).
#
# CHANGES FROM THE PARENT-CCPM VERSION
#   The former get.par() migration back-solve and the 2020-anchored sf.wpp arrays
#   (mx/mig/asfr/srb/base.pop) are retired. get.par() only produced projection
#   PARAMETERS for a forward engine (engine.R) that is NOT present in this repo,
#   was anchored at 2020 (this task needs 2017-2100), and indexed pop[y-2024,,]
#   (negative indices for y<2024). Because the deliverables here are the two
#   population TABLES, the projection is taken directly from WPP2024's own
#   medium-variant series — internally consistent with the observed series and
#   requiring no absent engine. get.lt() and `locations` are retained per spec.
#
# Run ONCE per WPP version. R >= 4.1.0 (native pipe).
################################################################################

if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required. Install with: install.packages('here')", call. = FALSE)
source(here::here("R", "packages.R"))
library(here)

# NOTE: load data.table BEFORE wpp2024. Loading wpp2024 first and data.table
# second triggers an OpenMP-runtime conflict on Windows that segfaults the first
# [.data.table call. (00_run_all.R already loads data.table first.)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

if (!requireNamespace("wpp2024", quietly = TRUE))
  stop("wpp2024 not installed.\n  Run: devtools::install_github('PPgp/wpp2024')", call. = FALSE)
if (!requireNamespace("countrycode", quietly = TRUE))
  stop("countrycode not installed.\n  Run: install.packages('countrycode')", call. = FALSE)

suppressPackageStartupMessages(library(wpp2024))

# Output directory: honour the wd_data global from 00_run_all.R when present,
# else resolve via here() for standalone sourcing.
if (!exists("wd_data")) wd_data <- paste0(here::here("data"), "/")
dir.create(file.path(wd_data, "wpp"), recursive = TRUE, showWarnings = FALSE)
OUT_RDA <- file.path(wd_data, "wpp", "indonesia_rhd_demography.Rda")

################################################################################
# 1  CONSTANTS
################################################################################

MODEL_COUNTRIES <- tribble(
  ~iso3,  ~wpp_name,
  "IDN",  "Indonesia"
)

OBS_YEARS  <- 1990:2023      # observed segment
PROJ_YEARS <- 2026:2100      # projection segment
TERMINAL_AGE <- 100L         # WPP single-year terminal OPEN age group (100+)
AGES <- 0:TERMINAL_AGE       # 101 single-year ages

################################################################################
# 2  LIFE TABLE FUNCTION  (retained from parent CCPM; saved for downstream use)
################################################################################

get.lt <- function(mx, parm = NULL, dws = NULL, qual.adj = NULL) {
  mx    <- mx[!is.na(mx)]
  n_age <- length(mx)
  nx    <- rep(1, n_age)
  qx    <- 1 - exp(-nx * mx); qx[n_age] <- 1
  ax    <- (nx + 1 / mx - nx / qx)
  lx    <- c(1, cumprod(1 - qx)[1:(n_age - 1)])
  dx    <- c(rev(diff(rev(lx))), lx[1] - sum(rev(diff(rev(lx)))))
  nLx   <- nx * lx - (nx - ax) * dx
  Tx    <- rev(cumsum(rev(nLx)))
  ex    <- Tx / lx

  if (is.null(qual.adj)) {
    Sx        <- nLx / c(1, nLx)[1:n_age]
    Sx[n_age] <- nLx[n_age] / (nLx[n_age - 1] + nLx[n_age])
    df.lt <- data.table(age = 0:(n_age - 1), ax, mx, lx, qx, dx, nLx, Tx, ex, Sx)
    if (is.null(parm)) df.lt else df.lt |> pull(parm)
  } else {
    qual  <- 1 - dws; qual <- qual[!is.na(qual)]
    nLxh  <- nLx * qual
    Txh   <- rev(cumsum(rev(nLxh)))
    exh   <- Txh / lx
    exh[1]
  }
}

################################################################################
# 3  LOAD WPP2024 PACKAGE DATA + LOCATION LOOKUP
################################################################################

message("\n── Loading WPP2024 package data ─────────────────────")

data(popAge1dt)
data(popprojAge1dt)
data(mx1dt)          # loaded so get.lt has a companion source if needed downstream
data(UNlocations)

locations <- suppressWarnings(
  UNlocations |>
    mutate(iso3 = countrycode::countrycode(name, "country.name", "iso3c")) |>
    rename(location_name = name) |>
    filter(!is.na(iso3),
           location_name != "Less developed regions, excluding China",
           country_code %in% unique(popAge1dt$country_code)) |>
    select(location_name, country_code, iso3)
) |>
  filter(iso3 %in% MODEL_COUNTRIES$iso3)

if (!"IDN" %in% locations$iso3)
  stop("Indonesia (IDN) not resolved in UNlocations. Check wpp2024 install.")

idn_name <- MODEL_COUNTRIES$wpp_name[MODEL_COUNTRIES$iso3 == "IDN"]

message("  Locations resolved: ", paste(locations$iso3, collapse = ", "))

################################################################################
# 4  OBSERVED single-year population, 1990-2023  (spline-smoothed, totals kept)
################################################################################

# 4a  raw observed long table (persons)
obs_raw <- as.data.frame(popAge1dt) |>
  # wpp2024 stores popAge1dt$year as integer but popprojAge1dt$year as character;
  # coerce both to integer age/year so the saved tables have consistent types.
  mutate(year = as.integer(year), age = as.integer(age)) |>
  filter(name == idn_name, year %in% OBS_YEARS, age %in% AGES) |>
  transmute(location = name, year, age,
            Female = popF * 1e3, Male = popM * 1e3) |>
  pivot_longer(c(Female, Male), names_to = "sex", values_to = "Nx") |>
  arrange(location, year, sex, age)

# 4b  spline-smooth the single-year age distribution within each (location,year,sex),
#     then rescale so the smoothed counts sum to the observed all-age total.
smooth_age_distribution <- function(df) {
  df <- df[order(df$age), ]
  raw_total <- sum(df$Nx)
  # smooth.spline needs >= 4 distinct x; ages 0..100 always satisfy this
  fit <- stats::smooth.spline(x = df$age, y = df$Nx)
  sm  <- stats::predict(fit, df$age)$y
  sm[sm < 0 | !is.finite(sm)] <- 0             # age heaping can push tails <0
  if (sum(sm) > 0) sm <- sm * (raw_total / sum(sm))   # PRESERVE the sex-year total
  df$Nx <- sm
  df
}

pop_observed <- obs_raw |>
  group_by(location, year, sex) |>
  group_modify(~ smooth_age_distribution(.x)) |>
  ungroup() |>
  select(location, year, sex, age, Nx) |>
  arrange(location, year, sex, age)

# diagnostic: how much did smoothing move any single-age cell?
smooth_shift <- obs_raw |>
  rename(Nx_raw = Nx) |>
  inner_join(pop_observed, by = c("location", "year", "sex", "age")) |>
  mutate(rel = abs(Nx - Nx_raw) / pmax(Nx_raw, 1))
message(sprintf("  Observed 1990-2023 smoothed | max single-age shift = %.1f%% | mean = %.2f%%",
                100 * max(smooth_shift$rel), 100 * mean(smooth_shift$rel)))

################################################################################
# 5  PROJECTION single-year population, 2026-2100  (WPP medium variant, direct)
################################################################################

pop_projection <- as.data.frame(popprojAge1dt) |>
  mutate(year = as.integer(year), age = as.integer(age)) |>   # year is character in wpp2024 proj
  filter(name == idn_name, year %in% PROJ_YEARS, age %in% AGES) |>
  transmute(location = name, year, age,
            Female = popF * 1e3, Male = popM * 1e3) |>
  pivot_longer(c(Female, Male), names_to = "sex", values_to = "Nx") |>
  select(location, year, sex, age, Nx) |>
  arrange(location, year, sex, age)

################################################################################
# 6  VALIDATION  (fail loudly — stop() — before writing anything)
################################################################################

message("\n── Validation ──────────────────────────────")

validate_pop <- function(dt, years_expected, label) {
  # completeness: every (year, sex) has all 101 single ages 0..100
  ages_ok <- dt |>
    group_by(year, sex) |>
    summarise(n_age = n(),
              amin = min(age), amax = max(age), .groups = "drop")
  if (any(ages_ok$n_age != length(AGES)) ||
      any(ages_ok$amin != 0) || any(ages_ok$amax != TERMINAL_AGE))
    stop(label, ": age grid is not a complete 0..", TERMINAL_AGE,
         " for every sex-year.", call. = FALSE)

  yrs <- sort(unique(dt$year))
  if (!identical(as.integer(yrs), as.integer(years_expected)))
    stop(label, ": year coverage is ", min(yrs), "-", max(yrs),
         " (n=", length(yrs), "), expected ",
         min(years_expected), "-", max(years_expected),
         " complete.", call. = FALSE)

  if (any(is.na(dt$Nx)))  stop(label, ": Nx contains NA.", call. = FALSE)
  if (any(dt$Nx < 0))     stop(label, ": Nx contains negative values.", call. = FALSE)
  if (!all(c("Female", "Male") %in% unique(dt$sex)))
    stop(label, ": both sexes (Female, Male) must be present.", call. = FALSE)
  invisible(TRUE)
}

validate_pop(pop_observed,   OBS_YEARS,  "pop_observed_1990_2023")
validate_pop(pop_projection, PROJ_YEARS, "pop_projection_2026_2100")

# totals preserved by smoothing (observed): each sex-year smoothed sum == raw sum
tot_check <- obs_raw |>
  group_by(location, year, sex) |>
  summarise(raw = sum(Nx), .groups = "drop") |>
  inner_join(pop_observed |>
               group_by(location, year, sex) |>
               summarise(sm = sum(Nx), .groups = "drop"),
             by = c("location", "year", "sex")) |>
  mutate(rel = abs(sm - raw) / raw)
if (max(tot_check$rel) > 1e-6)
  stop("Smoothing did not preserve sex-year totals (max rel err = ",
       signif(max(tot_check$rel), 3), ").", call. = FALSE)

# sane national totals: Indonesia ~275M around 2020
tot2020 <- pop_observed |> filter(year == 2020) |> summarise(t = sum(Nx)) |> pull(t)
message(sprintf("  Total population 2020 (observed): %s",
                formatC(round(tot2020), format = "d", big.mark = ",")))
if (tot2020 < 2e8 || tot2020 > 3.5e8)
  stop("Indonesia 2020 total = ", round(tot2020 / 1e6),
       "M, outside the sane band 200-350M.", call. = FALSE)

# no observed->projection level discontinuity across the 2023/2026 join
tot2023 <- pop_observed   |> filter(year == 2023) |> summarise(t = sum(Nx)) |> pull(t)
tot2026 <- pop_projection |> filter(year == 2026) |> summarise(t = sum(Nx)) |> pull(t)
join_ratio <- tot2026 / tot2023
message(sprintf("  Join check: total 2023 obs = %s -> 2026 proj = %s (ratio %.3f)",
                formatC(round(tot2023), format = "d", big.mark = ","),
                formatC(round(tot2026), format = "d", big.mark = ","),
                join_ratio))
if (join_ratio < 0.9 || join_ratio > 1.2)
  stop("Observed/projection join looks discontinuous (2026/2023 ratio = ",
       round(join_ratio, 3), ", expected ~1.0-1.1).", call. = FALSE)

message("  All demography validation checks passed ✓")

################################################################################
# 7  SAVE
################################################################################

saveRDS(pop_observed,   file = paste0(wd_data, "pop_observed_1990_2023.rds"))
saveRDS(pop_projection, file = paste0(wd_data, "pop_projection_2026_2100.rds"))
readr::write_csv(pop_observed,   paste0(wd_data, "pop_observed_1990_2023.csv"))
readr::write_csv(pop_projection, paste0(wd_data, "pop_projection_2026_2100.csv"))

save(get.lt, locations, pop_observed, pop_projection, file = OUT_RDA)

message("\n── Saved ──────────────────────────────────")
message("  pop_observed_1990_2023.rds/.csv   [", nrow(pop_observed),  " rows]  years ",
        min(pop_observed$year),   "-", max(pop_observed$year))
message("  pop_projection_2026_2100.rds/.csv [", nrow(pop_projection), " rows]  years ",
        min(pop_projection$year), "-", max(pop_projection$year))
message("  wpp/indonesia_rhd_demography.Rda  (get.lt, locations, both tables)")
message("\n── 02_build_demography.R complete ────────────────────")
message("  Next: 03_build_disease_model.R")
