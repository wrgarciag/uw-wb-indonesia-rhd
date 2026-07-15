################################################################################
# INDONESIA RHD MODEL — POPULATION BACKBONE FROM WPP2024
# scripts/02_build_demography.R
# ------------------------------------------------------------------------------
# Builds the single-year-of-age population backbone that the RHD disease-model
# input builder (03_build_disease_model.R), the calibration (04) and the initial
# state assembler (05) consume. Two tidy tables plus a life-table helper and a
# country lookup.
#
# INPUTS  (wpp2024 R package datasets — data only; devtools::install_github('PPgp/wpp2024'))
#   popAge1dt      "observed" single-year population, 1949-2023, ages 0-100 (100 = 100+)
#   popprojAge1dt  projected single-year population (medium variant), 2024-2100, ages 0-100
#   mx1dt          single-year mortality rates (kept only for the get.lt helper)
#   UNlocations    country lookup (-> `locations`)
#
# OUTPUTS (written to wd_data = data/):
#   pop_observed_1990_2024.rds / .csv
#       tidy long {location, year, sex, age, Nx}; 35 yr x 2 sex x 101 age = 7,070 rows.
#       Single-year age distribution SPLINE-SMOOTHED (removes single-year sawtooth /
#       age-heaping), then rescaled so each (sex, year) sums EXACTLY to the source
#       WPP all-age total for that sex-year (totals preserved).
#   pop_projection_2025_2100.rds / .csv
#       tidy long {location, year, sex, age, Nx}; 76 yr x 2 sex x 101 age = 15,352 rows.
#       WPP2024 medium-variant projection, taken directly (UN single-year projections
#       are already graduated, so no extra smoothing is applied).
#   wpp/indonesia_rhd_demography.Rda
#       get.lt (life-table fn), locations (lookup), and copies of both tables —
#       for any downstream script that still expects the .Rda bundle.
#
# YEAR RANGES : observed 1990-2024 ; projection 2025-2100.  This CLOSES the former
#   2024-2025 gap (previously observed 1990-2023 + projection 2026-2100).  The split
#   mirrors the UNWPP2024 structure: WPP2024 estimates run through 2024 and the
#   medium-variant projection proper begins 2025.
#
#   IMPORTANT provenance note: in the wpp2024 R PACKAGE, popAge1dt (the "observed"
#   table) ends at 2023, and 2024 is stored as the FIRST row of popprojAge1dt (the
#   projection jump-off / base year). So the 2024 "last observed" slice is taken
#   from popprojAge1dt, and the projection proper is 2025-2100 (also popprojAge1dt).
#   Because 2024 and 2025 are consecutive rows of the SAME series, the observed ->
#   projection join (2024 -> 2025) is internally continuous (validated below).
# AGE         : single year 0..100; 100 is the WPP terminal OPEN age group (100+).
# Nx UNITS    : persons (WPP reports thousands; multiplied by 1e3 here).
#
# CHANGES FROM THE PARENT-CCPM VERSION
#   The former get.par() migration back-solve and the 2020-anchored sf.wpp arrays
#   (mx/mig/asfr/srb/base.pop) are retired. get.par() only produced projection
#   PARAMETERS for a forward engine (engine.R) that is NOT present in this repo,
#   was anchored at 2020, and indexed pop[y-2024,,] (negative indices for y<2024).
#   Because the deliverables here are the two population TABLES, the projection is
#   taken directly from WPP2024's own medium-variant series — internally consistent
#   with the observed series and requiring no absent engine. get.lt() and
#   `locations` are retained per spec.
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

# Year windows (honour globals from 00_run_all.R when present, else defaults).
# OBS_YEARS closes the former 2024-2025 gap by extending observed through 2024;
# PROJ_YEARS is the medium-variant projection proper, 2025-2100.
if (!exists("OBS_YEARS"))  OBS_YEARS  <- 1990:2024   # observed segment (2024 = WPP jump-off)
if (!exists("PROJ_YEARS")) PROJ_YEARS <- 2025:2100   # medium-variant projection segment

# In the wpp2024 package, popAge1dt (observed) ends at 2023 and popprojAge1dt
# begins at 2024. Split the requested observed window accordingly.
OBS_YEARS_FROM_OBS  <- OBS_YEARS[OBS_YEARS <= 2023L]         # taken from popAge1dt
OBS_YEARS_FROM_PROJ <- OBS_YEARS[OBS_YEARS >= 2024L]         # taken from popprojAge1dt (jump-off)

TERMINAL_AGE <- 100L         # WPP single-year terminal OPEN age group (100+)
AGES <- 0:TERMINAL_AGE       # 101 single-year ages

# self-describing labels used in filenames / object names / messages
OBS_LABEL  <- sprintf("%d_%d", min(OBS_YEARS),  max(OBS_YEARS))    # "1990_2024"
PROJ_LABEL <- sprintf("%d_%d", min(PROJ_YEARS), max(PROJ_YEARS))   # "2025_2100"

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

# helper: pull a single-year raw slice (persons) from a wpp2024 table
pull_wpp_slice <- function(tbl, years_keep) {
  as.data.frame(tbl) |>
    # popAge1dt$year is integer but popprojAge1dt$year is character; coerce both.
    mutate(year = as.integer(year), age = as.integer(age)) |>
    filter(name == idn_name, year %in% years_keep, age %in% AGES) |>
    transmute(location = name, year, age,
              Female = popF * 1e3, Male = popM * 1e3)
}

################################################################################
# 4  OBSERVED single-year population, 1990-2024  (spline-smoothed, totals kept)
#    1990-2023 from popAge1dt; 2024 from popprojAge1dt (WPP jump-off / base year).
################################################################################

obs_raw <- bind_rows(
  pull_wpp_slice(popAge1dt,     OBS_YEARS_FROM_OBS),   # 1990-2023
  pull_wpp_slice(popprojAge1dt, OBS_YEARS_FROM_PROJ)   # 2024 (jump-off in proj series)
) |>
  pivot_longer(c(Female, Male), names_to = "sex", values_to = "Nx") |>
  arrange(location, year, sex, age)

if (length(OBS_YEARS_FROM_PROJ))
  message("  Observed 2024 slice sourced from popprojAge1dt (WPP jump-off year).")

# 4b  spline-smooth the single-year age distribution within each (location,year,sex),
#     then rescale so the smoothed counts sum to the source all-age total.
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

pop_observed_1990_2024 <- obs_raw |>
  group_by(location, year, sex) |>
  group_modify(~ smooth_age_distribution(.x)) |>
  ungroup() |>
  select(location, year, sex, age, Nx) |>
  arrange(location, year, sex, age)

# diagnostic: how much did smoothing move any single-age cell?
smooth_shift <- obs_raw |>
  rename(Nx_raw = Nx) |>
  inner_join(pop_observed_1990_2024, by = c("location", "year", "sex", "age")) |>
  mutate(rel = abs(Nx - Nx_raw) / pmax(Nx_raw, 1))
message(sprintf("  Observed %s smoothed | max single-age shift = %.1f%% | mean = %.2f%%",
                OBS_LABEL, 100 * max(smooth_shift$rel), 100 * mean(smooth_shift$rel)))

################################################################################
# 5  PROJECTION single-year population, 2025-2100  (WPP medium variant, direct)
################################################################################

pop_projection_2025_2100 <- pull_wpp_slice(popprojAge1dt, PROJ_YEARS) |>
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

validate_pop(pop_observed_1990_2024,   OBS_YEARS,  paste0("pop_observed_", OBS_LABEL))
validate_pop(pop_projection_2025_2100, PROJ_YEARS, paste0("pop_projection_", PROJ_LABEL))

# totals preserved by smoothing (observed): each sex-year smoothed sum == raw sum
tot_check <- obs_raw |>
  group_by(location, year, sex) |>
  summarise(raw = sum(Nx), .groups = "drop") |>
  inner_join(pop_observed_1990_2024 |>
               group_by(location, year, sex) |>
               summarise(sm = sum(Nx), .groups = "drop"),
             by = c("location", "year", "sex")) |>
  mutate(rel = abs(sm - raw) / raw)
if (max(tot_check$rel) > 1e-6)
  stop("Smoothing did not preserve sex-year totals (max rel err = ",
       signif(max(tot_check$rel), 3), ").", call. = FALSE)

# sane national totals: Indonesia ~275M around 2020
tot2020 <- pop_observed_1990_2024 |> filter(year == 2020) |> summarise(t = sum(Nx)) |> pull(t)
message(sprintf("  Total population 2020 (observed): %s",
                formatC(round(tot2020), format = "d", big.mark = ",")))
if (tot2020 < 2e8 || tot2020 > 3.5e8)
  stop("Indonesia 2020 total = ", round(tot2020 / 1e6),
       "M, outside the sane band 200-350M.", call. = FALSE)

# no observed->projection level discontinuity across the 2024/2025 join
join_obs_yr  <- max(OBS_YEARS)     # 2024
join_proj_yr <- min(PROJ_YEARS)    # 2025
totA <- pop_observed_1990_2024   |> filter(year == join_obs_yr)  |> summarise(t = sum(Nx)) |> pull(t)
totB <- pop_projection_2025_2100 |> filter(year == join_proj_yr) |> summarise(t = sum(Nx)) |> pull(t)
join_ratio <- totB / totA
message(sprintf("  Join check: total %d obs = %s -> %d proj = %s (ratio %.4f)",
                join_obs_yr,  formatC(round(totA), format = "d", big.mark = ","),
                join_proj_yr, formatC(round(totB), format = "d", big.mark = ","),
                join_ratio))
# consecutive WPP years: expect ~1 year of growth (tight band since same series)
if (join_ratio < 0.97 || join_ratio > 1.05)
  stop("Observed/projection join looks discontinuous (", join_proj_yr, "/", join_obs_yr,
       " ratio = ", round(join_ratio, 4), ", expected ~1.00-1.02).", call. = FALSE)

message("  All demography validation checks passed ✓")

################################################################################
# 7  SAVE
################################################################################

f_obs  <- paste0(wd_data, "pop_observed_",   OBS_LABEL)
f_proj <- paste0(wd_data, "pop_projection_", PROJ_LABEL)

saveRDS(pop_observed_1990_2024,   file = paste0(f_obs,  ".rds"))
saveRDS(pop_projection_2025_2100, file = paste0(f_proj, ".rds"))
readr::write_csv(pop_observed_1990_2024,   paste0(f_obs,  ".csv"))
readr::write_csv(pop_projection_2025_2100, paste0(f_proj, ".csv"))

save(get.lt, locations, pop_observed_1990_2024, pop_projection_2025_2100, file = OUT_RDA)

message("\n── Saved ──────────────────────────────────")
message("  pop_observed_",   OBS_LABEL,  ".rds/.csv   [", nrow(pop_observed_1990_2024),  " rows]  years ",
        min(pop_observed_1990_2024$year),   "-", max(pop_observed_1990_2024$year))
message("  pop_projection_", PROJ_LABEL, ".rds/.csv [", nrow(pop_projection_2025_2100), " rows]  years ",
        min(pop_projection_2025_2100$year), "-", max(pop_projection_2025_2100$year))
message("  wpp/indonesia_rhd_demography.Rda  (get.lt, locations, both tables)")
message("\n── 02_build_demography.R complete ────────────────────")
message("  Next: 03_build_disease_model.R")
