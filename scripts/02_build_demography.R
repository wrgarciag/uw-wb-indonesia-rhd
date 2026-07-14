################################################################################
# INDONESIA INTEGRATED NCD MODEL — BASELINE DEMOGRAPHY PREPARATION
# scripts/02_build_demography.R
# ─────────────────────────────────────────────────────────────────────────────
# Builds sf.wpp — the WPP2024-based demographic backbone for the NCD model.
# Uses the same get.par() logic as the parent CCPM, which back-solves migration
# as the residual between WPP projected population and a full CCPM forward
# projection (including births).
#
# NOTE: get.par() and get.lt() are defined within this script. They are NOT
# from any external package — they are custom functions from the NCD Countdown
# parent model, adapted here for the Indonesia V1 pipeline.
#
# SOURCE: wpp2024 R package (data only)
#   devtools::install_github("PPgp/wpp2024")
# Also requires: countrycode (install.packages("countrycode"))
#
# INPUT:  wpp2024 package datasets
#
# OUTPUT: data/wpp/indonesia_ncd_demography.Rda
#   sf.wpp     — named list, primary entry IDN (Indonesia).
#                Each element produced by get.par():
#                  $mx       [131 × 2 × 101]  all-cause mortality rates
#                  $base.pop [2 × 101]         population at 2020
#                  $mig      [131 × 2 × 101]  net migration rates
#                  $years    [131]             2020:2155
#                  $asfr     [131 × 45]        age-specific fertility rates
#                  $srb      [131]             sex ratio at birth
#   get.lt     — life table function needed by the projection engine at runtime
#   locations  — country lookup table needed by the projection engine at runtime
#
# Run ONCE. Re-run only if WPP version or country set changes.
################################################################################

if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required. Install with: install.packages('here')", call. = FALSE)
source(here::here("R", "packages.R"))
library(here)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(data.table)
})

if (!requireNamespace("wpp2024", quietly = TRUE))
  stop("wpp2024 not installed.\n  Run: devtools::install_github('PPgp/wpp2024')", call. = FALSE)
if (!requireNamespace("countrycode", quietly = TRUE))
  stop("countrycode not installed.\n  Run: install.packages('countrycode')", call. = FALSE)

suppressPackageStartupMessages(library(wpp2024))

OUT_FILE <- here("data", "wpp", "indonesia_ncd_demography.Rda")
dir.create(here("data", "wpp"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1  CONSTANTS
################################################################################

MODEL_COUNTRIES <- tribble(
  ~iso3,  ~wpp_name,
  "IDN",  "Indonesia"
)

GBD_IDS <- c(IDN = 11)

################################################################################
# 2  LIFE TABLE FUNCTION
# Unchanged from the parent CCPM. Saved alongside sf.wpp so the projection
# engine (engine.R) can use it at runtime without re-sourcing this script.
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
# 3  LOAD WPP2024 PACKAGE DATA
################################################################################

message("\n── Loading WPP2024 package data ─────────────────────────────────────────")

data(popAge1dt)
data(popprojAge1dt)
data(mx1dt)
data(tfr1dt)
data(tfrproj1dt)
data(percentASFR1dt)
data(sexRatio1dt)
data(UNlocations)
data(misc1dt)
data(miscproj1dt)

locations <- suppressWarnings(
  UNlocations |>
    mutate(iso3 = countrycode::countrycode(name, "country.name", "iso3c")) |>
    rename(location_name = name) |>
    filter(!is.na(iso3),
           location_name != "Less developed regions, excluding China",
           country_code %in% unique(popprojAge1dt$country_code)) |>
    select(location_name, country_code, iso3)
) |>
  filter(iso3 %in% MODEL_COUNTRIES$iso3)

message("  Locations resolved: ", paste(locations$iso3, collapse = ", "))

if (!"IDN" %in% locations$iso3)
  stop("Indonesia (IDN) not resolved in UNlocations. Check wpp2024 install.")

pop.dt <- rbind(
  popAge1dt     |> select(country_code, name, year, age, popF, popM),
  popprojAge1dt |> select(country_code, name, year, age, popF, popM)
) |>
  right_join(locations, by = "country_code") |>
  filter(year >= 2020, !is.na(age)) |>
  rename(Female = popF, Male = popM) |>
  select(-c(name, country_code)) |>
  gather(sex, Nx, Female, Male) |>
  mutate(Nx = Nx * 1e3) |>
  arrange(location_name, year, sex, age)

mort.dt <- mx1dt |>
  right_join(locations, by = "country_code") |>
  filter(year >= 2020, !is.na(age)) |>
  rename(Female = mxF, Male = mxM) |>
  select(-c(mxB, name, country_code)) |>
  gather(sex, mx, Female, Male) |>
  arrange(location_name, year, sex, age)

tfr.dt <- rbind(
  tfr1dt     |> select(country_code, name, year, tfr),
  tfrproj1dt |> select(country_code, name, year, tfr)
) |>
  right_join(locations, by = "country_code") |>
  filter(year >= 2020, !is.na(tfr)) |>
  select(location_name, iso3, year, tfr) |>
  arrange(location_name, year)

pasfr.dt <- percentASFR1dt |>
  right_join(locations, by = "country_code") |>
  filter(year >= 2020, !is.na(age)) |>
  select(location_name, iso3, year, age, pasfr) |>
  arrange(location_name, year, age)

srb.dt <- sexRatio1dt |>
  right_join(locations, by = "country_code") |>
  filter(year >= 2020, !is.na(srb)) |>
  select(location_name, iso3, year, srb) |>
  arrange(location_name, year)

births.dt <- rbind(
  misc1dt     |> select(country_code, name, year, births),
  miscproj1dt |> select(country_code, name, year, births)
) |>
  right_join(locations, by = "country_code") |>
  filter(year >= 2020, !is.na(births)) |>
  select(location_name, iso3, year, births) |>
  arrange(location_name, year)

message("  All WPP2024 data tables built \u2713")
message("  Pop years  : ", paste(range(pop.dt$year),   collapse = "\u2013"))
message("  Mort years : ", paste(range(mort.dt$year),  collapse = "\u2013"))
message("  Fert years : ", paste(range(tfr.dt$year),   collapse = "\u2013"))

################################################################################
# 4  get.par() — BACK-SOLVE MIGRATION AND EXTEND TO 2155
#
# For each country, back-solves net migration as the CCPM residual:
#   migration[k] = WPP_population[k] / CCPM_projection[k] - 1
# Then extends mx and mig from 2100 to 2155 (131 years total).
################################################################################

get.par <- function(loc) {
  years   <- 2020:2100
  n_yrs   <- length(years)
  lab_sex <- c("Female", "Male")
  lab_age <- c(0:99, "100+")
  n_sex   <- length(lab_sex)
  n_age   <- length(lab_age)
  
  deaths <- migr.m <- mx.m <- popin <- pop <-
    array(NA_real_, dim = c(n_yrs, n_sex, n_age),
          dimnames = list(Year = years, Sex = lab_sex, Age = lab_age))
  births <- array(NA_real_, dim = c(n_yrs, n_sex),
                  dimnames = list(Year = years, Sex = lab_sex))
  
  p.df  <- pop.dt  |> filter(location_name == loc, sex %in% lab_sex) |>
    spread(age, Nx) |> arrange(year, sex)
  mx.df <- mort.dt |> filter(location_name == loc, sex %in% lab_sex) |>
    spread(age, mx) |> arrange(year, sex)
  
  for (y in years) {
    pop[y - 2024, , ]  <- p.df  |> filter(year == y) |>
      select(-c(year, location_name, iso3, sex)) |> as.matrix()
    mx.m[y - 2024, , ] <- mx.df |> filter(year == y) |>
      select(-c(year, location_name, iso3, sex)) |> as.matrix()
  }
  
  pasfrm <- pasfr.dt |> filter(location_name == loc) |>
    spread(age, pasfr) |>
    select(-c(year, location_name, iso3)) |> as.matrix()
  pasfrm <- t(apply(pasfrm, 1, function(x) x / sum(x, na.rm = TRUE)))
  
  tfrv <- tfr.dt |> filter(location_name == loc) |> pull(tfr)
  srbv <- srb.dt |> filter(location_name == loc) |> pull(srb)
  
  Sx.m <- aperm(apply(mx.m, 1:2, get.lt, "Sx"), c(2, 3, 1))
  
  obsbir <- 1e3 * (births.dt |> filter(location_name == loc) |> pull(births))
  
  asfrm      <- pasfrm
  asfrm[1, ] <- NA
  
  for (k in 2:n_yrs) {
    popin[k - 1, , ]    <- pop[k - 1, , ]
    popin[k, , 2:n_age] <- popin[k - 1, , 1:(n_age - 1)] * Sx.m[k, , 2:n_age]
    popin[k, , n_age]   <- popin[k, , n_age] + popin[k - 1, , n_age] * Sx.m[k, , n_age]
    
    tbirths1   <- sum(0.5 * (popin[k, 1, 11:55] + popin[k - 1, 1, 11:55]) *
                        tfrv[k] * pasfrm[k, ])
    asfrm[k, ] <- obsbir[k] / tbirths1 * tfrv[k] * pasfrm[k, ]
    
    tbirths      <- sum(0.5 * (popin[k, 1, 11:55] + popin[k - 1, 1, 11:55]) *
                          asfrm[k, ])
    births[k, 2] <- tbirths * srbv[k] / (1 + srbv[k])
    births[k, 1] <- tbirths - births[k, 2]
    
    popin[k, , 1]  <- births[k, ] * Sx.m[k, , 1]
    migr.m[k, , ]  <- pop[k, , ] / popin[k, , ] - 1
  }
  
  # Extend from 2100 to 2155 (131 years total)
  n_max <- 131L
  mxm <- mig <- array(dim = c(n_max, 2, 101))
  
  mig[1:n_yrs, , ] <- migr.m
  mig[(n_yrs + 1):n_max, , ] <-
    aperm(replicate(n_max - n_yrs, migr.m[n_yrs, , ]), c(3, 1, 2))
  
  # BUG FIX: pmax guard prevents log(0) = -Inf at young ages with mx = 0
  delta.mxm <- log(pmax(mx.m[n_yrs, , ], 1e-10)) -
    log(pmax(mx.m[n_yrs - 1, , ], 1e-10))
  mxm[1:n_yrs, , ] <- mx.m
  for (y in n_yrs:(n_max - 1))
    mxm[y + 1, , ] <- exp(log(pmax(mxm[y, , ], 1e-10)) + delta.mxm)
  
  srb  <- c(srbv,  replicate(n_max - n_yrs, srbv[n_yrs]))
  asfr <- rbind(asfrm, t(replicate(n_max - n_yrs, asfrm[n_yrs, ])))
  
  list(
    base.pop = pop[1, , ],
    mig      = mig[1:n_max, , ],
    years    = 2020 + 0:(n_max - 1),
    srb      = srb[1:n_max],
    asfr     = asfr[1:n_max, ],
    mx       = mxm[1:n_max, , ]
  )
}

################################################################################
# 5  BUILD sf.wpp
################################################################################

message("\n── Building sf.wpp ──────────────────────────────────────────────────────")

sf.wpp <- list()

for (i in seq_len(nrow(MODEL_COUNTRIES))) {
  iso <- MODEL_COUNTRIES$iso3[i]
  loc <- MODEL_COUNTRIES$wpp_name[i]
  message(sprintf("  [%d/%d]  %s [%s] ...", i, nrow(MODEL_COUNTRIES), loc, iso))
  sf.wpp[[iso]] <- get.par(loc)
  p <- sf.wpp[[iso]]
  message(sprintf("    mx [%s] | mig [%s] | base.pop [%s] | Female mx(60, 2020) = %.4f",
                  paste(dim(p$mx),       collapse = "\u00d7"),
                  paste(dim(p$mig),      collapse = "\u00d7"),
                  paste(dim(p$base.pop), collapse = "\u00d7"),
                  p$mx[1L, 1L, 61L]))
}

################################################################################
# 6  VALIDATION
################################################################################

message("\n── Validation ───────────────────────────────────────────────────────────")

for (iso in names(sf.wpp)) {
  p <- sf.wpp[[iso]]
  stopifnot(
    "mx wrong dims"       = identical(dim(p$mx),       c(131L, 2L, 101L)),
    "mig wrong dims"      = identical(dim(p$mig),      c(131L, 2L, 101L)),
    "base.pop wrong dims" = identical(dim(p$base.pop), c(2L, 101L)),
    "years wrong length"  = length(p$years) == 131L,
    "years start at 2020" = p$years[1]   == 2020L,
    "years end at 2155"   = p$years[131] == 2150L
  )
  if (any(is.na(p$mx[1:76, , ])))  warning(iso, ": NAs in mx (2020-2100)")
  if (any(is.na(p$base.pop)))       warning(iso, ": NAs in base.pop")
  
  total_pop <- sum(p$base.pop)
  message(sprintf("  %s \u2713  years %d\u2013%d | total pop 2020: %s",
                  iso, p$years[1], p$years[131],
                  formatC(round(total_pop), format = "d", big.mark = ",")))
  if (iso == "IDN" && (total_pop < 2e8 || total_pop > 3.5e8))
    warning("IDN total pop 2020 = ", round(total_pop / 1e6), "M — expected ~275M")
}

################################################################################
# 7  CROSS-CHECK: WPP(2020) vs GBD(2023) all-cause mx, ages 30-69
################################################################################
# 
# message("\n── Cross-check: WPP(2020) vs GBD(2023) all-cause mx, ages 30\u201369 ─────────")
# 
# GBD_AC_FILE <- here("data", "gbd", "gbd_allcause_mx.csv")
# if (!file.exists(GBD_AC_FILE)) {
#   message("  GBD file not found \u2014 run 01_prepare_gbd_inputs.R first.")
# } else {
#   gbd_ac <- read_csv(GBD_AC_FILE, show_col_types = FALSE)
#   age_idx <- 31:70
#   
#   cal <- map_dfr(names(sf.wpp), function(iso) {
#     p   <- sf.wpp[[iso]]
#     gid <- GBD_IDS[iso]
#     wf  <- mean(p$mx[1L, 1L, age_idx]) * 1e5
#     wm  <- mean(p$mx[1L, 2L, age_idx]) * 1e5
#     gs  <- gbd_ac |>
#       filter(location_id == gid, year == 2023, age_mid >= 30, age_mid <= 67) |>
#       group_by(sex) |>
#       summarise(m = mean(mx_all, na.rm = TRUE), .groups = "drop")
#     gf <- gs$m[gs$sex == "Female"]; gm <- gs$m[gs$sex == "Male"]
#     tibble(iso3    = iso,
#            wpp_f   = round(wf, 1), gbd_f = round(if (length(gf)) gf else NA_real_, 1),
#            ratio_f = round(wf / max(if (length(gf)) gf else 1, 0.001), 3),
#            wpp_m   = round(wm, 1), gbd_m = round(if (length(gm)) gm else NA_real_, 1),
#            ratio_m = round(wm / max(if (length(gm)) gm else 1, 0.001), 3))
#   })
#   print(cal)
#   bad <- filter(cal, abs(ratio_f - 1) > 0.25 | abs(ratio_m - 1) > 0.25)
#   if (nrow(bad) > 0) {
#     message("  >25% discrepancy: ", paste(bad$iso3, collapse = ", "),
#             " — expected; cause fractions absorb this by design.")
#   } else {
#     message("  All WPP/GBD ratios within 25% \u2713")
#   }
# }

################################################################################
# 8  SAVE
################################################################################

save(sf.wpp, get.lt, locations, file = OUT_FILE)

message("\n── Saved: ", normalizePath(OUT_FILE))
message("  Objects: sf.wpp, get.lt, locations")
message("  sf.wpp[[\"IDN\"]] dims: mx/mig [131\u00d72\u00d7101] | base.pop [2\u00d7101]")
message("\n── 02_build_demography.R complete ───────────────────────────────────────")
message("  Next: 03_build_cause_fractions.R")






