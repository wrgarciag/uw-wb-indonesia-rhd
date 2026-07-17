#...........................................................
# GBD 2023 Data ----
#...........................................................

# Standalone-sourcing guard (honour 00_run_all.R globals when present).
library(data.table)
if (!exists("wd_raw"))   wd_raw   <- paste0(here::here("data-raw"), "/")
if (!exists("wd_data"))  wd_data  <- paste0(here::here("data"), "/")
if (!exists("LOCATION")) LOCATION <- "Indonesia"   # 00_run_all.R sets this per COUNTRY
dir.create(wd_data, recursive = TRUE, showWarnings = FALSE)

# GBD Level-2 country = LOCATION

# https://collab2023.healthdata.org/gbd-results?params=gbd-api-2023-permalink/e23ae880499d3ae1d643dca738057425

# load 1990-2023

# List all CSV files
files <- list.files(paste0(wd_raw,"epidemiology/"), pattern = "\\.csv$", full.names = TRUE)

# Read and combine using rbindlist
dt_23 <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)

dt<-data.table(dt_23)

# RETAIN upper/lower (GBD 95% uncertainty interval) on the epidemiology file: the
# calibration (04) uses them for OPTIONAL inverse-variance weighting of the
# log-scale residuals (Part D). They are dropped from the population file below
# (not needed there). Guarded so 01 still runs if a source CSV lacks them.
if (!all(c("upper", "lower") %in% names(dt))) {
  if (!"upper" %in% names(dt)) dt[, upper := NA_real_]
  if (!"lower" %in% names(dt)) dt[, lower := NA_real_]
}

unique(dt$year)
unique(dt$location_name)
unique(dt$cause_name)
unique(dt$age_name)

# Fix countries names

# # Remove unnecessary dx
dx_include <- c("All causes","Rheumatic heart disease")

cause_map <- c(
  rhd      = "Rheumatic heart disease",
  all      = "All causes"
)

# AFTER  – define the vector once, reuse it
cause_cols <- names(cause_map)

# Filter the data to include only the specified causes
dt <- dt[cause_name %in% dx_include,]

# Filter only the target country (00_run_all.R sets LOCATION)
dt <- dt[location_name == LOCATION,]

# save temp baseline rates from gbd 2023 (COUNTRY-scoped wd_data)
saveRDS(dt, file = paste0(wd_data,"temp_baseline_rates_gbd.rds"))

#...........................................................
# GBD 2023 Population ----
#...........................................................

# load 1990-2023

# List all CSV files
files <- list.files(paste0(wd_raw,"population/"), pattern = "\\.csv$", full.names = TRUE)

# Read and combine using rbindlist
pop_23 <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)

dt_pop<-data.table(pop_23)

dt_pop[, upper:=NULL]
dt_pop[, lower:=NULL]

dt_pop <- dt_pop[age_name!="80+ years",]

# Filter only the target country
dt_pop <- dt_pop[location_name == LOCATION,]

# save temp population from gbd 2023 (was mistakenly saving `dt` = epi rates)
saveRDS(dt_pop, file = paste0(wd_data,"temp_population_gbd.rds"))

# Here make population projecions based on gbd population

# Here conduct a GLM projection model of all cause, background and disease specific
# Rates
