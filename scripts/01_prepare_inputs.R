#...........................................................
# GBD 2023 Data ----
#...........................................................

# Indonesia Level 2

# https://collab2023.healthdata.org/gbd-results?params=gbd-api-2023-permalink/e23ae880499d3ae1d643dca738057425

# load 1990-2023

# List all CSV files
files <- list.files(paste0(wd_raw,"epidemiology/"), pattern = "\\.csv$", full.names = TRUE)

# Read and combine using rbindlist
dt_23 <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)

dt<-data.table(dt_23)

dt[, upper:=NULL]
dt[, lower:=NULL]

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

# Filter only Indonesia
#dt <- dt[location_name == "Indonesia",]

# save temp baseline rates from gbd 2023
saveRDS(dt, file = paste0(wd_raw,"temp_baseline_rates_gbd.rds"))

