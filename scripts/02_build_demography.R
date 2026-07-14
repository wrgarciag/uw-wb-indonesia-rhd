

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

# save temp baseline rates from gbd 2023
saveRDS(dt, file = paste0(wd_raw,"temp_population_gbd.rds"))

# Here make population projecions based on gbd population

# Here conduct a GLM projection model of all cause, background and disease specific
# Rates




