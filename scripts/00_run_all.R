rm(list=ls()) 

#libraries
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)   
library(countrycode)
library(stringr)
library(parallel)
library(doParallel)
library(foreach)
library(gmodels)

# For forecasting mortality
library(forecast)

wd <- "C:/Users/wrgar/OneDrive - UW/02Work/WorldBank-Indonesia/uw-wb-indonesia-rhd/"

wd_code <- paste0(wd,"scripts/")

# Raw data not available on GitHub
wd_raw <- paste0(wd,"data-raw/")

# Processed data (from base rates and tps)
wd_data <- paste0(wd,"data/")
wd_outp <- paste0(wd,"output/")

# Create a temporary directory for the processing data change to wd in final version
wd_temp <- paste0("C:/Users/wrgar/OneDrive - UW/02Work/WorldBank-Indonesia/","temp-rhd/")
if (!dir.exists(wd_temp)) {
  dir.create(wd_temp, recursive = TRUE)
}

setwd(paste0(wd_code))

#...........................................................
# 0. Functions and parameters-----
#...........................................................

#source("01_utils.R")

run_calibration_par <- TRUE # set to TRUE to run parallel calibration

run_adjustment_model <- FALSE # set to TRUE to run adjustment model

run_aod_par <- FALSE # set to TRUE to run model with dementia

run_adjustments_inputs <- FALSE

run_bgmx_trend <- TRUE

run_CF_trend   <- TRUE

# Baseline scenario 80% of secular trend. 20% historically explained by
# HTN control improvements

run_CF_trend_80   <- TRUE

run_CF_trend_ihme  <- FALSE

# Remove unnecessary dx

dx_include <- c("All causes",
                "Rheumatic heart disease"
                )

cause_map <- c(
  rhd     = "Rheumatic heart disease",
  all      = "All causes"
)

# AFTER  – define the vector once, reuse it
cause_cols <- names(cause_map)

#...........................................................
# NOTE: sourcing updated to the ACTUAL script filenames in scripts/.
# The previous list (03_calibration.R, 04_define_interventions.R,
# cvd/05..09_*.R) referred to parent-NCD files that do not exist in this repo.
#...........................................................

setwd(wd_code)   # scripts are sourced from scripts/ ; each uses absolute wd_* / here() for I/O

#...........................................................
# 01. Load & filter GBD inputs -> data-raw/temp_baseline_rates_gbd.rds
#...........................................................
source("01_prepare_inputs.R")

#...........................................................
# 02. Build WPP2024 population backbone
#     -> data/pop_observed_1990_2023.rds, data/pop_projection_2026_2100.rds
#...........................................................
source("02_build_demography.R")

#...........................................................
# 03. RHD disease model (data-fed) -> in-memory `ref`/`sap` + economic globals
#...........................................................
source("03_build_disease_model.R")

#...........................................................
# 04. Random-search TP calibration (RHD-native) -> data/adjusted_searo_part*.rds
#     Keep run_adjustment_model == FALSE (this bakes multipliers into IR/CF).
#...........................................................
source("04_calibration_random_tp.R")

#...........................................................
# 05-06. Parent-NCD baseline (05_build_baseline.R) and the toy prevention model
#     (06_run_prevention_model.R, a duplicate of 03) are NOT yet ported to the
#     RHD pipeline and reference inputs absent from this repo. Skipped by default;
#     flip run_downstream_ncd <- TRUE once they are ported.
#...........................................................
run_downstream_ncd <- FALSE
if (run_downstream_ncd) {
  source("05_build_baseline.R")
  source("06_run_prevention_model.R")
}

#...........................................................
# 07. Make outputs (budget impact, BCR, DALYs) from 03's `ref`/`sap`
#     -> written into output/
#...........................................................
setwd(wd_outp)
source(paste0(wd_code, "07_make_outputs.R"))
