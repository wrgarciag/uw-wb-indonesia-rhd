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

wd <- "C:/Users/wrgar/OneDrive - UW/02Work/WorldBank-Indonesia/uw-wb-indonesia-ncd/"

wd_code <- paste0(wd,"code/")

# Raw data not available on GitHub
wd_raw <- paste0(wd,"data/raw/")

# Processed data (from base rates and tps)
wd_data <- paste0(wd,"data/processed/")
wd_outp <- paste0(wd,"output/")

# Create a temporary directory for the processing data change to wd in final version
wd_temp <- paste0("C:/Users/wrgar/OneDrive - UW/02Work/WorldBank-Indonesia/","temp/")
if (!dir.exists(wd_temp)) {
  dir.create(wd_temp, recursive = TRUE)
}

setwd(paste0(wd_code,"cvd/"))

#...........................................................
# 0. Functions and parameters-----
#...........................................................

source("01_utils.R")

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
                "Ischemic heart disease",
                "Ischemic stroke",
                "Intracerebral hemorrhage",
                "Hypertensive heart disease",
                "Rheumatic heart disease",
                "Cardiomyopathy and myocarditis"                
                )

cause_map <- c(
  ihd      = "Ischemic heart disease",
  istroke  = "Ischemic stroke",
  hstroke  = "Intracerebral hemorrhage",
  hhd      = "Hypertensive heart disease",
  rhd     = "Rheumatic heart disease",
  cmd     = "Cardiomyopathy and myocarditis",
  all      = "All causes"
)

# AFTER  – define the vector once, reuse it
cause_cols <- names(cause_map)

#...........................................................
# 02. Load inputs-----
#...........................................................

source("02_load_inputs_indonesia.R")

#...........................................................
# 03. Clean and process inputs-----
#...........................................................

source("03_calibration.R")

#...........................................................
# 04. define interventions ----
#...........................................................

source("04_define_interventions.R")

#...........................................................
# 05. build baseline ----
#...........................................................

# Run CVD multiple interventions
setwd(paste0(wd_code,"cvd/"))
source("05_build_baseline_indonesia.R")

#...........................................................
# 06. Run model ----
#...........................................................

# Run CVD multiple interventions
setwd(paste0(wd_code,"cvd/"))
source("06_run_scenarios_indonesia_fair.R")

#...........................................................
# 07. Run Burden of Disease ----
#...........................................................

setwd(paste0(wd_code,"cvd/"))
source("07_output_dalys.R")

#...........................................................
# 08. Run Economic Value ----
#...........................................................
setwd(paste0(wd_code,"cvd/"))
source("08_economic_value_calculation.R")


#...........................................................
# 09. Run Validation ----
#...........................................................
setwd(paste0(wd_code,"cvd/"))
source("09_validation_indonesia.R")
