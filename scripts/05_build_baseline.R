# Required inputs----

#...........................................................
#add covid mx data ----
#...........................................................

setwd(wd_raw)
load(paste0(wd_data,"wpp.adj.Rda"))

wpp.adj<-wpp.adj%>%
  mutate(location_name = ifelse(location_name=="North Korea", "Democratic People's Republic of Korea", location_name))
#Covid mx ~= excess mortality

# check for names as GBD and UWPP 2024
locs_wpp.adj <- unique(wpp.adj$location_name)

#baseline rates calculated in file calibration:

files <- list.files(
  path       = wd_data, 
  pattern    = "adjusted", 
  full.names = TRUE
)

dt_list <- lapply(files, function(f) {
  dt <- readRDS(f)
  setDT(dt)  # convert to data.table by reference if it isn't already
  dt
})

# Bind them all together, matching columns by name and filling missing ones
b_rates <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

rm(dt_list, files)

locs_b_rates <- unique(b_rates$location)
b_rates[location=="United States of America",location:="United States"]
b_rates[location=="Bolivia (Plurinational State of)",location:="Bolivia"]
b_rates[location=="United Republic of Tanzania",location:="Tanzania"]

b_rates <- b_rates[!is.na(location),]
b_rates[,c("percent_lag","percent_diff"):=NULL]

b_rates<-left_join(b_rates, wpp.adj%>%
                     rename(location = location_name)%>%
                     select( -Nx, -mx, -iso3))

# Update to UNWPP 2024
dt_pop_unwpp <- as.data.table(readRDS(paste0(wd_data,"PopulationsSingleAge0050.rds")))

dt_pop_unwpp[age>=95, age:= 95]

setnames(dt_pop_unwpp, c("year_id"), c("year"))

dt_pop_unwpp <- dt_pop_unwpp[, .(Nx = sum(Nx)), by = .(location, year, sex, age)]

b_rates <- merge(b_rates,
                 dt_pop_unwpp[, .(location,year,age,sex,Nx2=Nx)],
                 by = c("location", "year", "age","sex"),
                 all.x = TRUE
)

# replace Nx with Nx2
b_rates <- b_rates[, Nx := ifelse(is.na(Nx2), Nx, Nx2)]

b_rates[,Nx2 := NULL]

locs <- unique(b_rates$location)

#...........................................................
# Population data from UNWPP
pop20 <- read.csv(paste0(wd_data,"PopulationsAge20_2050.csv"), stringsAsFactors = F)

b_rates<-left_join(b_rates, pop20%>%rename(Nx2=Nx, year=year_id)%>%filter(year>=2017), 
                   by=c("location", "year", "sex", "age"))%>%
  mutate(Nx = ifelse(is.na(Nx2), Nx, Nx2), pop=Nx)%>%
  select(-c(Nx2))

#...........................................................
# Blood Pressure data ----
# blood pressure data calculated in file: "Blood pressure.R"
#...........................................................

data.in<-fread(paste0(wd_data,"bp_data6.csv"))%>%rename(location = location_gbd)%>%select(-Year, -Country)
# 
# data.in$salt[data.in$location=="China"]<-4.83*2.54
# length(unique(data.in$location))

#...........................................................
# HTN add scale-up data ----
#...........................................................

inc <- read.csv(paste0(wd_data,"covfxn2.csv"), stringsAsFactors = F)%>%
  select(iso3, location, Year, aroc, p_change, a_change, refwsalt, aspwsalt, reach_base,
         aroc2, p_change2, a_change2, ideal)

bpcats<-c("<120", "120-129", "130-139", 
          "140-149", "150-159", "160-169", 
          "170-179", "180+")

data.in<-merge(bpcats, data.in)%>%rename(bp_cat = x)

data.in <- as.data.table(data.in)
# Fixes location names

name_map <- c(
  "Brunei"                            = "Brunei Darussalam",
  "Cape Verde"                        = "Cabo Verde",
  "Cote d'Ivoire"                     = "Ivory Coast",
  "Czech Republic"                    = "Czechia",
  "Federated States of Micronesia"    = "Micronesia (Federated States of)",
  "Iran"                              = "Iran (Islamic Republic of)",
  "Laos"                              = "Lao People's Democratic Republic",
  "Macedonia"                         = "North Macedonia",
  "Moldova"                           = "Republic of Moldova",
  "South Korea"                       = "Republic of Korea",
  "Swaziland"                         = "Eswatini",
  "Syria"                             = "Syrian Arab Republic",
  "The Bahamas"                       = "Bahamas",
  "The Gambia"                        = "Gambia",
  "Venezuela"                         = "Venezuela (Bolivarian Republic of)",
  "Vietnam"                           = "Viet Nam",
  "North Korea"                       = "Democratic People's Republic of Korea"
)

# 3. update your data.in in place, using fcoalesce() so that
#    any location not in name_map stays unchanged
data.in[, location := fcoalesce(name_map[location], location)]

inc <- as.data.table(inc)
inc[, location := fcoalesce(name_map[location], location)]

unique(data.in$location)
any(is.na(data.in))

locs_data.in <- unique(data.in$location)

#? testing covid.x =0
#b_rates[covid.mx==0, covid.mx:=0]

#rebalance TPs w/ covid such that they sum to less than 1
#especially @ old ages where covid deaths are high
b_rates[,check_well := BG.mx+covid.mx+IR]
b_rates[,check_sick := BG.mx+covid.mx+CF]

#first ensure that background mortality + covid <1
b_rates[check_well>1 | check_sick>1, covid.mx:=ifelse(1-BG.mx<covid.mx, 1-BG.mx, covid.mx)]
#then proportionally reduce rates by check_well
b_rates[check_well>1, covid.mx:= covid.mx - covid.mx*(check_well-1)/(covid.mx+BG.mx+IR)]
b_rates[check_well>1, BG.mx   := BG.mx    - BG.mx*   (check_well-1)/(covid.mx+BG.mx+IR)]
b_rates[check_well>1, IR      := IR       - IR*      (check_well-1)/(covid.mx+BG.mx+IR)]

b_rates[,check_well := BG.mx+covid.mx+IR]
b_rates[check_well>1]

#same process for check_sick
b_rates[check_sick>1, covid.mx:= covid.mx - covid.mx*(check_sick-1)/(covid.mx+BG.mx+CF)]
b_rates[check_sick>1, BG.mx   := BG.mx    - BG.mx*   (check_sick-1)/(covid.mx+BG.mx+CF)]
b_rates[check_sick>1, CF      := CF       - CF*      (check_sick-1)/(covid.mx+BG.mx+CF)]

b_rates[,check_sick := BG.mx+covid.mx+CF]
b_rates[check_sick>1]

#check that no BG.mx.all+covid>1
b_rates[covid.mx+BG.mx.all>1]

#...........................................................
###fxn ----
#...........................................................

repYear<-function(row){
  2017+floor((row-1)/224)
}

data.in<-data.table(data.in%>%select(-age)%>%rename(age=Age.group))
b_rates[, newcases:=0]

##repeat rates for years 2020-2050
rep<-b_rates%>%filter(year==2019)

for (i in 2020:2050){
  b_rates<-bind_rows(b_rates, rep%>%mutate(year=i))
}

# # rename causes to match abbreviated names
# b_rates[,cause:=ifelse(cause=="Ischemic heart disease", "ihd",
#                        ifelse(cause=="Ischemic stroke", "istroke",
#                               ifelse(cause=="Intracerebral hemorrhage", "hstroke",
#                                      ifelse(cause=="Hypertensive heart disease", "hhd",
#                                             ifelse(cause=="Alzheimer's disease and other dementias", "aod",
#                                                    cause)))))]

# Build reverse lookup: full name -> abbreviation
cause_lookup <- setNames(names(cause_map), cause_map)

# Rename in place; unmapped causes keep their original name
b_rates[, cause := fcoalesce(cause_lookup[cause], cause)]

# #...........................................................
# # Adjustments ----
# #...........................................................
# ?? Adjustment of incidence rates and CF 

if(run_adjustment_model == TRUE) {
  
  adjustments <- fread(file = paste0(wd_data,"adjustments2023_age.csv"))
  
  adjustments <- adjustments[,c("location","sex","cause","age_group","IRadjust", "CFadjust"),with=FALSE]
  
  gbd_breaks <- c(seq(20, 95, by = 5), Inf)
  gbd_labels <- c(
    paste0(seq(20, 90, by = 5), "-", seq(24, 94, by = 5)),
    "95+"
  )
  
  # 2) (Optionally) wrap in a helper
  create_gbd_age_group <- function(age) {
    cut(
      age,
      breaks        = gbd_breaks,
      labels        = gbd_labels,
      right         = FALSE,      # [20,25), [25,30), …, [95,Inf)
      include.lowest = TRUE
    )
  }
  
  b_rates[,age_group := create_gbd_age_group(age)]
  # Adjustments for age group
  #b_rates <- merge(b_rates,adjustments,by=c("location","sex","cause"),all.x = T)
  b_rates <- merge(b_rates,adjustments,by=c("location","sex","cause","age_group"),all.x = T)
  
  b_rates[ , age_group:=NULL]
  
  b_rates[!is.na(IRadjust), IR:=IR * IRadjust]
  b_rates[!is.na(CFadjust), CF:=CF * CFadjust]
  
  b_rates[,c("IRadjust", "CFadjust"):=NULL]
  
}

# #...........................................................
# # UNWPP 2024 Pop ----
# #...........................................................
# Adjust pop 20 to unwpp

b_rates<-left_join(b_rates, pop20%>%rename(Nx2=Nx, year=year_id)%>%filter(year>=2017), 
                   by=c("location", "year", "sex", "age"))%>%
  mutate(Nx = ifelse(is.na(Nx2), Nx, Nx2), pop=Nx)%>%
  select(-c(Nx2))

# #...........................................................
# # Covid 2020/2021 ----
# #...........................................................

b_rates[,covid.mx:=NULL]
b_rates <- merge(b_rates,wpp.adj[,c("location_name","year","sex","age","covid.mx"),with=F],
                 by.x=c("location","year","sex","age"),
                 by.y=c("location_name","year","sex","age"),all.x=T)

b_rates[is.na(covid.mx), covid.mx:=0]
b_rates[covid.mx>=1, covid.mx:=0.9]


#...........................................................
# Mortality downward trends ----
#...........................................................

if(run_bgmx_trend == TRUE){
  
  bgmx_fcst <- readRDS(file = paste0(wd_data,"tps_bgmx_forecasted.rds"))
  
  bgmx_fcst[,BG.mx.all:=NULL]
  
  bgmx_fcst <- bgmx_fcst[year>2019,]
  
  bgmx_fcst <- unique(bgmx_fcst,by=c("age","sex","cause","year"))
  
  bgmx_fcst[, cause := fcase(
    cause == "Ischemic heart disease", "ihd",
    cause == "Ischemic stroke", "istroke",
    cause == "Intracerebral hemorrhage", "hstroke",
    cause == "Hypertensive heart disease", "hhd",
    cause == "Alzheimer's disease and other dementias", "aod",
    default = cause
  )]
  
  summary(b_rates$BG.mx)
  
  b_rates <- merge(b_rates,bgmx_fcst,,by=c("age","sex","cause","year"),all.x = T)
  
  b_rates[year>2019 & !is.na(percent_diff),BG.mx:=BG.mx*(1+percent_diff)]
  b_rates[,c("percent_lag","percent_diff"):=NULL]
  
  summary(b_rates$BG.mx)
  
  # All dead envelope
  bgmx_fcst <- readRDS(file = paste0(wd_data,"tps_bgmx_all_forecasted.rds"))
  
  bgmx_fcst <- bgmx_fcst[year>2019,]
  
  bgmx_fcst[,BG.mx.all:=NULL]
  
  bgmx_fcst <- unique(bgmx_fcst,by=c("age","sex","cause","year"))
  
  bgmx_fcst[, cause := fcase(
    cause == "Ischemic heart disease", "ihd",
    cause == "Ischemic stroke", "istroke",
    cause == "Intracerebral hemorrhage", "hstroke",
    cause == "Hypertensive heart disease", "hhd",
    cause == "Alzheimer's disease and other dementias", "aod",
    default = cause
  )]
  
  summary(b_rates$BG.mx.all)
  
  b_rates <- merge(b_rates,bgmx_fcst,by=c("age","sex","cause","year"),all.x = T)
  
  b_rates[year>2019 & !is.na(percent_diff),BG.mx.all:=BG.mx.all*(1+percent_diff)]
  b_rates[,c("percent_lag","percent_diff"):=NULL]
  
  summary(b_rates$BG.mx.all)
}

## Adjusting also CF with downward trend

if(run_CF_trend== TRUE){
  
  if(run_CF_trend_ihme== TRUE){
    
    bgmx_fcst <- readRDS(file = paste0(wd_data,"tps_bgmx_cvd_ihme.rds"))
    
    bgmx_fcst <- bgmx_fcst[year>2019,]
    
    bgmx_fcst <- unique(bgmx_fcst,by=c("cause","year"))
    
    b_rates <- merge(b_rates,bgmx_fcst,by=c("cause","year"),all.x = T)
    
    b_rates[year>2019 & !is.na(percent_diff),CF:=CF*(1+percent_diff)]
    b_rates[,c("percent_diff"):=NULL]
    
  }else{
    
    # All dead envelope
    #bgmx_fcst <- readRDS(file = paste0(wd_data,"tps_bgmx_all_forecasted.rds"))
    bgmx_fcst <- readRDS(file = paste0(wd_data,"tps_bgmx_cvd_forecasted.rds"))
    
    bgmx_fcst <- bgmx_fcst[year>2019,]
    
    bgmx_fcst[,BG.mx.all:=NULL]
    
    bgmx_fcst <- unique(bgmx_fcst,by=c("age","sex","cause","year"))
    
    bgmx_fcst[, cause := fcase(
      cause == "Ischemic heart disease", "ihd",
      cause == "Ischemic stroke", "istroke",
      cause == "Intracerebral hemorrhage", "hstroke",
      cause == "Hypertensive heart disease", "hhd",
      cause == "Alzheimer's disease and other dementias", "aod",
      default = cause
    )]
    
    
    b_rates <- merge(b_rates,bgmx_fcst,by=c("age","sex","cause","year"),all.x = T)
    
    if(run_CF_trend_80 == TRUE){
      b_rates[year>2019 & !is.na(percent_diff),CF:=CF*(1+percent_diff*0.8)]
    }else{
      b_rates[year>2019 & !is.na(percent_diff),CF:=CF*(1+percent_diff)]
    }
    
    b_rates[,c("percent_lag","percent_diff"):=NULL]
    
  }
  
}

# Clean up environment
rm("adjustments","bgmx_fcst","dt_pop_unwpp","wpp.adj","rep","pop20")
