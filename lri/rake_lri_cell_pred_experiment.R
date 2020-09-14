#load packages
#load all the functions
repo <- '<<<< FILEPATH REDACTED >>>>'
funks = list.files(paste0(repo,'mbg_central'), pattern = 'functions',full.names = T)
for(fff in funks){
  source(fff)
}

## Step 2: Rake all new cell preds to 2017 GBD values for prevalence

#load in required functions
source('<<<< FILEPATH REDACTED >>>>')
source('<<<< FILEPATH REDACTED >>>>')

library(plyr)
library(dplyr)
library(mgcv)
library(raster)
library(sf)
library(fasterize)
library(rgdal)
library(rgeos)
library(data.table)
library(INLA)
library(seegSDM)
library(seegMBG)
library(dismo)
library(gbm)
library(foreign)
library(parallel)
library(doParallel)
library(grid)
library(gridExtra)
library(pacman)
library(gtools)
library(glmnet)
library(ggplot2)
library(RMySQL)
library(tictoc)
library(magrittr)

#recover variables from qsub call
raketo <- as.character(commandArgs()[6])
reg <- as.character(commandArgs()[7])
outputdir <- as.character(commandArgs()[8])
makeholdouts <- as.character(commandArgs()[9])
year_list <- c(2000:2017)
indicator <- as.character(commandArgs()[11])
modeling_shapefile_version <- as.character(commandArgs()[12])
raking_shapefile_version <- as.character(commandArgs()[13])
fun_tol <- as.numeric(commandArgs()[14])

print(raketo)
print(reg)
print(outputdir)
print(makeholdouts)
print(year_list)
print(indicator)
print(modeling_shapefile_version)
print(raking_shapefile_version)

if (raketo == 'incidence'){
  measure_id <- 6
} else if (raketo == 'prevalence'){
  measure_id <- 5
} else if (raketo == 'mortality'){
  measure_id <- 1
}

cp_dir <- outputdir

if (makeholdouts == F){
  holdout = 0
} else if (makeholdouts == T){
  holdout = 0
  message('Note! Only raking holdout = 0')
}

# pull gbd results
# GBD Location ID Table
locs <- read.csv('<<<< FILEPATH REDACTED >>>>')

# Subset locations to National and Subnational estimates
country <- subset(locs, level >=3)$location_id

gbd_lri <- get_outputs(topic= 'cause',
                       cause_id=322,
                       location_id=country,
                       year_id=year_list,
                       age_group_id=1,
                       gbd_round_id=5,
                       metric_id=3,
                       measure_id = measure_id,
                       sex_id = 3,
                       version='latest')%>%
  
  as.data.frame() %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)), newdata = data.frame(year_id = year_list)))))

gbd_lri$val <- ifelse(is.na(gbd_lri$val), gbd_lri$inter_val, gbd_lri$val)
gbd_lri <- gbd_lri[complete.cases(gbd_lri),]

#rename columns in gbd_lri to be read correctly by rake_cell pred
#gaul_code column named 'nqme'
#year column named 'year'
#value column named 'mean'
gbd_lri$name <- gbd_lri$location_id 
gbd_lri$year <- gbd_lri$year_id
gbd_lri$mean <- gbd_lri$val

#load in cell pred object
cell_pred <- readRDS(paste0(cp_dir, indicator,'_cell_draws_eb_bin', '0_', reg,'_', holdout, '.rds'))
gc()

#rake to gbd
raked_cell_pred <- rake_cell_pred(cell_pred = cell_pred,
                                  rake_to = gbd_lri,
                                  reg = reg,
                                  year_list = year_list,
                                  pop_measure = 'a0004t',
                                  rake_method = 'logit',
                                  rake_subnational = T,
                                  crosswalk = T,
                                  shapefile_path = get_admin_shapefile(admin_level = 0, raking = T, version = 'current'),
                                  field = 'loc_id',
                                  zero_heuristic = F,
                                  approx_0_1 = F,
                                  simple_raster = NULL,
                                  simple_polygon = NULL,
                                  pop_raster = NULL,
                                  modeling_shapefile_version = modeling_shapefile_version,
                                  raking_shapefile_version = raking_shapefile_version,
                                  if_no_gbd = 'return_unraked',
                                  FunTol = fun_tol)

#save as rds file 
saveRDS(raked_cell_pred, paste0(cp_dir, indicator,'_raked_cell_draws_eb_bin', '0_', reg,'_', holdout,'_',raketo, '_', fun_tol, '.RData'))
