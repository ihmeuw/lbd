## Step 2: Rake all new cell preds to 2017 GBD values for INCIDENCE

#load in required functions
source('<<<< FILEPATH REDACTED >>>>') #raking functions
source('<<<< FILEPATH REDACTED >>>>') #prep fucntions
commondir      <- '<<<< FILEPATH REDACTED >>>>'
core_repo = '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>', header = FALSE)))
source('<<<< FILEPATH REDACTED >>>>') #setup
mbg_setup(package_list = package_list, repos = core_repo)
library(mgcv)
library(raster)
library(sf)
library(fasterize)
library(dplyr)

#prepare to load in raked cell preds
cp_dir <- '<<<< FILEPATH REDACTED >>>>'
tag <- '<<<< FILEPATH REDACTED >>>>'

#get gbd incidence estimates
#load function
source('<<<< FILEPATH REDACTED >>>>') #get outputs

# GBD Location ID Table
locs <- read.csv('<<<< FILEPATH REDACTED >>>>')

# Subset locations to National and Subnational estimates
country <- subset(locs, level >=3)$location_id

gbd_lri_inc <- get_outputs(topic= 'cause',
                           cause_id=322,
                           location_id=country,
                           year_id=2000:2017,
                           age_group_id=1,
                           gbd_round_id=5,
                           metric_id=3,
                           measure_id = 6,
                           sex_id = 3,
                           version='latest')%>%
  
  as.data.frame() %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)),
                                            newdata = data.frame(year_id = 2000:2017))))) 

gbd_lri_inc$val <- ifelse(is.na(gbd_lri_inc$val), gbd_lri_inc$inter_val, gbd_lri_inc$val)
gbd_lri_inc <- gbd_lri_inc[complete.cases(gbd_lri_inc),]

#rename columns in gbd_lri_inc to be read correctly by rake_cell pred
#gaul_code column named "nqme"
#year column named "year"
#value column named "mean"
gbd_lri_inc$name <- gbd_lri_inc$location_id 
gbd_lri_inc$year <- gbd_lri_inc$year_id
gbd_lri_inc$mean <- gbd_lri_inc$val


##cssa
cssa_cell_pred <- readRDS('<<<< FILEPATH REDACTED >>>>')
gc()

#rake cssa cell pred
cssa_raked_cell_pred <- rake_cell_pred(cell_pred = cssa_cell_pred,
                                       rake_to = gbd_lri_inc,
                                       reg = 'cssa',
                                       year_list = c(2000:2017),
                                       pop_measure = 'a0004t',
                                       rake_method = "logit",
                                       rake_subnational = F,
                                       crosswalk = T,
                                       shapefile_path = get_admin_shapefile(admin_level = 0, raking = T, version = 'current'),
                                       field = "loc_id",
                                       zero_heuristic = F,
                                       approx_0_1 = F,
                                       simple_raster = NULL,
                                       simple_polygon = NULL,
                                       pop_raster = NULL,
                                       modeling_shapefile_version = '2018_08_28',
                                       raking_shapefile_version = '2018_12_04',
                                       if_no_gbd = 'return_unraked')
#save this raked cell pred
saveRDS(cssa_raked_cell_pred, '<<<< FILEPATH REDACTED >>>>')

#remove the cell pred given memory concerns
rm(cssa_raked_cell_pred)
rm(cssa_cell_pred)
gc()

##sssa
sssa_cell_pred <- readRDS('<<<< FILEPATH REDACTED >>>>')
gc()

#rake sssa cell pred
sssa_raked_cell_pred <- rake_cell_pred(cell_pred = sssa_cell_pred,
                                       rake_to = gbd_lri_inc,
                                       reg = 'sssa',
                                       year_list = c(2000:2017),
                                       pop_measure = 'a0004t',
                                       rake_method = "logit",
                                       rake_subnational = F,
                                       crosswalk = T,
                                       shapefile_path = get_admin_shapefile(admin_level = 0, raking = T, version = 'current'),
                                       field = "loc_id",
                                       zero_heuristic = F,
                                       approx_0_1 = F,
                                       simple_raster = NULL,
                                       simple_polygon = NULL,
                                       pop_raster = NULL,
                                       modeling_shapefile_version = '2018_08_28',
                                       raking_shapefile_version = '2018_12_04',
                                       if_no_gbd = 'return_unraked')
gc()

#save this raked cell pred
saveRDS(sssa_raked_cell_pred, '<<<< FILEPATH REDACTED >>>>')

#remove the cell pred given memory concerns
rm(sssa_raked_cell_pred)
rm(sssa_cell_pred)
gc()


##wssa
wssa_cell_pred <- readRDS('<<<< FILEPATH REDACTED >>>>')
gc()

#rake wssa cell_pred
wssa_raked_cell_pred <- rake_cell_pred(cell_pred = wssa_cell_pred,
                                       rake_to = gbd_lri_inc,
                                       reg = 'wssa',
                                       year_list = c(2000:2017),
                                       pop_measure = 'a0004t',
                                       rake_method = "logit",
                                       rake_subnational = F,
                                       crosswalk = T,
                                       shapefile_path = get_admin_shapefile(admin_level = 0, raking = T, version = 'current'),
                                       field = "loc_id",
                                       zero_heuristic = F,
                                       approx_0_1 = F,
                                       simple_raster = NULL,
                                       simple_polygon = NULL,
                                       pop_raster = NULL,
                                       modeling_shapefile_version = '2018_08_28',
                                       raking_shapefile_version = '2018_12_04',
                                       if_no_gbd = 'return_unraked')
gc()

#save this raked cell pred
saveRDS(wssa_raked_cell_pred, '<<<< FILEPATH REDACTED >>>>')

#remove the cell pred given memory concerns
rm(wssa_raked_cell_pred)
rm(wssa_cell_pred)
gc()

##name
name_cell_pred <- readRDS('<<<< FILEPATH REDACTED >>>>')
gc()

#rake cssa cell pred
name_raked_cell_pred <- rake_cell_pred(cell_pred = name_cell_pred,
                                       rake_to = gbd_lri_inc,
                                       reg = 'name',
                                       year_list = c(2000:2017),
                                       pop_measure = 'a0004t',
                                       rake_method = "logit",
                                       rake_subnational = F,
                                       crosswalk = T,
                                       shapefile_path = get_admin_shapefile(admin_level = 0, raking = T, version = 'current'),
                                       field = "loc_id",
                                       zero_heuristic = F,
                                       approx_0_1 = F,
                                       simple_raster = NULL,
                                       simple_polygon = NULL,
                                       pop_raster = NULL,
                                       modmodeling_shapefile_version = '2018_08_28',
                                       raking_shapefile_version = '2018_12_04',
                                       if_no_gbd = 'return_unraked')
gc()

#save this raked cell pred
saveRDS(name_raked_cell_pred, '<<<< FILEPATH REDACTED >>>>')

#remove the cell pred given memory concerns
rm(name_raked_cell_pred)
rm(name_cell_pred)
gc()

##essa
essa_cell_pred <- readRDS('<<<< FILEPATH REDACTED >>>>')
gc()

#rake essa cell pred
essa_raked_cell_pred <- rake_cell_pred(cell_pred = essa_cell_pred,
                                       rake_to = gbd_lri_inc,
                                       reg = 'essa',
                                       year_list = c(2000:2017),
                                       pop_measure = 'a0004t',
                                       rake_method = "logit",
                                       rake_subnational = F,
                                       crosswalk = T,
                                       shapefile_path = get_admin_shapefile(admin_level = 0, raking = T, version = 'current'),
                                       field = "loc_id",
                                       zero_heuristic = F,
                                       approx_0_1 = F,
                                       simple_raster = NULL,
                                       simple_polygon = NULL,
                                       pop_raster = NULL,
                                       modeling_shapefile_version = '2018_08_28',
                                       raking_shapefile_version = '2018_08_28',
                                       if_no_gbd = 'return_unraked')
gc()

#save this raked cell pred
saveRDS(essa_raked_cell_pred, '<<<< FILEPATH REDACTED >>>>')