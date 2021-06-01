#########
# Fractional raking and aggregation implementation
#########

## SETUP ################################################################################
rm(list=ls())
stratum=as.character(commandArgs()[4])
run_date=as.character(commandArgs()[5])
indicator=as.character(commandArgs()[6])
indicator_group=as.character(commandArgs()[7])
geos_node =as.logical(commandArgs()[8])
modeling_shapefile_version=as.character(commandArgs()[9])
raking_shapefile_version=as.character(commandArgs()[9])
interval_mo <- 12
age <- 0
holdout <- 0

rake_list <- c("unraked", "raked", "raked_c")

# Define directories
main_dir <- paste0('<<<<FILEPATH REDACTED>>>>>', indicator_group, '/', indicator, '/output/', run_date, '/')
temp_dir <- paste0(main_dir, "temp_post_est/")

# Load objects from convenience temp file
load(paste0(temp_dir, "post_est_temp_objs.RData"))
summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list)
# Print some settings to console
message('indicator')
message(indicator)
message('indicator_group')
message(indicator_group)
message('run_date')
message(run_date)
message('stratum')
message(stratum)
message('pop_measure')
message(pop_measure)
message('age_group')
message(age_group)
message('sex_id')
message(sex_id)
message('shapefile_version')
message(shapefile_version)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

# For now just assign stratum to reg (will need to modify the below for strata beyond reg)
reg <- stratum
region <- stratum
# do we really need this again?
sharedir       <- sprintf('<<<<FILEPATH REDACTED>>>>>%s/%s',indicator_group,indicator)
commondir      <- sprintf('<<<<FILEPATH REDACTED>>>>>common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Set repo
core_repo  <- paste0("<<<<FILEPATH REDACTED>>>>>", "/lbd_core/")
indic_repo <- paste0("<<<<FILEPATH REDACTED>>>>>", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)


pathaddin<-paste0('_bin',age,'_',reg,'_',holdout)

require(rgdal)

# loading some bespoke functions developed for this process
source(paste0(indic_repo, "mbg/functions/fractional_raking_functions_hiv.r"))

#####################################################################
# Loading in the raw draws data
message("Loading Cell Draws that are unraked")
# Load cell draws
message('Loading Data...')
load(paste0(main_dir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))

# setting a reference for the number of dreaws
ndraws = ncol(cell_pred)

#####################################################################
# Getting the simple polygon and simple raster objects for this region alone
message("Getting the spatial objects associated with this region.")

#load regional simple raster
simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = shapefile_version),
                                           buffer = 0.4, subset_only = FALSE,
                                           shapefile_version = shapefile_version)


subset_shape   <- simple_polygon_list[[1]]
simple_polygon = simple_polygon_list[[2]]
message("Building simple raster from subset_shape")
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
rm(raster_list,simple_polygon_list);gc()

message("All done loading spatial template files (subset_shape,simple_polygon,simple_raster,pop_raster, etc)")

#####################################################################
# Determining a list of the valid pixel indices based on the simple raster template

pixel_id <- seegSDM:::notMissingIdx(simple_raster)

#####################################################################
#load the cell id to admin units link

link_table <- get_link_table_hiv(simple_raster, shapefile_version = shapefile_version)

#####################################################################
# collect and load the population data
covdt <- load_populations_cov(reg, pop_measure, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo)

#####################################################################
# Prepping the cell_pred and link table to be linked
link <- prep_link_table(link_table = link_table,
                        simple_raster = simple_raster,
                        pixel_id = pixel_id)

cell_ids <- link_table[[2]]

# getting the connector for sub-national raking
connector <- get_gbd_locs(rake_subnational,
                          reg, shapefile_version=shapefile_version)

# merge the connector on to the link table

link <- sub_nat_link_merge(rake_subnational,
                           link,
                           connector)

#set cell pred as a data table, and rename things
cell_pred <- prep_cell_pred(cell_pred=cell_pred,
                            cell_ids = cell_ids,
                            pixel_id = pixel_id,
                            covdt = covdt)

#merge on the link
cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian=TRUE)

# space
link <- NULL

message("raking population")
#convert to fractional population
cell_pred = cell_pred[,pop:= pop * area_fraction]

#NA out population where the pixel value is NA (to prevent weirdness with denominators)
cell_pred = cell_pred[is.na(V1), pop := NA]

scalars <- calculate_pop_scalars(cell_pred = cell_pred,
                                 age_group = age_group,
                                 connector = connector,
                                 sex = sex_id,
                                 sharedir = sharedir,
                                 run_date = run_date,
                                 indicator = indicator,
                                 stratum = stratum)
# add back to the cell_pred as a population rf
cell_pred <- merge(cell_pred, scalars, by = c("location_id", "year"))

# rake fractional populations
cell_pred$pop_raked <- 0
cell_pred = cell_pred[,pop_raked := pop * pop_scalar]
cell_pred$pop <- NULL

message("raking rates")
# Calculate Fractional Raking Factors
fractional_rf <- calculate_fractional_rfs(ndraws=ndraws,
                                          cell_pred = cell_pred,
                                          gbd=gbd,
                                          sharedir = sharedir,
                                          run_date = run_date,
                                          indicator = indicator,
                                          stratum = stratum)
#########
# Clearing up after raking calculations, preparing to apply the raking factors and then aggregate
# This is nessecary because to get back to a rates cell pred we cannot simply divide by pop the
# same as we multiplied by pop because of 0 pop locations
#########

## SETUP ################################################################################
rm(list=ls())
stratum=as.character(commandArgs()[4])
run_date=as.character(commandArgs()[5])
indicator=as.character(commandArgs()[6])
indicator_group=as.character(commandArgs()[7])
geos_node =as.logical(commandArgs()[8])
modeling_shapefile_version=as.character(commandArgs()[9])
raking_shapefile_version=as.character(commandArgs()[9])
interval_mo <- 12
age <- 0
holdout <- 0

rake_list <- c("unraked", "raked", "raked_c")

# Define directories
main_dir <- paste0('<<<<FILEPATH REDACTED>>>>>', indicator_group, '/', indicator, '/output/', run_date, '/')
temp_dir <- paste0(main_dir, "temp_post_est/")

# Load objects from convenience temp file
load(paste0(temp_dir, "post_est_temp_objs.RData"))
summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list)
# Print some settings to console
message('indicator')
message(indicator)
message('indicator_group')
message(indicator_group)
message('run_date')
message(run_date)
message('stratum')
message(stratum)
message('pop_measure')
message(pop_measure)
message('age_group')
message(age_group)
message('sex_id')
message(sex_id)
message('shapefile_version')
message(shapefile_version)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

# For now just assign stratum to reg (will need to modify the below for strata beyond reg)
reg <- stratum
region <- stratum
# do we really need this again?
sharedir       <- sprintf('<<<<FILEPATH REDACTED>>>>>%s/%s',indicator_group,indicator)
commondir      <- sprintf('<<<<FILEPATH REDACTED>>>>>common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Set repo
core_repo  <- paste0("<<<<FILEPATH REDACTED>>>>>", "/lbd_core/")
indic_repo <- paste0("<<<<FILEPATH REDACTED>>>>>", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)


pathaddin<-paste0('_bin',age,'_',reg,'_',holdout)

require(rgdal)

# loading some bespoke functions developed for this process
source(paste0(indic_repo, "mbg/functions/fractional_raking_functions_hiv.r"))

#####################################################################
# Loading in the raw draws data
message("Loading Cell Draws that are unraked")
# Load cell draws
message('Loading Data...')
load(paste0(main_dir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))

# setting a reference for the number of dreaws
ndraws = ncol(cell_pred)

#####################################################################
# Getting the simple polygon and simple raster objects for this region alone
message("Getting the spatial objects associated with this region.")

#load regional simple raster
simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = shapefile_version),
                                           shapefile_version = shapefile_version,
                                           buffer = 0.4, subset_only = FALSE)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon = simple_polygon_list[[2]]
message("Building simple raster from subset_shape")
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
rm(raster_list,simple_polygon_list);gc()

message("All done loading spatial template files (subset_shape,simple_polygon,simple_raster,pop_raster, etc)")

#####################################################################
# Determining a list of the valid pixel indices based on the simple raster template

pixel_id <- seegSDM:::notMissingIdx(simple_raster)

#####################################################################
#load the cell id to admin units link

link_table <- get_link_table_hiv(simple_raster, shapefile_version = shapefile_version)

#####################################################################
# collect and load the population data
covdt <- load_populations_cov(reg, pop_measure, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo)

#####################################################################
# Prepping the cell_pred and link table to be linked and then merging them
link <- prep_link_table(link_table = link_table,
                        simple_raster = simple_raster,
                        pixel_id = pixel_id)

cell_ids <- link_table[[2]]

# getting the connector for sub-national raking
connector <- get_gbd_locs(rake_subnational,
                          reg, shapefile_version=shapefile_version)

# merge the connector on to the link table

link <- sub_nat_link_merge(rake_subnational,
                           link,
                           connector)

#set cell pred as a data table, and rename things
cell_pred <- prep_cell_pred(cell_pred=cell_pred,
                            cell_ids = cell_ids,
                            pixel_id = pixel_id,
                            covdt = covdt)

#merge on the link
cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian=TRUE)


############################################################
#adding the raking factors and scaling the populations

message("adding raking factors")
#convert to fractional population
cell_pred = cell_pred[,pop:= pop * area_fraction]

scalars <- read.csv(file = paste0(sharedir, "/output/", run_date, "/", indicator, "_", stratum, "_pop_rf.csv"))

# add back to the cell_pred as a population rf
cell_pred <- merge(cell_pred, scalars, by = c("location_id", "year"))

## load the raking factors
fractional_rf <- read.csv(file = paste0(sharedir, "/output/", run_date, "/", indicator, "_", stratum, "_rf.csv"))

## merge them onto the cell_pred
cell_pred <- merge(cell_pred, fractional_rf, by = c("location_id", "year"))


############################################################
# creating a raked rates cell_pred (this happens first becasue once we go to counts in the cell_pred we can't do back with out loading a fresh copy)
message("creating a raked rates cell_pred")
# rake rates
overs = paste0('V',1:ndraws)
cell_pred = cell_pred[, (overs) := lapply(overs, function(x) get(x) * rf )]

# multiply the cell_pred by the area fraction for the dedupe function (so that each cell will add to 1 and the constituent rates are weighted by area)
cell_pred = cell_pred[, (overs) := lapply(overs, function(x) get(x) * area_fraction) ]

raked_cell_pred <- dedupe_linked_cell_pred(cell_pred, overs)

# Save raked rates cell_pred object
save(raked_cell_pred, file = paste0(sharedir, "/output/", run_date, "/",
                                       indicator, "_raked_cell_draws_eb_bin0_", reg, "_0.RData" ))
# SPACE!!!!!
raked_cell_pred <- NULL

# un do the area fraction (so that you can use this cell pred again)
overs = paste0('V',1:ndraws)
cell_pred = cell_pred[, (overs) := lapply(overs, function(x) get(x) / area_fraction) ]

############################################################
# creating a raked counts cell_pred
message("creating a raked counts cell_pred")
# rake fractional populations
cell_pred$pop_raked <- 0
cell_pred = cell_pred[,pop_raked := pop * pop_scalar]
cell_pred$pop <- NULL

message("converting from prevalence to counts")
#set the variables to aggregate
overs = paste0('V',1:ndraws)
cell_pred = cell_pred[, (overs) := lapply(overs, function(x) get(x) * pop_raked )]

#NA out population where the pixel value is NA (to prevent weirdness with denominators)
cell_pred = cell_pred[is.na(V1), pop_raked := NA]
raked_cell_pred_c <- dedupe_linked_cell_pred(cell_pred, overs)
save(raked_cell_pred_c, file = paste0(sharedir, "/output/", run_date, "/",
                                       indicator, "_raked_c_cell_draws_eb_bin0_", reg, "_0.RData" ))
# SPACE
raked_cell_pred_c <- NULL

############################################################
# creating a raked counts aggregations
message("creating a raked counts aggregations")

#do calculations!
admin_2 = cell_pred[, lapply(c(overs,'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c('year', 'ADM2_CODE', 'ADM0_CODE')]
admin_1 = cell_pred[, lapply(c(overs,'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c('year','ADM1_CODE', 'ADM0_CODE')]
admin_0 = cell_pred[, lapply(c(overs,'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c('year','ADM0_CODE')]


setnames(admin_2, grep('V[0-9]', names(admin_2),value = T),c(overs, 'pop_raked'))
setnames(admin_1, grep('V[0-9]', names(admin_1),value = T),c(overs, 'pop_raked'))
setnames(admin_0, grep('V[0-9]', names(admin_0),value = T),c(overs, 'pop_raked'))

#create the spatial hierarchy
sp_hierarchy_list = unique(link[ADM0_CODE %in% unique(admin_0[,ADM0_CODE]),.(ADM0_CODE, ADM1_CODE, ADM2_CODE, ADM0_NAME, ADM1_NAME, ADM2_NAME, region)])
sp_hierarchy_list[ ,region:= reg]

# cleaning raked admin draws in count space
admin_0$pop <- admin_0$pop_raked
admin_0$pop_raked <- NULL
admin_1$ADM0_CODE <- NULL
admin_1$pop <- admin_1$pop_raked
admin_1$pop_raked <- NULL
admin_2$ADM0_CODE <- NULL
admin_2$pop <- admin_2$pop_raked
admin_2$pop_raked <- NULL

## save raked counts aggregations
raked <- TRUE
save(admin_0,admin_1,admin_2,sp_hierarchy_list,
     file=paste0(main_dir, indicator, "_", ifelse(raked, "raked", "unraked"),"_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))

############################################################
# creating raked rates aggregations (you can work back at the admin level because there are no admin's with pop = 0)
message("creating a raked rates aggregations")

#convert back to rates/prevalence
overs = paste0('V',1:ndraws)
admin_0 = admin_0[, (overs) := lapply(overs, function(x) get(x)/pop) ]
admin_1 = admin_1[, (overs) := lapply(overs, function(x) get(x)/pop) ]
admin_2 = admin_2[, (overs) := lapply(overs, function(x) get(x)/pop) ]

## save raked rates aggregations
raked <- TRUE
save(admin_0,admin_1,admin_2,sp_hierarchy_list,
     file=paste0(main_dir, indicator, "_", ifelse(raked, "raked", "unraked"), "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))


############################################################
# creating unraked counts aggregations (can be done two ways, un raking the counts cell_pred or reloading the unraked cell_pred and multiplying by the pop.  I chose this to avoid reloading cell_preds)
message("creating a unraked counts aggregations")

#unrake all draws
overs = paste0('V',1:ndraws)
cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x)/rf) ]

## Create unraked counts agregation
admin_2 = cell_pred[, lapply(c(overs,'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c('year', 'ADM2_CODE', 'ADM0_CODE')]
admin_1 = cell_pred[, lapply(c(overs,'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c('year','ADM1_CODE', 'ADM0_CODE')]
admin_0 = cell_pred[, lapply(c(overs,'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c('year','ADM0_CODE')]

setnames(admin_2, grep('V[0-9]', names(admin_2),value = T),c(overs, 'pop_raked'))
setnames(admin_1, grep('V[0-9]', names(admin_1),value = T),c(overs, 'pop_raked'))
setnames(admin_0, grep('V[0-9]', names(admin_0),value = T),c(overs, 'pop_raked'))

admin_0$pop <- admin_0$pop_raked
admin_0$pop_raked <- NULL
admin_1$pop <- admin_1$pop_raked
admin_1$pop_raked <- NULL
admin_1$ADM0_CODE <- NULL
admin_2$pop <- admin_2$pop_raked
admin_2$pop_raked <- NULL
admin_2$ADM0_CODE <- NULL

# SPACE!!
cell_pred <- NULL

## save unraked counts aggregations
raked <- FALSE
save(admin_0,admin_1,admin_2,sp_hierarchy_list,
     file=paste0(main_dir, indicator, "_", ifelse(raked, "raked", "unraked"),"_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))


############################################################
# creating unraked rates aggregations
message("creating a unraked rates aggregations")

#convert back to rates
overs = paste0('V',1:ndraws)
admin_0 = admin_0[, (overs) := lapply(overs, function(x) get(x)/pop) ]
admin_1 = admin_1[, (overs) := lapply(overs, function(x) get(x)/pop) ]
admin_2 = admin_2[, (overs) := lapply(overs, function(x) get(x)/pop) ]

## save unraked rates aggregations
raked <- FALSE
save(admin_0,admin_1,admin_2,sp_hierarchy_list,
     file=paste0(main_dir, indicator, "_", ifelse(raked, "raked", "unraked"), "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))



############################################################
## save raster objects

save_cell_pred_summary <- function(summstat, raked, ...) {
  message(paste0('Making unraked summmary raster for: ',summstat, " (", raked, ")"))
  if (raked=="unraked") {
    load(paste0(main_dir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))
    cpred <- "cell_pred"
  } else if (raked=="raked") {
    load(paste0(sharedir, "/output/", run_date, "/", indicator, "_raked_cell_draws_eb_bin0_", reg, "_0.RData" ))
    cpred <- "raked_cell_pred"
  } else if (raked=="raked_c") {
    load(paste0(sharedir, "/output/", run_date, "/", indicator, "_raked_c_cell_draws_eb_bin0_", reg, "_0.RData" ))
    cpred <- "raked_cell_pred_c"
  }
  ras <- make_cell_pred_summary(
    draw_level_cell_pred = get(cpred),
    mask                 = simple_raster,
    return_as_raster     = TRUE,
    summary_stat         = summstat,
    ...)
  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ifelse(raked == 'raked_c', '_raked_c', '')),'_',summstat,'_raster'))
}


lapply(1:nrow(summ_list), function(i) {
  summstat <- as.character(summ_list[i, 1])
  raked <- as.character(summ_list[i, 2])
  save_cell_pred_summary(summstat, raked)
})


## Can't pass additional params in the above framework, so will code by hand here
for (r in rake_list) {
  if ("p_below" %in% summstats) {
    save_cell_pred_summary(summstat = "p_below",
                           raked = r,
                           value = 0.8,
                           equal_to = F)
  }
}

# Write a file to mark done
output_dir <- paste0('<<<<FILEPATH REDACTED>>>>>', indicator_group, '/', indicator, '/output/', run_date)
pathaddin <- paste0('_bin0_',reg,'_0') # To allow us to use waitformodelstofinish()
write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

# All done
message(paste0("Done with post-estimation for ", stratum))

