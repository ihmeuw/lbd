####################################################################################################
## Launch HIV prevalence (hiv_prev_disagg) models for Africa
####################################################################################################

###############################################################################
## SETUP
###############################################################################

## Clear environment
rm(list=ls())

## Set indicator
indicator_group <- 'hiv'
indicator       <- 'hiv_prev_disagg'

## Set repo
core_repo  <- paste0(<<<< FILEPATH REDACTED >>>>"/lbd_core/")
if (!dir.exists(core_repo)) core_repo <- "<<<< FILEPATH REDACTED >>>>/lbd_core/"
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>"/lbd_hiv/")

setwd(core_repo)

## Load libraries and  MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = core_repo)
source(paste0(indic_repo, 'mbg/', indicator, '/3_functions/waitformodelstofinish.r'))

## Load passed arguments
config_name <- commandArgs()[4]
covs_name <- commandArgs()[5]
run_date <- commandArgs()[6]
test<- F

print(config_name)
print(covs_name)
print(run_date)

## Record git status info
outputdir <- paste(<<<< FILEPATH REDACTED >>>>, run_date, '', sep='/')

tmb_agg_run_date<-NULL

## Load and check config
config <- set_up_config(repo            = indic_repo,
                        core_repo       = core_repo,
                        indicator       = "",
                        indicator_group = "",
                        config_file     = paste0(indic_repo, 'mbg/', indicator, '/2_modeling/', config_name,   '.csv'),
                        covs_file       = paste0(indic_repo, 'mbg/', indicator, '/2_modeling/', covs_name,   '.csv'),
                        run_tests       = FALSE)

if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (length(summstats) == 1 & grepl("c(", summstats, fixed =T)) summstats <- eval(parse(text=summstats))
if (class(z_list) == "character") z_list <- eval(parse(text=z_list))
if (exists('sex_list')) sex_list <- eval(parse(text=sex_list))

## Make sure pc prior options are correct class
if (is.character(pc_prior)) pc_prior <- as.logical(pc_prior)
if (pc_prior & is.character(sigma_upper_tail) | pc_prior & is.character(range_lower_tail)){
  range_lower_tail <- as.numeric(range_lower_tail)
  sigma_upper_tail <- as.numeric(sigma_upper_tail)
}

##Set cluster environment for fair cluster
Sys.setenv(SGE_ENV=SGE_ENV)

## Get cluster project-related arguments
stopifnot(queue %in% c('long.q', 'geospatial.q'))
if (queue == 'geospatial.q') {
  project <- 'proj_geo_nodes_hiv'
  use_geos_nodes <- T
}
if (queue == 'long.q') {
  project <- 'proj_geospatial_hiv'
  use_geos_nodes <- F
}

###############################################################################
## Make Holdouts
###############################################################################
if(!exists('overwrite_holdout')) overwrite_holdout = F
if(!exists('create_holdout_files')) create_holdout_files = F
if(as.logical(makeholdouts) & as.logical(create_holdout_files)){
  'creating holdout files'
  for(reg in Regions){
    #name the holdout file
    holdout_file = paste0(outputdir, '/holdout_file_', n_ho_folds, 'f_', reg, '.rda')
  # load the full input data
  df <- load_input_data(indicator   = indicator,
                        removeyemen = TRUE,
                        years       = yearload,
                        yl          = year_list,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share),
                        region = reg)

# restrict to just age-specific points
df_point <- df[df$point==1 & !(is.na(agebin)) & type!='ANC', ]

# add in location information
df_point <- merge_with_ihme_loc(df_point, shapefile_version = modeling_shapefile_version)

#Set target number of holdout observations
ho_ts <- c(cssa = 425, 'essa_sdn-COM' = 500, sssa = 425, 'wssa-CPV-STP-MRT' = 500)[reg]

# make 5 qt folds
                         stratum_ho <- make_folds(data       = df_point,
                                                  n_folds    = as.numeric(n_ho_folds),
                                                  spat_strat = 'qt',
                                                  temp_strat = 'prop',
                                                  strat_cols = 'region',
                                                  ts         = as.numeric(ho_ts),
                                                  mb         = as.numeric(ho_mb),
                                                  save.file = paste0(<<<< FILEPATH REDACTED >>>, '/output/',
                                                                     run_date, '/stratum_', reg, '.rds'))
                         
                         df.points.tmp <- stratum_ho[[1]]
                         save(df.points.tmp, file=holdout_file)

  }
}

###############################################################################
## Launch Parallel Script
###############################################################################

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if (makeholdouts) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)
loopvars$submitted <- F

if (setequal(c('essa_sdn-COM', 'wssa-CPV-STP-MRT', 'cssa', 'sssa'), Regions)) {
  loopvars$Var1 <- factor(loopvars$Var1, levels = c('essa_sdn-COM', 'wssa-CPV-STP-MRT', 'cssa', 'sssa')) # temporary hack to make sure bigger models get submitted first
  loopvars <- loopvars[order(loopvars$Var1, loopvars$Var3), ]
  loopvars$Var1 <- as.character(loopvars$Var1)
}

## Loop over them, save images and submit qsubs
for (i in 1:nrow(loopvars)) {

  message(paste(loopvars[i, 2], as.character(loopvars[i, 1]), loopvars[i, 3]))

  # set memory
  memory <- c('cssa' = 30, 'essa_sdn-COM' = 60, 'sssa' = 50, 'wssa-CPV-STP-MRT' = 50)[as.character(loopvars[i, 1])]

  # make a qsub string
  qsub <- make_qsub_share(age            = loopvars[i, 2],
                          reg            = as.character(loopvars[i, 1]),
                          holdout        = loopvars[i, 3],
                          test           = as.logical(test),
                          indic          = indicator,
                          saveimage      = TRUE,
                          memory         = memory,
                          cores          = 5,
                          proj           = 'proj_geo_nodes',
                          geo_nodes      = FALSE,
                          use_c2_nodes   = FALSE,
                          singularity    = '<<<< FILEPATH REDACTED >>>>/lbd_full_20200128.simg',
                          singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1),
                          queue          = 'geospatial.q',
                          run_time       = '05:00:00:00')

  system(qsub)
  loopvars$submitted[i] <- T
  Sys.sleep(10)
}

## Check to make sure models are done before continuing
loopvars <- loopvars[loopvars$submitted == T, ]
waitformodelstofinish(lv = cbind(as.character(loopvars[, 1]), loopvars[, 3], 0, 0), sleeptime = 60)

###############################################################################
#Predictive validity
##############################################################################
#Combine the regional input data files
data_master<-rbindlist(lapply(Regions, function(r){
  load(paste0(outputdir, 'processed_input_data_bin0_', r, '_0.rda'))
  df$region <- r
  return(df)
}), use.names=TRUE)
write.csv(data_master, file=paste0(outputdir, '/processed_input_data.csv'))

saveRDS(data_master, file=paste0(outputdir, '/processed_input_data.RDS'))
csvs <- list.files(outputdir, pattern = "input_data_(.*).csv", full.names = T)
csv_master <- rbindlist(lapply(csvs, fread))
csv_master[, V1 := NULL]
write.csv(csv_master, file=paste0(outputdir, '/input_data.csv')) #Will now include the data as age/polygon disaggregated

#Submit job for model validation
qsub <- paste('qsub -e', paste0(outputdir, 'errors'), '-o', paste0(outputdir, 'output'), '-P', project, paste0('-q ', queue), '-l fthread=5 -l m_mem_free=150G -l h_rt=4:08:00:00 -N', paste0('pv_', run_date), '-l archive=TRUE -v sing_image= <<<< FILEPATH REDACTED >>>>lbd_full_20200128.simg -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',paste0('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/share_scripts/shell_sing.sh'),paste0(indic_repo, '/mbg/functions/predictive_validity.r'), indicator, indicator_group, run_date, core_repo, modeling_shapefile_version)
system(qsub)


###############################################################################
## Post-Estimation, Aggregate, and PLHIV
###############################################################################
## Save strata for Shiny to use in producing aggregated fit statistics
strata <- unique(as.character(loopvars[, 1]))
dir.create(paste0(outputdir, '/fit_stats'))
save(strata, file = paste0(outputdir, '/fit_stats/strata.RData'))

#Set up GBD age groups
  gbd_ages<-c(8:16)
  
  #Set WorldPop male & female age groups (to be summed)
  pop_measures_m<-c('a1519m',
                    'a2024m',
                    'a2529m',
                    'a3034m',
                    'a3539m',
                    'a4044m',
                    'a4549m',
                    'a5054m',
                    'a5559m')
  pop_measures_f<-c('a1519f',
                    'a2024f',
                    'a2529f',
                    'a3034f',
                    'a3539f',
                    'a4044f',
                    'a4549f',
                    'a5054f',
                    'a5559f')

pop_measure<-'a1549t'
## Parallelized post-estimation over region, agebin ###########################

for (s in strata) {
  for(agebin in 1:length(z_list)) {
        for(sex in 1:length(sex_list)) {
## Load GBD Estimates for this indicator which will be used in raking
rake_subnational <- TRUE
age_group        <- gbd_ages[agebin]
pop_measurem     <- pop_measures_m[agebin]
pop_measuref     <- pop_measures_f[agebin]
sex_id           <- sex_list[sex]
connector        <- get_gbd_locs(rake_subnational = rake_subnational, reg = s, shapefile_version=shapefile_version)

source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
gbd <- get_outputs(topic = "cause",
                  version = "latest", 
                  gbd_round_id = 6,
                  cause_id = 298,
                  measure_id = 5,
                  metric_id = 3,
                  sex_id = sex,
                  age_group_id = age_group,
                  location_id = connector$location_id,
                  year_id = year_list,
                  decomp_step = 'step5')
gbd <- gbd[, list(name = location_id, year = year_id, mean = val)]


## Prepare for parallel post-estimation - save file with objs to re-load in child processes
source(paste0(indic_repo, 'mbg/', indicator, '/3_functions/prep_postest.r'))
prep_postest(indicator       = indicator,
             indicator_group = indicator_group,
             run_date        = run_date,
             reg             = s,
             age             = z_list[agebin],
             sex             = sex_list[sex],
             save_objs       = c("core_repo", "indic_repo", "gbd", "year_list",
                                 "summstats", "rake_transform", "pop_measurem", "pop_measuref",
                                 "rake_subnational", "age_group", "sex_id",
                                 "shapefile_version", "config", "pop_release"))

  # set memory (based on observed memory use)
  memory <- c(cssa = 500, essa_sdn-COM = 700, sssa = 450, wssa-CPV-STP-MRT = 750)[s]
  
  # make qsub string
  source(paste0(indic_repo, 'mbg/', indicator, '/3_functions/make_qsub_postest.r'))
  qsub <- make_qsub_postest(code           = 'postest_frax_script',
                            stratum        = s,
                            agebin         = z_list[agebin],
                            sex            = sex_list[sex],
                            log_location   = 'sharedir',
                            memory         = memory,
                            cores          = 5,
                            run_time       = '1:05:00:00',
                            proj           ='proj_geo_nodes_hiv',
                            geo_nodes      = FALSE,
                            use_c2_nodes   = FALSE,
                            singularity    = '<<<< FILEPATH REDACTED >>>>/lbd_full_20200128.simg',
                            queue          = 'geospatial.q',
                            modeling_shapefile_version = shapefile_version,
                            raking_shapefile_version = shapefile_version)
  system(qsub)
        }
  }
}

## Check to make sure post-est done before continuing
source(paste0(indic_repo, 'mbg/', indicator, '/3_functions/waitformodelstofinish.r'))
waitformodelstofinish(lv = expand.grid(strata, 0, z_list, sex_list), sleeptime=60, age=z_list, sex=sex_list)


######################################################################################################
## Combine post est stuff across regions and save needed outputs

#Set up to change sex list to be a 3 to call files that got named differently during postestimation
if(sex_list == 0) sex_list <- 3

for(agebin in 1:length(z_list)) {
  for(sex in 1:length(sex_list)) {
post_load_combine_save(regions    = strata,
                              summstats  = c( "mean", "lower", "upper", "cirange"),
                              raked      = c('unraked', "raked", "raked_c"),
                              rf_table   = TRUE,
                              run_summ   = FALSE,
                              indic      = indicator,
                              ig         = indicator_group,
                              sdir       = sharedir,
                              age        = z_list[agebin],
                              sex        = sex_list[sex])

## Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)
}
}
combine_aggregation(rd                   = run_date,
                           indic                = indicator,
                           ig                   = indicator_group,
                           ages                 = z_list,
                           sexes                = sex_list,
                           regions              = strata,
                           holdouts             = 0,
                           raked                = c("unraked", "raked", "raked_c"),
                           dir_to_search        = NULL,
                           delete_region_files  = F,
                           merge_hierarchy_list = F)

summarize_admins(ind = indicator,
                        ig = indicator_group,
                        summstats = c("mean", "lower", "upper", "cirange"),
                        raked    = c( "unraked", "raked", "raked_c"),
                        ad_levels = c(0,1,2),
                        file_addin = NULL,
                        ages = z_list,
                        sexes = sex_list)


###############################################################################
## Make maps and plots
###############################################################################
## Create change uncertainty estimates
source(paste0(indic_repo, "mbg/functions/admin_uncertainty_draws.r"))
for (age in z_list) {
  for (sex in sex_list) {
    try(make_admin_change_uncertainty_tables(run_date, indicator, indicator_group, strata, core_repo, outputdir, age=age, sex = sex))
    try(make_pixel_difference_uncertainty_rasters(run_date, indicator, indicator_group, strata, summstats, core_repo, outputdir, age=age, sex = sex))
  }
}




############################################################################################################################
#Reaggregate estimates across age--run this section manually to create estimates for adults ages 15-49 year and 15-59 years
############################################################################################################################
## SETUP ######################################################################
for(reg in Regions){
  for(agebin in 1:9){
    for(sex in 1:2){
stratum <- reg
holdout <-0

interval_mo <- 12

gbd_ages<-c(8:16)


pathaddin <- paste0('_bin0_',reg,'_',holdout)
load(paste0('<<<< FILEPATH REDACTED >>>>/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))

sex_list<-1:2
z_list <-1:9

#Set WorldPop male & female age groups (to be summed)
pop_measures_m<-c('a1519m',
                  'a2024m',
                  'a2529m',
                  'a3034m',
                  'a3539m',
                  'a4044m',
                  'a4549m',
                  'a5054m',
                  'a5559m')
pop_measures_f<-c('a1519f',
                  'a2024f',
                  'a2529f',
                  'a3034f',
                  'a3539f',
                  'a4044f',
                  'a4549f',
                  'a5054f',
                  'a5559f')

## Parallelized post-estimation over region, agebin ###########################
###########################################################################################################
## Load GBD Estimates for this indicator which will be used in raking
rake_subnational <- TRUE
age_group        <- gbd_ages[agebin]
pop_measurem     <- pop_measures_m[agebin]
pop_measuref     <- pop_measures_f[agebin]
sex_id           <- sex_list[sex]

connector        <- get_gbd_locs(rake_subnational = rake_subnational, reg=reg, shapefile_version=shapefile_version)
year_list <- 2000:2018
shapefile_version <- modeling_shapefile_version <- raking_shapefile_version <- '2020_02_20'

# Define directories
main_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
temp_dir <- paste0(main_dir, "temp_post_est/")
aam_dam_dir <- paste0(main_dir, "reaggregated_results/")

# Get the necessary variables out from the config object into global env
rake_countries <- TRUE
rake_subnational <- TRUE
countries_not_to_rake <- 'ESH+GUF'
countries_not_to_subnat_rake <-'KEN'
year_list <- 2000:2018
metric_space <- 'rates'
pop_release <-'2020_03_20'


# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(stratum)
message(agebin)
message(sex)
message(pop_measurem)
message(pop_measuref)



# For now just assign stratum to reg (will need to modify the below for strata beyond reg)

stratum <- reg<- as.character(commandArgs()[4])
#age <- 0
holdout <- 0


if(sex==1) pop_measure = pop_measurem else if (sex ==2) pop_measure = pop_measuref

## PREPARE RASTERS, ETC. ################################################################
crosswalk <- F

# Assume linear raking unless specified as logit
rake_method <- "linear"

##### Prep input data into raking:

## Get the simple and new_simple rasters prepped up for us
print("Getting simple and prepped rasters")
raster_outputs <- prep_shapes_for_raking(
  reg = reg,
  modeling_shapefile_version = modeling_shapefile_version,
  raking_shapefile_version = raking_shapefile_version,
  field = "loc_id"
)

## Take out the objects from the list that actually matters to us:
simple_raster <- raster_outputs[["simple_raster"]]
new_simple_raster <- raster_outputs[["new_simple_raster"]]

simple_polygon <- raster_outputs[["simple_polygon"]]
new_simple_polygon <- raster_outputs[["new_simple_polygon"]]

pixel_id <- raster_outputs[["pixel_id"]]


#####################################################################
#load the cell id to admin units link
link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)

#####################################################################
# collect and load the population data from the WorldPop rasters
print(reg)
print(pop_measure)
print(year_list)
print(interval_mo)
covdt <- load_populations_cov(reg, pop_measure, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
print(nrow(covdt))
#####################################################################
# Prepping the cell_pred and link table to be linked by making sure they have the appropriate identifiers.  Also performs a 
# zippering at the region boundary where cells that have a portion of their area outside of the modeling region are reintegrated
# as a whole cell and considered to be only in one region.  This works becasue a given cell is only modeled in one region.

connector <- get_gbd_locs(rake_subnational = rake_subnational,
                          reg = reg,
                          shapefile_version = shapefile_version)

# getting the connector for sub-national raking - used to implement countries_not_to_subnat_rake
nat_connector <- get_gbd_locs(rake_subnational = F,
                              reg = reg,
                              shapefile_version = shapefile_version)

cell_ids <- link_table[[2]]

link <- prep_link_table(link_table = link_table,
                        simple_raster = simple_raster,
                        pixel_id = pixel_id)


# getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
# gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary 

# merge the connectors on to the link table
link <- sub_nat_link_merge(rake_subnational,
                           link,
                           connector,
                           nat_connector,
                           countries_not_to_subnat_rake)

#############################################################################################
# Load cell draws
message("Loading Data...")
load(paste0(main_dir, indicator, "_cell_draws_eb_bin", agebin, "_sex", sex, "_", reg, "_0.RData"))
##### Using fractional raking #####

## First, create all the fractional rake factors
# setting a reference for the number of draws
ndraws = ncol(cell_pred)

# set cell pred as a data table, and rename things
cell_pred <- prep_cell_pred(cell_pred = cell_pred,
                            cell_ids  = cell_ids,
                            pixel_id  = pixel_id,
                            covdt     = covdt)

# merge cell_pred on the link
cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)


## Raking Population ###################################################################
# This is done to ensure that the total pop in each raking geography is the same as GBD
message("raking population")

#convert to fractional population 
cell_pred = cell_pred[,pop_af := pop * area_fraction] 

#NA out population where the pixel value is NA (to prevent weirdness with denominators)
cell_pred = cell_pred[is.na(V1), pop_af := NA]


scalars<-setDT(read.csv(paste0(outputdir, 'hiv_prev_disagg_bin',agebin, '_sex', sex, '_', reg, '_pop_rf.csv')))

# add back to the cell_pred as a population rf
cell_pred <- merge(cell_pred, scalars, by = c("location_id", "year"))

# rake fractional populations
cell_pred$pop_raked <- 0
cell_pred = cell_pred[,pop_raked := pop_af * pop_scalar]
cell_pred$pop <- NULL
cell_pred$pop_af <- NULL

## load the raking factors
fractional_rf <- setDT(read.csv(file =  paste0(sharedir, "/output/", run_date, "/", indicator, "_bin", agebin, "_sex", sex, "_", reg, "_rf.csv") ))

## merge them onto the cell_pred
cell_pred <- merge(cell_pred, fractional_rf, by = c("location_id", "year"))

############################################################
# creating a raked rates cell_pred (this happens first becasue once we go to counts in the cell_pred we can't do back without loading a fresh copy)
message("creating a raked rates cell_pred")
# rake rates
overs <- paste0("V", 1:ndraws)

cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * rf)]

#######
message("converting from prevalence to counts")
# set the variables to aggregate
overs <- paste0("V", 1:ndraws)

cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * pop_raked)]
# NA out population where the pixel value is NA (to prevent weirdness with denominators)
cell_pred <- cell_pred[is.na(V1), pop_raked := NA]

#
cell_pred[ , plhiv_mean := apply(.SD, 1, mean, na.rm=TRUE), .SDcols=overs]
cell_pred[ , plhiv_upper := apply(.SD, 1, quantile, c(.975), na.rm=TRUE), .SDcols=overs]
cell_pred[ , plhiv_lower := apply(.SD, 1, quantile, c(.025), na.rm=TRUE), .SDcols=overs]


saveRDS(cell_pred[,list(pop_raked,plhiv_mean,plhiv_lower,plhiv_upper,area_fraction,ADM0_CODE,ADM1_CODE,ADM2_CODE,year,cell_pred_id)], 
        paste0(aam_dam_dir, 'fractionalized_cps_age', agebin, '_sex', sex, '_', reg, '.RDS'))
    }
  }
}

message('done with fractionalizing plhiv!')

###########################################################################################
##Sum fractionalized PLHIV to ages 15-49 and calculate associated 15-49 prevalence
############################################################################################
for(reg in Regions){
  sex_list<-1:2
  z_list <-1:7

  ## Parallelized post-estimation over region, agebin ###########################
  ###########################################################################################################
  ## Load GBD Estimates for this indicator which will be used in raking
  rake_subnational <- TRUE
  connector        <- get_gbd_locs(rake_subnational = rake_subnational, reg=reg, shapefile_version=shapefile_version)
  year_list <- 2000:2018
  
  
  
  # Define directories
  main_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  temp_dir <- paste0(main_dir, "temp_post_est/")
  aam_dam_dir <- paste0(main_dir, "reaggregated_results/")
  
  pop_measure <- 'a1549t'
  stratum <- reg
  # Print some settings to console
  message(indicator)
  message(indicator_group)
  message(run_date)
  message(stratum)
  message(pop_release)
  
  # Print raking info
  print(paste0("Metric Space                       : ", metric_space))
  print(paste0("Subnational raking                 : ", rake_subnational))
  print(paste0("Countries not to rake at all       : ", countries_not_to_rake))
  print(paste0("Countries not to rake subnationally: ", countries_not_to_subnat_rake))
  
  
  
  # For now just assign stratum to reg (will need to modify the below for strata beyond reg)
  reg <- stratum
  #age <- 0
  holdout <- 0
  
  ## Get the simple and new_simple rasters prepped up for us
  print("Getting simple and prepped rasters")
  raster_outputs <- prep_shapes_for_raking(
    reg = reg,
    modeling_shapefile_version = modeling_shapefile_version,
    raking_shapefile_version = raking_shapefile_version,
    field = "loc_id"
  )
  
  ## Take out the objects from the list that actually matters to us:
  simple_raster <- raster_outputs[["simple_raster"]]
  ######################################################
  #Read in and add all the plhiv and pop values together
  df_all<-readRDS(paste0(aam_dam_dir, 'fractionalized_cps_age1_sex1_', reg, '.RDS'))
  df_all[,plhiv_mean:=0]
  df_all[,plhiv_upper:=0]
  df_all[,plhiv_lower:=0]
  df_all[,pop_raked:=0]
  
  for(sex in 1:2){
    for(agebin in 1:7){
      df<-readRDS(paste0(aam_dam_dir, 'fractionalized_cps_age', agebin, '_sex', sex, '_', reg, '.RDS'))
      df_all[,pop_raked := pop_raked + df$pop_raked]
      df_all[,plhiv_mean := plhiv_mean + df$plhiv_mean]
      df_all[,plhiv_lower := plhiv_lower + df$plhiv_lower]
      df_all[,plhiv_upper := plhiv_upper + df$plhiv_upper]
      remove(df)
    }
  }
  
  #Estimate prevalence where pop_raked==0
  df_all[,prev_mean:=0]
  df_all[,prev_upper:=0]
  df_all[,prev_lower:=0]
  
  adm2<-df_all[!is.na(pop_raked),list(prev_mean=sum(plhiv_mean)/sum(pop_raked),
                                      prev_upper=sum(plhiv_upper)/sum(pop_raked),
                                      prev_lower=sum(plhiv_lower)/sum(pop_raked)), by='year,ADM2_CODE']
  
  #Set prevalence for areas with pop==0, by the prevalence at the admin2 level (fractionalized)
  lapply(unique(adm2$ADM2_CODE), function(adm2_code){
    lapply(unique(2000:2018), function(y){
      df_all[pop_raked==0 & year==y & ADM2_CODE==adm2_code, prev_mean:=adm2[year==y & ADM2_CODE==adm2_code]$prev_mean]
      df_all[pop_raked==0 & year==y & ADM2_CODE==adm2_code, prev_upper:=adm2[year==y & ADM2_CODE==adm2_code]$prev_upper]
      df_all[pop_raked==0 & year==y & ADM2_CODE==adm2_code, prev_lower:=adm2[year==y & ADM2_CODE==adm2_code]$prev_lower]
    })
  })
  
  df_all[pop_raked==0, prev_mean:=prev_mean*area_fraction]
  df_all[pop_raked==0, prev_lower:=prev_lower*area_fraction]
  df_all[pop_raked==0, prev_upper:=prev_upper*area_fraction]
  #############################################################################
  #Admin summaries
  admin_0 <- df_all[, lapply(c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM0_CODE")]
  names(admin_0)[3:6] <- c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked')
  admin_0[,prev_mean := plhiv_mean/pop_raked]
  admin_0[,prev_lower := plhiv_lower/pop_raked]
  admin_0[,prev_upper := plhiv_upper/pop_raked]
  admin_0[,ADM1_CODE:=NA]
  
  admin_1 <- df_all[, lapply(c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM0_CODE", "ADM1_CODE")]
  names(admin_1)[4:7] <- c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked')
  admin_1[,prev_mean := plhiv_mean/pop_raked]
  admin_1[,prev_lower := plhiv_lower/pop_raked]
  admin_1[,prev_upper := plhiv_upper/pop_raked]
  admins<-rbind(admin_0, admin_1)
  saveRDS(admins, paste0(aam_dam_dir, 'admin_summaries_1549_', reg, '.RDS'))
  
  #Rasters
  unfrac_ras <- df_all[, lapply(c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked', 'prev_mean','prev_upper','prev_lower'), function(x) sum(get(x), na.rm = T)), by = c("year", "cell_pred_id")]
  setorder(unfrac_ras, cell_pred_id)
  names(unfrac_ras)[3:9] <- c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked','prev_mean','prev_upper','prev_lower')
  unfrac_ras[pop_raked!=0,prev_mean := plhiv_mean/pop_raked]
  unfrac_ras[pop_raked!=0,prev_lower := plhiv_lower/pop_raked]
  unfrac_ras[pop_raked!=0,prev_upper := plhiv_upper/pop_raked]
  
  ######Make the rasters######################
  plhiv_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$plhiv_mean), ncol=19))
  writeRaster(plhiv_ras, paste0(outputdir, 'hiv_prev_disagg_bin1549_sex0_', reg, '_raked_c_mean_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  prev_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$prev_mean), ncol=19))
  writeRaster(prev_ras, paste0(outputdir, 'hiv_prev_disagg_bin1549_sex0_', reg, '_raked_mean_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  
  
  
  plhiv_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$plhiv_upper), ncol=19))
  writeRaster(plhiv_ras, paste0(outputdir, 'hiv_prev_disagg_bin1549_sex0_', reg, '_raked_c_upper_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  prev_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$prev_upper), ncol=19))
  writeRaster(prev_ras, paste0(outputdir, 'hiv_prev_disagg_bin1549_sex0_', reg, '_raked_upper_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  
  
  
  plhiv_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$plhiv_lower), ncol=19))
  writeRaster(plhiv_ras, paste0(outputdir, 'hiv_prev_disagg_bin1549_sex0_', reg, '_raked_c_lower_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  prev_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$prev_lower), ncol=19))
  writeRaster(prev_ras, paste0(outputdir, 'hiv_prev_disagg_bin1549_sex0_', reg, '_raked_lower_raster.tif'),
              format='GTiff', overwrite=TRUE)
}
###########################################################################################
###########################################################################################
##Sum fractionalized PLHIV to ages 15-59 and calculate associated 15-59 prevalence
############################################################################################
## Parallelized over region, agebin, sex ###########################
for(reg in Regions){
  sex_list<-1:2
  z_list <-1:9
  ###########################################################################################################
  ## Load GBD Estimates for this indicator which will be used in raking
  rake_subnational <- TRUE
  connector        <- get_gbd_locs(rake_subnational = rake_subnational, reg=reg, shapefile_version=shapefile_version)
  year_list <- 2000:2018
  
  
  
  # Define directories
  main_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  temp_dir <- paste0(main_dir, "temp_post_est/")
  aam_dam_dir <- paste0(main_dir, "reaggregated_results/")
  
  pop_measure <- 'a1549t'
  stratum <- reg
  # Print some settings to console
  message(indicator)
  message(indicator_group)
  message(run_date)
  message(stratum)
  message(pop_release)
  
  # Print raking info
  print(paste0("Metric Space                       : ", metric_space))
  print(paste0("Subnational raking                 : ", rake_subnational))
  print(paste0("Countries not to rake at all       : ", countries_not_to_rake))
  print(paste0("Countries not to rake subnationally: ", countries_not_to_subnat_rake))
  
  
  
  # For now just assign stratum to reg (will need to modify the below for strata beyond reg)
  reg <- stratum
  #age <- 0
  holdout <- 0
  
  ## Get the simple and new_simple rasters prepped up for us
  print("Getting simple and prepped rasters")
  raster_outputs <- prep_shapes_for_raking(
    reg = reg,
    modeling_shapefile_version = modeling_shapefile_version,
    raking_shapefile_version = raking_shapefile_version,
    field = "loc_id"
  )
  
  ## Take out the objects from the list that actually matters to us:
  simple_raster <- raster_outputs[["simple_raster"]]
  ######################################################
  #Read in and add all the plhiv and pop values together
  df_all<-readRDS(paste0(aam_dam_dir, 'fractionalized_cps_age1_sex1_', reg, '.RDS'))
  df_all[,plhiv_mean:=0]
  df_all[,plhiv_upper:=0]
  df_all[,plhiv_lower:=0]
  df_all[,pop_raked:=0]
  
  for(sex in 1:2){
    for(agebin in 1:9){
      df<-readRDS(paste0(aam_dam_dir, 'fractionalized_cps_age', agebin, '_sex', sex, '_', reg, '.RDS'))
      df_all[,pop_raked := pop_raked + df$pop_raked]
      df_all[,plhiv_mean := plhiv_mean + df$plhiv_mean]
      df_all[,plhiv_lower := plhiv_lower + df$plhiv_lower]
      df_all[,plhiv_upper := plhiv_upper + df$plhiv_upper]
      remove(df)
    }
  }
  
  #Estimate prevalence where pop_raked==0
  df_all[,prev_mean:=0]
  df_all[,prev_upper:=0]
  df_all[,prev_lower:=0]
  
  adm2<-df_all[!is.na(pop_raked),list(prev_mean=sum(plhiv_mean)/sum(pop_raked),
                                      prev_upper=sum(plhiv_upper)/sum(pop_raked),
                                      prev_lower=sum(plhiv_lower)/sum(pop_raked)), by='year,ADM2_CODE']
  
  #Set prevalence for areas with pop==0, by the prevalence at the admin2 level (fractionalized)
  lapply(unique(adm2$ADM2_CODE), function(adm2_code){
    lapply(unique(2000:2018), function(y){
      df_all[pop_raked==0 & year==y & ADM2_CODE==adm2_code, prev_mean:=adm2[year==y & ADM2_CODE==adm2_code]$prev_mean]
      df_all[pop_raked==0 & year==y & ADM2_CODE==adm2_code, prev_upper:=adm2[year==y & ADM2_CODE==adm2_code]$prev_upper]
      df_all[pop_raked==0 & year==y & ADM2_CODE==adm2_code, prev_lower:=adm2[year==y & ADM2_CODE==adm2_code]$prev_lower]
    })
  })
  
  df_all[pop_raked==0, prev_mean:=prev_mean*area_fraction]
  df_all[pop_raked==0, prev_lower:=prev_lower*area_fraction]
  df_all[pop_raked==0, prev_upper:=prev_upper*area_fraction]
  #############################################################################
  #Admin summaries
  admin_0 <- df_all[, lapply(c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM0_CODE")]
  names(admin_0)[3:6] <- c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked')
  admin_0[,prev_mean := plhiv_mean/pop_raked]
  admin_0[,prev_lower := plhiv_lower/pop_raked]
  admin_0[,prev_upper := plhiv_upper/pop_raked]
  admin_0[,ADM1_CODE:=NA]
  
  admin_1 <- df_all[, lapply(c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM0_CODE", "ADM1_CODE")]
  names(admin_1)[4:7] <- c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked')
  admin_1[,prev_mean := plhiv_mean/pop_raked]
  admin_1[,prev_lower := plhiv_lower/pop_raked]
  admin_1[,prev_upper := plhiv_upper/pop_raked]
  admins<-rbind(admin_0, admin_1)
  saveRDS(admins, paste0(aam_dam_dir, 'admin_summaries_1559_', reg, '.RDS'))
  ##############################################################################
  #Rasters
  unfrac_ras <- df_all[, lapply(c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked', 'prev_mean','prev_upper','prev_lower'), function(x) sum(get(x), na.rm = T)), by = c("year", "cell_pred_id")]
  setorder(unfrac_ras, cell_pred_id)
  names(unfrac_ras)[3:9] <- c('plhiv_mean', 'plhiv_lower', 'plhiv_upper', 'pop_raked','prev_mean','prev_upper','prev_lower')
  unfrac_ras[pop_raked!=0,prev_mean := plhiv_mean/pop_raked]
  unfrac_ras[pop_raked!=0,prev_lower := plhiv_lower/pop_raked]
  unfrac_ras[pop_raked!=0,prev_upper := plhiv_upper/pop_raked]
  
  ######Make the rasters######################
  plhiv_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$plhiv_mean), ncol=19))
  writeRaster(plhiv_ras, paste0(outputdir, 'hiv_prev_disagg_bin1559_sex0_', reg, '_raked_c_mean_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  prev_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$prev_mean), ncol=19))
  writeRaster(prev_ras, paste0(outputdir, 'hiv_prev_disagg_bin1559_sex0_', reg, '_raked_mean_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  
  
  
  plhiv_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$plhiv_upper), ncol=19))
  writeRaster(plhiv_ras, paste0(outputdir, 'hiv_prev_disagg_bin1559_sex0_', reg, '_raked_c_upper_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  prev_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$prev_upper), ncol=19))
  writeRaster(prev_ras, paste0(outputdir, 'hiv_prev_disagg_bin1559_sex0_', reg, '_raked_upper_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  
  
  
  plhiv_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$plhiv_lower), ncol=19))
  writeRaster(plhiv_ras, paste0(outputdir, 'hiv_prev_disagg_bin1559_sex0_', reg, '_raked_c_lower_raster.tif'),
              format='GTiff', overwrite=TRUE)
  
  prev_ras <- insertRaster(simple_raster, matrix(as.matrix(unfrac_ras$prev_lower), ncol=19))
  writeRaster(prev_ras, paste0(outputdir, 'hiv_prev_disagg_bin1559_sex0_', reg, '_raked_lower_raster.tif'),
              format='GTiff', overwrite=TRUE)
}
###########################################################################################



###############################################################################
#Cap upper prevalence estimates at 100% prevalence
#Raking to GBD estimates causes a negligible number of draws to be estimated as >100% prevalence. 
#Therefore need to cap our upper (97.5% upper CI) estimates at 100% prevalence, and adjust all subsequent dependent estimates
###############################################################################
#Capping upper and mean raster estimates at 1
for(z in 1:9){
  for(x in 1:2){
    s<-stack(paste0('<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_upper_raked_bin',z,'_sex',x,'_raster.tif'))
    for(i in 1:19){
      if(s[[i]]@data@max >1){ s[[i]] <- min(s[[i]],1)
      }
    }
    
    writeRaster(s, paste0('<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_upper_raked_bin',z,'_sex',x,'_raster.tif'))
    
  }
}


for(z in 1:9){
  for(x in 1:2){
    s<-stack(paste0('<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_mean_raked_bin',z,'_sex',x,'_raster.tif'))
    for(i in 1:19){
      if(s[[i]]@data@max >1){ s[[i]] <- min(s[[i]],1)
      }
    }
    
    writeRaster(s, paste0('<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_mean_raked_bin',z,'_sex',x,'_raster.tif'))
    
  }
}


#Correspondingly, adjust PLHIV estimates
for(z in 1:9){
  for(x in 1:2){
    count_upper<-stack(paste0('<<<< FILEPATH REDACTED >>>>//hiv_prev_disagg_upper_raked_c_bin',z,'_sex',x,'_raster.tif'))
    count_mean<-stack(paste0('<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_mean_raked_c_bin',z,'_sex',x,'_raster.tif'))
    
    
    old_mean <-stack(paste0('<<<< FILEPATH REDACTED >>>>hiv_prev_disagg_mean_raked_bin',z,'_sex',x,'_raster.tif'))
    old_upper <-stack(paste0('<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_upper_raked_bin',z,'_sex',x,'_raster.tif'))
    
    new_mean <-stack(paste0('<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_mean_raked_bin',z,'_sex',x,'_raster.tif'))
    new_upper <-stack(paste0('<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_upper_raked_bin',z,'_sex',x,'_raster.tif'))
    
    
    old_upper <-raster::mask(old_upper, count_upper)
    old_mean <-raster::mask(old_mean, count_mean)  
    
    
    #calculate proportion of plhiv in capped vs uncapped age aggregations
    count_upper_df <-data.table(rasterToPoints(count_upper))
    count_mean_df <-data.table(rasterToPoints(count_mean))
    old_upper_df <-data.table(rasterToPoints(old_upper))
    old_mean_df <-data.table(rasterToPoints(old_mean))
    new_upper_df <-data.table(rasterToPoints(new_upper))
    new_mean_df <-data.table(rasterToPoints(new_mean))
    
    
    upper_prop <- count_upper_df[,1:2]  
    mean_prop <- count_mean_df[,1:2]  
    
    for(a in 3:21){
      
      upper_prop[[a]] <- new_upper_df[[a]] / old_upper_df[[a]]
      mean_prop[[a]] <- new_mean_df[[a]] / old_mean_df[[a]]
      if(any(is.na(upper_prop[[a]]))) upper_prop[which(is.na(upper_prop[[a]]))][[a]]<-1
      if(any(is.na(mean_prop[[a]])))  mean_prop[which(is.na(mean_prop[[a]]))][[a]]<-1
      
    }
    
    upper_prop_ras <- insertRaster(count_upper, matrix(as.matrix(upper_prop[,3:21]), ncol=19))
    mean_prop_ras <- insertRaster(count_mean, matrix(as.matrix(mean_prop[,3:21]), ncol=19))
    
    new_count_mean <- count_mean * mean_prop_ras
    new_count_upper <- count_upper * upper_prop_ras
    
    writeRaster(new_count_mean, paste0( '<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_mean_raked_c_bin',z,'_sex',x,'_raster.tif'),
                format='GTiff', overwrite=TRUE)
    writeRaster(new_count_upper, paste0( '<<<< FILEPATH REDACTED >>>>/hiv_prev_disagg_upper_raked_c_bin',z,'_sex',x,'_raster.tif'),
                format='GTiff', overwrite=TRUE)
    
    
  }
}

#########################################################################################################
#Combine regional 15-49 & 15-59 raster estimates
post_load_combine_save(regions    = strata,
                       summstats  = c( "mean", "lower", "upper"),
                       raked      = c( "raked", "raked_c"),
                       rf_table   = FALSE,
                       run_summ   = FALSE,
                       indic      = indicator,
                       ig         = indicator_group,
                       sdir       = sharedir,
                       age        = '1549',
                       sex        = as.character(0))

post_load_combine_save(regions    = strata,
                       summstats  = c( "mean", "lower", "upper"),
                       raked      = c( "raked", "raked_c"),
                       rf_table   = FALSE,
                       run_summ   = FALSE,
                       indic      = indicator,
                       ig         = indicator_group,
                       sdir       = sharedir,
                       age        = '1559',
                       sex        = as.character(0))


###############################################################################################################################################
#Capping upper admin estimates at 100% prevalence

#Cap upper prevalence estimates for admin2
for(x in 1:2){
  for(z in 1:9){
    s<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_2_raked_bin',z,'_sex',x,'_summary.csv'))
    
    print(nrow(s[upper>1]))
    summary(s$upper)
    
    s[upper>1, upper:=1]
    summary(s$upper)
    
    write.csv(s,file=paste0('<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_2_raked_bin',z,'_sex',x,'_summary.csv'))
    
  }
}

#Correspondingly adjust upper PLHIV estimates for admin2
for(z in 1:9){
  for(x in 1:2){
    count<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_2_raked_c_bin',z,'_sex',x,'_summary.csv'))
    
    
    old <-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_2_raked_bin',z,'_sex',x,'_summary.csv'))
    
    new <-fread(paste0('<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_2_raked_bin',z,'_sex',x,'_summary.csv'))
    
    upper_prop <- new$upper / old$upper
    if(any(is.na(upper_prop))) upper_prop[which(is.na(upper_prop))]<-1
    
    
    count$upper <- count$upper * upper_prop
    
    write.csv(count, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_2_raked_c_bin',z,'_sex',x,'_summary.csv'))
    
  }
}


########################################################################
#Save the admin1s and admin0s in the same place to be more easily adjusted to correspond with capped admin2 estimates
for(x in 1:2){
  for(z in 1:9){
    for(a in 0:1){
      for(c in c('', 'c_')){
        s<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_',a,'_raked_',c,'bin',z,'_sex',x,'_summary.csv'))
        
        write.csv(s,file=paste0('<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_',a,'_raked_',c,'bin',z,'_sex',x,'_summary.csv'))
      }
    }
  }
}

##Identify which admin1s and admin0s will need to also have upper estimates correspondingly adjusted
for(z in 1:9){
  for(x in 1:2){
    
    old <-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_2_raked_c_bin',z,'_sex',x,'_summary.csv'))
    
    new <-fread(paste0('<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_2_raked_c_bin',z,'_sex',x,'_summary.csv'))
    
    if(any(old$upper != new$upper)){
      
      print(paste0('age ', z))
      print(paste0('sex ', x))
      print(old[which(old$upper != new$upper)]$year)
      print(old[which(old$upper != new$upper)]$ADM1_CODE)
      print(old[which(old$upper != new$upper)]$ADM0_CODE)
      print(old[which(old$upper != new$upper)]$upper - new[which(old$upper != new$upper)]$upper)
    }
  }
}

######################################################################################
#Subtracting the plhiv difference in those >100% admin2s from admin1s and admin0s
####Bin 5 sex 1#######################################
count1<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_1_raked_c_bin5_sex1_summary.csv'))
ad1_c<-copy(count1)
count0<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_0_raked_c_bin5_sex1_summary.csv'))
ad0_c<-copy(count0)


prev1<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_1_raked_bin5_sex1_summary.csv'))
ad1<-copy(prev1)
prev0<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_0_raked_bin5_sex1_summary.csv'))
ad0<-copy(prev0)


ad1_c[year==2000 & ADM1_CODE==27118, upper:= upper - 155.6818]
ad0_c[year==2000 & ADM0_CODE==118,   upper:= upper - 155.6818]
ad1_c[year==2000 & ADM1_CODE==27118, cirange:= cirange - 155.6818]
ad0_c[year==2000 & ADM0_CODE==118,   cirange:= cirange - 155.6818]

prop1<-ad1_c$upper/count1$upper
prop0<-ad0_c$upper/count0$upper

prev1[, upper := upper * prop1]
prev1[, cirange := upper- lower]
prev0[, upper := upper * prop0]
prev0[, cirange := upper- lower]


write.csv(ad1_c, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_1_raked_c_bin5_sex1_summary.csv'))
write.csv(ad0_c, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_0_raked_c_bin5_sex1_summary.csv'))
write.csv(prev1, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_1_raked_bin5_sex1_summary.csv'))
write.csv(prev0, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_0_raked_bin5_sex1_summary.csv'))



####################################
#Adjusting the identified estimates: 
#Bin 6 sex1#########################
count1<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_1_raked_c_bin6_sex1_summary.csv'))
ad1_c<-copy(count1)
count0<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_0_raked_c_bin6_sex1_summary.csv'))
ad0_c<-copy(count0)


prev1<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_1_raked_bin6_sex1_summary.csv'))
ad1<-copy(prev1)
prev0<-fread(paste0('<<<< FILEPATH REDACTED >>>>/pred_derivatives/admin_summaries/hiv_prev_disagg_admin_0_raked_bin6_sex1_summary.csv'))
ad0<-copy(prev0)


ad1_c[year==2000 & ADM1_CODE==27118, upper:= upper - 146.84307]
ad0_c[year==2000 & ADM0_CODE==118,   upper:= upper - 146.84307]
ad1_c[year==2000 & ADM1_CODE==27118, cirange:= cirange - 146.84307]
ad0_c[year==2000 & ADM0_CODE==118,   cirange:= cirange - 146.84307]

ad1_c[year==2001 & ADM1_CODE==27118, upper:= upper - 112.02862]
ad0_c[year==2001 & ADM0_CODE==118,   upper:= upper - 112.02862]
ad1_c[year==2001 & ADM1_CODE==27118, cirange:= cirange - 112.02862]
ad0_c[year==2001 & ADM0_CODE==118,   cirange:= cirange - 112.02862]


ad1_c[year==2002 & ADM1_CODE==27118, upper:= upper - 19.45221]
ad0_c[year==2002 & ADM0_CODE==118,   upper:= upper - 19.45221]
ad1_c[year==2002 & ADM1_CODE==27118, cirange:= cirange - 19.45221]
ad0_c[year==2002 & ADM0_CODE==118,   cirange:= cirange - 19.45221]


prop1<-ad1_c$upper/count1$upper
prop0<-ad0_c$upper/count0$upper

prev1[, upper := upper * prop1]
prev1[, cirange := upper- lower]
prev0[, upper := upper * prop0]
prev0[, cirange := upper- lower]


write.csv(ad1_c, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_1_raked_c_bin6_sex1_summary.csv'))
write.csv(ad0_c, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_0_raked_c_bin6_sex1_summary.csv'))
write.csv(prev1, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_1_raked_bin6_sex1_summary.csv'))
write.csv(prev0, paste0( '<<<< FILEPATH REDACTED >>>>/capped_admins/hiv_prev_disagg_admin_0_raked_bin6_sex1_summary.csv'))







##################################################################################
#Aggregate estimates to 15-49 & 15-59 at the admin0,1, & 2 levels
##################################################################################
disagg_maps_path <- paste0("<<<< FILEPATH REDACTED >>>>/output/",run_date)

outdir <- paste0("<<<< FILEPATH REDACTED >>>>/output/",run_date, '/capped_admins/')

#rbindlist(lapply(c('raked', 'unraked'), function(rake_level){
pred_reagged <-rbindlist(lapply(c('to49','to59'), function(age_level){
  pred_reagged <-rbindlist(lapply(c(0,1,2), function(ad_level){
    if (ad_level==0){
      agg_levels<-'ADM0_CODE,ADM0_NAME,region,year,adm_level,agegroup'
    } else if(ad_level==1){
      agg_levels<-'ADM0_CODE,ADM0_NAME,ADM1_CODE,ADM1_NAME,region,year,adm_level,agegroup'
    } else if(ad_level==2){
      agg_levels<-'ADM0_CODE,ADM0_NAME,ADM1_CODE,ADM1_NAME,ADM2_CODE,ADM2_NAME,region,year,adm_level,agegroup'
    }
    
    if (age_level=='to49'){
      ages <- 1:7
    } else if(age_level=='to59'){
      ages <- 1:9
    }
    
    sexes <- 1:2
    
    # # pred_disagg
    pred_disagg <- rbindlist(lapply(ages, function(bin) {
      pred_disagg <- rbindlist(lapply(sexes, function(mf) {
        pred_disagg <- fread(paste0(disagg_maps_path, '/capped_admins/hiv_prev_disagg_admin_', ad_level, '_raked_bin',bin,"_sex", mf,"_" ,'summary.csv'))
        pred_disagg$agebin = bin
        pred_disagg$sex_id = mf
        pred_disagg$adm_level = ad_level
        pred_disagg$agegroup = age_level
        return(pred_disagg)
      }))
    }))
    
    
    pred_disagg_c <- rbindlist(lapply(ages, function(bin) {
      pred_disagg_c <- rbindlist(lapply(sexes, function(mf) {
        pred_disagg_c <- fread(paste0(disagg_maps_path, '/capped_admins/hiv_prev_disagg_admin_', ad_level, '_raked_c_bin',bin,"_sex", mf,"_" ,'summary.csv'))
        pred_disagg_c$agebin = bin
        pred_disagg_c$sex_id = mf
        pred_disagg_c$adm_level = ad_level
        pred_disagg_c$agegroup = age_level
        return(pred_disagg_c)
      }))
    }))
    
    
    #Weight by pop & re-aggregate
    pred_disagg[,plhiv_mean:=pred_disagg_c$mean]
    pred_disagg[,plhiv_upper:=pred_disagg_c$upper]
    pred_disagg[,plhiv_lower:=pred_disagg_c$lower]
    pred_disagg[,pop_mean:=plhiv_mean/mean]
    pred_disagg[,pop_upper:=plhiv_upper/upper]
    pred_disagg[,pop_lower:=plhiv_upper/lower]
    pred_disagg[,crude_wt_mean := pop_mean/sum(pop_mean), by=agg_levels]
    pred_disagg[,crude_wt_upper := pop_upper/sum(pop_upper), by=agg_levels]
    pred_disagg[,crude_wt_lower := pop_lower/sum(pop_lower), by=agg_levels]
    pred_disagg[,reagged_prev_mean:=mean*crude_wt_mean]
    pred_disagg[,reagged_prev_upper:=upper*crude_wt_upper]
    pred_disagg[,reagged_prev_lower:=lower*crude_wt_lower]
    
    pred_reagged<- pred_disagg[,list(prev_mean=sum(reagged_prev_mean), prev_upper=sum(reagged_prev_upper), prev_lower=sum(reagged_prev_lower),
                                     plhiv_mean=sum(plhiv_mean), plhiv_upper=sum(plhiv_upper), plhiv_lower=sum(plhiv_lower),
                                     pop_mean=sum(pop_mean), pop_upper=sum(pop_upper), pop_lower=sum(pop_lower)), 
                               by=agg_levels]
    
    if(!('ADM1_CODE' %in% names(pred_reagged))) pred_reagged$ADM1_CODE<-NA
    if(!('ADM1_NAME' %in% names(pred_reagged))) pred_reagged$ADM1_NAME<-NA
    if(!('ADM2_CODE' %in% names(pred_reagged))) pred_reagged$ADM2_CODE<-NA
    if(!('ADM2_NAME' %in% names(pred_reagged))) pred_reagged$ADM2_NAME<-NA
    
    return(pred_reagged)
  }), use.names = TRUE)
}), use.names = TRUE)


write.csv(pred_reagged[adm_level==0 & agegroup=='to49', ], paste0(disagg_maps_path, '/capped_admins/hiv_prev_disagg_admin_0_raked_bin','1549',"_sex0_" ,'summary.csv'))
write.csv(pred_reagged[adm_level==1 & agegroup=='to49', ], paste0(disagg_maps_path, '/capped_admins/hiv_prev_disagg_admin_1_raked_bin','1549',"_sex0_" ,'summary.csv'))
write.csv(pred_reagged[adm_level==2 & agegroup=='to49', ], paste0(disagg_maps_path, '/capped_admins/hiv_prev_disagg_admin_2_raked_bin','1549',"_sex0_" ,'summary.csv'))

write.csv(pred_reagged[adm_level==0 & agegroup=='to59', ], paste0(disagg_maps_path, '/capped_admins/hiv_prev_disagg_admin_0_raked_bin','1559',"_sex0_" ,'summary.csv'))
write.csv(pred_reagged[adm_level==1 & agegroup=='to59', ], paste0(disagg_maps_path, '/capped_admins/hiv_prev_disagg_admin_1_raked_bin','1559',"_sex0_" ,'summary.csv'))
write.csv(pred_reagged[adm_level==2 & agegroup=='to59', ], paste0(disagg_maps_path, '/capped_admins/hiv_prev_disagg_admin_2_raked_bin','1559',"_sex0_" ,'summary.csv'))

###############################################################################
## Close out
###############################################################################

run_date_file <- paste0('<<<< FILEPATH REDACTED >>>>/models_submitted_', substr(run_date, 1, 10), '.csv')
if (file.exists(run_date_file)) {
  temp <- fread(run_date_file)
  temp[run_date == get('run_date', .GlobalEnv), done := 1]
  write.table(temp, file = run_date_file, sep = ",", row.names = F, na = '')
}

message('Done!')
