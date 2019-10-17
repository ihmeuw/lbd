# load in the arguments
user            <- "<<<< FILEPATH REDACTED >>>>"
core_repo       <- sprintf("<<<< FILEPATH REDACTED >>>>") # main MBG repo
ig_repo         <- sprintf("<<<< FILEPATH REDACTED >>>>")      # u5m specific repo
indicator       <- 'died'
indicator_group <- 'u5m'

# load central lbd code stuff
source(paste0("<<<< FILEPATH REDACTED >>>>"))
pl <- c("rgeos",      "data.table", "raster",     "rgdal",      "INLA",
        "seegSDM",    "seegMBG",    "dismo",      "gbm",        "foreign",
        "parallel",   "doParallel", "grid",       "gridExtra",  "pacman",
        "gtools",     "glmnet",     "ggplot2",    "RMySQL",     "plyr",
        "tictoc",     "dplyr",      "magrittr",   "tidyr",      "sp",
        "sf",         "matrixStats", "fasterize")
mbg_setup(package_list = pl, repos = core_repo)

# load in the region and run_date
load_from_parallelize()
reg             <- region
pathaddin <- paste0("<<<< FILEPATH REDACTED >>>>")

# load the image from launch
load("<<<< FILEPATH REDACTED >>>>")
load_from_parallelize() # in case stuff was in the image that shouldnt be

#reassigning because these get overwritten by the image loaded above
user            <- Sys.info()['user']
core_repo       <- sprintf("<<<< FILEPATH REDACTED >>>>") # main MBG repo
ig_repo         <- sprintf("<<<< FILEPATH REDACTED >>>>")

# run the u5m specific setup
preset_run_date <- run_date
source("<<<< FILEPATH REDACTED >>>>") 

load_from_parallelize() # in case stuff was in the image that shouldnt be

# clean year list
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

# Print some settings to console
for(thing in c('indicator','indicator_group','run_date','reg','pop_measure', 'ageasindic','abs','agegroups','admin_levels_to_aggregate',
               'get_summary_parameter_table'))
  message(paste0(thing,': ',paste0(get(thing),collapse=', ')))

#pull config file from model run in separate environment to avoid namespace pollution
e1 <- new.env()
load(paste0("<<<< FILEPATH REDACTED >>>>"), e1)
config <- get('config', e1)
modeling_shapefile_version <- get('modeling_shapefile_version', e1)
raking_shapefile_version <- get('raking_shapefile_version', e1)
rm(e1)
message(sprintf("Modeling shapefile version: %s", modeling_shapefile_version))
message(sprintf("Raking shapefile version: %s", raking_shapefile_version))

#get config arguments (character(0) if not in config, setting to NULL)
countries_not_to_rake <- config[V1 == "countries_not_to_rake", V2]
countries_not_to_subnat_rake <- config[V1 == "countries_not_to_subnat_rake", V2] 
rake_subnational <- config[V1 == "rake_subnational", V2] 
if(length(countries_not_to_rake) == 0) countries_not_to_rake <- c("ESH+GUF")
if(length(countries_not_to_subnat_rake) == 0) countries_not_to_subnat_rake <- c("PHL+NGA+PAK")
if(length(rake_subnational) == 0) rake_subnational <- T

## U5M HOTFIX
countries_not_to_rake <-'ESH+GUF'

################################################################################
## save nice table of fitted parameter values
if(get_summary_parameter_table==TRUE){
  mf <- readRDS("<<<< FILEPATH REDACTED >>>>")
  paramsum <- fitted_param_table_tmb(mf)
  write.csv(paramsum, file="<<<< FILEPATH REDACTED >>>>")
}

################################################################################
## Combine age bins

# Combine cell_pred age groups in each region
# NOTE: this is obviously sensitive to the age bins as data are prepped
# returns a list object, one for each age bin
ab_fracs_sub <- ab_fracs[agegroups]
abs_required <- unique(unlist(ab_fracs_sub))

cell_pred <- combine_agebins(age_bins    = abs_required,
                             reg         = reg,
                             ageasindic  = ageasindic,    
                             monthlyprob = FALSE, # doing binly now (did monthly in stage 1)
                             fractions   = ab_fracs_sub)
gc(full=TRUE)
# Note look in future package to make some this go quicker


################################################################################
## Load spatial info for this region

message('Loading spatial data')
adm0_list      <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
## Load simple raster from file
e2 <- new.env()
load("<<<< FILEPATH REDACTED >>>>", e2)
simple_raster <- get('simple_raster', e2)
rm(e2)

pixel_id <- seegSDM:::notMissingIdx(simple_raster)
pixel_spatial<-data.table(pixel_id=pixel_id)

################################################################################
## Pull population brick using new covariates function
pop_raster_annual <- load_and_crop_covariates_annual(covs           = 'worldpop_raked',                
                                                     measures       = pop_measure,     
                                                     simple_polygon = simple_raster,
                                                     start_year     = min(year_list),
                                                     end_year       = max(year_list),
                                                     interval_mo    = as.numeric(interval_mo),
                                                     agebin=1)[[1]]

pop_raster_annual <- crop(pop_raster_annual, extent(simple_raster))
pop_raster_annual <- setExtent(pop_raster_annual, simple_raster)
pop_raster_annual <- mask(pop_raster_annual, simple_raster)

# getting the pixels from the population raster that correspond to the valid pixels in the pop-raster
pop <- data.table(raster::extract(pop_raster_annual, pixel_id)) 
pop[,pixel_id:=pixel_id]
# Melting the dataframe such that it is long () and should match the number of rows in cell_pred
pop<-melt(pop,id.vars="pixel_id")
# Converting "worldpop.1" variables to actual years.
pop[, year := (min(year_list) - 1) + as.numeric(gsub("worldpop_raked.", "", variable))]
pop<-pop[,list(pixel_id,year,pop=value)]
# Setting values where pop is NA or 0 to 0.01 to avoid NAs and NaNs in aggregation
pop[is.na(pop),pop:=0.01]
pop[pop == 0 ,pop:=0.01]

################################################################################
## Loading and Assigning Admin Units to Pixels
message("Getting the spatial location (admin unit) of each of the pixel locations.")

adm0_list<-get_adm0_codes(reg, shapefile_version = modeling_shapefile_version) # Getting the adm0 GAUL codes, we can use this to make sure we don't accidentally include countries from buffers that shouldn't be in this region

message("Rasterizing shapefiles; this may take a while.")
admin_levels<-list() # Emtpy list of levels that will be filled with admin levels
for(lvl in c(0,1,2)){
  fieldname<-paste0("ADM",lvl,"_CODE")
  
  #load in link table to deal with pixel assignment to countries outside the region in rasterize_check_coverage
  #This is due to an edge case where a pixel is in 3 admin units between two countries, where the largest area fraction is in a country outside the region
  link <- get_link_table(simple_raster, modeling_shapefile_version)
  link_table <- link[["link_table"]]
  id_ras <- link[["pixel_ids"]]  
  
  link_table <- link_table[ADM0_CODE %in% adm0_list]
  
  sub_link_table <- link_table[, c("pixel_id", "area_fraction", fieldname), with = F]
  aggregated_link_table <- sub_link_table[, .(area_fraction = sum(area_fraction)), by = c("pixel_id", fieldname)]
  
  #read in rds version of admin shapefile to speed things up
  custom_shapefile <- readRDS(get_admin_shapefile(admin_level = lvl, version = modeling_shapefile_version, suffix = ".rds"))
  custom_shapefile_sub <- custom_shapefile[custom_shapefile$ADM0_CODE %in% adm0_list,]
  
  simple_polygon <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, tolerance = 0.4,
                                        shapefile_version = modeling_shapefile_version, 
                                        custom_shapefile = custom_shapefile_sub)
  subset_shape   <- simple_polygon[['subset_shape']]
  
  raster_list    <- build_simple_raster_pop(subset_shape, field = fieldname, link_table = aggregated_link_table)
  simple_raster_lvl  <- raster_list[['simple_raster']]
  
  pixel_spatial[[fieldname]]<-raster::extract(simple_raster_lvl,pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
  if(sum(is.na(pixel_spatial[[fieldname]]))>0){ # Check to see if any of the pixels don't have a location assigned
    message(paste0("   Whoah, there are some pixels that are NA, and have not been assigned a location for level ",lvl))
  }
}

pop<-merge(pop,pixel_spatial,by="pixel_id",all.x=T) # Merging on the spatial information to the population information.
pop<-pop[order(year,pixel_id)] # Re-ordering the pop object by year such that pixels ascent, and years ascend (same as cell_pred)


################################################################################
## Paralell over agebins, do raking and save some summaries

# make a function to aprallelize over age bins afterwards
abrake <- function(group){
  message(paste0(group,group,group,group,group))
  #aboutputlist <- list()
  
  # update the output directory for current age bin
  outdir <- sprintf("<<<< FILEPATH REDACTED >>>>")
  
  rake_to <- get_gbd_q(group, year_list)
  
  #get draw column names
  ndraws <- ncol(cell_pred[[sprintf('%s_all',group)]])
  overs <- paste0('V',1:ndraws)
  
  #if column names in draws are repeated, reset column names based on column number
  if(length(unique(colnames(cell_pred[[sprintf('%s_all',group)]]))) != length(colnames(cell_pred[[sprintf('%s_all',group)]]))){
    colnames(cell_pred[[sprintf('%s_all',group)]]) <- overs
  }
  
  #if countries not to rake are specified, sets "name" in gbd (raking_targets) to -1. 
  #the if_no_gbd parameter in rake_cell_pred will return unraked estimates if no targets are found for that country
  if(!is.null(countries_not_to_rake)) {
    connector <- get_gbd_locs(
      rake_subnational = F,
      reg = reg,
      shapefile_version = raking_shapefile_version
    )
    not_to_rake_adm_codes <- get_adm0_codes(countries_not_to_rake, shapefile_version=raking_shapefile_version)
    not_to_rake_gbd_loc_ids <- connector[ADM_CODE %in% not_to_rake_adm_codes, location_id]
    if(length(not_to_rake_gbd_loc_ids) > 0){
      setDT(rake_to)
      rake_to[name %in% not_to_rake_gbd_loc_ids, name := -1]
    }
  }
  
  #make custom raking shapefile if countries_not_to_subnat_rake is set in config
  if (rake_subnational & !is.null("countries_not_to_subnat_rake")){
    custom_raking_shapefile <- make_custom_raking_shapefile(countries_not_to_subnat_rake, raking_shapefile_version)
  } else {
    custom_raking_shapefile <- NULL
  }

  #run raking
  outputs <- rake_cell_pred(cell_pred = cell_pred[[sprintf('%s_all',group)]],
                            rake_to = rake_to,                 #dataframe or table with columns name, year, mean where name is the admin0/1 gbd location id
                            reg              = reg,
                            year_list        = year_list,
                            pop_measure      = "a0004t",             
                            rake_method      = "linear",   #linear or logit
                            rake_subnational = rake_subnational,    #T or F
                            field            = 'loc_id',
                            simple_raster    = simple_raster,
                            simple_polygon   = NULL, # Don't need this is pop_raster is passed
                            pop_raster       = pop_raster_annual,
                            if_no_gbd        = "return_unraked",
                            shapefile_path = get_admin_shapefile(admin_level = 0, raking = T, version = modeling_shapefile_version), 
                            modeling_shapefile_version = modeling_shapefile_version,
                            raking_shapefile_version   = raking_shapefile_version,
                            custom_raking_shapefile = custom_raking_shapefile,
                            countries_not_to_subnat_rake = countries_not_to_subnat_rake)
  
  raked_cell_pred <- outputs[["raked_cell_pred"]]
  rf <- outputs[["raking_factors"]]
  
  ## summarize raked predictions for each cell
  message('  Summarizing unraked cell preds:\nmean')
  unraked_raster <- 
    make_cell_pred_summary( draw_level_cell_pred = cell_pred[[sprintf('%s_all',group)]],
                            mask                 = simple_raster,
                            return_as_raster     = TRUE,
                            summary_stat         = 'mean')
  save_post_est(unraked_raster, 'raster', sprintf('%s_mean_unraked_%s_raster',reg,group), indic = sprintf('died_%s',group))
  
  
  message('  Summarizing raked cell preds:')
  for(summeasure in c('mean','lower','upper')){ # 'cirange',
    message(sprintf('    %s',summeasure))
    tmp <- 
      make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                              mask                 = simple_raster,
                              return_as_raster     = TRUE,
                              summary_stat         = summeasure)
    save_post_est(tmp, 'raster', sprintf('%s_%s_raked_%s_raster',reg,summeasure,group), indic = sprintf('died_%s',group))
    
  }
  
  #add on year, admin codes, and population to raked and unraked cell draws

  meta_raked_cell_pred <- cbind(raked_cell_pred, pop)
  meta_unraked_cell_pred <- cbind(cell_pred[[sprintf('%s_all',group)]], pop)
  
  for(ad in admin_levels_to_aggregate){
    message(sprintf('    Admin Level %i',ad))

    bycol <- c(paste0("ADM", ad, "_CODE"), "year")
    
    missing_ads_raked <- find_missing_adms(lvl = ad, simple_raster = simple_raster, reg = reg, adm0_list = adm0_list, raked_cell_pred, year_list)
    missing_ads_unraked <- find_missing_adms(lvl = ad, simple_raster = simple_raster, reg = reg, adm0_list = adm0_list, cell_pred[[sprintf('%s_all',group)]], year_list)
    
    #aggregate raked cell pred
    raked_adm <- meta_raked_cell_pred[ADM0_CODE %in% adm0_list,lapply(.SD,weighted.mean,w=pop, na.rm=T), by=bycol, .SDcols=overs]
    setnames(raked_adm,  paste0("ADM", ad, "_CODE"), "name")
  
    #aggregate unraked cell pred
    unraked_adm <- meta_unraked_cell_pred[ADM0_CODE %in% adm0_list,lapply(.SD,weighted.mean,w=pop, na.rm=T), by=bycol, .SDcols=overs]
    setnames(unraked_adm, paste0("ADM", ad, "_CODE"), "name")
    
    #add on missing admins if they are available to raked cell pred
    raked_missing <- get_missing_admin_cell_pred(cp = raked_cell_pred, 
                                                 mo = missing_ads_raked,
                                                 yl = year_list,
                                                 ad = ad)
    if(!is.null(raked_missing)){
      rownames <- rownames(raked_missing)
      raked_missing <- data.table(raked_missing,rownames)
      raked_missing[, c("name", "year") := tstrsplit(rownames, "_", fixed=TRUE)]
      raked_missing[,rownames := NULL]
      
      #drop aggregate-admins that were picked up by the find_missing_adms function above 
      raked_adm <- raked_adm[!(name %in% raked_missing$name),]
      raked_adm <- rbind(raked_adm, raked_missing)
    }
    #add on missing admins if they are available to unraked cell pred
    unraked_missing <- get_missing_admin_cell_pred(cp = cell_pred[[sprintf('%s_all',group)]], 
                                                   mo = missing_ads_unraked,
                                                   yl = year_list,
                                                   ad = ad)
    if(!is.null(unraked_missing)){
      rownames <- rownames(unraked_missing)
      unraked_missing <- data.table(unraked_missing,rownames)
      unraked_missing[, c("name", "year") := tstrsplit(rownames, "_", fixed=TRUE)]
      unraked_missing[,rownames := NULL]
      
      #drop aggregate-admins that were picked up by the find_missing_adms function above 
      unraked_adm <- unraked_adm[!(name %in% unraked_missing$name),]
      unraked_adm <- rbind(unraked_adm, unraked_missing)
    }
    
    save_post_est(raked_adm, 'csv', sprintf('%s_draws_raked_ad%i',reg,ad),   indic = sprintf('died_%s',group))
    save_post_est(unraked_adm,'csv', sprintf('%s_draws_unraked_ad%i',reg,ad), indic = sprintf('died_%s',group))

  }
  
  # save cell preds and rf file
  message('  saving cell preds')
  saveRDS(raked_cell_pred, file=sprintf("<<<< FILEPATH REDACTED >>>>")
  cell_draws <- cell_pred[[sprintf('%s_all',group)]]
  saveRDS(cell_draws, file=sprintf("<<<< FILEPATH REDACTED >>>>")
  rm(cell_draws)
  
  write.csv(rf, file=sprintf("<<<< FILEPATH REDACTED >>>>")

}
  

# if using singularity image, make sure MKL CORES is set to 1 for mclapply to work
lapply(
  as.list(agegroups),
  function(x) {
    abrake(x)
    gc(full=TRUE)
  }
)


################################################################################
## run death crosswalk

lapply(
  as.list(agegroups),
  function(x){
    run_q_to_d_crosswalk(run_date          = run_date,
                         region            = reg,
                         ages              = x,
                         shapefile_version = modeling_shapefile_version,
                         simple_raster     = simple_raster, # This function uses UNRAKED, so use old SR
                         pop_raster_brick  = pop_raster_annual)    
  }
)
