###############################################################################
## MBG launch script
###############################################################################
## Setup -------------------------------------------------------------------------

## Clear environment
rm(list=ls())

#REDACTED

#load external packages
pacman::p_load(assertthat, magrittr)

#detect if running interactively
interactive <- F  %>% #manual override
  ifelse(., T, !length(commandArgs())>2) %>%  #check length of arguments being passed in
  ifelse(., T, !(is.na(Sys.getenv("RSTUDIO", unset = NA)))) #check if IDE

## if running interactively, set arguments
if (interactive) {
  warning('interactive is set to TRUE - if you did not mean to run MBG interactively then kill the model and set interactive to FALSE in parallel script')

  #REDACTED

  
} else {

  ## Set repo location, indicator group, and some arguments
  user            <- commandArgs()[4]
  core_repo       <- commandArgs()[5]
  indicator_group <- commandArgs()[6]
  indicator       <- commandArgs()[7]
  config_par      <- commandArgs()[8]
  cov_par         <- commandArgs()[9]
  reg         <- commandArgs()[10]
  parallel_script <- commandArgs()[11]
  
  ## Set some covariate options
  plot_covariates <- commandArgs()[12]
  covariate_plotting_only <- commandArgs()[13]
  
  ## Set to prod or geos nodes
  proj_arg <- commandArgs()[14]
  use_geos_nodes <- commandArgs()[15]
  
  ## Set run date
  run_date <- commandArgs()[16]
  
  ## Define personal repo
  my_repo <- commandArgs()[17]

}
  
## Singularity version
sing_dir <- 'REDACTED'
which_sing <- file.path(sing_dir, '#REDACTED')

## Released Singularity imgs
sing_dir <- 'REDACTED'
which_sing <- file.path(sing_dir, '#REDACTED')

## Load MBG packages
package_list <- c(t(read.csv(paste0(core_repo, '/#REDACTED/package_list.csv'), header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
  mbg_setup(package_list = package_list, repos = core_repo)
  
fix_diacritics <<- function(x) {
  
  require(mgsub)
  
  #first define replacement patterns as a named list
  defs <-
    list('??'='S', '??'='s', '??'='Z', '??'='z', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', 
         '??'='C', '??'='E', '??'='E','??'='E', '??'='E', '??'='I', '??'='I', '??'='I', '??'='I', '??'='N', '??'='O', 
         '??'='O', '??'='O', '??'='O', '??'='O', '??'='O', '??'='U','??'='U', '??'='U', '??'='U', '??'='Y', '??'='B', 
         '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='c','??'='e', '??'='e', '??'='e', 
         '??'='e', '??'='i', '??'='i', '??'='i', '??'='i', '??'='o', '??'='n', '??'='o', '??'='o', '??'='o', '??'='o',
         '??'='o', '??'='o', '??'='u', '??'='u', '??'='u', '??'='y', '??'='y', '??'='b', '??'='y', '??'='Ss')
  
  #then force conversion to UTF-8 and replace with non-diacritic character
  enc2utf8(x) %>% 
    mgsub(., pattern=enc2utf8(names(defs)), replacement = defs) %>% 
    return
  
}

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')
config <- set_up_config(repo            = my_repo,
                        indicator_group = indicator_group,
                        indicator       = indicator,
                        config_name     = paste0('/model/configs/config_', config_par),
                        covs_name       = paste0('/model/configs/covs_', cov_par)
                        )

## Create output folder with the run_date
outputdir      <- paste('#REDACTED', indicator_group, indicator, 'output', run_date, '', sep='/')
  dir.create(outputdir, recursive=T)

## save config type
file.create(paste0(outputdir, '/000_', config_par, '.note'))

## Create directory structure for this model run
create_dirs(indicator_group = indicator_group, indicator = indicator)

## Create proper year list object
if (class(year_list) == 'character') year_list <- eval(parse(text=year_list))

## Ensure you have defined all necessary settings in your config
check_config()

## If running individual countries make sure all country FEs and REs off
individual_countries <- ifelse(nchar(reg) == 3, TRUE, FALSE)
if (individual_countries) {
  use_child_country_fes <- FALSE
  use_inla_country_fes  <- FALSE
  use_inla_country_res  <- FALSE
}

## If not running a model with India make sure all subnational REs off
if (reg != 'dia_south_asia' & reg != 'IND') {
  use_subnat_res <- FALSE
}

## Make holdouts -------------------------------------------------------------------------

if(as.logical(makeholdouts) & !as.logical(skiptoinla)){
  message('Making holdouts')
  
  set.seed(98118)
  
  # load the full input data
  df <- load_input_data(indicator   = indicator,
                        years       = yearload,
                        yl          = year_list,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share))
  
  # add in location information
  df <- merge_with_ihme_loc(df, re = reg, shapefile_version = modeling_shapefile_version)
  
  # remove data for countries outside of region
  df <- df[!is.na(region)]
  
  # make holdouts based on nids, if selected
  if (holdout_strategy == 'nids') {
    
    # make sure we have enough nids to do holdouts
    if (length(unique(df$nid)) < 5) n_ho_folds <- length(unique(df$nid))

    # make a list of dfs for each region, with 5 nid folds identified in each
    stratum_ho <-   make_folds(data       = df,
                               n_folds    = as.numeric(n_ho_folds),
                               spte_strat = 'nids',
                               strat_cols = 'region',
                               ss_col     = ss_col,
                               yr_col     = yr_col,
                               seed       = 98112,
                               save.file = paste0('REDACTED',
                                                  indicator_group, '/',
                                                  indicator, '/output/',
                                                  run_date, '/stratum_', reg, '.rds'))
  } 
  
  # make holdouts based on admin 1, if selected
  if (holdout_strategy == 'admin') {
    
    # Load simple polygon template by region
    gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, shapefile_version = modeling_shapefile_version)
    subset_shape        <- simple_polygon_list[[1]]
    
    # Load simple raster by region
    raster_list        <- build_simple_raster_pop(subset_shape)
    simple_raster      <- raster_list[['simple_raster']]
    
    # load admin 2 shapefile and crop by region
    shapefile_admin <- shapefile(get_admin_shapefile(admin_level = 2, version = modeling_shapefile_version))
    shapefile_admin <- gBuffer(shapefile_admin, byid = TRUE, width = 0)
    shapefile_admin <- crop(shapefile_admin, extent(subset_shape))
    
    # load mask raster and crop by region
    raster_mask <- raster('REDACTED')
    raster_mask <- crop(raster_mask, extent(simple_raster))
    
    # make admin 2 holdouts
    stratum_ho <- make_folds(data = df,
                             n_folds = as.numeric(n_ho_folds),
                             spat_strat = 'poly',
                             temp_strat = 'prop',
                             strat_cols = 'region',
                             ts = as.numeric(ho_ts),
                             mb = as.numeric(ho_mb),
                             admin_shps = shapefile_admin,
                             admin_raster = simple_raster,
                             mask_shape = subset_shape,
                             mask_raster = raster_mask,
                             shape_ident = 'ADM2_CODE',
                             lat_col = lat_col,
                             long_col = long_col,
                             ss_col = ss_col,
                             yr_col = yr_col,
                             seed = 98112,
                             save.file = paste0('#REDACTED/',
                                                indicator_group, '/',
                                                indicator, '/output/',
                                                run_date, '/stratum_', reg, '.rds'))
  } 
  
}

## Launch parallel script -------------------------------------------------------------------------

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(as.logical(makeholdouts)) loopvars <- expand.grid(reg, 0, 0:as.numeric(n_ho_folds)) else loopvars <- expand.grid(reg, 0, 0)

## loop over them, save images and submit qsubs
for(i in 1:nrow(loopvars)){
  
  message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))
  
  # check that holdouts exist
  if (!as.logical(skiptoinla)) {
    if (loopvars[i, 3] > 0) {
      count <- sum(stratum_ho[[paste0('region__', loopvars[i, 1])]]$fold %in% loopvars[i, 3])
      if (count == 0) next
    }
  }
  
  # set cores by region
  r <- as.character(loopvars[i,1])
  region_cores <- 8
  if(r == 'NGA' | r == 'PAK' | r == 'KEN' | r == 'ZWE' | r == 'VNM' | r == 'THA' | r == 'dia_central_asia') region_cores <- 4
  if(r == 'dia_se_asia-VNM-THA' | r == 'MNG' | r == 'COD') region_cores <- 6
  if(r == 'dia_malay' | r == 'dia_name' | r == ' '| r == 'dia_cssa'| r == 'dia_se_asia') region_cores <- 10
  if(r %like% 'dia_south_asia' | r %like% 'dia_essa' | r=='CHN' | r=='dia_s_america-BRA' | r == 'dia_afr_horn') region_cores <- 12
  if(r == 'dia_s_america' | r == 'dia_wssa'| r == 'dia_chn_mng' | r == 'dia_wssa-NGA' | r=='BRA') region_cores <- 25
  if(loopvars[i, 3] > 0) region_cores <- round(region_cores*0.8)
  
  # convert cores approximately to memory and run time
  region_rt <- '06:00:00:00'
  if (region_cores < 9) region_rt <- '04:00:00:00'
  if (region_cores < 6) region_rt <- '03:12:00:00'
  if (region_cores > 9) region_rt <- '16:00:00:00'
  region_mem <- region_cores*35
  
  # set thread options
  threads <- ifelse(region_cores>10, 3, 6)
  
  # make a qsub string
  qsub <- make_qsub_share(age           = loopvars[i,2],
                          reg           = as.character(loopvars[i,1]),
                          holdout       = loopvars[i,3],
                          test          = F,
                          indic         = indicator,
                          saveimage     = TRUE,
                          memory        = region_mem,
                          cores         = threads,
                          run_time      = region_rt,
                          proj          = proj_arg,
                          geo_nodes     = as.logical(use_geos_nodes),
                          corerepo      = core_repo,
                          code          = parallel_script,
                          addl_job_name = paste0('orbit_',indicator_group, '_', as.character(loopvars[i,1])),
                          singularity   = which_sing,
                          singularity_opts = list(SET_OMP_THREADS=threads, SET_MKL_THREADS=threads))
  
  message(qsub) #for posterity
  
  # submit job
  system(qsub)
  
}
