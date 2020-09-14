##############################################################################
## Clean Year 4 deliverables outputs and save in separate folder
##############################################################################

# (1) Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# load packages
library(data.table)
library(raster)

# set general arguments
core_repo       <- '<<< FILEPATH REDACTED >>>'
use_inla_country_fes <- FALSE
get_dalys <- F 

#indicate which group we are working on (lri | ort=diarrhea)
indicator_group <- 'lri'


if(indicator_group=='lri') {
  # set MBG arguments for LRI
  indicator       <- 'has_lri'
  run_date        <- '2020_06_11_11_19_26'
  regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name-ESH', 'dia_sssa',
               'dia_mcaca', 'dia_s_america_n', 'dia_s_america_s', 'dia_central_asia',
               'dia_se_asia', 'dia_malay', 'dia_south_asia-IND', 'dia_mid_east', 'dia_essa', 'IND','MNG')
  #regions <- 'IND'
  etiologies      <- c('lri_pneumo')
  group_measure <- c('prevalence', 'incidence', 'mortality')
  
} else {
  # set MBG arguments for diarrhea
  indicator       <- 'had_diarrhea'
  run_date        <- '2019_09_17_14_12_54'
  regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_sssa', 'dia_name', 'dia_s_america',
               'dia_mcaca', 'dia_central_asia', 'MNG', 'IND',
               'dia_se_asia', 'dia_malay', 'dia_south_asia-ind', 'dia_mid_east', 'dia_essa')
  #regions <- 'IND'
  group_measure <- c('prevalence', 'incidence', 'deaths')
}

#set up dirs
base_dir        <- '<<< FILEPATH REDACTED >>>'
share_dir       <- file.path(base_dir, '/output/', run_date)
save_dir        <- file.path(base_dir, 'deliverables/')
splits_dir      <- file.path(save_dir, 'age_sex_splits/')

# create directory to save outputs in
dir.create(save_dir, showWarnings = FALSE)

# set age and sex groups
group_ids <- c('a2_s1', 'a3_s1', 'a4_s1', 'a5_s1', 
               'a2_s2', 'a3_s2', 'a4_s2', 'a5_s2', 
               'a2_s3', 'a3_s3', 'a4_s3', 'a5_s3', 
               'a1_s1', 'a1_s2', 'a1_s3')

# Load MBG packages
package_list <- c(t(read.csv('<<< FILEPATH REDACTED >>>',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Load custom lbd core functions
lapply(paste0('<<< FILEPATH REDACTED >>>',
              list.files('<<< FILEPATH REDACTED >>>', pattern = 'functions')),
       source)

# Load custom post-estimation functions
lapply(paste0('<<< FILEPATH REDACTED >>>',
              list.files('<<< FILEPATH REDACTED >>>')),
       source)



# (2) Merge and save global rasters for age and sex groups ----------------------------------------------------

saveSplitRasters <- function(s, g, gm) {

  # read in mean raster files
  rasters <- list.files(splits_dir,
                          pattern = glob2rx(paste0('*', gm, '*', s, '*', g, '*.grd*')),
                          full.names=T) %>% 
    lapply(brick) %>% 
    do.call(raster::merge, .) %>% 
    writeRaster(.,
                paste0(save_dir, indicator, '_', s, '_', gm, '_', g, '_raster.tif'),
                overwrite = TRUE)
}

grid <- expand.grid(s=c('mean', 'upper', 'lower'), g=group_ids, gm=group_measure)
tmp <- mcmapply(saveSplitRasters, grid$s, grid$g, grid$gm, mc.cores=15)

# (3) Merge and save global rasters for under 5, both sex inc, mort, prev ---------------------------------------

saveRasters <- function(s, measure) {

# read in mean raster files
rasters <- list.files(share_dir,
                      pattern = glob2rx(paste0('*', s, '*', measure, '*.grd*')),
                      full.names=T) %>% 
  lapply(brick) %>% 
  do.call(raster::merge, .) %>% 
  writeRaster(.,
              paste0(save_dir, indicator, '_',  ifelse(s == 'prediction', 'mean', s), '_', measure, '_raster.tif'),
              overwrite = TRUE)

}

grid <- expand.grid(s=c('prediction', 'upper', 'lower'), measure=group_measure)
tmp <- mcmapply(saveRasters, grid$s, grid$measure, mc.cores=9)


# (4) Merge and save dalys rasters for under 5, both sex ---------------------------------------
#TODO we dont need etis/dalys anymore?
if(get_dalys) {
  # loop over stats
  saveRasters <- function(s) {

    # read in mean raster files
    rasters <- list.files(share_dir,
                          pattern = glob2rx(paste0('*', s, '*daly*.grd*')),
                          full.names=T) %>%
      lapply(brick) %>%
      do.call(raster::merge, .) %>%
      writeRaster(.,
                  paste0(save_dir, indicator, '_',  ifelse(s == 'prediction', 'mean', s),  '_daly_raster.tif'),
                  overwrite = TRUE)
  }

  tmp <- mclapply(c('prediction', 'upper', 'lower'), saveRasters, mc.cores=3)

  # (5) Merge and save global rasters for under 5, both sex etiologies ---------------------------------------
  
  # loop over etiologies
  saveEtiRasters <- function(s, eti) {

      # read in mean raster files
      rasters <- list.files(share_dir,
                            pattern = glob2rx(paste0('*', s, '*mortality*.grd*')),
                            full.names=T) %>%
        lapply(brick) %>%
        do.call(raster::merge, .) %>%
        writeRaster(.,
                    paste0(save_dir, indicator, '_',  ifelse(s == 'prediction', 'mean', s),  '_mortality_raster.tif'),
                    overwrite = TRUE)

  }
  grid <- expand.grid(s=c('prediction', 'upper', 'lower'), eti=etiologies)
  mcmapply(saveEtiRasters, grid$s, grid$eti, mc.cores=3)
  
  # (6) Clean and move admin csvs ---------------------------------------
  
  #prev
  #inc
  #mort
  #daly
  #etiologies
  
  #combine and summarize daly admin draws
  message(paste0('Combining raked aggregated results for ', measure))
  combine_aggregation(rd       = run_date,
                      indic    = indicator,
                      ig       = indicator_group,
                      ages     = 0,
                      regions  = regions,
                      holdouts = 0,
                      raked    = T,
                      measure  = 'daly',
                      delete_region_files = F,
                      metrics = c('rates','counts'))

  # summarize admins
  summarize_admins(ad_levels = c(0,1,2), raked = T, measure = 'daly', metrics = c('rates','counts'))

  #combine and summarize etiology admin draws
  message(paste0('Combining raked aggregated results for ', paste(etiologies, collapse = ', '), ' mortality'))
  use_inla_country_fes <- FALSE
  combine_aggregation(rd       = run_date,
                      indic    = indicator,
                      ig       = indicator_group,
                      ages     = 0,
                      regions  = regions,
                      holdouts = 0,
                      raked    = T,
                      measure  = 'mortality',
                      delete_region_files = F,
                      metrics = c('rates','counts'),
                      eti = etiologies)

  # summarize admins
  summarize_admins(ad_levels = c(0,1,2), raked = T, measure = 'mortality', metrics = c('rates','counts'), eti = etiologies)
  
}

#move files to y4 deliverables folder
files_to_copy <- list.files(paste0(share_dir, 'pred_derivatives/admin_summaries/'), pattern = '_raked_', full.names=T)
file.copy(from = files_to_copy, to = save_dir, overwrite = TRUE)

#age and sex
for (a in 0:2) {
  files_to_copy <- list.files(paste0(splits_dir), pattern = paste0('adm', a, '_summary'), full.names=T)
  file.copy(from = files_to_copy, to = save_dir, overwrite = TRUE)
}
