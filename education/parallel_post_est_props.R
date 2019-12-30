core_repo <- as.character(commandArgs()[4])
indic_repo <- as.character(commandArgs()[5])
age_sex <- as.character(commandArgs()[6])
rd<- as.character(commandArgs()[7])
use_raked<- as.numeric(commandArgs()[8])
r<-as.character(commandArgs()[9])

reg <- r

commondir      <- sprintf("<<<< FILEPATH REDACTED >>>>")

## Load libraries and  MBG project functions.
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

setwd(core_repo)


# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)
load_R_packages('assertthat')


## Set up rescaled output directories (new indicators)
for(p in c('zero','no_primary','primary','secondary')) {
  dir <- paste0("<<<< FILEPATH REDACTED >>>>", 'edu_', p, '_prop_', age_sex, '_rescaled/')
  dir.create(paste0(dir), showWarnings = FALSE)
  dir.create(paste0(dir, '/output'), showWarnings = FALSE)
  dir.create(paste0("<<<< FILEPATH REDACTED >>>>"), showWarnings = FALSE)
}
  
## Get region_list for this model set.
config <- load_config(repo = indic_repo,
                      core_repo = core_repo,
                      indicator_group = 'education',
                      indicator = paste0('edu_zero_prop_', age_sex),
                      post_est_only = TRUE,
                      run_date = paste0("<<<< DATE REDACTED >>>>"))

## Load simple polygon template to model over
gaul_list           <- get_adm0_codes(r, shapefile_version = "<<<< DATE REDACTED >>>>")
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T, shapefile_version = '2019_02_27')
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
raster_list        <- build_simple_raster_pop(subset_shape)
simple_raster      <- raster_list[['simple_raster']]
pop_raster         <- raster_list[['pop_raster']]

## Load zero proportion draws for this region.
message(paste0('Rescaling draws for ', r, '...'))
if(use_raked==1) {
  load"<<<< FILEPATH REDACTED >>>>")
  zero_prop_draws <- raked_cell_pred
  raked_cell_pred <- NULL
}
if(use_raked==0) {
  load("<<<< FILEPATH REDACTED >>>>")
  zero_prop_draws <- cell_pred
  cell_pred <- NULL
}
gc(T)

## Load no primary proportion draws for this region (proportion of those with no primary OF THOSE who don't have zero)
message('Rescaling no primary...')
if(use_raked==1) {
  load("<<<< FILEPATH REDACTED >>>>")
  no_primary_prop_draws <- raked_cell_pred
  raked_cell_pred <- NULL
}
if(use_raked==0) {
  load("<<<< FILEPATH REDACTED >>>>")
  no_primary_prop_draws <- cell_pred
  cell_pred <- NULL
}
gc(T)

## Rescale no primary proportion from out of those who don't have zero to be out of TOTAL POPULATION.
no_primary_prop_draws <- no_primary_prop_draws * (1 - zero_prop_draws)

## Load primry proportions draws for this region (proportion of those with only primary school OF THOSE who don't have zero or no primary).
message('Rescaling primary...')
if(use_raked==1) {
  load("<<<< FILEPATH REDACTED >>>>")
  primary_prop_draws <- raked_cell_pred
  raked_cell_pred <- NULL
}
if(use_raked==0) {
  load("<<<< FILEPATH REDACTED >>>>")
  primary_prop_draws <- cell_pred
  cell_pred <- NULL
}
gc(T)

## Rescale no primary proportion from out of those who don't have zero or no primary to be out of TOTAL POPULATION.
primary_prop_draws = primary_prop_draws * (1 - zero_prop_draws - no_primary_prop_draws)

## Calculated secondary proportion out of total population as the complement of the sum of all of our rescaled proportions.
message('Rescaling secondary...')
secondary_prop_draws = (1 - zero_prop_draws - no_primary_prop_draws - primary_prop_draws)


create_dirs('education', paste0('edu_zero_no_primary_prop_', age_sex))
sharedir <- paste0("<<<< FILEPATH REDACTED >>>>")
dir.create(sharedir)

zero_no_primary_prop_draws <- zero_prop_draws + no_primary_prop_draws

## Save draws for the three proportions we had to rescale (zero_prop is already in terms of the total population, so we can just pull from that folder).
message(paste0('Saving rescaled draws for ', r, '...'))
if(use_raked==1) raked_file <- 'raked'
if(use_raked==0) raked_file <- 'unraked'
save(
  no_primary_prop_draws,
  file = ("<<<< FILEPATH REDACTED >>>>"),
  compress = TRUE
)

save(
  primary_prop_draws,
  file = ("<<<< FILEPATH REDACTED >>>>"),
  compress = TRUE
)

save(
  secondary_prop_draws,
  file = ("<<<< FILEPATH REDACTED >>>>"),
  compress = TRUE
)

save(
  zero_no_primary_prop_draws,
  file = ("<<<< FILEPATH REDACTED >>>>"),
  compress = TRUE
)


summstats <- c('mean', 'upper', 'lower')

raked <- 'raked'
indicator_group <- 'education'
run_date <- rd
save_cell_pred_summary <- function(summstat, cpred, indic, ...) {
  message(paste0('Making unraked summmary raster for: ',summstat, " (", raked, ")"))
  ras <- make_cell_pred_summary(
    draw_level_cell_pred = get(cpred),
    mask                 = simple_raster,
    return_as_raster     = TRUE,
    summary_stat         = summstat,
    ...)
  indicator <- indic
  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ""),'_',summstat,'_raster'), indic = indicator)
}

# Do this as lapply to not fill up memory in global env with big obs
cpred_list <- c('zero_prop_draws', 'no_primary_prop_draws', 'primary_prop_draws', 'secondary_prop_draws', 'zero_no_primary_prop_draws')
summ_list <- expand.grid(summstats[summstats != "p_below"], cpred_list)
summ_list[,3] <-  rep(c(paste0(c('edu_zero_prop_', 'edu_no_primary_prop_', 'edu_primary_prop_', 'edu_secondary_prop_', 'edu_zero_no_primary_prop_'), age_sex)), each =3)
lapply(1:nrow(summ_list), function(i) {
  summstat <- as.character(summ_list[i, 1])
  cpred <- as.character(summ_list[i, 2])
  indic <- as.character(summ_list[i, 3])
  save_cell_pred_summary(summstat, cpred, indic)
})


