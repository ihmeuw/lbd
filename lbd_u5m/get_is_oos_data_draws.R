# load in the arguments
user            <- Sys.info()['user']
core_repo       <- '<<<< FILEPATH REDACTED >>>>' # main MBG repo
ig_repo         <- '<<<< FILEPATH REDACTED >>>>'      # u5m specific repo
indicator       <- 'died'
indicator_group <- 'u5m'

# load central lbd code stuff
source(paste0(core_repo, '/mbg_central/setup.R'))
pl <- c("rgeos",      "data.table", "raster",     "rgdal",      "INLA",
        "seegSDM",    "seegMBG",    "dismo",      "gbm",        "foreign",
        "parallel",   "doParallel", "grid",       "gridExtra",  "pacman",
        "gtools",     "glmnet",     "ggplot2",    "RMySQL",     "plyr",
        "tictoc",     "dplyr",      "magrittr",   "tidyr",      "sp",
        "sf",         "matrixStats")
mbg_setup(package_list = pl, repos = core_repo)

# load in the region and run_date
load_from_parallelize()

# load the image from launch
load('<<<< FILEPATH REDACTED >>>>')
load_from_parallelize() # in case stuff was in the image that shouldnt be

# run the u5m specific setup
preset_run_date <- run_date
source(sprintf('%s/setup.R',ig_repo)) 

load_from_parallelize() # in case stuff was in the image that shouldnt be

# clean year list
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

# Print some settings to console
for(thing in c('region','agebin','run_date','year_list','ab_fracs','makeholdouts'))
  message(paste0(thing,': ',paste0(get(thing),collapse=', ')))

# TESTING ONLY!
#modeling_shapefile_version = 'current'

run_in_oos <- get_is_oos_draws_u5m(ind_gp        = indicator_group,
                                   ind           = indicator,
                                   rd            = run_date,
                                   region        = region,
                                 #  regions       = region,
                                   ind_fm        = 'binomial',
                                   abfracs       = ab_fracs,
                                   agebin        = agebin,
                                   yrs           = year_list,
                                   fold_table    = oos_fold_table,
                                   data_sample   = data_sample,
                                  # get.oos       = FALSE,
                                   year_col      = 'year',
                                   write.to.file = TRUE)

for(i in 1:10)
  message('~~~~~~~~~~~~~~~~~~~~~~~~ FINISHED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
