#####################################################################
## Generic parallel script for running MBG models                  ##
#####################################################################


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
indis <- c('s_piped', 's_imp_cr', 's_unimp_cr',
  'w_piped', 'w_imp_cr', 'w_unimp_cr')

for (ii in indis) {
  message(ii)
  # grab arguments
  # note this requires a shell script with "<$1 --no-save $@", because its starting at 4
  reg <- 'ind'
  age <- 0
  run_date <- 'ind_custom'
  test <- 0
  holdout <- 0
  indicator <- ii
  indicator_group <- 'wash'

  # make a pathaddin that get used widely
  pathaddin <- paste0("_bin", age, "_", reg, "_", holdout)

  # load an image of the main environment
  load("<<<< FILEPATH REDACTED >>>>")

  # In case anything got overwritten in the load, reload args
  reg <- 'ind'
  age <- 0
  run_date <- 'ind_custom'
  test <- 0
  holdout <- 0
  indicator <- ii
  indicator_group <- 'wash'
  pathaddin <- paste0("_bin", age, "_", reg, "_", holdout)
  outputdir <- "<<<< FILEPATH REDACTED >>>>"
  dir.create(outputdir, showWarnings = FALSE)

  # print run options
  message("options for this run:\n")
  for (arg in c(
    "reg", "age", "run_date", "test", "holdout",
    "indicator", "indicator_group", "pathaddin", "outputdir"
  )) {
    message(paste0(arg, ":  ", get(arg), "   // type: ", class(get(arg))))
  }

  ## load packages and custom functions
  ## drive locations
  root <- "<<<< FILEPATH REDACTED >>>>"
  sharedir <- "<<<< FILEPATH REDACTED >>>>"
  commondir <- "<<<< FILEPATH REDACTED >>>>"
  package_list <- c(t(read.csv(sprintf("%s/package_list.csv", commondir), header = FALSE)))

  setwd(repo)
  core_repo <- repo

  for (p in package_list) {
    try(library(p, character.only = T))
  }

  library(seegSDM)
  library(seegMBG)
  library(mgcv)
  library(sf)
  library(raster)

  # Load MBG packages and functions
  message("Loading in required R packages and MBG functions")
  source(paste0(repo, "/mbg_central/setup.R"))
  mbg_setup(package_list = package_list, repos = repo)

  ## Throw a check for things that are going to be needed later
  message("Looking for things in the config that will be needed for this script to run properly")

  # Ensure raster library is loaded and configure it.
  # print out session info so we have it on record
  sessionInfo()
  setwd(repo)

  # Throw a check for things that are going to be needed later
  message("Looking for things in the config that will be needed for this script to run properly")
  holdout_strategy <- "prop"
  individual_countries <- "NA"
  jn <- "NA"
  metric_space <- "NA"
  rake_countries <- FALSE
  rake_transform <- "logit"
  skip.inla <- FALSE
  skip.stacking <- FALSE
  spat_strat <- "prop"
  temp_strat <- "rand"
  use_geos_nodes <- TRUE
  test <- FALSE
  check_config()
  shapefile_version <- modeling_shapefile_version
  # cores to use
  cores_to_use <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
  reg_list <- readRDS(paste0(repo, "wash/00_reg_list.rds"))


  ## Load simple polygon template to model over
  simple_polygon_list <- load_simple_polygon(
    buffer = 1, tolerance = 0.4,
    shapefile_version = modeling_shapefile_version,
    custom_shapefile_path = "<<<< FILEPATH REDACTED >>>>"
  )
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]

  ## Load list of raster inputs (pop and simple)
  raster_list <- build_simple_raster_pop(subset_shape, pop_release = pop_release)
  simple_raster <- raster_list[["simple_raster"]]
  pop_raster <- raster_list[["pop_raster"]]

  year_list <- 2000:2018

  ## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
  ## For holdouts, we have depreciated val level, so each val level must be recorded in a different run date
    df <- load_input_data(
      indicator = gsub(paste0("_age", age), "", indicator),
      simple = subset_shape,
      agebin = age,
      removeyemen = FALSE,
      pathaddin = pathaddin,
      years = yearload,
      withtag = as.logical(withtag),
      datatag = datatag,
      use_share = as.logical(use_share),
      yl = year_list
    )

  year_list <- 2000:2018

  # Filter data.frame based on quality checks
  if (use_mean_year) {
    message("Used mean year per NID!")
    df <- as.data.frame(df)
    df <- df %>%
      dplyr::group_by(nid) %>%
      dplyr::mutate(year = weighted.mean(x = year, w = weight*N)) %>%
      dplyr::ungroup()
  }

  reg_list <- readRDS(paste0(repo, "wash/00_reg_list.rds"))
  if (length(reg_list[[reg]]) == 1) {
    use_inla_country_res <- FALSE
    message("set INLA country REs to FALSE since region only has 1 country!")
  }
  df <- as.data.frame(df)
  df$N <- round(df$N, digits = 7)
  df[, indicator] <- round(df[, indicator], digits = 7)
  df <- filter(df, N > 0)
  df <- filter(df, year >= min(year_list), year <= max(year_list))
  df <- filter(df, country %in% c('IND', 'PAK'))
  df <- filter(df, !is.na(latitude), !is.na(longitude))
  df$year <- round(df$year)

  # set 0s to 0.001 and 1s to 0.999
  df[,indicator] <- ifelse(df[,indicator]/df$N == 0, df$N*rnorm(length(df$N), 0.001, 0.0001), df[,indicator])
  df[,indicator] <- ifelse(df[,indicator]/df$N == 1, df$N*rnorm(length(df$N), 0.999, 0.0001), df[,indicator])

  # recalculate prop
  df$prop <- df[, indicator] / df$N

  ## Remove outliers
  df <- df[which(df[, indicator] <= df$N), ]

  df <- as.data.table(df)
  print(head(df))
  write.csv(df, file = "<<<< FILEPATH REDACTED >>>>", row.names = FALSE)
  print('########################################################################')
  print('########################################################################')
}
