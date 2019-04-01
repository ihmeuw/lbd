# HEADER ------------------------------------------------------------------
# misc_vaccine_functions.R
# Purpose: Various functions for vaccine coverage modeling
#**************************************************************************


# FUNCTIONS ---------------------------------------------------------------

# Timer Functions ---------------------------------------------------------

## Functions to time within the child stacking regions
##
## General usage:
##    require(tictoc)
## 
##    tic("Step 1")
##    **your code here**
##    toc(log = T)
##
##    ticlog <- tic.log(format = F)
##    generate_time_log(ticlog)   
##
##  Returns: data table with two columns
##     "step": names of events (e.g. "Step 1")
##     "time": time elapsed (as text: Xh Xm Xs)
##
##  Note: can nest tic/toc pairs

generate_time_log <- function(ticlog) {
  
  # Set up packages
  require(magrittr)
  require(data.table)
  
  # Functions in functions
  strip_time <- function(x) {
    sec <- as.numeric(x$toc - x$tic)
    time <- format_time(sec)
    name <- x$msg
    
    df <- c(name, time) %>%
            t %>%
            as.data.table
    
    names(df) <- c("step", "time")
    
    return(df)
  }
  
  format_time <- function(run_time) {
    run_time <- round(as.numeric(run_time),0)
    
    hours <- run_time %/% 3600
    remainder <- run_time %% 3600
    minutes <- remainder %/% 60 
    seconds <- remainder %% 60
    
    run_time <- paste0(hours, "h ", minutes, "m ", seconds, "s")
    return(run_time)
  }
  
  df_out <- lapply(ticlog, strip_time) %>% rbindlist
  
  return(df_out)

}

# Wrapper for easy plotting of run

plot_run_raster <- function(raster, 
                           vax_title,
                           summary_stat,
                           stat_title,
                           indicator = indicator, 
                           indicator_group = indicator_group, 
                           run_date = run_date,
                           layer_names = c(2000:2016)) {

  out_dir <- "<<<< FILEPATH REDACTED >>>>"
  dir.create(out_dir)

  plot_raster(ras = raster,
              indicator = paste0(indicator, "_", summary_stat),
              region = 'africa', 
              return_map = F,
              out_dir = out_dir,
              highisbad = F,
              min_value = 0,
              mid_value = 0.5,
              max_value = 1,
              legend_title = paste0(vax_title, "\n", stat_title),
              plot_title = paste0(vax_title, " ", stat_title),
              layer_names = as.character(layer_names),
              cores = 16,
              individual_layers = T)

}

# Create raster of GBD estimates

make_gbd_rasterbrick <- function(gbd_data, 
                                 indicator, 
                                 gaul_list,
                                 years = c(2000, 2005, 2010, 2015)) {

  message(paste0("Creating GBD RasterBrick for ", indicator))

  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                        buffer = 0.4)

  subset_shape   <- simple_polygon_list[[1]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]

  # create template raster
  ext <- extent(simple_raster)
  ncol <- ncol(simple_raster)
  nrow <- nrow(simple_raster)
  r <- raster(ext, ncol, nrow)

  gbd_brick_list <- lapply(years, function(yr) {

    my_spdf <- merge(subset_shape, gbd_data[year == yr], by.x = "GAUL_CODE", by.y = "name", all.x = T, all.y = F)
    ras <- rasterize(my_spdf, r, "mean")

  })

  gbd_brick <- brick(gbd_brick_list)

  return(gbd_brick) 
  
}

# Plot raking factors
make_rf_plot <- function(indicator_group, indicator, run_date) {

  require(data.table)
  require(ggplot2)
  require(magrittr)

  in_dir <- "<<<< FILEPATH REDACTED >>>>"
  out_dir <- paste0(in_dir, 'plots/')
  rf_file <- paste0(in_dir, indicator, "_rf.csv")

  table_file <- "<<<< FILEPATH REDACTED >>>>/gaul_list.csv"
  gaul_table <- read.csv(table_file) %>% data.table
  gaul_table <- subset(gaul_table, select = c("short_name", "gaul"))

  # fix Sudan
  gaul_table[short_name == "Sudan", gaul := 6 ]

  rf <- read.csv(rf_file, stringsAsFactors = F) %>% as.data.table
  rf[, X := NULL]
  setnames(rf, "name", "gaul")

  rf <- merge(rf, gaul_table, all.x = T)

  regions <- c("cssa", "essa", "sssa", "wssa", "name")

  for (reg in regions) {
    
  rf_reg <- rf[gaul %in% get_gaul_codes(reg)]

  gg <- ggplot(rf_reg, aes(x = rake_to_mean, y = geo_mean)) +
    geom_point(aes(color = year)) +
    xlim(c(0,1)) + ylim(c(0,1)) +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face="bold", size=10)) +
    labs(x = "GBD predicted", y = "MBG predicted", color = "Year") +
    facet_wrap(~short_name)

  png(filename = paste0(out_dir, "rf_scatter_", reg, ".png"),
      width = 10, 
      height = 6,
      units = "in", 
      res = 300,
      type = "cairo")
    
  print(gg)
    
  dev.off()

  }

}

# Transform spatial effects to more meaningful scales
transform_spatial_effects <- function(params.mat,age){

    nr=nrow(params.mat)

    tau_mean =exp(params.mat[(nr-2),1])
    kappa_mean=exp(params.mat[(nr-1),1])

    tau_lb =exp(params.mat[(nr-2),3])
    kappa_lb=exp(params.mat[(nr-1),3])

    tau_ub =exp(params.mat[(nr-2),5])
    kappa_ub=exp(params.mat[(nr-1),5])

    res=data.table('parameter'=c("GPRandom Range in Decimal Degrees","GPRandom Nominal Variance"))
        res[,paste0('a',age,'_mean'):=c(sqrt(8)/kappa_mean,1/(4*pi*kappa_mean^2*tau_mean^2))]
        res[,paste0('a',age,'_lb')  :=c(sqrt(8)/kappa_ub,1/(4*pi*kappa_ub^2*tau_ub^2))]
        res[,paste0('a',age,'_ub')  :=c(sqrt(8)/kappa_lb,1/(4*pi*kappa_lb^2*tau_lb^2))]

    return(res)

  }

# Create a clean table of model results (hyperparameter posteriors, beta coefficients, etc.)
clean_model_results <- function(rd   = run_date,
                                regs = Regions,
                                ages = 0,
                                nm   = '',
                                indicator = indicator, 
                                indicator_group = indicator_group){

  require(magrittr)

  sharedir <- paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', indicator, '/output/', rd, '/')

  # make loopvars
  lv <- expand.grid(regs,ages)

  # grab model fit objects
  mods <- model_fit_table(lv=lv,rd=rd,nullmodel=nm, indicator = indicator, indicator_group = indicator_group)

  # one table for each region
  tables=list()
  for(r in regs){

    # get parameter names
    tempbig <- lapply(ages, function(age) row.names(mods[[paste0(r, '_', age)]])) %>%
                  unlist %>%
                  unique %>%
                  data.table(parameter = .)

     tempbig = tempbig[!parameter%in%c('Theta1 for space','Theta2 for space','GroupRho for space'),]
     assign(r,tempbig)
      for(a in ages){
          tmp=data.table(cbind(parameter=row.names(mods[[paste0(r,'_',a)]]),
                                         mods[[paste0(r,'_',a)]]))[,c(1,2,4,6)]
          setnames(tmp, c("parameter",paste0("a",a,"_mean"),paste0("a",a,"_lb"),paste0("a",a,"_ub") ))

          # tranform spatial covariates to be more readbale
          tmp$parameter=as.character(tmp$parameter)
          tmp<-tmp[!parameter%in%c('Theta1 for space','Theta2 for space'),]
          tmp$parameter[tmp$parameter=="GroupRho for space"]="GPRandom Rho for time"
          tmp <- rbind(tmp,transform_spatial_effects(params.mat=mods[[paste0(r,'_',a)]],age=a))

          # round
          cols <- names(tmp)[2:4]
          tmp[,(cols) := round(.SD,4), .SDcols=cols]

          assign(r,merge(get(r),tmp,by='parameter',all=TRUE)) #,all.x=TRUE))

        #  get(r)$parameter[get(r)$parameter=="GroupRho for space"]
      }
      tables[[r]]<-get(r)

  }

  for(r in regs){
    write.csv(tables[[r]],sprintf('%s/model_results_table_%s.csv',sharedir,r))
  }
  return(tables)

}

# Load raking targets from GBD for vaccines

load_newest_gbd_vax  <- function(vaccine,
                                 gaul_list,
                                 return_cis = F,
                                 return_subnationals = F,
                                 gbd_year = 2017,
                                 years = c(2000:2016),
                                 gbd_date = "best") {

  str_match <- stringr::str_match

  # Convert between gaul code & IHME loc_id
  gaul_to_loc_id <- fread("<<<< FILEPATH REDACTED >>>>/gaul_to_loc_id.csv")
  names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"

  # Load table of GBD estimates to rake to
  gbd_file <- paste0('<<<< FILEPATH REDACTED >>>>/vacc_',vaccine,'.rds')
  gbd <- readRDS(gbd_file)
  gbd <- merge(gaul_to_loc_id, gbd, by="location_id")  

  if (return_subnationals) {
    # Use IHME shared functions to get location heirarchy
    source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
    source("<<<< FILEPATH REDACTED >>>>/get_location_metadata.R" )
    loc_hierarchy <- get_location_metadata(location_set_version = 280, gbd_round = 4)
    loc_hierarchy <- subset(loc_hierarchy, select = c("location_id", "parent_id", "location_name", 
                                                      "location_ascii_name", "location_name_short", 
                                                      "ihme_loc_id"))
    ihme_loc_ids <- gaul_to_loc_id[GAUL_CODE %in% gaul_list, ihme_lc_id]
    a_regexp <- paste(paste0(ihme_loc_ids, "_"), collapse = "|")
    subnat_location_ids <- loc_hierarchy[grepl(a_regexp, ihme_loc_id), location_id]

    gbdsubnat <- gbd[location_id %in% subnat_location_ids]
    gbdsubnat <- subset(gbdsubnat, year_id %in% years)
  }

  # subset to year & gaul lists
  gbd <- gbd %>%
            subset(., GAUL_CODE %in% gaul_list) %>%
            subset(., year_id %in% years)

  if (return_subnationals) {
    gbd <- rbind(gbd, gbdsubnat, use.names = T)
  }
  
  setnames(gbd, 
           c("GAUL_CODE", "year_id", "gpr_mean", "gpr_lower", "gpr_upper"),
           c("name", "year", "mean", "lower", "upper"))
  if (return_cis == F)  return_vars <- c("name", "year", "mean")
  if (return_cis == T)  return_vars <- c("name", "year", "mean", "lower", "upper")
  if (return_subnationals == T) {
    gbd[, parent_ihme_lc_id := str_match(ihme_lc_id,"(.*)_")[,2]]
    lc_id_table <- unique(subset(gbd, select = c("location_id", "loc_name", "ihme_lc_id")))
    setnames(lc_id_table, 
             c("location_id", "loc_name", "ihme_lc_id"), 
             c("parent_location_id", "parent_loc_name", "parent_ihme_lc_id"))
    gbd <- merge(gbd, lc_id_table, all.x=T, all.y=F)
    return_vars <- c(return_vars, c("loc_name", "location_id", "ihme_lc_id", 
                                    "parent_ihme_lc_id", "parent_location_id", "parent_loc_name"))
  }

  gbd <- subset(gbd, select = return_vars)

  return(gbd)
}

load_dropout_for_raking <- function(vaccine,
                                    gaul_list,
                                    years = c(2000:2016),
                                    abs_rel,
                                    gbd_date = "best") {

  vax1_gbd <- load_newest_gbd_vax(vaccine = paste0(vaccine, "1"),
                                  gaul_list = gaul_list,
                                  years = years,
                                  gbd_date = gbd_date)

  vax3_gbd <- load_newest_gbd_vax(vaccine = paste0(vaccine, "3"),
                                  gaul_list = gaul_list,
                                  years = years,
                                  gbd_date = gbd_date)

  setnames(vax1_gbd, "mean", "vax1")
  setnames(vax3_gbd, "mean", "vax3")
  vax1_3_dropout_gbd <- merge(vax1_gbd, vax3_gbd)
  if (abs_rel == "absolute") vax1_3_dropout_gbd[, mean := (vax1 - vax3)]
  if (abs_rel == "relative") vax1_3_dropout_gbd[, mean := ((vax1 - vax3)/vax1)]

  vax1_3_dropout_gbd <- subset(vax1_3_dropout_gbd, select = c("name", "year", "mean"))

   if (nrow(vax1_3_dropout_gbd[mean < 0,]) > 0) {
    n_fix <- nrow(vax1_3_dropout_gbd[mean < 0,])
    warning(paste0(n_fix, " country-year pairs have vax3 > vax1. Setting to 1e-4..."))
    message("Affected country-years: ")
    print(vax1_3_dropout_gbd[mean < 0,])
    vax1_3_dropout_gbd[mean < 0, mean := 1e-4]

  }

  return(vax1_3_dropout_gbd)

}

## check_input_data ################################################

#' Checks to see if input data has been saved into sharedir
#' If not, saves an input data object for the given region
#' This is useful if post-estimating indicators that are calculated
#' from other modeled indicators (e.g. if using continuation-ratio
#" approach for ordinal indicators)
#'
#' @param indicator Indicator
#' @param indicator_group Indicator group
#' @param age Age group
#' @param reg Region
#' @param run_date Run date
#' @param use_share Should the file look in share for master input_data csv?
#'
#' @return Does not return anything; saves input_data files by region
#'         into standard directory for this run
#'
#' @examples
#' 

check_input_data <- function(indicator,
                             indicator_group,
                             age,
                             reg,
                             holdout = 0,
                             run_date,
                             use_share = F,
                             ylist = year_list) {

  message(paste0("Checking input data for indicator: ", indicator,
                 " | region: ", reg, 
                 " | holdout: ", holdout))
  # make a pathaddin
  pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
  input_file <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", 
                       indicator, "/output/", run_date, "/input_data", pathaddin, ".csv")
  
  if (!file.exists(input_file)) {
    ## Load simple polygon template to model over
    gaul_list           <- get_gaul_codes(reg)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
    subset_shape        <- simple_polygon_list[[1]]
    simple_polygon      <- simple_polygon_list[[2]]

    message('Loading in full dataset using load_input_data() and saving...')

    df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                          simple      = simple_polygon,
                          agebin      = age,
                          removeyemen = TRUE,
                          pathaddin   = pathaddin,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = ylist)
  }
}

## distribute_config ################################################

#' Distribute config file to directories for a set of indicators with the same run date
#'
#' After running this, use load_config() with post_est_only == T to pull this config
#' Allows the config to be loaded once only from the master script, e.g. in case
#' you want to change it and launch another model shortly afterwards
#'
#' @param cfg config object (from `load_config()`)
#' @param indicators vector of indicator names to distribute config to
#' @param indicator_group indicator group
#' @param run_date shared run date for all indicators
#'
#' @return does not return any objects; writes CSV with config to each of the 
#'         indicator directories
#' @examples
#' 
# Distribute config file to directories for a set of indicators with the same run date.
# After running this, use `load_config()` with `post_est_only == T` to pull this config.
# Allows the config to be loaded once only from the master script, e.g. in case
#   you want to change it and launch another model shortly afterwards

distribute_config <- function(cfg = config,
                              indicators,
                              ig = indicator_group,
                              rd = run_date) {
  for (ind in indicators) {
    ind_dir <- paste0("<<<< FILEPATH REDACTED >>>>/",
                      ig, "/", ind,
                      "/output/", rd, "/")
    dir.create(ind_dir, recursive = T, showWarnings = F)
    write.csv(cfg, paste0(ind_dir,'/config.csv'), row.names = FALSE)
  }
}

## save_stratum_ho() ################################################

#' Save a stratum_ho (input data stratified into holdouts) 
#' object to a standardized location
#'
#' @param indic indicator for which holdouts have been generated
#' @param ig indicator group for `indic`
#' @param rd run_date for `indic`
#' @param stratum_obj the `stratum_ho` object to save
#' @return Saves an .rds file to the indicator's directory
#' @examples
#' 

save_stratum_ho <- function(indic = indicator,
                            ig = indicator_group,
                            rd = run_date,
                            stratum_obj = stratum_ho) {

  out_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", indic, "/output/", rd, "/")
  saveRDS(stratum_obj, file = paste0(out_dir, "stratum.rds"))
}

# Recreate stratum_ho object given a data frame, row_ids, and indicator for which holdouts have been generated
# This will assign row_ids the same holdout as in the load_from_indic
# Note that df (for your current vaccine) must be the same or a subset compared to the load_from_indic

## recreate_holdouts() ################################################

#' For indicator Y, recreate holdouts generated for indicator X using a 
#' row ID column. This is useful, for instance, in ordinal regression
#' settings if you want to generate holdouts once for the entire data
#' set, and apply those holdouts to each derivative data object (eg
#' the conditional data objects).  
#' 
#' In order for this to work, the holdouts must be generated on the 
#' most comprehensive data set (the one with all of the row_ids).  Those
#' holdouts can then be applied to other data sets with either the same
#' row_ids or a subset of the row_ids from the original data set
#'
#' @param data new data that you'd like to apply the holdouts to
#' @param row_id_col column in both `data` and the original data set
#'                   used to generate the master set of holdouts that
#'                   links the two data sets. Must be unique within
#'                   each data set
#' @param load_from_indic which indicator (within the same run_date and
#'                        indicator group) would you like to load the master
#'                        set of holdouts from?
#' @param rd run_date
#' @param ig indicator group
#' @return a `stratum_ho` style object containing the data from `data` but 
#'         divided up by the holdouts from `load_from_indic`
#' @examples
#' 

recreate_holdouts <- function(data = df,
                              row_id_col = "row_id",
                              load_from_indic,
                              rd = run_date,
                              ig = indicator_group) {
  data <- copy(data)
  n_before <- nrow(data)
  s_ho <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", load_from_indic, 
                         "/output/", rd, "/", "stratum.rds"))

  stratum_df <- rbindlist(s_ho)

  # Peel off just the columns of interest
  stratum_df <- subset(stratum_df, select = c(row_id_col, "t_fold", "fold", "region"))

  # Replace any region columns
  if ("region" %in% names(data)) data[, region := NULL]

  # Merge with df; drop if no holdout assigned (buffer points)
  data <- merge(data, stratum_df, all.x = F, all.y = F, by = row_id_col)
  n_after <- nrow(data)
  
  # get regions
  regions <- unique(data$region)

  out_list <- lapply(regions, function(r) {
    return(data[region == r,])
  })

  names(out_list) <- paste0("region__", regions)

  return(out_list)

}

# Check for .tif files for a given set of indicators, run_dates, and summary stats

## check_for_tifs ################################################

#' Check for .tif files for a given set of indicators, summary stats, and (optional) regions
#' within a particular run date.
#'
#' @param indicators vector of indicators
#' @param rd run date
#' @param ig indicator group
#' @param regs vector of regions. If NULL will search for the already-combined files
#' @param rake either T/F (for all to take that value) or vector of T/F indicating whether indicator was raked or not
#' 
#' @return either gives a success message and returns NULL or gives a data table of missing tifs
#' @examples
#' 
#' check_for_tifs(indicators = c("dpt3_cov", "dpt3_cov", dpt1_cov"),
#'                rd = run_date,
#'                regs = Regions,
#'                rake = c(T,F,F),
#'                summstats = c("mean", "upper", "lower"))

check_for_tifs <- function(indicators,
                           rd   = run_date,
                           ig   = indicator_group,
                           regs = NULL,
                           rake = F,
                           summstats) {

  library(data.table)
  str_match <- stringr::str_match

  # Set up raked params
  if (length(rake) == 1) {
    rake <- rep(rake, length(indicators)) 
  } else if ((length(rake) > 1) & (length(rake) != length(indicators))) {
    stop("length(rake) must equal length(indicators)")
  }

  pasted_indicators <- paste(indicators, rake, sep = "_-_")

  if (!is.null(regs)) {
    tif_table <- data.table(expand.grid(regs, summstats, pasted_indicators), stringsAsFactors = F)
    names(tif_table) <- c("region", "summstat", "pasted_indicator")
  } else if (is.null(regs)) {
    tif_table <- data.table(expand.grid(summstats, pasted_indicators), stringsAsFactors = F)
    names(tif_table) <- c("summstat", "pasted_indicator")
    
    # Fix cirange for naming convention for merged rasters
    tif_table[summstat == "cirange", summstat := "range"]
  }

  tif_table[, indicator := str_match(pasted_indicator, "(.*)_-_(.*)")[,2]]
  tif_table[, raked := str_match(pasted_indicator, "(.*)_-_(.*)")[,3]]
  tif_table[, pasted_indicator := NULL]
  tif_table[raked == T, raked_addin := "_raked"]
  tif_table[raked == F, raked_addin := ""]

  # Make filenames
  if (!is.null(regs)) {
    tif_table[, filename := paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", indicator, "/output/", run_date, "/",
                                   indicator, "_", region, raked_addin, "_", summstat, "_raster.tif")]
  } else if (is.null(regs)) {
    tif_table[, filename := paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", indicator, "/output/", run_date, "/",
                               indicator, "_", summstat, raked_addin, "_raster.tif")]
  }

  tif_table[,raked_addin := NULL]

  # Check if files exist
  tif_table[, exists := file.exists(filename)]
  missing_files <- subset(tif_table, exists == F)

  # Return
  if (nrow(missing_files) == 0) {
    message("All files are present.")
    return(NULL)
  } else if (nrow(missing_files) > 0) {
    message(paste0(nrow(missing_files), " out of ", nrow(tif_table), " files are missing."))
    message("Returning data table of missing files.")
    return(missing_files)
  }
  
}

# Save TIFs to the standard directory / format for mapping

copy_tif_to_map_input_dir <- function(ind,
                                      ig = indicator_group, 
                                      measure, 
                                      rd = run_date, 
                                      raked, 
                                      yl = year_list) {

  # Fix measures
  if (measure == "cirange") measure <- "range"

  # Fix raked if boolean
  if (raked == T) raked <- "raked"
  if (raked == F) raked <- "unraked"

  # Define input dirs / files
  in_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", ind, "/output/", rd, "/")
  in_file <- paste0(in_dir, ind, "_", measure, "_", ifelse(raked == "raked", "raked_", ""), "raster.tif")
  
  # Define output dirs / files
  out_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/", rd, "/inputs/")
  dir.create(out_dir, recursive = T, showWarnings =F)

  out_filename <- paste0(ind, "_", measure, "_", raked, "_", min(yl), "_", max(yl), ".tif")
  out_file <- paste0(out_dir, out_filename)

  success <- file.copy(in_file, out_file, overwrite = T)

  if (success) message(paste0(out_filename, " copied successfully!"))
  if (!success) warning(paste0("Failed to copy ", out_filename, "!")) 

  return(invisible(success))

}

# Save CSVs to the standard directory / format for mapping

copy_admins_to_map_input_dir <- function(ind,
                                         ig = indicator_group, 
                                         measure, 
                                         rd = run_date, 
                                         raked, 
                                         yl = year_list, 
                                         ad_level) {

  str_match <- stringr::str_match

  # Fix raked if boolean
  if (raked == T) raked <- "raked"
  if (raked == F) raked <- "unraked"

  # Define basic dirs / files
  in_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", rd, "/")
  out_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/", rd, "/inputs/")
  dir.create(out_dir, recursive = T, showWarnings =F)

  if (measure %in% c("mean", "lower", "upper", "cirange", "cfb")) {
    in_dir <- paste0(in_dir, "pred_derivatives/admin_summaries/")
    in_file <- paste0(in_dir, ind, "_admin_", ad_level, "_", raked, "_summary.csv")  
    in_df <- as.data.table(read.csv(in_file, stringsAsFactors = F))
    setnames(in_df, measure, "value")

    if (measure == "cirange") measure <- "range" #naming convention
    out_df <- subset(in_df, select = c(paste0("ADM", ad_level, "_CODE"), "year", "value")) 
    out_filename <- paste0(ind, "_", measure, "_", raked, "_ad", ad_level, ".csv")
    message(paste0("Writing ", out_filename, "..."))
    write.csv(out_df, file = paste0(out_dir, out_filename))
  
  } else if  (grepl("psup", measure)) {

    pct <- str_match(measure, "psup(.*)")[,2]
    pct_div100 <- as.numeric(pct)/100

    in_dir <- paste0(in_dir, "pred_derivatives/admin_summaries/")
    in_file <- paste0(in_dir, ind, "_admin_", ad_level, "_", raked, "_p_", pct_div100, "_or_better_summary.csv")  
    in_df <- as.data.table(read.csv(in_file, stringsAsFactors = F))
    setnames(in_df, "p_above", "value")

    out_df <- subset(in_df, select = c(paste0("ADM", ad_level, "_CODE"), "year", "value"))
    out_filename <- paste0(ind, "_psup", pct, "_", raked, "_ad", ad_level, ".csv")

    message(paste0("Writing ", out_filename, "..."))
    write.csv(out_df, file = paste0(out_dir, out_filename))

  } else if (grepl("diff", measure)) {
  
    in_dir <- paste0(in_dir, "pred_derivatives/admin_summaries/")
    in_file <- paste0(in_dir, ind, "_admin_", ad_level, "_", raked, "_", measure, ".csv")

    in_df <- as.data.table(read.csv(in_file, stringsAsFactors = F))
    setnames(in_df, "mean", "value")

    out_df <- subset(in_df, select = c(paste0("ADM", ad_level, "_CODE"), "year", "value"))
    out_filename <- paste0(ind, "_diff_2000-2016_", raked, "_ad", ad_level, ".csv")

    message(paste0("Writing ", out_filename, "..."))
    write.csv(out_df, file = paste0(out_dir, out_filename))

  } else {

    warning(paste0("Measure ", measure, " not found for indicator ", ind, " | ", raked, " | admin level ", ad_level))

  }
}


## plot_mbg_vs_gbd ################################################
  
#' Plot MBG vs GBD estimates along with data
#'
#' @param ind indicator
#' @param ig indicator_group
#' @param rd run_date
#' @param gbd gbd_indicator
#' @param new_vax use new vaccines from load_newest_gbd_vax? 
#'          (only T for now; need to build in basic GBD defaults)   
#' @param yl year_list
#' @param vax_title title of vaccine for display
#' @param raked use raked estimates? (T/F)
#' 
#' @return creates png images in a /comparisons/ folder in output directory
#' 
#' @examples
#' plot_mbg_vs_gbd(ind = "dpt3_cov",
#'                 ig = "vaccine",
#'                 rd = run_date,
#'                 gbd_ind = "dpt3",
#'                 yl = year_list,
#'                 vax_title = "DPT3")

plot_mbg_vs_gbd <- function(ind = indicator,
                            ig = indicator_group,
                            rd = run_date,
                            gbd_ind,
                            new_vax = T,
                            yl = year_list,
                            vax_title,
                            raked = F,
                            gbd_date = "best") {
  # Compare admin0s

  # Set up files and directories
  rd_dir <-paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", rd, "/")
  adm_dir <- paste0(rd_dir, "pred_derivatives/admin_summaries/")
  adm0_file <- paste0(adm_dir, ind, "_admin_0_", ifelse(raked, "raked", "unraked"), "_summary.csv")

  out_dir <- paste0(rd_dir, "comparisons/")
  dir.create(out_dir, showWarnings = F)

  adm0_df <- fread(adm0_file)

  gbd <- load_newest_gbd_vax(vaccine = gbd_ind,
                         gaul_list = unique(adm0_df$ADM0_CODE),
                         years = yl,
                         return_cis = T,
                         gbd_date = gbd_date)

  setnames(gbd, c("name", "mean", "lower", "upper"), c("ADM0_CODE", "gbd", "gbd_lower", "gbd_upper"))
  adm0_df <- merge(adm0_df, gbd, by = c("ADM0_CODE", "year"), all.x = T)

  missing_countries <- unique(adm0_df[is.na(gbd), ADM0_NAME])
  if (length(missing_countries) > 0) {
    warning(paste0("The following countries have no GBD value and will be dropped: \n",
            paste(missing_countries, collapse = ", ")))
  }

  adm0_df <- subset(adm0_df, !(ADM0_NAME %in% missing_countries))

  # Truncate long country names
  adm0_df[ADM0_NAME == "Democratic Republic of the Congo", ADM0_NAME := "DRC"]
  adm0_df[ADM0_NAME == "Sao Tome and Principe", ADM0_NAME := "STP"]
  adm0_df[ADM0_NAME == "United Republic of Tanzania", ADM0_NAME := "Tanzania"]
  adm0_df[ADM0_NAME == "Central African Republic", ADM0_NAME := "CAR"]
  adm0_df[ADM0_NAME == "Equatorial Guinea", ADM0_NAME := "Eq. Guinea"]

  # Now pull in the input data
  output_draws <- fread(paste0(rd_dir, "output_draws_data.csv"))

  # Drop draw columns
  output_draws <- subset(output_draws, select = names(output_draws)[!grepl("draw[0-9]+", names(output_draws))])
  output_draws[, V1 := NULL]
  ad0_lookup <- data.table(ad0 = unique(output_draws$ad0))
  ad0_lookup[, ADM0_CODE := get_gaul_codes(ad0)]
  output_draws <- merge(output_draws, ad0_lookup)
  output_draws[, outcome := get(ind) / N]

  group_by <- c("svy_id", "source", "country", "ADM0_CODE", "year")
  data_est <- output_draws[, .(mean = weighted.mean(outcome, w = weighted_n),
                               data_tot_n = sum(weighted_n)), 
                 by = group_by]
  data_est[, estimate := "Data"]

  # Make plots
  # Max = 16 per page
  n_pages <- ceiling(length(unique(adm0_df$ADM0_NAME))/16)
  adm0_df <- adm0_df[order(ADM0_NAME, year)]
  assign_table <- data.table(ADM0_NAME = unique(adm0_df$ADM0_NAME),
                             page = rep(1:n_pages, each=16, length.out = length(unique(adm0_df$ADM0_NAME))))
  adm0_df <- merge(adm0_df, assign_table)                           

  for (pg in 1:n_pages) {

    # Prepare time series data
    plot_df <- adm0_df[page == pg]
    plot_df_gbd <- subset(plot_df, select = !(names(plot_df) %in% c("mean", "lower", "upper")))
    plot_df_gbd[, estimate := "GBD Estimate"]
    setnames(plot_df_gbd, c("gbd", "gbd_lower", "gbd_upper"), c("mean", "lower", "upper"))
    plot_df_mbg <- subset(plot_df, select = !(names(plot_df) %in% c("gbd", "gbd_lower", "gbd_upper")))
    plot_df_mbg[, estimate := paste0(ifelse(raked, "Raked", "Unraked"), " MBG Estimate")]
    plot_df <- rbind(plot_df_mbg, plot_df_gbd)
    
    # Add data estimates
    data_df <- subset(data_est, ADM0_CODE %in% unique(plot_df$ADM0_CODE))
    setnames(data_df, "mean", "data_mean")
    data_df <- merge(data_df, unique(subset(plot_df, select = c("ADM0_CODE", "ADM0_NAME"))))

    all_df <- rbind(data_df, plot_df, fill=T, use.names=T)

    out_file = paste0(out_dir, ind, "_compare_", ifelse(raked, "raked", "unraked"),"_mbg_to_gbd_ad0_pg", pg, ".png")

    png(filename=out_file, 
        type="cairo",
        units="in", 
        width=14, 
        height=8, 
        pointsize=12, 
        res=300)

    gg <- ggplot(data = all_df, 
                 aes(x = year, color = estimate, group = estimate)) +
    geom_point(aes(y = mean), alpha = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0, alpha = 0.2) +
    geom_line(aes(y = mean), alpha = 0.2) +
    geom_point(aes(y = data_mean, size = data_tot_n), alpha = 0.5, shape = 4) +
    facet_wrap(~ADM0_NAME) +
    theme_bw() +
    labs(x = "Year", 
         y = paste0(vax_title, " Coverage"),
         color = "Estimate",
         size = "Data: Weighted N") +
    scale_color_manual(values = c("black", "red", "blue")) +
    ylim(c(0,1))

    plot(gg)

    dev.off()
  }
}

## plot_summary_metrics ################################################

#' Plot summary metrics by year and by region
#'
#' Function to plot summary metrics by year and by region, depending on
#' what options are selected.  Produces time trends.  Useful as diagnostic
#' plots, for instance, or for publications
#'
#' @param indic indicator
#' @param ig indicator group
#' @param rd run date
#' @param use_oos use OOS predictions? T, F, or c(T,F)
#' @param by_regions plot by regions? Requires that summary
#'                   metric tables have been generated by region
#' @param ad_levels country, ad1, or ad2 - or a vector of these
#' @return writes a series of png files to [run date dir]/summary_metrics/
#' @examples
#' plot_summary_metrics("dpt3_cov", "vaccine", "2018_02_10_30_20")

plot_summary_metrics <- function(indic,
                                 ig,
                                 rd,
                                 use_oos = c(T,F),
                                 by_regions = c(T,F),
                                 ad_levels = c("country", "ad1", "ad2")) {
  
  
  # Define a color scheme for the regions
  get_color_scheme <- function(theme){
    
    # Set up categorical colors
    carto_discrete <- c("#7F3C8D","#11A579","#3969AC","#F2B701",
                        "#E73F74","#80BA5A","#E68310","#008695",
                        "#CF1C90","#f97b72","#4b4b8f","#A5AA99")
    return(get(theme))
  }
  
  in_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", indic, "/output/", rd, "/summary_metrics/")
  
  for (level in ad_levels) {
    
  message(paste0("Working on the ", level, " level..."))
  # Read in csvs
  df_reg <- fread(paste0(in_dir, level, "_metrics_by_region.csv"))
  df_reg$region <- toupper(df_reg$region)
  df_all <- fread(paste0(in_dir, level, "_metrics.csv"))
  
  # Combine data tables
  df_all[, region := "All"]
  df_all <- rbind(df_reg, df_all, use.names = T)
  rename_table <- data.table(rbind(c("Year", "year", "Year"),
                                   c("OOS", "oos", "OOS"),
                                   c("Mean Err.", "me", "Mean Error"),
                                   c("Mean Abs Err.", "mae", "Mean Absolute Error"),
                                   c("Mean Pred.", "mean_pred", "Mean Prediction"),
                                   c("Mean Obs.", "mean_obs", "Mean Observed value"),
                                   c("RMSE", "rmse", "RMSE"),
                                   c("Median SS", "median_ss", "Median Sample Size"),
                                   c("Corr.", "corr", "Correlation"), 
                                   c("95% Cov.", "cov_95", "95% Coverage")))
  
  names(rename_table) <- c("old", "new", "title")
  setnames(df_all, rename_table$old, rename_table$new)
  
  message(paste0("Writing images to ", in_dir))
  # Loop over variables & make plots
    for (is_oos in use_oos) {
      for (by_region in by_regions) {
        for (plot_obj in c("me", "mae", "rmse", "corr", "cov_95")) {
          
          if (level == "ad2") ad_caption = "admin 2"
          if (level == "ad1") ad_caption = "admin 1"
          if (level == "country") ad_caption = "country"
            
          plot_title <- rename_table[new == plot_obj, title]
          
          gg <- ggplot()
          
          if (by_region == T) {
            gg <- gg + geom_line(data = df_all[region != "All" & oos == is_oos,],
                                 aes(x = year, y = get(plot_obj), group = region, color = region), alpha = 0.35) +
                       scale_color_manual(values = c("#000000", get_color_scheme("carto_discrete"))) +
                       labs(color = "Region") +
                       geom_line(data = df_all[region == "All" & oos == is_oos],
                                 aes(x = year, y = get(plot_obj), group = region, color = region))
          } else {
            gg <- gg + geom_line(data = df_all[region == "All" & oos == is_oos],
                                 aes(x = year, y = get(plot_obj)), color = "black")
          }
            
          gg <- gg +
            labs(x = "Year", y = plot_title, title = plot_title, 
                subtitle = paste0(ifelse(is_oos, "Out-of-sample", "In-sample"), " | Aggregated to the ", ad_caption, " level")) +
            theme_bw()
          
          # Get range
          range_obs <- df_all[oos == is_oos, get(plot_obj)]
          
          if (plot_obj == "me") gg <- gg + ylim(-max(abs(range_obs)), max(abs(range_obs))) + geom_abline(slope = 0, intercept = 0, color = "red", linetype = 2)
          if (plot_obj == "mae" | plot_obj == "rmse") gg <- gg + ylim(0, max(abs(range_obs))) + geom_abline(slope = 0, intercept = 0, color = "red", linetype = 2)
          if (plot_obj == "cov_95") gg <- gg + ylim(0,1) + geom_abline(slope = 0, intercept = 0.95, color = "red", linetype = 2)
          if (plot_obj == "corr") gg <- gg + ylim(0,1) + geom_abline(slope = 0, intercept = 1, color = "red", linetype = 2)
          
          message(paste0(" Writing file for ", plot_title, " | by_region = ", by_region, " | oos = ", is_oos))
          png_filename <- paste0(indic, "_summary_metric_plot_", level, "_", 
                                 ifelse(is_oos, "OOS" ,"IS"), "_", ifelse(by_region, "by_region", "all"), 
                                 "_", plot_obj, ".png")
        
          png(file = paste0(in_dir, png_filename), 
              width = 8, 
              height = 4, 
              units = "in",
              res = 300)
          plot(gg)
          dev.off()
          
          message(paste0("  File written to ", png_filename))
        }
      }
    }
  }
}
