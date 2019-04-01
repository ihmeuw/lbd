# HEADER ------------------------------------------------------------------
# 00_master.R
# Purpose: Set up options & launch child scripts that control the 
#          data preparation process
#**************************************************************************

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<<<< FILEPATH REDACTED >>>>'
indic_repo         <- '<<<< FILEPATH REDACTED >>>>'

## sort some directory stuff
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

# Specific function loads (avoid namespace conflicts)
str_match <- stringr::str_match

### II: DEFINE VARIABLES, ETC. #####################################################

# Variables ------------------------------------------------------------------------
  run_date     <- NULL # Specify here as "YYYY_MM_DD" if running on specific date
                       # Otherwise will use today if re-reading data
                       # Or the most recent existing if not
  cores 			 <- 30                    

# Vaccines & vaccine titles --------------------------------------------------------

  vaccine_list <- add_vaccine(prefix = "dpt", 
                              title = "DPT",
                              doses = 3,
                              age_min = 12,
                              age_max = 59)

  vaccines <- names(vaccine_list)

  graph_vaccines <- c("dpt3_cov")
  final_model_vaccines <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout", "dpt1_3_rel_dropout")
  graph_vaccine_titles <- c("DPT3")
  log_dir <- "<<<< FILEPATH REDACTED >>>>"
  dir.create(log_dir, showWarnings = F, recursive = T)
  
# Toggles --------------------------------------------------------------------------
  
  drop_pre_2000   <- T     # Drop surveys and individuals before 2000

  read_data       <- T     # Read in data and combine
  merge_data      <- T     # Merge all the files together to vaccine_clean_list_file (below)
  make_graphs     <- T     # Make data coverage graphs
  do_pointpoly    <- T     # Run point/poly processing
  make_csvs       <- T     # Generate CSVs

# Directories & files --------------------------------------------------------------
  
  # directory for data prep code
  script_dir <- paste0(repo, "vaccine/data_prep/")

  # filepath holding extracted survey data in .csv format
  in_dir <- "<<<< FILEPATH REDACTED >>>>"

  # output directory
  out_dir <- "<<<< FILEPATH REDACTED >>>>"

# Prep a run date-------------------------------------------------------------------

 if (is.null(run_date) & read_data == F) {
  # make list of dates already  generated
    dates <- list.files(in_dir)
    dates <- sort(dates[grep(".*_.*_.*", dates)], decreasing = T)
  # pull most recent date folder & use that
    run_date <- dates[1]
 } else if (is.null(run_date) & read_data == T) {
 	run_date <- Sys.Date() #date-stamp 
	run_date <- gsub("-", "_", run_date)
 }

 in_dir <- paste0(in_dir_stem, run_date, "/")
 dir.create(in_dir, showWarnings = F)

 out_dir <- paste0(out_dir_stem, run_date, "/")
 dir.create(out_dir, showWarnings = F)  

# File for temporary, cleaned data set of processed vaccines - in same data as input data
# Will have the following added to it:  "..._[vaccine_prefix].rds"
 vaccine_cleaned_file_prefix <- paste0(in_dir, "combined_vaccine_clean_")

# File for temporary, *resampled* data set of processed vaccines - in same data as input data
# Will have the following added to it:  "..._[vaccine_prefix].rds"
 vaccine_resampled_file_prefix <- paste0(in_dir, "combined_vaccine_resampled_")

### III: RUN BITS PER TOGGLES ######################################################

if (read_data == T) {
  # This will read in all data and geoposition it
  # Output: combined CSV files, one per survey series
  print_header()
  message("\nReading data...\n")
  source(paste0(script_dir, "01_highspeed_cluster_post_extraction_v2.R"))
}

if (merge_data == T) {
 # This will combine the survey-specific CSV files, then merge and apply case definitions
 # Output: an RData file (one line per individual)
 print_header()
 message("\nMerging vaccine data...\n")
 source(paste0(script_dir, "02_merge_vaccine.R"))
}

if (make_graphs == T) {
 # Creates graphs (maps and scatter plots) from the file made in merge_data
 # Outputs: vaccine-specific graphs and scatter plots
 print_header()
 message("\nGraphing coverage...\n")
 source(paste0(script_dir, "03_graph_coverage.R"))
}

if (do_pointpoly == T) {
 # Creates graphs from the file made in merge_data
 # Outputs: pointpoly-ed individual csvs in out_dir
 print_header()
 message("\nRunning pointpoly...\n")
 source(paste0(script_dir, "04_pointpoly.R"))
}

if (make_csvs == T) {
 # Creates graphs from the file made in merge_data
 # Outputs: pointpoly-ed individual csvs in out_dir
 print_header()
 message("\nMaking CSVs...\n")
 source(paste0(script_dir, "05_generate_vaccine_csvs.R"))
}
