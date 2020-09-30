## #############################################################################
## TEST THE NEW VERSION OF GET_GAUL_CODES AGAINST THE PREVIOUS VERSION
##
## Purpose: Test that the new get_gaul_codes() function is working the same as
##   the previous version for all previously defined regions.
##
## #############################################################################

library(data.table)


## HELPER FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load the GAUL lookup table
load_gaul_lookup_table <- function(){
  # Define consistent path to the GAUL lookup table
  lookup_table_filepath <-"<<<< FILEPATH REDACTED >>>>"
  # If data.table is loaded, read in as a data.table
  if ('package:data.table' %in% search()){
    lookup_table <- fread(lookup_table_filepath)
  } else {
     # Otherwise, read in as a data.frame
     lookup_table <- read.csv(lookup_table_filepath)
  }
  # Set ISO codes as lowercase for easy lookup
  lookup_table$iso3 <- tolower(lookup_table$iso3)
  return (lookup_table)
}


## Check the ISO code related to each GAUL code
gauls_to_isos <- function(gauls, lookup_table = load_gaul_lookup_table()){
  return(sapply(gauls,
                function(x) {paste0(x,'(',lookup_table[GAUL_CODE==x,iso3],')')}))
}



## Test that two lists ('old' and 'new') contain the same values
## Prints a short line if the lists contain the same values; prints a longer
##  error message if the lists contain different values
test_for_same_values <- function(old, new, region_name){
  missing_from_old <- new[ !(new %in% old) ]
  missing_from_new <- old[ !(old %in% new) ]

  if ( (length(missing_from_old)==0) & (length(missing_from_new)==0) ){
    message(sprintf("Test passed for %s.",region_name))
  } else {
    message(sprintf("**** DIFFERENT VALUES FOR REGION %s ****",region_name))
    message(sprintf("  Present in new, not in old: %s",
                    paste(gauls_to_isos(missing_from_old),collapse=', ')))
    message(sprintf("  Present in old, not in new: %s\n",
                    paste(gauls_to_isos(missing_from_new), collapse=', ')))
  }
}



## Run through a set of character vectors defined in both the old and new
##  versions of get_gaul_codes to make sure that all returned codes are the same
##  in both versions
compare_get_gaul_codes <- function(old_repo_path, new_repo_path, test_regions){
  # Load the old version of get_gaul_codes and change the function name
  source(paste0(old_repo_path,'mbg_central/prep_functions.R'))
  get_gaul_func_old <- get_gaul_codes
  # Load the new version of the repo
  source(paste0(new_repo_path,'mbg_central/prep_functions.R'))
  # Iterate through each test region, comparing output codes for each
  for (reg in test_regions){
    # Old function results
    old_results <- get_gaul_func_old(gaul=reg)
    # New function results
    new_results <- get_gaul_codes(gaul=reg, core_repo=new_repo_path)
    # Compare
    test_for_same_values(old_results, new_results, region_name=reg)
  }
}


## MAIN SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set repository paths
old_repo_path <- "<<<< FILEPATH REDACTED >>>>"
new_repo_path <- "<<<< FILEPATH REDACTED >>>>"

# Run test 1 (all user-defined modeling regions)
modeler_custom_regions <- c(
  'ind_all_territories','cssa_diarrhea','sssa_diarrhea','essa_diarrhea',
  'cssa_diarrhea2','essa_diarrhea2','name_diarrhea2','sssa_diarrhea2',
  'wssa_diarrhea2','essa_edu','cssa_edu','name_edu','sssa_edu','wssa_edu',
  'essa_sdn','cessa','cwssa','namelite','cessa2','sssa2','cssa_cam',
  'wssa_nocam','name_hi','essa_hi','essa_lo','cssa_hi','cssa_lo','wssa_hi',
  'wssa_lo','sssa_hi','sssa_lo','essa_hilo'
  )

# All of these should pass except for 'cessa2', because American Samoa is
#  definitely not in Africa
compare_get_gaul_codes(old_repo_path, new_repo_path, modeler_custom_regions)


# Run test 2 (large regions and continents)
large_regions <- c('africa','middle_east','eastern_europe','latin_america',
                   'south_asia','central_america','south_america','se_asia',
                   'stage1','stage2')
# Many of these are different, but the differences seem logical: either region
#  definitions have been updated recently, or the locations should have never
#  been in the given region in the first place. All stage 3 countries are also
#  now excluded.
compare_get_gaul_codes(old_repo_path, new_repo_path, large_regions)

# Run test 3 (selected countries in Africa)
afr_isos <- c('nam','egy','zaf','dji','eth','ken','gha')
compare_get_gaul_codes(old_repo_path, new_repo_path, afr_isos)
