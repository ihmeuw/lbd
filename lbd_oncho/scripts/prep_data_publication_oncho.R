### Prep Oncho MBG data for publication
### Purpose: Collapse data sources (by NID) and generate metadata. Also, collapse admins in aggregated estimate CSVs for making codebooks
##########################################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Script Prep #############################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Define scipt options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
indicator_group = "oncho"
indicator_poly = "had_oncho_poly"
indicator_w_resamp = "had_oncho_w_resamp"
data_tag = "_2020_08_03"
out_dir = <<<< FILEPATH REDACTED >>>>
pathadd = "_Y2020M08D18"
use_consideration_filename = "use_considerations_081820" #filename for use_consideration file (used to check for private sources to remove from output). See last code section.

### Set repo locations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get current user name
user <- Sys.info()[["user"]] ## Get current user name

## Set repo location and indicator group
core_repo <- <<<< FILEPATH REDACTED >>>>
indic_repo <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

# NOTE: package load errors likely due to R version problems
#       check .libPaths() and R.version. Change sing_imagename in config to singularity image with correct internals.
mbg_setup(package_list = package_list, repos = core_repo)

source(<<<< FILEPATH REDACTED >>>>)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Get Data ########################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
filename <- <<<< FILEPATH REDACTED >>>>
mbg_data <- fread(<<<< FILEPATH REDACTED >>>>)

################################################### COllAPSE DATA ###########################################################


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Get location names #############################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
locs <- get_location_metadata(gbd_round_id = 7, location_set_id = 35) %>% as.data.table %>% .[, .(ihme_loc_id, location_name)]
mbg_data <- merge(mbg_data, locs, by.x = "country", by.y = "ihme_loc_id", all.x = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Collapse #############################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rename var for readability in output
mbg_data[diagnostic == "ss", diagnostic := "skin snip"]
mbg_data[diagnostic == "nod", diagnostic := "nodule"]

# create data frame of data input records for output
end_frame <- mbg_data[is.na(nid) == F & nid != "", .(nid, location_name, diagnostic, point)][, lapply(.SD, FUN = function(x) paste0(unique(x), collapse = ", ")), by = nid]
end_frame[grep(",", end_frame$point), 'Geographic Detail' := "Intervention Unit, GPS"]
end_frame[point == 0, 'Geographic Detail' := "Intervention Unit"]
end_frame[point == 1 | is.na(point) == T, 'Geographic Detail' := "GPS"]
end_frame$point <- NULL

# summarize age_start by nid
end_frame$age_start <- lapply(end_frame$nid, FUN = function(x){
  begin = suppressWarnings(min(as.integer(mbg_data[nid == x, age_start]), na.rm = T))
  
  if(begin == 0|begin == Inf|begin == -Inf) begin <- 1
  return(begin)
}) %>% unlist

# summarize age_end by nid
end_frame$age_end <- lapply(end_frame$nid, FUN = function(x){
  end = suppressWarnings(max(as.integer(mbg_data[nid == x, age_end]), na.rm = T))
  if(end == 0|end == Inf|end == -Inf) end <- 99
  
  return(end)
}) %>% unlist

# summarize survey data years by nid
end_frame$years <- lapply(end_frame$nid, FUN = function(x){
  begin = suppressWarnings(min(as.integer(mbg_data[nid == x, year]), na.rm = T))
  end = suppressWarnings(max(as.integer(mbg_data[nid == x, year]), na.rm = T))
  
  if(begin == end){
    return(begin)
  }else{
    return(paste0(begin, "-", end))
  }
}) %>% unlist

# sum up sample size
end_frame$sample_size <- lapply(end_frame$nid, FUN = function(x){
  return(sum(mbg_data[nid == x, N]))
}) %>% unlist

# summarize # of points/polygons
end_frame$'Number of Geo-positioned Clusters' <- lapply(end_frame$nid, FUN = function(x){
  return(nrow(mbg_data[nid == x & point == 1,]))
}) %>% unlist

end_frame$'Number of Polygons (Areal)' <- lapply(end_frame$nid, FUN = function(x){
  return(nrow(mbg_data[nid == x & point == 0,]))
}) %>% unlist

# create fields that are consistent across all records
end_frame[, sex := "Both sexes"]
end_frame[, representativeness := "Subnationally representative only"]

#check for placeholder NIDs
if(length(grep("temp", end_frame$nid)) > 0 | length(grep(" ", end_frame$nid)) > 0) message(paste0(length(grep("temp", end_frame$nid) + length(grep(" ", end_frame$nid))), " NIDs are placeholders. Check on these."))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Make sure restricted access source remain private ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use the online tool for checking use consideration
all_nids_string <- paste0(end_frame[!(nid %in% grep("temp", end_frame$nid, value = T)) & !(nid %in% grep(" ", end_frame$nid, value = T)), nid], collapse = ",")
url <- <<<< FILEPATH REDACTED >>>>
message("Paste the following url to query use consideration by NID: ", url)

use_considerations <- fread(<<<< FILEPATH REDACTED >>>>)
restricted_nids <- use_considerations[Confidentiality == "Private", Nid]

end_frame_copy <- copy(end_frame)
end_frame <- end_frame[!(nid %in% restricted_nids),]

write.csv(end_frame, <<<< FILEPATH REDACTED >>>>, row.names = F)

################################################### COLLAPSE ADMINS FOR CODEBOOKS ###########################################################
if(TRUE){
  run_dates <- <<<< FILEPATH REDACTED >>>> # run_date for the MBG model; used for creating Codebooks
  
  aggregated_data <- lapply(run_dates, FUN = function(r){
    d_paths <- list.files(<<<< FILEPATH REDACTED >>>>, full.names = T)
    d <- lapply(d_paths, fread)
    return(d)
  })
  
  adm0s <- lapply(aggregated_data, FUN = function(d){
    d[[1]][, .(ADM0_CODE, ADM0_NAME)] %>% unique %>% return
  }) %>% rbindlist
  
  adm1s <- lapply(aggregated_data, FUN = function(d){
    d[[2]][, .(ADM1_CODE, ADM1_NAME)] %>% unique %>% return
  }) %>% rbindlist
  
  adm2s <- lapply(aggregated_data, FUN = function(d){
    d[[3]][, .(ADM2_CODE, ADM2_NAME)] %>% unique %>% return
  }) %>% rbindlist
  
  write.csv(adm0s, paste0(out_dir, "adm0s", pathadd, ".csv"), row.names = F)
  write.csv(adm1s, paste0(out_dir, "adm1s", pathadd, ".csv"), row.names = F)
  write.csv(adm2s[!(ADM2_NAME == "")], paste0(out_dir, "adm2s", pathadd, ".csv"), row.names = F)
  
}

################################################### PREP AGGREGATE FILES FOR GHDX ###########################################################
### Setup
run_date <- <<<< FILEPATH REDACTED >>>>
input_dir <- <<<< FILEPATH REDACTED >>>>
out_dir <- <<<< FILEPATH REDACTED >>>>

### Read results
adm0_results <- fread(<<<< FILEPATH REDACTED >>>>)
adm1_results <- fread(<<<< FILEPATH REDACTED >>>>)
adm2_results <- fread(<<<< FILEPATH REDACTED >>>>)

### Add fields
adm0_results[, c("age_group_id", "age_group_name", "sex_id", "sex", "measure", "metric") := list(22, "All Ages", 3, "Both", "Prevalence", "Rate")]
adm1_results[, c("age_group_id", "age_group_name", "sex_id", "sex", "measure", "metric") := list(22, "All Ages", 3, "Both", "Prevalence", "Rate")]
adm2_results[, c("age_group_id", "age_group_name", "sex_id", "sex", "measure", "metric") := list(22, "All Ages", 3, "Both", "Prevalence", "Rate")]

### Drop pre-2000 results
adm0_results <- adm0_results[year >= 2000]
adm1_results <- adm1_results[year >= 2000]
adm2_results <- adm2_results[year >= 2000]

### Drop results with NA ADM codes
adm0_results <- adm0_results[!is.na(ADM0_CODE)]
adm1_results <- adm1_results[!is.na(ADM0_CODE)]
adm2_results <- adm2_results[!is.na(ADM0_CODE)]

### Reorder and drop columns
adm0_results <- adm0_results[, c("ADM0_CODE", "ADM0_NAME", "year", "age_group_id", "age_group_name", "sex_id", "sex", "measure", "metric", "mean", "lower", "upper")]
adm1_results <- adm1_results[, c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "year", "age_group_id", "age_group_name", "sex_id", "sex", "measure", "metric", "mean", "lower", "upper")]
adm2_results <- adm2_results[, c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "ADM2_CODE", "ADM2_NAME", "year", "age_group_id", "age_group_name", "sex_id", "sex", "measure", "metric", "mean", "lower", "upper")]

### Save files
write.csv(adm0_results, <<<< FILEPATH REDACTED >>>>)
write.csv(adm1_results, <<<< FILEPATH REDACTED >>>>)
write.csv(adm2_results, <<<< FILEPATH REDACTED >>>>)

################################################### PREP RASTERS FOR GHDX ###########################################################
### Setup
run_date <- <<<< FILEPATH REDACTED >>>>
input_dir <- <<<< FILEPATH REDACTED >>>>
out_dir <- <<<< FILEPATH REDACTED >>>>

for (year in 2000:2018) {
  print(year)
  year_mean_current <- raster(<<<< FILEPATH REDACTED >>>>)
  year_lower_current <- raster(<<<< FILEPATH REDACTED >>>>)
  year_upper_current <- raster(<<<< FILEPATH REDACTED >>>>)
  
  writeRaster(year_mean_current, file = (<<<< FILEPATH REDACTED >>>>), overwrite = TRUE, format = "GTiff")
  writeRaster(year_lower_current, file = (<<<< FILEPATH REDACTED >>>>), overwrite = TRUE, format = "GTiff")
  writeRaster(year_upper_current, file = (<<<< FILEPATH REDACTED >>>>), overwrite = TRUE, format = "GTiff")
  
  # Capitalize file extension (writeRaster seems to automatically save files with lowercase file extensions)
  file.rename(<<<< FILEPATH REDACTED >>>>)
  file.rename(<<<< FILEPATH REDACTED >>>>)
  file.rename(<<<< FILEPATH REDACTED >>>>)
}
