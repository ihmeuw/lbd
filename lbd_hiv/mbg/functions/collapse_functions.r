####################################################################################################
## Description: Common functions used by multiple collapse codes
####################################################################################################


####################################################################################################
### PURPOSE: find most recent versions of dated files and return most recent date
### PARAMS: dir: directory to look for dated files in
###         date_format: format of dates used in filenames (default: YYYY_MM_DD)
###         out_format: format to return most recent date in (default: YYYY_MM_DD)
###         file_pattern: other pattern to look for in filenames
### RETURNS: date of most recent file in dir
####################################################################################################
most_recent_date <- function(dir, date_format = "%Y_%m_%d", out_format = "%Y_%m_%d", file_pattern = NULL) {
  date_pattern <- gsub("%y|%m|%d", "[[:digit:]]{2}", gsub("%Y", "[[:digit:]]{4}", date_format))
  dates <- dir(dir, pattern = date_pattern)
  if (!is.null(file_pattern)) dates <- grep(file_pattern, dates, value = T)
  dates <- gsub(paste0("(.*)(", date_pattern, ")(.*)"), "\\2", dates)
  dates <- as.Date(dates, date_format)
  format(max(dates), out_format)
}

####################################################################################################
### PURPOSE: save data for mbg with dates of geomatched & collapse. Add dates to log file.
### PARAMS: data: data.table in mbg ready format (only lacking dates)
###         geomatched_version: date of geomatched microdata transformed into data to be saved
###         collapse_version: date collapse was run to create data
####################################################################################################
save_data <- function(data, indicator, geomatched_version, collapse_version) {
  # add dates of geomatched data (geomatched_version) & collapsed data (collapse_version)
  data$geomatched_version <- geomatched_version
  data$collapse_version <- collapse_version

  # write dates to log file
  log_file <- paste0("<<<< FILEPATH REDACTED >>>>")
  if (!file.exists(log_file)) cat("geomatched_version, collapse_version", file=log_file, append=T)
  cat("\r\n", file=log_file, append=T) # make sure you are on a new line before appending dates
  cat(paste(geomatched_version, collapse_version, sep=","), file = log_file, append = TRUE)

  write.csv(data, file = paste0("<<<< FILEPATH REDACTED >>>>")
  saveRDS(data, file = paste0("<<<< FILEPATH REDACTED >>>>"))
}

####################################################################################################
### PURPOSE: check if sample sizes in passed data differ from passed initial start_n sample sizes.
###          Message warnings if so.
### PARAMS: data: data.table containing columns: N_obs, cluster_id, N, nid, country
###         start_n: data.table containing columns: country, nid
###         start (initial n value)
####################################################################################################
check_sample_sizes <- function(data, start_n) {
  data[, replicates := .N, cluster_id]
  end_n <- data[, list(end = sum(N_obs / replicates), end_wt = sum(N)), by='nid,country']
  sample_check <- merge(start_n, end_n, by=c("nid", "country"), all.x=T)
  sample_tmp <- sample_check[start != end,]
  if (sample_check[, sum(start != end)] != 0) message(paste0("sample sizes have changed for ", nrow(sample_tmp), " surveys. Here are the first few: "))
  print(head(sample_tmp))
  if (sample_check[, sum(end_wt > start)] != 0) message("effective sample size is larger than the actual sample size")
  if (sample_check[, sum(end_wt/end < 0.1)] != 0) message("effective sample size is less than 1/10th the actual sample size")
  data[, c('replicates') := NULL]
}

####################################################################################################
### PURPOSE: compare passed data for indicator to the previously saved version. Report differences.
### PARAMS: data: data.table containing columns: nid, country, source, year, latitude, longitude,
###               N, <indicator>, weight
###         indicator: name of the indicator that was used to save previous versions
###         prev_version: date of the version to compare to
####################################################################################################
compare_most_recent <- function(data, indicator, collapse_version, prev_version) {
  
  
  # Compare collapse version to previous version
  new_data <- copy(data)
  old_data <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
  
  # Remove indicator when it has BOTH, MN or WN
  if (indicator == "condom_last_time_BOTH") { 
    indicator_abrev <- "condom_last_time"
  } else if (indicator %in% c("multiple_partners_year_MN", "multiple_partners_year_WN")) {
    indicator_abrev <- "multiple_partners_year"
  } else if (indicator == "in_union_conservative") {
    indicator_abrev <- "in_union"
  } else if (indicator == "had_intercourse_WN") {
    indicator_abrev <- "had_intercourse"
  } else {
    indicator_abrev <- indicator
  }
  
  # Rename to indicator for generalizability 
  setnames(new_data, indicator_abrev, "indicator")
  setnames(old_data, indicator_abrev, "indicator")

  # Write file of changes
  dir.create(paste0("<<<< FILEPATH REDACTED >>>>"), showWarnings = F)
  sink(paste0("<<<< FILEPATH REDACTED >>>>"))
  
  cat(paste0("\n\n Checking differences between geomatched data from current date (", collapse_version, ") to data from ", prev_version)) 
  cat(paste0("\n\n*****NIDs that have been added:\n",   paste(sort(setdiff(new_data$nid, old_data$nid)), collapse=","), "\n\n"))
  cat(paste0("\n\n*****NIDs that have been dropped:\n", paste(sort(setdiff(old_data$nid, new_data$nid)), collapse=","), "\n\n"))
  
  cat("\n\n*****NIDs that have changed:\n")
  for (nn in sort(intersect(old_data$nid, new_data$nid))) {
    temp_old <- old_data[nid == nn, mget(c("nid", "country", "source", "year",
                                           "latitude","longitude","N","indicator","weight", "point"))] %>%
      arrange(nid, country, year, latitude, longitude, N, indicator, weight, point)
    
    temp_new <- new_data[nid == nn, mget(c("nid", "country", "source", "year",
                                           "latitude","longitude","N","indicator","weight", "point"))] %>%
      arrange(nid, country, year, latitude, longitude, N, indicator, weight, point)
    
    # If not equal, run some data checks
    if (!are_equal(temp_old, temp_new)) {
      cat(paste("\n\n", nn, "- data changed\n"))
      
      # Prev old; get rid of missing data
      old_nid <-
        temp_old %>%
        dplyr::summarize(ind_prev = round(100 * weighted.mean(indicator/N, N), 3),
                         sample_size = round(sum(N)))
      
      # Prev new; get rid of missingness for comparison
      new_nid <-
        temp_new %>%
        dplyr::summarize(ind_prev = round(100 * weighted.mean(indicator/N, N), 3),
                         sample_size = round(sum(N)))
      
      if (old_nid$ind_prev != new_nid$ind_prev) {
        cat(paste0("\n\n !! Weighted prevalence changed from ", old_nid$ind_prev, "% to ", new_nid$ind_prev, "%, this should not happen"))
      } else {
        cat("\n\nPrevalence remained the same")
      }
      
      if (old_nid$sample_size != new_nid$sample_size) {
        cat(paste0("\n\n !! Sample size changed from ", old_nid$sample_size, " (old data) to ", new_nid$sample_size, "; (new data) a difference of ", 
                   new_nid$sample_size - old_nid$sample_size, " people (", round(100 * abs(new_nid$sample_size - old_nid$sample_size) / nrow(temp_new), 1), "%)",
                   "\nthis should not happen"))
      } else {
        cat("\n\nSample size remained the same")
      }
      
      # Point/polygon breakdown
      old_points <- data.table(point = c(0, 1), n = c(nrow(filter(temp_old, point == 0)), nrow(filter(temp_old, point == 1))))
      new_points <- data.table(point = c(0, 1), n = c(nrow(filter(temp_new, point == 0)), nrow(filter(temp_new, point == 1))))
      
      if (are_equal(old_points, new_points)) {
        cat(paste0("\n\nGeographic information remained the same: ", new_points[point == 1, n], " people mapped to points and ", new_points[point == 0, n], " mapped to polygons\n\n"))
      } else {
        cat(paste0("\n\nGeographic information changed:\nthere were ", old_points[point == 1, n], " people mapped to points and ", old_points[point == 0, n],
                   " resampled polygons in the old data\nand now there are ", new_points[point == 1, n], " people mapped to points and ", new_points[point == 0, n], " resampled polygons in the new data. Likely due to population raster changes\n\n"))
      }
      
      # Additional checks if lengths are the same but still differences
      if (are_equal(dim(temp_old), dim(temp_new))){
        cat("all.equal function run for additional diagonstic information (current = new data, target = old data)\n\n")
        print(all.equal(temp_old, temp_new))
      }
    }
  }
  sink()
  message(paste0("Done, information saved to <<<< FILEPATH REDACTED >>>>"))
  return()
}