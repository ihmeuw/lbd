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
  log_file <- paste0("<<<< FILEPATH REDACTED >>>>", indicator, "_date_log.csv")
  if (!file.exists(log_file)) cat("geomatched_version, collapse_version", file=log_file, append=T)
  cat("\r\n", file=log_file, append=T) # make sure you are on a new line before appending dates
  cat(paste(geomatched_version, collapse_version, sep=","), file = log_file, append = TRUE)

  write.csv(data, file = paste0("<<<< FILEPATH REDACTED >>>>", indicator, ".csv"), row.names=FALSE)
  saveRDS(data, file = paste0("<<<< FILEPATH REDACTED >>>>", indicator, "_", format(Sys.Date(), "%Y_%m_%d"), ".rds"))
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
  data[, c('replicates', 'N_obs') := NULL]
}

####################################################################################################
### PURPOSE: compare passed data for indicator to the previously saved version. Report differences.
### PARAMS: data: data.table containing columns: nid, country, source, year, latitude, longitude,
###               N, <indicator>, weight
###         indicator: name of the indicator that was used to save previous versions
###         prev_version: date of the version to compare to
####################################################################################################
compare_most_recent <- function(data, indicator, prev_version) {
  old <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>", indicator, "_", prev_version, ".rds"))
  both_names <- intersect(names(old), names(data))
  cat(paste0("\n\n*****Columns have been added: \n", paste(setdiff(names(data),names(old)), collapse=", "), "\n\n"))

  old <- old[, both_names, with=F]
  setkeyv(old, key(data))

  cat(paste0("\n\n*****NIDs that have been added:\n", paste(setdiff(data$nid, old$nid), collapse=","), "\n\n"))
  cat(paste0("\n\n*****NIDs that have been dropped:\n", paste(setdiff(old$nid, data$nid), collapse=","), "\n\n"))

  indicator <- gsub("_WN|_MN|_BOTH", "", indicator)

  cat("\n\n*****NIDs that have changed:\n")
  for (nn in intersect(data$nid, old$nid)) {
    temp_old <- old[nid == nn, mget(c("nid", "country", "source", "year", "latitude", "longitude", "N", indicator, "weight"))]
    temp_old[, PREV := get(indicator) / N]
    temp_new <- data[nid == nn, mget(c("nid", "country", "source", "year", "latitude", "longitude", "N", indicator, "weight"))]
    temp_new[, PREV := get(indicator) / N]
    if (all.equal(temp_old, temp_new)[1] != "TRUE") {
      cat(paste("\n\n", nn, "- data changed\n"))
      print(all.equal(data.frame(temp_old), data.frame(temp_new)))
    }
  }

}
