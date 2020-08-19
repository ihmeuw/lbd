# ------------ function to load polygon data
# needs to:
# - go to proper directory
# - load dataset
# - subset to polygons only (using 'point' field)

load_polygon_data <- function(use_share = FALSE, indicator, withdate = FALSE, withtag = FALSE, datatag, region = region_list, adm0 = NULL, date = "", update_run_date = FALSE, pathaddin = "_polygons_only") {
  # date stuff if necessary

  # Load input data by indicator
  root <- ifelse(Sys.info()[1] == "Windows", <<<< FILEPATH REDACTED >>>>, <<<< FILEPATH REDACTED >>>>)
  if (use_share == FALSE) load_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  if (use_share == TRUE) load_dir <- <<<< FILEPATH REDACTED >>>>
  if (!withdate & !withtag) d <- fread(file = paste0(<<<< FILEPATH REDACTED >>>>), stringsAsFactors = FALSE)
  if (withtag) d <- fread(file = paste0(<<<< FILEPATH REDACTED >>>>), stringsAsFactors = FALSE)
  if (withdate) d <- fread(file = paste0(<<<< FILEPATH REDACTED >>>>), stringsAsFactors == FALSE)

  # Getting polygons
  d$point <- as.integer(d$point)
  d <- d[point == 0, ]
  d <- d[!is.na(shapefile), ]
  message(paste0("You have ", nrow(d), " rows of data in your dataset."))

  # Subset to region of interest
  d1 <- d[tolower(country) %in% tolower(adm0), ]
  message(paste0("There are ", nrow(d1), " polygons in ", region))
  
  # Save a copy
  if (is.null(update_run_date) == F) {
    if (update_run_date == TRUE) {
      if (dir.exists(paste0(<<<< FILEPATH REDACTED >>>>)) == TRUE) {
        existing_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
        new_try <- existing_dir
        index <- 0
        while (dir.exists(new_try)) {
          index <- index + 1
          new_try <- paste0(existing_dir, "_", index)
        }
        run_date <- paste0(run_date, "_", index)
        dir.create(new_try, showWarnings = FALSE)
        run_date_dir <- new_try
      }
      if (dir.exists(paste0(<<<< FILEPATH REDACTED >>>>)) == FALSE) {
        run_date_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
        dir.create(run_date_dir, showWarnings = FALSE)
      }
      write.csv(d, file = paste0(<<<< FILEPATH REDACTED >>>>))
      return(list(d, run_date))
    }

    if (update_run_date == FALSE) {
      run_date_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
      dir.create(run_date_dir, showWarnings = FALSE)
      write.csv(d1, file = paste0(<<<< FILEPATH REDACTED >>>>))
      return(d1)
    }
  } else {
    return(d1)
  }
}
