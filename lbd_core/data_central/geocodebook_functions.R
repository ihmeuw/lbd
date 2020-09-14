## get_geocodebooks
#'
#' @title Get geocodebooks
#' @description Pull geocodebooks from database and csvs
#'
#' @details Given a vector of NIDs, pull data from the geocodebooks.
#' Pulls first from the database, then from the old csvs if some
#' nids are unaccounted for. Prints a warning message for surveys not
#' in the database and in neither the database nor the csvs. Requires
#' the RPostgres, DBI, data.table, and dplyr libraries.
#'
#' @param nids A vector of NIDs. If pull_all_nids is true, this value is not used.
#' @param keep_all_columns boolean. Should all rows in the column be kept?
#' if F, returns only necessary columns ("nid", "iso3", "geospatial_id",
#' "point", "lat", "long", "shapefile", "location_code",
#' "location_name", "admin_level", "survey_series")
#' @param pull_all_nids boolean. Pull all nids?
#'
#' @return A data.table of geocodebook data for the nids requested
#'
get_geocodebooks <- function(nids = NULL, keep_all_columns = F, pull_all_nids = F){

  if(is.null(nids) & (!pull_all_nids)){
    stop("You must pass in a vector of nids or set pull_all_nids to True to continue.")
  }

  j <- "<<<< FILEPATH REDACTED >>>>"

  #if RPostgres isn't loaded in, try to load
  if(!("RPostgres" %in% loadedNamespaces())){
    tryCatch(
      {library(RPostgres)},
      error = function(err) {
        tryCatch(
          {library(RPostgres, lib.loc=paste0(j, "<<<< FILEPATH REDACTED >>>>"))},
          error = function(e) {
            stop("Could not source RPostgres library - install and source outside of function")
          })
      }
    )
  }

  if(pull_all_nids){
    nids <- c(-1)
    message("Pulling all nids from the codebooks")
  }

  #query database
  db_table <- get_geocodebook_db(nids, pull_all_nids)
  db_table <- cast_codebook_types(db_table)

  #leftover nids that were not in db
  nids_not_in_db <- nids[!(nids %in% unique(db_table$nid))]

  #if some nids were not found, look in old geocodebook csvs
  if(length(nids_not_in_db) != 0) {
    if(!pull_all_nids){
      message("The following nids are not in the geocodebook database: ", paste(as.character(nids_not_in_db), collapse = ", "))
    }
    csv_table <- get_geocodebook_csvs(nids_not_in_db, j, pull_all_nids)

    #if some nids were in the csvs, match db columns and datatypes
    if(nrow(csv_table) != 0) {
      csv_table <- csv_table[, colnames(db_table), with = FALSE]
      csv_table <- cast_codebook_types(csv_table)
    }

    #check if any surveys were not in either the db or csv
    nids_not_in_csv <- nids_not_in_db[!(nids_not_in_db %in% unique(csv_table$nid))]
    if((length(nids_not_in_csv) != 0) & (!pull_all_nids)) {
      message("The following nids are not in the geocodebook csvs: ", paste(as.character(nids_not_in_csv), collapse = ", "))
    }

    if(pull_all_nids){
      db_nids <- unique(db_table$nid)
      csv_table <- csv_table[!(nid %in% db_nids)]
    }

    #no data in db_table, data in csv
    if(nrow(db_table) == 0){
      combined_table <- csv_table
    #data in both
    } else {
      combined_table <- rbind(db_table, csv_table)
    }
  #data in db
  } else {
    combined_table <- db_table
  }

  #subset to necessary columns
  if(!keep_all_columns){
    geo_keep <- c("nid", "iso3", "geospatial_id", "point", "lat", "long", "shapefile", "location_code", "location_name", "admin_level", "survey_series")
    combined_table <- combined_table[, geo_keep, with=F]
  }

  return(combined_table)
}

## get_geocodebook_db
#'
#' @title Get geocodebooks from geocodebook database
#' @description Pull geocodebooks from database
#'
#' @details Given a vector of NIDs, pull data from the geocodebook database.
#' Used internally in `get_geocodebook`
#'
#' @param nids A vector of NIDs
#' @param pull_all_nids boolean. Pull all nids?
#'
#' @return A data.table of geocodebook data for the nids requested
#'
get_geocodebook_db <- function(nids, pull_all_nids = F) {

  #convert list of nids to string to append to db query
  nids_string = paste(nids, collapse = ", ")

  #make string of columns to keep to append to query
  columns_to_keep <- paste('survey_series_name', 'nid', 'file_path', 'file_name', 'location_id', 'survey_module',
                           'survey_year_start', 'survey_year_end', 'geospatial_id', 'svy_area1', 'svy_area2', 'svy_area3',
                           'point', 'ST_AsText(point_location.snapped_geocoord) AS snapped_geocoord', 'source', 'location_name', 'location_code',  'admin_level',
                           'polygon_location.shapefile_name', 'user_insert', 'user_note', 'loc_name1', 'loc_code1', 'admin_level1',
                           'shapefile1', 'loc_name2', 'loc_code2', 'admin_level2', 'shapefile2', 'loc_name3',
                           'loc_code3', 'admin_level3', 'shapefile3','ST_AsText(point_location.original_geocoord) AS original_geocoord', 'snapped', sep = ", ")

  query = paste("SELECT", columns_to_keep,
                "FROM survey_entry",
                "FULL OUTER JOIN survey_location ON",
                "survey_entry.survey_entry_id = survey_location.survey_entry_id",
                "FULL OUTER JOIN geo_location ON",
                "survey_location.geo_location_id = geo_location.geo_location_id",
                "FULL OUTER JOIN polygon_location ON",
                "geo_location.geo_location_id = polygon_location.geo_location_id",
                "FULL OUTER JOIN shapefile ON",
                "polygon_location.shapefile_name = shapefile.shapefile_name",
                "FULL OUTER JOIN point_location ON",
                "geo_location.geo_location_id = point_location.geo_location_id",
                "FULL OUTER JOIN admin_location ON",
                "survey_location.survey_location_id = admin_location.survey_location_id")

  if(!pull_all_nids){
    query <- paste(query, "Where survey_entry.nid IN (", nids_string, ")")
  }

  #names of columns pulled from db
  column_names = c('survey_series', 'nid', 'file_path', 'file_name', 'iso3', 'survey_module',
                   'start_year', 'end_year', 'geospatial_id', 'svy_area1', 'svy_area2', 'svy_area3',
                   'point', "snapped_geocoord", 'coords_source', 'location_name', 'location_code',  'admin_level',
                   'shapefile', 'user', 'notes', 'loc_name1', 'loc_code1',
                   'admin_level1', 'shapefile1', 'loc_name2', 'loc_code2', 'admin_level2',
                   'shapefile2', 'loc_name3', 'loc_code3', 'admin_level3', 'shapefile3', 'original_geocoord',
                   'snapped')

  # tries to create a connection to the postgres database. If it fails, make an empty data.table to return
  con <- tryCatch({
           con <- dbConnect(RPostgres::Postgres(), dbname = '<<< DBNAME REDACTED >>>',
                            host = '<<< HOST REDACTED >>>', port = 5432,
                            user = '<<< USER REDACTED >>>', password = '<<< PASSWORD REDACTED >>>')},
           error = function(err) {
             message("Could not connect to geocodebook database - skipping")
             return_table <- data.table(matrix(ncol = length(column_names), nrow = 0))
             colnames(return_table) <- column_names
             return_table[, names(return_table) := lapply(.SD, as.character)]
             fix_geocodebook_db_columns(return_table)
           }
         )

  #if there was an error with the db connection, return an empty data.table
  if(is.data.table(con)){
    return(con)
  }

  #query the database, save data as object, then close connection
  res <- dbSendQuery(con, query)
  db_result_dt <- data.table(dbFetch(res))
  dbClearResult(res)
  dbDisconnect(con)

  colnames(db_result_dt) <- column_names

  #if no rows returned from query, setup empty data table to return
  if(nrow(db_result_dt) == 0){
    message("None of the requested nids were in the geocodebook database")
    db_result_dt <- fix_geocodebook_db_columns(db_result_dt)
    return(db_result_dt)
  }

  #convert snapped_geocoord and original_geocoord back into numeric lat/longs
  db_result_dt[, c("lat", "long") := tstrsplit(snapped_geocoord, " ")]
  db_result_dt[, c("original_lat", "original_long") := tstrsplit(original_geocoord, " ")]
  db_result_dt[, lat := gsub("POINT\\(", "", lat)]
  db_result_dt[, original_lat := gsub("POINT\\(", "", original_lat)]
  db_result_dt[, long := gsub("\\)", "", long)]
  db_result_dt[, original_long := gsub("\\)", "", original_long)]
  db_result_dt[, c("original_geocoord", "snapped_geocoord") := NULL]

  point_cols <- c("lat", "long", "original_lat", "original_long")
  db_result_dt[, (point_cols) := lapply(.SD, as.numeric), .SDcols = point_cols]

  #DB encodes rows missing geographic information as points at (0,0) - convert back to NA
  db_result_dt[(lat == 0) & (long == 0), c("lat","long","original_lat", "original_long", "point")] = NA
  return(db_result_dt)
}

## fix_geocodebook_db_columns
#'
#' @title Fix geocodebook db table columns
#' @description drops the snapped_geocoord and original_geocoord columns and replaces them with lat, long, original_lat, and original_long. Used in cases where the database call breaks or returns no data. `get_geocodebook_db` has to return a data.table with the correct column names in order for `get_geocodebooks` to work properly.
#'
#'
#' @param dt A datatable with columns matching the query from the geocodebook db
#'
#' @return A data.table
#'
fix_geocodebook_db_columns <- function(dt){
  dt[, "lat" := numeric()]
  dt[, "long" := numeric()]
  dt[, "original_lat" := numeric()]
  dt[, "original_long" := numeric()]
  dt[, "snapped_geocoord" := NULL]
  dt[, "original_geocoord" := NULL]

  return(dt)
}

## get_geocodebook_csvs
#'
#' @title Get geocodebooks from csvs
#' @description Pull geocodebooks from csvs
#'
#' @details Given a vector of NIDs, pull data from the geocodebook csvs.
#' Used internally in `get_geocodebook`
#'
#' @param nids A vector of NIDs
#' @param j "<<<< FILEPATH REDACTED >>>>"
#' @param pull_all_nids boolean. Pull all nids?
#'
#' @return A data.table of geocodebook data for the nids requested
#'
get_geocodebook_csvs <- function(nids, j, pull_all_nids = F){
  #Read in geographies
  message("Pulling in all geography codebook csvs")
  files <- list.files("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T)
  geogs <- lapply(files, read_add_name_col)
  geo <- rbindlist(geogs, fill=T, use.names=T)

  #dedupe the geography codebook by geospatial_id and nid
  setkey(geo, nid, geospatial_id)
  geo <- unique(geo, use.key=T)

  if(!(pull_all_nids)) {
    geo <- geo[nid %in% nids,]
  }

  return(geo)
}

## read_add_name_col
#'
#' @title Read in and process geocodebook csv
#'
#' @details Given a filepath to the geocodebook csvs, read in file
#' and add a survey_series column based on the name of the file. Casts
#' all columns to character for consistency.
#'
#' @param file A filepath to one of the geocodebook csvs
#'
#' @return A data.table of geocodebook data for the nids requested
#'
read_add_name_col <- function(file){
  rn <- gsub(".csv", "", file, ignore.case=T)
  spl <- strsplit(rn, "/") %>% unlist()
  svy <- spl[length(spl)]
  df <- fread(file, integer64="character", data.table=T, showProgress = F)
  df[, survey_series := svy]
  df <- lapply(df, as.character)
  return(df)
}

## cast_codebook_types
#'
#' @title Cast the type of certain codebook columns
#'
#' @details Given a codebook data.table, cast certain columns
#' to integer or numeric.
#'
#' @param dt A codebook data table prepared from `get_geocodebook_csvs`
#' or `get_geocodebook_db`
#'
#' @return A data.table of geocodebook data with updated types
#'
cast_codebook_types <- function(dt){
  dt[, nid := as.integer(nid)]
  dt[, point := as.integer(point)]
  dt[, location_code := as.integer(location_code)]
  dt[, lat := as.numeric(lat)]
  dt[, long := as.numeric(long)]
  dt[, original_lat := as.numeric(original_lat)]
  dt[, original_long := as.numeric(original_long)]

  return(dt)
}
