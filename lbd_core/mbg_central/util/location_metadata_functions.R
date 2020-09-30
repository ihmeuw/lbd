library(data.table)
library(DBI)
library(RMySQL)


run_sql_query <- function(conn, query) {
    rows <- DBI::dbGetQuery(conn, query)
    # recommended by "R Cookbook" as a bengin defensive measure, in case MySQL
    # returns an additional result set with status information
    if (RMySQL::dbMoreResults(conn)) DBI::dbNextResults(conn)
    return(rows)
}

get_shared_db_conn <- function() {
    conn <- RMySQL::dbConnect(RMySQL::MySQL(), user = "<<<< USER REDACTED >>>>",
                              password = "<<<< PASSWORD REDACTED >>>>",
                              host = "<<<< HOST REDACTED >>>>",
                              client.flag = RMySQL::CLIENT_MULTI_RESULTS)
    return(conn)
}

# SQL Statements
# Retrieve core fields and path_to_top_parent (used to calculate ihme_lc_id)
all_locations_sql <- paste("SELECT",
                           "location_id         AS loc_id,",
                           "location_name       AS loc_name,",
                           "location_name_short AS loc_nm_sh,",
                           "path_to_top_parent  AS path_to_top_parent",
                           "FROM shared.location")

# Get GAUL_CODE metadata
gaul_code_sql <- paste("SELECT",
                       "location_id             AS loc_id,",
                       "location_metadata_value AS GAUL_CODE",
                       "FROM shared.location_metadata_history",
                       # id 26 is GAUL codes
                       "WHERE location_metadata_type_id = 26",
                       # 19 is 2017.g - the latest version of the metadata before
                       # the location_metadata table was blown away ~25 June
                       "  AND location_metadata_version_id = 19")

# Get ISO3 names (part of ihme_lc_id)
iso_lookup_sql <- paste("SELECT",
                        "location_id             AS loc_id,",
                        "location_metadata_value AS iso_name",
                        "FROM shared.location_metadata_history",
                        # id 1 is "short identifier"/ISO3 code
                        "WHERE location_metadata_type_id = 1",
                        "  AND location_metadata_version_id = 19")

#' Generate location ids
#'
#' Generate ihme_lc_id for a national or subnational location
#'
#' @param iso_lookup a data.table with loc_id and iso_name columns.
#' @param loc_id_hierarchy a comma-delimited string of loc_id values
#'   representing the logical path from EARTH through each admin unit to the
#'   desired location.
#' @return String representing the location id.
get_ihme_lc_id <- function(iso_lookup, loc_id_hierarchy) {
    path <- lapply(strsplit(loc_id_hierarchy, ","), as.numeric)[[1]]
    # first value is always EARTH id; second value is the country id
    nation_id <- path[2]
    nation_name <- iso_lookup[iso_lookup$loc_id == nation_id, "iso_name"]
    # Nations have a location code of their ISO name
    if (length(path) == 2) {
        return(nation_name)
    }
    # Subnationals have a location code of ISO_LOCATIONID
    id <- sprintf("%s_%i", nation_name, path[length(path)])
    # Special case - Puerto Rico
    if (id == "USA_385") {
        id <- "PRI"
    }
    return(id)
}

#' Return location metadata
#'
#' @param shapefile_version string indicating version of shapefile to use. 'gaul' or 'gadm' is inferred from this
#' @param fix_diacritics (default TRUE) whether to strip all diacritic characters from result.
#'
#' Get location metadata for:
#'  IHME location ids (loc_id)
#'  Admin codes (ADM_CODE)
#'  IHME loc ids (ihme_lc_id)
#'  Long location names (loc_name)
#'  Short location names (loc_nm_sh)
#'
#' Note: includes a 'GAUL_CODE' column (same values as ADM_CODE) for compatibility.
#'
#' @return data.frame with above columns
get_location_code_mapping <- function(shapefile_version, remove_diacritics = T) {
    shapefile_type = detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = shapefile_version))$shpfile_type

    if (shapefile_type == 'gaul') {
        data <- get_location_code_mapping_GAUL(remove_diacritics = remove_diacritics)
        data[['ADM_CODE']] <- data[['GAUL_CODE']]
        return(data)
    } else if (shapefile_type == 'gadm') {
        data <- fread("<<<< FILEPATH REDACTED >>>>")
        data[['GAUL_CODE']] <- data[['ADM_CODE']]

        if (remove_diacritics) {
          data$loc_name <- fix_diacritics(data$loc_name)
          data$loc_nm_sh <- fix_diacritics(data$loc_nm_sh)
        }

        return(data)
    } else {
        stop(paste0("Must provide gaul or gadm as shapefile_type, not ", shapefile_type))
    }
}


#' Return location metadata
#'
#' Get location metadata for:
#'  IHME location ids (loc_id)
#'  GAUL codes (GAUL_CODE)
#'  IHME loc ids (ihme_lc_id)
#'  Long location names (loc_name)
#'  Short location names (loc_nm_sh)
#' for locations with available GAUL codes.
#'
#' @return data.frame with above columns
get_location_code_mapping_GAUL <- function(remove_diacritics) {
    conn <- get_shared_db_conn()
    locs <- run_sql_query(conn, all_locations_sql)
    gaul_codes <- run_sql_query(conn, gaul_code_sql)
    gaul_codes['GAUL_CODE'] <- as.numeric(gaul_codes$GAUL_CODE)

    iso_lookup <- run_sql_query(conn, iso_lookup_sql)

    get_id <- function(path_to_parent_str) {
        get_ihme_lc_id(loc_id_helper, path_to_parent_str)  # nolint
    }

    data <- merge(locs, gaul_codes, by = "loc_id")
    data["ihme_lc_id"] <- sapply(data$path_to_top_parent, function(path_str) {
                          get_ihme_lc_id(iso_lookup, path_str)
                        })

    # remove path_to_top_parent
    data <- subset(data, select = c("loc_id",
                                    "loc_name",
                                    "loc_nm_sh",
                                    "ihme_lc_id",
                                    "GAUL_CODE"))
    data <- data.table(data)

    DBI::dbDisconnect(conn)

    # Fix diacritics
    if (remove_diacritics) {
      data$loc_name <- fix_diacritics(data$loc_name)
      data$loc_nm_sh <- fix_diacritics(data$loc_nm_sh)
    }

    # cast to data.table for usability/backwards compatibility
    return(data)
}

## fix_diacritics ################################################

#' Fixes diacritic encoding
#'
#' @param x String / character vector to fix
#' @return String/character vector with replaced diacritics
#' @examples
#' fix_diacritics(c("Côte", "Sør-Trøndela", "Löwenbräu"))

fix_diacritics <- function(x) {

  x <- iconv(x, to = "UTF-8")

  replacement_chars = list('Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 'Ç'='C', 'È'='E', 'É'='E',
                           'Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U',
                           'Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Þ'='B', 'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'ç'='c',
                           'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
                           'ö'='o', 'ø'='o', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y')

  replace_me <- paste(names(replacement_chars), collapse='')
  replace_with <- paste(replacement_chars, collapse = '')

  replaced <- chartr(replace_me, replace_with, x)

  replaced <- gsub('ß', 'Ss', replaced)

  return(replaced)
}
