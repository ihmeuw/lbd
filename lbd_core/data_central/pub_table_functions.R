################################################################################
##
## GENERATE STANDARD PUBLICATION TABLES
##
## Created: Sept 12, 2018
## Purpose: Functions for generating standard LBD tables that can be easily
##   formatted for publication. These functions use standard LBD datasets as
##   their primary inputs. As of this writing, the following functions are
##   available:
##    * pubtab_source_meta: Generate a source citation table for a publication SI.
##
################################################################################



## pubtab_source_meta --------------------------------------------------------->
#'
#' @title Generate a source metadata table
#' @description Generate an LBD source citation table for a publication SI.
#'
#' @details Given a model-ready MBG dataset, generate a source citation table
#'   that contains critical metadata from the GHDx and sample size information
#'   gleaned from the input dataset.
#'
#' @param indicator Name of the LBD indicator. This name should be related to
#'   a prepped MBG dataset located in the following filepath:
#'   "<<<< FILEPATH REDACTED >>>>".
#' @param core_repo [default="<<<< FILEPATH REDACTED >>>>"] Path to
#'   the core repository, used to load database reading packages.
#' @param resolve_urls [default=TRUE] Should GHDx links to each source be pulled?
#'   (Specifying `TRUE` in this argument will slow down the function a bit)
#' @param out_file [default=NULL] If not `NULL`, saves the prepped source
#'   metadata table to the specified filepath.
#' @param return_output [default=TRUE] If `TRUE`, returns the prepped source
#'   metadata table as a data.table
#'
#' @return If `return_output==TRUE`, returns the prepped source metadata table
#'   as a data.table.
#'
pubtab_source_meta <- function(
  indicator,
  core_repo     = "<<<< FILEPATH REDACTED >>>>",
  resolve_urls  = TRUE,
  out_file      = NULL,
  return_output = TRUE
  ){
  ## Import needed functions
  supp_source <- function(fp) suppressMessages(source(fp))
  supp_source( paste0(core_repo,'/mbg_central/setup.R') )
  suppressMessages( load_R_packages(c('data.table')) )
  supp_source( paste0(core_repo,'/data_central/query_dbs/ghdx_query_functions.R') )

  ## Validate, load, and clean input file
  input_fp <- "<<<< FILEPATH REDACTED >>>>"
  if( !file.exists(input_fp) ){
    stop(paste("File",input_fp,"not found. Check that your indicator is correct!"))
  }
  message("Valid indicator identified. Loading dataset...")
  prepped_dt <- fread(input_fp)
  message("  ...dataset loaded.")
  # Ensure that column names are standardized
  if(!('nid' %in% names(prepped_dt)))       setnames(prepped_dt, 'svy_id',  'nid')
  if(!('latitude' %in% names(prepped_dt)))  setnames(prepped_dt, 'latnum',  'latitude')
  if(!('longitude' %in% names(prepped_dt))) setnames(prepped_dt, 'longnum', 'longitude')
  if(!('location_code' %in% names(prepped_dt))){
    setnames(prepped_dt, 'GAUL_CODE', 'location_code')
  }
  # Get list of all NIDs
  all_nids <- unique(prepped_dt[,nid])

  message("Pulling GHDx metadata for all NIDs...")
  ## Get GHDx metadata for all sources
  ghdx_meta <- ghdx_construct_pub_table(
    nids         = all_nids,
    core_repo    = core_repo,
    resolve_urls = resolve_urls
  )
  message("  ...GHDx data loaded.")

  ## Get the total number of points and polygons by survey
  count_points_polys <- function(in_data){
    # Subset to valid points and polygons
    points <- in_data[(point==1) & !is.null(latitude) & !is.null(longitude),
                      .(nid, latitude, longitude)]
    polys  <- in_data[(point==0) & !is.null(location_code) & !is.null(shapefile),
                      .(nid, location_code, shapefile)]
    # Get counts of unique points and polygons by NID (survey)
    points <- unique(points)
    polys  <- unique(polys)
    points_count <- points[, .(num_points = .N), by=nid]
    polys_count  <-  polys[, .(num_polys  = .N), by=nid]
    # Merge surveys back together and fill NA values with zeroes
    combined <- merge(
      x   = points_count,
      y   = polys_count,
      by  = c('nid'),
      all = T
    )
    combined[is.na(num_points), num_points:=0 ]
    combined[is.na(num_polys),  num_polys :=0 ]
    # Return merged dataset
    return(combined)
  }
  message("Getting geometry data by NID and merging everything together...")
  num_points_polys <- count_points_polys( prepped_dt )

  ## Make sure input data nid field is numeric
  if(is.character(num_points_polys$nid)){
    num_points_polys$nid <- as.numeric(num_points_polys$nid)
  }

  ## Merge valid source data with point/polygon data, keeping only NIDs that are
  ##  PUBLIC and MOST DETAILED (not child NIDs)
  all_combined <- merge(
    x     = ghdx_meta,
    y     = num_points_polys,
    by    = c('nid'),
    all.x = TRUE
  )

  ## CLEANUP
  # Order by country, year, and title
  all_combined <- all_combined[order(iso, years, title)]
  # Add informative message to check when geopositioning has not been added to the GHDx
  all_combined[ is.na(geoprecision), geoprecision := '** UNKNOWN - CHECK **']

  ## Return data in the preferred format (saved or as the function result)
  message("  ...function complete.")
  if(!is.null(out_file)){
    message(paste("Writing results to",out_file))
    write.csv(all_combined, file=out_file, fileEncoding='latin1', row.names=FALSE)
  }
  if(return_output){
    return(all_combined)
  } else {
    return(NULL)
  }
}
