## these functions are used for annual rate of change and projections
## annualized rate of change and projections

## get_cell_pred_for_aroc ################################################

#' Pull cell pred object and format for aroc or projection calculations
#'
#' \code{get_cell_pred_for_aroc()} is a helper function that pulls a 
#' cell pred object and then formats it for use in other functions
#'
#' @param ind_gp indicator group
#' @param ind indicator
#' @param rd run_date
#' @param reg region
#' @param measure prevalence, incidence, mortality, etc
#' @param matrix_pred_name In \code{sprintf} notation. The one object passed into
#'   the string should will be a region name. this allows different regions to be 
#'   passed to different named matrix_preds (pixel level, ad0, ad1, ad2, ...)
#'   e.g. 'had_diarrhea_cell_draws_eb_bin0_%s_diarrhea2_0.RData' which
#'   will be passed to sprintf('had_diarrhea_cell_draws_eb_bin0_%s_0.RData', reg)
#' @param skip_cols columns to skip when reading in the cell preds
#'   For example, if the first two columns store non-pred information in your
#'   file format, \code{skip_cols = 2} will read in all columns from 3 onwards
#' @return a formated \code{cell_pred} object
#' @examples
#' 

get_cell_pred_for_aroc <- function(ind_gp,
                                   ind,
                                   rd,
                                   reg, 
                                   measure,
                                   matrix_pred_name = NULL,
                                   skip_cols = NULL, rk = T) {
  # Load the relevant pred object - loads an object named cell_pred      
  # Try both rds file and rdata file until we standardize this
  rds_file <- sprintf(paste0(<<<< FILEPATH REDACTED >>>>> '_raked_cell_pred.RDs'), 
                               ind_gp, ind, rd)

  if (rk) {
    rdata_file <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/', 
                       ind, "_raked_cell_draws_eb_bin0_", reg, "_0.RData")
  } else {
    rdata_file <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/', 
                       ind, "_cell_draws_eb_bin0_", reg, "_0.RData")
  }
  
  if (!is.null(matrix_pred_name) & file.exists(rds_file)) {
    cell_pred <- readRDS(rds_file)
  }
  
  if (rk) {
    if (file.exists(rdata_file)) {
      load(rdata_file)
      cell_pred <- raked_cell_pred
      rm(raked_cell_pred)
    }
  } else {
    if (file.exists(rdata_file)) {
      load(rdata_file)
    }
  }

  # Check to make sure loaded correctly
  if (!exists("cell_pred")) stop("Unable to load raked cell pred object!")

  # If extra columns at front of cell_pred, can skip here
  if(!(is.null(skip_cols))) cell_pred <- as.matric(cell_pred[, (skip_cols+1):ncol(cell_pred)])

  ## then we need to load the simple raster to determine the
  ## indexing between pixel-years, and matrix rows
  message('-- making simple raster')
  gaul_list <- get_gaul_codes(reg)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                             buffer = 0.4)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']] ## this is what we really need

  ## then we add two columns to the cell_preds. an spatial index
  ## and a year col
  cell_idx <- seegSDM:::notMissingIdx(simple_raster)
  num_yrs  <- nrow(cell_pred) / length(cell_idx)
  ## this should be an integer year. if it isn't, something isn't
  ## synced correctly between preds and simple_rasters
  if(round(num_yrs, 0) != num_yrs){
    print(sprintf("Error! The number of rows in the matrix_pred object for region %s is not divisible by ",
                  "the number of modeling cells in the associated simple raster object that was loaded for ",
                  "the region. Some dimension mismatch has occured and needs to be investigated.",
                  reg))
    break
  }

  new_cols <- cbind(rep(1:num_yrs, each = length(cell_idx)),
                    rep(cell_idx, num_yrs))
  colnames(new_cols) <- c('year', 'idx')
  cell_pred <- cbind(new_cols, cell_pred)

  return(cell_pred)
}

## get_sp_hierarchy ######################################################

#' Pull the hierarchy of admin codes
#'
#' \code{get_sp_hierarchy()} takes the admin2 shapefile and pulls a list of
#' admin level codes to be used in merges, etc. 
#'
#' @return a list of 3 data tables (adm0, adm1, adm2). Each data table includes
#'  the location codes of its administrative level and columns for the higher-
#'  order administrative levels as well.
#' @examples
#' 

get_sp_hierarchy <- function() {

  admin_level <- 2

  # Use most up to date admin2 files - in rds format for faster loading
  admin_shp <- readRDS("<<<< FILEPATH REDACTED >>>>>/g_2015_2014_2_modified.rds")

  # Assigning disputed areas to countries for admin 2 units
  admin_shp[[paste0('ADM', admin_level,'_CODE')]][admin_shp[[paste0('ADM', admin_level,'_CODE')]]==61013]=133
  admin_shp[[paste0('ADM', admin_level,'_CODE')]][admin_shp[[paste0('ADM', admin_level,'_CODE')]]==40760]=40765
  admin_shp[[paste0('ADM', admin_level,'_CODE')]][admin_shp[[paste0('ADM', admin_level,'_CODE')]]==40762]=145
  admin_shp[[paste0('ADM', admin_level,'_CODE')]][admin_shp[[paste0('ADM', admin_level,'_CODE')]]==102  ]=6

  # Pull admin codes
  admin_data <- copy(data.table(admin_shp@data))
  admin_codes <- lapply(0:2, function(aa) {
    col_list <- sapply(0:aa, function(adm) {
      return(c(paste0("ADM", adm, "_CODE"),
               paste0("ADM", adm, "_NAME")))
      })
    col_list <- as.vector(col_list)
    data_subset <- subset(admin_data, select = col_list)
    data_subset <- unique(data_subset)
    return(data_subset)
  })

  names(admin_codes) <- c("ADM0", "ADM1", "ADM2")

  return(admin_codes)

}

## merge_sp_hierarchy ####################################################

#' Takes a data frame or data table with a column of administrative codes
#' and merges on the accompanying admin code names, as well as any higher
#' order administrative codes within which that admin code is nested
#'
#' @param df Data frame or data table with a column of admin codes
#' @param admin_level The administrative level of the codes in the
#'   \code{df} object.  Integer; 0, 1, or 2. 
#' @param idx_col which column in \code{df} contains the administrative
#'   codes?  A character.
#' @param sp_h Spatial hierarchy object generated with 
#'   \code{get_sp_hierarchy()}. If \code{NULL}, then will be generated
#'   within the function. If using this function several times, it may 
#'   save time to generate the sp_h object once, then pass it to this 
#'   function
#' @return Data table with sp_hierarchy added on.  Note that the
#'   original admin_code column will be renamed to conform with
#'   the "ADMX_CODE" convention.
#' @examples
#' # add spatial hierarchy information to df, which contains
#' # admin 2 codes in the column "spatial_idx", using a pre-made
#' # spatial hierarchy object stored in sp_h
#' df2 <- merge_sp_hierarchy(df = df,
#'                           admin_level = 2,
#'                           idx_col = "spatial_idx",
#'                           sp_h = sp_h)

merge_sp_hierarchy <- function(df, admin_level, idx_col, sp_h = NULL) {

  # Grab the spatial hierarchy
  if (is.null(sp_h)) sp_h <- get_spatial_hierarchy()

  df <- as.data.table(df) %>%
        setnames(., idx_col, paste0("ADM", admin_level, "_CODE"))

  # grab df names but not the spatial idx
  df_names <- names(df)[names(df) != paste0("ADM", admin_level, "_CODE")]

  sp_idx_table <- sp_h[[paste0("ADM", admin_level)]]

  df <- merge(sp_idx_table, df, 
               by = c(paste0("ADM", admin_level, "_CODE")), 
               all.y = T, all.x = F) 

  idx_names <- sort(names(sp_idx_table))

  setcolorder(df, c(idx_names, df_names))
  setorderv(df, c(paste0("ADM", 0:admin_level, "_CODE")))

  return(df)
}

## make_aroc #############################################################
#' Generates a set of aroc objects from cell_preds or admin_preds, for use in 
#' projections or in other analyses that require AROC
#'
#' @param ind_gp indicator group
#' @param ind indicator
#' @param rd run_date
#' @param matrix_pred_name In \code{sprintf} notation. The one object passed into
#'   the string should will be a region name. this allows different regions to be 
#'   passed to different named matrix_preds (pixel level, ad0, ad1, ad2, ...)
#'   e.g. 'had_diarrhea_cell_draws_eb_bin0_%s_diarrhea2_0.RData' which
#'   will be passed to sprintf('had_diarrhea_cell_draws_eb_bin0_%s_0.RData', reg)
#' @param type Type of aroc to create. Options include \code{cell}, \code{admin},
#'   or \code{c('cell', 'admin')} for both.
#' @param measure prevalence, incidence, mortality, etc
#' @param skip_cols columns to skip when reading in the cell preds
#'   For example, if the first two columns store non-pred information in your
#'   file format, \code{skip_cols = 2} will read in all columns from 3 onwards
#' @param year_list Vector (integer) of years included in the model run / cell pred object
#' @param uselogit Should this be done in logit space?
#' @return writes cell- and admin- level aroc objects to standard directories and formats
#'   in the 'pred_derivatives' folder of the model run.  cell-level objects are in the
#'   cell_pred indexed format. Both cell- and admin- aroc objects are matrices wide by draw.
#' @examples
#' make_aroc(ind_gp = indicator_group, 
#'           ind = indicator,
#'           rd = run_date,
#'           matrix_pred_name = NULL,
#'           type = c("cell", "admin"),
#'           measure = "prevalence",
#'           year_list = c(2000:2015),
#'           uselogit = TRUE)


make_aroc <- function(ind_gp, ind, rd, 
                      matrix_pred_name = NULL, 
                      type, 
                      measure="mortality",
                      skip_cols = NULL,
                      year_list = c(2000:2015),
                      uselogit = FALSE,
                      raked) {

  # define directories
  share_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/')
  output_dir <- paste0(share_dir, "/pred_derivatives/aroc/")
  dir.create(output_dir, showWarnings=F, recursive = T)

  ## load in regions used in this model and run_date
  regions <- get_output_regions(in_dir =share_dir)

  ## get number of years
  num_yrs <- length(year_list)

  if('cell' %in% type){
    message('Working on CELL level')
    for(ii in 1:length(regions)){
      message(sprintf('- On region: %s', regions[ii]))
      
      cell_pred <- get_cell_pred_for_aroc(ind_gp,
                                          ind,
                                          rd, 
                                          regions[ii], 
                                          measure,
                                          matrix_pred_name,
                                          skip_cols, rk = raked)
                  
      ## now we have some year column and some id column and the rest are draw columns.
      ## we can now calculate AROC by idx across time
      num_draws <- ncol(cell_pred) - 2
      num_idx <- length(unique(cell_pred[, 2])) ## second col is 'idx'
      
      ## for each draw, calculate the rates of change between years and
      ## then the weighted total AROC across all years
      message('-- making AROC between years for each draw')
      aroc_draws <- matrix(ncol = num_draws, nrow = num_idx)
      #goal_draws <- copy(aroc_draws)
      if(uselogit){
        cell_pred <- as.data.table(cell_pred)
        for(jj in names(cell_pred)[-(1:2)]){ set(cell_pred,which(cell_pred[[jj]] >= 1),jj,1-1e-9)} ## for stuff raked abonve 1!
        ## cell_pred[which(cell_pred >= 1), ] <- 1 - 1e-9 ## for stuff raked abonve 1!
        ## cell_pred_logit[, 3:ncol(cell_pred)] <- log(cell_pred[, 3:ncol(cell_pred)] / (1 - cell_pred[, 3:ncol(cell_pred)]))
        cell_pred[, names(cell_pred)[-(1:2)] := lapply(.SD, logit), .SDcols = 3:1002]
        cell_pred_logit <- as.matrix(cell_pred)
      }
      for(dd in 1:(ncol(cell_pred) - 2)){
        if(dd %% 50 == 1) message(sprintf('---- on draw %i out of %i', dd, num_draws))
        aroc_mat <- matrix(ncol = num_yrs - 1, nrow = num_idx)
        for(yy in 1:(num_yrs - 1)){
          if(uselogit){
            aroc_mat[, yy] <- cell_pred_logit[1:num_idx + (yy) * num_idx, dd + 2] -  cell_pred_logit[1:num_idx + (yy - 1) * num_idx, dd + 2]
          }else{
            aroc_mat[, yy] <- log(cell_pred[1:num_idx + (yy) * num_idx, dd + 2]) -  log(cell_pred[1:num_idx + (yy - 1) * num_idx, dd + 2])
          }
        }
        year_wt <- 1:(num_yrs - 1)/ sum(1:(num_yrs - 1)) ## linear weighting
        aroc_vec <- aroc_mat %*% year_wt
        aroc_draws[, dd] <- aroc_vec
        #goal_draws[, dd] <- ifelse(aroc_vec < aroc_goal, 1, 0)  
      }
      # short_goal <- ifelse(aroc_draws < aroc_goal_2015, 1, 0)
      # goal_draws[is.na(goal_draws)] <- 0
      # short_goal[is.na(short_goal)] <- 0
      # relative_goal_prob <- rowMeans(goal_draws)
      # achieved_relative_prob <- rowMeans(short_goal)
      
      message('-- finished making AROC across draws. now saving')
      saveRDS(object = aroc_draws,
              file = sprintf('%s/%s_%s_aroc_cell_draw_matrix_%s%s.RDs',
                             output_dir, ind, measure, ifelse(uselogit, "logit_", ""), regions[ii]))
    } # Close region loop
  } # if ('cell' %in% type)

  if('admin' %in% type){
    
    message('Working on ADMIN level')
    ## load the admin objects
    ## try two different locations until we standardize
    file_1 <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s_%s_admin_draws_raked.Rdata',
                 ind_gp, ind, rd, ind, measure)
    file_2 <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/', 
                     ind, '_raked_admin_draws_eb_bin0_0.RData')
 
    if (file.exists(file_1)) {
      load(file_1)
    } else if (file.exists(file_2)) {
      load(file_2)
    } else {
      stop("Cannot load admin pred object!")
    }

    ## this contains admin_0, admin_1, and admin_2 matrix_draw objects
    ## col1 is year, col2 is ADM*_CODE, final column is pop. all other columns are draws
    
    for(aa in 0:2){ ## for each admin type
      message(sprintf('- On admin: %i', aa))
      
      cell_pred <- get(sprintf('admin_%i', aa))
      cell_pred <- as.data.table(cell_pred)
      
      ## format this object to look like the one we used for cell level (easier to copy code from above)
      ## revised now to work with a bigger array of cell pred objects
      str_match <- stringr::str_match
      draw_cols <- names(cell_pred)[grep("V[0-9]*", names(cell_pred))]
      keep_cols <- c("year", paste0("ADM", aa, "_CODE"), draw_cols)
      cell_pred <- subset(cell_pred, select = keep_cols)
      setnames(cell_pred, paste0("ADM", aa, "_CODE"), "idx")
      cell_pred <- as.matrix(cell_pred)
      
      num_draws <- ncol(cell_pred) - 2

      # Redid the below line to warn if mismatches in case of non-unique admin codes - JM
      num_idx <- nrow(cell_pred) / length(year_list)
      if (num_idx != length(unique(cell_pred[,2]))) {
        warning(paste0("The number of unique cell_pred indices does not equal the number of ",
                       "rows of cell_pred divided by the number of years.",
                       "\nYou may have duplicate admin codes at the admin ", aa, 
                       " level - check your aggregated objects!"))
      }
      num_yrs <- length(unique(cell_pred[, 1]))
      idx <- cell_pred[which(cell_pred[, 1] == min(cell_pred[, 1])), 2] ## idx in year1
      
      ## for each draw, calculate the rates of change between years and
      ## then the weighted total AROC across all years
      message('-- making AROC between years for each draw')
      aroc_draws <- matrix(ncol = num_draws, nrow = num_idx)

      if(uselogit){
        cell_pred_logit <- cell_pred
        cell_pred_logit[, 3:ncol(cell_pred)] <- log(cell_pred[, 3:ncol(cell_pred)] / (1 - cell_pred[, 3:ncol(cell_pred)]))
      }
      
      for(dd in 1:(ncol(cell_pred) - 2)){
        if(dd %% 50 == 1) message(sprintf('---- on draw %i out of %i', dd, num_draws))
        aroc_mat <- matrix(ncol = num_yrs - 1, nrow = num_idx)
        for(yy in 1:(num_yrs - 1)){
          if(uselogit){
            aroc_mat[, yy] <- cell_pred_logit[1:num_idx + (yy) * num_idx, dd + 2] -  cell_pred_logit[1:num_idx + (yy - 1) * num_idx, dd + 2]
          }else{
            aroc_mat[, yy] <- log(cell_pred[1:num_idx + (yy) * num_idx, dd + 2]) -  log(cell_pred[1:num_idx + (yy - 1) * num_idx, dd + 2])
          }
        }
        year_wt <- 1:(num_yrs - 1)/ sum(1:(num_yrs - 1)) ## linear weighting
        aroc_vec <- aroc_mat %*% year_wt
        aroc_draws[, dd] <- aroc_vec
      }
      
      message('-- finished making AROC across draws. now saving')
      final_aroc <- cbind(idx, aroc_draws)
      colnames(final_aroc)[1] <- sprintf('ADM%i_CODE', aa)
      
      saveRDS(object = final_aroc,
              file = sprintf('%s/%s_%s_aroc_adm%i_draw_matrix%s.RDs',
                             output_dir, ind, measure, aa, ifelse(uselogit, "_logit", "")))

    } # For aa in admin
  } # if ('admin' %in% type)...
}

## make_proj #############################################################

#' Generates a set of draw-level projection objects from aroc objects for 
#' a given set of target years
#'
#' @param ind_gp indicator group
#' @param ind indicator
#' @param rd run_date
#' @param type Type of aroc to create. Options include \code{cell}, \code{admin},
#'   or \code{c('cell', 'admin')} for both.
#' @param proj_years Vector (integer) of years that you want to project to.  Note 
#'   that this is different from \code{year_list}, which is the list of years that 
#'   were included in the model run / are included in the aroc object.
#' @param measure prevalence, incidence, mortality, etc
#' @param skip_cols columns to skip when reading in the cell preds
#'   For example, if the first two columns store non-pred information in your
#'   file format, \code{skip_cols = 2} will read in all columns from 3 onwards
#' @param matrix_pred_name In \code{sprintf} notation. The one object passed into
#'   the string should will be a region name. this allows different regions to be 
#'   passed to different named matrix_preds (pixel level, ad0, ad1, ad2, ...)
#'   e.g. 'had_diarrhea_cell_draws_eb_bin0_%s_diarrhea2_0.RData' which
#'   will be passed to sprintf('had_diarrhea_cell_draws_eb_bin0_%s_0.RData', reg)
#' @param year_list Vector (integer) of years included in the model run / cell pred object
#' @param uselogit Should this be done in logit space?
#' @return writes cell- and admin- level projection objects to standard directories and formats
#'   in the 'pred_derivatives' folder of the model run.  cell-level objects are in the
#'   cell_pred indexed format. Both cell- and admin- projection objects are matrices wide by draw.
#' @examples
#' make_proj(ind_gp = indicator_group, 
#'           ind = indicator, 
#'           rd = run_date, 
#'           type = c("cell", "admin"),
#'           proj_years = c(2020, 2025, 2030),
#'           measure = "prevalence",
#'           skip_cols = NULL,
#'           year_list = c(2000:2015), 
#'           uselogit = TRUE)

make_proj <- function(ind_gp, ind, rd, 
                      type, 
                      proj_years = c(2020, 2025), 
                      measure="mortality", 
                      skip_cols = NULL,
                      matrix_pred_name = NULL,
                      year_list = c(2000:2015),
                      uselogit = FALSE) {

  ## load in regions used in this model and run_date
  regions <- get_output_regions(in_dir = paste0('<<<< FILEPATH REDACTED >>>>>',
                                                ind_gp, '/',
                                                ind, '/output/',
                                                rd))

  # define directories
  share_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/')
  aroc_dir <- paste0(share_dir, "/pred_derivatives/aroc/")
  output_dir <- paste0(share_dir, "/pred_derivatives/proj/")

  dir.create(output_dir, recursive = T, showWarnings = F)

  ## make projections at the cell level
  if ('cell' %in% type) {
    message('Working on CELL level')
    for(ii in 1:length(regions)){
      message(sprintf('- On region: %s', regions[ii]))
      
      # Pull cell pred
      cell_pred <- get_cell_pred_for_aroc(ind_gp,
                                     ind,
                                     rd, 
                                     regions[ii], 
                                     measure,
                                     matrix_pred_name,
                                     skip_cols)
      
      # Load aroc object
      message('-- loading aroc')
      
      aroc_draws <- readRDS(sprintf('%s/%s_%s_aroc_cell_draw_matrix_%s%s.RDs',
                             aroc_dir, ind, measure, ifelse(uselogit, "logit_", ""), regions[ii]))

      # Get number of draws
      num_draws <- ncol(aroc_draws)
      if (ncol(aroc_draws) != ncol(cell_pred) - 2) stop("cell_pred & aroc draw #s do not match!")

      # Make projections -----------------------------------
      ## also make projections and save them
      message('-- making projections')
      
      ## grab last year of modeled estiamtes
      final_year <- cell_pred[which(cell_pred[, 1] == max(cell_pred[, 1])), ] ## first col is year
      last_year <- final_year[, -(1:2)]
      
      ## grab idx 
      idx <- final_year[,2]
      
      ## unlist all draws into vector, apply forecast, and convert back to matrix
      final_year <- as.vector(final_year[, -(1:2)]) ## unlist only the draw columns
      final_year[which(final_year) >= 1] <- 1 - 1e-9 ## for pixels raked above 1
      
      ## unlist aroc draw matrix
      aroc_draws <- as.vector(aroc_draws)

      ## set up a list to capture the projection output
      proj_draws_list <- list()

      for (yr in proj_years) {

        message(paste0('--- year: ', yr))

        ## figure out how many years to project
        proj_dur <- as.numeric(yr) - max(year_list)

        ## make projection from final_yr out proj_dur years
        if(uselogit){
          proj_draws_logit <- log(final_year / (1 - final_year)) + (aroc_draws * proj_dur)
          proj_draws <- exp(proj_draws_logit) / (1 + exp(proj_draws_logit))
        }else{
          proj_draws <- final_year * exp(aroc_draws * proj_dur)
        }
        
        ## convert back to matrix
        proj_draws <- matrix(proj_draws, ncol = num_draws)

        ## insert into list
        proj_draws_list[[as.character(yr)]] <- proj_draws
        rm(proj_draws)
         
      }
       
      ## save
      message('-- saving projections')

      lapply(1:length(proj_draws_list), function(i) {
        saveRDS(object = proj_draws_list[[i]],
                file = sprintf('%s/%s_%s_%s_projections_cell_draw_matrix_%s%s.RDs',
                               output_dir, ind, measure, names(proj_draws_list)[i],
                               ifelse(uselogit, "logit_", ""), regions[ii]))
        }) 
    } # close regions loop
  } # if ('cell' %in% type)

  # make projections as the admin level
  if ('admin' %in% type) {

    message('Working on ADMIN level')
    ## load the admin objects
    ## try two different locations until we standardize
    message('- loading admin objects')
    file_1 <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s_%s_admin_draws_raked.Rdata',
                 ind_gp, ind, rd, ind, measure)
    file_2 <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/', 
                     ind, '_raked_admin_draws_eb_bin0_0.RData')
 
    if (file.exists(file_1)) {
      load(file_1)
    } else if (file.exists(file_2)) {
      load(file_2)
    } else {
      stop("Cannot load admin pred object!")
    }

    for(aa in 0:2){ ## for each admin type 

      message(paste0("- admin level: ", aa))

      # create pseudo cell-pred object
      cell_pred <- get(sprintf('admin_%i', aa))
      cell_pred <- as.data.table(cell_pred)
      
      ## format this object to look like the one we used for cell level (easier to copy code from above)
      ## revised now to work with a bigger array of cell pred objects
      str_match <- stringr::str_match
      draw_cols <- names(cell_pred)[grep("V[0-9]*", names(cell_pred))]
      keep_cols <- c("year", paste0("ADM", aa, "_CODE"), draw_cols)
      cell_pred <- subset(cell_pred, select = keep_cols)
      setnames(cell_pred, paste0("ADM", aa, "_CODE"), "idx")
      cell_pred <- as.matrix(cell_pred)
    
      # Redid the below line to warn if mismatches in case of non-unique admin codes - JM
      num_idx <- nrow(cell_pred) / length(year_list)
      if (num_idx != length(unique(cell_pred[,2]))) {
        warning(paste0("The number of unique cell_pred indices does not equal the number of ",
                       "rows of cell_pred divided by the number of years.",
                       "\nYou may have duplicate admin codes at the admin ", aa, 
                       " level - check your aggregated objects!"))
      }
      num_yrs <- length(unique(cell_pred[, 1]))
      idx <- cell_pred[which(cell_pred[, 1] == min(cell_pred[, 1])), 2] ## idx in year1
      
      # load aroc
      message("-- loading aroc") 
      aroc_draws <- readRDS(sprintf('%s/%s_%s_aroc_adm%i_draw_matrix%s.RDs',
                             aroc_dir, ind, measure, aa, ifelse(uselogit, "_logit", "")))

      # remove first column (spatial index) and store separately
      spatial_idx <- aroc_draws[, 1]
      aroc_draws <- aroc_draws[, 2:ncol(aroc_draws)]

       # Get number of draws
      num_draws <- ncol(aroc_draws)
      if (ncol(aroc_draws) != ncol(cell_pred) - 2) stop("cell_pred & aroc draw #s do not match!")

      # make projections
      message('-- making projections')
      
      ## grab last year of modeled estiamtes
      final_year <- cell_pred[which(cell_pred[, 1] == max(cell_pred[, 1])), ] ## first col is year
      last_year <- final_year[, -(1:2)]
      
      ## grab idx 
      idx <- final_year[,2]
      
      ## unlist all draws into vector, apply forecast, and convert back to matrix
      final_year <- as.vector(final_year[, -(1:2)]) ## unlist only the draw columns
      
      ## unlist aroc draw matrix
      aroc_draws <- as.vector(aroc_draws)

      proj_draws_list <- list()
      
      ## make projection from final_yr out proj_dur years
      for (yr in proj_years) {

        message(paste0('--- year: ', yr))

        ## figure out how many years to project
        proj_dur <- as.numeric(yr) - max(year_list)

        if(uselogit){
          proj_draws_logit <- log(final_year / (1 - final_year)) + (aroc_draws * proj_dur)
          proj_draws <- exp(proj_draws_logit) / (1 + exp(proj_draws_logit))
        }else{
          proj_draws <- final_year * exp(aroc_draws * proj_dur)
        }
          
        ## convert back to matrix
        proj_draws <- matrix(proj_draws, ncol = num_draws)

        ## append spatial index
        proj_draws <- cbind(spatial_idx, proj_draws)

        ## insert into list
        proj_draws_list[[as.character(yr)]] <- proj_draws
        rm(proj_draws)
         
      }   
    
      ## save
      message('-- saving projections')
      lapply(1:length(proj_draws_list), function(i) {
        saveRDS(object = proj_draws_list[[i]],
                file = sprintf('%s/%s_%s_%s_projections_adm%i_draw_matrix%s.RDs',
                               output_dir, ind, measure, names(proj_draws_list)[i],
                               aa, ifelse(uselogit, "_logit", "")))
        })
    } # close admin levels loop 
   } # if ('admin' %in% type)

}

## add_goal ##############################################################

#' Adds a goal to a goal object, or creates a new one
#'
#' This function creates or modifies a "goal object" - which is just a data table
#' that contains all of the information in each row to check a given target.  This
#' setup is useful in case you want to check a series of different targets, or 
#' run the same analysis for a series of target years.  (Note that all target
#' years must have projection objects created by first running \code{make_proj()})
#'
#' @param goal_obj An existing goal object. If \{NULL}, will create a new goal object.
#' @param target_year Year for which the goal is meant. Must be a future (projected) year
#'   and the year projections must already be generated for this target_year using 
#'   \code{make_proj()}
#' @param target Target (e.g. 0.8 for 80% vaccine coverage); numeric
#' @param target_type Type of target ('greater' or 'less'); character
#' @param abs_rel Is this an 'absolute' or 'relative' goal? 
#' @param baseline_year For relative goals, what year should we compare to?
#' @param proj Placeholder - for now, always use projected years.  In the future this
#'   may support modeled years as well, but still need to code
#' @param pred_type Type of aroc to create. Options include \code{cell}, \code{admin},
#'   or \code{c('cell', 'admin')} for both.
#' @return A "goal object" (data.table in a specific format for use in 
#'   \code{compare_to_target()})
#' @examples
#' # Define goals: start by initializing goal object
#' goals <- add_goal(target_year = 2030, 
#'                   target = 0.8,
#'                   target_type = "greater",
#'                   abs_rel = "absolute",
#'                   pred_type = c("cell", "admin"))
#'
#' # Add goals to existing goal object by specifying goal_obj
#' goals <- add_goal(goal_obj = goals,
#'                   target_year = 2020, 
#'                   target = 0.8,
#'                   target_type = "greater",
#'                   abs_rel = "absolute",
#'                   pred_type = c("cell", "admin"))

add_goal <- function(goal_obj = NULL,
                     target_year,
                     target,
                     target_type,
                     abs_rel,
                     baseline_year = NA,
                     proj = TRUE,
                     pred_type = c("cell", "admin")) {

  if (!(abs_rel %in% c("absolute", "relative"))) {
    stop("abs_rel must be either \"absolute\" or \"relative\"")
  }

  if (!(target_type %in% c("greater", "less"))) {
    stop("target_type must be either \"greater\" or \"less\"")
  }

  if (abs_rel == "relative" & is.na(baseline_year)) {
    stop("Must specify baseline_year if abs_rel == relative")
  } else if (abs_rel == "absolute" & !is.na(baseline_year)) {
    warning("For absolute targets, ignoring baseline_year")
  }

  if (length(pred_type) > 1) {
    pred_type <- paste0("c(\"",
                        paste(pred_type, collapse = "\", \""),
                        "\")")
  }

  new_goal <- data.table(
                         target_year = target_year,
                         target = target,
                         target_type = target_type,
                         abs_rel = abs_rel,
                         baseline_year = baseline_year,
                         proj = proj, 
                         pred_type = as.character(pred_type))

  if (!is.null(goal_obj)) new_goal <- rbind(goal_obj, new_goal)
  return(new_goal)
}

## compare_to_target #####################################################

#' Runs comparison of projected years against a set of goals defined in a goal object
#'
#' @param ind_gp indicator group
#' @param ind indicator
#' @param rd run_date
#' @param goal_obj An existing goal object made with \code{add_goal()}
#' @param measure prevalence, incidence, mortality, etc
#' @param year_list Vector (integer) of years included in the model run / cell pred object
#' @param uselogit Should this be done in logit space?
#' @param skip_cols columns to skip when reading in the cell preds
#'   For example, if the first two columns store non-pred information in your
#'   file format, \code{skip_cols = 2} will read in all columns from 3 onwards
#' @param matrix_pred_name In \code{sprintf} notation. The one object passed into
#'   the string should will be a region name. this allows different regions to be 
#'   passed to different named matrix_preds (pixel level, ad0, ad1, ad2, ...)
#'   e.g. 'had_diarrhea_cell_draws_eb_bin0_%s_diarrhea2_0.RData' which
#'   will be passed to sprintf('had_diarrhea_cell_draws_eb_bin0_%s_0.RData', reg)
#' @return generates cell- and admin- level probabilities of meeting the specified goal,
#'   according to what is specified in \code{goal_obj}.  Objects are written to 
#'   standard directories and formats in the 'pred_derivatives' folder of the model run.  
#'   cell-level objects are in the cell_pred indexed format, but are no longer
#'   wide by draw (just a single column as they are probabilities).  admin objects are
#'   saved both as rds files and as .csv files with admin hierarchy appended. 
#' @examples
#' # Define goals: start by initializing goal object
#' goals <- add_goal(target_year = 2030, 
#'                   target = 0.8,
#'                   target_type = "greater",
#'                   abs_rel = "absolute",
#'                   pred_type = c("cell", "admin"))
#'
#' # Run comparison to target
#' compare_to_target(ind_gp = indicator_group, 
#'                   ind = indicator, 
#'                   rd = run_date, 
#'                   goal_obj = goals,
#'                   measure = "prevalence", 
#'                   year_list = c(2000:2015),
#'                   uselogit = T)

compare_to_target <- function(ind_gp, 
                              ind, 
                              rd,
                              goal_obj,
                              measure,
                              year_list = c(2000:2015),
                              uselogit,
                              skip_cols = NULL,
                              matrix_pred_name = NULL){

  ## load in regions used in this model and run_date
  regions <- get_output_regions(in_dir = paste0('<<<< FILEPATH REDACTED >>>>>',
                                                ind_gp, '/',
                                                ind, '/output/',
                                                rd))


  # define directories
  share_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/')
  aroc_dir <- paste0(share_dir, "/pred_derivatives/aroc/")
  proj_dir <- paste0(share_dir, "/pred_derivatives/proj/")
  output_dir <- paste0(share_dir, "/pred_derivatives/target_probs/")

  dir.create(output_dir, showWarnings = F, recursive = T)


  for (n in 1:nrow(goal_obj)) {

    # Load parameters for this row
    target_year <- as.integer(goal_obj[n, target_year])
    target <- as.numeric(goal_obj[n, target])
    target_type <- as.character(goal_obj[n, target_type])
    abs_rel <- as.character(goal_obj[n, abs_rel])
    if (abs_rel == "relative") baseline_year <- as.integer(goal_obj[n,baseline_year])
    proj <- as.logical(goal_obj[n, proj])
    pred_type <- as.character(goal_obj[n, pred_type])
    goal_type <- as.character(goal_obj[n, target_type])
    message(paste0("\n------------------------------------------------",
                   "\nWorking on target comparison for:", 
                   "\n  target_year: ", target_year, 
                   "\n  target: ", target, 
                   "\n  target_type: ", target_type, 
                   "\n  abs_rel: ", abs_rel, 
                   ifelse(abs_rel == "relative", paste0("\n  baseline_year: ",baseline_year), ""), 
                   "\n  proj: ", proj, 
                   "\n  pred_type: ", goal_obj[n, pred_type]), "\n") # character version for pred_type

    if(proj == F) stop("For now, can only use projected values")

    if ('cell' %in% pred_type) {
      message('Working on CELL level')

      for(ii in 1:length(regions)){
        message(sprintf('- On region: %s', regions[ii]))

        # Load proj_draws object
        message("-- Loading proj_draws object...")
        proj_draws <- readRDS(sprintf('%s/%s_%s_%s_projections_cell_draw_matrix%s_%s.RDs',
                                      proj_dir, ind, measure, target_year, 
                                      ifelse(uselogit, "_logit",""), regions[ii]))

        # Generate probabilities
        if (abs_rel == "absolute") {

          # Calculate cell-wise probability of meeting the goal
          if (target_type == "greater") {
            absolute_goal_draws <- ifelse(proj_draws >= target, 1, 0)
          } else if (target_type == "less") {
            absolute_goal_draws <- ifelse(proj_draws <= target, 1, 0)
          }

          absolute_goal_draws[is.na(absolute_goal_draws)] <- 0
          absolute_goal_prob <- rowMeans(absolute_goal_draws)

          # Save
          message("-- Saving probabilities...")
          saveRDS(object = absolute_goal_prob, 
                  file = paste0(output_dir, ind, "_", measure, "_", target_year, "_", abs_rel, "_",
                                target_type, '_', target, '_cell_target_probs_', regions[ii], '.RDs'))

          rm(absolute_goal_prob, absolute_goal_draws)

        } else if (abs_rel == "relative") {

          # Need to load the pred files for this one...
          # Pull cell pred
          cell_pred <- get_cell_pred_for_aroc(ind_gp,
                                              ind,
                                              rd, 
                                              regions[ii], 
                                              measure,
                                              matrix_pred_name,
                                              skip_cols)
                
          
          ## grab baseline year preds
          year_idx = which(year_list == baseline_year)
          baseline_year_draws <- cell_pred[which(cell_pred[, 1] == year_idx), ] ## first col is year 
          baseline_year_draws <- baseline_year_draws[, -(1:2)]

          # Now do the comparisons
          if (target_type == "greater") {
              relative_proj_draws <- ifelse(proj_draws / baseline_year_draws >= 1 - target, 1, 0)
            } else if (target_type == "less") {
              relative_proj_draws <- ifelse(proj_draws / baseline_year_draws <= 1 - target, 1, 0)
            }
          relative_proj_draws[is.na(relative_proj_draws)] <- 0
          relative_goal_prob <- rowMeans(relative_proj_draws)

          # Save
          message("-- Saving probabilities...")
          saveRDS(object = relative_goal_prob, 
                  file = paste0(output_dir, ind, "_", measure, "_", target_year, "_vs_", 
                                baseline_year, "_", abs_rel, "_", target_type, '_', target, 
                                '_cell_target_probs_', regions[ii], '.RDs'))

          rm(relative_goal_prob, relative_proj_draws)
        
        } # close abs/relative if/then/else 
      } #close regions loop
    } # close cell loop

    if ('admin' %in% pred_type) {
      message('Working on ADMIN level')

      # Load sp hierarchy
      sp_h <- get_sp_hierarchy()

      for(aa in 0:2) { ## for each admin type 

        # Load proj_draws object
        message("-- Loading proj_draws object...")
        proj_draws <- readRDS(sprintf('%s/%s_%s_%s_projections_adm%i_draw_matrix%s.RDs',
                                      proj_dir, ind, measure, target_year, aa,
                                      ifelse(uselogit, "_logit","")))

        # Split off spatial index
        spatial_idx <- proj_draws[,1]
        proj_draws <- proj_draws[, 2:ncol(proj_draws)]

        # Generate probabilities
        if (abs_rel == "absolute") {
          if (target_type == "greater") {
            absolute_goal_draws <- ifelse(proj_draws >= target, 1, 0)
          } else if (target_type == "less") {
            absolute_goal_draws <- ifelse(proj_draws <= target, 1, 0)
          }

          absolute_goal_draws[is.na(absolute_goal_draws)] <- 0
          absolute_goal_prob <- rowMeans(absolute_goal_draws)

          # Add spatial index
          absolute_goal_prob <- cbind(spatial_idx, absolute_goal_prob)
          
          # Save
          message("-- Saving probabilities...")

          # RDS
          saveRDS(object = absolute_goal_prob, 
                  file = paste0(output_dir, ind, "_", measure, "_", target_year, "_", abs_rel, "_",
                                target_type, '_', target, '_adm_', aa, '_target_probs.RDs'))
          
          # CSV
          absolute_goal_prob <- merge_sp_hierarchy(df = absolute_goal_prob,
                                                   admin_level = aa,
                                                   idx_col = "spatial_idx",
                                                   sp_h = sp_h)

          write.csv(absolute_goal_prob,
                    file = paste0(output_dir, ind, "_", measure, "_", target_year, "_", abs_rel, "_",
                                target_type, '_', target, '_adm_', aa, '_target_probs.csv'),
                    row.names = F)

        }
        
        ## Did it meet relative goal? 
        if (abs_rel == "relative") {
          ## load the admin objects
          ## try two different locations until we standardize
          message('- loading admin objects')
          file_1 <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/%s_%s_admin_draws_raked.Rdata',
                       ind_gp, ind, rd, ind, measure)
          file_2 <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/', 
                           ind, '_raked_admin_draws_eb_bin0_0.RData')
       
          if (file.exists(file_1)) {
            load(file_1)
          } else if (file.exists(file_2)) {
            load(file_2)
          } else {
            stop("Cannot load admin pred object!")
          }

          # Need to load "pseudo cell_pred" admin object to get baseline year
          cell_pred <- get(sprintf('admin_%i', aa))
          cell_pred <- as.data.table(cell_pred)
          
          ## format this object to look like the one we used for cell level (easier to copy code from above)
          ## revised now to work with a bigger array of cell pred objects
          str_match <- stringr::str_match
          draw_cols <- names(cell_pred)[grep("V[0-9]*", names(cell_pred))]
          keep_cols <- c("year", paste0("ADM", aa, "_CODE"), draw_cols)
          cell_pred <- subset(cell_pred, select = keep_cols)
          setnames(cell_pred, paste0("ADM", aa, "_CODE"), "idx")
          cell_pred <- as.matrix(cell_pred)
          
          num_draws <- ncol(cell_pred) - 2
          
          num_idx <- nrow(cell_pred) / length(year_list)
          if (num_idx != length(unique(cell_pred[,2]))) {
            warning(paste0("The number of unique cell_pred indices does not equal the number of ",
                           "rows of cell_pred divided by the number of years.",
                           "\nYou may have duplicate admin codes at the admin ", aa, 
                           " level - check your aggregated objects!"))
          }

          # Assess relative goals
          if (target_type == "greater") {
              relative_proj_draws <- ifelse(proj_draws / baseline_year_draws >= 1 - relative_goal, 1, 0)
            } else if (target_type == "less") {
              relative_proj_draws <- ifelse(proj_draws / baseline_year_draws <= 1 - relative_goal, 1, 0)
            }
          relative_proj_draws[is.na(relative_proj_draws)] <- 0
          relative_goal_prob <- rowMeans(relative_proj_draws)

          # Add spatial index
          relative_goal_prob <- cbind(spatial_idx, relative_goal_prob)

          message("-- Saving probabilities...")
          saveRDS(object = relative_goal_prob, 
                  file = paste0(output_dir, ind, "_", measure, "_", target_year, "_vs_", 
                                baseline_year, "_", abs_rel, "_", target_type, '_', target, 
                                '_adm_', aa, '_target_probs.RDs'))

          # CSV
          relative_goal_prob <- merge_sp_hierarchy(df = relative_goal_prob,
                                                   admin_level = aa,
                                                   idx_col = "spatial_idx",
                                                   sp_h = sp_h)

          write.csv(relative_goal_prob,
                    file = paste0(output_dir, ind, "_", measure, "_", target_year, "_vs_", 
                                baseline_year, "_", abs_rel, "_", target_type, '_', target, 
                                '_adm_', aa, '_target_probs.csv'),
                    row.names = F)

        } # close relative loop
      } # close admin loop
    } # close if_admin section
  } # close goal_obj row loop
}

## compare_all_admin_target #####################################################

#' Runs comparison of projected years against a goal that is meant to be
#' applied to all admins within a given country.
#'
#' This is somewhat of a special use scenario but does have some specific applications.
#' For instance, the GVAP goal for geographic parity in vaccine coverage is "All
#' districts with coverage > 80% within a country by 2020.  The proper way to 
#' calculate the probability of achieving this goal is to calculate by draw 
#' whether the country achieved the goal in the given draw, then
#' see in how many draws the country achieved the goal.  That is what this
#' function does.
#'
#' @param ind_gp indicator group
#' @param ind indicator
#' @param rd run_date
#' @param measure prevalence, incidence, mortality, etc
#' @param target_year Year for which the goal is meant. Can be either a future (projected) year
#'   with the year projections  already  generated for this target_year using \code{make_proj()},
#'   or can be a modeled year (depending on the \code{proj} parameter below.
#' @param target Target (e.g. 0.8 for 80% vaccine coverage); numeric
#' @param target_type Type of target ('greater' or 'less'); character
#' @param admin_level What admin level should the target be evaluated at? For instance, use
#'   \code{1} if the goal applies to all admin1s within a given country, or \code{2} if
#'   all the goal applies to all admin2s within a given country
#' @param uselogit Should this be done in logit space?
#' @param proj Is \code{target_year} projected (\code{proj = T}) or modeled (\code{proj = F})?
#' @return A data table by country with probability of meeting the specified goal
#' @examples
#' # Compare against all-admin2 > 80% target for GVAP
#' compare_all_admin_target(ind_gp = indicator_group, 
#'                          ind = indicator, 
#'                          rd = run_date,
#'                          measure = "prevalence",
#'                          target_year = 2020,
#'                          target = 0.8,
#'                          target_type = "greater",
#'                          admin_level = 2,
#'                          uselogit = T,
#'                          proj = T) 

compare_all_admin_target <- function(ind_gp, 
                                     ind, 
                                     rd,
                                     measure,
                                     target_year,
                                     target, 
                                     target_type,
                                     admin_level,
                                     uselogit,
                                     proj) {

  # Load sp hierarchy
  sp_h <- get_sp_hierarchy()

  # Define directories
  share_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/')
  proj_dir <- paste0(share_dir, "/pred_derivatives/proj/")

  # Convenience
  aa <- admin_level 

  if (proj == T) {

    # Load proj_draws object if working with projected data
    message("-- Loading proj_draws object...")
    proj_draws <- readRDS(sprintf('%s/%s_%s_%s_projections_adm%i_draw_matrix%s.RDs',
                                  proj_dir, ind, measure, target_year, aa,
                                  ifelse(uselogit, "_logit","")))
  } else if (proj == F) {

    # Load admin preds object if working with non-projected data
    # That is, if year is within the modeling frame

    message('Working on ADMIN level')
    ## load the admin objects
    ## try two different locations until we standardize
    file_1 <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/%s_%s_admin_draws_raked.Rdata',
                 ind_gp, ind, rd, ind, measure)
    file_2 <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/', 
                     ind, '_raked_admin_draws_eb_bin0_0.RData')
 
    if (file.exists(file_1)) {
      load(file_1)
    } else if (file.exists(file_2)) {
      load(file_2)
    } else {
      stop("Cannot load admin pred object!")
    }

    # Create a proj_draws object (same name for convenience)
    proj_draws <- get(paste0("admin_", aa))
    proj_draws <- subset(proj_draws, year == target_year)

    # Format to be the same as the other proj_draws object
    proj_draws <- subset(proj_draws, select = names(proj_draws)[!(names(proj_draws) %in% c("year", "pop"))])
    setnames(proj_draws, paste0("ADM", aa, "_CODE"), "spatial_idx")

  }

  # Split off spatial index
  spatial_idx <- proj_draws[,1]
  proj_draws <- proj_draws[, 2:ncol(proj_draws)]

  # Generate admin_level goals
  if (target_type == "greater") {
    absolute_goal_draws <- ifelse(proj_draws >= target, 1, 0)
  } else if (target_type == "less") {
    absolute_goal_draws <- ifelse(proj_draws <= target, 1, 0)
  }

  absolute_goal_draws[is.na(absolute_goal_draws)] <- 0

  # Switch to data table
  absolute_goal_draws <- as.data.table(absolute_goal_draws)

  # Add spatial index back on
  absolute_goal_draws <- cbind(spatial_idx, absolute_goal_draws)

  # Merge the rest of the info on
  absolute_goal_prob <- merge_sp_hierarchy(df = absolute_goal_draws,
                                           admin_level = aa,
                                           idx_col = "spatial_idx",
                                           sp_h = sp_h)
  
  # First split into a list by country, keeping only admin0 levels
  absolute_goal_list <- absolute_goal_prob %>%
    subset(., select = c("ADM0_CODE", "ADM0_NAME", names(.)[grep("V.*", names(.))])) %>%
    split(., by = "ADM0_CODE")

  # Check to see if, within draw, meeting goal for all admins in country
  adm_table_list <- lapply(absolute_goal_list, function(adm_table) {
     
      n <- nrow(adm_table)

      # Grab admin info
      adm_name <- adm_table[,1:2] %>% unique

      # Check if meets criteria in each draw
      adm_table <- adm_table[, 3:ncol(adm_table)] %>%
        summarize_all(sum)

      adm_table <- ifelse(adm_table == n, 1, 0)
      adm_table <- rowMeans(adm_table)
      adm_table <- cbind(adm_name, as.data.table(adm_table))

    })

  adm_table_summary <- do.call(rbind, adm_table_list) %>% 
                        as.data.table %>%
                        setnames(., "adm_table", "probability")

}
