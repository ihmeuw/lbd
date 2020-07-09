# -------------------------------------------------------------------------------------------------------------
# format cell pred function
format_cell_pred <- function(ind_gp,
                             ind,
                             rd,
                             reg,
                             measure,
                             pop_measure,
                             year_start,
                             year_end,
                             var_names = ind, # name using ind by default, but can pass custom name
                             matrix_pred_name = NULL,
                             skip_cols = NULL,
                             rk = T,
                             rk_measure = 'prevalence',
                             coastal_fix = T, # if model wasn't run w new coastal rasterization, force old simple raster process 
                             rake_subnational = rk,
                             shapefile_version = 'current') {
  
  
  message('loading simple raster & populations')
  
  ## Load simple polygon template to model over
  gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)
  simple_polygon_list <- load_simple_polygon(
    gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
    shapefile_version = shapefile_version
  )
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  
  ## Load list of raster inputs (pop and simple)
  if (coastal_fix) raster_list <- build_simple_raster_pop(subset_shape, link_table = shapefile_version) #uses new simple_raster 
  else raster_list <- build_simple_raster_pop(subset_shape, link_table = NULL) #uses old rasterize
  
  simple_raster <- raster_list[["simple_raster"]]
  pop_raster <- raster_list[["pop_raster"]]
  pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  
  ## get number of years
  year_list <- c(year_start:year_end)
  num_yrs <- length(year_list)
  
  message('loading links')
  #####################################################################
  # load the cell id to admin units link
  link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)
  
  #####################################################################
  # collect and load the population data from the WorldPop rasters
  covdt <- load_populations_cov(reg, pop_measure=pop_measure, measure = measure, simple_polygon, 
                                simple_raster, year_list, interval_mo=12, pixel_id = pixel_id)
  #####################################################################
  # Prepping the cell_pred and link table to be linked by making sure they have the appropriate identifiers.  Also performs a
  # zippering at the region boundary where cells that have a portion of their area outside of the modeling region are reintegrated
  # as a whole cell and considered to be only in one region.  This works becasue a given cell is only modeled in one region.
  link <- prep_link_table(
    link_table = link_table,
    simple_raster = simple_raster,
    pixel_id = pixel_id
  )
  
  # getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
  # gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary
  connector <- get_gbd_locs(
    rake_subnational, 
    reg = reg,
    shapefile_version = shapefile_version
  )
  
  # merge the connector on to the link table by making sure that each cell fragment gets connected to the appropriate
  # raking geography
  link <- sub_nat_link_merge(
    rake_subnational,
    link,
    connector
  )
  
  #helper function to load a list of cell preds and then merge them together
  loadCellPreds <- function(i) { 
    
    message('~>loading cell pred for: ', ind[[i]])
    
    # Load the relevant pred object - loads an object named cell_pred
    rdata_file <- paste0('<<<< FILEPATH REDACTED >>>>',
                         ind[[i]], 
                         ifelse(rk, paste0("_cell_draws_raked_", rk_measure, "_eb_bin0_"), "_cell_draws_eb_bin0_"), 
                         reg, "_0.RData")
    
    #TODO improve logic
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
    
    # Verify alignment  
    if(nrow(cell_pred)!=nrow(covdt)) stop('Dimensional mismatch between cell pred and simple raster!!')
    
    # set cell pred as a data table, and rename things
    cell_pred <- prep_cell_pred(
      cell_pred = cell_pred,
      cell_ids = link_table[[2]],
      pixel_id = pixel_id,
      covdt = covdt
    )
    
    message('~~>reshape long and formatting')
    #reshape long draws and key on the pixel ID/draw/year
    dt <- melt(cell_pred,
               measure = patterns("V"),
               variable.name = "draw",
               value.name = var_names[[i]]) %>% 
      #toconvert draw col to int instead of V1-250 as a factor
      setkey(., pixel_id, draw, year, cell_pred_id, cell_id, pop) %>%
      return
    
  }
  
  #load/format all the cell preds and then merge them together
  cell_pred <- lapply(1:length(ind), loadCellPreds) %>% 
    Reduce(function(...) merge(..., all = TRUE), .)
  
  # merge cell_pred on the link
  cell_pred <- merge(link[, -c('pixel_id')], cell_pred, by.x = "ID", by.y = "cell_id", allow.cartesian = TRUE)
  
  # space
  link <- NULL
  
  #subset to relevant columns and return
  keep_vars <- c('ADM0_CODE', 'ADM1_CODE', 'ADM2_CODE', 
                 'pixel_id', 'year', 'pop', 'area_fraction', 'draw', unlist(var_names))
  cell_pred[, (keep_vars), with=F] %>% 
    return
  
}
# -------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------
# format admins function
format_admin_results <- function(ind_gp,
                                 ind,
                                 rd,
                                 measure,
                                 suffix,
                                 var_names = ind, # name using ind by default, but can pass custom name
                                 rk) {
  
  #helper function to load a list of admin results and format them into a single DT
  loadAdmins <- function(i) { 
    
    message('~>loading results for: ', ind[[i]])
    
    # def directory
    dir <- '<<<< FILEPATH REDACTED >>>>'
    
    # read in all the admin objects from a single RData file
    paste0(dir, '/', ind[[i]], ifelse(rk[[i]], '_raked', '_unraked'), measure[[i]], '_admin_draws', suffix[[i]], '.RData') %>% 
      load(envir = globalenv(), verbose=T)
    
    # harmonize the hierarchy lists (some are using factor variables which don't merge well)
    factor_vars <- names(sp_hierarchy_list)[vapply(sp_hierarchy_list, is.factor, c(is.factor=FALSE))]
    if (length(factor_vars) > 0) {
      
      message('formatting spatial hierarchy table to num/chr - input table uses factor variables')
      
      #build helper functions
      formatHierarchy <- function(dt) {
        
        out <- copy(dt)
        
        #helper fx to convert factors without information loss
        #https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
        facToNum <- function(f) as.numeric(as.character(f))
        
        #we want to convert the codes to num and the names to chr
        cols_to_num <- factor_vars[factor_vars %like% 'CODE']
        cols_to_chr <- factor_vars[factor_vars %like% 'NAME']
        cols <- list(cols_to_num, cols_to_chr)
        funs <- rep(c(facToNum, as.character), lengths(cols))
        
        # convert class based on column name
        out <- out[, unlist(cols) := Map(function(f, x) f(x), funs, .SD), .SDcols = unlist(cols)] %>% 
          return
        
      }
      
      #convert vars
      sp_hierarchy_list <- formatHierarchy(sp_hierarchy_list)
      
      #reassess and test
      if(vapply(sp_hierarchy_list, is.factor, c(is.factor=FALSE)) %>% any) stop('Failed to convert sp hierarchy!')
      
    }
    
    # format and append
    bindResults <- function(input_dt, info) {
      
      dt <- copy(input_dt)
      
      #pull out the level of aggregation, rename the variable to harmonize and record the value for later
      lvl_str <- names(dt)[names(dt) %like% 'CODE']
      lvl <- substr(lvl_str, start=1, stop=4) #extract level value
      setnames(dt, lvl_str, 'code') # rename
      dt[, agg_level := lvl] #record the level
      dt[, c('pop', 'region') := NULL] #remove unecessary vars
      
      #merge on location names while simultaneously formatting them
      out <- merge(dt,
                   info[, .(code=get(lvl_str), 
                            name=paste0(lvl, '_NAME') %>% get)] %>% unique, 
                   by='code',
                   all.x=T)
      
      #melt and output
      melt(out,
           measure = patterns("V"),
           variable.name = "draw",
           value.name = var_names[[i]]) %>% 
        return
      
    }
    
    dt <- list(admin_0, admin_1, admin_2) %>% 
      lapply(., bindResults, info=sp_hierarchy_list) %>% 
      rbindlist %>% 
      return
    
  }
  
  #load/format all the admin results and then merge them together
  dt <- lapply(1:length(ind), loadAdmins) %>% 
    Reduce(function(...) merge(..., all = TRUE), .) %>% 
    return
  
}
# -------------------------------------------------------------------------------------------------------------
