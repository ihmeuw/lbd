# Link table generation (parallel) functions

# file name of final link table/id raster
LT_FILENAME <- 'lbd_standard_link.rds'
IR_FILENAME <- 'lbd_standard_id_raster.rds'

#' @title Partition Shapefile for Parallel Link Table Generation
#' 
#' @description partition the dataset by ADM0_CODE based on area
#' 
#' @param N (integer) number of pieces to break job into
#' @param tmpdir (string) temporary directory to save files
#' @param shapedir (string) shapefile directory
#' @param region (string) ISO region(s)/stages for subsetting shapefile

partition_link_table <- function(N, tmpdir, shapedir, region = 'stage1+stage2') {
  # column with which to group data on, will be divided into N groups
  group_col <- 'ADM0_CODE'
  
  # admin 2 shapefile
  polys <- st_read(get_admin_shapefile(admin_level = 2, version = shapedir), stringsAsFactors = FALSE)
  
  # subset shapefile codes
  s1 <- get_adm0_codes(region)
  polys <- polys[polys$ADM0_CODE %in% c(s1),]
  
  # get area of polygon geometry
  polys$area <- st_area(polys$geometry)
  
  ## greedy partitioning algorithm - partition adm0 codes based on area of polygons
  ## input should be of type "sf"
  greedy_partition <- function(sfdf, by, on){
    # get area for each adm0 code 
    areas <- aggregate(sfdf[[on]], by = list(sfdf[[by]]), FUN = sum)
    areas <- areas[with(areas, order(-areas[, 2])),]  # decreasing sort
  
    # keep track of totals for each partition
    totals <- double(length = N)
    
    # set up vector to assign task ids
    tasks <- vector('list', length = N)
    
    # assign each adm0 code to the current smallest (by area) partition
    for (i in 1:nrow(areas)) {
      row <- as.numeric(areas[i,])
      code <- row[1]
      area <- row[2]  
      smallest <- which.min(totals)  # get current smallest bucket
      totals[smallest] <- totals[smallest] + area  # add area to bucket total
      tasks[[smallest]] <- c(tasks[[smallest]], code)  # add adm0 to bucket
    }
    
    return(tasks)
  }
  
  tasks <- greedy_partition(polys, group_col, 'area')
  saveRDS(tasks, file.path(tmpdir, 'task_lookup.rds'))
  message(sprintf("File saved into '%s'.", tmpdir))
}


#' @title Build Sub Link Tables
#' 
#' @description Build link table for specific regions
#'
#' @param tmpdir (string) temporary directory to open/save files
#' @param shapedir (string) shapefile directory 
#' @param task_id (integer) task id gathered from the environment for array job
#' @param cores (integer) number of cores for local parallelization while computing intersections

build_sub_link_tables <- function(tmpdir, shapedir, task_id, cores) {
  # transform adm0 codes -> ISO3 (iso+iso notation) for build_link_table function
  tasks <- readRDS(file.path(tmpdir, 'task_lookup.rds'))
  adm0_codes <- tasks[[task_id]]
  lookup_table <- load_adm0_lookup_table()[, c('gadm_geoid', 'iso3')] 
  iso_codes <- lookup_table[lookup_table$gadm_geoid %in% adm0_codes]$iso3
  region <- paste(iso_codes, collapse = '+')
  
  # build link table for specified regions
  message(sprintf("Dispatching job #%i for region(s): %s.", task_id, region))
  link_table <- build_link_table(shapedir, cores, region, sample_ids = TRUE)
  
  message("Saving temporary files.")
  
  save_file <- sprintf('%s/link_table_%s.rds', tmpdir, task_id)
  saveRDS(link_table$link_table, save_file)
  
  save_file <- sprintf('%s/id_raster_%s.rds', tmpdir, task_id)
  saveRDS(link_table$id_raster, save_file)
}


#' @title Combine Sub Link Tables
#' 
#' @description combine all sub tables from temporary directory and save final link table
#' 
#' @param N (integer) number of pieces to break job into
#' @param tmpdir (string) temporary directory where link tables are stored
#' @param shapedir (string) shapefile directory
#' @param cores (integer) number of cores for local parallelization when loading sub tables

combine_link_tables <- function(N, tmpdir, shapedir, cores) {
  # read in partitioned tables
  registerDoParallel(cores = cores)
  link_table <- foreach(i = 1:N, .final = rbindlist) %dopar%
    readRDS(sprintf('%s/link_table_%i.rds', tmpdir, i))
  
  # fix area fractions
  link_table[, area_fraction := end_area / start_area]
  
  # fix grouping information (total area and n)
  link_table[, n := .N, by = ID]
  link_table[, total_area := sum(area_fraction), by = ID]
  link_table$area_fraction <- as.numeric(link_table$area_fraction)
  link_table$total_area <- as.numeric(link_table$total_area)
  
  ## fix area fractions w/ water
  # if there is only one admin area assigned and area fractions is less than 1, set area fraction to 1
  link_table[n == 1 & area_fraction < 1, area_fraction := 1]
  
  # recompute area for pixels in multiple admin units and water proportionally
  link_table[n > 1 & total_area < 1, area_fraction := area_fraction / total_area]
  
  save_file <- file.path(get_admin_shape_dir(version = shapedir), LT_FILENAME)
  message(sprintf("Saving link table to %s", save_file))
  saveRDS(link_table, save_file)
  Sys.chmod(save_file, '0744')

  # id raster
  id_raster <- empty_world_raster()
  id_raster[] <- 1:length(id_raster)
  
  save_file <- file.path(get_admin_shape_dir(version = shapedir), IR_FILENAME)
  message(sprintf("Saving id_raster to %s", save_file))
  saveRDS(id_raster, save_file)
  Sys.chmod(save_file, '0744')
}


#' @title Test and Clean Up Link Table After Build
#' 
#' @description check if link table was created correctly
#'
#' @param N (integer) number of nodes
#' @param tmpdir (string) temporary directory where link tables are stored
#' @param shapedir (string) shapefile directory

cleanup_link_tables <- function(N, tmpdir, shapedir) {
  # check if partition job failed
  if (!file.exists(file.path(tmpdir, 'task_lookup.rds'))) {
    msg <- 'Link Table job FAILED.\n\nCould not partition properly.'
  } else {
    # check if all sub link tables were created properly
    missing <- c()
    for (i in 1:N) {
      if (!file.exists(sprintf('%s/link_table_%i.rds', tmpdir, i))) {
        missing <- c(missing, i)
      }
    }
    
    # some link tables are missing
    if (length(missing) > 0) {
      # get adm0 codes that are missing
      task_list <- readRDS(file.path(tmpdir, 'task_lookup.rds'))
      missing_codes <- unlist(sapply(missing, function(i) {
        task_list[[i]]
      }))
      
      msg <- sprintf('Link table job FAILED. \n\nThe following admin codes are missing (any or all of them could be the problem): %s.',
                     paste(sort(as.integer(missing_codes)), collapse = ', '),
                     tmpdir)
    } else if (!file.exists(file.path(get_admin_shape_dir(version = shapedir), LT_FILENAME))) {
      msg <- 'Link table job FAILED. Link table is missing (most likely directory or permissions issues)'
    } else {
      msg <- sprintf('Link table job SUCCESS. \n\nNew table saved to: %s/%s.', shapedir, LT_FILENAME)
    }
  }
  
  # output result
  print(msg)
}