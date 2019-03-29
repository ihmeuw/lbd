# These functions get called in the fractional raking script itsself and as such this script is called from inside the fractional raking script.
# This will definetly get moved to a better location when we roll out fractional raking for everyone.


####### Some function definitons
fetch_from_rdata = function(file_location, item_name, use_grep = F){
  load(file_location)
  if (use_grep) {
    return(mget(grep(item_name, ls(), value = T)))
  }else{
    return(get(item_name))
  }
}


get_link_table_hiv <- function(simple_raster, shapefile_version) {
  #read in link_table
  global_link <- data.table(readRDS(paste0('<<<< FILEPATH REDACTED >>>>')))
  #read in id_raster
  global_raster <- readRDS(paste0('<<<< FILEPATH REDACTED >>>>'))
  #crop id_raster to simple_raster
  cropped_raster <- crop(global_raster, simple_raster)
  #mask id_raster with simple_raster
  masked_raster <- mask(cropped_raster, simple_raster)
  #extract id_raster
  pixel_ids <- raster::extract(masked_raster, extent(masked_raster), na.rm = T)
  pixel_ids <- pixel_ids[!is.na(pixel_ids)]
  #subset link table by pixel ids from extract
  link_table <- global_link[ID %in% pixel_ids,]
  #make sure list of pixel ids is same number of nonna cells in simple raster
  if (length(pixel_ids) != length(which(!is.na(getValues(simple_raster))))) {
    stop("The number of pixel_ids does not match the number of non-na values in simple raster. \nPart of the simple raster may be outside the extent of the global id raster.")
  }
  #return subsetted link table and list of pixel ids
  return(list("link_table" = link_table, "pixel_ids" = pixel_ids))
}


load_populations_cov <- function(reg,
                              pop_measure,
                              measure = 'prevalence',
                              simple_polygon,
                              simple_raster,
                              year_list,
                              interval_mo,
                              outputdir){
  message("Loading the world pop rasters for this region")
  ## Pull 2000-2015 annual population brick using new covariates function
  pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',
                                                       measures = pop_measure,
                                                       simple_polygon = simple_polygon,
                                                       start_year  = min(year_list),
                                                       end_year    = max(year_list),
                                                       interval_mo = as.numeric(interval_mo),
                                                       agebin = 1)

  ## extend and crop pop raster to ensure it matches the simple raster #not convinced this is totally needed
  pop <- pop_raster_annual[[1]]
  pop  <- extend(pop, simple_raster, values = NA)
  pop  <- crop(pop, extent(simple_raster))
  pop  <- setExtent(pop, simple_raster)
  pop  <- raster::mask(pop, simple_raster)

  ## check to ensure the pop raster matches the simple raster in extent and resolution
  if (extent(pop) != extent(simple_raster)) {
    stop("population raster extent does not match simple raster")
  }
  if (any(res(pop) != res(simple_raster))) {
    stop("population raster resolution does not match simple raster")
  }
  #writeRaster(pop, paste0(outputdir, indicator,"_", reg, "_pop_rasters.tif"), format = "GTiff", overwrite = TRUE)

  message("loading the covariate stakers for this model")
  if (measure == 'prevalence') {
    covs = fetch_from_rdata(paste0('<<<< FILEPATH REDACTED >>>>', pathaddin, '.RData'), 'cov_list')
    fes = fetch_from_rdata(paste0('<<<< FILEPATH REDACTED >>>>', pathaddin, '.RData'), 'all_fixed_effects')

    submodels = trimws(strsplit(fes, '+', fixed = T)[[1]])
    covs = covs[submodels]
    #make sure spatial extent is the same
    covs = lapply(covs, function(x) invlogit(crop(x, simple_raster)))
  }else{
    covs = list()
  }

  #bring everything into one place
  covs$pop = crop(pop, simple_raster)
  covnames = names(covs)

  #ensure the dimensions are the same
  for (ccc in covs) {
    stopifnot(dim(ccc)[1:2] == dim(simple_raster)[1:2])
  }

  message("converting the covariate stackers and pop data in to a data table")
  #convert to datables, reshape and stuff
  brick_to_dt = function(bbb){
    dt = setDT(as.data.frame(bbb))
    dt[, pxid := .I] #probably uncessary

    #drop rows now in cellIdx
    dt = dt[pixel_id,]

    dt = melt(dt, id.vars = 'pxid', variable.factor = F)
    dt = dt[,.(value)]
    return(dt)
  }

  covdt = lapply(covs, brick_to_dt)
  covdt = do.call(what = cbind, covdt)
  setnames(covdt, names(covs))
  covdt[,pixel_id := pixel_id]

  #set pop to 0 when pop is na
  covdt[is.na(pop), pop := 0]

  #add year to covdt
  yyy = as.vector(unlist(lapply(min(year_list):max(year_list), function(x) rep.int(x, times = length(pixel_id)))))
  covdt[,year := yyy]

  #free up a bit of space
  rm(covs)

  return(covdt)
}


get_gbd_locs <- function(rake_subnational = T,
                         reg){
  library(sf)
  if (rake_subnational == T) {
    connector <-
      st_read(get_admin_shapefile(admin_level = 0, raking = T), quiet = T) %>%
      st_set_geometry(NULL) %>%
      mutate(ADM0_CODE = as.numeric(as.character(ADM0_CODE))) %>%
      mutate(ADM1_CODE = as.numeric(as.character(ADM1_CODE))) %>%
      filter(ADM0_CODE %in% get_gaul_codes(reg)) %>%
      dplyr::select(ADM0_CODE, ADM1_CODE, loc_id, rak_level) %>%
      mutate(rak_level = as.character(rak_level)) %>%
      dplyr::rename(location_id = loc_id) %>%
      mutate(location_id = as.numeric(as.character(location_id)))
  } else {
    connector <-
      get_location_code_mapping()[GAUL_CODE %in% get_gaul_codes(reg), list(location_id = loc_id, GAUL_CODE)]
  }

}



sub_nat_link_merge <- function(rake_subnational,
                               link,
                               connector){
  if (rake_subnational == T) {
    connector0 <- connector[which(connector$rak_level == 0), ]
    connector0$ADM1_CODE <- NULL
    connector0$loc_type <- NULL
    connector1 <- connector[which(connector$rak_level == 1), ]
    connector1$ADM0_CODE <- NULL
    connector1$loc_type <- NULL
    link0 <- merge(link, connector0, by = c("ADM0_CODE"))
    link1 <- merge(link, connector1, by = c("ADM1_CODE"))
    link <- rbind(link0, link1)
  } else {
    link <- merge(link, connector, by.x = c("ADM0_CODE"), by.y = c("GAUL_CODE"))
  }
  return(link)
}



prep_link_table <- function(link_table = link_table,
                            simple_raster = simple_raster,
                            pixel_id = pixel_id){
  ####################################################################

  #keep only locations in the region
  link <- link_table[[1]]
  link = link[ADM0_CODE %in% unique(simple_raster),]# This means only cells and cell fractions that are in the region are included in the link table

  #keep only cells where simple raster identifies it as belonging to the country or region
  link_map = data.table(af_id = link_table[[2]], reg_id = pixel_id, simp_ras_val = simple_raster[pixel_id])#creates table relating pixel number in the africa simple raster (id in link table) to pixel number in simple raster (id in simple raster/cell pred), to the value in the simple raster(ADM0_CODE).
  link = merge(link, link_map, by.x = 'ID', by.y = 'af_id')#merges that relational table to the link table

  #scale up such that the sum of area fractions for each pixel is equal to 1 so pieces are not lost
  link[, total_frac := sum(area_fraction), by = "ID"]
  link[, c('area_fraction','old_frac') := list(area_fraction/total_frac, area_fraction)]

  return(link)
}



prep_cell_pred <- function(cell_pred=cell_pred,
                           cell_ids = cell_ids,
                           pixel_id = pixel_id,
                           covdt = covdt){
  #set cell pred as a data table, and rename things
  if (!inherits(x = cell_pred, 'data.table')) {
    cell_pred = as.data.table(cell_pred)
    names(cell_pred) = paste0('V',1:ncol(cell_pred))
  }

  cell_pred[, cell_pred_id := .I] #cell_pred object ID
  cell_pred[,cell_id := cell_ids] #cell id references the africa map
  cell_pred[,pixel_id := pixel_id] #pixel id references the regional map

  #add population, year and potentially the stackers
  cell_pred = cbind(cell_pred, covdt[,c('year', 'pop'),with = F])

  #make sure it behaved
  stopifnot(any(!(cell_pred[,pixel_id] != rep.int(covdt[,pixel_id], 18))))

  return(cell_pred)
}


calculate_pop_scalars <- function(cell_pred = cell_pred,
                                  age_group = age_group,
                                  connector = connector,
                                  sex = sex_id,
                                  sharedir = sharedir,
                                  run_date = run_date,
                                  indicator = indicator,
                                  stratum = stratum){

  #do calculations!
  rake_geo_pop = cell_pred[, lapply(c('pop'), function(x) sum(get(x), na.rm = T)), by = c('year','location_id')]
  rake_geo_pop$world_pop_unraked <- rake_geo_pop$V1
  loc_ids = unique(connector$location_id)

  #adjust to be GBD populations
  source(paste0('<<<< FILEPATH REDACTED >>>>',"/get_population.R"))

  gbd_pops = get_population(age_group_id = age_group ,location_id = loc_ids, sex = sex_id, year_id = unique(rake_geo_pop[,year]), gbd_round_id = 5)
  scalars = merge(rake_geo_pop, gbd_pops, by.x = c('location_id','year'), by.y = c('location_id','year_id'))
  scalars[,pop_scalar := population/world_pop_unraked]

  # for records
  scalars$gbd_pop <- scalars$population
  scalars$population <- NULL

  #trim Scalars object
  scalars$age_group_id <- NULL
  scalars$sex_id <- NULL
  scalars$run_id <- NULL
  scalars$V1 <- NULL

  write.csv(scalars, file = paste0('<<<< FILEPATH REDACTED >>>>', indicator, "_", stratum, "_pop_rf.csv"))

  return(scalars)
}



calculate_fractional_rfs <- function(ndraws=ndraws,
                                     cell_pred = cell_pred,
                                     gbd=gbd,
                                     sharedir = sharedir,
                                     run_date = run_date,
                                     indicator = indicator,
                                     stratum = stratum){
  message("converting from pervalence to counts")
  #set the variables to aggregate
  overs = paste0('V',1:ndraws)

  #covert to counts
  cell_pred = cell_pred[, (overs) := lapply(overs, function(x) get(x) * pop_raked )]

  #do calculations!
  rake_geo = cell_pred[, lapply(c(overs,'pop_raked'), function(x) sum(get(x), na.rm = T)), by = c('year', 'location_id')]
  setnames(rake_geo, grep('V[0-9]', names(rake_geo),value = T),c(overs, 'pop_raked'))

  #convert back to rates/prevalence
  rake_geo = rake_geo[, (overs) := lapply(overs, function(x) get(x)/pop_raked) ]

  # merge to admin estimates
  rake_geo <- merge(rake_geo, gbd, by.x = c('location_id', 'year'), by.y = c('name', 'year'))

  # finding the mean of the draws at the raking geographies
  rake_geo$gbd_hiv_prev <- rake_geo$mean
  if (ndraws == 100) {
    rake_geo$mbg_rate <-  rowMeans(rake_geo[ ,V1:V100])
  } else if (ndraws == 1000) {
    rake_geo$mbg_rate <-  rowMeans(rake_geo[ ,V1:V1000])
  } else if (ndraws == 15) {
    rake_geo$mbg_rate <-  rowMeans(rake_geo[ ,V1:V15])
  } else {
    message("the number of draws is not 15, 100, or 1000, critical failure.  see raking code")
    stop()
  }
  rake_geo$rf <- rake_geo$gbd_hiv_prev/rake_geo$mbg_rate

  # Clear Out
  message("creating fractional raking factors table")
  fractional_rf <- rake_geo
  fractional_rf$mbg_hiv_prev <- fractional_rf$mbg_rate
  fractional_rf <- fractional_rf[ , c("location_id", "year", "mbg_hiv_prev", "gbd_hiv_prev", "rf")]

  # saving the raking factors
  write.csv(fractional_rf, file = paste0('<<<< FILEPATH REDACTED >>>>', indicator, "_", stratum, "_rf.csv"))

  return(fractional_rf)
}










