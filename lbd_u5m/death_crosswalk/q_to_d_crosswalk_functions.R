#' Run q to raked d crosswalk functions
#' 
#' Description: Takes a cell pred or raster object with mortality probability (q) and
#' converts it to deaths. Fractionally rakes deaths to gbd
#'
#' @param run_date model run date
#' @param region character string, modelling region used
#' @param ages character string, only `under5` currently tested
#' @param shapefile_version string, default `2018_08_01`, which is the gadm shapefile. `current` or other shapefile dates in "<<<< FILEPATH REDACTED >>>>" also acceptable. 
#' @param simple_raster simple raster loaded from load_simple_polygon and build_simple_raster_pop. Optional to save time
#' @param pop_raster_brick pop_raster brick loaded from load_and_crop_covariates_annual. Optional to save time. 
#' 
#' @return saves death draws, raking factors, adm0, 1, and 2 aggregates, and mean, upper and lower rasters into `run_date` folders for `ages`.
#'
run_q_to_d_crosswalk <- function(run_date,
                                 region,
                                 ages = c("under5", "infant", "neonatal"),
                                 shapefile_version = "2018_08_01",
                                 simple_raster = NULL,
                                 pop_raster_brick = NULL){
  
  #pull config file from model run in separate environment to avoid namespace pollution
  e1 <- new.env()
  load("<<<< FILEPATH REDACTED >>>>", e1)
  config <- get('config', e1)
  rm(e1)
  
  #pull pop_measure, year_list, and interval_mo from config file from model run
  pop_measure <- config[V1 == "pop_measure", V2]
  year_list <- config[V1 == "year_list", V2]
  year_list <- eval(parse(text = year_list))
  interval_mo <- as.numeric(config[V1 == "interval_mo", V2])
  
  if(is.null(simple_raster)) {
    #load simple raster
    #shapefile version here is different from the shapefile version passed to run_q_to_d_crosswalk
    adm0_list      <- get_adm0_codes(region, shapefile_version = shapefile_version)
    simple_polygon <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, tolerance = 0.4,
                                               shapefile_version = shapefile_version)
    subset_shape   <- simple_polygon[['subset_shape']]
    simple_polygon <- simple_polygon[['spoly_spdf']]
    
    message('Loading simple raster')
    raster_list    <- build_simple_raster_pop(subset_shape) #,u5m=TRUE)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]
  }
  
  if(is.null(pop_raster_brick)){
    #get world_pop raster
    pop <- load_and_crop_covariates_annual(covs = 'worldpop',
                                           measures = pop_measure, # Defined above
                                           simple_polygon = simple_polygon,
                                           start_year  = min(year_list),
                                           end_year    = max(year_list),
                                           interval_mo = interval_mo,
                                           agebin=1)[[1]]
    
    pop <- extend(pop, simple_raster, values=NA)
    pop <- crop(pop, extent(simple_raster))
    pop <- setExtent(pop, simple_raster)
    pop <- mask(pop, simple_raster)
  } else {
    pop <- pop_raster_brick
  }
    
  for(age in ages) {
    message("~~~~~~~~~~~~~ Starting crosswalk for ", age, " ~~~~~~~~~~~~~")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # prepare GBD values for raking
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #get gbd values for q to m to d crosswalks(location hierarchy, life table, population, deaths)
    gbd <- get_gbd_crosswalk_values(age, year_list, shapefile_version = shapefile_version)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Crosswalk
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #load cell_pred in separate environment to avoid namespace pollution
    e1 <- new.env()
    cell_draws <- readRDS("<<<< FILEPATH REDACTED >>>>")
    cell_draws <- get('cell_draws', e1)
    rm(e1)
    
    #get draw column names
    ndraws <- ncol(cell_draws)
    overs <- paste0('V',1:ndraws)
    
    #if column names in draws are repeated, reset column names based on column number
    if(length(unique(colnames(cell_draws))) != length(colnames(cell_draws))){
      colnames(cell_draws) <- overs
    }
    
    #crosswalk raster from q to m
    cell_pred <- crosswalk_qm(cell_draws, "cell_pred", region, year_list, age, gbd, simple_raster, shapefile_version)
    rm(cell_draws)
    gc(full=TRUE)
    
    #go from m to d by multiplying by population
    message("calculating raw deaths from m")
    unraked_deaths_cell_pred <- calculate_deaths_cell_pred(cell_pred, simple_raster, year_list, pop, overs,age)
    rm(cell_pred)
    gc(full=TRUE)
    
    setnames(gbd, c("location_id", "year_id", "deaths"), c("name", "year", "mean"))
    gbd <- gbd[,c("name", "year", "mean")]
    count_cell_pred <- unraked_deaths_cell_pred[,overs, with = F]
    
    output <- fractionally_rake_counts(count_cell_pred = count_cell_pred,
                                       rake_to = gbd,
                                       reg = reg,
                                       year_list = year_list,
                                       rake_subnational = F,
                                       countries_not_to_subnat_rake = NULL,
                                       countries_not_to_rake = c("GUF+ESH"),
                                       simple_raster = simple_raster,
                                       modeling_shapefile_version = modeling_shapefile_version,
                                       raking_shapefile_version = raking_shapefile_version)
    
    
    #make mean, upper and lower raster bricks
    raster_bricks_list <- list()
    for(func in c("mean", "upper", "lower")) {
      raster_bricks_list[[func]] <- cell_pred_to_raster_brick(output[["raked_cell_pred"]], simple_raster, year_list, func)
    }
    
    #compile outputs
    output <- c(output, raster_bricks_list)
    
    #save outputs
    message("saving crosswalk outputs \n")
    save_death_crosswalk_outputs(
      output_list=output, age=age, run_date=run_date, region=region
    )
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# q to m crosswalk functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Crosswalk Cell Pred or Raster from q to m
#' 
#' Description: Takes a cell pred or raster object with mortality probability (q) and
#' converts it to mortality rate (m) using a from gbd life tables. Uses formula 
#' m = 1 / (-n + a  + n/q)
#'
#' @param data a cell pred or raster object
#' @param data_type character string, `raster` or `cell_pred` acceptable
#' @param region character string, modelling region used
#' @param year_list list of years in `data`
#' @param age character string, `under5`, `4q1`, `infant`, `neonatal` acceptable
#' @param gbd a table from `get_gbd_crosswalk_values()` with a values for Adm0 gaul codes
#' @param simple_raster default `NULL`, can pass in a simple raster object to speed things up,
#' otherwise pulls simple raster using `region`.`
#' @param shapefile_version string, date of shapefile. `current` also acceptable
#' 
#' @return a cell pred or raster object of type `data_type`
#'
crosswalk_qm <- function(data, data_type, region, year_list, age, gbd, simple_raster = NULL, shapefile_version) {
  
  #set n depending on age for q to m formula
  if(age == "under5") {
    n <- 5
  } else if (age == "4q1") {
    n <- 4
  } else if(age == "infant") {
    n <- 1
  } else if(age == "neonatal") {
    n <- 1/12
  }
  
  #skip if simple raster is passed in to function
  if(is.null(simple_raster)) {
    ## get simple polygon and simple raster used to produce cell pred
    message('Loading simple polygon')
    adm0_list           <- get_adm0_codes(reg, shapefile_version = shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, tolerance = 0.4,
                                               shapefile_version = shapefile_version)
    subset_shape   <- simple_polygon[['subset_shape']]
    simple_polygon <- simple_polygon[['spoly_spdf']]
    
    message('Loading simple raster')
    raster_list    <- build_simple_raster_pop(subset_shape) #,u5m=TRUE)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]
  }
  
  if(!(data_type %in% c("raster", "cell_pred"))) {
    stop("data_type must be either raster or cell_pred")
  }
  
  if(data_type == "raster") {
    
    ras_list <- unstack(data)
    
    if(!compareRaster(ras_list[[1]], simple_raster, stopiffalse = F)) {
      stop("There is a mismatch in the raster extent, dimensions, projection, or resolution")
    }
    
    message("crosswalking raster object from q to m")
    # unstack raster into list and extract
    ras <- lapply(ras_list, function(i){raster::extract(i, extent(i))})
    
    # extract simple raster to list
    simple <- raster::extract(simple_raster, extent(simple_raster))
    
    # combine raster layers with simple raster into table
    dt <- lapply(ras, function(i){data.table(i, simple)})
    dt <- lapply(dt, function(i){setnames(i, "i", "ras")})
    
    #get years in there somehow
    for(i in 1:length(year_list)) {
      dt[[i]][, year := year_list[[i]]]
      dt[[i]][, id := .I]
    }
    
    # merge on gbd to get a
    fi <- lapply(dt, function(i){merge(i, gbd, by.x = c("simple", "year"), by.y = c("gaul", "year_id"), all.x = T)
    })
    
    fi <- lapply(fi, function(i){
      i[order(id),]
    })
    
    # calculate m using function m = 1 / (-n + a  + n/q), or m = -12 * log(1 - q) for neonatal
    fin <- lapply(fi, function(i){
      i[,nmx := -1.]
      if(age == "neonatal"){
        i[!is.na(ras), nmx := as.double((-12 * log(1-ras)))]
      } else {
        i[!is.na(ras), nmx := as.double((1 / (-n + a + (n/ras))))]
      }
      i[nmx == -1, nmx := NA]
    })
    
    # assign new m values to rasters and rebrick
    new_ras <- ras_list
    for(i in 1:length(ras_list)) {
      values(new_ras[[i]]) <- 0
      new_ras[[i]] <- insertRaster(new_ras[[i]], cbind(fin[[i]]$nmx))
    }
    
    output <- brick(new_ras)
    return(output)
    
  } else {
    if(age == "neonatal") {
      m <- -12 * log(1 - data)
    } else {
      # extract simple raster and remove NAs
      simple <- raster::extract(simple_raster, extent(simple_raster))
      simple <- simple[!is.na(simple)]
      
      # calculate length of one year of cell pred object
      cell_pred_len <- dim(data)[1] / length(year_list)
      
      # check that the simple_raster pairs with the cell pred object
      if(length(simple) != cell_pred_len) {
        stop("the simple_raster does not have the same number of non-NA cells as the cell pred")
      }
      
      message("crosswalking cell_pred object from q to m")
      # duplicate years and gaul codes to match with length of cell pred and combine into data.table
      years <- rep(year_list, each=cell_pred_len)
      simple_list <- rep(simple, times = length(year_list))
      combined <- data.table(years, simple_list)
      combined[,id := .I]
      
      # merge gbd info on year and gaul_code
      merged <- merge(combined, gbd, by.x = c("years", "simple_list"), by.y = c("year_id", "gaul"), all.x=T, all.y = F)
      
      merged <- merged[order(id)]
      
      #fix for french guiana and western sahara with no gbd a values
      merged[is.na(a), a := 1]
      
      # make matrix of a values for crosswalk calculation
      a_matrix <- matrix(nrow = nrow(data), ncol = ncol(data))
      for(i in 1:ncol(a_matrix)) {
        a_matrix[,i] <- merged$a
      }
      
      # create new cell pred object with m
      # calculate m using function m = 1 / (-n + a  + n/q)
      m <- 1 / (-n + a_matrix + (n / data))
    }
    return(m)
  }
}


#' Pull nax from gbd life tables
#' 
#' a is not necessary for the neonatal calculation, so this function is only needed for under5 and infant
#' 
#' @param age character string, `under5`, `infant` acceptable
#' @param year_list list of years in `data`
#' @param shapefile_version string, date of shapefile. `current` also acceptable
#' 
#' @return a data.table with columns "age_group_id", "year_id", "location_id", "a", "gaul",
#'
get_gbd_a <- function(age,
                      year_list,
                      shapefile_version) {
  
  message("Pulling GBD nax from get_life_table")
  
  source("<<<< FILEPATH REDACTED >>>>")

  if(age == "under5"){
    age_id <- c("5", "28")
  } else if(age == "infant") {
    age_id <- c("28")
  }
  
  locs <- get_location_code_mapping(shapefile_version = shapefile_version) 
  loc_ids <- locs$loc_id
  
  lt <- data.table(get_life_table(age_group_id = age_id, 
                                  year_id = year_list, 
                                  sex_id = 3, 
                                  location_id = loc_ids, 
                                  life_table_parameter_id=c(2), 
                                  gbd_round_id = 5))
  
  combined <- merge(lt, locs, by.x = "location_id", by.y = "loc_id", all.x = T)
  setnames(combined, c("mean", "GAUL_CODE"), c("a", "gaul"))
  combined <- combined[,c("age_group_id", "year_id", "location_id", "a", "gaul")]
  
  combined[is.na(gaul), gaul := -1]
  
  return(combined)
}


#' Pull ndx from gbd get envelope
#' 
#' specified to pull death with shocks and hiv
#' 
#' @param age character string, `under5`, `infant`, `neonatal` acceptable
#' @param year_list list of years in `data`
#' @param shapefile_version string, date of shapefile. `current` also acceptable
#' 
#' @return a data.table with columns "age_group_id", "year_id", "location_id", "deaths", "gaul",
#'
get_gbd_d <- function(age,
                      year_list,
                      shapefile_version) {
  
  message("Pulling GBD ndx from get_envelope")
  
  source("<<<< FILEPATH REDACTED >>>>")
  
  if(age == "under5"){
    age_id <- c("5", "28")
  } else if(age == "infant") {
    age_id <- c("28")
  } else if(age == "neonatal"){
    age_id <- c("2", "3")
  }
  
  locs <- get_location_code_mapping(shapefile_version = shapefile_version) 
  loc_ids <- locs$loc_id
  
  lt <- data.table(get_envelope(age_group_id = age_id, 
                                year_id = year_list, 
                                sex_id = 3, 
                                location_id = loc_ids, 
                                gbd_round_id = 5,
                                with_shock = 1,
                                with_hiv = 1))
  
  combined <- merge(lt, locs, by.x = "location_id", by.y = "loc_id", all.x = T)
  setnames(combined, c("mean", "GAUL_CODE"), c("deaths", "gaul"))
  combined <- combined[,c("age_group_id", "year_id", "location_id", "deaths", "gaul")]
  
  combined[is.na(gaul), gaul := -1]
  
  if(age == "neonatal"){
    combined <- combined[, .(deaths = sum(deaths, na.rm=T)), by = c("year_id", "location_id", "gaul")]
    combined[, age_group_id := 42]
  }
  
  return(combined)
}  


#' pulls nax and ndx for q to d crosswalk
#' 
#' @param age character string, `under5`, `infant`, `neonatal` acceptable
#' @param year_list list of years in `data`
#' @param shapefile_version string, date of shapefile. `current` also acceptable
#' 
#' @return a data.table with columns "age_group_id", "year_id", "location_id", "deaths", "gaul", "a", and "deaths"
#'
get_gbd_crosswalk_values <- function(age, 
                                     year_list,
                                     shapefile_version){
  
  if(!(age %in% c("under5", "infant", "neonatal"))) {
    stop("age must be under5, infant, neonatal")
  }
  
  if(age == "under5") {
    a <- get_gbd_a(age, year_list, shapefile_version)
    d <- get_gbd_d(age, year_list, shapefile_version)
    gbd <- merge(a, d , by = c("age_group_id", "year_id", "location_id", "gaul"))
    gbd <- combine_to_5a0(gbd)
    
  } else if(age == "infant"){
    a <- get_gbd_a(age, year_list, shapefile_version)
    d <- get_gbd_d(age, year_list, shapefile_version)
    gbd <- merge(a, d , by = c("age_group_id", "year_id", "location_id", "gaul"))
    
  } else if(age == "neonatal"){
    gbd <- get_gbd_d(age, year_list, shapefile_version)
    gbd[, a := NA]
  }
  
  return(gbd)
}


#' Get 4a0 and 10 weighted with deaths to get 5a0
#' 
#' Description: GBD life tables only have a for 4a1 and 1a0, so they need to be combined.
#' Calculates a, population, and deaths for under-5.
#'
#' @param gbd output of `get_gbd_crosswalk_values()`
#' 
#' @return an updated gbd table with 5x0 values
#'
combine_to_5a0 <- function(gbd) {
  
  message("Combining 4a1 and 1a0 to 5a0")
  
  #split by age
  under1 <- gbd[age_group_id == 28,]
  child <- gbd[age_group_id == 5,]
  
  setnames(under1, c("deaths", "a"), c("under1_deaths", "under1_a"))
  setnames(child, c("deaths", "a"), c("child_deaths", "child_a"))
  #add 1 to 4a1 so it represents average age of death rather than time of death after entering age bin
  child$child_a <- child$child_a + 1
  
  #merge them back together
  under5 <- merge(under1, child, by = c("year_id", "location_id"))
  under5 <- under5[,c("year_id", "location_id", "under1_a", "under1_deaths", "child_a", "child_deaths", "gaul.x")]
  setnames(under5, c("gaul.x"), c("gaul"))
  
  #calculate a, deaths, and population for 5x0
  under5[, a := (((under1_deaths * under1_a) + (child_deaths * child_a)) / (under1_deaths + child_deaths))][,deaths := under1_deaths + child_deaths]
  under5[,c("year_id", "location_id", "gaul", "a", "deaths")]
  return(under5)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert m cell pred to deaths
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Convert m cell pred to deaths
#' 
#' Description: Convert a cell pred with mortality rates to deaths. 
#'
#' @param cell_pred a cell pred with mortality rates
#' @param simple_raster simple raster for the corresponding region
#' @param year_list list of years in the cell_pred
#' @param pop a population raster that matches the simple raster
#' @param overs a list of the column names in the cell pred corresponding to draws
#' @param age character string, `under5`, `infant`, `neonatal` acceptable
#' (typically V1:Vndraws)
#' 
#' @return a cell pred of deaths, includes a column for years, adm0 gaul code from 
#' simple raster, and population from pop.
#' 
calculate_deaths_cell_pred <- function(cell_pred, simple_raster, year_list, pop, overs,age) {
  
  cell_pred <- as.data.table(cell_pred)
  #get ids corresponding to rows in cell pred
  simple_pixel_id <- which(!is.na(getValues(simple_raster)))
  
  #prepare ADM0 gaul code
  simple_list <- raster::extract(simple_raster, extent(simple_raster))
  simple_list <- na.omit(simple_list)
  simple_list <- rep(simple_list, times = length(year_list))
  #prepare years
  year = as.vector(unlist(lapply(year_list, function(x) rep.int(x, times = length(simple_pixel_id)))))
  #prepare population
  pop_list <- raster::extract(pop, extent(pop))
  pop_list <- pop_list[simple_pixel_id,]
  pop_list <- as.vector(pop_list)
  
  #add on gaul code, years, and pop to cell_pred
  cell_pred_tag <- cbind(cell_pred, simple_list, year, pop_list)
  cell_pred_tag <- as.data.table(cell_pred_tag)
  
  #multiply all draws by population to get deaths
  cell_pred_deaths <- cell_pred_tag[, (overs) := lapply(overs, function(x) get(x)* pop_list)]
  
  ##temporary fix for inflated border pixel deaths.  
  if(age == "infant"){
    cell_pred_deaths / 5
  } else if(age == "neonatal"){
    cell_pred_deaths / 60
  }
  
  return(cell_pred_deaths)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Post raking death functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Save death crosswalk outputs
#'
#' @param output_list a named list with 8 objects -
#'          -"deaths cell pred"
#'          -"mean": mean raster of deaths
#'          -"upper": upper raster of deaths
#'          -"lower": lower raster of deaths
#'          -"admin_0": data.table with draws aggregated to admin 0
#'          -"admin_1": data.table with draws aggregated to admin 1
#'          -"admin_2": data.table with draws aggregated to admin 2
#'          -"raking_factors": raking factors used in fractional raking
#' @param age "under5", "infant", or "neonatal"
#' @param run_date run date of folder to save in
#' @param region name of modelling region
#' 
#' @return NULL
#'
save_death_crosswalk_outputs <- function(output_list, 
                                         age, 
                                         run_date, 
                                         region) {

  out_dir <- sprintf("<<<< FILEPATH REDACTED >>>>")
  saveRDS(output_list[["raked_cell_pred"]], "<<<< FILEPATH REDACTED >>>>")
  writeRaster(output_list[["mean"]], file= "<<<< FILEPATH REDACTED >>>>", format = "GTiff", overwrite=T)
  writeRaster(output_list[["upper"]], file="<<<< FILEPATH REDACTED >>>>", format = "GTiff", overwrite=T)
  writeRaster(output_list[["lower"]], file= "<<<< FILEPATH REDACTED >>>>", format = "GTiff", overwrite=T)
  saveRDS(output_list[["raked_adm0_draws"]], "<<<< FILEPATH REDACTED >>>>")
  saveRDS(output_list[["raked_adm1_draws"]], "<<<< FILEPATH REDACTED >>>>")
  saveRDS(output_list[["raked_adm2_draws"]], "<<<< FILEPATH REDACTED >>>>")
  write.csv(output_list[["raking_factors"]], "<<<< FILEPATH REDACTED >>>>")
}