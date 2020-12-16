#----HEADER-------------------------------------------------------------------------------------------------------------
# Date:    November 2018
# Purpose: For LBD - reshape data from long to wide, resample polygons to points, and save by vaccine.
# Run:     source("FILEPATH/resample.r")
# Inputs:  nid  ->  NID from dataset to tabulate; saved population raster
#                                       vaccine prefixes to resample data for    
#***********************************************************************************************************************

message("Beginning resample")
require(fasterize, lib.loc = "/FILEPATH") #fasterize is used for a modification in the resampling process.
#Brandon's modifications to run the popraster separately and load from and RDs may make this unnecessary
require(sf, lib.loc = "FILEPATH" )      #pairs with fasterize

build_simple_raster_pop_fast <- function(subset_shape,u5m=FALSE,field="ADM0_CODE") {
  master_pop <- brick('FILEPATH/WorldPop_total_global_stack.tif') #WorldPop_allStages_stack.tif')
  
  tifs <- list.files(path = "FILEPATH", pattern = ".tif$", full.names = TRUE)
  tifs <- tifs[32:51]
  tifs <- tifs[c(1, 6, 11, 16)]
  
  master_pop <- do.call(brick, lapply(tifs, raster))
  cropped_pop <- crop(master_pop, extent(subset_shape), snap="out")
  
  ## Using fasterize
  initial_raster <- fasterize(st_as_sf(subset_shape), cropped_pop[[1]], field = field)
  
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster))))!=0) {
    rasterized_shape <- merge(fasterize(st_as_sf(subset(subset_shape, !(get(field) %in% unique(initial_raster)))), cropped_pop[[1]], field = field), initial_raster)
  }
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster))))==0) {
    rasterized_shape <- initial_raster
  }
  #rasterized_shape <- rasterize(subset_shape, cropped_pop, field='GAUL_CODE')
  masked_pop <- raster::mask(x=cropped_pop, mask=rasterized_shape)
  
  raster_list <- list()
  raster_list[['simple_raster']] <- rasterized_shape
  raster_list[['pop_raster']] <- masked_pop
  
  return(raster_list)
  
}

message('gaul_list = get_adm0_codes("all", shapefile_version = "current") ')
gaul_list = get_adm0_codes('all', shapefile_version = 'current')
shapefile_version = "current"
use_1k_popraster = TRUE


#### loads popraster. Takes a long time, so saves to file location which can be read from in future runs.
if(!exists("popraster")){
 simpp <- suppressMessages(suppressWarnings(load_simple_polygon(gaul_list = gaul_list,buffer=0.4,subset_only=TRUE,
                                                              shapefile_version = shapefile_version)))
 
  raster_list<-suppressMessages(suppressWarnings(build_simple_raster_pop_fast(simpp$subset_shape)))
  if(use_1k_popraster){
    popraster <- disaggregate(raster_list[['pop_raster']],5) # needs to be 1km for density 0.001
  } else {
    popraster = raster_list[['pop_raster']]
  }
 }

message("resample_geo_polys")
resample_geo_polys <- function(nid, extraction_root = "FILEPATH", 
                               prefixes = c("mcv", "dpt", "pcv", "rota", "polio", "hepb", "hib", "rcv", "bcg", "yfv", "hib3_dpt3_ratio", "hepb3_dpt3_ratio", "dpt3_timeliness_ratio"),
                               crop_to_region="all", collapse_method="kish") {
  
  ######################################################
  # Setup ##############################################
  ######################################################
  
  fn <- paste0(extraction_root, "/tabulated/lbd/", nid, ".csv")
  
  if (file.exists(fn)) { #Display console message if nid does not have an LBD tabulated dataset
    df_tabulated <- fread(paste0(extraction_root, "/tabulated/lbd/", nid, ".csv"))
    for (n in c("SSW", "N_obs")) {  
      if (!n %in% names(df_tabulated)) {
        df_tabulated[ ,(n):=NA]
      }
    }
    
    if (record_data_drop) {
      filter_table <- load_filter_table()
    }
    
    df_tabulated <- df_tabulated[me_name != "vacc_rcv2"] 
    
    if (max(df_tabulated$year_id) >= 1999){  
      for (vax_prefix in prefixes) {  ###Vaccine loop
        message(paste0("Beginning ", vax_prefix, " in Vaccine loop" ))
        
        df_input <- df_tabulated
        #Clean up data types
        df_input <- df_input[, latitude := as.numeric(latitude)]
        df_input <- df_input[, longitude := as.numeric(longitude)]
        
        df_input$survey_name <- gsub("/", "_", df_input$survey_name) #Clean up country specific survey names
        
        ### Create vaccine variable
        
        # Figure out which doses are present in the data & grab max
        dose_vars <- unique(df_input$me_name)[grepl(vax_prefix, unique(df_input$me_name))] #ID dose variables
        
        if (!grepl("ratio", vax_prefix)) {
          dose_vars <- dose_vars[!grepl("ratio", dose_vars)]
        }
        
        if (length(dose_vars) > 0) { #Displays message to console and ends loop iteration if no dose variables for given prefix
          max_doses <- str_match(dose_vars, paste0("vacc_", vax_prefix, "([0-9])"))[,2] %>%
            as.numeric %>% max ###IDs max dose for numeric vaccines - adjust for rotac
          
          #subset vaccine of interest
          df_input <- df_input[me_name %in% dose_vars]
          #reshape wide
          idvars <- c(names(df_input)[!names(df_input) %in% c("me_name", "value")]) #me_name is vaccine name, will use to create columns
          #test <- reshape(df_input, timevar = "me_name", idvar = idvars, direction="wide" )
          df_input <- data.table(data.table::dcast(df_input, svy_id + country + point + latitude + longitude + location_code +
                                                     shapefile + survey_name + year_id + psu + weight ~ me_name, sum,value.var = c("value", "N", "N_obs", "SSW"), fill = NA))
          
          
          val_vars <- paste0("value_", dose_vars) #existing names for new columns to be changed
          N_vars <- paste0("N_", dose_vars)
          N_obs_vars <- paste0("N_obs_", dose_vars)
          SSW_vars <- paste0("SSW_", dose_vars)
          names(df_input)[names(df_input) %in% val_vars] <- dose_vars #change dose column names
          
          ###ID variables by prefix, rename and drop as nec. 
          
          if(vax_prefix == "hib3_dpt3_ratio"){
            names(df_input)[names(df_input) == "N_hib3_dpt3_ratio"]      <- "N"
            names(df_input)[names(df_input) == "N_obs_hib3_dpt3_ratio"]  <- "N_obs"
            names(df_input)[names(df_input) == "SSW_hib3_dpt3_ratio"]    <- "SSW"
          }
          
          if(vax_prefix == "hepb3_dpt3_ratio"){
            names(df_input)[names(df_input) == "N_hepb3_dpt3_ratio"]      <- "N"
            names(df_input)[names(df_input) == "N_obs_hepb3_dpt3_ratio"]  <- "N_obs"
            names(df_input)[names(df_input) == "SSW_hepb3_dpt3_ratio"]    <- "SSW"
          }
          
          if(vax_prefix == "mcv2_mcv1_ratio"){
            names(df_input)[names(df_input) == "N_mcv2_mcv1_ratio"]      <- "N"
            names(df_input)[names(df_input) == "N_obs_mcv2_mcv1_ratio"]  <- "N_obs"
            names(df_input)[names(df_input) == "SSW_mcv2_mcv1_ratio"]    <- "SSW"
          }
          
          if(vax_prefix == "pcv3_dpt3_ratio"){
            names(df_input)[names(df_input) == "N_pcv3_dpt3_ratio"]      <- "N"
            names(df_input)[names(df_input) == "N_obs_pcv3_dpt3_ratio"]  <- "N_obs"
            names(df_input)[names(df_input) == "SSW_pcv3_dpt3_ratio"]    <- "SSW"
          }
          
          if(vax_prefix == "rotac_dpt3_ratio"){
            names(df_input)[names(df_input) == "N_rotac_dpt3_ratio"]      <- "N"
            names(df_input)[names(df_input) == "N_obs_rotac_dpt3_ratio"]  <- "N_obs"
            names(df_input)[names(df_input) == "SSW_rotac_dpt3_ratio"]    <- "SSW"
          }
          
          if(vax_prefix == "dpt3_timeliness_ratio"){
            names(df_input)[names(df_input) == "vacc_dpt3_timeliness_ratio"]       <- "dpt3_timeliness_ratio"
            names(df_input)[names(df_input) == "N_vacc_dpt3_timeliness_ratio"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_dpt3_timeliness_ratio"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_dpt3_timeliness_ratio"]   <- "SSW"
          }
          
          if(vax_prefix == "mcv"){
            names(df_input)[names(df_input) == "vacc_mcv1"]       <- "mcv_dose_1"
            names(df_input)[names(df_input) == "vacc_mcv2"]       <- "mcv_dose_2"
            
            names(df_input)[names(df_input) == "N_vacc_mcv1"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_mcv1"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_mcv1"]   <- "SSW"
            # Rename
            df_input[ , mcv_dose_0 := N - mcv_dose_1]
            if ("mcv_dose_2" %in% names(df_input)) df_input[!is.na(mcv_dose_2), mcv_dose_1 := mcv_dose_1 - mcv_dose_2]
          }
          
          if(vax_prefix == "dpt"){
            names(df_input)[names(df_input) == "vacc_dpt1"]       <- "dpt_dose_1"
            names(df_input)[names(df_input) == "vacc_dpt2"]       <- "dpt_dose_2"
            names(df_input)[names(df_input) == "vacc_dpt3"]       <- "dpt_dose_3"
            names(df_input)[names(df_input) == "N_vacc_dpt1"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_dpt1"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_dpt1"]        <- "SSW"
            df_input[ , dpt_dose_0 := N - dpt_dose_1]
            df_input[ , dpt_dose_1 := dpt_dose_1 - dpt_dose_2]
            df_input[ , dpt_dose_2 := dpt_dose_2 - dpt_dose_3]
          }
          
          if(vax_prefix == "pcv"){
            names(df_input)[names(df_input) == "vacc_pcv1"]       <- "pcv_dose_1"
            names(df_input)[names(df_input) == "vacc_pcv2"]       <- "pcv_dose_2"
            names(df_input)[names(df_input) == "vacc_pcv3"]       <- "pcv_dose_3"
            names(df_input)[names(df_input) == "N_vacc_pcv1"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_pcv1"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_pcv1"]        <- "SSW"
            df_input[ , pcv_dose_0 := N - pcv_dose_1]
            df_input[ , pcv_dose_1 := pcv_dose_1 - pcv_dose_2]
            df_input[ , pcv_dose_2 := pcv_dose_2 - pcv_dose_3]
          }
          
          if(vax_prefix == "rota"){
            names(df_input)[names(df_input) == "vacc_rota1"]       <- "rota_dose_1"
            names(df_input)[names(df_input) == "vacc_rota2"]       <- "rota_dose_2"
            names(df_input)[names(df_input) == "vacc_rota3"]       <- "rota_dose_3"
            names(df_input)[names(df_input) == "vacc_rotac"]       <- "rota_dose_c"
            if("rota_dose_c" %in% names (df_input)) {  
              names(df_input)[names(df_input) == "N_vacc_rotac"]     <- "N"
              names(df_input)[names(df_input) == "N_obs_vacc_rotac"] <- "N_obs"
              names(df_input)[names(df_input) == "SSW_vacc_rotac"]        <- "SSW"
            } else {
              names(df_input)[names(df_input) == "N_vacc_rota1"]     <- "N"
              names(df_input)[names(df_input) == "N_obs_vacc_rota1"] <- "N_obs"
              names(df_input)[names(df_input) == "SSW_vacc_rota1"]        <- "SSW"
              df_input[ , rota_dose_0 := N - rota_dose_1]
              df_input[ , rota_dose_1 := rota_dose_1 - rota_dose_2]
              df_input[ , rota_dose_2 := rota_dose_2 - rota_dose_3]
            }
          }
          
          if(vax_prefix == "polio"){
            names(df_input)[names(df_input) == "vacc_polio1"]       <- "polio_dose_1"
            names(df_input)[names(df_input) == "vacc_polio2"]       <- "polio_dose_2"
            names(df_input)[names(df_input) == "vacc_polio3"]       <- "polio_dose_3"
            names(df_input)[names(df_input) == "N_vacc_polio1"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_polio1"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_polio1"]        <- "SSW"
            df_input[ , polio_dose_0 := N - polio_dose_1]
            df_input[ , polio_dose_1 := polio_dose_1 - polio_dose_2]
            df_input[ , polio_dose_2 := polio_dose_2 - polio_dose_3]
          }
          if(vax_prefix == "hepb"){
            names(df_input)[names(df_input) == "vacc_hepb1"]       <- "hepb_dose_1"
            names(df_input)[names(df_input) == "vacc_hepb2"]       <- "hepb_dose_2"
            names(df_input)[names(df_input) == "vacc_hepb3"]       <- "hepb_dose_3"
            names(df_input)[names(df_input) == "N_vacc_hepb1"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_hepb1"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_hepb1"]   <- "SSW"
            df_input[ , hepb_dose_0 := N - hepb_dose_1]
            df_input[ , hepb_dose_1 := hepb_dose_1 - hepb_dose_2]
            df_input[ , hepb_dose_2 := hepb_dose_2 - hepb_dose_3]
          }
          if(vax_prefix == "hib"){
            names(df_input)[names(df_input) == "vacc_hib1"]       <- "hib_dose_1"
            names(df_input)[names(df_input) == "vacc_hib2"]       <- "hib_dose_2"
            names(df_input)[names(df_input) == "vacc_hib3"]       <- "hib_dose_3"
            names(df_input)[names(df_input) == "N_vacc_hib1"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_hib1"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_hib1"]   <- "SSW"
            df_input[ , hib_dose_0 := N - hib_dose_1]
            df_input[ , hib_dose_1 := hib_dose_1 - hib_dose_2]
            df_input[ , hib_dose_2 := hib_dose_2 - hib_dose_3]
          }
          if(vax_prefix == "rcv"){
            names(df_input)[names(df_input) == "vacc_rcv1"]       <- "rcv_dose_1"
            names(df_input)[names(df_input) == "N_vacc_rcv1"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_rcv1"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_rcv1"]   <- "SSW"
          }
          if(vax_prefix == "bcg"){
            names(df_input)[names(df_input) == "vacc_bcg"]       <- "bcg_dose_1"
            names(df_input)[names(df_input) == "N_vacc_bcg"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_bcg"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_bcg"]        <- "SSW"
            df_input[ , bcg_dose_0 := N - bcg_dose_1]
          }
          if(vax_prefix == "yfv"){
            names(df_input)[names(df_input) == "vacc_yfv"]       <- "yfv_dose_1"
            names(df_input)[names(df_input) == "N_vacc_yfv"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_yfv"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_yfv"]   <- "SSW"
            df_input[ , yfv_dose_0 := N - yfv_dose_1]
          }
          
          
          # Remove other vaccine dose variables
          message("Remove other vaccine dose variables")
          
          
          if (grepl("ratio", vax_prefix)) {
            df_input <- subset(df_input, select = !(names(df_input) %in% c(paste0("N_", dose_vars),
                                                                           paste0("N_obs_", dose_vars),
                                                                           paste0("SSW_", dose_vars))))
          } else {
            df_input <- subset(df_input, select = !(names(df_input) %in% c(dose_vars, 
                                                                           paste0("N_", dose_vars),
                                                                           paste0("N_obs_", dose_vars),
                                                                           paste0("SSW_", dose_vars))))
          }
          
          ####### Save Design Effect information
          message(paste0("Save design effect information for ", vax_prefix ))
          
          if(nrow(df_input[point==0, ]) > 0){  #DE for polygons only
            de_data <- data.table(read.csv("FILEPATH/design_effect_nids.csv"))  #existing DE data to append to/modify
            sid = nid  ### so we don't have column name/object name confusion later
            N_sum <- sum(df_input[point == 0, N]) #Sum of weighted N in polygon data
            Nobs_sum <- sum(df_input[point == 0, N_obs]) #Sum of unweighted N in polygon data
            de_calc <- N_sum/Nobs_sum ###Survey level DE is WtN/UnWtN
            if(paste0(nid, vax_prefix) %in% unique(paste0(de_data$nid, de_data$vacc))){ ###If NID/vaccine exists in sheet, modify that row
              de_data[ nid == sid  & vacc == vax_prefix, N := N_sum]
              de_data[ nid == sid  & vacc == vax_prefix, N_obs := Nobs_sum]
              de_data[ nid == sid  & vacc == vax_prefix, de := de_calc]
            } else {  ### IF NID/vaccine pair does not exists in DE data, append it
              de_data <- rbind(de_data, data.table("nid" = sid, "vacc" = vax_prefix, "N" = N_sum, "N_obs" = Nobs_sum, "de" = de_calc)) 
            }
            write.csv(de_data, "FILEPATH/design_effect_nids.csv", row.names=FALSE) ###Save the DE data back to the same place with modification/append
          }
          #######
          
          
          if (nrow(df_input[point == 0]) > 0) {
            # Run the time-intensive point/poly solution 
            source(paste0(core_repo, "FILEPATH/polygon_functions.R"))
            source(paste0(core_repo, "FILEPATH/shapefile_functions.R"))
            source(paste0(core_repo, "FILEPATH/prep_functions.R"))
            
            df_pointpoly <- resample_polygons(data = df_input, 
                                              cores = 10, 
                                              indic = vax_prefix, 
                                              density = 0.001,
                                              gaul_list = get_adm0_codes('all', shapefile_version = 'current'),
                                              pull_poly_method = "fast")
            
          } else {
            
            # If no polys, don't resample!
            
            df_pointpoly <- df_input
          }
          
          #Save resampled data
          if (!dir.exists( file.path(paste0(extraction_root, "/FILEPATH/"), vax_prefix)) ) dir.create( file.path(paste0(extraction_root, "/FILEPATH/"), vax_prefix))
          
          write.csv(df_pointpoly, paste0(extraction_root, "/FILEPATH/", vax_prefix, "/", nid, ".csv"), row.names=FALSE)
          
          message(paste0(vax_prefix, " resampled"))
          
          
        } else {
          message(paste0(nid, " ", vax_prefix, " has no dose variables, skipping resample"))
          if(record_data_drop & vax_prefix == "mcv"){
            message("Filter Table: No MCV data, recording in missing_indicator, setting all else to NA")
            
            filter_table[, missing_indicator := TRUE]
            filter_table[, record_filter_successful := TRUE]
            filter_table[, setdiff(names(filter_table), c("survey_id", "missing_indicator", "survey_in_year_range", "outlier", "record_filter_successful")) := NA]
            
            record_data_drop <- FALSE
            write.csv(filter_table, file = paste0(filter_table_path, nid, ".csv"), row.names = FALSE)
          }
        }
      } #close vaccine loop
    } else {
      message(paste0(nid, " all pre-2000, skipping resample"))
      
      if(record_data_drop){
        
        filter_table[, survey_in_year_range     := FALSE]
        filter_table[, record_filter_successful := TRUE]
        filter_table[, setdiff(names(filter_table), c("survey_id", "missing_indicator", "survey_in_year_range", "outlier", "record_filter_successful")) := NA]
        record_data_drop <- FALSE
        
        write.csv(filter_table, file = paste0(filter_table_path, nid, ".csv"), row.names = FALSE)
      }
    }
  } else {
    message(paste0(nid, " not tabulated for lbd. skipping resampling step."))
  }
} 

