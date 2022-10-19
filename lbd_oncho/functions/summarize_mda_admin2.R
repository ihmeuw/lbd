## summarize_mda_admin2 function ###################################################

# This function generates a summarization of the MDA data
# to be used in the visualization function and add programmatic timelines to plots
#
# This function pulls the mdavector.numrounds data from the covariate directory. It then summarizes whether MDA occurred each year
# for each adm2 unit (based off of central shapefile)
#
# NOTE: each polygon is stripped of its outer layer after crop and mask in order to ensure that boundary pixels are not captured.
#
# @param cov What covariate is being summarized? Function currently geared towards MDA.
# @param measure What measure of the covariate should be used? Function currently geared towards "numrounds" (Number of rounds).
# @param reg Region to feed into simple_polygon & covariate loading.
# @param convert How to summarize covariate? Currently set to "boolean" (For a given admin2 & year, was MDA implemented?)
# @param yl param feeding into simple polygon & covariate laoding
# @param interval_mo
# @param update
# @param save
# @param save_filepath
#
# @return a single data table which is the input data with spatial characteristics assigned
#

summarize_mda_admin2 <- function(cov = "allmda",
                                 measure = "numrounds",
                                 reg = "oncho_endem_afr",
                                 convert = "boolean",
                                 yl = c(2000:2015),
                                 interval_mo = 12,
                                 update,
                                 save = FALSE,
                                 save_filepath = NULL,
                                 load_filepath = NULL,
                                 ind_repo = indic_repo,
                                 shapefile_version = NULL,
                                 ind_gp = indicator_group,
                                 indic = indicator,
                                 rd = run_date,
                                 age = 0, holdout = 0) {
  if (update) {
    configurations <- fread(<<<< FILEPATH REDACTED >>>>)
    config_sh_version <- configurations[V1 == "modeling_shapefile_version", V2]

    if (is.null(shapefile_version) == T) {
      shapefile_version <- config_sh_version
      message(paste0("Using modeling_shapefile_version specified in config: ", shapefile_version))
    }

    ### Function Prep ################################################################

    if (!(file.exists(<<<< FILEPATH REDACTED >>>>))) {
      gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)
      simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = F, shapefile_version = shapefile_version)
      subset_shape <- simple_polygon_list[[1]]
      simple_polygon <- simple_polygon_list[[2]]
      message("Building simple raster from subset_shape")
      raster_list <- build_simple_raster_pop(subset_shape)
      simple_raster <- raster_list[["simple_raster"]]
      pop_raster <- raster_list[["pop_raster"]]
    } else {
      message("-- loading spatial templates")

      load(<<<< FILEPATH REDACTED >>>>)
    }

    temp_config <- fixed_effects_config[covariate == cov]
    temp_config$year_lag <- 0
    
    cov_layers <- load_lagged_covariates(covariate_config = temp_config,
                                         template = simple_raster_raf1,
                                         start_year = min(predict_years),
                                         end_year = max(predict_years),
                                         raster_agg = as.integer(raster_agg_factor))
    covariate <- cov_layers[[cov]]
    
    subset_codes <- get_adm0_codes(reg, shapefile_version = shapefile_version)
    ad2_shape <- readOGR(get_admin_shapefile(admin_level = 2, raking = F, version = shapefile_version)) %>% subset(ADM0_CODE %in% subset_codes)
    ad2_data <- as.data.table(ad2_shape@data)
    covariate <- raster::crop(covariate, extent(simple_raster))
    covariate <- setExtent(covariate, simple_raster)
    covariate <- raster::mask(covariate, simple_raster)

    if (nlayers(covariate) < length(yl)) warning("The specified annual covariate does not match up with the number of years specified in year_list. Function will break. Please fix.")

    crop_mask_getvalue <- function(code, shape = ad2_shape, raster = covariate, conv = convert, yl = yl) {
      masked <- raster::crop(raster, extent(shape[shape@data$ADM2_CODE == code, ]))
      masked <- raster::mask(masked, shape[shape@data$ADM2_CODE == code, ])

      # remove border which could contain values from adjacent admin units
      inner <- raster::mask(masked[[length(yl)]], boundaries(masked[[length(yl)]]), maskvalue = 1)
      if (length(unique(raster::values(inner))) > 1) masked <- raster::brick(lapply(c(1:length(yl)), FUN = function(i) return(raster::mask(masked[[i]], boundaries(masked[[i]]), maskvalue = 1))))

      if (conv == "boolean") {
        vals <- lapply(c(1:length(yl)), FUN = function(i) {
          if (sum(raster::values(masked[[i]]), na.rm = T) > 0) {
            bool <- 1
          }
          else {
            bool <- 0
          }
          return(bool)
        })
      }
      return(unlist(vals))
    }

    ad2_list <- ad2_data$ADM2_CODE
    cov_values <- lapply(ad2_list, FUN = crop_mask_getvalue, yl = yl)
    cov_table <- data.table(matrix(unlist(cov_values), nrow = length(cov_values), byrow = T), stringsAsFactors = FALSE)
    names(cov_table) <- paste0("year_", yl)
    cov_table <- cbind(ad2_data, cov_table)
    year_start_index <- grep("year_2000", names(cov_table))
    cov_table$mda_start <- apply(cov_table, 1, FUN = function(v) {
      if (sum(as.numeric(v[year_start_index:ncol(cov_table)])) == 0) {
        return(NA)
      } else {
        return(which(v[year_start_index:ncol(cov_table)] == 1)[1] + 1999)
      }
    })

    if (save) {
      if (is.null(save_filepath)) {
        warning("Filepath for saving updated MDA admin 2 table was not specified. Please not that saving will not completed.")
      }
      else {
        write.csv(cov_table, save_filepath, row.names = F)
      }
    }

    return(cov_table)
  } else {
    if (save) "Saving unnecessary since update was not made"
    if (is.null(load_filepath) == T) {
      load_filepath <- <<<< FILEPATH REDACTED >>>>
      warning(paste0("Using default loading path: ", load_filepath))
    }
    cov_table <- fread(load_filepath)
    return(cov_table)
  }
}
