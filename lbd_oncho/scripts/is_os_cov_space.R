## is_os_cov_space ###################################################

#' Diagnostic plots digging into the covariate range of data input locations vs. the covariate range of the modeling region
#'
#' @param effects_l Either character string with '+' separating the covariate names or as a vector
#' @param measures_l Either character string with '+' separating the covariate measure names or as a vector
#' @param year_list integer vector of years to load for covariates
#' @param geo string of geographic extent. Fed into get_adm0_codes
#' @param data data.frame/data.table; in general MBG data input format
#' @param year_load specifying options in covariate loading
#' @param interval_mo specifying optios in covariate loading
#' @param centrescale should covariates be standardized during loading?
#' @param num_bins number of bins for histogram plotting and for determining gaps in data covariate range v. region covariate range
#' @param out_dir either file_path or "default" which is specified as "REDACTED"
#' @param plot_pdf boolean; whether or not a pdf output should be created
#'
#' @return returns a list; [[1]] = raster summarizing number of out of range covariates per pixel in specified region; [[2]] = comparison histograms; [[3]] = binary rasters (out of range or not) per covariate

user <- Sys.info()[["user"]]

## Set repo locations
core_repo <- <<<< FILEPATH REDACTED >>>>
indic_repo <- <<<< FILEPATH REDACTED >>>>

path <- paste0("/homes/", user, "/rlibs/")

indicator_group <- "oncho"
indicator <- "had_oncho_w_resamp"
run_date <- "2020_08_08_15_01_40"
region <- "oncho_endem_afr"
holdout <- 0

outputdir <- file.path(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

mbg_setup(package_list = package_list, repos = core_repo)

library(fasterize)
library(matrixStats)
library(sf)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(<<<< FILEPATH REDACTED >>>>)
}

config <- set_up_config_focal_3(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator,
                                config_file = <<<< FILEPATH REDACTED >>>>, run_tests = FALSE)

effects_l <- fixed_effects
measures_l <- fixed_effects_measures
if (class(year_list) == "character") year_list <- eval(parse(text = year_list))
geo <- region

data <- fread(<<<< FILEPATH REDACTED >>>>)

num_bins <- 100
out_dir <- "default"
plot_pdf <- TRUE

is_os_cov_space <- function(fixed_effects_config, year_list, geo, data, yearload = "annual", interval_mo = "12", centrescale = FALSE, num_bins = 250, out_dir = "default", plot_pdf = TRUE, use_river_size = TRUE, use_oncho_suitability = TRUE) {

  message(paste0("Covariates = ", paste(fixed_effects, collapse = " + ")))
  message(paste0("Measures = ", paste(fixed_effects_measures[[2]], collapse = " + ")))

  if (yearload == "annual") {
    period_map <-
      make_period_map(modeling_periods = c(min(year_list):max(year_list)))
  }
  if (yearload == "five-year") {
    period_map <-
      make_period_map(modeling_periods = seq(min(year_list), max(year_list), by = 5))
  }

  ## Load simple polygon template to model over
  gaul_list <- get_adm0_codes(geo, shapefile_version = modeling_shapefile_version)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = F)
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]

  ## Build administrative and population rasters
  raster_list <- build_simple_raster_pop(subset_shape)
  simple_raster <- raster_list[["simple_raster"]]
  
  cov_layers <- load_lagged_covariates(covariate_config = fixed_effects_config, template = simple_raster, start_year = min(year_list), end_year = max(year_list), raster_agg = 1)

  if (as.logical(use_river_size)) {
    river_size <- raster(<<<< FILEPATH REDACTED >>>>)
    river_size <- projectRaster(river_size, crs = crs(cov_layers[[1]][[1]]))
    river_size <- raster::resample(river_size, cov_layers[[1]][[1]])
    river_size <- raster::crop(river_size, cov_layers[[1]][[1]])
    river_size[is.na(river_size[])] <- 0
    river_size <- raster::mask(river_size, cov_layers[[1]][[1]])
    names(river_size) <- "river_size"
    
    cov_layers[["river_size"]] <- river_size
  }
  
  if (as.logical(use_oncho_suitability)) {
    oncho_suitability <- stack(<<<< FILEPATH REDACTED >>>>)
    oncho_suitability <- subset(oncho_suitability, 1)
    oncho_suitability <- projectRaster(oncho_suitability, crs = crs(cov_layers[[1]][[1]]))
    oncho_suitability <- raster::resample(oncho_suitability, cov_layers[[1]][[1]])
    oncho_suitability <- raster::crop(oncho_suitability, cov_layers[[1]][[1]])
    oncho_suitability[is.na(oncho_suitability[])] <- 0
    oncho_suitability <- raster::mask(oncho_suitability, cov_layers[[1]][[1]])
    names(oncho_suitability) <- "oncho_suitability"
    
    cov_layers[["oncho_suitability"]] <- oncho_suitability
  }
  
  # copy the dataset to avoid unintended namespace conflicts
  the_data <- copy(data) %>% as.data.table()

  # add a row id column
  the_data[, a_rowid := seq(1:nrow(the_data))]

  # extract covariates to the points and subset data where its missing covariate values
  cs_covs <- extract_covariates(the_data,
    cov_layers,
    id_col = "a_rowid",
    return_only_results = TRUE,
    centre_scale = centrescale,
    period_var = "year",
    period_map = period_map
  )

  ## Check for data where covariate extraction failed
  rows_missing_covs <- nrow(the_data) - nrow(cs_covs)
  if (centrescale) rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]])
  if (rows_missing_covs > 0) {
    pct_missing_covs <- round((rows_missing_covs / nrow(the_data)) * 100, 2)
    warning(paste0(
      rows_missing_covs, " out of ", nrow(the_data), " rows of data ",
      "(", pct_missing_covs, "%) do not have corresponding ",
      "covariate values and will be dropped from child models..."
    ))
    if (rows_missing_covs / nrow(the_data) > 0.1) {
      stop(paste0(
        "Something has gone quite wrong: more than 10% of your data does not have ",
        "corresponding covariates.  You should investigate this before proceeding."
      ))
    }
  }

  if (centrescale) {
    the_data <- merge(the_data, cs_covs[[1]], by = "a_rowid", all.x = F, all.y = F)
    covs <- copy(cs_covs[[1]])
  } else {
    the_data <- merge(the_data, cs_covs, by = "a_rowid", all.x = F, all.y = F)
    covs <- copy(cs_covs)
  }

  # store the centre scaling mapping
  if (centrescale) covs_cs_df <- cs_covs[[2]]

  the_data <- as.data.table(the_data)

  ### Make histograms and map output objects############################
  bicolor_hist <- vector("list", length = length(cov_layers))
  summarized_layers <- vector("list", length = length(cov_layers))
  summary <- copy(cov_layers[[1]][[1]])
  values(summary) <- 0
  summary <- raster::mask(summary, cov_layers[[1]][[1]])
  
  # these objects are mostly for debugging purposes
  d <- vector("list", length = length(cov_layers))
  matrix_l <- vector("list", length = length(cov_layers))

  for (i in 1:length(cov_layers)) {
    print(i)
    nlay <- nlayers(cov_layers[[i]])

    if (nlay == 1) {
      cov_value <- data.table(value = values(cov_layers[[i]]))
      summarized_layers[[i]] <- copy(cov_layers[[i]])
    } else {
      cov_value <- lapply(1:nlay, function(l) data.table(value = values(cov_layers[[i]][[l]]))) %>% rbindlist()
      summarized_layers[[i]] <- calc(cov_layers[[i]], fun = mean)
    }
    cov_value <- cov_value[is.na(value) == F, ]
    d[[i]] <- the_data[, names(cov_layers)[i], with = F]

    bicolor_hist[[i]] <- ggplot() +
      stat_bin(data = cov_value, aes(x = value, y = ..count.. / sum(..count..), fill = "In Region"), bins = num_bins) + theme_classic() +
      stat_bin(data = d[[i]], aes(x = get(names(cov_layers)[i]), y = ..count.. / sum(..count..), fill = "In Data"), bins = num_bins, alpha = .75) +
      labs(title = paste0("IS/OS covariate values: ", names(cov_layers)[i]), x = "Values", y = "Proportion in data source")# + guides(fill = "none")

    # converting continuous to discrete
    a_step <- (max(cov_value$value, na.rm = T) - min(cov_value$value, na.rm = T)) / num_bins
    lower <- seq(from = min(cov_value$value, na.rm = T), to = max(cov_value$value, na.rm = T) - a_step, length.out = num_bins)
    upper <- c(lower[2:length(lower)], max(cov_value$value, na.rm = T))
    matrix_l[[i]] <- cbind(lower, upper, id = c(1:num_bins)) %>% as.data.table()
    matrix_l[[i]]$cat <- cut(x = upper, breaks = c(lower[1], upper), include.lowest = TRUE)

    d[[i]]$cat <- cut(unlist(d[[i]]), breaks = c(lower[1], upper), include.lowest = TRUE)
    d[[i]] <- merge(d[[i]], matrix_l[[i]], by = "cat")

    matrix_l[[i]][(id %in% d[[i]]$id) == T, out_data := 0]
    matrix_l[[i]][is.na(out_data), out_data := 1]
    
    lower_quantile <- 0.025 #0.0
    upper_quantile <- 0.975 #1.0
    
    quants <- quantile(as.data.frame(the_data)[, names(cov_layers)[i]], probs = c(lower_quantile, upper_quantile), na.rm = TRUE)
    
    summary[values(summarized_layers[[i]]) < quants[1]] <- summary[values(summarized_layers[[i]]) < quants[1]] + 1
    summary[values(summarized_layers[[i]]) > quants[2]] <- summary[values(summarized_layers[[i]]) > quants[2]] + 1
    
    summarized_layers[[i]] <- reclassify(summarized_layers[[i]], matrix_l[[i]][, .(lower, upper, out_data)], include.lowest = T)
  }

  names(bicolor_hist) <- names(cov_layers)
  names(summarized_layers) <- names(cov_layers)
  names(d) <- names(cov_layers)
  names(matrix_l) <- names(cov_layers)

  output_names <- names(cov_layers)
  sum_raster <- brick(summarized_layers[output_names])
  sum_raster <- calc(sum_raster, fun = sum)

  plot(summary) # This is the figure we want
  writeRaster(summary, <<<< FILEPATH REDACTED >>>>, format = "GTiff")
  
  ### Plotting and output ###################################
  if (plot_pdf) {
    library(viridis)
    library(rasterVis)
    library(grid)
    library(gridExtra)
    library(RColorBrewer)

    theme_empty <- theme_classic() +
      theme(
        axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 15)
      )

    colr <- brewer.pal(9, "RdYlGn")

    if (out_dir == "default") {
      outdir <- <<<< FILEPATH REDACTED >>>>
    } else {
      outdir <- out_dir
    }

    max_v <- max(values(sum_raster$layer), na.rm = T)

    pdf(file = <<<< FILEPATH REDACTED >>>>)
    levelplot(sum_raster,
      xlab = NULL, ylab = NULL, main = list("Extrapolation - # of covariates out of range", color = "black"),
      margin = FALSE, # suppress marginal graphics
      colorkey = list(
        space = "bottom", # plot legend at bottom
        labels = list(at = 0:max_v, font = 4) # legend ticks and labels
      ),
      par.settings = list(
        axis.line = list(col = "transparent") # suppress axes and legend outline
      ),
      scales = list(draw = FALSE), # suppress axis labels
      col.regions = colr, # colour ramp
      at = seq(0, max_v, by = 1)
    ) + # colour ramp breaks
      layer(sp.polygons(subset_shape, lwd = 1))

    for (n in output_names) {
      i <- which(names(cov_layers) == n)

      p <- levelplot(summarized_layers[[i]],
        xlab = NULL, ylab = NULL, main = list(paste0("Out of range: ", names(cov_layers)[i]), color = "black"),
        margin = FALSE, # suppress marginal graphics
        colorkey = list(
          space = "bottom", # plot legend at bottom
          labels = list(at = 0:1, font = 4) # legend ticks and labels
        ),
        par.settings = list(
          axis.line = list(col = "transparent") # suppress axes and legend outline
        ),
        scales = list(draw = FALSE), # suppress axis labels
        col.regions = colr, # colour ramp
        at = seq(0, 1, len = 3)
      ) + # colour ramp breaks
        layer(sp.polygons(subset_shape, lwd = 1))

      # arrange plots side by side
      grid.arrange(bicolor_hist[[i]], p, ncol = 2)
    }
    dev.off()
  }

  pdf(file = <<<< FILEPATH REDACTED >>>>)
  plot(sum_raster)
  dev.off()
  
  return(list(sum_raster, bicolor_hist, summarized_layers))
}

oncho_endem_afr <- is_os_cov_space(fixed_effects_config = fixed_effects_config, year_list = year_list, geo = region, data = data, yearload = yearload, interval_mo = interval_mo, centrescale = FALSE, num_bins = num_bins, out_dir = out_dir, plot_pdf = plot_pdf)
