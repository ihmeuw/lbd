############################################################################################################################
## Crosswalking function for report data in seperate age ranges ############################################################
############################################################################################################################

#' crosswalk_prepare function ######################################
#'
#' This function takes the microdata collected for an indicator and then prepares it to be used in the crosswalk_report function below. It prepares the
#' data by first seperating out the point data and then assigning it to the correct admin1 unit by its lat/long using a GADM admin1 shapefile. It
#' then subsets the polygon data to those that have a recorded admin1 unit, and appends the point and polygon data together at the end of the function.
#'
#'
#' @param data  This is the microdata that is assigned to the correct admin1 level
#' @param indicator Model indicator name
#' @param weight Sample weight column name, typically pweight
#'
#'
#' @return
#' This function returns the microdata that is ready to be used in crosswalk_reports(). The microdata has the appropriate variables and
#' has been assigned an admin1, where the indicator is now assigned name indicator, and weight function is now called weight


crosswalk_data_prepare <- function(data, indicator = "hiv_test", weight = "hiv_weight", shapefile_version=shapefile_version){

  # Load required packages
  library(sf)

  # Subset data to only required variables, make sure none are missing
  data <- data.table(data)
  setnames(data, c(indicator, weight), c("indicator", "weight"))

  data <- data[, .(nid, country, survey_series, survey_name, year, strata, psu, point, shapefile, location_code,
                 latitude, longitude, sex_id, age_year, int_year, indicator, weight, admin_1, admin_level)]
  for (i in 1:ncol(data)) attributes(data[[i]]) <- NULL

  #correct for encoding errors in admin_1
  Encoding(data$admin_1) <- 'latin1'
  data$admin_1 <- stringi::stri_trans_general(data$admin_1, 'Latin-ASCII') %>%
    tolower()

  data <-
    data %>%
    mutate_at(c("point", "year", "latitude", "longitude"), as.numeric) %>%
    filter(!is.na(age_year), !is.na(indicator), !is.na(year), !is.na(nid), !is.na(weight)) %>%
    dplyr::rename(source = survey_series)

  # Drop countries outside of Africa
  source("<<<< FILEPATH REDACTED >>>>")
  loc <- data.table(get_location_metadata(location_set_id = 2, gbd_round_id = 4))
  loc <- loc[grepl("Africa", region_name) & location_type == "admin0", ihme_loc_id]
  data <- filter(data, country %in% loc)

  # Subset point data and assign to admin 1 ---------------------------------
  point_data <- filter(data, point == 1, !is.na(latitude),!is.na(longitude))

  # Read in admin 1 shapefile
  admin_shp <- st_read(get_admin_shapefile(admin_level = 1), quiet = T)

  # Get location name to ihme_lc_id to subset down shapefile for speed
  gadm_to_loc_id <-
    get_location_code_mapping("current") %>%
    dplyr::select(ADM_CODE, loc_name, ihme_lc_id) %>%
    dplyr::rename(location_name = loc_name) %>%
    filter(ADM_CODE %in% get_adm0_codes(c("cssa", "wssa", "essa_sdn", "sssa")))

  point_data <-
    point_data %>%
    left_join(gadm_to_loc_id, by = c("country" = "ihme_lc_id")) %>%
    mutate(location_name = ifelse(location_name== "Tanzania", "United Republic of Tanzania", location_name)) %>%
    mutate(location_name = ifelse(location_name == "Cote d'Ivoire", "Côte d'Ivoire", location_name)) %>% data.table()

  # Subset shapefile to include countries we have data on
  countries <- unique(point_data$location_name)

  admin_shp <-
    admin_shp %>%
    mutate(ADM0_NAME = as.character(ADM0_NAME)) %>%
    filter(ADM0_NAME %in% countries)

  # Assign each point to correct admin unit
  message("Assigning point data to the correct admin1 unit")

  point_data <-
    point_data %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(admin_shp)) %>%
    st_join(admin_shp) %>%
    st_set_geometry(NULL)

  # Again subset and rename some variables. Note there still are a couple (30 ish) cases where the lat long are assigned
  # to a different country then before
  point_data <-
    point_data %>%
    filter(location_name == ADM0_NAME) %>% # Drop where countries do not align (very few do this)
    dplyr::rename(admin1 = ADM1_NAME, country_long = location_name) %>%
    mutate(admin1 = tolower(admin1)) %>%
    dplyr::select(nid, country_long, country, source, year, strata,
                  psu, point, shapefile, location_code, sex_id, age_year,
                  int_year, indicator, weight, admin1) %>% data.table()

  # Make sure polygon data is assigned to the correct admin unit ---------------------------------------------------
  polygon_data <-
    data %>%
    filter(point == 0, !is.na(shapefile), !is.na(admin_1)) %>% # This drops the ZAF ACDIS survey
    left_join(gadm_to_loc_id, by = c("country" = "ihme_lc_id")) %>%
    mutate(admin1 = tolower(admin_1)) %>%
    dplyr::rename(country_long = location_name) %>%
    dplyr::select(nid, country_long, country, source, year, strata,
                  psu, point, shapefile, location_code, sex_id, age_year,
                  int_year, indicator, weight, admin1) %>% data.table()

  micro_data <- rbind(point_data, polygon_data)

  # Reset names to specified indicator, weight column
  return(micro_data)
}



#' crosswalk_reports function ####################################################################################
#'
#' This function takes report data that has HIV prevalence of age ranges that are not the standard 15-49 age range,
#' and maps the recorded prevalence to the standard 15-49 age range by comparing micro data we have form surveys that cover
#' both the standard and unstandard age ranges. This function then downsizes the sample size N to reflect the relative uncertainty
#' inherent in mapping the prevalence from a different age range to 15-49, renaming the sample size the N_eff or effective sample size
#'
#'
#' @param crosswalk_data This is the microdata used to fit a linear model that compares the HIV prevalence of the unstandard age range
#'                       to that of our standard 15-49 age range at the admin1 level. This data is processed using the above crosswalk_prepare() function.
#' @param report_data    Report data that has an unstandard age range (for example, HIV prevalence of 15-54 year olds) which will be mapped to HIV prevalence
#'                       of 15-49 year olds using the crosswalked data
#'
#'
#' @return
#' This function returns the report_data data table, with new columns hiv_old (the old hiv prevalence, named hiv_test), hiv_new (the new prevlance of 15-49)
#' hiv_var (the variance associated with the uncertainty in our estimate of hiv prevalence), and N_eff (the downsampled HIV prevalence to reflect uncertainty)
#'
#' PDF of linear modle in log space fitting microdata hiv prevalence in seperate age ranges ('J/WORK/11_geospatial/10_mbg/hiv/hiv_test/crosswalk_output/fiting_standard_[age_range])
#' PDF of new to old hiv prevalences for the crosswalked report data ('J/WORK/11_geospatial/10_mbg/hiv/hiv_test/crosswalk_output/Adjusted_hiv_test_[age_range])



crosswalk_reports <- function(crosswalk_data, 
                              report_data, 
                              indicator = "hiv_test", 
                              weight = "hiv_weight",
                              filepath_indicator = filepath_indicator) {


  # Load required packages
  library(MASS)
  library(dplyr)
  summarize <- dplyr::summarize
  set.seed(130)

  # Subset to age range
  age_range <- c(unique(report_data$start_age), unique(report_data$end_age))
  lower_range <- min(c(15, age_range[1]))
  upper_range <- max(c(49, age_range[2]))

  message("Age crosswalk", " ",age_range[1] ,"-",age_range[2])

  # Make sure survey covers at least upper and lower bounds and then subset to those bounds
  nids <- 
    crosswalk_data %>%
    group_by(nid) %>%
    summarize(low = min(age_year), high = max(age_year))

  #This is becuase surveys that are lower than 14 and higher than 64 presumably sample all ages
  if (lower_range < 1)  nids <- filter(nids, low < 14)  else nids <- filter(nids, low <= lower_range)
  if (upper_range > 65) nids <- filter(nids, high > 64) else nids <- filter(nids, high >= upper_range)

  # Subset crosswalk data to surveys that cover entire range, subset to the correct year range
  crosswalk_data <- crosswalk_data[nid %in% nids$nid & between(age_year, lower_range, upper_range)]

  # Drop any surveys that are sex-specific if hiv, not for circumcision
  if (indicator == "hiv_test") {
    crosswalk_data <-
      crosswalk_data %>%
      group_by(nid) %>%
      summarize(sexes = max(sex_id, na.rm = T) - min(sex_id, na.rm = T)) %>%
      filter(sexes == 1) %>%
      semi_join(crosswalk_data, ., by = "nid")

    # Drop sex specific surveys that only ask males or females older than 49
    if (age_range[2] > 49) {
      crosswalk_data <-
        crosswalk_data %>%
        filter(age_year > 49) %>%
        group_by(nid) %>%
        summarize(sexes = max(sex_id, na.rm = T) - min(sex_id, na.rm = T)) %>%
        filter(sexes == 1) %>%
        semi_join(crosswalk_data, ., by = "nid")
    }
  }

  # Data_standard has prevalence for 15-49 age range to crosswalk, data_crosswalked has entire age range
  data_standard <- filter(crosswalk_data, between(age_year, 15, 49))
  data_crosswalked <- filter(crosswalk_data, between(age_year, age_range[1], age_range[2]))

  # See how much more crosswalk data is there for 15-49 vs the other age_ranges
  adjustment_factor <- min(c(nrow(data_standard), nrow(data_crosswalked))) / max(c(nrow(data_standard), nrow(data_crosswalked)))

  #data collapsed to admin1 level
  data_standard <- 
    data_standard %>%
    group_by(nid, country, source, year, point, shapefile, admin1) %>%
    summarize(
      prev_standard = weighted.mean(indicator, weight),
      N = sum(weight) ^ 2 / sum(weight ^ 2),
      N_obs = n()
    ) %>%
    ungroup() %>%
    dplyr::select(nid, country, year, point, admin1, prev_standard, N, N_obs)

  data_crosswalked <- 
    data_crosswalked %>%
    group_by(nid, country, source, year, point, shapefile, admin1) %>%
    summarize(
      prev_crosswalked = weighted.mean(indicator, weight),
      N = sum(weight) ^ 2 / sum(weight ^ 2),
      N_obs = n()
    ) %>%
    ungroup() %>%
    dplyr::select(nid, country, year, point, admin1, prev_crosswalked, N, N_obs)

  # In some cases the weight of summing N is zero, and in this small cases N is missing
  data_crosswalked <- data_crosswalked %>% filter(!is.na(N))

  #add country gadm codes, then add region in seperate column to divide admin1 by region
  data_countries <- data_standard$country %>% unique
  data_gadm <- get_adm0_codes(data_countries, shapefile_version = shapefile_version)

  data_standard <-
    data.table(country = data_countries, gadm = data_gadm) %>%
    left_join(data_standard, ., by = "country")

  #Matrix of gadm codes matched with region
  regions <- c("cssa", "essa_sdn", "sssa", "wssa")
  region_gadm <- data.table(matrix(ncol = 2, nrow = 0))
  colnames(region_gadm) <- c('region', 'gadm')
  for (rr in regions) {
    rr_gadm <- get_adm0_codes(rr, shapefile_version = shapefile_version)
    region_gadm <- rbind(region_gadm,
                         data.table(region = rep(rr, length(rr_gadm)),
                                    gadm = rr_gadm))
  }

  region_gadm$region <- as.character(region_gadm$region)

  ## merge on region
  data_standard <- left_join(data_standard, region_gadm, by = "gadm")

  #Merge datasets so we have seperate estimates for hiv prevalence in standard and crosswalking age groups
  data_compare <- merge(data_standard, data_crosswalked, by = c('nid', 'country', 'admin1', 'year', 'point'), all.x = TRUE, all.y=FALSE)
  data_compare <- filter(data_compare, !is.na(prev_crosswalked)) # Missing for same reason as stated above, if N is missing

  #Make point a factor to subset plots
  data_compare$region <- as.factor(data_compare$region)

  #Rename the sample sizes of both sample and crosswalked data
  data_compare <- data_compare %>%
    setnames(c('N.x', 'N.y'), c("sample_standard", "sample_crosswalked"))

  #passing years to four year bins
  data_compare$years <- data_compare$year
  data_compare$years[between(data_compare$years, 1998, 2002)] <- 2000
  data_compare$years[between(data_compare$years, 2003, 2007)] <- 2005
  data_compare$years[between(data_compare$years, 2008, 2012)] <- 2010
  data_compare$years[between(data_compare$years, 2013, 2017)] <- 2015
  data_compare$years <- as.factor(data_compare$years)

  # Emp Logit tranform the prevalence ------------------------------------------------------------------------------------

  epsilon1 <- data_compare %>%
    filter(prev_standard > 0) %>%
    pull(prev_standard) %>%
    min()/2

  epsilon2 <- data_compare %>%
    filter(prev_crosswalked > 0) %>%
    pull(prev_crosswalked) %>%
    min()/2

  epsilon <- min(epsilon1, epsilon2)

  # Add epsilon or subtract epsilon such that we can logit transform
  data_compare <-
    data_compare %>%
    mutate(prev_standard = ifelse(prev_standard == 1, prev_standard - epsilon, prev_standard)) %>%
    mutate(prev_standard = ifelse(prev_standard == 0, prev_standard + epsilon, prev_standard)) %>%
    mutate(prev_crosswalked = ifelse(prev_crosswalked == 1, prev_crosswalked - epsilon, prev_crosswalked)) %>%
    mutate(prev_crosswalked = ifelse(prev_crosswalked == 0, prev_crosswalked + epsilon, prev_crosswalked)) %>%
    mutate(prev_crosswalked = logit(prev_crosswalked)) %>%
    mutate(prev_standard = logit(prev_standard))

  # Grab r squared for later
  r_squared <- summary(lm(prev_standard ~ prev_crosswalked, data = data_compare))$r.squared
  
  # Visual comparison of prevalences in logit space --------------------------------------------------------------------------------
  gg1 <-
    ggplot(data_compare, aes(x = prev_crosswalked, y = prev_standard)) +
    geom_point(aes(color = years, size = sample_standard), alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE, fullrange = T) +
    facet_wrap( ~ region) +
    scale_color_discrete("Year", drop = F) +
    scale_size_continuous("Sample size", range = c(0,9)) +
    # geom_text(x = quantile(data_compare$prev_crosswalked, .01),
    #           y = max(data_compare$prev_standard),
    #           label = lm_eqn(data_compare), parse = TRUE) +
    stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse = TRUE,
      xpos = quantile(data_compare$prev_crosswalked, .01),
      ypos = max(data_compare$prev_standard)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    coord_equal() +
    theme_bw()
  if (indicator == "hiv_test") {
    gg1 <- gg1 +
      labs(x = paste0("Logit HIV prevalence ", age_range[1], '-', age_range[2]),
           y = "Logit HIV prevalence 15-49",
           title = "Relationship between Logit HIV prevalence at different age ranges by year (Admin1)")
  } else {
    gg1 <- gg1 +
      labs(x = paste0("Logit circumcision prevalence ", age_range[1], '-', age_range[2]),
           y = "Logit circumcision prevalence 15-49",
           title = "Relationship between logit circumcision prevalence at different age ranges by year (Admin1)")
  }


  #Crosswalking results and adjusting N effective size and  prevalence ----------------------------------------------------

  # Inspecting report data ------------------------
  setnames(report_data, indicator, "indicator")

  #remove observations where observed prevalence is zero or 1, these values of N are adjusted and then added on in the end so they will not break N adjustment
  if (nrow(report_data[indicator == 0]) != 0) {
    report_zero <-
      report_data %>%
      filter(indicator == 0) %>%
      mutate(N_eff = N * adjustment_factor) %>%
      dplyr::rename(prev_old = indicator) %>%
      mutate(prev_new = prev_old) %>%
      mutate(prev_var = 0)

    #remove when hiv_test is 0 for subsequent analysis
    report_data <- report_data[indicator != 0]
  }

  if (nrow(report_data[indicator == 1]) != 0) {
    report_one <-
      report_data %>%
      filter(indicator == 1) %>%
      mutate(N_eff = N * adjustment_factor) %>%
      dplyr::rename(prev_old = indicator) %>%
      mutate(prev_new = prev_old) %>%
      mutate(prev_var = 1)

    #remove when hiv_test is 0 for subsequent analysis
    report_data <- report_data[indicator != 1]
  }



  # Begin generating uncertainty -------------------------------------------------------------------------------------------------------------
  # Estimates for covariates in linear model, as well as the covariance matrix
  linear <- lm(prev_standard ~ prev_crosswalked, data = data_compare)
  cov <- vcov(linear)
  mean_coef <- linear$coefficients
  sigma <- summary(linear)$sigma

  # Draw 1000 estimates for intercept and slope of linear regression model
  draws <- 1000
  intercept_slope <- mvrnorm(n = draws, mean_coef, cov)

  # Draw 1000 estimates for standard error of linear model
  stand_error <- rnorm(draws, 0, sigma)

  # Convert to data table for ease
  fit <- data.table(Beta_0 = intercept_slope[, 1],
                    Beta_1 = intercept_slope[, 2],
                    Error = stand_error)

  # Generate 1000 draws, will be used for 1000 draws of the HIV prev 15-49
  draw_mat <- matrix(data = NA,
                     nrow = nrow(report_data),
                     ncol = draws)

  # Make sure report is in logit space, and then show uncertainty in estimate
  report_data[, indicator := logit(indicator)]

  for (i in 1:nrow(report_data)) {
    draw_mat[i,] <-
      fit$Beta_0 + fit$Beta_1 * report_data[i, indicator] + fit$Error
  }

  # Return to prevalence space
  draw_mat <- inv.logit(draw_mat)
  report_data[, indicator := inv.logit(indicator)]


  # Propogate uncertainty inherent in binomial model
  for (i in 1:nrow(report_data)) {
    for (j in 1:draws){
      draw_mat[i, j] <- rbinom(1, prob = draw_mat[i, j],
                               size = round(report_data$N[i])) / round(report_data$N[i])
    }
  }

  # New N effective is calculating using the binomial theorem
  report_data <- report_data %>%
    mutate(prev_new = rowMeans(draw_mat)) %>%
    mutate(prev_var = apply(draw_mat, 1, var)) %>%
    mutate(N_eff = ((prev_new * (1 - prev_new)) / prev_var)) %>%
    dplyr::rename(prev_old = indicator)

  # # As hiv prevalence increaes, so does the Variance
  # ggplot(report_data,aes(x = hiv_old, y = hiv_var)) +
  #   geom_point(aes(size =N)) +
  #   geom_smooth(method = "lm")
  #
  # # As N increaes, the variance increases
  # ggplot(report_data, aes(x = N, y = hiv_var)) +
  #   geom_point(aes(color = hiv_old)) +
  #   geom_smooth(method = "lm")
  #
  # # Relationship between N_eff and N
  #Adjustment factor
  # mutate(report, adjust = N_eff / N) %>%
  #   dplyr::select(adjust) %>%
  #   boxplot()

  # Add back on when prevalence was zero or 1
  if (exists("report_zero")){
    report_data <- rbind(report_data, report_zero)
    rm(report_zero)
  }

  if (exists("report_one")){
    report_data <- rbind(report_data, report_one)
    rm(report_one)
  }

  #making plots of results ----------------------------------------------------------------------------------

  # Map linetype to linear fit or equivalence plot
  if (nrow(report_data) > 1) {
  linear_coef <- coefficients(lm(prev_new ~ prev_old, data = report_data))
  report_line <- data.table(slope = c(1, linear_coef[[2]]),
                            intercept = c(0, linear_coef[[1]]),
                            linetype = as.factor(c(1L, 3L)))
  }

  max_data <- max(report_data$prev_new, report_data$prev_old)
  gg2 <-
    ggplot(report_data, aes(x = prev_old, y = prev_new)) +
    geom_point(aes(color = N_eff / N, size = N_eff), alpha = 0.9) +
    stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse = TRUE,
                     xpos = 0.05,
                     ypos = max_data + 0.075) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_size_continuous("N effective") +
    scale_color_distiller("Ratio of N effective\nto sample size",
                          palette = "RdYlBu",
                          direction = 1) +
    theme_bw() +
    coord_fixed(ratio = 1,
                xlim = c(0, max_data + .1),
                ylim = c(0, max_data + .1)) +
    guides(size = guide_legend(order = 1),
           color = guide_colorbar(order = 2))

  if (nrow(report_data) > 1) {

    gg2 <- gg2 +
      geom_abline(data = report_line,
                  aes(slope = slope,
                      intercept = intercept,
                      linetype = linetype),
                  alpha = 0.7, color = "grey") +
      scale_linetype_manual(
        "Line",
        values = c(1, 3),
        breaks = c(1L, 3L),
        labels = c("Equivalence", "Linear\nFit")
      )
  } else {
    gg2 <- gg2 +
      geom_abline(slope = 1, intercept = 0, alpha = 0.7, color = "grey")
  }

  if(indicator == "hiv_test"){
    gg2 <- gg2 +
      labs(y = "HIV prevalence 15-49",
           x = paste0("HIV prevalence ", age_range[1], '-', age_range[2]),
           title = "Adjusted HIV prevalence")
  } else {
    gg2 <- gg2 +
      labs(y = "Male Circumcision Prevalence 15-49",
           x = paste0("Male Circumcision Prevalence ", age_range[1], '-', age_range[2]),
           title = "Adjusted Male Circumcision Prevalence")
  }

  # Save ggplot objects to a pdf
  root <- "<<<< FILEPATH REDACTED >>>>"
  pdf(paste0(root, 'WORK/11_geospatial/10_mbg/hiv/', filepath_indicator,
             '/crosswalk_output/Adjusted_', indicator, '_', age_range[1], '_', age_range[2], '.pdf'),
      width = 16, height = 10)
  print(gg1)
  print(gg2)

  dev.off()

  # add in r_squared to report data
  report_data$r_square = r_squared
  return(report_data)
}

## Utility functions ###############################################################################################
# stat_smooth_func --------------------------------------------------------------------------------------------
# This function is used to add geom text of linear fit onto graphs
stat_smooth_func <- function(mapping = NULL,
                             data = NULL,
                             geom = "smooth",
                             position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto(
  "StatSmooth",
  Stat,

  setup_params = function(data, params) {
    # Figure out what type of smoothing to do: loess for small datasets,
    # gam with a cubic regression basis for large data
    # This is based on the size of the _largest_ group.
    if (identical(params$method, "auto")) {
      max_group <- max(table(data$group))

      if (max_group < 1000) {
        params$method <- "loess"
      } else {
        params$method <- "gam"
        params$formula <-
          y ~ s(x, bs = "cs")
      }
    }
    if (identical(params$method, "gam")) {
      params$method <- mgcv::gam
    }

    params
  },

  compute_group = function(data,
                           scales,
                           method = "auto",
                           formula = y ~ x,
                           se = TRUE,
                           n = 80,
                           span = 0.75,
                           fullrange = FALSE,
                           xseq = NULL,
                           level = 0.95,
                           method.args = list(),
                           na.rm = FALSE,
                           xpos = NULL,
                           ypos = NULL) {
    if (length(unique(data$x)) < 2) {
      # Not enough data to perform fit
      return(data.frame())
    }

    if (is.null(data$weight))
      data$weight <- 1

    if (is.null(xseq)) {
      if (is.integer(data$x)) {
        if (fullrange) {
          xseq <- scales$x$dimension()
        } else {
          xseq <- sort(unique(data$x))
        }
      } else {
        if (fullrange) {
          range <- scales$x$dimension()
        } else {
          range <- range(data$x, na.rm = TRUE)
        }
        xseq <-
          seq(range[1], range[2], length.out = n)
      }
    }
    # Special case span because it's the most commonly used model argument
    if (identical(method, "loess")) {
      method.args$span <- span
    }

    if (is.character(method))
      method <- match.fun(method)

    base.args <-
      list(quote(formula),
           data = quote(data),
           weights = quote(weight))
    model <-
      do.call(method, c(base.args, method.args))

    m = model
    eq <-
      substitute(
        italic(y) == a + b %.% italic(x) * "," ~  ~ italic(r) ^ 2 ~ "=" ~ r2,
        list(
          a = format(coef(m)[1], digits = 3),
          b = format(coef(m)[2], digits = 3),
          r2 = format(summary(m)$r.squared, digits = 3)
        )
      )
    func_string = as.character(as.expression(eq))

    if (is.null(xpos))
      xpos = min(data$x) * 0.9
    if (is.null(ypos))
      ypos = max(data$y) * 0.9
    data.frame(x = xpos, y = ypos, label = func_string)

  },

  required_aes = c("x", "y")
)