## estimate_parallel_model_mem() ----------------------------------------------
#'
#' @title Estimate parallel_model.R memory requirements
#'
#' @description Based off of profiled model runs from 2018, AOZ
#'   performed some nonlinear regression to estimate the relationship
#'   between peak RAM usage and number of pixels, number of years, and
#'   the number of posterior samples in prediction
#'
#' @details This function is meant as a STARTING POINT to determine how
#'   much memory a job requires. For each region (or pixels in a
#'   region) submitted to the function, the function will estimate
#'   suggested memory in GB that your run might need. It gives 4
#'   options per region that are increasingly 'safe' in the sense that
#'   they are increasingly conservative.
#'
#' The regression that is used for estimation uses ACE/AVAS
#' methodology to perform the nonlinear regression. We assume that
#' peak memory usage occurs in genrating the large cell_pred object
#' and thus is a function of non-NA pixels in the region, the number
#' of years in the prediction, and the number of posterior samples.
#'
#' In the returned suggestions, the first one is the best fit linear
#' relationship - and it goes through the middle of the dataset
#' meaning that it often would not estimate large enough memory
#' requirements. Suggestions 2, 3, and 4 add 1, 2, and 3 times
#' (respectively) the residual standard deviation from the nonlinear
#' fit. Basically this just adds some buffers of safety to the
#' regression. The outputs (printed and the output list) also store
#' the number of times that level of suggestion would have resulted in
#' the model using more memory than the prediction suggests. As you
#' choose a safer level to start, you are less likely to have you job
#' killed, but you eat up more resources... We leave it to the users
#' to determine how much safety vs how many jobs they would initially
#' like to run.
#'
#' REMEMBER, this is only a starting point. Once you have runs
#' complete, look at the results from the run (qacct -j <jobid>) to
#' check what the actual mem max use was and then adjust your requests
#' accordingly.
#'
#' ** ## PARAMETERS ## **
#'
#' @param regions string or string vector naming regions that you need
#'   RAM estimates for. Either regions or npixels must be given.
#' @param npixels number of pixels in a region. if you already know
#'   how many pixels are in a region, enter it directly to save the
#'   function the time required to load and count the number of pixels
#'   in a region
#' @param nyears number of years in model output. Integer
#' @param ndraws number of posterior draws the run will make. Integer
#' @param core_repo path to core_repo
#' @param shapefile_version string specifying version of shapefile to
#'   use. This is needed if regions are given instead of npixels
#' @param print.tables logical. if true, the function will print the
#'   RAM results tables in addition to returning a named list of the
#'   tables
#'
#' @return A named list (names are entries in regions if supplied,
#'   otherwise names are ingteger value corresponding to the order of
#'   entries in npixels). Each list entry is a 2-column matrix
#'   containing a column of RAM suggestions and an estiamte of the
#'   change (in percantage points) that the RAM suggestion may not be
#'   large enough.
#'
#' @examples
#' estimate_parallel_model_mem(regions=c('essa', 'name'), nyears=17, ndraws=500)
#' $`Region: essa`
#'     SUGGESTED RAM (GB) CHANCE OF EXCEEDENCE (%)
#' 1:           308.3816               42.0287253
#' 2:           424.5765               15.7450628
#' 3:           584.5523                3.3752244
#' 4:           804.8052                0.9335727
#'
#' $`Region: name`
#'     SUGGESTED RAM (GB) CHANCE OF EXCEEDENCE (%)
#' 1:           311.0521               42.0287253
#' 2:           428.2532               15.7450628
#' 3:           589.6143                3.3752244
#' 4:           811.7745                0.9335727

estimate_parallel_model_mem <- function(regions = NULL,
                                        npixels = NULL,
                                        nyears,
                                        ndraws,
                                        core_repo = "<<<< FILEPATH REDACTED >>>>",
                                        shapefile_version = 'current',
                                        print.tables = FALSE){

  ## basic checks
  if(is.null(npixels) & is.null(regions)){
    stop('You must supply either a region name or the number pixels in the region to proceed')
  }

  if(!is.null(npixels) & !is.null(regions)){
    message('You have supplied both region names and pixel counts. To save time we will proceed by using the supplied pixel counts')
    regions <- NULL
  }

  if(length(nyears) != 1 | length(ndraws) != 1){
    warning('This function was only tested using integer values for ndraws and nyears. You appear to have entered in vectors for at least one of these arguments and we do not promise that this function will work as expected with your inputs!')
  }

  ## determine number of regions
  nreg <- ifelse(is.null(npixels), length(regions), length(npixels))

  ## make a dataset full of relevant predictors
  train.dat <- data.table(matrix(0.0, ncol = 3, nrow = nreg))
  colnames(train.dat) <- c('region', 'pixels', 'samples')
  if(is.null(regions)){
    train.dat[, region := 1:nreg]
  }else{
    train.dat[, region := regions]
  }
  train.dat[, samples := rep(ndraws, nreg)]

  if(is.null(regions)){
    train.dat[, pixels := npixels * nyears]
  }else{

    ## determine pixel counts if only region names were give
    simple_raster_list <- simple_polygon_list <- NULL
    missing.regs <- NULL
    for(ii in 1:nreg){

      reg <- regions[ii]

      message(sprintf('Getting pixel count from region: %s', reg))

      adm0_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)

      if(length(adm0_list) == 0){

        missing.regs <- c(missing.regs, reg)
        message(sprintf('\n\n--could not find adm0 list for region: %s\n\n', reg))

      }else{

        polygon_list <- suppressMessages(load_simple_polygon(gaul_list = adm0_list,
                                                             buffer = 1, tolerance = 0.4,
                                                             shapefile_version = shapefile_version))
        subset_shape       <- polygon_list[[1]]
        raster_list        <- suppressMessages(build_simple_raster_pop(subset_shape))
        simple_raster      <- raster_list[['simple_raster']] ## get pixels from this

        ## get pixel count
        train.dat[ii, pixels := ncell(simple_raster) * nyears]

        ## save things b/c it's slow
        simple_raster_list[[reg]]  <- raster_list
        simple_polygon_list[[reg]] <- simple_polygon_list
      }
    }
  }

  ## load in the ACE/AVAS model fits that have been estimated
  data(slot_estimation_fits)

  ## make predictions

  ## now we can predict RAM use in the test data and compare against requested RAM
  train.dat[, ram.pred0sd := exp(f$log_ram_gb(predict(ace.r,train.dat) + 0 * res.sd))] ## no sd buffer
  train.dat[, ram.pred1sd := exp(f$log_ram_gb(predict(ace.r,train.dat) + 1 * res.sd))] ## 1 sd buffer
  train.dat[, ram.pred2sd := exp(f$log_ram_gb(predict(ace.r,train.dat) + 2 * res.sd))] ## 2 sd buffer
  train.dat[, ram.pred3sd := exp(f$log_ram_gb(predict(ace.r,train.dat) + 3 * res.sd))] ## 3 sd buffer

  ## now we print out the results

  ram.list <- NULL
  for(rr in 1:nreg){

    dt <- data.table(RAM = t(train.dat[rr, .( ram.pred0sd,  ram.pred1sd,  ram.pred2sd,  ram.pred3sd)]) / 1000,
                     EXCEEDENCE = percent.over)
    colnames(dt) <- c('SUGGESTED RAM (GB)', 'CHANCE OF EXCEEDENCE (%)')

    if(print.tables){
      message(sprintf('FOR REGION: %s', train.dat[rr, region]))
      message('Here are 4 starting suggestions for RAM use and the associated chance \n (from the training data we used) that you will need more RAM than the suggested amount')
      print(dt)
    }

    ram.list[[paste0("Region: ", train.dat[rr, region])]] <- copy(round(dt, 1))
  }

  return(ram.list)

}

