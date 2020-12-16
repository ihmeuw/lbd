# HEADER ------------------------------------------------------------------
# Project: Geospatial - general
# Purpose: Diagnostic functions to assess various aspects of your model
# Details: Pairs well with "prep_diagnostic_shiny.R" script
#**************************************************************************


#' @title An empty `ggplot()` theme
#' @description An empty `ggplot()` theme, to be added to a \code{ggplot} object
#'
#' @return theme object for ggplot

theme_empty <- function() {
  theme_empty <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
  return(theme_empty)
}


#' @title Plot stacker betas
#' @description  Create a plot of stacker betas
#'
#' @param model_results table from [whatever]_model_results_table.csv
#' @param stackers vector containing names of stackers
#' @return ggplot object plotting stackers by region

plot_stacker_betas <- function(model_results, stackers, xaxis = "stacker") {
  stacker_beta_table <- subset(model_results, parameter %in% stackers)
  dodge <- position_dodge(width = 0.4)
  if (xaxis == "stacker") {

    gg <- ggplot(data = stacker_beta_table, aes(x = parameter, y = q0.5, color = region)) +
      geom_point(position = dodge) +
      geom_pointrange(aes(ymin=q0.025, ymax=q0.975), position = dodge) +
      labs(x = "Stacker", y = "Beta", color = "Region") +
      scale_color_manual(values = get_color_scheme("carto_discrete")) +
      theme_classic()
  } else if (xaxis == "region") {
    gg <- ggplot(data = stacker_beta_table, aes(x = region, y = q0.5, color = parameter)) +
      geom_point(position = dodge) +
      geom_pointrange(aes(ymin=q0.025, ymax=q0.975), position = dodge) +
      labs(x = "Region", y = "Beta", color = "Stacker") +
      scale_color_manual(values = get_color_scheme("carto_discrete")) +
      theme_classic()
  }
  return(gg)
}


#' @title Plot other parameters
#' @description Create a plot of other parameters from INLA fit
#'
#' @param model_results table from [whatever]_model_results_table.csv
#' @param other_params names of other parameters (not stackers)
#' @return ggplot object plotting parametersby region

plot_other_params <- function(model_results, other_params) {

  other_param_table <- subset(model_results, parameter %in% other_params)
  dodge <- position_dodge(width = 0.4)
  ggplot(data = other_param_table, aes(x = region, y = q0.5, color = region)) +
    geom_point(position = dodge) +
    geom_pointrange(aes(ymin=q0.025, ymax=q0.975), position = dodge) +
    labs(x = "Region", y = "Value", color = "Region") +
    facet_wrap(~parameter, scales = "free_y") +
    scale_color_manual(values = get_color_scheme("carto_discrete")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#' @title Plot Stackers
#' @description Plot maps of stackers, mean raster of covariate, and data informing model
#' @param reg region
#' @param ig indicator group
#' @param ind indicator
#' @param rd run date
#' @param ss vector of stacker names
#' @param yl year list (in vector form, e.g. `c(2000:2015)`)
#' @param zmin minimum value for color scheme (if NULL will calculate from data)
#' @param zmax maximum value for color scheme (if NULL will calculate from data)
#' @param mo_dir model output dir (`FILEPATH` directory, including run date)
#' @param highisbad should high values be colored in red ("bad")? Logical.
#' @param o_dir output directory
#' @param individual_countries should individual countries be graphed as well? Logical.
#' @param shapefile_version string specifies shapefile version to be used
#' @return Writes a series of image fileswith maps of each stacker, mean covariate raster, and
#'         a map of input data for each year-region combination (and year-country if individual_countries = T)
#'         in prespecified folder structure

# plot stackers & mean outcome
plot_stackers <- function(reg,
                          ig = indicator_group,
                          ind = indicator,
                          rd = run_date,
                          ss = stackers,
                          yl = year_list,
                          zmin = NULL, zmax = NULL,
                          mo_dir = modir,
                          highisbad = F,
                          o_dir = out_dir,
                          individual_countries = T,
                          shapefile_version = 'current') {

  # Load master shape for outlines
  # master_shape <- readRDS("FILEPATH")
  master_shape <- readRDS(get_admin_shapefile(admin_level = 0, version = shapefile_version, suffix = ".rds"))
  master_shape <- subset(master_shape, ADM0_CODE %in% get_adm0_codes(reg, shapefile_version = shapefile_version))

  # Set up output dir
  o_dir <- paste0(o_dir, "/stacker_maps/")
  dir.create(o_dir, recursive = T, showWarnings = FALSE)

  # Load stacker objects
  load(paste0("FILEPATH"))
  # subset cov_list to stackers
  stacker_list <- cov_list[which(names(cov_list) %in% ss)]


  # load mean unraked raster for estimates. Try region rasters first; whole raster & crop if not.
  if (file.exists(paste0(mo_dir, ind, "_", reg, "_unraked_mean_raster.tif"))) {
    result_brick <- brick(paste0(mo_dir, ind, "_", reg, "_unraked_mean_raster.tif"))
  } else if (file.exists(paste0(mo_dir, ind, "_", reg, "_mean_raster.tif"))) {
    result_brick <- brick(paste0(mo_dir, ind, "_", reg, "_mean_raster.tif"))
  } else if (file.exists(paste0(mo_dir, ind, "_unraked_mean_raster.tif"))) {
    result_brick <- brick(paste0(mo_dir, ind, "_unraked_mean_raster.tif"))
    result_brick <- crop(result_brick, master_shape)
    result_brick <- mask(result_brick, master_shape)
  } else if (file.exists(paste0(mo_dir, ind, "_mean_raster.tif"))) {
    result_brick <- brick(paste0(mo_dir, ind, "_mean_raster.tif"))
    result_brick <- crop(result_brick, master_shape)
    result_brick <- mask(result_brick, master_shape)
  } else {
    stop("Could not find unraked raster .tif.")
  }

  # add mean raster to stacker list
  stacker_list[[ind]] <- result_brick

  # load input data from csv if present; re-generate if needed
  if (file.exists(paste0(mo_dir, "input_data_bin0_", reg, "_0.csv"))) {
    input_df <- read.csv(paste0(mo_dir, "input_data_bin0_", reg, "_0.csv"),
                         stringsAsFactors=F) %>%
                as.data.table
  } else {
    gaul_list <- get_adm0_codes(reg,
                                shapefile_version = shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                               buffer = 0.4,
                                               shapefile_version = shapefile_version)
    subset_shape   <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
    raster_list    <- build_simple_raster_pop(subset_shape)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]

    input_df <-  load_input_data(indicator = indicator,
                                 removeyemen = TRUE,
                                 yl = yl,
                                 region = reg)
  }

  input_df[, outcome := get(indicator) / N]

  # Define function to save maps
  save_maps <- function(input_df, stacker_list, master_shape,
                        result_brick, zmin, zmax, yl, ind, ig, mo_dir,
                        highisbad, o_dir, ctry = NULL,
                        shapefile_version) {

    if (!is.null(ctry)) {
      input_df <- subset(input_df, country == ctry)
      master_shape <- subset(master_shape, ADM0_CODE == get_adm0_codes(ctry, shapefile_version = shapefile_version))
      stacker_list <- lapply(stacker_list, function(x) {
        x <- suppressMessages(crop(x, extent(master_shape)))
        x <- suppressMessages(mask(x, master_shape))
        extent(x) <- extent(master_shape)
        return(x)
      })
    }

    master_shape_df <- suppressMessages(fortify(master_shape))

    # Cap N at 90% for plot interpretability
    input_df[, cap_N := ifelse(N >= quantile(N, 0.9), quantile(N, 0.9), N)]

    # define mins / maxes from global values across all years in stacker_list & input data
    if (is.null(zmin)) zmin <- min(sapply(stacker_list, function(x) min(minValue(x))),
                                   min(input_df$outcome))
    if (is.null(zmax)) zmax <- max(sapply(stacker_list, function(x) max(maxValue(x))),
                                   max(input_df$outcome))

    # Check # years correct
    if (nlayers(result_brick) != length(yl)) stop("Number of mean raster brick layers does not equal length of year list")

    # rearrange
    for (i in 1:length(yl)) {

      message(paste0("   year ", yl[i], "..."))

      make_gg_map <- function(rbrick, title, i) {
        r_df <- rbrick[[i]] %>% as("SpatialPixelsDataFrame") %>% as.data.frame
        names(r_df) <- c("value", "x", "y")
        gg_result <-ggplot() +
          geom_raster(data = r_df, aes(x = x, y = y, fill = value)) +
          coord_equal() +
          theme_empty() +
          labs(title = title) +
          scale_fill_distiller(palette = "RdYlBu",
                               direction = ifelse(highisbad, -1, 1),
                               limits = c(zmin,zmax)) +
          geom_path(data = master_shape_df,
                    aes(x = long, y = lat, group = group))
        return(gg_result)
      }

      gg_stackers_results <- lapply(1:length(stacker_list), function(n) {
        the_title <- paste0(names(stacker_list)[n], ": ", yl[i])
        the_rbrick <- stacker_list[[n]]
        return(make_gg_map(the_rbrick, the_title, i))
      })

      # Make data plot
      gg_data <- ggplot() +
        geom_point(data = subset(input_df, year == yl[i] & weight < 1),
                   aes(x = longitude, y = latitude, size = cap_N, alpha = weight, color = outcome)) +
        geom_point(data = subset(input_df, year == yl[i] & weight == 1),
                   aes(x = longitude, y = latitude, size = cap_N, color = outcome)) +
        coord_equal() +
        theme_empty() +
        scale_color_distiller(palette = "RdYlBu",
                              direction = ifelse(highisbad, -1, 1),
                              limits = c(zmin,zmax))+
        geom_path(data = master_shape_df,
                  aes(x = long, y = lat, group = group)) +
        #  scale_alpha(range = c(0,1)) +
        scale_size_continuous(limits = c(NA, max(input_df$cap_N))) +
        labs(title = paste0("data: ", yl[i]))

      gg_stackers_results[[length(gg_stackers_results) + 1]] <- gg_data

      # Use first legend only
      the_legend <- g_legend(gg_stackers_results[[1]])
      gg_stackers_results <- lapply(gg_stackers_results, function(x) return(x + theme(legend.position="none")))

      if (is.null(ctry)) {
        reg_dir <- paste0(o_dir, reg, "/")
        dir.create(reg_dir, recursive = T, showWarnings = FALSE)
        fn <-paste0(reg_dir, "stacker_map_", reg, "_", yl[i], ".png")
      } else if (!is.null(ctry)) {
        ctry_dir <- paste0(o_dir, ctry, "/")
        dir.create(ctry_dir, recursive = T, showWarnings = FALSE)
        fn <- paste0(ctry_dir, "stacker_map_", ctry, "_", yl[i], ".png")
      }

      png(filename = fn,
          width = 16, height = 9, units = "in",
          res = 200, pointsize = 10,
          type = "cairo-png")

      multiplot(plotlist = gg_stackers_results,
                cols = ceiling(length(gg_stackers_results) / 2),
                legend = the_legend)

      dev.off()
    }
  }

  # Save map for entire region
  save_maps(input_df, stacker_list, master_shape, result_brick, zmin,
            zmax, yl, ind, ig, mo_dir, highisbad, o_dir, ctry = NULL,
            shapefile_version = shapefile_version)

  # Save maps for individual countries
  if (individual_countries == T) {
    gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
    gaul_list <- data.table(GAUL_CODE = get_adm0_codes(reg, shapefile_version = shapefile_version))
    gaul_list <- merge(gaul_list, gaul_to_loc_id, by = "GAUL_CODE")
     for (c in unique(gaul_list$ihme_lc_id)) {
      if (!(get_adm0_codes(c, shapefile_version = shapefile_version) %in% master_shape$ADM0_CODE)) {
        message(paste0("No shapes in master_shape corresponding to admin 0 code for ", c, " - skipping..."))
      } else {
        message(paste0("Saving stacker maps for country: ", c))
        save_maps(input_df, stacker_list, master_shape, result_brick,
                  zmin, zmax, yl, ind, ig, mo_dir, highisbad, o_dir,
                  ctry = c, shapefile_version = shapefile_version)
      }
    }
  } 
}

#' @title Pulls legend out of a ggplot object for use later
#' @description Pulls legend out of a ggplot object for use later, based off of function originally written by Hadley Wickham
#' @param a.gplot ggplot() object
#' @return grob containing legend
g_legend<-function(a.gplot){
  pdf(NULL) # Workaround for bug in ggplot_gtable causing empty Rplots.pdf to be created
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    legend <- tmp$grobs[[leg]]
  } else {
    legend <- NULL
  }
  graphics.off()
  return(legend)
}


#' @title Plot gam models
#' @description Creates plots of component smooth objects from fitted gam model
#' This function creates plots in a standardized 4x3 format, saved to
#' standardized directories and filenames
#' @param child_model_list Nested list of child model objects (models within regions)
#'                         such that \code{child_model_list[["cssa"]][["gam"]]} returns the
#'                         appropriate gam model object for region cssa
#' @param regions character vector of regions
#' @param o_dir output directory
#' @return writes png files in standardized format to `o_dir\\gam\\[region]\\`

plot_gam_models <- function(child_model_list, regions, o_dir) {

  # Create gam plots in a standardized output format (3x3 grid) for each
  # of the child models (for each region)

  str_match <- stringr::str_match

  for (reg in regions) {

    # Set up directory
    dirname <- paste0(o_dir, "gam/", reg, "/")
    dir.create(dirname, recursive = T, showWarnings = FALSE)

    # Load child model
    gam_child <- child_model_list[[reg]][["gam"]]
    labs <- sapply(gam_child$smooth, function(x) return(x$label))
    n_plots <- length(labs)

    # Create list of which plots go in which pages
    x <- 1:n_plots
    length(x) <- suppressWarnings(prod(dim(matrix(x, ncol = 12, byrow = T))))
    x <- matrix(x, ncol = 12, byrow = T)

    for (i in 1:nrow(x)) {

      # Set up each page (row in layout matrix) & plot

      terms <- x[i,]
      terms <- terms[!is.na(terms)]

      filename <- paste0(dirname, "gam_plot_", reg, "_page_", i, ".png")

      png(filename=filename,
          type="cairo",
          units="in",
          width=12,
          height=8,
          pointsize=14,
          res=200)

      par(mfrow = c(3,4))
      for (z in terms) {
        plot(gam_child, select = z)
      }

      dev.off()
    }
  }
}


#' @title Multiple ggplots
#' @description Function to plot multiple ggplot objects together
#' @source Adapted from \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#' @param ... Can pass in ggplot objects in \code{...} or in plotlist
#' @param plotlist Can pass in ggplot objects as a list with this argument instead of in \code{...}
#' @param cols Number of columns in layout
#' @param layout Matrix specifying the layout; i.e. \code{matrix(c(1,2,3,3), nrow = 2, byrow = T)}.
#'               If \code{layout} is specified, then \code{cols} is ignored
#' @param legend A legend object.  If legend is passed, then this will add an extra cell at the
#'               end of the grid layout and insert the legend there (good, for instance, if you
#'               have common legends for all of your plots and only want to show it once).
#' @return Prints a gridded output of your ggplot objects to the active graphical device
#' @note gg_stackers_results is a list of ggplot objects

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, legend = NULL) {

  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)
  if (!is.null(legend)) numPlots <- numPlots + 1 # add legend

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols),
                     byrow = T)
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:(numPlots-1)) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }

    # Print the legend if present
    if (!is.null(legend)) {
      matchidx <- as.data.frame(which(layout == numPlots, arr.ind = T))
      legend$vp <- viewport(layout.pos.row = matchidx$row,
                            layout.pos.col = matchidx$col)
      grid.draw(legend)
    }
  }
}


#' @title Function to submit a qsub to run the model_diagnostics.R script
#' @description \code{make_model_diagnostics} creates a qsub string to run the '
#'   model_diagnostics.R script on the cluster and does the system call to actually
#'   submit the job
#' @param code Name of script, with relative path if desired.
#' @param code_path Full path to R script. Overrides \code{code} and \code{script_dir}
#' @param cores Number of threads. Default: 5
#' @param memory RAM to be reserved, in GBs
#' @param geo_nodes If TRUE, your job will be submitted to the geos (LBD)
#'   cluster, if FALSE, it will be submitted to the prod cluster. Note that if
#'   using the 'proj' argument, make sure to use project name which is valid on
#'   the cluster you are submitting to. [default = FALSE]
#' @param use_c2_nodes If TRUE, your job will be submitted to the C2 nodes on
#'   the prod cluster, if FALSE, the C2 nodes are not specified. Note that if
#'   FALSE, your job may land on a node with much less memory or your node may
#'   still land on a C2 node anyway. If both the 'use_c2_nodes' and 'geo_nodes'
#'   arguments are set to TRUE, then the code will issue a warning and default
#'   to the geos nodes. [default = FALSE]
#' @param proj Can pass in a project name to submit your job under. If default
#'   and the 'geo_nodes' argument is left as its default of 'FALSE', jobs
#'   will be submitted to the prod cluster under the default project
#'   'proj_geospatial'. If default and with 'geos_nodes = TRUE', jobs will be
#'   submitted to the geos (LBD) nodes under the default project
#'   'proj_geo_nodes'. If a project name is passed in for 'proj' the job will
#'   be submitted under that project. Note that this function does not check for
#'   valid project names since these are likely to change often and likely
#'   valid project names are different on each cluster. [default = NULL]
#' @param queue Queue to be used on the fair cluster.
#' @param run_time Run-time to be used on the fair cluster.
#' @param priority Job priority that can be deprioritized if needed, and can only be used for values in [-1023,0]. Default = 0.
#' This value will get bounded to 0 or -1023 if the user supplies a value outside those bounds.
#' @param singularity Instead of using the default R installation on the geos
#'   or prod nodes, launch R from a Singularity image. This arg currently takes
#'   three options: the default is NULL, indicating not to launch a Singularity
#'   container, 'default' if you wish to launch a Singularity container from the
#'   default image, or you can provide a string which can be either a complete
#'   path to a Singularity image that is not located at the default image
#'   location, or just the name of the Singularity image that is assumed located
#'   at the default image location.  If 'default' is chosen, the default image
#'   is defined in the shell script executed by this R script ('shell_sing.sh')
#'   so that no R code need be updated when the default image is updated.
#'   Different versions of a Singularity image or test versions may be specified
#'   by providing the name or path of the image. Currently, all standard images
#'   for LBD are kept at the default location of "FILEPATH".
#'   [default = NULL]
#' @param singularity_opts pass in a named list of environmental variables.
#'   \code{qsub_sing_envs} will check that the names of the list members passed
#'   in match the environmental variables that the shell_sing.sh script knows
#'   about: 'SET_OMP_THREADS' and/or 'SET_MKL_THREADS'. Passing in other
#'   environmental names in the list will result in an error. If this is left
#'   as 'NULL' and a Singularity image is used, SET_OMP_THREADS and
#'   SET_MKL_THREADS will remain unset and the shell_sing.sh script will use
#'   the default setting of SET_OMP_THREADS=1 and SET_MKL_THREADS={max_threads}
#'   (see shell_sing.sh comments). For example SET_OMP_THREADS=1 and
#'   SET_MKL_THREADS=4 can be achieved by passing in
#'     \code{envs = list(SET_OMP_THREADS=1, SET_MKL_THREADS=4)}
#'   [default = list(SET_OMP_THREADS = cores, SET_MKL_THREADS = cores)]
#' @seealso This function uses the following functions found in
#'   mbg_central/misc_functions.R:
#'   \code{\link{get_singularity}}
#'   \code{\link{qsub_sing_envs}}
#' @export
#'
make_model_diagnostics <- function(user = 'USERNAME'
                                   code_path = NULL,
                                   cores = 5,
                                   memory = 10,
                                   proj = NULL,
                                   ig = indicator_group,
                                   corerepo = core_repo,
                                   indic = indicator,
                                   rd = run_date,
                                   log_location = "sgeoutput",
                                   code = "model_diagnostics",
                                   script_dir = "FILEPATH",
                                   keepimage = FALSE,
                                   shell = "r_shell.sh",
                                   geo_nodes = FALSE,
                                   use_c2_nodes = FALSE,
                                   queue = NULL,
                                   run_time = NULL,
                                   priority = 0,
                                   singularity = singularity_version,
                                   singularity_opts = list(SET_OMP_THREADS = cores, SET_MKL_THREADS = cores)) {
  # Define project first (necessary to validate node options)
  proj <- get_project(proj, use_geo_nodes = geo_nodes)

  # Validate arguments
  validate_singularity_options(singularity, singularity_opts)
  validate_node_option(geo_nodes, use_c2_nodes, proj)

  temp_dir <- path_join(get_model_output_dir(ig, indic, rd), "temp_post_est")
  dir.create(temp_dir, showWarnings = FALSE)

  # Determine where stdout and stderr files will go
  output_err <- setup_log_location(log_location, user, indic, ig, rd)
  output_log_dir <- output_err[[1]]
  error_log_dir <- output_err[[2]]

  # Since we no longer rely on `setwd()`'s, we need to construct a sensible
  # "script_dir". If someone wants to use a special script, we assume it is
  # somewhere in their corerepo here:
  script_dir <- path_join(corerepo, script_dir)
  code <- path_join(script_dir, paste0(code, ".R"))

  # If code_path is not NULL, then override `code`
  if(!is.null(code_path)) {
    code <- code_path
  }

  # Define remaining job attributes
  job_name <- paste0("job_dx_", indic)
  run_time <- get_run_time(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  queue <- get_queue(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  shell <- paste0("FILEPATH")
  sing_image <- get_singularity(image = singularity)

  # resources are all the -l qsub arguments
  resources <- get_resources(use_geo_nodes = geo_nodes, cores = cores, ram_gb = memory, runtime = run_time)

  qsub <- generate_qsub_command(
    # qsub-specific arguments
    stderr_log = error_log_dir,
    stdout_log = output_log_dir,
    project = proj,
    resources = resources,
    job_name = job_name,
    singularity_str = qsub_sing_envs("", singularity_opts, sing_image),
    cores = cores,
    queue = queue,
    priority = priority,
    # Command to qsub
    shell, code, indic, ig, rd, cores, corerepo
  )

  # make the qsub call
  system(qsub)
}

#' @title extract_hyperpar_prior_posteriors
#' @description extract hyperparameter pior and posterior marginals from inla
#' @param res inla result object
#' @param transform_prec transform precisions to std.deviations
#' @return data.table with hyper parameter prior and posterior marginals
#' @rdname extract_hyperpar_prior_posteriors
#' @import dplyr
#' @import tibble
#' @export

extract_hyperpar_prior_posteriors <- function(res, transform_prec=TRUE){
  all.hyper <- inla_all_hyper_postprocess(res$all.hyper)
  hyper <- res$marginals.hyperpar
  nhyper <- length(hyper)
  # loop through all hyperpriors
  hyperDF <- bind_rows(lapply(1:nhyper, function(i){
    hh <- hyper[[i]]
    label <- INLA:::inla.nameunfix(names(hyper)[i])
    # Extract the posterior valuse into m
    m <- as_tibble(INLA::inla.smarginal(hh)) %>%
      mutate(Distribution="Posterior", label=label)
    id <- unlist(strsplit(attr(hyper[[i]], "hyperid"), "\\|"))
    prior <- INLA:::inla.extract.prior(tolower(id[2]), id[1], all.hyper)
    # Extract the prior values, if else statements required mostly for SPDE
    if(grepl("the Gaussian observations", label)){
      label <- gsub("the Gaussian observations", "Normal obs", label)
    }

    if(prior$prior == "none"){
      df <- m
    }

    else if(prior$prior == "pcmatern"){
      # extract lambda1, lambda2 and d from inla
      l1 <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[1]
      l2 <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[2]
      d <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[3]

      # calculate pc priors as done in 'Constructing Priors that
      # Penalize the Complexity of Gaussian Random Fields' by Fuglstad
      # et al, 2019.  The joint PC prior on top of pg 448 (RHS)
      # implies sigma and rho are independent with sigma ~ expo(lambda2)
      calc_pc_prior_ranges <- function(ranges){
        d/2 * l1 * ranges^(-d/2 - 1) * exp(-l1 * ranges^(-d/2))
      }

      calc_pc_prior_sigmas <- function(sigmas){
        l2 * exp(- l2 * sigmas)
      }

      m2 <- as_tibble(INLA::inla.smarginal(hyper[[i+1]]))
      offset_ <- m2$x[2] - m2$x[1]
      df <- bind_rows(
        m,
        m %>%
          mutate(y=calc_pc_prior_ranges(x)) %>%
          mutate(Distribution="Prior", label=label),
        tibble(x=seq(offset_, max(m2$x), by=offset_)) %>%
          mutate(y=calc_pc_prior_sigmas(x)) %>%
          mutate(
            Distribution="Prior",
            label=gsub("Range", "Stdev", label))
      )
    }

    # id 23001 gives us the spde sigma parameter untransformed if we
    # also have a normal prior we backtransfrom for interpretability
    else if((prior$prior == "normal") & (id[1] == "23001")){
      
      kapp <- exp(res$summary.hyperpar[i+1,"0.5quant"])
      fx <- function(x) (1 / (4 * pi * kapp^2 * exp(x) ^ 2))^.5
      label <- gsub("Theta1", "Stdev", label)
      malt <- INLA::inla.smarginal(hh) %>%
        {INLA::inla.tmarginal(fx, .)} %>%
        as_tibble() %>%
        mutate(Distribution="Posterior", label=label)
      r_ <- range(m$x)
      df <- bind_rows(
        malt,
        INLA:::inla.get.prior.xy(
          tolower(id[2]), id[1], all.hyper, range = r_) %>%
          {INLA::inla.tmarginal(fx, .)} %>%
          as_tibble() %>%
          mutate(Distribution="Prior", label=label) %>%
          filter(x < max(malt$x))
      )
    }

    # id 23002 gives us the spde kappa parameter untransformed if we
    # also have a normal prior we backtransfrom for interpretability
    else if((prior$prior == "normal") & (id[1] == "23002")){
      label <- gsub("Theta2", "Range", label)
      malt <- INLA::inla.smarginal(hh) %>%
        {INLA::inla.tmarginal(function(x) sqrt(8) / exp(x), .)} %>%
        as_tibble() %>%
        mutate(Distribution="Posterior", label=label)
      r_ <- range(m$x)
      df <- bind_rows(
        malt,
        INLA:::inla.get.prior.xy(
          tolower(id[2]), id[1], all.hyper, range = r_) %>%
          {INLA::inla.tmarginal(function(x) sqrt(8) / exp(x), .)} %>%
          as_tibble() %>%
          mutate(Distribution="Prior", label=label) %>%
          filter(x < max(malt$x))
      )
    }

    # transformation for precisions to be in stdev space for interpretability
    else if(grepl("Precision", label) & transform_prec){
      label <- gsub("Precision", "Stdev", label)

      malt <- INLA::inla.smarginal(hh) %>%
        {INLA::inla.tmarginal(function(x) x^-.5, .)} %>%
        as_tibble() %>%
        mutate(Distribution="Posterior", label=label)
      r_ <- range(m$x)
      df <- bind_rows(
        malt,
        tibble(x=malt$x) %>%
          mutate(y=dnorm(x, mean=prior$param[1], sd=prior$param[2]^-.5)) %>%
          mutate(Distribution="Prior", label=label) %>%
          filter(x < max(malt$x)))
    }

    else{
      # portion for any prior except pcmatern prior
      r_ <- range(m$x)
      df <- bind_rows(
        m,
        INLA:::inla.get.prior.xy(
          tolower(id[2]), id[1], all.hyper, range = r_) %>%
          as_tibble() %>%
          mutate(Distribution="Prior", label=label)
      )
    }

    filter(df, !is.na(x) & !is.na(y)) %>%
      mutate(prior=as.character(prior$prior))
  }))
  return(data.table::as.data.table(hyperDF))
}

#' @title get_data_estimate_df
#' @description Constructs the modeling dir
#' @param indicator_group Indicator Group
#' @param indicator Indicator String
#' @param run_date The run date
#' @param region a region of input data to pull, may be multiple
#' @param yrs vector of ordered integres for years modeled
#' @param cache_results logical, save results to disk for easy access later
#' @param use_cache logical, if cache exists use it to load data
#' @param save_cols character, additional columns to save for comparison
#' @param age_inputs_separated logical, are age inputs grouped together ex: u5m
#' @param age integer, age groups to pull may be multiple
#' @param sex integer, future proof input for sex categories
#' @param shapefile_version character, string indicating linktable version
#' @param holdout integer, holdout number to extract, may be multiple
#' @return character vector of length 1 for the modeling directory
#' @rdname get_data_estimate_df
#' @import data.table
#' @import sp
#' @export

get_data_estimate_df <- function(
  indicator_group,
  indicator,
  run_date,
  region,
  yrs                         = 2000:2015,
  cache_results               = TRUE,
  use_cache                   = TRUE,
  save_cols                   = c("data_type", "nid"),
  age_inputs_separated        = TRUE,
  age                         = 0,
  sex                         = 0,
  shapefile_version           = "current",
  holdout                     = 0){


  if((length(region) > 1) | (length(age) > 1) | (length(holdout) > 1)){
    demogDF <- expand.grid(ages=age, regions=region, holdouts=holdout)
    df <- as.data.table(bind_rows(lapply(1:nrow(demogDF), function(i){
      get_data_estimate_df(
        indicator_group             = indicator_group,
        indicator                   = indicator,
        run_date                    = run_date,
        region                      = as.character(demogDF$regions[i]),
        yrs                         = yrs,
        cache_results               = cache_results,
        use_cache                   = use_cache,
        save_cols                   = save_cols,
        age_inputs_separated        = age_inputs_separated,
        age                         = demogDF$ages[i],
        sex                         = sex,
        shapefile_version           = shapefile_version,
        holdout                     = demogDF$holdouts[i]
      )
    })))

    return(df)
  }

  mod_dir <- get_model_output_dir(indicator_group, indicator, run_date)
  viz_dir <- paste0(mod_dir, "/vizualization_outputs/")

  if(!dir.exists(viz_dir)){
    dir.create(viz_dir)
  }

  save_file <- sprintf(
    paste0(viz_dir, "data_estimate_bin%s_%s_%s.csv"),
    age, region, holdout)

  if(use_cache & file.exists(save_file)){
    df.r <- fread(save_file)
    return(df.r)
  }

  save_cols_add <- c(
    c("ihme_loc_id", "year", "N", "weight", "age", "sex"),
    save_cols, indicator,
    c("holdout", "latitude", "longitude", "ADM_CODE", "region"))


  nperiod <- length(yrs)

  mod_dir <- get_model_output_dir(indicator_group, indicator, run_date)
  adm_to_loc_id <- get_location_code_mapping(
    shapefile_version=shapefile_version)

  df <- data.table()

  message("LOADING IN-SAMPLE RAW DATA")

  base_path <- "%s/%s_trainingdata_bin%s_%s_%s%s"

  f_ <- sprintf(base_path, mod_dir, indicator, age, region, holdout, ".csv")
  if(!file.exists(f_)){
    f_ <- sprintf(base_path, mod_dir, indicator, age, region, holdout, "")
  }

  if(!age_inputs_separated){
    f_ <- sprintf(base_path, mod_dir, indicator, 0, region, holdout, ".csv")
    if(!file.exists(f_)){
      f_ <- sprintf(base_path, mod_dir, indicator, 0, region, holdout, "")
    }
  }

  tmp <- fread(f_)
  tmp$holdout <- holdout

  tmp$region <- region
  if(!("age" %in% names(tmp))){
    tmp$age <- 0
  }
  if(!("sex" %in% names(tmp))){
    tmp$sex <- 0
  }

  # subset to the sex and age that is specified in the function call
  aa <- age
  ss <- sex
  tmp <- tmp[(sex == ss) & (age == aa),]


  # removing countries not in region (BORDER DATA)
  adm0_codes <- get_adm0_codes(region, shapefile_version=shapefile_version)
  tmp$ihme_lc_id <- as.character(tmp$country)
  tmp <- merge(tmp, adm_to_loc_id, by="ihme_lc_id",all.x=T)
  tmp <- tmp[ADM_CODE %in% adm0_codes]

  # rbind it on to the full dataset
  save_cols_present <- save_cols_add[save_cols_add %in% names(tmp)]
  tmp <- tmp[, save_cols_present, with=FALSE]

  if(holdout!=0){

    f_2 <- sprintf(base_path, mod_dir, indicator, age, region, ".csv")
    if(!file.exists(f_)){
      f_2 <- sprintf(base_path, mod_dir, indicator, age, region, "")
    }

    if(!age_inputs_separated){
      f_2 <- sprintf(base_path, mod_dir, indicator, 0, region, ".csv")
      if(!file.exists(f_)){
        f_2 <- sprintf(base_path, mod_dir, indicator, 0, region, "")
      }
    }

    tmp2 <- fread(f_2)

    tmp2$region <- region
    if(!("age" %in% names(tmp2))){
      tmp2$age <- 0
    }
    if(!("sex" %in% names(tmp2))){
      tmp2$sex <- 0
    }

    # subset to the sex and age that is specified in the function call
    tmp2 <- tmp2[(sex == ss) & (age == aa),]

    # removing countries not in region (BORDER DATA)
    tmp2 <- tmp2[, ihme_lc_id := as.character(country)]
    tmp2 <- merge(tmp2, adm_to_loc_id, by="ihme_lc_id",all.x=T)
    tmp2 <- tmp2[ADM_CODE %in% adm0_codes]

    # rbind it on to the full dataset
    save_cols_present <- save_cols_add[save_cols_add %in% names(tmp2)]
    tmp2 <- tmp2[, save_cols_present, with=FALSE]

    tmp <- as.data.table(left_join(tmp2, tmp))
    tmp$in.sample <- !is.na(tmp$holdout)
    tmp$holdout <- holdout
  }

  df  <- rbind(df, tmp)

  message("Identify ad1 and ad2 membership")

  # load admin2 shapefile (also contains admin1)
  admin2_shapefile <- rgdal::readOGR(get_admin_shapefile(
    admin_level=2,
    version=shapefile_version))

  for (v in grep("CODE", names(admin2_shapefile@data))){
    admin2_shapefile@data[[v]] <- as.numeric(as.character(
      admin2_shapefile@data[[v]]))
  }

  # identify the admin2 each unique point belongs to
  ucoords  <- unique(data.table(
    longitude=df$longitude,
    latitude=df$latitude))[,ucid:=1:.N] # unique for sppedups on sp::over
  locs     <- SpatialPointsDataFrame(
    cbind(ucoords$longitude, ucoords$latitude),
    data        = ucoords,
    proj4string = CRS(proj4string(admin2_shapefile)))
  adm.df   <- sp::over(locs, admin2_shapefile)

  # for those that don't fall inside a polygon, assign them to the nearest
  # polygon (this tends to happen on coastlines) do this by country to speed
  # things up and to make sure the point ends up at least in the correct admin0
  ii          <- which(is.na(adm.df[,1]))
  temp_shape  <- admin2_shapefile[
    admin2_shapefile@data$ADM0_CODE %in% unique(df$ADM_CODE),] # only region
  distmat     <- rgeos::gDistance(locs[ii,], temp_shape, byid = TRUE)
  jj          <- apply(distmat, 2, which.min)
  adm.df[ii,] <- temp_shape@data[jj,]
  rm("distmat")

  # merge admins back on to df based on coordinates
  adm.df     <- cbind(adm.df, ucoords)
  adcols <- c(
    paste0("ADM", 0:2, "_CODE"), "ADM1_NAME", "ADM2_NAME", "longitude", "latitude")

  adm.df.tmp <- data.table(adm.df)[,adcols,with=FALSE]
  df <- merge(df,adm.df.tmp,by=c("latitude","longitude"),all.x=TRUE)

  # copy admin information into df
  setnames(df,paste0("ADM",0:2,"_CODE"),paste0("ad",0:2))

  message('DONE GETTING RAW DATA, NOW GETTING DRAWS')

  message(paste('...Region:', region))

  # load the simple raster template
  message("......load simple raster template")
  gaul_list <- get_adm0_codes(region, shapefile_version = shapefile_version)
  simple_polygon_list <- load_simple_polygon(
    gaul_list = gaul_list, buffer = 0.4, subset_only = TRUE,
    shapefile_version = shapefile_version
  )
  subset_shape <- simple_polygon_list[[1]]
  raster_list <- build_simple_raster_pop(
    subset_shape, link_table = shapefile_version)
  template <- raster_list[["simple_raster"]]

  # subset data to this region
  rr <- region
  df.r <- df[region == rr, ]
  loc.r <- as.matrix(df.r[, list(longitude, latitude)])

  # 'nudge' coordinates to make them fall inside the nearest cell non-NA cell
  # in template if they are in an NA cell (this can happen near coastlines)
  check <- raster::extract(template, loc.r)
  miss <- which(is.na(check))
  for (ii in miss) {
    dist <- replace(raster::distanceFromPoints(
      template, loc.r[ii, ]), is.na(template), NA)
    loc.r[ii, ] <- raster::xyFromCell(template, raster::which.min(dist)[1])
  }

  # make a raster of cell_pred row ids
  id_raster <- seegMBG::insertRaster(
    template, matrix(1:(length(seegMBG::cellIdx(template)) * nperiod), ncol = nperiod))

  # get cell_pred row ids for each data point
  for (yy in 1:nperiod) {
    this_year <- which(df.r$year == yrs[yy])
    if (length(this_year) == 0) next
    df.r[
      this_year,
      cell_pred_id := raster::extract(
        id_raster[[yy]], matrix(loc.r[this_year, ], ncol=2))]
  }

  df.r$run_date <- run_date


  # load cell pred objects
  message(paste("......load cell preds for holdout", holdout))
  load(paste0(mod_dir, sprintf(
    "/%s_cell_draws_eb_bin%i_%s_%i.RData", indicator, age, region, holdout)))
  if(!("cell_pred" %in% ls())){
    cell_pred <- cptmp
  }

  df.r[, paste0("draw", 1:ncol(cell_pred)) := as.data.table(cell_pred[cell_pred_id, ])]

  if(cache_results){
    write.csv(df.r, save_file, row.names=FALSE)
  }

  return(df.r)
}

#' @title inla_all_hyper_postprocess
#' @description Postprocess multivariate normal hyperparameters to get thier
#' normal marginals
#' @param all.hyper hyper parameters from inla fit ex: inla.result$all.hyper
#' @return updated hyperparameter object
#' @rdname inla_all_hyper_postprocess

inla_all_hyper_postprocess <- function(all.hyper){
  ## postprocess all.hyper, by converting and replacing prior = 'mvnorm' into
  ## p marginals. this is for the spde-models

  len.n <- function(param, max.dim = 10000){
    len <- function(n) {
      return (n + n^2)
    }

    len.target <- length(param)
    for(n in 1:max.dim) {
      if (len(n) == len.target) {
        return(n)
      }
    }
    stop(paste("length(param) is wrong:", len.target))
  }

  get.mvnorm.marginals = function(param){
    n <- len.n(param)
    mu <- param[1:n]
    Q <- matrix(param[-(1:n)], n, n)
    Sigma <- solve((Q + t(Q))/2.)
    return(list(mean = mu, prec = 1/diag(Sigma)))
  }

  for (i in seq_along(all.hyper$random)) {
    for(j in seq_along(all.hyper$random[[i]]$hyper)) {

      if (all.hyper$random[[i]]$hyper[[j]]$prior == "mvnorm") {
        ## replace this one, and the p-following ones, with its marginals
        m <- get.mvnorm.marginals(all.hyper$random[[i]]$hyper[[j]]$param)
        for(k in 1:length(m$mean)) {
          kk <- j + k - 1
          all.hyper$random[[i]]$hyper[[kk]]$prior <- "normal"
          all.hyper$random[[i]]$hyper[[kk]]$param <- c(m$mean[k], m$prec[k])
        }
      }
    }
  }

  return (all.hyper)
}

#' @title load_admin_draws
#' @description Load admin level indicator draws
#' @param indicator_group character Indicator Group
#' @param indicator character Indicator String
#' @param run_date character The run date
#' @param regions character string of regions to read NULL will return all regions
#' @param raked bool Get raked results
#' @param age integer, age groups to pull may be multiple
#' @param sex integer, future proof input for sex categories
#' @param level administrative level to return data for
#' @param holdout integer, holdout number to extract, may be multiple
#' @return data.frame like object with draw level indicator data
#' @details Aggregates input data usually created in launch script into a single
#' file to read later in the plotting process.
#' @rdname load_admin_draws
#' @export

load_admin_draws <- function(
  indicator_group, indicator, run_date, region,
  raked=FALSE,
  age=0,
  sex=0,
  level=0,
  holdout=0){

  ## model outputs
  modir <- get_model_output_dir(
    indicator_group, indicator, run_date)

  data_list <- list()
  e <- new.env()
  admCodeDF <- sf::read_sf(get_admin_shapefile(2)) %>%
    as_tibble() %>%
    select(c(
      paste0("ADM", 0:level, "_CODE"),
      paste0("ADM", 0:level, "_NAME"))) %>%
    unique

  ## load in all admin results
  for(r in region){
    for(a in age){
      for(h in holdout){
        paste0(
          modir, "/", indicator, ifelse(raked, "_raked_", "_unraked_"),
          "admin_draws_eb_bin", a,"_", r, "_", h, ".RData") %>%
        load(envir = e, verbose = F)

        data_list <- c(
          data_list,
          get(paste0("admin_", level), envir = e) %>%
          mutate(age=a, holdout=h, region=r, sex=sex, run_date=run_date) %>%
          mutate(raked=raked) %>%
          list())

        rm(list=ls(envir = e), envir=e)
        gc(full=TRUE)
      }
    }
  }

  ## combine and format
  df <- as.data.table(left_join(
    bind_rows(data_list), admCodeDF, by=paste0("ADM", level, "_CODE")))
  old_draw_names <- grep("V[0-9]*", names(df), value=T)
  new_draw_names <- sub("V", "draw", old_draw_names)
  setnames(df, old_draw_names, new_draw_names)
  
  return(df)
}

#' @title read_inla_model
#' @description Reads an inla model from disk
#' @param indicator_group Indicator Group
#' @param indicator Indicator String
#' @param run_date The run date
#' @param region region that was run
#' @param age age group run
#' @param sex sex run (not yet implemented)
#' @param holdout holdout number
#' @param custom a custom file path for home directory
#' @return completed INLA model
#' @rdname read_inla_model
#' @export
read_inla_model <- function(
  indicator_group, indicator, run_date, region, age=0, sex=0, holdout=0, custom=NULL){

  f_ <- paste0(
    get_model_output_dir(indicator_group, indicator, run_date),
    sprintf(
      "/inla_model_fit_pre_preds_%s_holdout_%s_agebin_%s.RDS",
      region, holdout, age))

  return(readRDS(f_))
}

#' @title save_plot_list
#' @description Saves plot list object in the same format as the list construction
#' @param plotList plot list with ggplot2 plot objects
#' @param indicator_group Indicator Group
#' @param indicator Indicator String
#' @param run_date The run date
#' @param custom a custom file path for the save directory directory
#' @param save_image_format whether to save image flat files
#' @param ext image extension name
#' @param ... additional arguments passed to ggsave
#' @rdname save_plot_list
#' @export

save_plot_list <- function(
  plotList, indicator_group=NULL, indicator=NULL, run_date=NULL, custom=NULL,
  save_image_format=FALSE, ext="eps", ...){
  mdir <- get_model_output_dir(indicator_group, indicator, run_date)
  if(!is.null(custom)){
    mdir <- custom
  }

  viz_dir <- paste0(mdir, "/vizualization_outputs/")

  if(!dir.exists(viz_dir)){
    dir.create(viz_dir, showWarnings = FALSE)
  }

  img_dir <- paste0(viz_dir, "figures/")
  saveRDS(plotList, file=paste0(viz_dir, "plotList.RDS"))

  if(save_image_format){
    save_plots(plotList, img_dir, ext, ...)
  }
}

#' @title save_plots
#' @description Saves plot list object in recursive fashion; INTERNAL ONLY
#' @param obj recursive plot list; MUST BE NAMED
#' @param dir current save directory
#' @param ext image extension name
#' @param ... additional arguments passed to ggsave
#' @rdname save_plot

save_plots <- function(obj, dir, ext, ...){
  dir.create(dir, showWarnings = FALSE)
  for(n in names(obj)){
    if(inherits(obj[[n]], "gg")){
      f_ <- paste0(dir, gsub("/", " ", n), ".", ext)
      ggsave(f_, obj[[n]], ...)
    }
    else{
      new_dir <- paste0(dir, gsub("/", " ", n), "/")
      save_plots(obj[[n]], new_dir, ext, ...)
    }
  }
}

#' @title summarize_draws
#' @description Summarize draw level data from outputs
#' @param df data frame from load_admin_draws or get_data_estimate_df
#' @param .width the confidence interval widths to calculate
#' @param .metric summary metric to calculate, usually mean or median
#' @return data.table with summary statistics calculated from draws
#' @details Aggregates input data usually created in launch script into a single
#' file to read later in the plotting process.
#' @rdname summarize_draws
#' @export

summarize_draws <- function(df, .width=c(.95), .metric=mean){
  draw_cols <- grep("draw[0-9]*", names(df), value=T)
  draws <- as.matrix(as.data.table(df)[,draw_cols,with=FALSE])

  v_ <- apply(draws, 1, .metric, na.rm=TRUE)

  data.table(bind_rows(lapply(.width, function(w){
    df %>%
      select(-!!draw_cols) %>%
      mutate(
        .value = v_,
        .lwr   = rcpp_row_quantile(draws, .5 - w/2),
        .upr   = rcpp_row_quantile(draws, .5 + w/2),
        .width = w
      )})))
}
