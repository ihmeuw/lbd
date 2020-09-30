combine_aggregation <- function(rd = run_date,
                                indic = indicator,
                                ig = indicator_group,
                                ages,
                                regions,
                                holdouts,
                                raked,
                                dir_to_search = NULL,
                                delete_region_files = T,
                                merge_hierarchy_list = F,
                                check_for_dupes = F,
                                measure = 'prevalence',
                                metrics = c('rates', 'counts'),
                                eti = 'none') {
  
  # Combine aggregation objects across region
  
  # Args:
  #   run_date, indicator, indicator_group for this run
  #   ages: single value or vector of ages
  #   regions: vector of regions used
  #   holdouts: vector of holdouts used, e.g. just 0 or c(1,2,3,4,5,0)
  #   raked: vector of raked values, e.g. just T, just F, or c(T,F)
  #   dir_to_search: which directory to search in (defaults to share directory)
  #   delete_region_files: logical. Should we delete the region-specific intermediate files?
  #   merge_hierarchy_list: logical. Do you want to merge the sp_hierarchy_list onto your admin tables?
  
  # Outputs:
  #   rdata files for each combo of age/holdout/raked
  #   each with admin_0, admin_1, admin_2 data table objects & the sp_hierarchy_list object
  #   that maps them to names of admin units
  
  if (is.null(dir_to_search)) {
    dir_to_search <- '<<<< FILEPATH REDACTED >>>>'
  }
  
  message("Combining aggregation results...")
  for (et in eti){
    for(rake in raked) {
      for (holdout in holdouts) {
        for (age in ages) {
          for (metric in metrics) {
            message(paste0("\nWorking on age: ", age, " | holdout: ", holdout, " | raked: ", rake, " | metric: ", metric, ' | etiology: ', et))
            
            # Set up lists
            ad0 <- list()
            ad1 <- list()
            ad2 <- list()
            sp_h <- list()
            
            
            for (reg in regions) {
              message(paste0("  Region: ", reg))
              
              load(paste0(dir_to_search, indic, "_", 
                          ifelse(et == 'none', '', paste0(et, '_')),
                          ifelse(raked, paste0("raked_", measure), "unraked"), 
                          ifelse(metric == "counts", "_c", ""),
                          "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))
              
              if(merge_hierarchy_list == T) {
                # Prepare hierarchy list for adm0
                ad0_list <- subset(sp_hierarchy_list, select = c("ADM0_CODE", "ADM0_NAME", "region")) %>% unique
                
                # Prepare hierarchy list for adm1
                ad1_list <- subset(sp_hierarchy_list,
                                   select = c("ADM0_CODE", "ADM1_CODE", "ADM0_NAME", "ADM1_NAME", "region")) %>%
                  unique
                
                # Merge
                admin_0 <- merge(ad0_list, admin_0, by = "ADM0_CODE", all.y = T)
                admin_1 <- merge(ad1_list, admin_1, by = "ADM1_CODE", all.y = T)
                admin_2 <- merge(sp_hierarchy_list, admin_2, by = "ADM2_CODE", all.y = T)
                rm(ad0_list, ad1_list)
              }
              if(check_for_dupes){
                adms <- get_adm0_codes(reg)
                sp_hier <- get_sp_hierarchy()
                include_ad0 <- sp_hier$ADM0[ADM0_CODE %in% adms, ADM0_CODE]
                include_ad1 <- sp_hier$ADM1[ADM0_CODE %in% adms, ADM1_CODE]
                include_ad2 <- sp_hier$ADM2[ADM0_CODE %in% adms, ADM2_CODE]
                
                ad0[[reg]] <- admin_0[ADM0_CODE %in% include_ad0]
                ad1[[reg]] <- admin_1[ADM1_CODE %in% include_ad1]
                ad2[[reg]] <- admin_2[ADM2_CODE %in% include_ad2]
                sp_h[[reg]] <- sp_hierarchy_list
              } else{
                ad0[[reg]] <- admin_0
                ad1[[reg]] <- admin_1
                ad2[[reg]] <- admin_2
                sp_h[[reg]] <- sp_hierarchy_list
              }
              
              rm(admin_0, admin_1, admin_2, sp_hierarchy_list)
              
              # remove any fixed effects collumns
              if (use_inla_country_fes) {
                if (nchar(reg) != 3) {
                  ad0[[reg]][, c(grep(pattern = 'gaul_code', colnames(admin_0))) := NULL]
                  ad1[[reg]][, c(grep(pattern = 'gaul_code', colnames(admin_1))) := NULL]
                  ad2[[reg]][, c(grep(pattern = 'gaul_code', colnames(admin_2))) := NULL]
                }
              }
              
              # add region names to file
              ad0[[reg]][, region := reg]
              ad1[[reg]][, region := reg]
              ad2[[reg]][, region := reg]
            }
            
            # Get to long format & save
            message("  Combining...")
            admin_0 <- rbindlist(ad0)
            admin_1 <- rbindlist(ad1)
            admin_2 <- rbindlist(ad2)
            sp_hierarchy_list <- rbindlist(sp_h)
            
            message("  Saving combined file...")
            save(admin_0, admin_1, admin_2, sp_hierarchy_list,
                 file = paste0(dir_to_search, indic, "_",
                               ifelse(et == 'none', '', paste0(et, '_')),
                               ifelse(rake, paste0("raked_", measure), "unraked"),
                               ifelse(metric == "counts", "_c", ""),
                               "_admin_draws_eb_bin", age, "_",
                               holdout, ".RData"))
          }
        }
      }
    }
  }
  
  
  if (delete_region_files == T) {
    # Make sure all full files are written
    combos <- expand.grid(ifelse(raked, "raked", "unraked"), ages, holdouts)
    files_to_check <- sapply(1:nrow(combos), function(i) {
      paste0(dir_to_search, indic, "_", combos[i, 1],
             "_admin_draws_eb_bin", combos[i, 2], "_", combos[i,3], ".RData")
    })
    
    if (all(file.exists(files_to_check))) {
      message("All anticipated combined files were created successfully.  Deleting intermediate files...")
      combos <- expand.grid(ifelse(raked, "raked", "unraked"), ages, regions, holdouts)
      files_to_delete <- sapply(1:nrow(combos), function(i) {
        paste0(dir_to_search, indic, "_", combos[i, 1],
               "_admin_draws_eb_bin", combos[i, 2], "_", combos[i,3], "_", combos[i, 4], ".RData")
      })
      unlink(files_to_delete)
    } else {
      warning("Did not delete intermediate files - not all output files created successfully!")
    }
  }
  
  # Finally, delete the "fin" files
  fin_files_to_delete <- list.files(dir_to_search, pattern = "fin_agg_", full.names=T)
  unlink(fin_files_to_delete)
}


## submit_aggregation_script ---------------------------------------------------
# This origianlly was in the post_estimation_functions.R script, but was moved
# to misc_functions.R since all of the other functions which construct qsub
# commands are here and it needed to be updated to take arguments for
# Singularity which requires some helper functions that are also defined in
# misc_functions.R.
#
#' Constructs a qsub string and executes it
#'
#' @param geo_nodes If TRUE, your job will be submitted to the geos (LBD)
#'   cluster, if FALSE, it will be submitted to the prod cluster. Note that if
#'   using the 'proj' argument, make sure to use project name which is valid on
#'   the cluster you are submitting to. [default = FALSE]
#'
#' @param use_c2_nodes If TRUE, your job will be submitted to the C2 nodes on
#'   the prod cluster, if FALSE, the C2 nodes are not specified. Note that if
#'   FALSE, your job may land on a node with much less memory or your node may
#'   still land on a C2 node anyway. If both the 'use_c2_nodes' and 'geo_nodes'
#'   arguments are set to TRUE, then the code will issue a warning and default
#'   to the geos nodes. [default = FALSE]
#'
#' @param proj Can pass in a project name to submit your job under. If default
#'   and the 'geo_nodes' argument is left as its default of 'FALSE', jobs
#'   will be submitted to the prod cluster under the default project
#'   'proj_geospatial'. If default and with 'geos_nodes = TRUE', jobs will be
#'   submitted to the geos (LBD) nodes under the default project
#'   'proj_geo_nodes'. If a project name is passed in for 'proj' the job will
#'   be submitted under that project. Note that this function does not check for
#'   valid project names since these are likely to change often and likely
#'   valid project names are different on each cluster. [default = NULL]
#' 
#' @param queue Queue to be used on the fair cluster.
#'
#' @param run_time Run-time to be used on the fair cluster.
#'
#' @param singularity Launch R from a Singularity image. The default is
#   'default' indicating that you wish to launch a Singularity container from
#'   the default image. You may also provide a string which can be either a complete
#'   path to a Singularity image that is not located at the default image
#'   location, or just the name of the Singularity image that is assumed located
#'   at the default image location. NULL is also accepted, which will launch R
#'   using the default R installation on the geos or prod nodes, but this is
#'   no longer recommended and will likely be deprecated at some point in the
#'   future.
#'
#'   If 'default' is chosen, the default image is defined in the shell script
#'   executed by this R script ('shell_sing.sh') so that no R code need be
#'   updated when the default image is updated. Different versions of a
#'   Singularity image or test versions may be specified by providing the name
#'   or path of the image. Currently, all standard images for LBD are kept at
#'   the default location of sing image.
#'   [default = 'default']
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
#'   [default = NULL]
#' @param modeling_shapefile_version character string specifying which shapefile version was used in modeling
#' @param raking_shapefile_version character string specifying which shapefile version was used in raking
#' @param cores specify number of cores to use, defaults to NULL. If this is provided by the user, it is used to assign resources in get_resources

submit_aggregation_script <- function(indicator,
                                      indicator_group,
                                      proj             = NULL,
                                      run_date,
                                      raked,
                                      pop_measure,
                                      overwrite,
                                      ages,
                                      holdouts,
                                      regions,
                                      corerepo         = core_repo,
                                      log_dir,
                                      geo_nodes        = FALSE,
                                      use_c2_nodes     = FALSE,
                                      queue            = NULL,
                                      run_time         = NULL,
                                      priority         = 0,
                                      slots            = cores,
                                      cores            = 8,
                                      memory           = 20,
                                      singularity      = singularity_version,
                                      singularity_opts = NULL,
                                      modeling_shapefile_version = 'current',
                                      raking_shapefile_version = 'current',
                                      measures,
                                      submit_qsubs     = TRUE) {

  # Takes vectors of regions, holdouts, and ages and submits qsubs for each of these
  # Define project first (necessary to validate node options)
  proj <- get_project(proj, use_geo_nodes=geo_nodes)
  
  # Validate arguments
  validate_singularity_options(singularity, singularity_opts)
  validate_node_option(geo_nodes, use_c2_nodes, proj)
  
  # Create sharedir
  sharedir = get_model_output_dir(indicator_group, indicator, run_date)
  dir.create(sharedir, showWarnings = FALSE)
  output_err = setup_log_location(log_dir, user, indic, ig, rd)
  output_log_dir = output_err[[1]]
  error_log_dir = output_err[[2]]

  # Define remaining job attributes
  
  run_time <- get_run_time(use_geo_nodes=geo_nodes, use_c2_nodes=use_c2_nodes, queue=queue, run_time=run_time)
  queue <- get_queue(use_geo_nodes=geo_nodes, use_c2_nodes=use_c2_nodes, queue=queue, run_time=run_time)
  
  shell <- paste0(corerepo, '/mbg_central/share_scripts/shell_sing.sh')
  sing_image <- get_singularity(image = singularity)
  singularity_str <- qsub_sing_envs("", singularity_opts, sing_image)
  # resources are all the -l qsub arguments
  if(!is.null(cores) & !is.null(slots)) warning("Slots and cores are both specified, cores will be used to assign resources")
  if(is.null(cores)) cores <- slots
  resources <- get_resources(use_geo_nodes=geo_nodes, cores=cores, ram_gb=memory, runtime=run_time)

  code <- path_join(corerepo, 'mbg_central', 'share_scripts', 'aggregate_results.R')

  qsubs_to_make <- expand.grid(regions, holdouts, ages, raked)
  
  aggregation_qsubs <- make_qsubs_aggregation(qsubs_to_make, error_log_dir, output_log_dir, proj, resources, singularity_str, queue, priority, slots, shell, code,
                                              indicator, indicator_group, run_date, pop_measure, overwrite, corerepo, raking_shapefile_version, modeling_shapefile_version)

  if (submit_qsubs) {
      for(qsub in aggregation_qsubs) {
        for (measure in measures){
          system(paste(qsub, measure))
        }
      }
  }
  return(aggregation_qsubs)
}

waitforaggregation <- function(rd = run_date,
                               indic = indicator,
                               ig = indicator_group,
                               ages,
                               regions,
                               holdouts,
                               raked,
                               measure = agg_measures,
                               dir_to_search = NULL) {
  
  # waitformodelstofinish() analog for aggregation
  
  if (is.null(dir_to_search)) {
    dir_to_search <- '<<<< FILEPATH REDACTED >>>>'
  }
  
  lv <- expand.grid(regions, holdouts, ages, raked, measure) %>% as.data.table
  names(lv) <- c("reg", "holdout", "age", "raked", "measure")
  lv[, file := paste0(dir_to_search, "fin_agg_", reg, "_", holdout, "_", age, "_", raked, "_", measure)]
  n_left <- nrow(lv)
  
  message("Waiting for aggregation to finish...")
  
  while(n_left > 0) {
    message(paste0("\n=====================================",
                   "\nCurrent time: ", Sys.time()))
    
    lv[, file_exists := file.exists(file)]
    n_left <- nrow(lv[file_exists == F])
    lv_left <- lv[file_exists==F]
    
    if (n_left > 0) {
      message(paste0("Aggregation jobs remaining: ", n_left))
      for (i in 1:nrow(lv_left)) {
        message(paste0("  ",
                       "Region: ", lv_left[i, reg], " | ",
                       "Measure: ", lv_left[i, measure], " | ",
                       "Holdout: ", lv_left[i, holdout], " | ",
                       "Age: ", lv_left[i, age], " | ",
                       "Raked: ", lv_left[i, raked], " | "))
      }
    }
    
    Sys.sleep(60)
    
  }
}
