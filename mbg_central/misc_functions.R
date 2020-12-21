
pullupstream<-function(){
  system('cd <<<< FILEPATH REDACTED >>>>\ngit pull origin master')
}

# miscellaneous functions
lstype<-function(type='closure'){
    inlist<-ls(.GlobalEnv)
    if (type=='function') type <-'closure'
    typelist<-sapply(sapply(inlist,get),typeof)
    return(names(typelist[typelist==type]))
}
# to listen in on parallel processes in R: http://stackoverflow.com/questions/10903787/how-can-i-print-when-using-dopar
Log <- function(text, ...) {
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  cat(msg)
  write.socket(log.socket, msg)
}


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
                                metrics = c('rates', 'counts')) {
  
  # Combine aggregation objects across region
  
  # Args:
  #   run_date, indicator, indicator_group for this run
  #   ages: single value or vector of ages
  #   regions: vector of regions used
  #   holdouts: vector of holdouts used, e.g. just 0 or c(1,2,3,4,5,0)
  #   raked: vector of raked values, e.g. just T, just F, or c(T,F)
  #   dir_to_search: which directory to search in
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
  
  for(rake in raked) {
    for (holdout in holdouts) {
      for (age in ages) {
        for (metric in metrics) {
          message(paste0("\nWorking on age: ", age, " | holdout: ", holdout, " | raked: ", rake, " | metric: ", metric))
          
          # Set up lists
          ad0 <- list()
          ad1 <- list()
          ad2 <- list()
          sp_h <- list()
          
          
          for (reg in regions) {
            message(paste0("  Region: ", reg))
            
            load(paste0(dir_to_search, indic, "_", 
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
                             ifelse(rake, paste0("raked_", measure), "unraked"),
                             ifelse(metric == "counts", "_c", ""),
                             "_admin_draws_eb_bin", age, "_",
                             holdout, ".RData"))
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

## get_singularity ------------------------------------------------------------
#' Which Singularity image to use
#'
#' \code{get_singularity} determines which Singularity image to use. The
#' default is the 'default' keyword. Image names without the full path to the
#' file defined are assumed to exist as the default location:
#'    <<<< FILEPATH REDACTED >>>>
#' In either case it will test to make sure the file exists and exit if it does
#' not. The default image is hardcoded into the shell script used to launch
#' Singularity containers:
#'   <<<< FILEPATH REDACTED >>>>
#'
#' @param image A string that defines which Singularity image to launch
#'   [default = 'default']. If the 'default' keyword is passed in or left blank,
#'   the default keyword will be returned. Either the full path to the image
#'   may be provided or only the Singularity image name. In the latter case,
#'   the image is assumed to live in the default image location:
#'   <<<< FILEPATH REDACTED >>>>.
#'
#' @return When image = 'default', 'default' is returned. When this keyword is
#'   passed to the shell_sing.sh script through `qsub` it will use the default
#'   Singularity image hardcoded into it. Otherwise, if the function is
#'   successful at verifying that the Singularity image file specified exists,
#'   it will return the full path to that image.
#'
#' @seealso This function is used by:
#'   \code{\link{parallelize}}
#'   \code{\link{make_qsub}}
#'   \code{\link{make_qsub_share}}
#'   \code{\link{make_qsub_postest}}
#'   \code{\link{submit_aggregation_script}}
#'
get_singularity <- function(image = 'default') {
  if(image == 'default') {        # use default image
    sing_image <- 'default'
  } else if(grepl('/', image)) {  # user supplied path to image
    sing_image <- image
  } else {                        # image at default location
    sing_image <- paste0('<<<< FILEPATH REDACTED >>>>', image)
  }
  # If something other than the default image is being used, let's make sure
  # the image file actually exists:
  if(!sing_image == 'default' & !file.exists(sing_image)) {
    stop(paste0("Could not locate Singularity image: ", sing_image))
  }
  return(sing_image)
}

## qsub_sing_envs -------------------------------------------------------------
#' Adds environmental variables to a qsub string
#'
#' \code{qsub_sing_envs} assumes that a qsub string is being built to launch a
#' Singularity container. It always adds in the '-v sing_image=sing_image' as
#' expected by lbd_core/mbg_central/shell_sing.sh script that ultimately
#' launches the container. Optionally, users may want to pass the additional
#' environmental variables 'SET_OMP_THREADS' and/or 'SET_MKL_THREADS' to
#' shell_sing.sh. If one or both of those are passed into \code{qsub_sing_envs}
#' it will add those environmental variables and their values as additional
#' '-v' flags in the construction of the qsub command.
#'
#' @param qsub_stub A string containing the initial qsub string to which
#'   environmental variables will be concatenated to in the form of
#'   '-v ENV=VALUE'.
#'
#' @param envs This should be a named list of environmental variables.
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
#'
#' @param image The keyword (e.g. 'default') or path to the Singularity image.
#'   This should have been defined by \code{get_singularity} so likely
#'   \code{get_singularity} should have been run on the 'singularity' argument
#'   (in \code{make_qsub_share} or \code{parallelize} for example) before this
#'   function is run.
#'
#' @return Returns a string with at least '-v sing_image=image' and possibly
#'   other environmental variables values if they were passed into the
#'   'singularity_opts' argument of functions like \code{make_qsub_share}.
#'
#' @seealso The function \code{\link{get_singularity}} should likely be run
#'   before this function is run. This function is used by:
#'     \code{\link{parallelize}}
#'     \code{\link{make_qsub}}
#'     \code{\link{make_qsub_share}}
#'     \code{\link{make_qsub_postest}}
#'     \code{\link{submit_aggregation_script}}
#'
qsub_sing_envs <- function(qsub_stub, envs, image) {
  # we always want to add this to our qsub to launch a Singularity container
  qsub_envs <- c(sing_image=image)
  # valid variables in the list passed into 'envs' argument
  valid_env_vars <- c("SET_OMP_THREADS", "SET_MKL_THREADS")
  if(!is.null(envs)){
    # verify first that what is passed in is a list since this is all this
    # function knows how to deal with
    if(!is.list(envs)) {
      message("Expected a list for 'envs' argument")
      stop("Exiting!")
    }
    # if invalid variables were passed in, give a message and exit
    else if(!all(names(envs) %in% valid_env_vars)) {
      message("The following variables were unexpected for the 'envs' argument:")
      for(env in names(envs)[!names(envs) %in% valid_env_vars]) {
        message(paste0("  ", env))
      }
      stop("Exiting!")
    } else qsub_envs <- c(qsub_envs, envs)
  }
  # Concatenate all of the variables to the qsub string with the '-v'
  qsub_stub <- paste0(
                 qsub_stub,
                 paste0(
                   " -v ", names(qsub_envs), "=", qsub_envs, collapse = ""))
  return(qsub_stub)
}

# make_qsub_share --------------------------------------------------------------
#'
#' @title Make qsub string (share)
#'
#' @description Constructs a qsub string and returns it
#'
#' @param user Username
#' 
#' @param code Name of script, with relative path if desired.
#' 
#' @param code_path Full path to R script. Override \code{code}
#' 
#' @param cores Number of threads
#' 
#' @param memory RAM to be reserved, in GBs
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
#' @param ig Indicator Group
#' 
#' @param indic Indicator
#' 
#' @param reg Region
#' 
#' @param age Age
#' 
#' @param rd Run date
#' 
#' @param log_location Location of logs
#' 
#' @param addl_job_name Additional name appended to end of job
#'
#' @param saveimage Save image of prerun image?
#' 
#' @param test Run test job?
#' 
#' @param holdout Holdout
#' 
#' @param corerepo Path to core repo
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

#' @param queue Queue to be used on the fair cluster.
#'
#' @param run_time Run-time to be used on the fair cluster.
#'
#' @param priority Job priority that can be deprioritized if needed, and can only be used for values in [-1023,0]. Default = 0.
#' This value will get bounded to 0 or -1023 if the user supplies a value outside those bounds.
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
#'   the default location of <<<< FILEPATH REDACTED >>>>.
#'   [default = 'default']
#' 
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
#'
#' @return Returns a qsub string
#' @export
#'
make_qsub_share <- function(user = Sys.info()["user"],
                            code = NULL,
                            code_path = NULL,
                            cores = slots,
                            memory = 100,
                            proj = NULL,
                            ig = indicator_group,
                            indic = indicator,
                            reg = "test",
                            age = 0,
                            rd = run_date,
                            log_location = "sharedir",
                            addl_job_name = "",
                            saveimage = FALSE,
                            test = FALSE,
                            holdout = 0,
                            corerepo = core_repo,
                            geo_nodes = FALSE,
                            use_c2_nodes = FALSE,
                            queue = NULL,
                            run_time = NULL,
                            priority = 0,
                            singularity = singularity_version,
                            singularity_opts = NULL) {
  
  # save an image
  if (saveimage == TRUE) save.image(pre_run_image_path(ig, indic, rd, age, reg, holdout))
  
  # Define project first (necessary to validate node options)
  proj <- get_project(proj, use_geo_nodes = geo_nodes)
  
  # Validate arguments
  validate_singularity_options(singularity, singularity_opts)
  validate_node_option(geo_nodes, use_c2_nodes, proj)
  
  # Create sharedir (TODO is this necessary?)
  sharedir <- get_model_output_dir(ig, indic, rd)
  dir.create(sharedir, showWarnings = FALSE)
  
  # Determine where stdout and stderr files will go
  output_err <- setup_log_location(log_location, user, indic, ig, rd)
  output_log_dir <- output_err[[1]]
  error_log_dir <- output_err[[2]]
  
  # Define remaining job attributes
  job_name <- paste0("job_", addl_job_name, "_", reg, "_", age, "_", holdout)
  run_time <- get_run_time(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  queue <- get_queue(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  shell <- paste0(corerepo, "<<<< FILEPATH REDACTED >>>>/shell_sing.sh")
  sing_image <- get_singularity(image = singularity)
  
  # resources are all the -l qsub arguments
  resources <- get_resources(use_geo_nodes = geo_nodes, cores = cores, ram_gb = memory, runtime = run_time)
  
  # set code to shared parallel if its NULL
  if (is.null(code)) {
    code <- sprintf("<<<< FILEPATH REDACTED >>>>/parallel_model.R", corerepo)
  } else {
    code <- sprintf("%s/%s.R", corerepo, code)
  }
  
  # If code_path is not NULL, then override `code`
  if(!is.null(code_path)) {
    code <- code_path
  }
  
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
    shell, code, reg, age, rd, as.numeric(test), holdout, indic, ig, "fin"
  )
  
  return(qsub)
}


# make_qsub_postest -------------------------------------------------------------
#'
#' @title Make qsub string for post estimation
#'
#' @description Constructs a qsub string for the post estimation script and returns it
#'
#' @param user Username
#'
#' @param code Name of script, with relative path if desired.
#'
#' @param code_path Full path to R script. Overrides \code{code} and \code{script_dir}
#'
#' @param cores Number of threads
#'
#' @param memory RAM to be reserved, in GBs
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
#' @param ig Indicator Group
#'
#' @param indic Indicator
#'
#' @param stratum "Strata" , usually region
#'
#' @param age Age
#'
#' @param rd Run date
#'
#' @param log_location Location of logs
#'
#' @param script_dir Location of post-estimation script
#'
#' @param addl_job_name Additional name appended to end of job
#'
#' @param saveimage Save image of prerun image?
#'
#' @param test Run test job?
#'
#' @param holdout Holdout
#'
#' @param corerepo Path to core repo
#'
#' @param modeling_shapefile_version Modeling shapefile version, defaults to "current"
#'
#' @param raking_shapefile_version Raking shapefile version, defaults to "current",
#'
#' @param subnat_raking Do subnational raking?
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

#' @param queue Queue to be used on the fair cluster.
#'
#' @param run_time Run-time to be used on the fair cluster.
#'
#' @param priority Job priority that can be deprioritized if needed, and can only be used for values in [-1023,0]. Default = 0.
#' This value will get bounded to 0 or -1023 if the user supplies a value outside those bounds.
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
#'   the default location of <<<< FILEPATH REDACTED >>>>.
#'   [default = 'default']
#'
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
#'
#' @return Returns a qsub string
#'
#' @export
make_qsub_postest <- function(user = Sys.info()["user"],
                              code,
                              code_path = NULL,
                              cores = slots,
                              memory = 100,
                              proj = NULL,
                              ig = indicator_group,
                              indic = indicator,
                              stratum = "test",
                              addl_job_name = "",
                              rd = run_date,
                              log_location = "sgeoutput",
                              script_dir = NULL,
                              keepimage = FALSE,
                              corerepo = core_repo,
                              modeling_shapefile_version = "current",
                              raking_shapefile_version = "current",
                              subnat_raking = TRUE,
                              geo_nodes = FALSE,
                              use_c2_nodes = FALSE,
                              queue = NULL,
                              run_time = NULL,
                              priority = 0,
                              singularity = singularity_version,
                              singularity_opts = NULL) {
  
  # Create test_post_est dir within model dir.
  temp_dir <- path_join(get_model_output_dir(ig, indic, rd), "test_post_est")
  dir.create(temp_dir, showWarnings = F)
  
  # Define project first (necessary to validate node options)
  proj <- get_project(proj, use_geo_nodes = geo_nodes)
  
  # Validate arguments
  validate_singularity_options(singularity, singularity_opts)
  validate_node_option(geo_nodes, use_c2_nodes, proj)
  
  # Create sharedir (TODO is this necessary?)
  sharedir <- get_model_output_dir(ig, indic, rd)
  dir.create(sharedir, showWarnings = FALSE)
  
  # Determine where stdout and stderr files will go
  output_err <- setup_log_location(log_location, user, indic, ig, rd)
  output_log_dir <- output_err[[1]]
  error_log_dir <- output_err[[2]]
  
  # Define remaining job attributes
  job_name <- paste("job_pe", indic, stratum, sep = "_")
  job_name <- paste0(job_name, addl_job_name)
  run_time <- get_run_time(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  queue <- get_queue(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  shell <- paste0(corerepo, "<<<< FILEPATH REDACTED >>>>/shell_sing.sh")
  sing_image <- get_singularity(image = singularity)
  
  # resources are all the -l qsub arguments
  resources <- get_resources(use_geo_nodes = geo_nodes, cores = cores, ram_gb = memory, runtime = run_time)
  
  # code is short-hand for a file in the script_dir directory; script_dir defaults to share_scripts
  if (is.null(script_dir)) script_dir <- path_join(corerepo, "mbg_central", "share_scripts")
  code <- path_join(script_dir, paste0(code, ".R"))
  
  # If code_path is not NULL, then override `code`
  if (!is.null(code_path)) {
    code <- code_path
  }
  
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
    shell, code, stratum, rd, indic, ig, as.character(geo_nodes),
    modeling_shapefile_version, raking_shapefile_version, subnat_raking
  )
  
  return(qsub)
}



# use my condSim update until forked repo is pulled by NG
#source('../seegMBG/R/gis_functions.R')
cellIdx <- function (x) which(!is.na(getValues(x[[1]])))
insertRaster <- function (raster, new_vals, idx = NULL) {


  # calculate cell index if not provided
  if (is.null(idx)) idx <- cellIdx(raster)

  # check the index makes superficial sense
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))

  # create results raster
  n <- ncol(new_vals)
  raster_new <- raster::brick(replicate(n,
                                        raster[[1]],
                                        simplify = FALSE))
  names(raster_new) <- colnames(new_vals)

  # update the values
  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }

  return (raster_new)

}
condSim <- function (vals, weights = NULL, group = NULL, fun = NULL, ...) {
  # given a matrix of pixel-level prevalence samples `prev`
  # where each rows are pixels and columns are draws, a vector
  # of corresponding pixel populations `pop`, and an optional pixel
  # grouping factor `group`, return draws for the total deaths in each
  # group, or overall if groups are not specified

  # get dimensions of vals
  ncell <- nrow(vals)
  ndraw <- ncol(vals)

  # capture function as a string
  fun_string <- deparse(substitute(fun))

  # check fun accepts a

  # check dimensions of weights and group, set to 1 if not specified
  if (is.null(weights)) {
    weights <- rep(1, ncell)
  } else {
    if (length(weights) != ncell) {
      stop (sprintf('number of elements in weights (%i) not equal to number of cells in vals (%i)',
                    length(weights),
                    ncell))
    }
  }

  if (is.null(group)) {
    group <- rep(1, length(weights))
  } else {
    if (length(group) != ncell) {
      stop (sprintf('number of elements in group (%i) not equal to number of cells in vals (%i)',
                    length(group),
                    ncell))
    }
  }

  # otherwise, get the levels in group and create a matrix of results
  levels <- unique(na.omit(group))
  nlevel <- length(levels)

  ans <- matrix(NA,
                ncol = ndraw,
                nrow = nlevel)
  rownames(ans) <- levels

  # loop through levels in group, getting the results
  for (lvl in 1:nlevel) {

    # get an index o pixels in the level
    idx <- which(group == levels[lvl])

    # by default, calculate a weighted sum
    if (is.null(fun)) {

      # get draws and add to results
      # exception for if area has 1 cell, transpose matrix so it conforms (RB)
      if(all(dim(t(vals[idx, ]))==c(1,ndraw))){
        ans[lvl, ] <- weights[idx] %*% t(vals[idx, ])
      } else {
        ans[lvl, ] <- weights[idx] %*% vals[idx, ]
      }

    } else {

      # otherwise, apply function to each column
      ans[lvl, ] <- apply(vals[idx, ], 2, fun, weights = weights[idx], ...)

    }

  }

  # if only one level, make this a vector
  if (nlevel == 1) ans <- as.vector(ans)

  # return result
  return (ans)

}

defaultOptions <- function (resolution = 5,
                            location = loc(),
                            cores = switch(loc(), oxford = 60, seattle = 60),
                            spacetime = TRUE,
                            start = Sys.time()) {

  # built in check for if they selected more than one country

  # if these options are not already set, set them at these default values

  # list all visible objects in this environment
  object_names <- ls()

  # get them in a named list
  objects <- lapply(object_names,
                    get,
                    envir = environment())
  names(objects) <- object_names

  # empty vector of the options to keep
  keep <- c()

  # loop through identifying those that don't already exist
  for (i in 1:length(objects)) {

    # if doesn't exist yet...
    if (is.null(options(names(objects)[i])[[1]])) {

      # add it to the list
      keep <- c(keep, i)

    }
  }

  # keep the new options
  objects <- objects[keep]

  if (length(objects) > 0) {
    # notify the user
    message('\nThe following options were not defined and have been set to their defaults:\n')
    for(i in 1:length(objects)) {
      message(sprintf('    %s = %s',
                      names(objects)[i],
                      paste(objects[[i]], collapse = ', ')))
    }
    message('\n')
  }

  # add to options
  options(objects)
}


# for cleaning clutter
"rmlike" <- function(...) {
  names <- sapply(
    match.call(expand.dots = FALSE)$..., as.character)
  names = paste(names,collapse="|")
  Vars <- ls(1)
  r <- Vars[grep(paste("^(",names,").*",sep=""),Vars)]
  rm(list=r,pos=1)
}



reportTime <- function () {

  # get start time stored in options
  start <- options()$start

  # if it exists
  if (!is.null(start)) {

    # get elapsed time
    now <- Sys.time()
    diff <- difftime(now, start)

    # format nicely and report
    diff_string <- capture.output(print(diff, digits = 3))
    diff_string <- gsub('Time difference of', 'Time elapsed:', diff_string)
    message (diff_string)
  }
}


tidyInstall <- function (packages, versions = NULL, lib = .libPaths()[1]) {
  # given a vector of package names and optional vector of version numbers
  # install them from CRAN if they aren't already loaded

  npkg <- length(packages)

  # check the arguments
  stopifnot(inherits(packages, 'character'))

  if (!is.null(versions)) {
    stopifnot(length(versions) == npkg)
    stopifnot(inherits(versions, 'character'))
  }

  # run a loop if multiple required
  if (npkg > 1) {

    for (i in 1:npkg) {

      if (is.null(versions))
        versions_tmp <- NULL
      else
        versions_tmp <- versions[i]

      tidyInstall(packages[i], versions_tmp, lib)
    }

  } else {
    # otherwise just for one

    # first make sure versions and devtools are installed
    for (pkg in c('versions', 'devtools')) {
      if (!(pkg %in% installed.packages(lib.loc = lib)))
        install.packages(pkg, lib = lib)
    }

    # check whether the package is installed
    installed <- packages %in% rownames(installed.packages(lib.loc = lib))

    if (installed && is.null(versions)) {

      # if it is installed the user doesn't care about the version exit
      return (invisible())

    } else {

      # otherwise, do some installing...

      # if it's installed, check the version
      if (installed) {

        installed_version <- versions::installed.versions(packages, lib = lib)

        # if that's the required version, exit
        if (versions == installed_version) {
          return (invisible())
        } else {
          # otherwise, install the correct version
          devtools::with_lib(new = lib,
                             devtools::install_version(packages, versions))
        }

      } else {
        # if it isn't installed ...

        if (!is.null(versions)) {
          # if a version is required, install it and exit
          devtools::with_lib(new = lib,
                             devtools::install_version(packages, versions))
          return (invisible())
        } else {
          # otherwise get the most recent version and exit
          install.packages(packages, lib = lib)
          return (invisible())

        }
      }
    }
  }
}

subtractYears <- function (Date, years) {
  # given a vector of dates in Date format,
  # and a number of years (positive numeric, can be non-integer)
  # subtract the required number of years and return

  stopifnot(years >= 0)
  stopifnot(inherits(Date, 'Date'))

  n_Date <- length(Date)
  n_years <- length(years)

  if (n_years != n_Date) {
    if (n_years == 1) {
      years <- rep(years, n_Date)
    } else {
      stop('Date and years have different lengths')
    }
  }

  Date <- as.POSIXlt(Date)

  month <- round(12 * years %% 1)
  year <- floor(years)

  Date$year <- Date$year - year
  Date$mon <- Date$mon - month

  Date <- as.Date(Date)

  return (Date)

}

toYear <- function(Date) format(Date, '%Y')

targetYear <- function (period, period_end, width = 60) {
  # get the target year, based on the period, period size (in months) and end
  # date of the final period

  # convert date to cmc
  period_end <- Date2cmc(period_end)

  # get number of months preceeding
  months <- width / 2 + width * (period - 1)

  # subtract
  cmc <- period_end - months + 1

  # format as a year
  ans <- toYear(cmc2Date(cmc))

  return (ans)

}

# convert strings or numerics to a sequential numeric index
# pass any number of vectors of the same length, get a numeric ID for the
# unique ones
idx <- function (...) {
  list <- list(...)

  # get their lengths
  lengths <- lapply(list, length)
  if (!isTRUE(do.call(all.equal, lengths))) {
    stop('vectors do not appear to have the same length')
  }
  # combine them
  x <- do.call(paste, list)
  #get numeric index
  match(x, unique(x))
}


# given the admin level, a GAUL code and the admin1 and 2 shapefiles,
# return an SPDF with the relevant area
getPoly <- function (level, code, admin1, admin2) {

  # get admin level
  if (level == 1) {
    admin <- admin1
  } else {
    admin <- admin2
  }

  # find the reight area
  idx <- match(code, admin$GAUL_CODE)

  # if it's valid
  if (length(idx) == 1) {
    ans <- admin[idx, ]
  } else {
    # handle errors
    warning (paste0("something's wrong on row ", i))
    ans <- NULL
  }

  # return result
  return (ans)
}

# get predictions from predictive INLA glm models
predGLM <- function (result,
                     intercept = '(Intercept)',
                     fixed_continuous = NULL,
                     fixed_group = NULL,
                     random_group = 'cluster_id') {

  # starting means and variances
  pred_mean <- pred_var <- 0

  # add intercept terms
  if (!is.null(intercept)) {
    pred_mean <- pred_mean + result$summary.fixed[intercept, 'mean']
    pred_var <- pred_var + result$summary.fixed[intercept, 'sd'] ^ 2
  }

  # add continuous fixed effects terms
  if (!is.null(fixed_continuous)) {
    for (name in fixed_continuous) {

      # get coefficients
      coef_mean <- result$summary.fixed[name, 'mean']
      coef_var <- result$summary.fixed[name, 'sd'] ^ 2

      # get covariate
      cov <- result$model.matrix[, name]

      pred_mean <- pred_mean + cov * coef_mean
      pred_var <- pred_var + (cov ^ 2) * coef_var

    }
  }

  # add discrete fixed effects terms
  if (!is.null(fixed_group)) {
    for (name in fixed_group) {

      # find group members
      idx <- grep(sprintf('^%s*', name), rownames(result$summary.fixed))

      # get clean names
      names <- rownames(result$summary.fixed)[idx]
      names_clean <- gsub(name, '', names)

      # get coefficients
      coef_mean <- result$summary.fixed[idx, 'mean']
      coef_var <- result$summary.fixed[idx, 'sd'] ^ 2

      # get covariate
      cov <- result$model.matrix[, names]

      pred_mean <- pred_mean + as.vector(cov %*% coef_mean)
      pred_var <- pred_var + as.vector((cov ^ 2) %*% coef_var)

    }
  }

  # add discrete random effects terms
  if (!is.null(random_group)) {
    for (name in random_group) {

      # get coefficients
      coef_level <- result$summary.random[[name]][, 'ID']
      coef_mean <- result$summary.random[[name]][, 'mean']
      coef_var <- result$summary.random[[name]][, 'sd'] ^ 2

      # get covariate
      cov <- result$.args$data[, name]

      # match up the levels
      length(cov)
      length(coef_level)
      idx <- match(cov, coef_level)

      pred_mean <- pred_mean + coef_mean[idx]
      pred_var <- pred_var + coef_var[idx]

    }
  }

  # return result
  ans <- data.frame(mean = pred_mean,
                    sd = sqrt(pred_var))

  return (ans)

}

getCols <- function (df, period = 1) {
  # subset results matrix
  df[, grep(sprintf('^p%s_', period),
            colnames(df))]
}

rnormMatrix <- function (n, mean, sd) {
  # sample random normals with matrix parameters
  # returna an array as a result

  # coerce to matrix
  mean <- as.matrix(mean)
  sd <- as.matrix(sd)

  # get & check dimensions
  ncol <- ncol(mean)
  nrow <- nrow(mean)
  stopifnot(ncol(sd) == ncol)
  stopifnot(nrow(sd) == nrow)

  # convert to vector
  mean <- as.vector(mean)
  sd <- as.vector(sd)

  # sample
  draws <- rnorm(n * length(mean), mean, sd)

  # reshape
  ans <- array(draws, dim = c(nrow, ncol, n))

  return (ans)
}

accumulate <- function (mean, sd, months = c(1, 11, 24, 24), nsample = 1000) {
  # given matrices of means and standard deviations of
  # logit component death probabilities, each column giving
  # consecutive and adjacent age bins and rows giving the records,
  # calculate the logit of the cumulative mortality probability across
  # the age bins by Monte Carlo sampling

  # generate random draws from these logits
  draws <- rnormMatrix(nsample, mean, sd)

  # convert to draws of survival probabilities
  draws_p <- 1 - plogis(draws)

  # raise to power of number of months per bin
  for (i in 1:dim(draws_p)[2]) {
    draws_p[, i, ] <- draws_p[, i, ] ^ months[i]
  }

  # loop through age bins accumulating them
  for (i in 2:dim(draws_p)[2]) {
    draws_p[, i, ] <- draws_p[, i, ] * draws_p[, i - 1, ]
  }

  # convert back to logit mortality probabilities
  draws_l <- qlogis(1 - draws_p)

  # calculate the means and standard deviations of these logits
  l_mean <- apply(draws_l, c(1, 2), mean)
  l_sd <- apply(draws_l, c(1, 2), sd)

  # return as a list
  return (list(mean = l_mean,
               sd = l_sd))

}


# functions to centre and scale matrices, columnwise

getCentreScale <- function (x, exclude = NULL, na.rm = TRUE) {
  # get dataframe of centreing and scaling values to convert x
  # to the standard normal. exclude is an optional character vector
  # giving column names to exclude from scaling

  # get means and SDs for all columns
  df <- data.frame(name = colnames(x),
                   mean = colMeans(x, na.rm = na.rm),
                   sd = apply(x, 2, sd, na.rm = na.rm))
  rownames(df) <- NULL

  # replace any zero standard deviations with 1
  # to avoid divide-by-zero errors
  df$sd[df$sd == 0] <- 1

  # if any named covariates are to be excluded, set mean to 0 and sd to 1
  if (!is.null(exclude)) {
    idx <- match(exclude, df$name)
    df$mean[idx] <- 0
    df$sd[idx] <- 1
  }

  return (df)
}

centreScale <- function (x, df, inverse = FALSE) {
  # apply pre-calculated centreing/scaling to matrix x,
  # with fixed dataframe of means/sds df
  # or uncentre/unscale if inverse = TRUE

  # get the centreing/scaling dataframe if not available
  if (is.null(df))
    df <- getCentreScale(x)

  # get index to match up values with column names
  names <- colnames(x)
  idx <- match(names, df$name)

  if (any(is.na(idx))) {
    stop ('could not match up column names with the values in df')
  }

  df <- df[idx, ]

  # apply transformations
  if (!inverse) {
    # move to standard normal

    # centre
    x <- sweep(x, MARGIN = 2, STATS = df$mean, FUN = '-')
    # scale
    x <- sweep(x, MARGIN = 2, STATS = df$sd, FUN = '/')

  } else {
    # inverse case, move from standard normal to original

    # unscale
    x <- sweep(x, MARGIN = 2, STATS = df$sd, FUN = '*')
    # uncentre
    x <- sweep(x, MARGIN = 2, STATS = df$mean, FUN = '+')

  }

  return (x)

}

startLog <- function(file = 'full_run.log') {
  # create a logfile and start logging
  con <- file(file)
  sink(con,
       split = TRUE)
  sink(con,
       type = 'message')

  # report session info
  message(sprintf('# starting log at %s\n', Sys.time()))
  message('# session info:\n')
  print(sessionInfo())
  message('\n# run log:\n')
}

stopLog <- function() {
  # report session info
  message('\n# session info:\n')
  print(sessionInfo())
  message(sprintf('\n# stopping log at %s', Sys.time()))
  # stop logging to the logfile
  sink()
  sink(type = 'message')
}

unpack <- function (tmp) {
  # unpack a nested list where the outer list is unnamed, but each inner
  # list contains a single named element.
  # Assign this element to the parent environment and delete the list
  # it is called on fromt he parent ennvironment

  # get name of list in calling environment
  tmp_name <- deparse(substitute(tmp))

  # get calling environment
  pf <- parent.frame()

  # unpack into a single list
  tmp <- unlist(tmp, recursive = FALSE)

  # loop through assigning to calling environment
  for (i in 1:length(tmp)) {
    assign(names(tmp)[i], tmp[[i]], envir = pf)
  }

  # remove object from calling environment
  rm(list = tmp_name, envir = pf)

  # return nothing
  return (invisible(NULL))
}

prepLines <- function (rate,
                       year,
                       ylab = '',
                       title = '',
                       xlim = c(1999, 2016),
                       ylim = c(0, max(rate)),
                       line_years = c(2000, 2005, 2010, 2015)) {
  # set up the base plot for the custom line chart

  plot(rate ~ year,
       type = 'n',
       ylab = '',
       xlab = '',
       axes = FALSE,
       xlim = xlim,
       ylim = ylim)

  for (ly in line_years) {
    lines(x = rep(ly, 2),
          y = ylim,
          lwd = 3,
          lty = 3,
          col = grey(0.8))
  }

  axis(1,
       lty = 0,
       col.axis= grey(0.4),
       line = -1)

  axis(2,
       las = 2,
       cex.axis = 0.8,
       col = grey(0.4),
       col.axis= grey(0.4))

  title(main = title,
        col.main = grey(0.35),
        cex.main = 1.2,
        line = 0.5)

  title(ylab = ylab,
        col.lab = grey(0.4),
        cex.lab = 1.2)
}

addLines <- function (rate,
                      year,
                      country,
                      countries = NULL,
                      col = grey(0.7),
                      size = 1,
                      border = grey(0.4)) {

  # given vectors: rate (y axis), year (x axis)
  # and country (grouping factor), make a nice line
  # plot with lollipop-like likes with border colour
  # 'border' and fill colour 'col' (repeated if length one)
  # If 'countries' is NULL, all countries are plotted,
  # otherwise only those named in this character vector

  # check inputs
  stopifnot(all.equal(length(rate),
                      length(year),
                      length(country)))


  # sort countries
  all_countries <- sort(unique(country))
  if (is.null(countries)) {
    countries <- all_countries
  } else {
    stopifnot(all(countries %in% all_countries))
  }
  n_ctry <- length(countries)

  # expand col and bg if needed
  if (length(col) == 1) {
    col <- rep(col, n_ctry)
  } else {
    stopifnot(length(col) == n_ctry)
  }

  # loop through countries
  for (i in 1:n_ctry) {

    ctry <- countries[i]

    idx_ctry <- which(country == ctry)

    # dark grey outline
    lines(rate[idx_ctry] ~ year[idx_ctry],
          col = border,
          lwd = 7.5 * size)
    points(rate[idx_ctry] ~ year[idx_ctry],
           col = border,
           cex = 1 * size,
           pch = 16)

    # coloured foreground
    lines(rate[idx_ctry] ~ year[idx_ctry],
          col = col[i],
          lwd = 6 * size)

    points(rate[idx_ctry] ~ year[idx_ctry],
           col = col[i],
           pch = 16,
           cex = 0.85 * size)

  }
}

addLabels <- function (rate,
                       year,
                       country,
                       countries = NULL,
                       col = grey(0.7),
                       gap = diff(range(rate)) / 60,
                       cex = 0.7,
                       adj = 0,
                       xpd = NA,
                       ...) {

  # add country names on RHS
  # arguments as before, with gap to define spacing between labels.
  # dots are passed to text
  require (plotrix)

  # check inputs
  stopifnot(all.equal(length(rate),
                      length(year),
                      length(country)))

  # sort countries
  all_countries <- sort(unique(country))
  if (is.null(countries)) {
    countries <- all_countries
  } else {
    stopifnot(all(countries %in% all_countries))
  }
  n_ctry <- length(countries)

  # expand col and bg if needed
  if (length(col) == 1) {
    col <- rep(col, n_ctry)
  } else {
    stopifnot(length(col) == n_ctry)
  }

  # keep only those for latest year
  max_year <- max(as.numeric(year))
  year_idx <- which(year == max_year)
  rate <- rate[year_idx]
  year <- year[year_idx]
  country <- country[year_idx]

  # keep only those for countries requires
  ctry_idx <- which(country %in% countries)
  rate <- rate[ctry_idx]
  year <- year[ctry_idx]
  country <- country[ctry_idx]

  # order by countries
  ctry_o <- match(countries, country)
  rate <- rate[ctry_o]
  year <- year[ctry_o]
  country <- country[ctry_o]

  # get a good gap
  y_pos <- plotrix::spreadout(rate, gap)
  text(x = max_year + 1,
       y = y_pos,
       labels = country,
       col = col,
       cex = cex,
       adj = adj,
       xpd = xpd,
       ...)

}

splitGeoNames <- function (geo) {
  # split country/years (in format 'country_year') out of rownames
  # of a geostatistical conditional simulation object and add as columns
  splits <- strsplit(rownames(geo), split = '_')
  ctry <- sapply(splits, '[', 1)
  year <- sapply(splits, '[', 2)
  geo <- data.frame(iso3 = ctry,
                    year = year,
                    geo,
                    stringsAsFactors = FALSE)
  rownames(geo) <- NULL
  return (geo)
}

getEst <- function (df, iso3, year) {
  # for each dataframe of estimates, line up the country
  # and year, add object name as prefix to the column nmes and
  # return as a dataframe

  # get object name
  prefix <- deparse(substitute(df))

  # get index
  iso3_year_target <- paste(iso3, year, sep = '_')
  iso3_year_df <- paste(df$iso3, df$year, sep = '_')
  idx <- match (iso3_year_target, iso3_year_df)

  # subset, add prefix to column names and return
  df <- df[idx, -(1:2)]
  colnames(df) <- paste(prefix, colnames(df), sep = '_')
  return (df)
}

elogit <- function (y, n) log ( (y + 0.5) / (n - y + 0.5) )

cirange = function(x){
          z=quantile(x,probs=c(.025,.975),na.rm=T)
          return(z[2]-z[1])
      }

lower = function(x) quantile(x,probs=.025,na.rm=T)

upper = function(x) quantile(x,probs=.975,na.rm=T)

cleanup_inla_scratch <- function(run_date) {

  if(keep_inla_files==FALSE) {

    # Clean up INLA intermediate directories unless user has specified to keep them.
    inla_working_dir <- '<<<< FILEPATH REDACTED >>>>'
    inla_working_dirs <- list.dirs(inla_working_dir, recursive = FALSE)
    inla_working_dirs <- inla_working_dirs[grepl(run_date, inla_working_dirs)]
    for(inla_dir in inla_working_dirs) {
      unlink(inla_dir, recursive=TRUE)
    }

  }

  if(keep_inla_files==TRUE) {

    message('Keeping INLA intermediate files because keep_inla_files==TRUE in config.')
    message(paste0('Files stored here: <<<< FILEPATH REDACTED >>>>'))

  }

}


# Timer Functions ---------------------------------------------------------

## Functions to time within the child stacking regions
##
## General usage:
##    require(tictoc)
##
##    tic("Step 1")
##    **your code here**
##    toc(log = T)
##
##    ticlog <- tic.log(format = F)
##    generate_time_log(ticlog)
##
##  Returns: data table with two columns
##     "step": names of events (e.g. "Step 1")
##     "time": time elapsed (as text: Xh Xm Xs)
##
##  Note: can nest tic/toc pairs

generate_time_log <- function(ticlog) {

  # Set up packages
  require(magrittr)
  require(data.table)

  # Functions in functions
  strip_time <- function(x) {
    sec <- as.numeric(x$toc - x$tic)
    time <- format_time(sec)
    name <- x$msg

    df <- c(name, time) %>%
            t %>%
            as.data.table

    names(df) <- c("step", "time")

    return(df)
  }

  format_time <- function(run_time) {
    run_time <- round(as.numeric(run_time),0)

    hours <- run_time %/% 3600
    remainder <- run_time %% 3600
    minutes <- remainder %/% 60
    seconds <- remainder %% 60

    run_time <- paste0(hours, "h ", minutes, "m ", seconds, "s")
    return(run_time)
  }

  df_out <- lapply(ticlog, strip_time) %>% rbindlist

  return(df_out)

}


clean_model_results_table <- function(rd   = run_date,
                                      regs = Regions,
                                      ages = 0,
                                      nm   = '',
                                      indic = indicator,
                                      ig = indicator_group,
                                      stackers = stacked_fixed_effects,
                                      coefs.sum1 = as.logical(coefs_sum1)){

  str_match <- stringr::str_match

  require(magrittr)

  sharedir <- paste0('<<<< FILEPATH REDACTED >>>>')

  # make loopvars
  lv <- expand.grid(regs,ages)

  # grab formatted model fit objects
  stacker_names <- strsplit(stackers, " + ", fixed=T)[[1]]
  mods <- model_fit_table(lv=lv,rd=rd,nullmodel=nm, indicator = indic, indicator_group = ig,
                          coefs.sum1 = coefs.sum1, stacker_name_vec = stacker_names,
                          use_stacking_covs = use_stacking_covs, use_gp = use_gp,
                          spde_prior = spde_prior)

  # add region column
  all_mods <- lapply(names(mods), function(n) {

                mod <- mods[[n]] %>% as.data.table

                r <- str_match(n,"(.*)_")[1,2]
                a <- str_match(n, "_(.*)")[1,2]
                mod[, region := r]
                mod[, age := a]
                return(mod)
               })

  all_mods <- rbindlist(all_mods)
  colorder <- c("region", "age",
                names(all_mods)[!(names(all_mods) %in% c("region", "age"))])
  setcolorder(all_mods, colorder)

  write.csv(all_mods, paste0(sharedir, indic, "_model_results_table.csv"),
            row.names = F)

  return(mods)

}

## ~~~~~~~~~~~ make table of model results ~~~~~~~~~~~~~~~~~

model_fit_table <- function(lv=loopvars[loopvars[,3]==0,],
                            rd=run_date,
                            nullmodel='',
                            indicator = indicator,
                            indicator_group = indicator_group,
                            holdout = 0,
                            coefs.sum1,
                            use_gp = use_gp,
                            spde_prior = spde_prior,
                            use_stacking_covs = use_stacking_covs,
                            stacker_name_vec){
  # load models
  require(INLA)
  message(sprintf('Pulling together table for %s models',rd))
  tlist=list()
  sharedir <- paste0('<<<< FILEPATH REDACTED >>>>')
  for(i in 1:nrow(lv)){
    reg <- lv[i,1]
    age <- lv[i,2]
    message(sprintf('%i of %i loopvars. %s %i',i,nrow(lv),lv[i,1],lv[i,2]))

    # Recreate the INLA data stack
    pathaddin <- paste0('_bin', age, '_', reg, '_', holdout)
    f=paste0('<<<< FILEPATH REDACTED >>>>/', rd, pathaddin, ".RData")

    if(!file.exists(f)){
      message('FAILED TO OPEN')
    } else {
      load(f)
    }

    stacker_name_vec <- intersect(stacker_name_vec, colnames(df))
    if(length(stacker_name_vec) > 0){
      df = df[,paste0(stacker_name_vec) := lapply(stacker_name_vec,
                                                   function(x) get(paste0(x,'_cv_pred')))]
    }

    # Get the spde
    if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
    input_data <- build_mbg_data_stack(df = df,
                                       fixed_effects = all_fixed_effects,
                                       mesh_s = mesh_s,
                                       mesh_t = mesh_t,
                                       use_ctry_res = FALSE,
                                       use_nugget = FALSE,
                                       stacker_names = stacker_name_vec,
                                       exclude_cs    = stacker_name_vec,
                                       spde_prior = spde_prior)

     spde <- input_data[[2]]

     # Now get & transform model fit
     message('::::loading in INLA fit\n')
     f=paste0(sharedir, indicator, '_model_eb_bin', age, "_", reg, "_0.RData")

    if(!file.exists(f)){
      message('FAILED TO OPEN')
    } else {
      load(f)
    }

    ## columns we'll show return
    keep.cols <- c('0.025quant', '0.5quant', '0.975quant')

    ## other hyperparmas
    hyps <- summary(res_fit)$hyperpar[-(1:2), keep.cols] ## first two rows are
                                                         ## theta1/range, theta2/sd

    if(as.logical(use_gp)){
      if(eval(parse(text=spde_prior))$type=="pc") {
        ## extract values from the fit directly
        range <- res_fit$summary.hyperpar[1,keep.cols]
        nom.var <- res_fit$summary.hyperpar[2,keep.cols]^2
      } else {
        ## now we extract what we need from the fit to get transformed spatial params
        res.field <- inla.spde2.result(res_fit, 'space', spde, do.transf=TRUE)
        
        ## nominal range at 0.025, 0.5, 0.975 quantiles
        range   <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.range.nominal[[1]])
        nom.var <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.variance.nominal[[1]])
      }
      spat.hyps <- rbind(range, nom.var)
      rownames(spat.hyps) <- c('Nominal Range', 'Nominal Variance')
      colnames(spat.hyps) <- keep.cols
    }

    ## fixed effects from coefs.sum1
    if(as.logical(coefs.sum1) & as.logical(use_stacking_covs)){
      fixed.sum1 <- res_fit$summary.random$covar
      fixed.sum1$ID <- NULL
      rownames(fixed.sum1) <- stacker_name_vec
      fixed.sum1 <- fixed.sum1[, keep.cols]
    }else{
      fixed.sum1 <- NULL
    }

    ## all other coefs (e.g. intercept and raw covs)
    fixed <- summary(res_fit)$fixed
    if(is.null(nrow(fixed))){
      fixed <- matrix(fixed, ncol = length(fixed)) ## in the event of one row, convert numeric back to data.frame
      rownames(fixed) <- rownames( summary(res_fit)$fixed )
      colnames(fixed) <- colnames( summary(res_fit)$fixed )
    }
    fixed <- fixed[, keep.cols]

    ## combine the two types of 'fixed' results
    fixed <- rbind(fixed, fixed.sum1)

    ## combine them all and just keep three quantiles
    all.res <- rbind(fixed,
                     if(use_gp){spat.hyps}else{NULL},
                     hyps)

    ## rename GPRandom rho for time
    all.res <- as.data.table(all.res, keep.rownames = T)
    setnames(all.res, "rn", "parameter")
    if(use_gp) all.res[parameter == "GroupRho for space", parameter := "GPRandom rho for time"]
    all.res

    tlist[[paste0(reg, "_", age)]] <- all.res

  }

  return(tlist)
}

check_config <- function(cr = core_repo) {
  must_haves <- read.csv(paste0(cr, '<<<< FILEPATH REDACTED >>>>/config_must_haves.csv'), header = F, stringsAsFactors = F)$V1
  
  message("\nRequired covariates: ")
  for(confs in must_haves){
    if (exists(confs)) {
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'use_gp') {
      message("You are missing a 'use_gp' argument in your config. Defaulting it to TRUE")
      use_gp <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'use_stacking_covs') {
      message("You are missing a 'use_stacking_covs' argument in your config. Defaulting it to TRUE")
      use_stacking_covs <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'use_raw_covs') {
      message("You are missing a 'use_raw_covs' argument in your config. Defaulting it to FALSE")
      use_raw_covs <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'fit_with_tmb') {
      message("You are missing a 'fit_with_tmb' argument in your config. Defaulting it to FALSE")
      fit_with_tmb <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'gbd_fixed_effects_measures') {
      message("You are missing a 'gbd_fixed_effects_measures' argument in your config. Defaulting it to 'covariate' for all elements of gbd_fixed_effects")
      gbd_fixed_effects_measures <<- paste(rep("covariate", length(strsplit(gbd_fixed_effects, split = " \\+ ")[[1]])), collapse = " + ")
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'gbd_fixed_effects_age') {
      message("You are missing a 'gbd_fixed_effects_age' argument in your config. Defaulting to '2 3 4 5'")
      gbd_fixed_effects_age <<- '2 3 4 5'
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'z_list') {
      message("You are missing a 'z_list' argument in your config. Defaulting it to 0")
      z_list <<- 0
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'zcol') {
      message("You are missing a 'zcol' argument in your config. Defaulting it to z_column_default_blank")
      zcol <<- 'z_column_default_blank'
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'summstats') {
      message("You are missing a 'summstats' argument in your config. Defaulting to c('mean','lower','upper','cirange')")
      summstats <<- c('mean','lower','upper','cirange')
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'scale_gaussian_variance_N') {
      message("You are missing a 'scale_gaussian_variance_N' argument in your config. Defaulting to TRUE")
      scale_gaussian_variance_N <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "nugget_prior") {
      message("You are missing a 'nugget_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      nugget_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "ctry_re_prior") {
      message("You are missing a 'ctry_re_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      ctry_re_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "nid_re_prior") {
      message("You are missing a 'nid_re_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      nid_re_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "use_nid_res") {
      message("You are missing a 'use_nid_res' argument in your config. Defaulting to FALSE")
      use_nid_res <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "rho_prior") {
      message("You are missing a 'rho_prior' argument in your config. Defaulting to 'list(prior = 'normal', param = c(0, 0.1502314))'")
      rho_prior <<- "list(prior = 'normal', param = c(0, 1/(2.58^2)))"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "gp_sum_to_zero") {
      message("You are missing a 'gp_sum_to_zero' argument in your config. Defaulting to FALSE")
      gp_sum_to_zero <<- FLASE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "use_s2_mesh") {
      message("You are missing a 'use_s2_mesh' argument in your config. Defaulting to FALSE")
      use_s2_mesh <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "s2_mesh_params") {
      message("You are missing a 's2_mesh_params' argument in your config. Defaulting to c(50, 500, 1000)")
      s2_mesh_params <<- "c(25, 500, 1000)"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "sparse_ordering") {
      message("You are missing a 'sparse_ordering' argument in your config. Defaulting to TRUE")
      sparse_ordering <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "modeling_shapefile_version") {
      message("You are missing a 'modeling_shapefile_version' argument in your config. Defaulting to 'current'")
      modeling_shapefile_version <<- "current"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "raking_shapefile_version") {
      message("You are missing a 'raking_shapefile_version' argument in your config. Defaulting to 'current'")
      raking_shapefile_version <<- "current"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "subnational_raking") {
      message("You are missing a 'subnational_raking' argument in your config. Defaulting to TRUE")
      subnational_raking <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "check_cov_pixelcount") {
      message("You are missing a 'check_cov_pixelcount' argument in your config. Defaulting to FALSE")
      check_cov_pixelcount <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "gbd_fixed_effects_constraints") {
      message("You are missing a 'gbd_fixed_effects_constraints' argument in your config. Defaulting to FALSE")
      gbd_fixed_effects_constraints <<- "c(0)"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "fixed_effects_constraints") {
      message("You are missing a 'fixed_effects_constraints' argument in your config. Defaulting to FALSE")
      fixed_effects_constraints <<- "c(0)"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "memory") {
      message("You are missing a 'memory' argument in your config. Defaulting to 10G")
      memory <<- 10
      
    } else if (confs == "singularity_version") {
      message("You are missing a 'singularity_version' argument in your config. Defaulting to 'default'")
      singularity_version <<- "default"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "queue") {
      message("You are missing a 'queue' argument in your config. Defaulting to 'long.q', unless you have use_geos_nodes to TRUE, which will override this to geospatial.q")
      queue <<- "long.q"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "run_time") {
      message("You are missing a 'run_time' argument in your config. Defaulting to 16 days ('16:00:00:00')")
      run_time <<- "16:00:00:00"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "countries_not_to_rake") {
      message("You are missing a 'countries_not_to_rake' argument in your config. Defaulting to ESH+GUF")
      countries_not_to_rake <<- "ESH+GUF"
      message(paste0("  ", confs, ": ", get(confs)))
      
    } else if (confs == "countries_not_to_subnat_rake") {
      message("You are missing a 'countries_not_to_subnat_rake' argument in your config. Defaulting to PHL+NGA+PAK+ETH+KEN")
      countries_not_to_subnat_rake <<- "PHL+NGA+PAK+ETH+KEN"
      message(paste0("  ", confs, ": ", get(confs)))
      
    } else if (confs == "rake_countries") {
      message("You are missing a 'rake_countries' argument in your config. Defaulting to TRUE")
      rake_countries <<- TRUE
      message(paste0("  ", confs, ": ", get(confs)))
      
    } else if (confs == "use_space_only_gp") {
      message("You are missing a 'use_space_only_gp' argument in your config. Defaulting to FALSE")
      use_space_only_gp <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "st_gp_int_zero") {
      message("You are missing a 'st_gp_int_zero' argument in your config. Defaulting to FALSE")
      st_gp_int_zero <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "s_gp_int_zero") {
      message("You are missing a 's_gp_int_zero' argument in your config. Defaulting to FALSE")
      s_gp_int_zero <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "use_time_only_gmrf") {
      message("You are missing a 'use_time_only_gmrf' argument in your config. Defaulting to FALSE")
      use_time_only_gmrf <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "time_only_gmrf_type") {
      message("You are missing a 'time_only_gmrf_type' argument in your config. Defaulting to 'rw2'")
      time_only_gmrf_type <<- "rw2"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "spde_prior") {
      message("You are missing a 'spde_prior' argument in your config. Defaulting to 'list(type='pc')'")
      spde_prior <<- "list(type='pc')"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else {
      stop(paste0(confs, " is missing, add it to your config"))
    }
  }


  ## Test for subnational random effect
  if(exists("use_subnat_res", envir = .GlobalEnv)) {
    stopifnot(exists("subnat_country_to_get", envir = .GlobalEnv))
    # stopifnot(length(eval(parse(text = subnat_country_to_get))) == 1)
  } else {
    use_subnat_res <<- FALSE
    subnat_country_to_get <<- FALSE
  }


  message("\nAdditional covariates: ")
  extras <- config$V1[!(config$V1 %in% must_haves)]
  for (extra in extras) message(paste0('  ', extra, ': ', get(extra)))

  ## print out shapeilfe info
  m.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = modeling_shapefile_version))
  r.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = raking_shapefile_version))
  message("\n\n\nSHAPEFILE VERSION INFORMATION: ")
  message(sprintf("\n--MODELING SHAPEFILE VERSION: %s -- which contains %s codes", m.sf.info$shpfile_date, toupper(m.sf.info$shpfile_type)))
  message(sprintf("\n--RAKING SHAPEFILE VERSION:   %s -- which contains %s codes\n", r.sf.info$shpfile_date, toupper(r.sf.info$shpfile_type)))
}


#' @title Easy eval-parse
#' @description Allows for easily eval-parsing through a config dataset
#' @param data The data.table 
#' @param column The column with the string call
#' @return Evaluated call
#' @export
#' @rdname ez_evparse
ez_evparse <- function(data, column) {
  return(eval(parse(text = data[, column, with = FALSE])))
}






#' @title Set up config
#' @description Setting up configuration variables for an MBG run
#' @param repo Location where you've cloned the MBG repository for your indicator.
#' @param core_repo Location where you've cloned the lbd_core repository. Not necessary in the package version.
#' @param indicator_group Category of indicator, i.e. "education"
#' @param indicator Specific outcome to be modeled within indicator category, i.e. "edu_0"
#' @param config_name Name of configuration file in the indicator folder, Default: NULL
#' @param config_file Full path to configuration file that overrides \code{config_name}, Default: NULL
#' @param covs_name Name of covariates configuration file, Default: NULL
#' @param covs_file Full path to covariates configuration file that overrides \code{covs_name}, Default: NULL
#' @param post_est_only Set up only for post estimation? Default: FALSE
#' @param run_date Run date, Default: ''
#' @param push_to_global_env Should the config parameters be pushed to the global environment? Default: TRUE
#' @param run_tests Run the assertion tests? This will run the tests and error out if there's an 
#' inconsistent config parameter. Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   config <- load_config(repo = core_repo,
#'     indicator_group = indicator_group,
#'     indicator = indicator,
#'     config_name = 'config_training',
#'     covs_name = 'covs_training')
#' }
#' }
#' @rdname set_up_config
#' @importFrom assertthat is.flag is.string is.number
#' @export
set_up_config <- function(repo, 
                          core_repo = repo,
                          indicator_group, 
                          indicator, 
                          config_name = NULL, 
                          config_file = NULL, 
                          covs_name = NULL, 
                          covs_file = NULL,
                          post_est_only = FALSE, 
                          run_date = "", 
                          push_to_global_env = TRUE,
                          run_tests = TRUE) {
  
  ###### Block 1: Equivalent to load_config ######
  
  print("[1/6] Load the configs")
  
  ####### Logic checking for model config ####### 
  ## Make sure only one of config_name or config_file are not null
  if (!is.null(config_name) & !is.null(config_file)) {
    stop("You must specify just one of config_name or config_file, not both", call. = FALSE)
  }
  
  ## Pull config from indicator repo
  if (is.null(config_name) & is.null(config_file)) {
    message("Pulling config from default folder, since config_name and config_file are NULL")
    ## If new model run, pull config from /share repo
    if (post_est_only == FALSE) 
      config <- fread(paste0(repo, "/<<<< FILEPATH REDACTED >>>>/config_", indicator, ".csv"), header = FALSE)
    ## If running analysis on existing model, use config from that model's outputs folder
    if (post_est_only == TRUE) 
      config <- fread(paste0("<<<< FILEPATH REDACTED >>>>/config.csv"))
  }
  
  ## Pull by specific config name
  if (!is.null(config_name) & is.null(config_file)) {
    message("Pulling config from specified name")
    config <- fread(paste0(repo, "/<<<< FILEPATH REDACTED >>>>/", config_name, ".csv"), header = FALSE)
  }
  ## Pull specified config file
  if (is.null(config_name) & !is.null(config_file)) {
    message("Pulling config from specified filepath")
    config <- fread(config_file, header = FALSE)
  }
  
  ####### Logic checking for covariates config ####### 
  ## Make sure only one of covs_name or covs_file are not null
  if (!is.null(covs_name) & !is.null(covs_file)) {
    stop("You must specify just one of covs_name or covs_file, not both", call. = FALSE)
  }
  
  ## Covs not pulled
  if (is.null(covs_name) & is.null(covs_file)) {
    message("Not pulling covs since covs_name and covs_file are NULL")
    covs <- NULL
  }
  
  ## Pull by specific covs name
  if (!is.null(covs_name) & is.null(covs_file)) {
    message("Pulling covs from specified name")
    covs <- fread(paste0(repo, "/<<<< FILEPATH REDACTED >>>>/", covs_name, ".csv"), header = TRUE)
  }
  
  ## Pull specified covs file
  if (is.null(covs_name) & !is.null(covs_file)) {
    message("Pulling covs from specified filepath")
    covs <- fread(covs_file, header = TRUE)
  }
  
  ## For parsimony, let's make sure that the config column names are V1 and V2
  config <- data.table(config)
  if(colnames(config)[1] != "V1" & colnames(config)[2] != "V2") {
    warning("Renaming config column names to V1 and V2. Please verify that 'config' is properly built")
    colnames(config) <- c("V1", "V2")
  }
  
  
  # If a covariate .csv file exists, use that instead
  if (!is.null(covs)) {
    
    # Grab fixed effects & measures (and gbd fixed effects & measures) from CSV if present
    
    # After update to data.table 1.11.4, 'T' and 'F' are not read in as logical, 
    ## but as characters, which we need to remedy here. 
    ## We are assuming that the 'covs.csv' has an 'include' and 'gbd' column here
    covs[, `:=`(gbd, as.logical(gbd))]
    covs[, `:=`(include, as.logical(include))]
    covs <- subset(covs, include == T)  # Use only those where include flag set to true
    fe <- subset(covs, gbd == F)
    gbd <- subset(covs, gbd == T)
    gbd[measure != "output", `:=`(measure, "covariate")]  # FIXME: This is a hack for backwards compatability -- basically it assumes you meant 'covariate' if you specified anything other than 'outcome' (eg, mean or NA)
    fixed_effects <- paste0(fe$covariate, collapse = " + ")
    fixed_effects_measures <- paste0(fe$measure, collapse = " + ")
    gbd_fixed_effects <- paste0(gbd$covariate, collapse = " + ")
    gbd_fixed_effects_measures <- paste0(gbd$measure, collapse = " + ")
    
    # Remove any other versions from original config and 
    # override with covariates config outputs
    for(varz in c("fixed_effects", "fixed_effects_measures", "gbd_fixed_effects", "gbd_fixed_effects_measures")) {
      if(!(varz %in% colnames(config))) {
        config <- rbindlist(list(config, data.table(V1 = varz, V2 = get(varz))))
      } else {
        config[V1 == varz, V2:= get(varz)]
      }
    }
    
  }
  
  
  ###### Block 2: Add fields in config that are not in the default set ######
  
  print("[2/6] Add fields in config that are not in the default set")
  
  ## Load in the default config dataset
  config_values <- fread(paste0(core_repo, '/<<<< FILEPATH REDACTED >>>>/config_values.csv'), header = TRUE, stringsAsFactors = FALSE)
  
  ## Now, go through each of the values in `config_values` and 
  ## add on all the fields that are not in the user-specified config
  for(conf_vars in names(config_values)) {
    if(!(conf_vars %in% config$V1)) {
      message(paste0(conf_vars, ' not found in config. Adding in default value of: ', config_values[1, conf_vars, with=F]))
      config <- rbindlist(list(config, data.table(V1 = conf_vars, V2 = config_values[1, conf_vars, with=F] )))
    }
  }
  
  
  ###### Block 3: Extra parameters in config ######
  
  print("[3/6] Add fields that are in in config but not in the default set")
  message("\nAdditional covariates: ")
  extras <- config$V1[!(config$V1 %in% names(config_values))]
  for (extra in extras) {
    message(paste0("  ", extra, ": ", config[V1 == extra, V2] ))
  }
  
  ###### Block 4: Print out shapefile info from config ######
  
  print("[4/6] Print out shapefile info from config")
  m.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = config[V1 == 'modeling_shapefile_version', V2]))
  r.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = config[V1 == 'raking_shapefile_version', V2]))
  message("\n\n\nSHAPEFILE VERSION INFORMATION: ")
  message(sprintf("\n--MODELING SHAPEFILE VERSION: %s -- which contains %s codes", m.sf.info$shpfile_date, toupper(m.sf.info$shpfile_type)))
  message(sprintf("\n--RAKING SHAPEFILE VERSION:   %s -- which contains %s codes\n", r.sf.info$shpfile_date, toupper(r.sf.info$shpfile_type)))
  
  
  ###### Block 5: Run tests on all the configuration variables loaded ######
  if(run_tests) {
    print("[5/6] Running simple type-assertion tests on config parameters")
    config_tests <- fread(paste0(core_repo, '/<<<< FILEPATH REDACTED >>>>/config_tests.csv'), header = TRUE, stringsAsFactors = FALSE)
    
    ## Test for params only in the config_tests list of params
    for (param in sort(config[, V1])) {
      cat(paste0("Testing config parameter: ", param))
      if(param %in% config_tests$variable) {
        test_call_1 <- config_tests[variable == param, test_call]
        test_call_2 <- config_tests[variable == param, extra_test1]
        test_call_3 <- config_tests[variable == param, extra_test2]
        
        if(test_call_1 != "") {
          ## For a string in the config file, the eval-parse combo will
          ## fail to evaluate it, and so we build in this exception for that
          tryCatch(
            get(test_call_1)(ez_evparse(config[V1 == param, ], "V2")),
            error = function(e) {
              if(attributes(e)$class[[1]] == 'simpleError' & test_call_1 == "is.string") {
                message(paste0("Assertion on ", param, " errored out because it's tested as a string. Please check for the real type manually"))
              }
            }
          )
        }
        if(test_call_2 != ""  ) {
          tryCatch(
            assertthat::assert_that(eval(parse(text = test_call_2))),
            error = function(e) {
              stop(paste0("The following test failed: ", test_call_2) )
            }
          )
        }
        if(test_call_3 != ""  ) {
          tryCatch(
            assertthat::assert_that(eval(parse(text = test_call_3))),
            error = function(e) {
              stop(paste0("The following test failed: ", test_call_3) )
            }
          )
        }
        cat(" OK. \n")
      }
    }
      
  } else {
    warning("[5/6] Skipping over type-assertion")
  }
  
  
  

  ###### Final Block :  Assign all the covariates to the environment if desired ######
  
  if(push_to_global_env) {
    print("[6/6] Pushing config parameters into global environment")
    for (param in config[, V1]) { 
      assign(param, config[V1 == param, V2], envir = globalenv())
    }
  } else {
    print("[6/6] Config parameters not passed into global environment")
  }
  
  ## Return the config data.table
  message("Saving out config...")
  return(config)
  
  
}





## proportional_area_map ################################################

#' Create proportional area maps for count data at various admin levels
#'
#' @param data data frame or data table with at least a value and ADMX_CODE column
#' @param ad_level admin level to map (0,1,2)
#' @param value_col column in `data` holding the value to plot (proportional to size)
#' @param main_title main title for the map (default = NULL)
#' @param legend_title legend title for the map
#' @param fill_color color to use to fill in the proportional area circles
#' @param alpha transparency (alpha) of the proportional area circles
#' @param scale_size_max max size of the proportional area circles
#'                       (passed to `scale_size_area()`)
#' @param scale_size_breaks breaks to specify for size scale legend
#'                          (passed to `scale_size_area()`)
#' @param scale_size_breaks labels to specify for size scale legend
#'                          (passed to `scale_size_area()`)
#' @param lake_river_color controls the color of the lakes/rivers
#' @param out_file output file name (png)
#' @param out_file_height height of output file (inches)
#' @param out_file_width width of ouput file (inches)
#' @param out_file_res resolution of output file (dpi)
#' @param out_file_pointsize pointsize of text in output file (generally would first
#'                           adjust legend text/title size as below with the relevant
#'                           options, and then use this to modify line spacing, etc.)
#' @param legend_text_size adjust the size of the legend text
#' @param legend_title_size adjust the size of the legend title
#'
#' @return ggplot object containing the map object; writes png to `out_file` if not NULL
#' @examples
#'
#' # Not run
#' a1_df <- fread("<<<< FILEPATH REDACTED >>>>/[indic]_admin_1_raked_summary.csv")
#' a1_df <- subset(a1_df, year == 2015)
#'
#' gg_obj <- proportional_area_map(data = a1_df,
#'                                 ad_level = 1,
#'                                 value_col = "mean",
#'                                 main_title = NULL,
#'                                 legend_title = "DPT3: 2015",
#'                                 fill_col = "red",
#'                                 alpha = 0.25,
#'                                 size_scale_max = 5,
#'                                 out_file = "/my/out/file.png")


proportional_area_map <- function(data, # Needs to have ADM0_CODE, ADM1_CODE, ADM2_CODE
                                  ad_level,
                                  value_col,
                                  main_title = NULL,
                                  legend_title,
                                  fill_col,
                                  alpha = 0.25,
                                  size_scale_max = 6,
                                  size_scale_breaks = NULL,
                                  size_scale_labels = NULL,
                                  lake_river_color = "lightblue",
                                  out_file = NULL,
                                  out_file_height=12,
                                  out_file_width=12,
                                  out_file_res=300,
                                  out_file_pointsize = 16,
                                  legend_text_size = 16,
                                  legend_title_size = 18,
                                  shapefile_version = 'current') {

  # Set up and define objects
  library("ggplot2")
  library("data.table")
  library("magrittr")
  library("raster")
  library("rgeos")

  # Set up admin column code
  adm_code_col <- paste0("ADM", ad_level, "_CODE")

  # Copy input data to avoid environment/reference issues
  df <- copy(as.data.table(data))

  message("Loading background map...")
  message("Note: background map only available for Africa for now")
  background_shp <- readRDS('/<<<< FILEPATH REDACTED >>>>/background_map_africa.rds')
  background_map <- fortify(background_shp)

  message("Loading master admin shape...")
  ad_shape <- rgdal::readOGR(get_admin_shapefile(ad_level, version = shapefile_version))
  if (ad_level == 2) {
    ad_shape <- subset(ad_shape, ADM2_CODE %in% unique(df[, get(adm_code_col)]))
  } else if (ad_level == 1) {
    ad_shape <- subset(ad_shape, ADM1_CODE %in% unique(df[, get(adm_code_col)]))
  } else if (ad_level == 0) {
    ad_shape <- subset(ad_shape, ADM0_CODE %in% unique(df[, get(adm_code_col)]))
  }

  message("Loading masks, lakes & rivers...")

  ### Load & process mask for ggplot2
  mask <- raster('/<<<< FILEPATH REDACTED >>>>/mask_master.tif')
  mask_df <- rasterToPoints(mask)
  mask_df <- data.frame(mask_df)
  colnames(mask_df) <- c("long", 'lat', 'mask')

  ### Load & process lakes & rivers for ggplot2
  lakes <- raster('/<<<< FILEPATH REDACTED >>>>/lakes_all_2.tif')
  lakes_df <- rasterToPoints(lakes)
  lakes_df <- data.frame(lakes_df)
  colnames(lakes_df) <- c("long", 'lat', 'lakes')

  message("Getting centroids...")
  if (ad_level == 2) centroids <- gCentroid(ad_shape, byid = T, id = ad_shape$ADM2_CODE)
  if (ad_level == 1) centroids <- gCentroid(ad_shape, byid = T, id = ad_shape$ADM1_CODE)
  if (ad_level == 0) centroids <- gCentroid(ad_shape, byid = T, id = ad_shape$ADM0_CODE)

  centroids <- as.data.frame(centroids) %>%
    cbind(rownames(.), .) %>%
    as.data.table(.) %>%
    setnames(., names(.), c("adm_code", "x", "y"))
  centroids$adm_code <- as.numeric(levels(centroids$adm_code))[centroids$adm_code]
  setnames(centroids, "adm_code", adm_code_col)

  message("Merging centroids to data")
  df <- merge(df, centroids, by = eval(adm_code_col), all.x = T, all.y = F)
  n_missing <- nrow(df[is.na(x) | is.na(y),])

  if (n_missing > 0) {
    missing_ad_shapes <- unique(df[is.na(x) | is.na(y), get(adm_code_col)])
    message(paste0("The following admin shapes are missing: ", missing_ad_shapes))
  }

  message("Preparing to plot...")
  ad_shape_df <- fortify(ad_shape, region = adm_code_col)

  # A custom, mostly blank theme
  theme_empty <- theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank())

  # Set up an object to hold the scale size options
  if (is.null(size_scale_breaks) & is.null(size_scale_labels)) {

    size_scale_object <- scale_size_area(max_size = size_scale_max)

  } else if (!is.null(size_scale_breaks) & is.null(size_scale_labels)) {

    size_scale_object <- scale_size_area(max_size = size_scale_max,
                                         breaks = size_scale_breaks)

  } else if (is.null(size_scale_breaks) & !is.null(size_scale_labels)) {

    size_scale_object <- scale_size_area(max_size = size_scale_max,
                                         labels = size_scale_labels)

  } else if (!is.null(size_scale_breaks) & !is.null(size_scale_labels)) {

    size_scale_object <- scale_size_area(max_size = size_scale_max,
                                         breaks = size_scale_breaks,
                                         labels = size_scale_labels)

  }

  # For convenience change df name of plotting column
  setnames(df, value_col, "plot_me")

  message("Plotting...")
  g_map <- ggplot() +

    # Plot background
    geom_polygon(data = background_map,
                 aes(x=long,
                     y=lat,
                     group = group),
                 fill = "white")  +

    # Plot national country outlines
    geom_path(data = background_map,
              aes(x = long,
                  y = lat,
                  group = group),
              size = 0.5,
              color = "black")   +

    # Plot admin unit outlines
    geom_path(data = ad_shape_df,
              aes(x = long,
                  y = lat,
                  group = group),
              size = 0.3,
              color = "black") +

    # Plot lakes & rivers
    annotate(geom = 'raster', x = lakes_df$long, y = lakes_df$lat,
             fill = lake_river_color,
             alpha = 0.75) +

    # Plot points
    geom_point(data = df,
               aes(x = x,
                   y = y,
                   size = plot_me),
               color = fill_col,
               alpha = alpha) +

    size_scale_object +

    # Theme/style additions & labels
    theme_empty +
    coord_equal(ratio = 1) +
    labs(size = legend_title,
         main = title) +
    guides(size = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.text=element_text(size=legend_text_size)) +
    theme(legend.title=element_text(size=legend_title_size))

  if (!is.null(out_file)) {
    message(paste0("Saving to ", out_file))
    png(filename = out_file,
        type = "cairo",
        units = "in",
        width = out_file_width,
        height = out_file_height,
        pointsize = out_file_pointsize,
        res = out_file_res)

    print(g_map)
    dev.off()
  }

  return(g_map)
}


## New suite of functions to use for qsubbing ############################

## parallelize ################################################

#' parallelize() is a versatile function to take an R script and run it in
#' parallel on the cluster combination of an arbitrary number of variables.
#'
#' This function is meant to replace the many qsub functions that are
#' floating around and provide a single new function that can be used in
#' almost all circumstances, with several additional features:
#'
#' - By pairing this function with a `load_from_parallelize()` call in the
#'   child script itself, objects are loaded into the child script's
#'   environment without having to ensure that a series of `commandArgs()`
#'   are in the appropriate order
#'
#' - `parallelize()` returns a list of job_ids and loop variables, along
#'   with the original qsub call.  This object - when paired with `monitor_jobs()`,
#'   which is a replacement for `waitformodelstofinish()` - allows closer tracking
#'   of the jobs with respect to their status on the cluster, and eliminates
#'   the need to write clunky empty files like `fin_[whatever]` to mark that
#'   jobs are done.  Finally, the `monitor_jobs()` function can automatically
#'   resubmit jobs and notify the user (in progress; Pushover notifications only
#'   currently supported) when jobs fail.
#'
#' @param user current user [default = user name]
#' @param slots number of slots to request for each job [required]
#' @param memory memory to request (in GB) [required]
#' @param script R script to be run in parallel. Assumes that script name
#'               ends with '.R'. This script should have a
#'               `load_from_parallelize()` call towards the top. [required]
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
#' @param ig indicator group [default = indicator_group]
#' @param indic indicator [default = indicator]
#' @param rd run date [default = run_date]
#' @param expand_vars a named list of objects to `grid.expand()` over.  One job
#'                    will be submitted for each named item in the list, using
#'                    the name of the item as the variable name.  For instance:
#'                    `expand_vars = list(region = c("cssa", "essa", "sssa"))`
#'                    would submit one job for each region with `region = "cssa"
#'                    in the first job, `region = "essa"` in the second job, and
#'                    so forth. If a second item were added to that list, then
#'                    all combinations will be submitted:
#'                    `expand_vars = list(region = c("cssa", "essa", "sssa"),
#'                                        raked = c(T,F))`
#'                    would submit cssa-raked, cssa-unraked, essa-raked,
#'                    essa-unraked, etc...
#'                    Only `expand_vars` or `lv_table` can be given, but not both.
#'                    [default = NULL]
#' @param save_objs character vector of objects that should be available to all
#'                  child scripts. Unlike `expand_vars`, these will be the *same*
#'                  throughout all of the child scripts - they are saved to a
#'                  temporary file and loaded by each child script. For instance,
#'                  `save_objs = c("run_date", "indicator_group")` would load the
#'                  `run_date` and `indicator_group` objects into the environment
#'                  of each child script. [default = NULL]
#' @param lv_table pass in the loop vars table? Pass in this table if you want to
#'                 supply something more fine-tuned than what can by given by
#'                 `grid.expand()` with `expand_vars`. Only `expand_vars` or
#'                 `lv_table` can be given, but not both. [default = NULL]
#' @param script_dir directory to look in for the R script that is to be run in
#'                   parallel. If this is NULL, then script will look in
#'                   'corerepo/mbg_central/share_scripts' (see `corerepo` below)
#'                   for the script given in `script`. [default = NULL]
#' @param prefix prefix to be appended to all jobs.  The jobs will have the naming
#'               convention of prefix_[first_expand_var]_[second_expand_var]_[etc]
#'               [default = 'job']
#' @param log_location where to save the logs? [default = 'sgeoutput']
#' @param corerepo location of the lbd_core repo [default = core_repo]
#' @param geos_nodes run on the geos nodes or not? Defaults to running on prod
#'                   with the 'proj_geospatial' project (see the 'proj' arg).
#'                   [default = FALSE]
#' @param runtime Run-time for usage in the fair cluster (unused with prod)
#' @param priority Job priority that can be deprioritized if needed, and can only be used for values in [-1023,0]. Default = 0. 
#' This value will get bounded to 0 or -1023 if the user supplies a value outside those bounds.
#' @param threads numeric number of threads to request on fair cluster (unused with prod)
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
#'   the default location of <<<< FILEPATH REDACTED >>>>.
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
#'
#' @return a list containing:
#'           -`lv`: data table of loop variables including qsub commands
#'           -`fname`: filename of the temporary file containing `save_objs`
#'
#' @examples
#'
#' # In master script:
#'
#' if (as.logical(makeholdouts) == T) holdout_vector <- 0:as.numeric(n_ho_folds)
#' if (as.logical(makeholdouts) == F) holdout_vector <- 0
#'
#' combine_lv<- list(region = Regions,
#'                   holdout = holdout_vector,
#'                   doses = 3)
#'
#' combination_output <- parallelize(script = "indicator_specific_script",
#'                                   script_dir = "<<<< FILEPATH REDACTED >>>>",
#'                                   expand_vars = combine_lv,
#'                                   save_objs = c("indicator_group", "run_date",
#'                                                 "vaccine"),
#'                                   prefix = "combine",
#'                                   cores = 20,
#'                                   memory = 100,
#'                                   geo_nodes = TRUE)
#'
#' monitor_jobs(combination_output)
#'
#' @seealso \code{\link{load_from_parallelize}}: Child scripts launched
#'   by parallelize should call this function near the top.
#'   \code{\link{get_singularity}}: Determines which Singularity to use.
#'   \code{get_singularity}: Which Singularity image to use
#'   \code{qsub_sing_envs}: Adds environmental variables to a qsub string
#'
parallelize <- function(user = Sys.info()['user'],
                        slots,
                        memory,
                        script,
                        proj             = NULL,
                        ig               = indicator_group,
                        indic            = indicator,
                        rd               = run_date,
                        expand_vars      = NULL,
                        save_objs        = NULL,
                        lv_table         = NULL,
                        script_dir       = NULL,
                        prefix           = 'job',
                        log_location     = 'sgeoutput',
                        corerepo         = core_repo,
                        geo_nodes        = FALSE,
                        use_c2_nodes     = FALSE,
                        queue            = NULL,
                        run_time         = NULL,
                        priority         = 0,
                        threads          = 2,
                        singularity      = singularity_version,
                        singularity_opts = NULL) {

  # Setup ---------------------------------------------------------------
  # allocate cores based on threads argument if on the new cluster. Otherwise respect historical "slot" usage
  cores <- ifelse(is_new_cluster(), threads, slots)
  # Attempt to locate the R script we intend to run and verify that it
  # actually exists
  # We assume that the script ends with '.R'. If the user supplied
  # the file name with a '.R' appended make sure not to append one again
  script <- ifelse(stringr::str_sub(script, -2) == '.R', script, paste0(script, '.R'))
  if(is.null(script_dir)) script_dir <- paste0(corerepo, '<<<< FILEPATH REDACTED >>>>')
  script_file <- paste0(script_dir, "/", script)
  if(!file.exists(script_file)) {
    stop(paste0("Could not locate R script: ", script_file))
  }

  # Define project first (necessary to validate node options)
  proj <- get_project(proj, use_geo_nodes=geo_nodes)

  # Validate arguments
  validate_singularity_options(singularity, singularity_opts)
  validate_node_option(geo_nodes, use_c2_nodes, proj)

  # Determine where stdout and stderr files will go
  output_err = setup_log_location(log_location, user, indic, ig, rd)
  output_log_dir = output_err[[1]]
  error_log_dir = output_err[[2]]

  # Define remaining attributes
  run_time <- get_run_time(use_geo_nodes=geo_nodes, use_c2_nodes=use_c2_nodes, queue=queue, run_time=run_time)
  queue <- get_queue(use_geo_nodes=geo_nodes, use_c2_nodes=use_c2_nodes, queue=queue, run_time=run_time)
  shell <- paste0(corerepo, '/<<<< FILEPATH REDACTED >>>>/shell_sing.sh')
  sing_image <- get_singularity(image = singularity)
  
  # resources are all the -l qsub arguments
  resources <- get_resources(use_geo_nodes=geo_nodes, cores=cores, ram_gb=memory, runtime=run_time)

  # Expand loop variables from expand_vars or lv_table (depending on what was provided)
  if (is.null(lv_table) & is.null(expand_vars)) {
    stop("Need to have one of either lv_table or expand_vars")
  }
  if (!is.null(lv_table) & !is.null(expand_vars)) {
    stop("Can only have one of either lv_table or expand_vars")
  }

  if (is.null(lv_table) & !is.null(expand_vars)) {
    lv <- data.table(expand.grid(expand_vars, stringsAsFactors = F))
  } else if (!is.null(lv_table) & is.null(expand_vars)) {
    lv <- copy(as.data.table(lv_table))
  }

  lv$jobname <- paste0(prefix, "_", apply(lv, 1, paste, collapse = "_"))

  # Ensure all whitespace collapsed
  lv[, jobname := gsub(" ", "", jobname, fixed = TRUE)]

  # Save objects ---------------------------------------------------------

  # Create filename using time stamp
  tmp_dir <- "<<<< FILEPATH REDACTED >>>>"
  # Append random string on the end to avoid overlapping filenames for
  # rapidly-submitted jobs
  fname <- paste0(user, "_", gsub("-|:| ","_",Sys.time()),  sample(1:100000, 1))

  values_to_save <- c("lv")
  if (!is.null(save_objs)) values_to_save <- c(values_to_save, save_objs)
  save(file = paste0(tmp_dir, fname, ".RData"),
       list = values_to_save)

  # Qsub over lv rows
  for (i in 1:nrow(lv)) {
    job_name <- lv[i, jobname]

    qsub <- generate_qsub_command(
      # qsub-specific arguments
      stderr_log=error_log_dir,
      stdout_log=output_log_dir,
      project=proj,
      resources=resources,
      job_name=job_name,
      singularity_str=qsub_sing_envs("", singularity_opts, sing_image),
      cores=cores,
      queue=queue,
      priority=priority,
      # Command to qsub
      shell, script_file, fname, i)

    returned <- system(qsub, intern = T)
    message(returned)
    job_id <- as.numeric(stringr::str_match(returned,"Your job ([0-9]*) ")[,2])
    lv[i, jobid := job_id]
    lv[i, the_qsub := qsub]

  }
  return(list(lv, fname, qsub))
}

## monitor_jobs() ################################################

#' Function to monitor & resubmit jobs if they fail
#'
#' This is intended as a replacement for `waitformodelstofinish()`
#' and its ilk.  This function takes output from `parallelize()`
#' and periodically submits qstat and qacct requests to see how
#' the jobs are running on the cluster.
#'
#' @param parallelize_output output from the `parallelize()` function, where
#'                           `parallelize_output[[1]] is a data.table of loopvars
#'                           and `parallelize_output[[2]] is the filename that
#'                           contains the `save_objs` from parallelize
#' @param sleeptime how long to sleep for between checking the status of the jobs
#'                  running on the cluster? (numeric, in seconds)
#' @param title title for the looping output
#' @param keep_temp_file keep the temp file after this function exits?
#'                       logical; if `keep_temp_file = F` then temp file is deleted
#' @param return_lv should this function return the loopvars?
#' @param max_tries maximum number of times to resubmit a job before giving up
#' @param notification how would you like to be notified when jobs fail?
#'                     The only current option is "pushover" (via `pushover_notify()`)
#'                     but could be expanded to include email, etc. if desired
#' @return loopvars in data table (if `return_lv = T`)
#' @examples
#' # see example for `parallelize()` for workflow
#'
#' # In master script:
#' monitor_jobs(output_from_parallelize,
#'              max_tries = 3,
#'              notification = "pushover")

monitor_jobs <- function(parallelize_output,
                         sleeptime = 100,
                         title = "Job Monitor",
                         keep_temp_file = F,
                         return_lv = F,
                         max_tries = 1,
                         notification = "none") {

  lv <- parallelize_output[[1]]
  fname <- parallelize_output[[2]]

  str_match <- stringr::str_match

  # Wait a minute to let all jobs be submitted
  Sys.sleep(60)

  # Add a column to lv to hold exit statuses
  if (!("exit_status" %in% names(lv))) lv[, exit_status := numeric()]
  if (!("tries" %in% names(lv))) lv[, tries := 1]
  if (!("give_up" %in% names(lv))) lv[, give_up := F]

  # Function for updating of loopvars table
  update_lv_table <- function(lv) {

    # Grab and parse qstat
    get_qstat_table <- function() {
      qs <- system("qstat", intern = T)
      qs <- qs[3:length(qs)] # Trim headers
      qs <- lapply(qs, function(x) gsub("^ *|(?<= ) | *$", "", x, perl = TRUE)) %>% unlist
      qs <- lapply(qs, function(x) unlist(strsplit(x, " ", fixed = T)))
      qs <- lapply(qs, function(x) return(x[1:5])) # Trim to just the useful stuff
      qs <- rbindlist(lapply(qs, function(x) setDT(as.list(x))[]))
      names(qs) <- c("jobid", "prior", "name", "user", "state")
      qs[, jobid:=as.numeric(jobid)]
      return(qs)
    }

    qstat <- get_qstat_table()
    if ("state" %in% names(lv)) lv[, state := NULL] # clear state if present
    lv <- merge(lv, subset(qstat, select = c("jobid", "state")), by = "jobid", all.x = T, all.y = F)

    # For any jobs without an exit status that have closed, grab exit status
    get_qacct_exit_status <- function(jobid) {
      qa <- system(paste0("qacct -j ", jobid), intern = T)
      qa <- str_match(qa, "exit_status\\s+([0-9]*)")[,2]
      qa <- as.numeric(qa[!is.na(qa)])
      return(qa)
    }

    get_qa_wrapper <- function(jobids) {
      Sys.sleep(30) # Give a bit of time to make sure that exit status generated
      return(sapply(jobids, get_qacct_exit_status))
    }

    lv[is.na(state) & is.na(exit_status),
       exit_status := get_qa_wrapper(jobid)]

    # update states
    lv[is.na(state), state := "x"]

    return(lv)
  }

  lv <- update_lv_table(lv)
  n_finished <- nrow(lv[state == "x" & exit_status == 0,])

  while(n_finished < nrow(lv)) {

    lv <- update_lv_table(lv)
    n_finished <- nrow(lv[state == "x" & exit_status == 0,])

    message('\n====================================================================================')
    message(sprintf('==============================      %s      ===============================',title))
    message(paste0('\nAt ',Sys.time(),' .... ',n_finished,' of ', nrow(lv),' jobs have finished.'))
    message("\nJob status:")
    for (i in 1:nrow(lv)) {
      message(paste0("  Job: ", lv[i,'jobname'],
                     " | ID: ", lv[i, 'jobid'],
                     " | Tries: ", lv[i, 'tries'],
                     " | State: ", lv[i, 'state'],
                     " | Exit status: ", lv[i, 'exit_status']))
    }
    message('\n====================================================================================')
    message('====================================================================================')
    message("\n")

    check_if_nonzero_exit <- function(lv, notification) {
      for (i in 1:nrow(lv)) {
        es <- lv[i, 'exit_status']
        if (!is.na(es) & (es != 0)) {

          if (lv[i, 'tries'] <= max_tries - 1) {
            qs <- lv[i, 'the_qsub']
            returned <- system(as.character(qs), intern = T)
            new_job_id <- as.numeric(str_match(returned,"Your job ([0-9]*) ")[,2])
            if (notification == "pushover") {
              pushover_notify(paste0("Resubmitted ", lv[i, 'jobname'], " with new job id ", new_job_id, ". ",
                                     "Failed job id: ", lv[i, 'jobid'], " | exit status: ", lv[i, 'exit_status'], "."),
                              title = paste0("Job failed: ", lv[i, 'jobname']))
            }
            lv[i, 'exit_status'] <- NA
            lv[i, 'jobid'] <- new_job_id
            lv[i, 'tries'] <- lv[i, 'tries'] + 1
          } else if (lv[i, tries] == max_tries) {
            if (lv[i, give_up] == F) {
              if (notification == "pushover") {
                pushover_notify(paste0("Job ", lv[i, 'jobname'], " was resubmitted ", max_tries, " times and will not be resubmitted again.",
                                       "Most recent job id was ", lv[i, 'jobid'], "."))
              }
              lv[i, give_up := TRUE]
            }
          }
        } # close if statement to catch non-zero exit statuses
      } # close for loop over lv rows
      return(lv)
    } # close check_if_nonzero_exit()
    lv <- check_if_nonzero_exit(lv, notification)
    Sys.sleep(sleeptime)
  } # close while loop
  # Exiting function...
  if (return_lv == T) return(lv)
  if (keep_temp_file == F) unlink(paste0("<<<< FILEPATH REDACTED >>>>", fname, ".RData"))
}

## load_from_parallelize() ################################################

#' This function takes the two things passed in a qsub created by
#' `parallelize()` - the temp file name and which row in loopvars the current iteration
#' of the child script should load from - and loads the appropriate `save_objs` and
#' `expand_vars` from `parallelize()` into the environment of the child script.
#'
#' Note that this is meant to be run from the `child script`; by default
#' both `fname` and `rownumber` should be loaded appropriately from
#' `commandArgs()`
#'
#' @param fname filename of the temp file created by `parallelize()`
#' @param rownumber which row of the `lv` object should this particular
#'                  instance of the child script load from
#' @return nothing; assigns objects to child script global environment.
#' @examples
#' # Note: this is within the CHILD SCRIPT, not the master script
#'
#' # Ensure that you've first loaded your functions
#' # (need `load_from_parallelize()` loaded before you can use this)
#' # A good place to put this is right after you've sourced all the
#' # mbg_central function scripts.  Then simply run the
#' # function to set up your environment:
#'
#' load_from_parallelize()
#'

load_from_parallelize <- function(fname = as.character(commandArgs()[4]),
                                    rownumber = as.numeric(commandArgs()[5])) {

  message(paste0("fname: ", fname))
  message(paste0("rownumber: ", rownumber))
  tmp_dir <- "<<<< FILEPATH REDACTED >>>>"

  # Load in the temporary object
  load(paste0(tmp_dir, fname, ".RData"), envir = .GlobalEnv, verbose = T)

  lv <- data.frame(lv)

  # Load loopvars
  this_lv <- lv[rownumber, -which(names(lv) == "jobname"), drop = FALSE]

  # Assign loopvars
  for (n in names(this_lv)) {
    assign(n, this_lv[1, which(names(this_lv) == n)], envir = .GlobalEnv)
    message(paste0(n, ": ", get(n)))
  }
}

## clean_path() --------------------------------------------------------------->
#'
#' @param messy_path a file path (character vector). This assumes that the path
#'    is constructed with one or more '/' as a separator and not '\'
#'
#' @return a 'clean' path with a single '/' as a separator and no trailing '/'
#'
clean_path <- function(messy_path) {
  # convert all '/' to white space separating the diretory names between
  dir_names  <- strsplit(messy_path, '/')[[1]]
  # paste together the diretory names with '/' separators ignoring empty space
  clean_path <- paste(dir_names[which(dir_names != "")], collapse = '/')
  # make sure the tack on any prepended '/' if the original path began with a '/'
  if(substr(messy_path, start = 1, stop = 1) == '/') clean_path <- paste0('/', clean_path)
  return(clean_path)
}

## get_git_status() ----------------------------------------------------------->
#'
#' Given a path to a directory, determine if it is a git repository, and if so,
#' return a string of git hash, branch, date of last commit, and status plus the
#' diff if requested
#'
#' @param repo The git repository directory
#' @param repo_name The name of the repo to print in the header
#' @param show_diff Print the diff
#'
#' @return character vector of a report intended to be printed with cat
#'
get_git_status <- function(repo, repo_name, show_diff = FALSE) {
  # start the report with the name of the repo in a header
  report <- paste("\n\n**********", repo_name, "**********\n", repo,
                  collapse = "\n ")
  # first check and see if this is a git repository and give a warning if it isn't
  if(!dir.exists(paste0(repo, '/.git'))) {
    report <- paste0(c(report,
                       "\n WARNING: 'repo' does not appear to be a git repository",
                       "\n          Cannot print git hash\n"))
  } else {
    # Collect some git commands for the report
    repo_git_cmd <- paste0("cd ", repo, "; git")
    branch       <- system(paste(repo_git_cmd, "branch | grep \'*\' | cut -d ' ' -f2-"),
                           intern = TRUE)
    commit       <- system(paste(repo_git_cmd, "log -1 --pretty=oneline"),
                           intern = TRUE)
    commit_date  <- system(paste(repo_git_cmd, "log -1 --format=%cd"),
                           intern = TRUE)
    status       <- system(paste(repo_git_cmd, "status --long"), intern = TRUE)
    if(show_diff) diff <- system(paste(repo_git_cmd, "diff HEAD"), intern = TRUE)
    # finish constructing the report
    report <- paste(c(report,
                      "\n******* Branch *******", branch,
                      "\n******* Commit *******", commit,
                      "\n** Commit Date/Time **", commit_date,
                      "\n******* Status *******", status, "\n"),
                    collapse = "\n ")
    if(show_diff) report <- paste(c(report,
                                    "********Diff********", diff, "\n"),
                                  collapse = "\n ")
  }
  return(report)
}

## record_git_status() ################################################

#' Retrieve information about current git status (eg, branch, commit,
#' uncommitted changes) for core repository and optionally a separate
#' indicator repository. A check can also be made to see if your core
#' repository is in sync with the LBD core repo
#'
#' Can be used at a minimum to print the git hash of the code being used
#' for posterity.
#'
#' @param core_repo file path to the lbd_core repo
#' @param indic_repo file path to an indicator-specific repo (optional)
#' @param show_diff logical. If there are uncomitted changes, should the
#'     output from git diff be shown?
#' @param check_core_repo logical. Will check whatever has been set as
#'     'core_repo' is the default LBD core code master repo and will give
#'     messages and warnings if not.
#' @param file file path to a text file where output should be saved. This
#'     is optional. If no file path is provided, the output will instead be
#'     printed to the screen.
#'

record_git_status <- function(core_repo,
                              indic_repo = NULL,
                              show_diff = FALSE,
                              check_core_repo = TRUE,
                              file = NULL) {

  # the core code repo directory
  lbd_core_repo <- '<<<< FILEPATH REDACTED >>>>'

  # if a file is specified, start a sink to record output
  if (!is.null(file)) sink(file)

  # print out the core_repo status
  cat(get_git_status(repo = core_repo, repo_name = 'Core repo', show_diff = FALSE))

  # if the user wants to make sure their repo is up-to-date with LBD master
  if(check_core_repo) {
    check_repo_report <- "\n********** REPO CHECK **********\n"
    # The two repo paths are the same
    if(clean_path(core_repo) == lbd_core_repo) {
      check_repo_report <- paste0(c(check_repo_report,
                                    "'core_repo' set to default LBD core code master repo: '",
                                    core_repo, "'\n"))
    } else {
      # Get the git hashes of the standard LBD core code repo and user repo and compare
      lbd_core_repo_hash <- system(paste0('cd ', lbd_core_repo, '; git rev-parse HEAD'),
                                   intern = TRUE)
      core_repo_hash     <- system(paste0('cd ', core_repo, '; git rev-parse HEAD'),
                                   intern = TRUE)
      # The two repos paths are not the same, but the hash matches (separate, up-to-date clones)
      if(lbd_core_repo_hash == core_repo_hash) {
        check_repo_report <- paste0(c(check_repo_report,
                                      "Current 'core_repo' clone is up-to-date with LBD core code master repo:\n'",
                                      lbd_core_repo, "' == '", core_repo, "'\n"))
      # The two repo paths are not the same and the hash doesn't match, repos are out of sync
      } else {
        # print out the LBD core code master repo information for reference against current
        # repo being used as a warning.
        check_repo_report <- paste0(c(check_repo_report,
                                      "WARNING: Current 'core_repo' clone is out of sync with the LBD core code master repo:\n'",
                                       lbd_core_repo, "' != '", core_repo, "'\n",
                                       "\n\n** LBD Core Code MASTER Repo Git Info **\n"))
        check_repo_report <- paste0(c(check_repo_report,
                                      get_git_status(repo = lbd_core_repo,
                                                     repo_name = 'LBD Core Repo',
                                                     show_diff = FALSE)))
      }
    }
    cat(check_repo_report)
  }
  # run git log and git status on the indicator repo
  if (!is.null(indic_repo)) {
    cat(get_git_status(repo = indic_repo, repo_name = 'Indicator repo', show_diff = FALSE))
  }
  # if a file is specified, end the sink recording output
  if (!is.null(file)) sink()
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
#'   the default location of /share/singularity-images/lbd.
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
#'
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
                                      submit_qsubs     = TRUE,
                                      measure          = 'prevalence') {
  
  # Define project first (necessary to validate node options)
  proj <- get_project(proj, use_geo_nodes=geo_nodes)
  
  # Validate arguments
  validate_singularity_options(singularity, singularity_opts)
  validate_node_option(geo_nodes, use_c2_nodes, proj)
  
  # Create sharedir (TODO is this necessary?)
  sharedir = get_model_output_dir(indicator_group, indicator, run_date)
  dir.create(sharedir, showWarnings = FALSE)
  
  # Determine where stdout and stderr files will go
  output_err = setup_log_location(log_dir, user, indic, ig, rd)
  output_log_dir = output_err[[1]]
  error_log_dir = output_err[[2]]
  
  # Define remaining job attributes
  
  run_time <- get_run_time(use_geo_nodes=geo_nodes, use_c2_nodes=use_c2_nodes, queue=queue, run_time=run_time)
  queue <- get_queue(use_geo_nodes=geo_nodes, use_c2_nodes=use_c2_nodes, queue=queue, run_time=run_time)
  
  shell <- paste0(corerepo, '/<<<< FILEPATH REDACTED >>>>/shell_sing.sh')
  sing_image <- get_singularity(image = singularity)
  singularity_str <- qsub_sing_envs("", singularity_opts, sing_image)
  # resources are all the -l qsub arguments
  if(!is.null(cores) & !is.null(slots)) warning("Slots and cores are both specified, cores will be used to assign resources")
  if(is.null(cores)) cores <- slots
  resources <- get_resources(use_geo_nodes=geo_nodes, cores=cores, ram_gb=memory, runtime=run_time)
  
  code <- path_join(corerepo, 'mbg_central', 'share_scripts', 'aggregate_results.R')
  
  qsubs_to_make <- expand.grid(regions, holdouts, ages, raked)
  
  aggregation_qsubs <- make_qsubs_aggregation(qsubs_to_make, error_log_dir, output_log_dir, proj, resources, singularity_str, queue, priority, slots, shell, code,
                                              indicator, indicator_group, run_date, pop_measure, overwrite, corerepo, raking_shapefile_version, modeling_shapefile_version, measure)
  
  if (submit_qsubs) {
    for(qsub in aggregation_qsubs) {
      system(qsub)
    }
  }
  return(aggregation_qsubs)
}

make_qsubs_aggregation <- function(qsubs_to_make, stderr_log, stdout_log, project, resources, singularity_str, queue = NULL, priority = 0, slots, shell, code,
                                   indicator, indicator_group, run_date, pop_measure, overwrite, corerepo, raking_shapefile_version, modeling_shapefile_version, measure) {
  qsubs <- c()
  for (i in 1:nrow(qsubs_to_make)) {
    region <- qsubs_to_make[i, 1]
    holdout <- qsubs_to_make[i, 2]
    age <- qsubs_to_make[i, 3]
    rake <- qsubs_to_make[i, 4]
    shapefile_version <- if(rake) raking_shapefile_version else modeling_shapefile_version
    job_name <- paste(indicator, region, 'aggregate', sep='_')

    
    qsub <- generate_qsub_command(stderr_log=stderr_log,
                                  stdout_log=stdout_log,
                                  project=project,
                                  resources=resources,
                                  job_name=job_name,
                                  singularity_str=singularity_str,
                                  queue=queue,
                                  cores=slots,
                                  priority=priority,
                                  shell, code,
                                  indicator,
                                  indicator_group,
                                  run_date,
                                  rake,
                                  pop_measure,
                                  overwrite,
                                  age,
                                  holdout,
                                  region,
                                  corerepo,
                                  shapefile_version,
                                  measure)
    qsubs <- c(qsubs, qsub)
  }
  return(qsubs)
}



## draw.mesh  ---------------------------------------------------
# this function allows 3d drawing and imaging of a mesh.inla
# object. it can plot either a mesh constructed on the R2 or S2
# manifold. NOTE that it requires the 'rgl' R package and it spawns
# interactive graphics windows. it is untested on the cluster and is
# meant for use on local machines
#
#' @description Plots, in 3d, a s2-manifold or r2-manifold mesh
#'
#' @param mesh an inla.mesh object
#'
#' @param draw.edges Logical. Draw the edges between the vertices?
#'
#' @param draw.segments Logical. Draw the segments that bound the mesh
#'   object?
#'
#' @param draw.plane Logical. Draw a planar shape to aid in displaying
#'   curvature of mesh?
#'
#' @param node.cols Numeric vector with length equal to the number of
#'   mesh vertices (mesh$n). e.g. pass in the posterior mean of the
#'   spatial random effects heights to visualize the fitted GP.
#'
#' @param col.type String taking value of either 'bw' or 'col' and
#'   determining whether the color of the surface should be drawn in
#'   black and white or in color. Only used if a node.cols vector is
#'   passed in.
#'
#' @param window.dims 2d numeric vector describing the width and
#'   height of the plotting in pixels
#'
#' @param plot.mirror Logical. Should a mirror image of the mesh be
#'   added to the plot (could help visualize mesh if printing to a
#'   static image)
#'
#' @param returns nothing but spawns an interactive plot window
#'
#' @examples
#'
#' ## plot a mesh. add color to the background that results from
#' ## the linear interpolation of randomly generated nodal (basis
#' ## height) values
#' draw.s2.mesh(mesh_s, draw.edges = T, draw.segments = T, col.type = 'col',
#'              node.cols = rnorm(n = mesh_s$n), draw.plane = F)
#'
#' ## take a snapshot to save to file
#' fig.path <- '<<<< FILEPATH REDACTED >>>>'
#' rgl.snapshot(file.path(fig.path, "mesh.png"), top=TRUE)
#'
#' ## shut down the graphics window
#' rgl.close()

draw.mesh <- function(mesh, draw.edges=TRUE, draw.segments=TRUE,
                      draw.plane = F, node.cols = NULL,
                      col.type = 'bw', window.dims = c(840, 400),
                      plot.mirror = FALSE){

  require('rgl') ## this is an interactive R graphics. won't work on the cluster

  window.dims = c(50, 50, 50 + window.dims[1], 50 + window.dims[2])

  if(is.null(node.cols)){
    node.cols <- rep(1, mesh$n)
  }

  if(col.type == 'bw')  cp <- function (n,...) { return (grey.colors(n,0.95,0.05,...))}
  if(col.type == 'col') cp <- colorRampPalette(c("darkblue", "blue", "cyan",
                                                 "yellow", "red", "darkred"))

  mesh0 = inla.mesh.create(loc=cbind(0,0), extend=list(offset=1.1,n=4))

  mesh01 = mesh0
  mesh02 = mesh0
  mesh1 = mesh
  mesh2 = mesh
  mesh02$loc[,1] = mesh02$loc[,1]*(-1)
  mesh02$loc[,3] = mesh02$loc[,3]*(-1)
  mesh2$loc[,1] = mesh2$loc[,1]*(-1)
  mesh2$loc[,3] = mesh2$loc[,3]*(-1)

  mesh01$loc[,1] = mesh01$loc[,1]-1.1
  mesh02$loc[,1] = mesh02$loc[,1]+1.1
  mesh1$loc[,1] = mesh1$loc[,1]-1.1
  mesh2$loc[,1] = mesh2$loc[,1]+1.1

  open3d(windowRect=window.dims)
  if(draw.plane){
    plot(mesh01, rgl=TRUE, col="white", color.palette=cp,
         draw.vertices=FALSE, draw.edges=FALSE, add=TRUE)
    if(plot.mirror){
      plot(mesh02, rgl=TRUE, col="white", color.palette=cp,
           draw.vertices=FALSE, draw.edges=FALSE, add=TRUE)
    }
  }
  plot(mesh1, rgl=TRUE, col=node.cols, color.palette=cp,
       draw.vertices=FALSE, draw.edges=draw.edges, add=TRUE,
       draw.segments=draw.segments)
  if(plot.mirror){
    plot(mesh2, rgl=TRUE, col=node.cols, color.palette=cp,
         draw.vertices=FALSE, draw.edges=draw.edges, add=TRUE,
         draw.segments=draw.segments)
  }

  view3d(0,0,fov=0,zoom=0.4)
  rgl.bringtotop()
}



## lonlat3D ---------------------------------------------------
# this function takes in a vector of longitude and a vector of
# latitude and returns coordinates on the S2 sphere (globe living in
# 3D) in (x, y, z) coords
#'
#' @description Converts from long-lat coords to corresponding 3D
#'   coords, (x, y, z), on a sphere with radius 1
#'
#' @param lon numeric vector of longitude coords
#' @param lat numeric vector of latitude coords
#'
#' @return 3 column numeric matrix where each row is a (x,y,z) of the
#'   transformed (long, lat) coords

lonlat3D <- function(lon,lat){
  cbind(cos((lon/180)*pi)*cos((lat/180)*pi),
        sin((lon/180)*pi)*cos((lat/180)*pi),
        sin((lat/180)*pi))
}


## get_output_regions  ---------------------------------------------------
# this function takes in a directory where mbg modeling has stored
# outputs (*cell_pred* objects) and infers the regions specirfied in
# the model
#'
#' @description Determines modeling regions from written output dir objects
#'
#' @param in_dir directory path containing completed mbg cell_pred objects
#'
#' @return A vector string of region names

get_output_regions <- function(in_dir) {
  return(unique(stringr::str_match(list.files(in_dir, pattern = paste0('_cell_draws_eb_')),
                                   '_cell_draws_[^_]+_[^_]+_(.*)_[0-9].RData')[,2]))
}


## delete_model_ouputs
# Quickly delete the model outputs of a number of model runs will create
# a unique combination for each indicator, indicator_group, run_date combination
#'
#' @description Deletes model outputs from share directories
#'
#' @param indicator_group The indicator group level example "u5m"
#' @param indicator The indicator level example "died"
#' @param run_date The name of your model but often its just a run date
#' @param dryrun Only print the directories that would be deleted dont actually
#' @param ... Additional arguments passed to unlink function
#' delete them

delete_model_outputs <- function(
  indicator_group, indicator, run_date, dryrun=FALSE, ...){
  require(dplyr)
  fileDF <- expand.grid(ind=indicator, ig=indicator_group, rd=run_date) %>%
    mutate(moddir='<<<< FILEPATH REDACTED >>>>')

  if(dryrun){
    print("Will Delete the following directories")
    print(fileDF$moddir)
  }

  else{
    for(i in 1:nrow(fileDF)){
      unlink(fileDF$moddir[i], recursive=T, ...)
      print("Deleted directory")
      print(fileDF$moddir[i])
    }
  }
}

#' @title rank_draws
#' @description function to transform values to ranks by draw
#' @param df data frame with columns for each draw
#' @param 'high' means a high value is rank 1, 'low' means a
#' low value should be rank 1
#' @param columns: vector of column names that contain the draws
#' @export

rank_draws <- function(df, ordr, columns) {
  # load library
  library(data.table)

  # Ensure data.table class
  df <- as.data.table(df)
  # Separate into df with draws and the rest
  df1 <- df[, setdiff(names(df), columns), with = FALSE]
  df2 <- df[, columns, with = FALSE]

  # Transform draws to ranks
  if (ordr == 'high') {
    descending <- TRUE
  } else {
    descending <- FALSE
  }
  df2 <- as.data.table(apply(df2, 2, order, decreasing = descending))

  # Summarize ranks
  df2 <- df2[, `:=`(median=apply(df2, 1, median),
                    lci=apply(df2, 1, quantile, 0.025),
                   uci=apply(df2, 1, quantile, 0.975))]

  # Bind to make a single df
  df <- cbind(df1, df2)
  return(df)

}
