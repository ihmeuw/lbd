# miscellaneous functions


waitformodelstofinish <- function(sleeptime=100,
                                    path =  "<<<< FILEPATH REDACTED >>>>",
                                    rd   = run_date,
                                    lv   = loopvars,
                                    showfiles = TRUE,
                                    showcluster = FALSE){

      n_finished <- length(grep("fin_", list.files(path)))

      lv <- data.table(lv)
      names(lv) <- c("reg", "holdout")
      lv$reg <- as.character(lv$reg)

      lv[, file := paste0(path, "/", "fin__bin0_", reg, "_", holdout)]

      while(n_finished != nrow(lv)){

            n_finished <- length(grep("fin_", list.files(path)))

            message('\n====================================================================================')
            message(sprintf('=====================      Run Date: %s      ======================',rd))
            message(paste0('\nAt ',Sys.time(),' .... ',n_finished,' Models have written output.'))
            if(showfiles){
              message('\nCurrently missing models:')
              for(i in 1:nrow(lv))
                if(file.exists(lv[i, file]) == F)
                  message(paste('Region =',lv[i,1],'| Holdout =',lv[i,2]))
            }
            n_cluster <- system("qstat | grep job_ | wc | awk '{print $1}'", intern = T)
            message(paste0('\nFuthermore, there are still ', n_cluster, ' jobs running on the cluster.'))
            if(showcluster){
              system("qstat -r | grep jobname | grep -oP \"(?<=job_)(.*)\"")
            }
            message('\n====================================================================================')
            message('====================================================================================')
            message("\n")
            Sys.sleep(sleeptime)
      }

      unlink(lv$file) # Clean up by deleting extra files once done with the loop
}

waitforpostesttofinish <- function(sleeptime = 100,
                                   indicator = indicator,
                                   indicator_group = indicator_group,
                                   run_date = run_date,
                                   strata,
                                   showfiles = TRUE) {

path <- "<<<< FILEPATH REDACTED >>>>/temp_post_est/"

lv <- data.table(strata)
names(lv) <- "stratum"

lv[, filename := paste0(path, stratum, "_post_est_list.RData")]

n_finished <- 0
n_needed <- length(strata)

 while (n_finished != n_needed) {
   lv[, exists := file.exists(filename)]
   n_finished <- nrow(lv[exists == T])

   message('\n====================================================================================')
   message('\nRunning post-estimation in parallel')
   message(paste0('Run Date: ', run_date))
   message(paste0('\nAt ',Sys.time(),' .... ',n_finished,' out of ', n_needed, ' strata have written output.'))
   if (showfiles) {
    message(paste0('\nCurrently missing strata:'))
    for(i in 1:nrow(lv[exists == F])) {
      message(subset(lv, exists == F)[i, stratum])
    }
   }
   message('\n====================================================================================')
   message('====================================================================================')
   Sys.sleep(sleeptime)
 }
}


waitforaggregation <- function(rd = run_date,
                               indic = indicator,
                               ig = indicator_group,
                               ages,
                               regions,
                               holdouts,
                               raked,
                               dir_to_search = NULL) {

  # waitformodelstofinish() analog for aggregation

  if (is.null(dir_to_search)) {
    dir_to_search <- "<<<< FILEPATH REDACTED >>>>"
  }

  lv <- expand.grid(regions, holdouts, ages, raked) %>% as.data.table
  names(lv) <- c("reg", "holdout", "age", "raked")
  lv[, file := paste0(dir_to_search, "fin_agg_", reg, "_", holdout, "_", age, "_", raked)]
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
                       "Holdout: ", lv_left[i, holdout], " | ",
                       "Age: ", lv_left[i, age], " | ",
                       "Raked: ", lv_left[i, raked]))
      }
    }

    Sys.sleep(60)

  }
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
                                merge_hierarchy_list = F) {

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
    dir_to_search <- "<<<< FILEPATH REDACTED >>>>"
  }

  message("Combining aggregation results...")

   for(rake in raked) {
    for (holdout in holdouts) {
      for (age in ages) {
        message(paste0("\nWorking on age: ", age, " | holdout: ", holdout, " | raked: ", rake))

        # Set up lists
        ad0 <- list()
        ad1 <- list()
        ad2 <- list()
        sp_h <- list()

        for (reg in regions) {
          message(paste0("  Region: ", reg))
          load(paste0(dir_to_search, indic, "_", ifelse(rake, "raked", "unraked"),
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

          ad0[[reg]] <- admin_0
          ad1[[reg]] <- admin_1
          ad2[[reg]] <- admin_2
          sp_h[[reg]] <- sp_hierarchy_list

          rm(admin_0, admin_1, admin_2, sp_hierarchy_list)
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
                           ifelse(rake, "raked", "unraked"),
                           "_admin_draws_eb_bin", age, "_",
                           holdout, ".RData"))
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
#' file defined are assumed to exist as the default location
#' In either case it will test to make sure the file exists and exit if it does
#' not. The default image is hardcoded into the shell script used to launch
#' Singularity containers:
#'   lbd_core/mbg_central/share_scripts/shell_sing.sh
#'
#' @param image A string that defines which Singularity image to launch
#'   [default = 'default']. If the 'default' keyword is passed in or left blank,
#'   the default keyword will be returned. Either the full path to the image
#'   may be provided or only the Singularity image name. In the latter case,
#'   the image is assumed to live in the default image location:
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
    sing_image <- "<<<< FILEPATH REDACTED >>>>"
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

# make_qsub --------------------------------------------------------------------
#'
#' Constructs a qsub string and returns it
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
#'   the default location
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
#'   [default = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1)] (TEMPORARY UNTIL
#'   MCLAPPLY WORKAROUND IN PLACE. ORIGINAL: [default = NULL])
#'
#' @return Returns a qsub string
#'
make_qsub <- function(user                     = Sys.info()['user'],
                      code,
                      cores                    = slots,
                      memory                   = 100,
                      proj                     = 'proj_geospatial',
                      ig                       = indicator_group,
                      indic                    = indicator,
                      reg                      = "test",
                      age                      = 0,
                      vallevel                 = '',
                      hard_set_sbh_wgt         = TRUE,
                      pct_sbh_wgt              = 100,
                      rd                       = run_date,
                      log_location             = 'sgeoutput',
                      saveimage                = FALSE,
                      test                     = FALSE,
                      holdout                  = 0,
                      use_base_covariates      = FALSE,
                      test_pv_on_child_models  = TRUE,
                      constrain_children_0_inf = FALSE,
                      child_cv_folds           = 5,
                      fit_stack_cv             = TRUE,
                      shell                    = "r_shell.sh",
                      modeltype                = 'full',
                      corerepo                 = core_repo,
                      geo_nodes                = FALSE,
                      use_c2_nodes             = FALSE,
                      singularity              = 'default',
                      singularity_opts         = NULL) {

  if(test) t=1 else t=0
  if(use_base_covariates) b=1 else b=0
  if(test_pv_on_child_models) pvc=1 else pvc=0
  if(constrain_children_0_inf) cc0i=1 else cc0i=0
  if(fit_stack_cv) fscv=1 else fscv=0
  if(hard_set_sbh_wgt) hssw=1 else hssw=0

  if(saveimage==TRUE)  save.image(paste0('<<<< FILEPATH REDACTED >>>>/pre_run_tempimage_', rd, '_bin',age,'_',reg,'_',holdout,'.RData'))

  dir.create("<<<< FILEPATH REDACTED >>>>")

  sharedir <- "<<<< FILEPATH REDACTED >>>>"

  if(log_location=='sgeoutput')
    logloc = "<<<< FILEPATH REDACTED >>>>"
  if(log_location=='sharedir'){
    logloc = "<<<< FILEPATH REDACTED >>>>"
    dir.create(sprintf('%s/errors',logloc), showWarnings = F)
    dir.create(sprintf('%s/output',logloc), showWarnings = F)
  }

 # Submit to prod or geo nodes. The geos nodes only ever use the project
 # "proj_geo_nodes" and require the "-l geos_node=TRUE". When submitting to the
 # "prod" cluster, we want to be able to use different projects with the default
 # being "project_geospatial".
  if(geo_nodes) {
    node.flag <- "-l geos_node=TRUE" ## send the job to geo nodes
    # "proj" is set to whatever is passed in
  } else {
    if(use_c2_nodes) node.flag <- "-q all.q@@c2-nodes" else node.flag <- ""
  }

  # If the user has passed in options for the Singularity container in the
  # 'singularity_opts' argument, but the 'singularity' argument is 'NULL' exit
  # with an error.
  if(is.null(singularity) & !is.null(singularity_opts)) {
    message("'singularity' argument is 'NULL' but options passed in for 'singularity_opts'")
    stop("Exiting!")
  }

  # If the script is to be run with R from a Singularity container, point to
  # the shell script to launch the container. Users can provide the 'default'
  # keyword to get the default Singulariy image, just the name of the image
  # located at the default path, or the full path to the image.
  if(!is.null(singularity)) {
    shell <- paste0(corerepo, '/mbg_central/share_scripts/shell_sing.sh')
    # Determine which Singularity image to use:
    sing_image <- get_singularity(image = singularity)
    # TEMPORARY DEFAULT: CAN REMOVE AFTER MCLAPPLY WORKAROUND IN PLACE
    if(is.null(singularity_opts)) singularity_opts <- list(SET_OMP_THREADS=1, SET_MKL_THREADS=1)
  } else {
    # if not, fire up the appropriate version of R depending on the cluster node
    shell <- paste0(corerepo, '/mbg_central/share_scripts/shell_cluster.sh')
  }

  # Piece together lengthy `qsub` command
  qsub <- paste0("qsub",
                 " -e ", logloc, "/errors",
                 " -o ", logloc, "/output",
                 " -cwd -l mem_free=", memory, "G",
                 " -pe multi_slot ", cores,
                 " -P ", proj, " ", node.flag)

  # If a Singularity image is being used, pass the name of the image from
  # `get_singularity` as well as any other environmentals the user asked for
  # from the 'singularity_opts' argument to the qsub job
  if(!is.null(singularity)) qsub <- qsub_sing_envs(qsub, singularity_opts,
                                                   sing_image)

  # append the rest of the qsub command
  qsub <- paste0(qsub,
                 " -N job_", reg, "_", age, "_", holdout, vallevel,
                 " ", ig, "/", shell, " ", ig, "/", code, ".R")
  qsub <- paste(qsub,
                reg, age, rd, t, holdout, indic, ig, b, pvc, cc0i,
                child_cv_folds, fscv, modeltype, vallevel, pct_sbh_wgt, hssw,
                sep = " ")

  return(qsub)
}

# make_qsub_share --------------------------------------------------------------
#'
#' Constructs a qsub string and returns it
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
#'   the default location
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
#'   [default = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1)] (TEMPORARY UNTIL
#'   MCLAPPLY WORKAROUND IN PLACE. ORIGINAL: [default = NULL])
#'
#' @return Returns a qsub string
#'
make_qsub_share <- function(user             = Sys.info()['user'],
                            code             = NULL,
                            cores            = slots,
                            memory           = 100,
                            proj             = 'proj_geospatial',
                            ig               = indicator_group,
                            indic            = indicator,
                            reg              = "test",
                            age              = 0,
                            rd               = run_date,
                            log_location     = 'sharedir',
                            addl_job_name    = '',
                            saveimage        = FALSE,
                            test             = FALSE,
                            holdout          = 0,
                            corerepo         = core_repo,
                            geo_nodes        = FALSE,
                            use_c2_nodes     = FALSE,
                            new_geos_nodes   = NULL,
                            singularity      = 'default',
                            singularity_opts = NULL) {

  # save an image
  if(saveimage==TRUE) save.image(paste0('<<<< FILEPATH REDACTED >>>>/pre_run_tempimage_', rd, '_bin',age,'_',reg,'_',holdout,'.RData'))

  # sort directories
  sharedir <- "<<<< FILEPATH REDACTED >>>>"

  if(log_location=='sgeoutput')
    logloc = "<<<< FILEPATH REDACTED >>>>"
  if(log_location=='sharedir'){
    logloc = "<<<< FILEPATH REDACTED >>>>"
    dir.create(sprintf('%s/errors',logloc), showWarnings = FALSE)
    dir.create(sprintf('%s/output',logloc), showWarnings = FALSE)
  }

  # Submit to prod or geo nodes. The geos nodes only ever use the project
  # "proj_geo_nodes" and require the "-l geos_node=TRUE". When submitting to the
  # "prod" cluster, we want to be able to use different projects with the default
  # being "project_geospatial".
  if(geo_nodes) {
    node.flag <- "-l geos_node=TRUE" ## send the job to geo nodes

    # Specify new/old geos nodes if desired
    if (!is.null(new_geos_nodes)) {
      if (new_geos_nodes) {
        # Specify "new" geos nodes, less memory faster computation
        node.flag <- paste0(node.flag, " -q geospatial.q@@intel_gold_6130_756gb")
      } else {
        # Specify "old" nodes
        node.flag <- paste0(node.flag, " -q geospatial.q@@intel_e5_2660v4_1024gb")
      }
    }
  } else {
    if(use_c2_nodes) node.flag <- "-q all.q@@c2-nodes" else node.flag <- ""
  }

  # If the user has passed in options for the Singularity container in the
  # 'singularity_opts' argument, but the 'singularity' argument is 'NULL' exit
  # with an error.
  if(is.null(singularity) & !is.null(singularity_opts)) {
    message("'singularity' argument is 'NULL' but options passed in for 'singularity_opts'")
    stop("Exiting!")
  }

  # If the script is to be run with R from a Singularity container, point to
  # the shell script to launch the container. Users can provide the 'default'
  # keyword to get the default Singulariy image, just the name of the image
  # located at the default path, or the full path to the image.
  if(!is.null(singularity)) {
    shell <- paste0(corerepo, '/mbg_central/share_scripts/shell_sing.sh')
    # Determine which Singularity image to use:
    sing_image <- get_singularity(image = singularity)
    # TEMPORARY DEFAULT: CAN REMOVE AFTER MCLAPPLY WORKAROUND IN PLACE
    if(is.null(singularity_opts)) singularity_opts <- list(SET_OMP_THREADS=1, SET_MKL_THREADS=1)
  } else {
    # if not, fire up the appropriate version of R depending on the cluster node
    shell <- paste0(corerepo, '/mbg_central/share_scripts/shell_cluster.sh')
  }

  # set code to shared paralel if its NULL
  if(is.null(code)) {
    code <- sprintf('%s/mbg_central/share_scripts/parallel_model.R',corerepo)
  } else {
    code <- sprintf('%s/%s/%s.R',corerepo,ig,code)
  }

  # Piece together lengthy `qsub` command
  qsub <- paste0("qsub",
                 " -e ", logloc, "/errors",
                 " -o ", logloc, "/output",
                 " -cwd -l mem_free=", memory, "G",
                 " -pe multi_slot ", cores,
                 " -P ", proj, " ", node.flag)

  # If a Singularity image is being used, pass the name of the image from
  # `get_singularity` as well as any other environmentals the user asked for
  # from the 'singularity_opts' argument to the qsub job
  if(!is.null(singularity)) qsub <- qsub_sing_envs(qsub, singularity_opts,
                                                   sing_image)

  # append the rest of the qsub command
  qsub <- paste0(qsub,
                 " -N job_", addl_job_name, "_", reg, "_", age, "_", holdout)
  qsub <- paste(qsub,
                shell, code, reg, age, rd, as.numeric(test), holdout, indic, ig, "fin",
                sep = " ")

  return(qsub)
}

# make_qsub_postest -------------------------------------------------------------
#'
#' Constructs a qsub string for the post_estimation script and returns it
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
#'   the default location
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
#'   [default = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1)] (TEMPORARY UNTIL
#'   MCLAPPLY WORKAROUND IN PLACE. ORIGINAL: [default = NULL])
#'
#' @return Returns a qsub string
#'
make_qsub_postest <- function(user             = Sys.info()['user'],
                              code,
                              cores            = slots,
                              memory           = 100,
                              proj             = 'proj_geospatial',
                              ig               = indicator_group,
                              indic            = indicator,
                              stratum          = "test",
                              rd               = run_date,
                              log_location     = 'sgeoutput',
                              script_dir       = NULL,
                              keepimage        = FALSE,
                              corerepo         = core_repo,
                              geo_nodes        = FALSE,
                              new_geos_nodes   = NULL,
                              use_c2_nodes     = FALSE,
                              singularity      = 'default',
                              singularity_opts = NULL) {

  sharedir <- "<<<< FILEPATH REDACTED >>>>"
  temp_dir <- paste0(sharedir, 'temp_post_est/')
  dir.create(temp_dir, showWarnings = F)

  if(is.null(script_dir)) script_dir <- paste0(corerepo, '/mbg_central/share_scripts')

  if(log_location=='sgeoutput')
    logloc = "<<<< FILEPATH REDACTED >>>>"
  if(log_location=='sharedir'){
    logloc = sharedir
    dir.create(sprintf('%serrors',logloc), showWarnings = F)
    dir.create(sprintf('%soutput',logloc), showWarnings = F)
  }

  # Submit to prod or geo nodes. The geos nodes only ever use the project
  # "proj_geo_nodes" and require the "-l geos_node=TRUE". When submitting to the
  # "prod" cluster, we want to be able to use different projects with the default
  # being "project_geospatial".
  if(geo_nodes) {
    node.flag <- "-l geos_node=TRUE" ## send the job to geo nodes
    # "proj" is set to whatever is passed in
    # Specify new/old geos nodes if desired
    if (!is.null(new_geos_nodes)) {
      if (new_geos_nodes) {
        # Specify "new" geos nodes, less memory faster computation
        node.flag <- paste0(node.flag, " -q geospatial.q@@intel_gold_6130_756gb")
      } else {
        # Specify "old" nodes
        node.flag <- paste0(node.flag, " -q geospatial.q@@intel_e5_2660v4_1024gb")
      }
    }
  } else {
    if(use_c2_nodes) node.flag <- "-q all.q@@c2-nodes" else node.flag <- ""
  }

  gn <- as.character(geo_nodes)

  # If the user has passed in options for the Singularity container in the
  # 'singularity_opts' argument, but the 'singularity' argument is 'NULL' exit
  # with an error.
  if(is.null(singularity) & !is.null(singularity_opts)) {
    message("'singularity' argument is 'NULL' but options passed in for 'singularity_opts'")
    stop("Exiting!")
  }

  # If the script is to be run with R from a Singularity container, point to
  # the shell script to launch the container. Users can provide the 'default'
  # keyword to get the default Singulariy image, just the name of the image
  # located at the default path, or the full path to the image.
  if(!is.null(singularity)) {
    shell <- paste0(corerepo, '/mbg_central/share_scripts/shell_sing.sh')
    # Determine which Singularity image to use:
    sing_image <- get_singularity(image = singularity)
    # TEMPORARY DEFAULT: CAN REMOVE AFTER MCLAPPLY WORKAROUND IN PLACE
    if(is.null(singularity_opts)) singularity_opts <- list(SET_OMP_THREADS=1, SET_MKL_THREADS=1)
  } else {
    # if not, fire up the appropriate version of R depending on the cluster node
    shell <- paste0(corerepo, '/mbg_central/share_scripts/shell_cluster.sh')
  }

  # Piece together lengthy `qsub` command
  # last line below sets args for parallel script
  qsub <- paste0("qsub",
                " -e ", logloc, "errors",
                " -o ", logloc, "output",
                " -cwd -l mem_free=", memory, "G",
                " -pe multi_slot ", cores,
                " -P ", proj, " ", node.flag)

  # If a Singularity image is being used, pass the name of the image from
  # `get_singularity` as well as any other environmentals the user asked for
  # from the 'singularity_opts' argument to the qsub job
  if(!is.null(singularity)) qsub <- qsub_sing_envs(qsub, singularity_opts,
                                                   sing_image)

  # append the rest of the qsub command
  qsub <- paste0(qsub,
                 " -N job_pe_", indic, "_", stratum, " ",
                shell, " ",
                script_dir, "/", code, ".R")
  qsub <- paste(qsub,
                stratum, rd, indic, ig, gn,
                sep = " ")

  return(qsub)
}

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

check_config <- function(cr = core_repo) {
  must_haves <- read.csv(paste0(cr, '/mbg_central/share_scripts/common_inputs/config_must_haves.csv'), header = F, stringsAsFactors = F)$V1

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

    } else if (confs == "rho_prior") {
      message("You are missing a 'rho_prior' argument in your config. Defaulting to 'list(prior = 'normal', param = c(0, 0.1502314))'")
      rho_prior <<- "list(prior = 'normal', param = c(0, 1/(2.58^2)))"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "use_nid_res") {
      message("You are missing a 'use_nid_res' argument in your config. Defaulting to FALSE")
      use_nid_res <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else {
      stop(paste0(confs,' is missing, add it to your config'))
    }
  }

  message("\nAdditional covariates: ")
  extras <- config$V1[!(config$V1 %in% must_haves)]
  for (extra in extras) message(paste0('  ', extra, ': ', get(extra)))
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
#' if a fork is being used.
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
  lbd_core_repo <- '<<<< FILEPATH REDACTED >>>>/lbd_core'

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
