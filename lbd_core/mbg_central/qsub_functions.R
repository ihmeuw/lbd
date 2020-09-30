# Constants. These should probably live elsewhere.
MBG_ROOT <- "<<<< FILEPATH REDACTED >>>>"
SGE_OUTPUT_ROOT <- "<<<< FILEPATH REDACTED >>>>"

## path_join
#'
#' @title join values together using filesystem path separator.
#'
#' @param ... values to join.
#'
#' @return path built from arguments.
#'
path_join <- function(...) {
  result <- paste(..., sep = .Platform$file.sep)
  if (length(result) == 0) stop("no arguments provided to path_join")
  return(result)
}

## get_indicator_dir
#'
#' @title compute the directory for outputs for a given indicator
#'
#' @param indicator_group the indicator group e.g., "hiv" or "u5m".
#'
#' @param indicator the indicator e.g., "wasting_mod_b" or "male_circumcision".
#'
#' @return filesystem path.
get_indicator_dir <- function(indicator_group, indicator) {
  path_join(MBG_ROOT, indicator_group, indicator)
}

## get_model_output_dir
#'
#' @title return path to model output directory.
#'
#' @param indicator_group the indicator group e.g., "hiv" or "u5m".
#'
#' @param indicator the indicator e.g., "wasting_mod_b" or "male_circumcision".
#'
#' @param run_date string run date in YYYY_MM_DD_HH_MM_SS format.
#'
#' @return the path to the models output directory.
#'
get_model_output_dir <- function(indicator_group, indicator, run_date) {
  path_join(get_indicator_dir(indicator_group, indicator), "output", run_date)
}

## sge_output_dir
#'
#' @title return the directory where sgeoutput is usually kept.
#'
#' @param user the name of the user.
#'
#' @return the root directory for sge error/output logs.
#'
get_sge_output_dir <- function(user) {
  path_join(SGE_OUTPUT_ROOT, user)
}

## pre_run_image_path
#'
#' @title return the path to save model data prior to running parallel script.
#'
#' @param indicator_group the indicator group e.g., "hiv" or "u5m".
#'
#' @param indicator the indicator e.g., "wasting_mod_b" or "male_circumcision".
#'
#' @param run_date string run date in YYYY_MM_DD_HH_MM_SS format.
#'
#' @param age the age for the model which uses this data.
#'
#' @param region str region name.
#'
#' @param holdout numeric holdout value (0 for no holdout).
#'
#' @return the path to save the file.
#'
pre_run_image_path <- function(indicator_group, indicator, run_date, age, region, holdout) {
  f <- sprintf("pre_run_tempimage_%s_bin%s_%s_%s.RData", run_date, age, region, holdout)
  path_join(get_indicator_dir(indicator_group, indicator), "model_image_history", f)
}


## setup_log_location
#'
# " @title Determines error/output file path and optionally ensures dirs exist.
#'
#' @param log_location the location of the log files OR 'sgeoutput' OR
#'  'sharedir'. 'sgeoutput' and 'sharedir' will appropriately substitute the
#'  user's default sgeoutput directory or the models default output directory
#'  (respectively).
#'
#' @param user the name of the user.
#'
#' @param indicator the indicator e.g., "wasting_mod_b" or "male_circumcision".
#'
#' @param run_date string run date in YYYY_MM_DD_HH_MM_SS format.
#'
#' @param age the age for the model which uses this data.
#'
#' @param make_dirs whether to make the directories. [default = TRUE]
#'
setup_log_location <- function(log_location, user, indicator, indicator_group, run_date, make_dirs = TRUE) {
  if (log_location == "sgeoutput") {
    log_root <- get_sge_output_dir(user)
  } else if (log_location == "sharedir") {
    log_root <- get_model_output_dir(indicator_group, indicator, run_date)
  } else {
    log_root <- log_location
  }
  output_err <- c(path_join(log_root, "output"), path_join(log_root, "errors"))
  if (make_dirs) sapply(output_err, dir.create, showWarnings = FALSE)
  return(output_err)
}

## validate_node_option
#'
#' @title Validates node selection emitting warning message if inconsistent.
#'
#' @param use_geo_nodes logical indicating whether to use geospatial nodes.
#'
#' @param use_c2_nodes logical indicating whether to use "c2-" prefixed nodes.
#'
validate_node_option <- function(use_geo_nodes, use_c2_nodes, project) {
  if (use_geo_nodes & use_c2_nodes) {
    message("WARNING: Both 'geo_nodes' and 'use_c2_nodes' arguments were set to TRUE")
    message(paste0("         Submitting job to LBD nodes under project: '", project, "'"))
  }
}

## validate_singularity_options
#'
#' @title Validates singularity-related options and errors if invalid.
#'
#' @param singularity the string name of the singularity image to use.
#'
#' @param singularity_opts named list of options.
#'
validate_singularity_options <- function(singularity, singularity_opts) {
  if (is.null(singularity) & !is.null(singularity_opts)) {
    message("'singularity' argument is 'NULL' but options passed in for 'singularity_opts'")
    stop("Exiting!")
  }
}

## get_project
#'
#' @title Return default project if not set.
#'
#' @param proj str project name (may be NULL)
#'
#' @param use_geo_nodes logical indicating whether to use geospatial nodes.
#'
#' @return proj (if not null) or an appropriate default project.
#'
get_project <- function(proj, use_geo_nodes) {
  if (!is.null(proj)) return(proj)
  if (use_geo_nodes) "proj_geo_nodes" else "proj_geospatial"
}

#' @title Complete user-supplied runtime string into DHMS format
#'
#' @description Check the input string of runtime and convert its values into a DD:HH:MM:SS format
#'
#' @param run_time String. User supplied run-time in format \code{[DD]:[HH]:[MM]:[SS]}
#'
#' @return The same run-time supplied by user once validations have passed,
#' but with days appended if not supplied, in format \code{DD:HH:MM:SS}
#'
#' @importFrom stringr str_split
#' @export
complete_runtime_string <- function(run_time) {

  ## Split runtime into its components
  ## We could have 3 (H:M:S) or 4 (D:H:M:S) components here
  rt_split <- as.numeric(stringr::str_split(run_time, ":")[[1]])


  ## Return a 'complete' timestamp
  ## with 'DD' added if the validations are correct

  ## DD:HH:MM:SS
  if (length(rt_split) == 4) {
    return(run_time)
  } else if (length(rt_split) == 3) {
    return(paste0("00:", run_time))
  } else if (length(rt_split) == 2) {
    return(paste0("00:00:", run_time))
  } else if (length(rt_split) == 1) {
    return(paste0("00:00:00:", run_time))
  } else {
    stop("Run-time is not in correct format. Exiting.")
  }
}


#' @title Run-time to hours
#'
#' @description Convert a validated \code{DD:HH:MM:SS} runtime string to hours
#'
#' @param run_time String. Run-time in format \code{DD:HH:MM:SS}.
#'
#' @return number of hours
#'
#' @export
dhms_to_hours <- function(run_time) {

  ## Split into numeric vector
  rt_split <- as.numeric(stringr::str_split(
    run_time, ":"
  )[[1]])

  ## Convert to hours space
  return(rt_split[1]*24 + rt_split[2]*1 + rt_split[3]/60 + rt_split[4]/3600)

}


#' @title Validate user-supplied queue string
#'
#' @description Check the input string of queue and validate whether its a legitimate queue or not.
#' This checks for whether the queue values are matching with \code{all}, \code{geospatial} or \code{long}
#'
#' @param queue String. User supplied queue.
#'
#' @return the queue post-validating
#'
#' @export
validate_queue_string <- function(queue) {
  stopifnot(queue %in% c("all.q", "geospatial.q", "long.q"))
  return(queue)
}



#' @title Get max runtime of queue
#'
#' @description Check the input string of queue and the maximum runtime available for that queue
#'
#' @param queue String. User supplied queue
#'
#' @return a string of run-time in format \code{DD:HH:MM:SS}
#'
#' @export
get_max_runtime_by_queue <- function(queue) {
  if (queue %like% "all") {
    return("03:00:00:00")
  } else if (queue %like% "geospatial") {
    return("25:00:00:00")
  } else {
    return("16:00:00:00")
  }
}


#' @title Get run-time
#'
#' @description Get run-time based on user supplied queue, run-time, and whether to use geos nodes or c2 nodes
#' If run-time is already supplied by the user, then we just validate that string and return it
#'
#' @param use_geo_nodes Boolean. Use geo nodes?
#' @param use_c2_nodes Boolean. Use c2 nodes? Only relevant to cluster-prod.
#' @param run_time String. User supplied run-time
#' @param queue String. User supplied queue.
#'
#' @return A run-time string in format \code{DD:HH:MM:SS} if using fair cluster. Returns \code{NULL} otherwise.
#' If run-time is empty, then it return the max runtime
#'
#' @export
get_run_time <- function(use_geo_nodes, use_c2_nodes, queue = NULL, run_time = NULL) {

  ## Get legacy cluster queue:
  if (!is_new_cluster()) {
    return(NULL)
  } else {
    ## If RT is null:
    if (is.null(run_time)) {
      if (is.null(queue)) {
        if(!use_geo_nodes | is.null(use_geo_nodes)) {
          message("No queue supplied. Reverting to long.q")
        }
        return(get_max_runtime_by_queue("long.q"))
      } else {
        if (use_geo_nodes) {
          if (!(queue %like% "geospatial")) {
            warning("You did not specify geospatial.q but supplied TRUE for use_geo_nodes. Reverting to the geospatial.q")
          }
          run_time <- get_max_runtime_by_queue("geospatial.q")
        } else {
          if (queue %in% "all.q") {
            run_time <- get_max_runtime_by_queue("all.q")
          } else {
            run_time <- get_max_runtime_by_queue("long.q")
          }
        }
      }
    } else {
      ## Nothing to do except just validate run_time and return below!
    }
    return(complete_runtime_string(run_time))
  }
}



#' @title Get queue
#'
#' @description Get queue based on user supplied queue, run-time, and whether to use geos nodes or c2 nodes
#' If queue is already supplied by the user, then we just validate that string and return it
#'
#' @param use_geo_nodes Boolean. Use geo nodes?
#' @param use_c2_nodes Boolean. Use c2 nodes? Only relevant to cluster-prod.
#' @param run_time String. User supplied run-time
#' @param queue String. User supplied queue.
#'
#' @return A queue string.
#'
#' @importFrom stringr str_split
#' @export
get_queue <- function(use_geo_nodes, use_c2_nodes, queue = NULL, run_time = NULL) {

  ## Legacy cluster chunk:
  if (!is_new_cluster()) {
    if (use_geo_nodes) {
      return("geospatial.q")
    } else if (use_c2_nodes) {
      return("all.q@@c2-nodes")
    } else {
      return("all.q")
    }
  } else {

    ## All about fair cluster now

    ## If queue is supplied, then validate that string and return it:
    if (!is.null(queue)) {

      ## If use_geo_nodes is TRUE, then we always want to return geospatial.q
      if(use_geo_nodes) {
        message("use_geo_nodes was found to be TRUE. Returning geospatial.q regardless of the queue supplied.")
        return("geospatial.q")
      }

      return(validate_queue_string(queue))
    } else {
      ## If queue is null, then infer the queue based on the booleans and
      ## supplied runtime


      ## Split out the validated runtime into hours
      rt_hrs <- dhms_to_hours(
        get_run_time(use_geo_nodes = use_geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
      )


      ## Use the booleans to derive the proper queue:
      if (use_geo_nodes) {
        if (rt_hrs > 25*24) {
          stop("No queues exist with run-time greater than 25 days. Exiting.")
        }
        return("geospatial.q")
      } else {
        if (rt_hrs > 16*24) {
          stop("You are asking for a run-time of greater than 16 days but not asking for geospatial.q. Exiting.")
        } else if (rt_hrs <= 3*24) {
          return("all.q")
        } else {
          return("long.q")
        }
      }
    }
  }
}



#' @title Is this the new cluster?
#'
#' @description Returns logical indicating if program is running on the new cluster.
#'
#' @return logical indicating if program is on the new cluster.
#'
#' @export
is_new_cluster <- function() {
  if (Sys.getenv("SGE_ENV") == "" & Sys.getenv("RSTUDIO_HTTP_REFERER") == "") {
    stop("Neither SGE_ENV nor RSTUDIO_HTTP_REFERER are found. Unidentifiable environment, exiting.")
  }

  if (Sys.getenv("SGE_ENV") == "") {
    warning("SGE_ENV is not set")

    # NOTE: If this is RStudio, then SGE_ENV isn't inherited, hence
    # we check using the URL
    is_new_hack <- grepl("-uge-", Sys.getenv("RSTUDIO_HTTP_REFERER"))
    is_old_hack <- grepl("cn|c2", Sys.getenv("RSTUDIO_HTTP_REFERER"))
    if (is_new_hack) {
      warning("Running RStudio on the fair cluster")
      return(TRUE)
    } else if (is_old_hack) {
      warning("Using... the old cluster?")
      return(FALSE)
    } else {
      stop("URL does not either have 'uge' or start with 'cn' or 'c2'")
    }
  } else {
    # per @matpp in #new-cluster-rollout at 9:02PM on 25 October 2018
    return(grepl("-el7$", Sys.getenv("SGE_ENV")))
  }
}

## get_resources
#'
#' @title Get list of resources to request for new cluster.
#'
#' @param use_geo_nodes logical in dicating if code should run on LBD nodes.
#'
#' @param cores numeric indicating the number of cores to request.
#'
#' @param ram_gb numeric indicating the number of GB of RAM to request.
#'
#' @param runtime (default NULL) runtime in 'DD:HH:MM:SS' format.
#'
#' @return named vector with resource name and string values.
#'
get_resources <- function(..., use_geo_nodes, cores, ram_gb, runtime = NULL) {
  result <- c()
  if (is_new_cluster()) {
    result <- c(result, c(
      archive = "TRUE",
      m_mem_free = paste0(ram_gb, "G"),
      fthread = cores
    ))
    if (is.null(runtime)) {
      if (use_geo_nodes) {
        runtime <- "25:00:00:00"
      } else {
        runtime <- "3:00:00:00"
      }
    } else {
      runtime <- get_run_time(use_geo_nodes, use_c2_nodes, run_time = runtime)
    }
    result <- c(result, c(h_rt = runtime))
  } else {
    result <- c(result, c(mem_free = paste0(ram_gb, "G")))
    if (use_geo_nodes) result <- c(result, c(geos_node = "TRUE"))
  }
  return(result)
}

## get_resource_str
#'
#' @title Compute a string of resources to request for a qsub command.
#'
#' @param resources named vector of resources.
#'  @seealso \code{\link{get_resources}}
#'
get_resource_str <- function(resources) {
  if (length(resources) == 0) return("")
  paste0("-l ", names(resources), "=", resources, collapse = " ")
}

## generate_qsub_command
#'
#' @title Generates a qsub command that should work on any IHME cluster.
#'
#' @param ... The command + arguments to be run as a hibh e.g., the
#'  singularity shell, the parallel script and its 8 arguments.
#'
#' @param stderr_log string indicating directory or file to log stderr to.
#'
#' @param stdout_log string indicating directory or file to log stdout to.
#'
#' @param project the project to run the job under.
#'
#' @param resources named vector of resources.
#'  @seealso \code{\link{get_resources}}
#'
#' @param job_name string name for the job.
#'
#' @param singularity_opts string options to pass to singularity.
#'  @seealso \code{\link{qsub_sing_envs}}
#'
#' @param cores numeric number of cores to use
#'
#' @param priority Job priority that can be deprioritized if needed, and can only be used for values in [-1023,0]. Default = 0.
#' This value will get bounded to 0 or -1023 if the user supplies a value outside those bounds.
#'
generate_qsub_command <- function(..., stderr_log, stdout_log, project, resources, job_name, singularity_str, queue, cores = NULL, priority = 0) {

  # Check on priority value
  if(priority < -1023) {
    warning("Priority value supplied is < -1023. Setting to -1023 (lowest possible value).")
    priority <- -1023
  } else if(priority > 0) {
    warning("Priority value supplied is > 0 Setting to 0 (highest possible value).")
    priority <- 0
  }

  new_cluster <- is_new_cluster()
  if (!new_cluster & is.null(cores)) stop("cores not specified, but is mandatory on old cluster.")
  resource_str <- get_resource_str(resources)
  paste(
    "qsub",
    # qsub arguments
    "-e", stderr_log, "-o", stdout_log, "-P", project, "-N", job_name, "-q", queue,
    "-cwd", resource_str, singularity_str,
    "-p", priority,
    # shim for old cluster slot
    if (!new_cluster) {
      paste("-pe multi_slot", cores)
    },
    # command + arguments
    ...
  )
}
