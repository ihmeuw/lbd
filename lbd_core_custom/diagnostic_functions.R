## make_model_diagnostics() --------------------------------------------------->
#'
#' @title Function to submit a qsub to run the model_diagnostics.R script
#'
#' @description \code{make_model_diagnostics} creates a qsub string to run the '
#'   model_diagnostics.R script on the cluster and does the system call to actually
#'   submit the job
#'
#' @param code Name of script, with relative path if desired.
#' 
#' @param code_path Full path to R script. Overrides \code{code} and \code{script_dir}
#' 
#' @param cores Number of threads. Default: 5.
#' 
#' @param memory RAM to be reserved, in GBs
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
#' @param priority Job priority that can be deprioritized if needed, and can only be used for values in [-1023,0]. Default = 0.
#' This value will get bounded to 0 or -1023 if the user supplies a value outside those bounds.
#'
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
#'   for LBD are kept at the default location of the sing image.
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
#'   [default = list(SET_OMP_THREADS = cores, SET_MKL_THREADS = cores)<Paste>]
#'
#' @return None
#'
#' @seealso This function uses the following functions found in
#'   mbg_central/misc_functions.R:
#'   \code{\link{get_singularity}}
#'   \code{\link{qsub_sing_envs}}
#'
#' @export
#'
make_model_diagnostics <- function(user = Sys.info()["user"],
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
                                   script_dir = "mbg_central/share_scripts",
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
  script_dir <- '<<<< FILEPATH REDACTED >>>>'
  code <- path_join(script_dir, paste0(code, ".R"))
  
  # If code_path is not NULL, then override `code`
  if(!is.null(code_path)) {
    code <- code_path
  }
  
  # Define remaining job attributes
  job_name <- paste0("job_dx_", indic)
  run_time <- get_run_time(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  queue <- get_queue(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  shell <- paste0(corerepo, "/mbg_central/share_scripts/shell_sing.sh")
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

