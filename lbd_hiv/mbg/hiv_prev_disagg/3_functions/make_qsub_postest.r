make_qsub_postest <- function(user = Sys.info()["user"],
                              code,
                              code_path = NULL,
                              cores = slots,
                              memory = 100,
                              proj = NULL,
                              ig = indicator_group,
                              indic = indicator,
                              stratum = "test",
                              agebin  = 0,
                              sex     = 0,
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
  job_name <- paste("job_pe", indic, agebin, sex, stratum, sep = "_")
  job_name <- paste0(job_name, addl_job_name)
  run_time <- get_run_time(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  queue <- get_queue(use_geo_nodes = geo_nodes, use_c2_nodes = use_c2_nodes, queue = queue, run_time = run_time)
  shell <- paste0(corerepo, "/mbg_central/share_scripts/shell_sing.sh")
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
    shell, code, stratum, agebin, sex, rd, indic, ig, as.character(geo_nodes),
    modeling_shapefile_version, raking_shapefile_version, subnat_raking
  )
  
  return(qsub)
}
