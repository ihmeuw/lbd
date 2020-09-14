submit_aggregation_script_experiment <- function(indicator,
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
                                      slots            = 8,
                                      singularity      = 'default',
                                      singularity_opts = NULL,
                                      modeling_shapefile_version = 'current',
                                      raking_shapefile_version = 'current',
                                      measures,
                                      fun_tol_try) {
  
  # Takes vectors of regions, holdouts, and ages and submits qsubs for each of these
  
  dir.create(log_dir)
  dir.create(paste0(log_dir, '/errors'))
  dir.create(paste0(log_dir, '/output'))
  
  # Submit to geo or prod nodes with different default projects. 
  if(geo_nodes) {
    # The geo nodes have default project 'proj_geo_nodes' if the 'proj' argument 
    # is left as NULL and also require the '-l geos_node=TRUE' complex for UGE
    if(is.null(proj)) proj <- "proj_geo_nodes"
    node.flag <- "-l geos_node=TRUE"
  } else {
    # The prod nodes have default project 'proj_geospatial' if the 'proj'
    # argument is left as NULL
    if(is.null(proj)) proj <- "proj_geospatial"
    if(use_c2_nodes) node.flag <- "-q all.q@@c2-nodes" else node.flag <- ""
  }
  
  # At least give a warning if both 'geo_nodes' and 'use_c2_nodes' are requested
  if(geo_nodes & use_c2_nodes) {
    message("WARNING: Both 'geo_nodes' and 'use_c2_nodes' arguments were set to TRUE")
    message(paste0("         Submitting job to LBD nodes under project: '", proj, "'"))
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
  } else {
    # if not, fire up the appropriate version of R depending on the cluster node
    shell <- paste0(corerepo, '/mbg_central/share_scripts/shell_cluster.sh')
  }
  proj.flag <- paste('-P', proj, node.flag, sep = ' ')
  
  for (fun_tol in fun_tol_try){
    
    qsubs_to_make <- expand.grid(regions, holdouts, ages, raked)
    
    for (i in 1:nrow(qsubs_to_make)) {
      region <- qsubs_to_make[i, 1]
      holdout <- qsubs_to_make[i, 2]
      age <- qsubs_to_make[i, 3]
      rake <- qsubs_to_make[i, 4]
      shapefile_version <- ifelse(rake, raking_shapefile_version, modeling_shapefile_version)
      
      for (measure in measures){
        qsub <- paste0('qsub -e ', log_dir, '/errors -o ', log_dir, '/output',
                       ' -cwd -pe multi_slot ', slots, ' ', proj.flag)
        
        # If a Singularity image is being used, pass the name of the image from
        # `get_singularity` as well as any other environmentals the user asked for
        # from the 'singularity_opts' argument to the qsub job
        if(!is.null(singularity)) qsub <- qsub_sing_envs(qsub, singularity_opts,
                                                         sing_image)
        if (rake){r_tag <- 'raked'} else{r_tag <- 'unraked'}
        
        # append the rest of the qsub command
        qsub <- paste0(qsub, ' -N ', indicator, '_', region, '_', measure,'_',r_tag, '_aggregate ',
                       shell, ' ', '<<<< FILEPATH REDACTED >>>>',
                       indicator, ' ',         # arg 4
                       indicator_group, ' ',   # arg 5
                       run_date, ' ',          # arg 6
                       rake, ' ',              # arg 7
                       pop_measure, ' ',       # arg 8
                       overwrite, ' ',         # arg 9
                       age, ' ',               # arg 10
                       holdout, ' ',           # arg 11
                       region, ' ',            # arg 12
                       corerepo, ' ',          # arg 13
                       shapefile_version, ' ', # arg 14
                       measure, ' ',           # arg 15
                       fun_tol)                # arg 16
        
        system(qsub)
      }
    }
  }
}

