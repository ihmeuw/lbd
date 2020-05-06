## Script to run raking and aggregating for diarrhea after cell preds have been generated
## Instructions: If not running interactively, set run_interactive to FALSE and submit a qsub using 01_launcher_dia.R

# Setup ------------------------------------------------------------------------------------------------------------------------------------------

run_interactive <- FALSE

if (run_interactive) {
  
  # set arguments
  indicator_group          <- 'ort'
  indicator                <- 'had_diarrhea'
  reg                      <- 'dia_sssa'
  run_date                 <- '2019_09_17_14_12_53'
  makeholdout              <- FALSE
  get_dalys                <- TRUE

} else {
  
  # Load arguments from qsub
  indicator_group          <- commandArgs()[6]
  indicator                <- commandArgs()[7]
  reg                      <- commandArgs()[10]
  run_date                 <- commandArgs()[16]
  makeholdout              <- commandArgs()[17]
  get_dalys                <- commandArgs()[19]
}

# set holdouts
if (as.logical(makeholdout)) {
  holdouts <- c(0,1,2,3,4,5) 
} else {
  holdouts <- 0
}


# loop over holdouts
for (ho in holdouts) {

  # Set filepaths
  age <- 0
  test <- 0
  pathaddin <- paste0('_bin', age, '_', reg, '_', ho)
  outputdir <- file.path('<<<< FILEPATH REDACTED >>>>')
  
  # Load files
  load(paste0('/<<<< FILEPATH REDACTED >>>>/pre_run_tempimage_', run_date, pathaddin,'.RData'))

  # Overwrite project for this re-launch 
  proj_arg <- commandArgs()[14]
  use_geos_nodes <- commandArgs()[15]
  indicator_group          <- commandArgs()[6]
  indicator                <- commandArgs()[7]
  reg                      <- commandArgs()[10]
  run_date                 <- commandArgs()[16]
  makeholdout              <- commandArgs()[17]
  
  
  ## Aggregate ----------------------------------------------------------------------------------------------------------
  
  # set measures
  measures <- c('deaths', 'incidence', 'prevalence')
  
  for (measure in measures) {
    
    # set specific arguments
    jname           <- paste('agg', reg, indicator, sep = '_')
    
    # set memory by region
    individual_countries <- ifelse(nchar(reg) == 3, TRUE, FALSE)
    mymem <- 120
    if(r == 'dia_malay' | r == 'dia_name') mymem <- 150
    if(r == 'dia_wssa' | r =='dia_south_asia') mymem <- 180
    if(r == 'dia_s_america') mymem <- 270
    
    # set up qsub
    sys.sub <- paste0('qsub -e ', outputdir, '/errors -o ', outputdir, '/output ', 
                      '-l m_mem_free=', mymem, 'G -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                      '-l fthread=2 -l h_rt=03:00:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE ')
    r_shell <- paste0(repo, 'mbg_central/share_scripts/shell_sing.sh')
    script <- '<<<< FILEPATH REDACTED >>>>/mbg_central/share_scripts/frax_script_ort.R'
    args <- paste(user, repo, indicator_group, indicator, config_par, cov_par, reg, proj_arg, 
                  use_geos_nodes, run_date, measure, ho)
    
    # print
    print(paste(sys.sub, r_shell, script, args))
    
    # submit qsub
    system(paste(sys.sub, r_shell, script, args))
    
  }

}