## Script to run raking and aggregating for diarrhea after cell preds have been generated


# Setup ------------------------------------------------------------------------------------------------------------------------------------------

run_interactive <- F

if (run_interactive) {
  
  # set arguments
  indicator_group          <- 'ort'
  indicator                <- 'any_ors'
  reg                      <- 'SLE'
  run_date                 <- '2019_12_13_12_00_46'
  makeholdout              <- TRUE
  
} else {
  
  # Load arguments from qsub
  indicator_group          <- commandArgs()[6]
  indicator                <- commandArgs()[7]
  reg                      <- commandArgs()[10]
  run_date                 <- commandArgs()[16]
  makeholdout              <- commandArgs()[17]
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
  pathaddin <- paste0('_bin',age,'_',reg,'_',ho)
  outputdir <- file.path('<<<< FILEPATH REDACTED >>>>')
  
  # Load files
  load(paste0('<<<< FILEPATH REDACTED >>>>/pre_run_tempimage_', run_date, pathaddin,'.RData'))
  
  # make sure we have repo argument
  repo <- core_repo
  
  # Overwrite project for this re-launch 
  proj_arg <- commandArgs()[14]
  use_geos_nodes <- commandArgs()[15]
  indicator_group          <- commandArgs()[6]
  indicator                <- commandArgs()[7]
  reg                      <- commandArgs()[10]
  run_date                 <- commandArgs()[16]
  makeholdout              <- commandArgs()[17]
  
  
  ## Aggregate ----------------------------------------------------------------------------------------------------------
  
  # set specific arguments
  measure         <- 'prevalence'
  jname           <- paste('agg', reg, indicator, sep = '_')
  
  # set memory by region
  individual_countries <- ifelse(nchar(reg) == 3, TRUE, FALSE)
  mymem <- 120
  if (as.logical(individual_countries) & reg != 'MLI') mymem <- 50
  
  # set up qsub
  sys.sub <- paste0('qsub -e ', outputdir, '/errors -o ', outputdir, '/output ', 
                    '-l m_mem_free=', mymem, 'G -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                    '-l fthread=2 -l h_rt=00:12:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE ')
  r_shell <- paste0(repo, 'mbg_central/share_scripts/shell_sing.sh')
  script <- '<<<< FILEPATH REDACTED >>>>/ort/mbg_central/share_scripts/frax_script_ort.R'
  args <- paste(user, repo, indicator_group, indicator, config_par, cov_par, reg, proj_arg, 
                use_geos_nodes, run_date, measure, ho)
  
  # submit qsub
  system(paste(sys.sub, r_shell, script, args))

}
