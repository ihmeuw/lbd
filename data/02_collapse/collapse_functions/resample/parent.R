source('<<<< FILEPATH REDACTED >>>>collapse_functions/general_functions.R')

proj <- '-P proj_geospatial'
user <- '<<<< USERNAME REDACTED >>>>'

indicator <- 'has_lri'
run_date <- Sys.Date()
folder_out <- '<<<< FILEPATH REDACTED >>>>'

latest_collapse <- get_latest_file(folder_out, '*.csv')

polydat <- fread(paste0(folder_out, latest_collapse))
polydat <- subset(polydat, point == 0)

agg_jobs = list()
  
for (shp in unique(polydat$shapefile)) { 
  jname <- paste(indicator, shp, sep = "_")
  agg_jobs <- append(agg_jobs, jname)
  mythreads <- 1
  mymem <- '4G'
  sys.sub <- paste0("qsub ",proj,paste0(" -e <<<< FILEPATH REDACTED >>>> -o <<<< FILEPATH REDACTED >>>> "),
                    "-cwd -N ", jname, " ", "-l fthread=", mythreads, " ", "-l m_mem_free=", mymem, ' -q all.q -l h_rt=00:25:00 -l archive=TRUE')
  script <- "<<<< FILEPATH REDACTED >>>>collapse_functions/resample/child.R"
  r_shell <- '<<<< FILEPATH REDACTED >>>>'
  
  args <- paste(shp, indicator, run_date, latest_collapse)
  system(paste(sys.sub, r_shell, script, args)) 
}
#stall code until qsubs finish
suppressWarnings(system(paste('qsub -b y -sync y -l m_mem_free=1G -l fthread=1 -q all.q -l h_rt=00:01:00 -P proj_geospatial -hold_jid', 
                              paste(agg_jobs, collapse =','), paste('-N ', 'lri_aggs'), shQuote('echo 1'))))