root <- '<<<< FILEPATH REDACTED >>>>' #filepath root
indicator <- "has_lri"
run_date <- Sys.Date()
folder_out <- '<<<< FILEPATH REDACTED >>>>' #file out path

latest_collapse <- '<<<< FILEPATH REDACTED >>>>' #latest collapse data ready to be resampled

polydat <- fread(paste0(folder_out, latest_collapse))

#subset to only polygon data
polydat <- subset(polydat, point == 0)

#create list of job names for hold_jid
agg_jobs = list()

#launch qsubs by shapefile
for (shp in unique(polydat$shapefile)) {
  jname <- paste(indicator, shp, sep = "_") #name of job
  agg_jobs <- append(agg_jobs, jname) #add name of job to job name list
  mythreads <- "" #threads to be requested
  mymem <- '8G' #memory to be requested
  proj <- '<<<< PROJECT NAME REDACTED >>>>' #cluster project name
  user <- Sys.info()[["user"]] # name of user
  error_dir <- '<<<< FILEPATH REDACTED >>>>' #error directory
  output_dir <- '<<<< FILEPATH REDACTED >>>>' #output directory
  run_time <- "" #amount of time per job
  sys.sub <- paste0("qsub ",proj,paste0(" -e ", error_dir," -o ", output_dir," "),
                    "-cwd -N ", jname, " ", "-l fthread=", mythreads, " ",
                    "-l m_mem_free=", mymem, " -q all.q -l h_rt=", run_time, " -l archive=TRUE")
  script <- "resample_child.R"
  r_shell <- '<<<< FILEPATH REDACTED >>>>' #path to R shell script

  #build qsub call
  args <- paste(shp, indicator, run_date, latest_collapse)
  #launch qsub
  system(paste(sys.sub, r_shell, script, args))
}
#stall code until qsubs finish
suppressWarnings(system(paste("qsub -b y -sync y -l m_mem_free=1G -l fthread=1 -q all.q -l h_rt=00:01:00 -P ", proj, " -hold_jid",
                              paste(agg_jobs, collapse =","), paste("-N ", "lri_aggs"), shQuote("echo 1"))))
