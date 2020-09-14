####################################################################################################
## Description: Submit jobs to run models, prediction, and post-estimation code based on the
##              specifications given in the 'settings' file in the specified directory.
##
## Passed args: main_dir [character] -- home directory for settings and final output.
##              type [character] -- what type of run is this? options:
##                'models' -- just fit models and generate initial area-sex-year-level rates.
##                'post_estimation' -- just run post-estimation aggregation, collapse, and
##                  compile code.
##                'all' -- run all steps (equivalent to running w/ 'models_only' and then
##                  'post_estimation').
##              resub [logical] -- should jobs where the relevant output file already exists be
##                skipped?
##              proj [character] -- cluster project to use when submitting jobs.
##
## Requires:    a fully specified settings file ([main_dir]/settings.csv).
##
## Outputs:     submitted jobs for all steps specified by the passed 'type' argument.
##              'runtime_info.txt' which includes information about provided settings, system
##                settings, and job IDs for submitted jobs.
##
## NOTES:       Must be run from within 'sae_central' directory with 2 slots!!
####################################################################################################

library(TMB)
library(data.table)

rm(list=ls())

## Get and check settings --------------------------------------------------------------------------
main_dir <- commandArgs()[4]
type <- commandArgs()[5]
resub <- as.logical(commandArgs()[6])
proj <- commandArgs()[7]

# set defaults for arguments not provided
if (is.na(type)) type <- "all"
if (is.na(resub)) resub <- F
if (is.na(proj)) proj <- NULL

# check settings (skip for resubmits since this will have already happened)
source("settings.r")
if (resub) {
  get_settings(main_dir)
} else {
  check_settings(main_dir)
}

if (!type %in% c("all", "models", "post_estimation")) stop("type must be 'all', 'models', or 'post_estimation'")

# create temp_dir
dir.create(temp_dir, recursive=T)
if(!file.exists(temp_dir)) stop("temp_dir does not exist and was not created")

# load qsub()
source("qsub.r")

# precompile TMB models
TMB::compile(paste0("models/mod_", model, ".cpp"))

# create a table for holding job ids
jids <- CJ(sex = c(sexes, 3), year = years, level = c(area_var, names(geoagg_files)))

if (!is.null(geoagg_files)) { # we may not be able to aggregate in all years, so year-level combinations with no crosswalk need to be removed.
  for (this_level in names(geoagg_files)) {
    load(geoagg_files[this_level])
    jids <- jids[level != this_level | year %in% unique(weights$year),]
    rm(weights)
  }
}

## Submit model jobs -------------------------------------------------------------------------------
if (type %in% c("models", "all")) {

  # Prep input data
  jids[, prep_inputs := qsub(code = "models/prep_inputs.r",
                             arguments = c(main_dir),
                             slots = 2,
                             sgeoutput = temp_dir,
                             proj = proj,
                             skip_if_exists = if(resub) paste0(temp_dir, "/data.rdata"))]

  # Fit models (by sex)
  jids[sex != 3,
       fit_mod := qsub(code = paste0("models/fit_mod_", model, ".r"),
                       arguments = c(main_dir, sex),
                       hold = unique(prep_inputs),
                       slots = 16,
                       sgeoutput = temp_dir,
                       proj = proj,
                       intel = T,
                       skip_if_exists = if(resub) paste0(temp_dir, "/model_fit_", sex, ".rdata")),
       by='sex']

  # Plot model fits
  jids[sex != 3,
       plot_mod := qsub(code = "models/plot_mod.r",
                        arguments = main_dir,
                        hold = unique(fit_mod),
                        slots = 4,
                        sgeoutput = temp_dir,
                        proj = proj,
                        skip_if_exists = if(resub) paste0(main_dir, "/model_fit.pdf"))]

  # Predict area-sex-level rates (by sex)
  jids[sex != 3,
       pred := qsub(code = "models/pred.r",
                    arguments = c(main_dir, sex),
                    hold = unique(fit_mod),
                    sgeoutput = temp_dir,
                    slots = 54,
                    proj = proj,
                    skip_if_exists = if(resub) paste0(temp_dir, "/est_", area_var, "_", unique(year), "_", sex, ".rdata")),
       by='sex']

} else {
  jids[, pred := 1] # 1 here is just a placeholder jobid

}

## Submit post-estimation jobs (aggregation, collapse, and compile) --------------------------------
if (type %in% c("post_estimation", "all")) {

# Aggregate geographies (by year, sex, level)
  if (!is.null(geoagg_files)) {
    jids[sex != 3 & level != area_var,
         agg_geos := qsub(code = "post_estimation/agg_geos.r",
                          arguments = c(main_dir, year, sex, level),
                          hold = unique(pred),
                          sgeoutput = temp_dir,
                          slots = 9,
                          proj=proj,
                          skip_if_exists = if(resub) paste0(temp_dir, "/est_", level, "_", year, "_", sex, ".rdata")),
         by='year,sex,level']
  }

# Aggregate sexes (by year)
  jids[, agg_sex := qsub(code = "post_estimation/agg_sex.r",
                         arguments = c(main_dir, year, level),
                         hold = na.omit(if (level != area_var) unique(agg_geos) else unique(pred)),
                         sgeoutput = temp_dir,
                         slots = if (level == area_var) 13 else 3,
                         proj=proj,
                         skip_if_exists = if(resub) paste0(temp_dir, "/est_", level, "_", year, "_3.rdata")),
       by='year,level']

# Compile estimates
  jids[, compile := qsub(code = "post_estimation/compile_estimates.r",
                         arguments = c(main_dir),
                         hold = unique(agg_sex),
                         sgeoutput = temp_dir,
                         slots = 2,
                         proj=proj)]
  
# Graph results
  jids[, final_graphs := qsub(code = "post_estimation/model_viz.r",
                         arguments = c(main_dir),
                         hold = unique(compile),
                         sgeoutput = temp_dir,
                         slots = 4,
                         proj=proj)]

}

## Save runtime info to main_dir -------------------------------------------------------------------
options(width=500)
sink(file=paste0(temp_dir, "/runtime_info.txt"), append=T)

cat("\n*********************************************************************************************\n")
cat(paste("\n******************\nTime Submitted\n******************\n", Sys.time(), "\n"))

cat("\n******************\nSettings\n******************\n")
cat(paste("main_dir:", main_dir, "\n"))
cat(paste("type:", type, "\n"))
cat(paste("resub:", resub, "\n\n"))
print(read.csv(paste0(main_dir, "/settings.csv"), header=F, col.names=c("", "Value"), row.names=1), right=F)

cat("\n******************\nJob Ids\n******************\n")
print(jids, nrows=nrow(jids), row.names=F)

cat("\n******************\nGit Status & Version Info\n******************\n")
cat("git status\n")
cat(paste(system("git status", intern=T), collapse="\n"))
cat("\n\ngit log -n 1\n")
cat(system("git log -n 1", intern=T)[1])

cat("\n\n******************\nSystem Info\n******************\n")
Sys.info()

cat("\n******************\nSession Info\n******************\n")
sessionInfo()
sink()
