# Parallel script for post-estimation
#
# Note that prep_postest is used to pass objects to the parallel script
# (saved in a temporarly RData file).  If you want to use logit raking,
# pass "rake_transform" as a string that equals "logit".  You can set this
# up easily in the config file with a row for rake_transform with value "logit"
#
# Sample use:
#
# # Prepare for parallel post-estimation - save file with objs to re-load
# prep_postest(indicator = indicator,
#              indicator_group = indicator_group,
#              run_date = run_date,
#              save_objs = c("core_repo", "gbd", "year_list", "summstats", "rake_transform"))
#
# ## Parallelized post-estimation over region
# postest_script <- "postest_script"
#
# for (s in strata) {
#  qsub <- make_qsub_postest(code = postest_script,
#                            stratum = s,
#                            log_location = 'sharedir',
#                            memory        = mem,
#                            cores         = cores)
#  system(qsub)
# }
#
# ## check to make sure post-est done before continuing
# waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)

## SETUP ################################################################################

stratum = as.character(commandArgs()[4])
run_date = as.character(commandArgs()[5])
indicator = as.character(commandArgs()[6])
indicator_group = as.character(commandArgs()[7])
geos_node = as.logical(commandArgs()[8])
interval_mo <- 12

# Define directories
main_dir <- paste0('<<<< FILEPATH REDACTED >>>>')
temp_dir <- paste0('<<<< FILEPATH REDACTED >>>>')

# Load objects from convenience temp file
load(paste0(temp_dir, "post_est_temp_objs.RData"))

# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(stratum)
message(pop_measure)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

# For now just assign stratum to reg (will need to modify the below for strata beyond reg)
reg <- stratum

# do we really need this again?
sharedir       <- sprintf('<<<< FILEPATH REDACTED >>>>')
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header = FALSE)))

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)


## PREPARE RASTERS, ETC. ################################################################

# Load cell draws
message('Loading Data...')
load(paste0(main_dir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))

# Check if load was successful; stop if not
if (!exists("cell_pred")) {
  message(filename_rds)
  stop("Unable to load cell_pred object! Check to make sure that the relevant object exists.")
}

## SAVE THE RESULTS #####################################################################
message('Saving results...')


# make and save summaries

save_cell_pred_summary <- function(summstat, raked, ...) {
  message(paste0('Making unraked summmary raster for: ',summstat, " (", raked, ")"))
  if (raked == "unraked") cpred <- "cell_pred"
  if (raked == "raked")   cpred <- "raked_cell_pred"
  ras <- make_cell_pred_summary(
    draw_level_cell_pred = get(cpred),
    mask                 = simple_raster,
    return_as_raster     = TRUE,
    summary_stat         = summstat,
    ...)
  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ""),'_',summstat,'_raster'))
}

# Do this as lapply to not fill up memory in global env with big obs
if (is.null(gbd))  rake_list <- c("unraked")
if (!is.null(gbd)) rake_list <- c("unraked", "raked")

summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list)

lapply(1:nrow(summ_list), function(i) {
  summstat <- as.character(summ_list[i, 1])
  raked <- as.character(summ_list[i, 2])
  save_cell_pred_summary(summstat, raked)
})

## Can't pass additional params in the above framework, so will code by hand here
for (r in rake_list) {
  if ("p_below" %in% summstats) {
    save_cell_pred_summary(summstat = "p_below",
                           raked = r,
                           value = 0.8,
                           equal_to = F)
  }
}

# Write a file to mark done
output_dir <- paste0('<<<< FILEPATH REDACTED >>>>')
pathaddin <- paste0('_bin0_',reg,'_0') # To allow us to use waitformodelstofinish()
write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

# All done
message(paste0("Done with post-estimation for ", stratum))
