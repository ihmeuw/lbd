####################################################################################################
## Description:   Get in- and (if required) out-of-sample predictive validity draws and then
##                calculate and plot summary statistics.
####################################################################################################

## Load passed arguments
indicator <- commandArgs()[4]
indicator_group <- commandArgs()[5]
run_date <- commandArgs()[6]
core_repo <- commandArgs()[7]

## Load libraries and project functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = core_repo)

## Load config
config <- load_config(repo            = core_repo,
                      indicator_group = indicator_group,
                      indicator       = indicator,
                      post_est_only   = TRUE,
                      run_date        = run_date)

if (class(year_list) == "character") year_list <- eval(parse(text = year_list))
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text = Regions))
outputdir <- paste("<<<< FILEPATH REDACTED >>>>")

## Get in- and out- of sample draws
run_in_oos <- get_is_oos_draws(ind_gp        = indicator_group,
                               ind           = indicator,
                               rd            = run_date,
                               ind_fm        = 'binomial',
                               model_domain  = 'argument_not_used',
                               age           = 0,
                               nperiod       = length(year_list),
                               yrs           = year_list,
                               get.oos       = as.logical(makeholdouts),
                               year_col      = "year",
                               write.to.file = TRUE)

## Load draws and subset to not include ANC
draws.df <- fread(paste0(outputdir, "/output_draws_data.csv"))
draws.df <- draws.df[type != "ANC",]

## Calculate and save PV summary statistics
dir.create(paste0(outputdir, "/summary_metrics/"), recursive = T, showWarnings = F)
pvtab <- rbindlist(lapply(list(c("oos"), c("oos", "year"), c("oos", "region"), c("oos", "point")), function(results_by) {
  pv <- get_pv_table(d               = draws.df,
                     indicator_group = indicator_group,
                     rd              = run_date,
                     indicator       = indicator,
                     aggregate_on    = c("country", "ad1", "ad2"),
                     result_agg_over = results_by,
                     coverage_probs  = c(25, 50, 75, 90, 95),
                     plot_ci         = TRUE,
                     draws           = as.numeric(samples),
                     save_csv        = FALSE,
                     out.dir         = paste0(outputdir, "/summary_metrics/"))
  ldply(pv, .id = "Level")
}), fill = T)
write.csv(pvtab, file = paste0(outputdir, "/summary_metrics/pv_metrics.csv"), row.names = F)
