####################################################################################################
## Description:   Get in- and (if required) out-of-sample predictive validity draws and then
##                calculate and plot summary statistics.
####################################################################################################

## Load passed arguments
indicator         <- commandArgs()[4]
indicator_group   <- commandArgs()[5]
run_date          <- commandArgs()[6]
core_repo         <- commandArgs()[7]
shapefile_version <- commandArgs()[8]

## Load libraries and project functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = core_repo)

indic_repo <- paste0("/lbd_hiv/")

outputdir <- paste('<<<< FILEPATH REDACTED >>>>')

config <- set_up_config(repo            = outputdir,
                        core_repo       = core_repo,
                        indicator       = "",
                        indicator_group = "",
                        config_file     = paste0(outputdir,'config.csv'),
                        covs_file       = paste0(indic_repo, 'mbg/covariates/2_modeling/cov_list.csv'),
                        post_est_only    = TRUE,
                        run_tests = FALSE)


if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(z_list) == "character") z_list <- eval(parse(text=z_list))
if(exists('sex_list')) if (class(sex_list) == "character") sex_list <- eval(parse(text=sex_list))

## Get in- and out- of sample draws
if(indicator!='hiv_prev_disagg') {
  run_in_oos <- get_is_oos_draws(ind_gp        = indicator_group,
                                 ind           = indicator,
                                 rd            = run_date,
                                 ind_fm        = 'binomial',
                                 model_domain  = 'argument_not_used',
                                 age           = 0,
                                 sex           = 0,
                                 nperiod       = length(year_list),
                                 yrs           = year_list,
                                 get.oos       = as.logical(makeholdouts),
                                 year_col      = "year",
                                 write.to.file = TRUE,
                                 shapefile_version = modeling_shapefile_version)
  } else {
  source(paste0(indic_repo, 'mbg/hiv_prev_disagg/3_functions/get_is_oos_draws_function.r'))
  run_in_oos <- get_is_oos_draws_agebins(ind_gp        = indicator_group,
                               ind           = indicator,
                               rd            = run_date,
                               ind_fm        = 'binomial',
                               model_domain  = 'argument_not_used',
                               agebins       = z_list,
                               sexes         = sex_list,
                               yrs           = year_list,
                               get.oos       = as.logical(makeholdouts),
                               year_col      = "year",
                               write.to.file = TRUE,
                               shapefile_version = modeling_shapefile_version)
  }

## Load draws and subset to not include ANC
draws.df <- fread(paste0(outputdir, "/output_draws_data.csv"))
if (str_detect(indicator, "hiv")) draws.df <- draws.df[type != "ANC",]

## Calculate and save PV summary statistics
dir.create(paste0(outputdir, "/summary_metrics/"), recursive = T, showWarnings = F)

if(!length(z_list)>1) {
  pvtab <- rbindlist(lapply(list(c("oos"), c("oos", "year"), c("oos", "region"), c("oos", "point")), function(results_by) {
  pv <- get_pv_table(d               = draws.df,
                     indicator_group = indicator_group,
                     rd              = run_date,
                     indicator       = indicator,
                     aggregate_on    = c("country", "ad1", "ad2", "cluster_id"),
                     result_agg_over = results_by,
                     coverage_probs  = c(25, 50, 75, 90, 95),
                     plot_ci         = TRUE,
                     draws           = as.numeric(samples),
                     save_csv        = FALSE,
                     out.dir         = paste0(outputdir, "/summary_metrics/"))
  ldply(pv, .id = "Level")
}), fill = T)
} else {
  pvtab <- rbindlist(lapply(list(c("oos"), c("oos", "year"), c("oos", "region"), c("oos", "point"), c("oos", "agebin")), function(results_by) {
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
}
write.csv(pvtab, file = paste0(outputdir, "/summary_metrics/pv_metrics.csv"), row.names = F)
