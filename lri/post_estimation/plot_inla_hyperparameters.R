#load the two functions
source('<<<< FILEPATH REDACTED >>>>/plot_hyperparameters.R')
source('<<<< FILEPATH REDACTED >>>>/plot_priors_posteriors.R')

indicator <- 'has_lri'
indicator_group <- 'lri'
run_dates <- c('2019_04_25_15_49_17', '2019_04_25_15_50_13', '2019_04_01_15_28_31')
age <- 0
holdout <- 0
region_list <- c('cssa', 'essa', 'wssa', 'name', 'sssa')


#make hyperparameter plots
for (run_date in run_dates){
  plot_hyperparameters(indicator = indicator,
                       indicator_group = indicator_group,
                       run_date = run_date,
                       age = age,
                       holdout = holdout,
                       save_file = NULL,
                       regs = region_list)
  
  for (region in region_list){
    plot_spatial_priors(indicator = indicator,
                        indicator_group,
                        run_date,
                        region,
                        holdout = 0,
                        age = 0)
  }
}
