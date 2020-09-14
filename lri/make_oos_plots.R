library('data.table')
library('ggplot2')

indicator = 'has_lri'
indicator_group = 'lri'
run_dates = c('2018_05_10_16_22_42_1000_draws')
samples = 1000
ho_strat = data.frame(run_date = run_dates)
ho_strat[,'Spatial Holdout'] = c('nid')
period = data.frame(ppp = 0:4, Period = c('2000 - 2016', '2000 - 2003', '2004 - 2007', '2008 - 2011', '2012 - 2017'))
source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/validation_functions.R')

for(run_date  in run_dates){
  diagfolder = '<<<< FILEPATH REDACTED >>>>'
  dir.create(diagfolder)
  draws.df = fread('<<<< FILEPATH REDACTED >>>>')
  
  if(!any(names(draws.df) %in% 'ho_id')){
    draws.df[, ho_id:=nid]
  }
  
  dim(draws.df)
  draws.df = merge(draws.df, period, by = 'ppp')
  dim(draws.df)
  drawcols <- grep("draw", names(draws.df), value = T)
  data <- dplyr::select(draws.df, c(-drawcols))
  
  periodpv = get_pv_table(d = draws.df,
                          indicator_group = indicator_group,
                          rd = run_date,
                          indicator=indicator,
                          aggregate_on='ho_id',
                          result_agg_over = c('oos', 'Period'),
                          draws = as.numeric(samples),
                          out.dir = paste0(diagfolder,'/'),
                          plot_ci = T,
                          plot_ncol = 2)
  
  allpv = get_pv_table(d = draws.df,
                       indicator_group = indicator_group,
                       rd = run_date,
                       indicator=indicator,
                       aggregate_on='ho_id',
                       result_agg_over = c('oos'),
                       draws = as.numeric(samples),
                       out.dir = paste0(diagfolder,'/'),
                       plot_ci = T,
                       plot_ncol = 2)
  
  allpv$Period = '2000 - 2016'
  
  res = rbind(periodpv, allpv, fill = T)
  res$run_date = run_date
  
  res[OOS == T, 'Sample' := 'Out of Sample']
  res[OOS == F, 'Sample' := 'In Sample']
  
  res = merge(res, ho_strat, by = 'run_date')
  
  tab = res[, c('Sample', 'Spatial Holdout', 'Period', 'Mean Err.', 'RMSE', 'Corr.', '95% Cov.'), with = F]
  tab[, '95% Cov.' := round(get('95% Cov.') * 100, 2)]
  tab[, 'Mean Err.' := round(get('Mean Err.'), 5)]
  tab[, 'RMSE' := round(get('RMSE'), 5)]
  tab[, 'Corr.' := round(get('Corr.'), 2)]
  
  write.csv(tab, paste0(diagfolder,'/formated_validation_table.csv'), row.names = F)
}

