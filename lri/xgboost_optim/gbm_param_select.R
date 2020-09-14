#collect arguments
indicator = as.character(commandArgs()[4])
indicator_group = as.character(commandArgs()[5])
run_date = as.character(commandArgs()[6])
region = as.character(commandArgs()[7])

for(aaa in c(indicator, indicator_group, run_date, region)){
  message(aaa)
}

pathaddin <- paste0('_bin',0,'_',region,'_',0)

library('mbgstacking')

#load from the first checkpoint: has input data and other things
load('<<< FILEPATH REDACTED >>>')

wf = '<<< FILEPATH REDACTED >>>'
dir.create(wf, recursive = T)
logloc = '<<< FILEPATH REDACTED >>>'
oot = paste0(logloc, 'output/')
eet = paste0(logloc, 'errors/')
r_path = '<<< FILEPATH REDACTED >>>'

sgeset =init_sge(working_folder = wf, r_path = r_path, repeat_iterations = 2,
                 project_name = 'proj_geospatial_ppp', other_options = NULL, slots_per_job = c(3,6),
                 package_location = '<<< FILEPATH REDACTED >>>')


#create a bunch of brt models varying the parameters
pg = expand.grid(eta = c(.005,.007,.01,.02,.03, .04, .05,.06),
                 max_depth = c(8, 10, 12, 14, 16),
                 min_child_weight = 1,
                 subsample = c(.5, .7, 1),
                 colsample_bytree = c(.5, .7, 1),
                 nrounds = c(500, 1000, 2000, 5000, 10000))

#make the various stacker objects
brtobjs = lapply(1:nrow(pg), function(x){
  
  init_xgboost(model_name = paste0('xgb_',x),
               params_arg = list(eta = pg$eta[x],
                                 max_depth = pg$max_depth[x],
                                 min_child_weight = pg$min_child_weight[x],
                                 subsample = pg$subsample[x],
                                 colsample_bytree = pg$colsample_bytree[x],
                                 nthread = 1),
               nrounds = pg$nrounds[x],
               binomial_evaluation = 'prev',
               weight_column = 'weighted_n')
  
})

#do hold out folds by dropping out datasets rather than random rows
#this is the input data for reg

nids = unique(df[,nid])


for(sfo in 1:3){
  #randomly rearrange
  nids = sample(nids, size = length(nids))
  
  #chop into 5 folds
  foldz = split(nids, ceiling(seq_along(nids)/(length(nids)/3)))
  
  iter = 1
  for(fdz in foldz){
    df[nid %in% fdz, paste0('sfold_', sfo):=iter]
    iter = iter + 1
  }
}

steak = init_stacker(brtobjs, #earth_model,  #em_def, em_chg, enet, gm_def,
                     inlist = T,
                     data = df,
                     indicator = indicator,
                     indicator_family = indicator_family,
                     fe_equation = paste(names(all_cov_layers), sep = ' + '),
                     covariate_layers = all_cov_layers,
                     centre_scale = T,
                     time_var = 'year',
                     time_scale = year_list,
                     weight_col = 'weight',
                     num_fold_cols = c('sfold_1','sfold_2', 'sfold_3'),
                     num_folds = 3,
                     cores = cores_to_use,
                     sge_parameters = sgeset)

#run child models
model_results = run_stacking_child_models(steak)

#make rasters
saveRDS(pg, '<<< FILEPATH REDACTED >>>')
saveRDS(model_results, '<<< FILEPATH REDACTED >>>')