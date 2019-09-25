rerun = F

#collect gbm params
indicator = 'has_lri'
indicator_group = 'lri'
run_date = '<<<< RUN DATE REDACTED >>>>'
region = 'essa'
for(aaa in c(indicator, indicator_group, run_date, region)){
  message(aaa)
}

pathaddin <- paste0('_bin',0,'_',region,'_',0)
package_lib    <- '<<<< FILEPATH REDACTED >>>>'

#load packages
library('ggplot2')
library('data.table')
.libPaths(package_lib)
library('mbgstacking')
library('wCorr')

#load from the first checkpoint
load('<<<< FILEPATH REDACTED >>>>')

weighted.rmse <- function(error, w) sqrt( sum(  (w/sum(w)) * ((error)^2) ) )

pg = readRDS('<<<< FILEPATH REDACTED >>>>')

calc_best_mods = function(rrr){
  wf = '<<<< FILEPATH REDACTED >>>>'

  #read the
  steak = readRDS('<<<< FILEPATH REDACTED >>>>')
  modres = readRDS('<<<< FILEPATH REDACTED >>>>')
  preds = modres[[1]]

  #clear up space
  rm(modres)

  #out of sample validity
  cv_names = grep('cv_pred', names(preds), fixed = T, value = T)
  preds = preds[, c('rid', cv_names),with = F]

  preds = cbind(steak$data, preds)
  preds[, region:= rrr]
  return(preds)
}

if(rerun){
  res = parallel::mclapply(region_list, calc_best_mods, mc.cores = 6)
  res = rbindlist(res)
  saveRDS(res, '<<<< FILEPATH REDACTED >>>>')
}else{
  res = readRDS('<<<< FILEPATH REDACTED >>>>')
}

setDT(pg)

#find the model that was eventually selected as "best"
pg[, id := .I]

#keep the 150 best models based on rmse
rmses = res[, lapply(paste0('xgb_',pg[,id],'_cv_pred'), function(cvn) weighted.rmse(has_lri/N - get(cvn), data_weight * N ))]
rmses = unlist(rmses)

besty = pg[eta==.007 & max_depth ==8 & subsample == .7 & colsample_bytree ==.5 & nrounds ==1000, id]
res[, best_model := get(paste0('xgb_',besty,'_cv_pred'))]

#for each model, calculate the weighted correlation
thecorrs = lapply(paste0('xgb_',pg[,id],'_cv_pred'), function(mod) weightedCorr(res[,best_model], res[,get(mod)], weights = res[,alt_weight], method = 'pearson'))

rmse_rel_best = rmses - rmses[besty]

#weighted mean lri prevalence
mean_prev = weighted.mean(res[,has_lri/N], res$alt_weight)

comparez = data.table(model = paste0('xgb_',pg[,id],'_cv_pred'), rmse_dif = rmse_rel_best, corr_with_best = unlist(thecorrs))

pdf('<<<< FILEPATH REDACTED >>>>', height = 10, width = 10)

#high level overview
g = ggplot(comparez, aes(x = rmse_dif, y = corr_with_best)) + geom_point() +theme_bw() + geom_vline(xintercept = 0, color = 'red') +
  labs(title = paste('Corr and RMSE Dif relative to "Best" BRT'), subtitle = paste('Mean Prev:', round(mean_prev,6)))
plot(g)


#focused
g2 = g + coord_cartesian(xlim = c(-2e-5,.001))
plot(g2)

#compare the models with rmse dif <.0006 & correlation of less than .7. Choose the top 10 lowest correlated
setorder(comparez, corr_with_best)
uncorr = comparez[rmse_dif <.0006,][1:10,model]

for(unc in uncorr){
  res[, compare_model := get(unc)]
  p = ggplot(res, aes(x = best_model, y = compare_model)) + geom_point() + theme_bw() + geom_abline(slope = 1, color = 'red') +ylab(unc) +
    ggtitle(paste0('Best Model vs. ', unc)) + xlim(c(0,max(c(res$best_model,res$compare_model)))) + ylim(c(0,max(c(res$best_model,res$compare_model))))
  plot(p)
}

dev.off()

#examine which models those are
pg[, name:= paste0('xgb_',id,'_cv_pred')]

#take the average of the parameters
params = c('eta','max_depth','min_child_weight','subsample','colsample_bytree','nrounds')
lapply(params, function(x) table(pg[name %in% comparez[corr_with_best<.7 & abs(rmse_dif) <.0006, model], get(x)]))
