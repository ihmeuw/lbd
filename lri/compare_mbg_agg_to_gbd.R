library('magrittr')
library('data.table')
library('rgeos')
library('sp')
library('raster')
library('ggplot2')
library('dplyr')

package_lib  = '<<<< FILEPATH REDACTED >>>>'
.libPaths(package_lib)
library('rgdal', lib = package_lib)

repo = '<<<< FILEPATH REDACTED >>>>'
setwd(repo)

#load all the functions
funks = list.files(paste0(repo,'mbg_central/'), pattern = 'functions',full.names = T)
for(fff in funks){
  source(fff)
}

source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
gaul_to_loc = '<<<< FILEPATH REDACTED >>>>'
g2l = readRDS(gaul_to_loc)
g2l <- unique(g2l[,.(ADM0_CODE, ADM0_NAME)])

source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
loc_ids <- read.csv('<<<< FILEPATH REDACTED >>>>')
country <- subset(loc_ids, level >=3)$location_id
loc_ids <- filter(loc_ids, location_type=='admin0')
loc_ids <- loc_ids[,c('location_name','location_id')]

#fix discrepant location names
loc_ids$location_name <- as.character(loc_ids$location_name)
loc_ids$location_name[loc_ids$location_name == "Cote d'Ivoire"] <- "Côte d'Ivoire"
loc_ids$location_name[loc_ids$location_name == 'Congo'] <- 'Republic of Congo'
loc_ids$location_name[loc_ids$location_name == 'Sao Tome and Principe'] <- "São Tomé and Príncipe"
loc_ids$location_name[loc_ids$location_name == 'The Gambia'] <- "Gambia"


g2l <- merge(g2l, loc_ids, by.x = 'ADM0_NAME', by.y = 'location_name',all.x = TRUE)
g2l$loc_id <- g2l$location_id
g2l <- g2l[,c('ADM0_NAME','ADM0_CODE','loc_id')]

indicator <- 'has_lri'
indicator_group <- 'lri'
run_dates <- c('2017_stage_1_updates')

summarize_admin0 = function(measure){
  admin = paste0(in_dir, indicator,'_raked_admin_draws_eb_bin0_0_',measure,'.RData')
  load(admin)
  
  admin_0 = copy(admin_0)
  
  draw_cols = grep('V[0-9]*', names(admin_0), value = T)
  
  #get rates
  admin_0[, paste0('mean_rate') := rowMeans(.SD), .SDcols = draw_cols]
  admin_0[, paste0('lower_rate') := apply(.SD[,draw_cols,with=F], 1, quantile, probs = 0.025, na.rm = T)]
  admin_0[, paste0('upper_rate') := apply(.SD[,draw_cols,with=F], 1, quantile, probs = 0.975, na.rm = T)]
  
  #remove draws
  admin_0 = admin_0[,!names(admin_0) %in% draw_cols, with = F]
  admin_0[, measure := measure]
  
  return(admin_0)
}
for(run_date in run_dates){
  
  out_dir = '<<<< FILEPATH REDACTED >>>>'
  in_dir = '<<<< FILEPATH REDACTED >>>>'
  dir.create(out_dir, recursive = T)
  
  mbg = rbindlist(lapply(c('prevalence','incidence','mortality'), summarize_admin0), fill = T)
  
  #add location id
  mbg = merge(mbg, g2l, by.x = 'ADM0_CODE', by.y = 'ADM0_CODE', all.x = T)
  
  #convert to measure id
  meas_id = data.table(measure = c('prevalence','incidence','mortality'), measure_id =c(5,6,1))
  mbg = merge(mbg, meas_id, by = 'measure')
  
  #get corrosponding outputs from GBD
  gbd = get_outputs(topic = 'cause', measure_id = c(1,5,6), location_id = country, year_id = 2000:2017, sex_id = 3, age_group_id = 1, metric_id = c(1,3), gbd_round_id = 5, cause_id = 322)
  
  gbd[metric_name == 'Number', metric_name := 'count']
  gbd[metric_name == 'Rate', metric_name := 'rate']
  
  setnames(gbd, c('val', 'upper','lower'), c('gbd_mean','gbd_lower','gbd_upper'))
  gbd = gbd[,.(location_id, measure_id, year_id, metric_name, gbd_mean, gbd_lower, gbd_upper)]
  gbd = dcast(gbd, location_id + measure_id + year_id ~ metric_name, value.var = c('gbd_mean','gbd_lower','gbd_upper'))
  
  #merge gbd and mbg
  mer = merge(mbg, gbd, by.x = c('loc_id', 'year', 'measure_id'), by.y = c('location_id', 'year_id', 'measure_id'))
  
  stopifnot(nrow(mer) == nrow(mbg))
  
  pdf(paste0(out_dir, '/compare_gbd_mbg.pdf'), width = 12, height = 10)
  overall_rate = ggplot(mer, aes(x = mean_rate, y = gbd_mean_rate)) + geom_point() + theme_bw() + facet_wrap(~measure, scales = 'free') + 
    geom_abline(slope = 1) + ggtitle("Rate compare, GBD vs. MBG")
  plot(overall_rate)
  
  dev.off()
  
  #by country, make a scatter of raked mbg prev by gbd prev
  prev = mer[measure == 'prevalence',]
  pdf(paste0(out_dir, '/examine_prev_raking.pdf'), width = 12, height = 10)
  for(loc in unique(prev[,ADM0_NAME])){
    dat = prev[ADM0_NAME == loc,]
    
    maxval = max(c(dat[,mean_rate], dat[,gbd_mean_rate]))
    
    g = ggplot(dat, aes(x = mean_rate, y = gbd_mean_rate)) + geom_point() + theme_bw() + xlim(0,maxval) + ylim(0, maxval) + 
      geom_abline(slope = 1) + ggtitle(paste('Prevalence raking for',loc))
    plot(g)
    
    
  }
  dev.off()
  
  #by country, make a scatter of raked mbg inc by gbd inc
  inc = mer[measure == 'incidence',]
  pdf(paste0(out_dir, '/examine_inc_raking.pdf'), width = 12, height = 10)
  for(loc in unique(inc[,ADM0_NAME])){
    dat = inc[ADM0_NAME == loc,]
    
    maxval = max(c(dat[,mean_rate], dat[,gbd_mean_rate]))
    
    g = ggplot(dat, aes(x = mean_rate, y = gbd_mean_rate)) + geom_point() + theme_bw() + xlim(0,maxval) + ylim(0, maxval) + 
      geom_abline(slope = 1) + ggtitle(paste('Incidence raking for',loc))
    plot(g)
    
  }
  dev.off()
  
  #by country, make a scatter of raked mbg mort by mort 
  mort = mer[measure == 'mortality',]
  pdf(paste0(out_dir, '/examine_mort_raking.pdf'), width = 12, height = 10)
  for(loc in unique(mort[,ADM0_NAME])){
    dat = mort[ADM0_NAME == loc,]
    
    maxval = max(c(dat[,mean_rate], dat[,gbd_mean_rate]))
    
    g = ggplot(dat, aes(x = mean_rate, y = gbd_mean_rate)) + geom_point() + theme_bw() + xlim(0,maxval) + ylim(0, maxval) + 
      geom_abline(slope = 1) + ggtitle(paste('Mortality raking for',loc))
    plot(g)
    }
  dev.off()
}
