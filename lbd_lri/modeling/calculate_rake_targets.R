library('data.table')
library('magrittr')

core_repo = '<<<< FILEPATH REDACTED >>>>'
#load get outputs
source(paste0(core_repo,'mbg_central/mbg_functions.R'))                  
source(paste0(user_repo,'mbg_central/prep_functions.R'))               
source(paste0(core_repo,'mbg_central/covariate_functions.R'))             
source(paste0(core_repo,'mbg_central/misc_functions.R'))                 
source(paste0(core_repo,'mbg_central/post_estimation_functions.R'))
source(paste0(core_repo,'mbg_central/gbd_functions.R'))
source(paste0(core_repo,'mbg_central/holdout_functions.R'))
source(paste0(core_repo,'mbg_central/validation_functions.R'))
source(paste0(core_repo,'mbg_central/validation_report_functions.R'))
source(paste0(core_repo,'mbg_central/seegMBG_transform_functions.R'))
source('<<<< FILEPATH REDACTED >>>>')
source('<<<< FILEPATH REDACTED >>>>')

#lri specific settings
severe = .151
cause_id = 322
rei_ids = 187:190
rei_names = c('flu','pneumo','hib','rsv')
year_list = 2000:2016
#location metadata
locs = get_location_metadata(location_set_id = 9, gbd_round_id = 4)
locs = locs[level>=3, location_id]

#a small function to fetch raking targets
fetch_rts = function(measure_id, g2l){
  #takes year list and gbd id from the global environment
  dat = get_outputs(topic = 'cause', measure_id = measure_id, location_id = locs, year_id = eval(parse(text = year_list)), 
                    sex_id = 3, age_group_id = 1, cause_id = cause_id,
                    gbd_round_id = 4, metric = 3)
  
  #convert to severe incidence; stolen from the global environment
  if(measure_id == 6){
    dat[,val:= val * severe]
  }
  
  dat = merge(dat, g2l, by.x= 'location_id', by.y = 'loc_id', all.x = T)
  setnames(dat, c('GAUL_CODE', 'year_id', 'val'), c('name','year', 'mean'))
  
  dat = dat[ ,.(name, year, mean, measure_id, metric_id, location_id)]
  
  return(dat)
  
}

#get the baseline rates
rake_targets = lapply(c(1:6), function(x) fetch_rts(x,g2l))
names(rake_targets) = c('mortality','daly','yld','yll', 'prevalence', 'incidence')

#get the pafs for etiologies
dat = get_outputs(topic = 'rei', measure_id = c(2,3,4), location_id = locs,
                  year_id = year_list, sex_id = 3, age_group_id = 1, 
                  rei_id = rei_ids, cause_id = cause_id,
                  gbd_round_id = 4, metric_id = 2)

#linearly interpolate things
dat[, new_val:= approx(x = year_id, y = val, xout = min(year_id):max(year_id))$y, by = c('rei_id', 'location_id', 'measure_id')]

#keep only required columns
dat = dat[, .(year_id, location_id, measure_id, rei_id, new_val)]

#reshape wide
dat = dcast(dat, formula = year_id + location_id ~ rei_id + measure_id, value.var = 'new_val')

#convert rake targets to etiology rake targets
convert_to_etiology = function(rake_target, fracs, eti_id){
  
  aaa = copy(rake_target)
  aaa = merge(rake_target, fracs, by.x = c('location_id', 'year'), by.y = c('location_id', 'year_id'))
  
  #frac_metric
  meas = unique(rake_target$measure_id)
  if(any(meas %in% c(1,4))){
    frac_metric = 4
  }else if( meas == 2){
    frac_metric =2
  }else{
    frac_metric = 3
  }
  
  print(frac_metric)
  
  aaa[, mean:= mean * get(paste0(eti_id,'_',frac_metric))]
  return(aaa[, .(name, year, mean, location_id)])
}

#create etiologies
flu = lapply(rake_targets,function(x) convert_to_etiology(x,  fracs = dat, 187))
names(flu) = paste0('flu',names(flu))

pneumo = lapply(rake_targets,function(x) convert_to_etiology(x,  fracs = dat, 188))
names(pneumo) = paste0('pneumo',names(pneumo))

hib = lapply(rake_targets,function(x) convert_to_etiology(x,  fracs = dat, 189))
names(hib) = paste0('hib',names(hib))

rsv = lapply(rake_targets,function(x) convert_to_etiology(x,  fracs = dat, 190))
names(rsv) = paste0('rsv',names(rsv))

for(bbb in list(flu, pneumo, hib, rsv)){
  rake_targets = append(rake_targets, bbb)
}

names(rake_targets)

saveRDS(rake_targets, file =  '<<<< FILEPATH REDACTED >>>>')