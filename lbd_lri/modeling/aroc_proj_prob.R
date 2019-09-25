library('data.table')

#indicator, indicator_group, run_date, mbg functions, datadir, diagfolder, region_list, measure_list, and year_list inherited from model_post_processing.R

#make folder structure
for(fff in c('proj','aroc','target_probs')){
  dir.create(paste0(datadir, 'pred_derivatives/',fff), recursive = T)
}


source('<<<< FILEPATH REDACTED >>>>') #aroc proj functions
source('<<<< FILEPATH REDACTED >>>>') #admin functions

#create LRI goals
if(indicator == 'has_lri'){
  g1 = add_goal(target_year = 2025,
                target = 0.75, #75% reduction since 2010
                baseline_year = 2010,
                target_type = "less",
                abs_rel = "relative",
                pred_type = 'admin')
  g1$measure ='incidence'

  g2 = add_goal(target_year = 2025,
                target = 0.75, #75% reduction since 2010
                baseline_year = 2010,
                target_type = "less",
                abs_rel = "relative",
                pred_type = 'cell')
  g2$measure ='incidence'

  g3 <- add_goal(target_year = 2025,
                 target = 3/1000, #1 per 1000
                 baseline_year = 2010,
                 target_type = "less",
                 abs_rel = "absolute",
                 pred_type = 'admin')
  g3$measure ='mortality'

  g4 <- add_goal(target_year = 2025,
                 target = 3/1000, #1 per 1000
                 baseline_year = 2010,
                 target_type = "less",
                 abs_rel = "absolute",
                 pred_type = 'cell')
  g4$measure ='mortality'

  admin_goals = list(g1, g3)
  cell_goals = list(g2, g4)
  matpredname = 'has_lri_%s_%s_raked_0.rds'
}

if(arocproj_results){
  #calculate AROC and project at the admin level
  for(mmm in c('incidence','mortality')){
    admin_aroc(admin_levels = 0:2,
               indicator_group = indicator_group,
               indicator = indicator,
               run_date = run_date,
               measure = mmm,
               year_list = year_list,

               uselogit = F,
               save = T)

    admin_proj(admin_levels = 0:2, indicator_group = indicator_group, indicator = indicator, run_date = run_date, proj_years = c(2025),
               measure = mmm, year_list = year_list, uselogit = F, save = T)

  }

  #for admin goals, compare to target
  for(ggg in admin_goals){
    compare_to_target_admin(admin_levels = 0:2, indicator_group = indicator_group, indicator = indicator, run_date = run_date,
                            goal_obj = ggg, measure = ggg$measure, year_list = year_list, uselogit = F, save = T)
  }

  #aroc and project at the pixel level
  for(mmm in c('incidence','mortality')){
    make_aroc(ind_gp = indicator_group,
              ind = indicator,
              rd = run_date,
              matrix_pred_name = matpredname,
              type = 'cell',
              raked = T,
              measure = mmm,
              year_list = year_list,
              uselogit = F) #set to false because incidence can be over 1

    make_proj(ind_gp = indicator_group,
              ind = indicator,
              rd = run_date,
              type = 'cell',
              proj_years = 2025,
              measure = mmm,
              skip_cols = NULL,
              matrix_pred_name = matpredname,
              year_list = year_list,
              uselogit = F)

  }


  #for cell goals, compare to target
  for(ggg in cell_goals){
    compare_to_target(ind_gp = indicator_group, ind = indicator, rd = run_date, goal_obj = ggg,
                      measure = ggg$measure, year_list = year_list, uselogit = F,
                      matrix_pred_name = matpredname)
  }
}

if(arocproj_graphics){
  #maps and tables
  #aroc for all the indicators
  for(measure in c('incidence','mortality')){
    aroc_ras = lapply(region_list, function(x) make_new_ras(x, '<<<< FILEPATH REDACTED >>>>'))
    names(aroc_ras) = NULL
    aroc_ras = do.call(raster::merge, aroc_ras)
    outj = '<<<< FILEPATH REDACTED >>>>'
    writeRaster(x = aroc_ras,filename = '<<<< FILEPATH REDACTED >>>>', overwrite = T)
  }

  #mortality goal in 2016
  gappd15 = lapply(region_list, function(x) make_new_ras(x, '<<<< FILEPATH REDACTED >>>>', T, .003))
  names(gappd15) = NULL
  gappd15 = do.call(raster::merge, gappd15)
  gappd15 = gappd15[[16]]
  outj = '<<<< FILEPATH REDACTED >>>>'
  writeRaster(gappd15, filename = '<<<< FILEPATH REDACTED >>>>', overwrite = T)
  #

  #incidence and mortality projections in 2025
  for(ggg in cell_goals){
    #get the file names

    if(ggg$abs_rel=='absolute'){
      prob_pre = '<<<< FILEPATH REDACTED >>>>'
    }else{
      prob_pre= '<<<< FILEPATH REDACTED >>>>'
    }

    projras = lapply(region_list, function(x) make_new_ras(x, '<<<< FILEPATH REDACTED >>>>'))
    names(projras) = NULL
    projras = do.call(raster::merge, projras)
    outj = '<<<< FILEPATH REDACTED >>>>'

    writeRaster(x = projras,filename = '<<<< FILEPATH REDACTED >>>>', overwrite = T)
  }

  #get the admin unit info
  #aroc
  for(measure in c('incidence', 'mortality')){

    aroc_ad = readRDS('<<<< FILEPATH REDACTED >>>>')
    aroc_ad[,value := rowMeans(.SD), .SDcols = grep('V[0-9]*', names(aroc_ad), value = T)]
    aroc_ad[,year := max(year_list)]
    aroc_ad = aroc_ad[,.(ADM1_CODE, year, value)]

    outj = '<<<< FILEPATH REDACTED >>>>'

    write.csv(aroc_ad, '<<<< FILEPATH REDACTED >>>>', row.names = F)

  }

  #2025 projections
  for(ggg in admin_goals){

    prob_pre= '<<<< FILEPATH REDACTED >>>>'


    prob_ad = readRDS('<<<< FILEPATH REDACTED >>>>')
    setnames(prob_ad, 'goal_prob', 'value')

    outj = '<<<< FILEPATH REDACTED >>>>'

    #load draws if mortality
    if(ggg$measure == 'mortality'){
      dat = summarize_admin(measure = ggg$measure,admin_level = 1,draws = T, raked = T)

      #check if the rates are below the target
      dat = dat[year == max(year_list),]
      dcols = grep('V[0-9]*\\_rate', names(dat), value = T)
      dat[, paste0('prob',seq(dcols)) := lapply(dcols, function(x) get(x) <= ggg$target)]
      dat[, value := rowMeans(.SD), .SDcols = paste0('prob',seq(dcols))]

      prob_ad = rbind(prob_ad, dat[, .(ADM1_CODE, year, value)])

    }

    write.csv(prob_ad, '<<<< FILEPATH REDACTED >>>>', row.names = F)

  }
}
