fetch_from_rdata = function(file_location, item_name, use_grep = F){
  load(file_location)
  
  if(use_grep){
    ret_obj = lapply(item_name, function(x) mget(grep(x, ls(), value=T)))
  }else{
    ret_obj = lapply(item_name, function(x) get(x))
  }
  
  if(length(item_name)==1){
    ret_obj = ret_obj[[1]]
  }
  
  return(ret_obj)
  
}

build_africa_raster = function(ras_files, save = T, rake = T){
  ras = lapply(ras_files, raster::brick)
  af_brick = do.call(raster::merge, ras)
  
  rakeaddin = '<<<< FILEPATH REDACTED >>>>'
  
  if(save){
    bn = basename(ras_files[[1]])
    if(rake){
      metric = substring(bn, regexpr('_raked_',bn, fixed = T)[1]+7, regexpr('_raster',bn, fixed = T)[1]-1)
    }else{
      metric = strsplit(bn, "_", fixed = T)[[1]]
      metric = metric[[length(metric)-1]]
    }
    pos_meas = c('prevalence','mortality','incidence',"daly","yll","yld")
    measure = pos_meas[sapply(pos_meas, function(x) grepl(x,bn))]
    
    
    outfolder_j = '<<<< FILEPATH REDACTED >>>>'
    dir.create(outfolder_j, recursive = T)
    
    outfolder = dirname(ras_files[[1]])
    
    
    dest_name = '<<<< FILEPATH REDACTED >>>>'
    dest_name_j = '<<<< FILEPATH REDACTED >>>>'
    writeRaster(af_brick, filename = dest_name, overwrite = T)
    writeRaster(af_brick, filename = dest_name_j, overwrite = T)
    return(dest_name)
  }else{
    return(af_brick)
  }
}


oos = function(run_date, folds = 0:5, post_proc = T, save = T){
  if(length(folds)>1){
    
    stratum_ho = fetch_from_rdata('<<<< FILEPATH REDACTED >>>>', 'stratum_ho')
    
    xyt = rbindlist(stratum_ho)
    
  } else{
    region_list = fetch_from_rdata('<<<< FILEPATH REDACTED >>>>', 'region_list')
    xyt = rbindlist(lapply(region_list, function(x) {
      pathaddin = paste0('_bin',0,'_',x,'_',0)
      df = fetch_from_rdata('<<<< FILEPATH REDACTED >>>>', 'df')
      df[, region := x]
    }))
    
    xyt[, fold:=0]
    
  }
  message(run_date)
  res = get_is_oos_draws(xyt = xyt , xytfr_cols = c('longitude','latitude','year','fold','region'),
                         indicator_group = indicator_group, indicator = indicator,run_date = run_date, age = 0,
                         year_list = year_list,folds = folds)
  
  dat = cbind(xyt, res)
  
  if(post_proc){
    #add admin 2
    admin_level <- 2
    #load the shape
    shapes <- readRDS('<<<< FILEPATH REDACTED >>>>')
    
    
    locs <- SpatialPoints(cbind(dat$longitude, dat$latitude), proj4string = CRS(proj4string(shapes)))
    adm.df <- sp::over(locs, shapes)
    
    dat = cbind(dat, adm.df[,names(adm.df)[!names(adm.df) %in% names(dat)]])
    
    insnames = grep('is\\.[0-9]', names(dat), value = T)
    newnames = paste0('draw',seq_along(insnames))
    
    insamp = copy(dat)
    insamp[, fold:=0]
    setnames(insamp, insnames, newnames)
    if(!all(folds == 0)){
      oosnames = grep('oos\\.[0-9]', names(dat), value = T)
      insamp[, (oosnames) := NULL]
      dat[, (insnames) := NULL]
      setnames(dat,oosnames, newnames)
      dat = rbind(dat, insamp)
    }else{
      dat = insamp
    }
    
  }
  
  if(save){
    saveRDS(dat, '<<<< FILEPATH REDACTED >>>>')
  }
  
  return(dat)
  
}

summarize_oos_results = function(draws.df){
  
  hold = copy(draws.df)
  
  country.pvtable <- get_pv_table(d = draws.df,
                                  indicator_group = indicator_group,
                                  rd = run_date,
                                  indicator=indicator,
                                  aggregate_on='ad0',
                                  draws = as.numeric(samples),
                                  out.dir = '<<<< FILEPATH REDACTED >>>>',
                                  plot_ci = T)
  
  write.csv(country.pvtable,
            file = '<<<< FILEPATH REDACTED >>>>')
  
  draws.df = copy(hold)
  
  ad1.pvtable <- get_pv_table(d = draws.df,
                              indicator_group = indicator_group,
                              rd = run_date,
                              indicator=indicator,
                              aggregate_on='ad1',
                              draws = as.numeric(samples),
                              out.dir = '<<<< FILEPATH REDACTED >>>>',
                              plot_ci = T)
  write.csv(ad1.pvtable,
            file = '<<<< FILEPATH REDACTED >>>>')
  
  draws.df = copy(hold)
  ad2.pvtable <- get_pv_table(d = draws.df,
                              indicator_group = indicator_group,
                              rd = run_date,
                              indicator=indicator,
                              aggregate_on='ad2',
                              draws = as.numeric(samples),
                              out.dir = '<<<< FILEPATH REDACTED >>>>',
                              plot_ci = T)
  write.csv(ad2.pvtable,
            file = '<<<< FILEPATH REDACTED >>>>')
  
  draws.df = copy(hold)
  ho.pvtable <- get_pv_table(d = draws.df,
                             indicator_group = indicator_group,
                             rd = run_date,
                             indicator=indicator,
                             aggregate_on='ho_id',
                             draws = as.numeric(samples),
                             out.dir = '<<<< FILEPATH REDACTED >>>>',
                             plot_ci = T)
  write.csv(ho.pvtable,
            file = '<<<< FILEPATH REDACTED >>>>')
  
  return(invisible())
  
}

summarize_admin = function(measure, admin_level, draws = F, raked = T){
  
  inme = '<<<< FILEPATH REDACTED >>>>'
  outme = '<<<< FILEPATH REDACTED >>>>'
  
  if(raked){
    admin = '<<<< FILEPATH REDACTED >>>>'
  } else{
    admin = '<<<< FILEPATH REDACTED >>>>'
  }
  
  load(admin)
  
  admin = copy(get(paste0('admin_',admin_level)))
  
  draw_cols = grep('V[0-9]*', names(admin), value = T)
  
  #create counts
  #get rates
  admin[, paste0('mean_rate') := rowMeans(.SD), .SDcols = draw_cols]
  admin[, paste0('lower_rate') := apply(.SD[,draw_cols,with=F], 1, quantile, probs = 0.025, na.rm = T)]
  admin[, paste0('upper_rate') := apply(.SD[,draw_cols,with=F], 1, quantile, probs = 0.975, na.rm = T)]
  
  setnames(admin,draw_cols, paste0(draw_cols,'_rate'))
  
  #get counts
  admin[, paste0(c('mean','lower','upper',draw_cols), '_count') := lapply(c('mean','lower','upper', draw_cols), function(x)
    get(paste0(x,'_rate')) * pop * pop_scalar)]
  
  admin[, measure := measure]
  
  dcs = c(paste0(draw_cols, '_rate'), paste0(draw_cols,'_count'))
  ndcs = names(admin)[!names(admin) %in% dcs]
  
  setcolorder(admin,c(ndcs, dcs))
  
  if(!draws){
    admin[,(dcs):=NULL] 
  }
  
  link = readRDS('<<<< FILEPATH REDACTED >>>>')
  kcols = expand.grid(adlev = 0:admin_level, type = c('NAME','CODE'))
  kcols = paste0('ADM',kcols$adlev,'_',kcols$type)
  link = unique(link[,kcols, with = F])
  
  merge_vars = intersect(names(admin), names(link))
  
  admin = merge(admin, link, by = merge_vars)
  
  
  if(nrow(admin[ADM0_NAME == 'Cape Verde' & year == 2016 & is.na(mean_rate),])>0){
    cpv = admin[ADM0_NAME == 'Cape Verde' & year ==2015,]
    cpv[,year:= 2016]
    admin = admin[!(ADM0_NAME == 'Cape Verde' & year ==2016),]
    admin = rbind(admin,cpv)
  }
  
  return(admin)
}

make_new_ras = function(reg, cell_pred_file, condition = F, condition_val){
  print(cell_pred_file)
  stopifnot(file.exists(cell_pred_file))
  
  #get simple raster
  gaul_list           <- get_gaul_codes(reg)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
  subset_shape        <- simple_polygon_list[[1]]
  simple_polygon      <- simple_polygon_list[[2]]
  
  ## Load list of raster inputs (pop and simple)
  raster_list        <- build_simple_raster_pop(subset_shape)
  simple_raster      <- raster_list[['simple_raster']]
  
  #load the aroc cell pred
  cell_pred = readRDS(cell_pred_file)
  
  if(condition){
    print('blarg')
    cell_pred = cell_pred <=condition_val
  }
  
  #summarize
  if(class(cell_pred) == 'matrix'){
    cell_pred = rowMeans(cell_pred)
  }
  
  ras = newInsertRaster(simple_raster, cell_pred)
  
  return(ras)
  
}

newInsertRaster <- function (ras, new_vals) {
  
  
  # calculate cell index if not provided
  idx <- cellIdx(ras)
  
  # check the index makes superficial sense
  stopifnot(max(idx) <= ncell(ras))
  
  # create results raster
  
  stopifnot(length(new_vals)%%length(idx) == 0 )
  
  nlay = length(new_vals)/length(idx)
  raster_new <- raster::brick(replicate(nlay,
                                        ras[[1]],
                                        simplify = FALSE))
  
  raster_new[] = NA
  
  # update the values
  #get the indices
  idxs = unlist(lapply(1:nlay, function(x) idx + ((x - 1) * length(ras))))
  
  stopifnot(length(idxs) == length(new_vals))
  
  values(raster_new)[idxs] = new_vals
  
  return(raster_new)
  
}

get_running_jobs = function(){
  qstat = system2('qstat', stdout = T)
  qstat = strsplit(qstat,split = '\\s+')
  
  #remove second row
  #remake the headings
  qstat = qstat[c(-1,-2)]
  
  qstat = lapply(qstat, function(x) data.table(t(x)))
  
  qstat = lapply(qstat, function(x){
    if(ncol(x) ==9){
      x[,V10:=V9]
      x[,V9:= ""]
    }
  })
  
  
  qstat = rbindlist(qstat)
  qstat[,V1:=NULL]
  
  setnames(qstat, c('jobid','prior','name','user','state', 'date_start', 'time_start', 'node','slots'))
  
  
  return(qstat)
  
}