# ---------------------------------------------------------------------------------------------
# plot_stackers_by_adm01

# Function that pulls in aggregated results, data, and stackers
# and plots them for admin 0 and admin 1 over time
#
# Inputs:
# indicator - modeled indicator
# indicator_group - modeled indicator group
# run_date - run date corresponding to model run
# reg - region modeling (NOTE: the function only takes one at a time)
# measure - the name you'd like for your measure, e.g. 'prevalence' or 'proportion'
# draws - whether or not to use draws
# raked - whether or not to plot raked estimates
# credible_interval - credible interval to plot for the mean unraked estimates
# N_breaks - breaks you'd like to use to plot points in proportion to N
# admin_data - optional
#              if NULL:
#                will aggregate data pulled from <<<< FILEPATH REDACTED >>>>
#              if proivded, must be:
#                list of aggregated admin 0 and admin 1 data with data.table/data.frame 
#                elements named ad0 and ad1, respectively, and the following collumns:
#                ad0: svy_id, ADM0_NAME, ADM0_CODE, year, outcome, N
#                ad1: svy_id, ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_CODE, year, outcome, N
#                which can be created using mbg_central function input_aggregate_admin()
#
# Outputs:
# plots of data, stackers, mean unraked estimates with upper and lower credible intervals,
# and mean raked estimates
# ---------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# Start function
plot_stackers_by_adm01 <- function(indicator = indicator,
                                   indicator_group = indicator_group,
                                   run_date = run_date,
                                   regions = Regions,
                                   measure,
                                   draws = T,
                                   raked = T,
                                   credible_interval = 0.95,
                                   N_breaks = c(0, 10, 50, 100, 500, 1000, 2000, 4000),
                                   admin_data = NULL,
                                   data_tag,
                                   rm_yemen = F) {
  # -----------------------------------------------------------------------------------
  
  
  # --------------------------
  # Load packages
  library('magrittr')
  library('data.table')
  library('rgeos')
  library('sp')
  library('raster')
  library('ggplot2')
  library('rgdal')
  library('ggplot2')
  library('ggthemes')
  # --------------------------
  
  # ------------------------------------------------------------------------------------------------------------------------
  # Load functions
  
  # Get outputs
  source('<<<< FILEPATH REDACTED >>>>')
  
  # Get model information from RData objects
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
  
  # Summarize model results by admin
  summarize_admin = function(outputdir, measure, admin_level, admin_data, draws, raked, credible_interval){
    
    if(raked){
      admin = paste0(outputdir, indicator,'_', measure, '_raked_admin_draws_eb_bin0_0.RData')
    } else{
      admin = paste0(outputdir, indicator, '_', measure, '_unraked_admin_draws_eb_bin0_0.RData')
    }
    
    load(admin)
    
    admin = copy(get(paste0('admin_',admin_level)))
    
    draw_cols = grep('V[0-9]*', names(admin), value = T)
    
    #create counts
    #get rates
    admin[, paste0('mean_rate') := rowMeans(.SD, na.rm = T), .SDcols = draw_cols]
    admin[, paste0('lower_rate') := apply(.SD[,draw_cols,with=F], 1, quantile, probs = (1-credible_interval)/2, na.rm = T)]
    admin[, paste0('upper_rate') := apply(.SD[,draw_cols,with=F], 1, quantile, probs = 1-(1-credible_interval)/2, na.rm = T)]
    
    setnames(admin,draw_cols, paste0(draw_cols,'_rate'))
    
    #get counts
    admin[, paste0(c('mean','lower','upper',draw_cols), '_count') := lapply(c('mean','lower','upper', draw_cols), function(x)
      get(paste0(x,'_rate')) * pop)]
    
    admin[, measure := measure]
    
    dcs = c(paste0(draw_cols, '_rate'), paste0(draw_cols,'_count'))
    ndcs = names(admin)[!names(admin) %in% dcs]
    
    setcolorder(admin,c(ndcs, dcs))
    
    if(!draws){
      admin[,(dcs):=NULL] 
    }
    
    # add admin names
    link <- admin_data[['ad1']]
    kcols = expand.grid(adlev = 0:admin_level, type = c('NAME','CODE'))
    kcols = paste0('ADM',kcols$adlev,'_',kcols$type)
    link = unique(link[,kcols, with = F])
    
    merge_vars = intersect(names(admin), names(link))
    
    admin = merge(admin, link, by = merge_vars)
    
    return(admin)
  }
  # ------------------------------------------------------------------------------------------------------------------------
  
  
  # ------------------------------------------------------------------------------------------------------------------------
  # Load and clean data and model information
  # load and aggregate data if it is not given
  if (is.null(admin_data)) {
    message(paste0('You have not provided admin data. Aggregating data from <<<< FILEPATH REDACTED >>>>'))
    input_data <- fread('<<<< FILEPATH REDACTED >>>>')
    admin_data <- input_aggregate_admin(indicator = indicator, 
                                        indicator_group = indicator_group, 
                                        reg = regions, 
                                        run_date = run_date,
                                        input_data = input_data,
                                        sample_column = 'sum_of_sample_weights')
  }
  
  # set model output directories
  message('Loading and cleaning model information')
  outputdir = paste0('<<<< FILEPATH REDACTED >>>>')
  imagedir = paste0('<<<< FILEPATH REDACTED >>>>')
  
  # loop over regions
  for (reg in regions) {
    message(reg)
    
    # load submodel names
    config = fetch_from_rdata('<<<< FILEPATH REDACTED >>>>', 'config')
    
    # check to see if new stacking was used
    contains_object = function(file, obj){
      load(file)
      return(exists(obj))
    }
    
    if(contains_object(paste0(imagedir, 'pre_run_tempimage_', run_date, '_bin0_', reg, '_0.RData'), 'stacking_models')){
      stack = fetch_from_rdata(paste0(imagedir, 'pre_run_tempimage_', run_date, '_bin0_', reg, '_0.RData'), 'stacking_models')
      submodels = unlist(lapply(stack, function(x) x$model_name))
    }else{
      submodels = trimws(strsplit(config[V1 == 'stacked_fixed_effects',V2], '+', fixed = T)[[1]])
    }
    
    # set naming
    if (raked) {
      varnames = data.table(variable = c(submodels, 'stacking', 'Unraked', 'Raked'),
                            var_name = c(submodels, 'stacking', 'Unraked', 'Raked'))
    } else {
      varnames = data.table(variable = c(submodels, 'stacking', 'Unraked'),
                            var_name = c(submodels, 'stacking','Unraked'))
    }
    
    # load coefficients
    mod = lapply(paste0(outputdir, indicator, '_model_eb_bin0_',reg,'_0.RData'), 
                 function(x) fetch_from_rdata(x, 'res_fit'))
    if (as.logical(use_inla_country_fes)) {
      # fixed effects
      mod = mod[[1]]
      coefs = data.table(child = submodels, coef = invlogit(mod$summary.fixed$mean[2:4]), region = reg)
    } else {
      # random effects
      mod = lapply(mod, function(x) data.table(child = submodels, coef = x$summary.random$covar$mean))
      coefs = rbindlist(lapply(seq(reg), function(x) mod[[x]][,region:=reg[x]]))
    }
    
    
    # get unraked admin1 results
    mbg = summarize_admin(outputdir, measure, 1, admin_data, draws, raked = F, credible_interval)
    
    #add stackers weighted by their coefficient
    mbg <- as.data.frame(mbg)
    mbg$stacking <- 0
    for (n in 1:length(submodels)){
      child_name <- submodels[n]
      child_col <- mbg[,child_name]
      child_coef <- coefs[child == submodels[n]]$coef
      mbg$stacking <- child_col*child_coef + mbg$stacking
    }
    
    mbg <- as.data.table(mbg)
    
    # reshape mbg long
    mbg = mbg[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'mean_rate', 'lower_rate', 'upper_rate', submodels, 'stacking'), with = F]
    mbg = melt(mbg, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'lower_rate', 'upper_rate'), variable.factor = F)
    mbg[variable == 'mean_rate', variable:= 'Unraked']
    
    if (raked) {
      # get raked admin1 results
      mbg_raked = summarize_admin(outputdir, measure, 1, admin_data, draws, raked = T, credible_interval)
      
      #add stackers
      mbg_raked <- as.data.frame(mbg_raked)
      mbg_raked$stacking <- 0
      for (n in 1:length(submodels)){
        child_name <- submodels[n]
        child_col <- mbg_raked[,child_name]
        child_coef <- coefs[child == submodels[n]]$coef
        mbg_raked$stacking <- child_col*child_coef + mbg_raked$stacking
      }
      
      mbg_raked <- as.data.table(mbg_raked)
      
      # reshape mbg_raked long
      mbg_raked = mbg_raked[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'mean_rate', 'lower_rate', 'upper_rate', submodels), with = F]
      mbg_raked = melt(mbg_raked, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'lower_rate', 'upper_rate'), variable.factor = F)
      mbg_raked[variable == 'mean_rate', variable:= 'Raked']
      
      # add to data
      mbg = rbindlist(list(mbg, mbg_raked), use.names = TRUE)
      rm(mbg_raked)
    }
    
    # get unraked admin0 results
    mbg0 = summarize_admin(outputdir, measure, 0, admin_data, draws, raked = F, credible_interval)
    
    #add stackers weighted by their coefficient
    mbg0 <- as.data.frame(mbg0)
    mbg0$stacking <- 0
    for (n in 1:length(submodels)){
      child_name <- submodels[n]
      child_col <- mbg0[,child_name]
      child_coef <- coefs[child == submodels[n]]$coef
      mbg0$stacking <- child_col*child_coef + mbg0$stacking
    }
    
    mbg0 <- as.data.table(mbg0)
    
    # reshape mbg long
    mbg0 = mbg0[, c('ADM0_CODE', 'ADM0_NAME', 'year', 'mean_rate', 'lower_rate', 'upper_rate', submodels, 'stacking'), with = F]
    mbg0 = melt(mbg0, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'year', 'lower_rate', 'upper_rate'), variable.factor = F)
    mbg0[variable == 'mean_rate', variable:= 'Unraked']
    
    if (raked) {
      # get raked admin1 results
      mbg0_raked = summarize_admin(outputdir, measure, 0, admin_data, draws, raked = T, credible_interval)
      
      #add stackers
      mbg0_raked <- as.data.frame(mbg0_raked)
      mbg0_raked$stacking <- 0
      for (n in 1:length(submodels)){
        child_name <- submodels[n]
        child_col <- mbg0_raked[,child_name]
        child_coef <- coefs[child == submodels[n]]$coef
        mbg0_raked$stacking <- child_col*child_coef + mbg0_raked$stacking
      }
      
      mbg0_raked <- as.data.table(mbg0_raked)
      
      # reshape mbg_raked long
      mbg0_raked = mbg0_raked[, c('ADM0_CODE', 'ADM0_NAME', 'year', 'mean_rate', 'lower_rate', 'upper_rate', submodels), with = F]
      mbg0_raked = melt(mbg0_raked, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'year', 'lower_rate', 'upper_rate'), variable.factor = F)
      mbg0_raked[variable == 'mean_rate', variable:= 'Raked']
      
      # add to data
      mbg0 = rbindlist(list(mbg0, mbg0_raked), use.names = TRUE)
      rm(mbg0_raked)
    }
    
    # get unraked admin0 results
    adm0 = summarize_admin(outputdir, measure, 0, admin_data, draws, raked = F, credible_interval)
    
    #add stackers
    adm0 <- as.data.frame(adm0)
    adm0$stacking <- 0
    for (n in 1:length(submodels)){
      child_name <- submodels[n]
      child_col <- adm0[,child_name]
      child_coef <- coefs[child == submodels[n]]$coef
      adm0$stacking <- child_col*child_coef + adm0$stacking
    }
    
    adm0 <- as.data.table(adm0)
    
    # reshape adm0 long
    adm0 = adm0[, c('ADM0_CODE', 'ADM0_NAME', 'year', 'mean_rate', 'lower_rate', 'upper_rate', submodels, 'stacking'), with = F]
    adm0 = melt(adm0, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'year', 'lower_rate', 'upper_rate'), variable.factor = F)
    adm0[variable == 'mean_rate', variable:= 'Unraked']
    
    if (raked) {
      # get raked admin0 results
      adm0_raked = summarize_admin(outputdir, measure, 0, admin_data, draws, raked = T, credible_interval)
      
      #add stackers
      adm0_raked <- as.data.frame(adm0_raked)
      adm0_raked$stacking <- 0
      for (n in 1:length(submodels)){
        child_name <- submodels[n]
        child_col <- adm0_raked[,child_name]
        child_coef <- coefs[child == submodels[n]]$coef
        adm0_raked$stacking <- child_col*child_coef + adm0_raked$stacking
      }
      
      adm0_raked <- as.data.table(adm0_raked)
      
      # reshape adm0_raked long
      adm0_raked = adm0_raked[, c('ADM0_CODE', 'ADM0_NAME', 'year', 'mean_rate', 'lower_rate', 'upper_rate', submodels, 'stacking'), with = F]
      adm0_raked = melt(adm0_raked, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'year', 'lower_rate', 'upper_rate'), variable.factor = F)
      adm0_raked[variable == 'mean_rate', variable:= 'Raked']
      
      # add to data
      adm0 = rbindlist(list(adm0, adm0_raked), use.names = TRUE)
      rm(adm0_raked)
    }
    # ------------------------------------------------------------------------------------------------------------------------
    
    
    # ------------------------------------------------------------------------------------------------------------------------
    # Make plots
    message('Making plots')
    
    # create pdf to write to
    dir.create(paste0(outputdir, 'diagnostic_plots/'))
    pdf(paste0(outputdir, 'diagnostic_plots/admin_stacker_line_plots_', reg, '.pdf'), height = 10, width =14)
    
    # load the aggregated data that went into the model
    ad0_dat = admin_data$ad0
    ad1_dat = admin_data$ad1
    
    # add breaks for plotting points in proportion to N
    minn = min(c(ad0_dat$N,ad1_dat$N))
    maxn = max(c(ad0_dat$N,ad1_dat$N))
    
    ad0_dat[, Ncut := cut(N, N_breaks)]
    ad1_dat[, Ncut := cut(N, N_breaks)]
    ad0_dat[, Nrd := round(N, 2)]
    ad1_dat[, Nrd := round(N,2)]
    if(nrow(ad0_dat) != nrow(unique(ad1_dat[,.(svy_id,ADM0_CODE, year)]))) warning('The number of years of admin 0 data does not match number of years of admin 1 data.')
    
    #for each region in Regions
    # make a 6 panel plot with all countries in the region
    # get a0 codes
    country_codes <- get_adm0_codes(reg)
    a0_gd = mbg0[ADM0_CODE %in% country_codes,]
    a0_gd = merge(a0_gd, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
    a0_gd = merge(a0_gd, varnames)
    a0_gd[, ADM0_NAME := as.character(ADM0_NAME)]
    a0_gd[round(coef, 3) == 0, var_name := paste(var_name, signif(coef, 3)) ]
    a0_gd[!is.na(coef), var_name := paste(var_name, round(coef, 3)) ]
    submodel_names <- unique(grep("[.]", a0_gd$var_name, value=TRUE))
    a0_gd$var_name <- factor(a0_gd$var_name, levels = c(submodel_names, 'stacking','Unraked','Raked'))
    
    maxval = max(c(a0_gd[ADM0_CODE %in% country_codes,upper_rate], adm0[ADM0_CODE %in% country_codes,value], ad0_dat[ADM0_CODE %in% country_codes,outcome]), na.rm = TRUE)
    if (maxval > 1){
      maxval <- 0.01
    }
    childplot = ggplot(a0_gd, aes(x = year, y = value, color = ADM0_NAME)) + geom_line() + theme_bw() +
      ggtitle(paste(toupper(indicator), " Model Results for", reg, 'at admin 0 level')) +ylim(0, maxval) +
      facet_wrap(~var_name) + theme(legend.position="bottom")
    
    if(nrow(ad0_dat[ADM0_CODE %in% country_codes,])>0){
      childplot = childplot + geom_point(data = ad0_dat[ADM0_CODE %in% country_codes,], mapping = aes(x = year, y = outcome, size = N)) +
        scale_size_continuous(breaks = N_breaks, limits = c(head(N_breaks, n = 1), tail(N_breaks, n = 1)))
    }
    
    plot(childplot)
  
    # for each country in the region
    for(admin_id in intersect(get_adm0_codes(reg), unique(mbg[,ADM0_CODE]))){
      
      adname = unique(mbg[ADM0_CODE == admin_id, ADM0_NAME])
      message(adname)
      
      if (adname == 'Maldives') {
        
        message('Skipping Maldives plot')
        
      } else if (adname == 'Yemen') {
        message('Skipping Yemen plot')
      } else {
        
        a0_gd = adm0[ADM0_CODE == admin_id,]
        a0_gd$ADM0_NAME <- NULL
        
        a0_gd = merge(a0_gd, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
        a0_gd = merge(a0_gd, varnames)
        
        a0_gd[round(coef, 3) == 0, var_name := paste(var_name, signif(coef, 3)) ]
        a0_gd[!is.na(coef), var_name := paste(var_name, round(coef, 3)) ]
        
        # make an admin 0 plot of the stackers and whatnot
        maxval = max(c(adm0[ADM0_CODE == admin_id,upper_rate], adm0[ADM0_CODE == admin_id,value], ad0_dat[ADM0_CODE == admin_id,outcome]), na.rm = TRUE)
        if (maxval > 1){
          maxval <- 0.01
        }
        
        adplot = ggplot() + geom_line(data = a0_gd, aes(x = year, y = value, color = var_name), size = 1.2) +
          geom_ribbon(data = a0_gd[variable == 'Unraked'], aes(x = year, ymin = lower_rate, ymax = upper_rate), na.rm = TRUE, alpha = 0.2, fill = 'purple') + 
          theme_bw() + ggtitle(paste(toupper(indicator), " Model Results for", adname, 'in', reg)) + ylim(0, maxval)
        
        if(nrow(ad0_dat[ADM0_CODE == admin_id,])>0){
          adplot = adplot + geom_text(data = ad0_dat[ADM0_CODE == admin_id,], mapping = aes(x = year, y = outcome, label = Nrd))
        }
        
        plot(adplot)
        
        # make an admin 1 plot by child model
        a1_gd = mbg[ADM0_CODE == admin_id,]
        
        a1_gd = merge(a1_gd, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
        a1_gd = merge(a1_gd, varnames)
        a1_gd[, ADM1_NAME := as.character(ADM1_NAME)]
        a1_gd[round(coef, 3) == 0, var_name := paste(var_name, signif(coef, 3)) ]
        a1_gd[!is.na(coef), var_name := paste(var_name, round(coef, 3)) ]
        submodel_names <- unique(grep("[.]", a1_gd$var_name, value=TRUE))
        a1_gd$var_name <- factor(a1_gd$var_name, levels = c(submodel_names, 'stacking','Unraked','Raked'))
        
        maxval = max(c(a1_gd[ADM0_CODE == admin_id,upper_rate], adm0[ADM0_CODE == admin_id,value], ad1_dat[ADM0_CODE == admin_id,outcome]), na.rm = TRUE)
        if (maxval > 1){
          maxval <- 0.01
        }
        childplot = ggplot(a1_gd, aes(x = year, y = value, color = ADM1_NAME)) + geom_line() + theme_bw() +
          ggtitle(paste(toupper(indicator), " Model Results for", adname, 'at admin 1 level')) +ylim(0, maxval) +
          facet_wrap(~var_name) + theme(legend.position="bottom")
        
        if(nrow(ad0_dat[ADM0_CODE == admin_id,])>0){
          childplot = childplot + geom_point(data = ad1_dat[ADM0_CODE == admin_id,], mapping = aes(x = year, y = outcome, size = N)) +
            scale_size_continuous(breaks = N_breaks, limits = c(head(N_breaks, n = 1), tail(N_breaks, n = 1)))
        }
        
        plot(childplot)
        
        # make a 6 panel plot by child model for each admin 1
        for(a1 in unique(a1_gd$ADM1_CODE)){
          ad1name = unique(a1_gd[ADM1_CODE==a1, ADM1_NAME])
          
          # check for NA's
          check <- mbg[ADM1_NAME == ad1name]
          ad1name = unique(a1_gd[ADM1_CODE==a1, ADM1_NAME])
          
          # check for NA's
          check <- mbg[ADM1_NAME == ad1name]
          
          if (all(is.na(mean(check$value)))) {
            
            write(paste0(ad1name, ' in ', adname, '; '), 
                  file = paste0(outputdir, 'diagnostic_plots/ad1_missing_estimates.txt'),
                  append = TRUE)
            
          } else {
            a1_data = mbg[ADM1_CODE == a1,]
            
            a1_data = merge(a1_data, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
            a1_data = merge(a1_data, varnames)
            a1_data[, ADM1_NAME := as.character(ADM1_NAME)]
            
            a1_data[round(coef, 3) == 0, var_name := paste(var_name, signif(coef, 3)) ]
            a1_data[!is.na(coef), var_name := paste(var_name, round(coef, 3)) ]
            submodel_names <- unique(grep("[.]", a1_data$var_name, value=TRUE))
            a1_data$var_name <- factor(a1_data$var_name, levels = c(submodel_names, 'stacking','Unraked','Raked'))
            
            maxval = max(c(a1_data[ADM0_CODE == admin_id,upper_rate], adm0[ADM0_CODE == admin_id,value], ad1_dat[ADM0_CODE == admin_id,outcome]), na.rm = TRUE)
            if (maxval > 1){
              maxval <- 0.01
            }
            childplot = ggplot(a1_data, aes(x = year, y = value, color = ADM1_NAME)) + geom_line() + theme_bw() +
              ggtitle(paste0(toupper(indicator), " Model Results for ", ad1name, " (", adname, ")")) +ylim(0, maxval) +
              facet_wrap(~var_name) + theme(legend.position="bottom")
            
            if(nrow(ad0_dat[ADM0_CODE == admin_id,])>0){
              childplot = childplot + geom_point(data = ad1_dat[ADM1_CODE == a1,], mapping = aes(x = year, y = outcome, size = N)) +
                scale_size_continuous(breaks = N_breaks, limits = c(head(N_breaks, n = 1), tail(N_breaks, n = 1)))
            }
            
            plot(childplot)
          }
          
          ad1name = unique(a1_gd[ADM1_CODE==a1, ADM1_NAME])
          
          # check for NA's
          check <- mbg[ADM1_NAME == ad1name]
          
          if (all(is.na(mean(check$value)))) {
            
            write(paste0(ad1name, ' in ', adname, '; '), 
                  file = paste0(outputdir, 'diagnostic_plots/ad1_missing_estimates.txt'),
                  append = TRUE)
            
          } else {
            
            # make an admin 1 plot of the stackers and whatnot
            maxval = max(c(a1_gd[ADM1_CODE==a1,upper_rate], adm0[ADM0_CODE == admin_id,value], ad1_dat[ADM1_CODE == a1,outcome]), na.rm = TRUE)
            if (maxval > 1){
              maxval <- 0.01
            }
            
            a1_plot = ggplot() + geom_line(data = a1_gd[ADM1_CODE==a1,], aes(x = year, y = value, color = var_name), size = 1.2) +
              geom_ribbon(data = a1_gd[ADM1_CODE==a1 & variable == 'Unraked',], aes(x = year, ymin = lower_rate, ymax = upper_rate), na.rm = TRUE, alpha = 0.2, fill = 'purple') + 
              theme_bw() + ggtitle(paste(toupper(indicator), ' Model Results for:',ad1name,'|',adname, '|', a1)) + ylim(0, maxval)
            
            if(nrow(ad1_dat[ADM1_CODE == a1,])>0){
              a1_plot = a1_plot + geom_text(data = ad1_dat[ADM1_CODE == a1,], mapping = aes(x = year, y = outcome, label = Nrd))
            }
            
            plot(a1_plot)
          }
          
        }
          
      }
      
    }
    
    dev.off()
    message(paste0('Plots saved at ', outputdir, 'diagnostic_plots/admin_stacker_line_plots_', reg, '.pdf'))
    # ------------------------------------------------------------------------------------------------------------------------
  } # end loop over regions
  
  
  # ------------------------------------------------------------------------------------------------------------------------
  # End function
  return('Stacker line plots complete!')
  
}
# ------------------------------------------------------------------------------------------------------------------------

