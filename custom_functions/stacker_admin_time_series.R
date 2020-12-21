# ---------------------------------------------------------------------------------------------
# plot_stackers_by_adm01
#
# Function that pulls in aggregated results, data, and stackers
# and plots them for admin 0 and admin 1 over time
#
# Inputs:
# admin_data - must be:
#                list of aggregated admin 0 and admin 1 data, stackers, and mbg results
#                with data.table/data.frame elements named ad0 and ad1, respectively, 
#                and the following collumns:
#                  ad0: ADM0_CODE, year, ADM0_NAME, region, mean      
#                       lower, upper, cirange, pop, gam, gbm, enet [or your aggregated stackers],    
#                       nid, input_mean, input_ss
#                  ad1: ADM0_CODE, ADM1_CODE, year, ADM0_NAME, ADM1_NAME, region, mean      
#                       lower, upper, cirange, pop, gam, gbm, enet [or your aggregated stackers],    
#                       nid, input_mean, input_ss
# indicator - modeled indicator
# indicator_group - modeled indicator group
# run_date - run date corresponding to model run
# reg - region modeling (NOTE: the function only takes one at a time)
# measure - the name you'd like for your measure, e.g. 'prevalence' or 'proportion'
# raked - whether or not to plot raked estimates
# credible_interval - credible interval to plot for the mean unraked estimates
# N_breaks - breaks you'd like to use to plot points in proportion to N
#
# Outputs:
# plots of data, stackers, mean unraked estimates with upper and lower credible intervals,
# and mean raked estimates (optionally) over time saved in /share/ mbg output folder
# corresponding to indicator and run date in file /diagnositc_plots/
# ---------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Start function
plot_stackers_by_adm01 <- function(admin_data,
                                   indicator = indicator,
                                   indicator_group = indicator_group,
                                   run_date = run_date,
                                   regions = Regions,
                                   measure = 'prevalence',
                                   raked = T,
                                   credible_interval = 0.95,
                                   N_breaks = c(0, 10, 50, 100, 500, 1000, 2000, 4000)) {
  # -----------------------------------------------------------------------------------
  
  
  # -----------------------------
  # Load packages and functions
  source('<<<< FILEPATH REDACTED >>>>/mbg_central/fractional_raking_functions.R')
  library('magrittr')
  library('data.table')
  library('rgeos')
  library('sp')
  library('raster')
  library('ggplot2')
  library('rgdal')
  library('ggplot2')
  library('ggthemes')
  # -----------------------------
  
  
  # ------------------------------------------------------------------------------------------------------------------------
  # Load and clean data and model information
  
  # set model output directories
  message('Loading and cleaning model information')
  outputdir = '<<<< FILEPATH REDACTED >>>>'
  imagedir = '<<<< FILEPATH REDACTED >>>>'
  
  # loop over regions
  for (reg in regions) {
    message(reg)
    
    # make sure input files are data tables
    admin_data[['ad0']] <- as.data.table(admin_data[['ad0']])
    admin_data[['ad1']] <- as.data.table(admin_data[['ad1']])
    
    # subset input data
    ad0_data <- admin_data[['ad0']][region == reg]
    ad1_data <- admin_data[['ad1']][region == reg]
    
    # load submodel names
    config = fetch_from_rdata(paste0(imagedir, 'pre_run_tempimage_', run_date, '_bin0_', reg, '_0.RData'), 'config')
    
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
      varnames = data.table(variable = c(submodels, 'Unraked', 'Raked'),
                            var_name = c(submodels, 'Unraked', 'Raked'))
    } else {
      varnames = data.table(variable = c(submodels, 'Unraked'),
                            var_name = c(submodels, 'Unraked'))
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
    
    # get unraked admin1 results and reshape long
    mbg = ad1_data[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'mean', 'lower', 'upper', submodels), with = F]
    mbg = melt(mbg, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'lower', 'upper'), variable.factor = F)
    mbg[variable == 'mean', variable:= 'Unraked']
    
    if (raked) {
      # get unraked admin1 results and reshape long
      mbg_raked = ad1_data[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'mean_raked', 'lower_raked', 'upper_raked', submodels), with = F]
      mbg_raked = melt(mbg_raked, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'lower_raked', 'upper_raked'), variable.factor = F)
      mbg_raked[variable == 'mean_raked', variable:= 'Raked']
      
      # add to data
      setnames(mbg_raked, c('lower_raked', 'upper_raked'), c('lower', 'upper'))
      mbg = rbindlist(list(mbg, mbg_raked), use.names = TRUE)
    }
    
    # get unraked admin0 results  and reshape long
    adm0 = ad0_data[, c('ADM0_CODE', 'ADM0_NAME', 'year', 'mean', 'lower', 'upper', submodels), with = F]
    adm0 = melt(adm0, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'year', 'lower', 'upper'), variable.factor = F)
    adm0[variable == 'mean', variable:= 'Unraked']
    
    if (raked) {
      # get raked admin0 results and reshape long
      adm0_raked = ad0_data[, c('ADM0_CODE', 'ADM0_NAME', 'year', 'mean_raked', 'lower_raked', 'upper_raked', submodels), with = F]
      adm0_raked = melt(adm0_raked, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'year', 'lower_raked', 'upper_raked'), variable.factor = F)
      adm0_raked[variable == 'mean_raked', variable:= 'Raked']
      
      # add to data
      setnames(adm0_raked, c('lower_raked', 'upper_raked'), c('lower', 'upper'))
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
    
    # load and clean the aggregated data that went into the model
    ad0_dat = ad0_data[!is.na(input_mean)]
    setnames(ad0_dat, c('input_mean', 'input_ss'), c('outcome', 'N'))
    ad1_dat = ad1_data[!is.na(input_mean)]
    setnames(ad1_dat, c('input_mean', 'input_ss'), c('outcome', 'N'))
    
    # add breaks for plotting points in proportion to N
    minn = min(c(ad0_dat$N,ad1_dat$N))
    maxn = max(c(ad0_dat$N,ad1_dat$N))
    
    ad0_dat[, Ncut := cut(N, N_breaks)]
    ad1_dat[, Ncut := cut(N, N_breaks)]
    ad0_dat[, Nrd := round(N, 2)]
    ad1_dat[, Nrd := round(N,2)]
    if(nrow(ad0_dat) != nrow(unique(ad1_dat[,.(nid,ADM0_CODE, year)]))) warning('The number of years of admin 0 data does not match number of years of admin 1 data.')
    
    # for each country in the region
    for(admin_id in intersect(get_adm0_codes(reg), unique(mbg[,ADM0_CODE]))){

      adname = unique(mbg[ADM0_CODE == admin_id, ADM0_NAME])
      message(adname)
      
      if (nrow(adm0[ADM0_CODE == admin_id & !is.na(value)]) == 0) {
        
        message(paste0('Skipping ', adname, ' plot because there are not estimates here'))
        
      } else {
        
        a0_gd = adm0[ADM0_CODE == admin_id,]
        a0_gd$ADM0_NAME <- NULL
        
        a0_gd = merge(a0_gd, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
        a0_gd = merge(a0_gd, varnames)
        
        a0_gd[!is.na(coef), var_name := paste(var_name, round(coef, 2)) ]
        
        # make an admin 0 plot of the stackers and whatnot
        maxval = max(c(adm0[ADM0_CODE == admin_id,upper], adm0[ADM0_CODE == admin_id,value], ad0_dat[ADM0_CODE == admin_id,outcome]), na.rm = TRUE)
        
        adplot = ggplot() + geom_line(data = a0_gd, aes(x = year, y = value, color = var_name), size = 1.2) +
          geom_ribbon(data = a0_gd[variable == 'Unraked'], aes(x = year, ymin = lower, ymax = upper), na.rm = TRUE, alpha = 0.2, fill = 'purple') + 
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
        
        a1_gd[!is.na(coef), var_name := paste(var_name, round(coef, 2)) ]
        maxval = max(c(a1_gd[ADM0_CODE == admin_id,upper], adm0[ADM0_CODE == admin_id,value], ad1_dat[ADM0_CODE == admin_id,outcome]), na.rm = TRUE)
        childplot = ggplot(a1_gd, aes(x = year, y = value, color = ADM1_NAME)) + geom_line() + theme_bw() +
          ggtitle(paste(toupper(indicator), " Model Results for", adname, 'at admin 1 level')) +ylim(0, max(a1_gd[ADM0_CODE == admin_id,value])) +
          facet_wrap(~var_name) + theme(legend.position="bottom")
        
        if(nrow(ad0_dat[ADM0_CODE == admin_id,])>0){
          childplot = childplot + geom_point(data = ad1_dat[ADM0_CODE == admin_id,], mapping = aes(x = year, y = outcome, size = N)) +
            scale_size_continuous(breaks = N_breaks, limits = c(head(N_breaks, n = 1), tail(N_breaks, n = 1)))
        }
        
        plot(childplot)
        
        # make a plot per admin 1
        for(a1 in unique(a1_gd$ADM1_CODE)){
          
          ad1name = unique(a1_gd[ADM1_CODE==a1, ADM1_NAME])
          
          # check for NA's
          check <- mbg[ADM1_NAME == ad1name]
          
          if (all(is.na(mean(check$value)))) {
            
            write(paste0(ad1name, ' in ', adname, '; '), 
                  file = paste0(outputdir, 'diagnostic_plots/ad1_missing_estimates.txt'),
                  append = TRUE)
            
          } else {
            
            # make an admin 1 plot of the stackers and whatnot
            maxval = max(c(a1_gd[ADM1_CODE==a1,upper], adm0[ADM0_CODE == admin_id,value], ad1_dat[ADM1_CODE == a1,outcome]), na.rm = TRUE)
            
            a1_plot = ggplot() + geom_line(data = a1_gd[ADM1_CODE==a1,], aes(x = year, y = value, color = var_name), size = 1.2) +
              geom_ribbon(data = a1_gd[ADM1_CODE==a1 & variable == 'Unraked',], aes(x = year, ymin = lower, ymax = upper), na.rm = TRUE, alpha = 0.2, fill = 'purple') + 
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

