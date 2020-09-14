# ##################################################################################################
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
# N_breaks - breaks you'd like to use to plot points in proportion to N, can be number of total breaks or vector
#
# Outputs:
# plots of data, stackers, mean unraked estimates with upper and lower credible intervals,
# and mean raked estimates (optionally) over time
#########################################################################################################

plot_stackers_by_adm01 <- function(admin_data,
                                   shapefile_version = raking_shapefile_version,
                                   indicator = indicator,
                                   indicator_group = indicator_group,
                                   run_date = run_date,
                                   regions = Regions,
                                   measure = 'prevalence',
                                   raked = T,
                                   credible_interval = 0.95,
                                   force_plot = FALSE,
                                   two_xgboost) {
  
  # (1) Setup -----------------------------------------------------------------------------------------------
  
  # packages
  library('magrittr')
  library('data.table')
  library('rgeos')
  library('sp')
  library('raster')
  library('ggplot2')
  library('rgdal')
  library('ggthemes')
  library('ggrepel')
  
  # functions
  source('<<<< FILEPATH REDACTED >>>>')
  source('<<<< FILEPATH REDACTED >>>>')
  
  # (2) Load and clean model information ---------------------------------------------------------------------
  
  # model information
  message('Loading and cleaning model information')
  outputdir = '<<<< FILEPATH REDACTED >>>>'
  imagedir = '<<<< FILEPATH REDACTED >>>>'
  
  # loop over regions
  for (reg in regions) {
    message(reg)
    
    # make sure input files are data tables
    admin_data[['ad0']] <- as.data.table(admin_data[['ad0']])
    admin_data[['ad1']] <- as.data.table(admin_data[['ad1']])
    
    # subset input data to region
    ad0_data <- admin_data[['ad0']][region == reg]
    ad1_data <- admin_data[['ad1']][region == reg]
    
    # load submodel names
    config = fetch_from_rdata(paste0(imagedir, 'pre_run_tempimage_', run_date, '_bin0_', reg, '_0.RData'), 'config')
    if(class(config) == 'list'){
      config <- config$config
    }
    submodels = trimws(strsplit(config[V1 == 'stacked_fixed_effects',V2], '+', fixed = T)[[1]])
    if (two_xgboost) submodels <- c(submodels, "xgboost2")
    
    # set variable names
    if (raked) {
      varnames = data.table(variable = c(submodels, 'Stacking', 'Unraked', 'Raked'),
                            var_name = c(submodels, 'Stacking', 'Unraked', 'Raked'))
    } else {
      varnames = data.table(variable = c(submodels, 'Stacking', 'Unraked'),
                            var_name = c(submodels, 'Stacking', 'Unraked'))
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
    
    if (force_plot == T){
      mod[[1]]$coef <- 1
      coefs$coef <- 1
    }
    
    # (3) Format admin 0 and 1 results -----------------------------------------------------------------------------------------------------
    
    # get unraked admin1 results and reshape long
    mbg = ad1_data[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'mean', 'lower', 'upper', submodels), with = F]
    
    # calculate "stacking only" values
    mbg$Stacking <- 0
    for (stacker in submodels){
      mbg$Stacking <- mbg[,..stacker]*coefs[child == stacker,coef] + mbg$Stacking
    }
    
    mbg = melt(mbg, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'lower', 'upper'), variable.factor = F)
    mbg[variable == 'mean', variable:= 'Unraked']
    
    if (raked) {
      
      # get unraked admin1 results and reshape long
      mbg_raked = ad1_data[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'mean_raked', 'lower_raked', 'upper_raked', submodels), with = F]
      mbg_raked = melt(mbg_raked, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME', 'year', 'lower_raked', 'upper_raked'), variable.factor = F)
      mbg_raked[variable == 'mean_raked', variable:= 'Raked']
      mbg_raked <- dplyr::rename(mbg_raked, lower = lower_raked, upper = upper_raked)
      
      # add to data
      mbg = rbindlist(list(mbg, mbg_raked), use.names = TRUE)
    }
    
    # get unraked admin0 results and reshape long
    adm0 = ad0_data[, c('ADM0_CODE', 'ADM0_NAME', 'year', 'mean', 'lower', 'upper', submodels), with = F]
    
    # calculate "stacking only" values
    adm0$Stacking <- 0
    for (stacker in submodels){
      adm0$Stacking <- adm0[,..stacker]*coefs[child == stacker,coef] + adm0$Stacking
    }
    
    adm0 = melt(adm0, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'year', 'lower', 'upper'), variable.factor = F)
    adm0[variable == 'mean', variable:= 'Unraked']
    
    if (raked) {
      
      # get raked admin0 results and reshape long
      adm0_raked = ad0_data[, c('ADM0_CODE', 'ADM0_NAME', 'year', 'mean_raked', 'lower_raked', 'upper_raked', submodels), with = F]
      adm0_raked = melt(adm0_raked, id.vars = c('ADM0_CODE', 'ADM0_NAME', 'year', 'lower_raked', 'upper_raked'), variable.factor = F)
      adm0_raked[variable == 'mean_raked', variable:= 'Raked']
      adm0_raked <- dplyr::rename(adm0_raked, lower = lower_raked, upper = upper_raked)
      
      # add to data
      adm0 = rbindlist(list(adm0, adm0_raked), use.names = TRUE)
      rm(adm0_raked)
    }
    
    # (4) filter out border pixel input data aggregated to wrong country --------------------------------------------------------------
    
    # read in input data
    nid_to_loc <- fread('<<<< FILEPATH REDACTED >>>>') %>%
      dplyr::select(country, nid) %>%
      unique.data.frame()
    
    iso3_link <- get_location_code_mapping(shapefile_version) %>%
      dplyr::select(ADM_CODE, ihme_lc_id)
    
    nid_to_loc <- merge(nid_to_loc, iso3_link, by.x = 'country', by.y = 'ihme_lc_id') %>%
      dplyr::rename(nid_ad0_code = ADM_CODE) %>%
      dplyr::select(-country)
    
    #merge onto data
    ad0_data <- merge(ad0_data, nid_to_loc, by = 'nid', all.x  = TRUE )
    ad0_data[,in_country_data := ifelse(ADM0_CODE == nid_ad0_code, 'yes', 'no')]
    
    ad1_data <- merge(ad1_data, nid_to_loc, by = 'nid', all.x  = TRUE )
    ad1_data[,in_country_data := ifelse(ADM0_CODE == nid_ad0_code, 'yes', 'no')]
    
    # filter data
    ad0_data <- filter(ad0_data, in_country_data == 'yes') %>%
      as.data.table()
    ad1_data <- filter(ad1_data, in_country_data == 'yes') %>%
      as.data.table()
    
    #temp fixes for NA predictions
    if (nrow(filter(adm0, is.na(lower))) == nrow(adm0)) adm0$lower <- 0
    if (nrow(filter(adm0, is.na(upper))) == nrow(adm0)) adm0$upper <- 0
    
    if (nrow(filter(mbg, is.na(lower))) == nrow(mbg)) mbg$lower <- 0
    if (nrow(filter(mbg, is.na(upper))) == nrow(mbg)) mbg$upper <- 0

    # (5) Make plots ------------------------------------------------------------------------------------------------------------------------

    message('Making plots')
    
    # create pdf to write to
    dir.create(paste0(outputdir, 'diagnostic_plots/'))
    pdf(paste0(outputdir, 'diagnostic_plots/admin_stacker_line_plots_', reg, '.pdf'),  height = 10, width =14)
    
    # load and clean the aggregated data that went into the model
    ad0_dat = ad0_data[!is.na(input_mean)]
    setnames(ad0_dat, c('input_mean', 'input_ss'), c('outcome', 'N'))
    ad1_dat = ad1_data[!is.na(input_mean)]
    setnames(ad1_dat, c('input_mean', 'input_ss'), c('outcome', 'N'))
    
    # round N for plots later
    ad0_dat[, Nrd := signif(N, 3)]
    ad1_dat[, Nrd := signif(N, 3)]
    
    if(nrow(ad0_dat) != nrow(unique(ad1_dat[,.(nid,ADM0_CODE, year)]))) warning('The number of years of admin 0 data does not match number of years of admin 1 data.')
    
    # Plot type (1): region-level, 6 panel -----------------------------------------------------------------------------------
    a0_gd = adm0
    
    a0_gd = merge(a0_gd, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
    a0_gd = merge(a0_gd, varnames)
    
    a0_gd[!is.na(coef), var_name := paste(var_name, round(coef, 2)) ]
    
    submodel_names <- unique(grep("[ ]", a0_gd$var_name, value=TRUE))
    a0_gd$var_name <- factor(a0_gd$var_name, levels = c(submodel_names, 'Stacking','Unraked','Raked'))
    
    childplot = ggplot(a0_gd, aes(x = year, y = value, color = ADM0_NAME)) + geom_line() + theme_bw() +
      ggtitle(paste(toupper(indicator), " Model Results for", reg, 'at admin 0 level')) +
      facet_wrap(~var_name) + theme(legend.position = "bottom") + ylim(0, NA)
    
    if(nrow(ad0_dat)>0){
      ad0_dat_breaks <-  seq(min(ad0_dat$N, na.rm = T), max(ad0_dat$N, na.rm = T), length.out = 5)
      childplot = childplot + geom_point(data = ad0_dat, alpha = 0.7, mapping = aes(x = year, y = outcome, size = N)) +
        scale_size_continuous(breaks = ad0_dat_breaks, name = "sample size", labels = signif(ad0_dat_breaks, digits = -3))                                                                                                                   #labels = c('in-country data', 'out-of-country data'))
    }
    
    plot(childplot)
    
    # loop over regions
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
        
        # Plot type (2): country-level, 1 panel ----------------------------------------------------------------------------------------------------------
        submodel_names <- unique(grep("[ ]", a0_gd$var_name, value=TRUE))
        a0_gd$var_name <- factor(a0_gd$var_name, levels = c(submodel_names, 'Stacking','Unraked','Raked'))
        
        adplot = ggplot() + geom_line(data = a0_gd, aes(x = year, y = value, color = var_name), size = 1.2) +
          geom_ribbon(data = a0_gd[variable == 'Unraked'], aes(x = year, ymin = lower, ymax = upper), na.rm = TRUE, alpha = 0.2, fill = 'deepskyblue') + 
          theme_bw() + ggtitle(paste(toupper(indicator), " Model Results for", adname, 'in', reg)) + ylim(0, NA) + xlim(1999, NA)
        
        if(nrow(ad0_dat[ADM0_CODE == admin_id,])>0){
          adplot = adplot + geom_label_repel(data = ad0_dat[ADM0_CODE == admin_id,], mapping = aes(x = year, y = outcome, label = Nrd))
        }
        
        plot(adplot)
        
        # Plot type (3): country-level, 6 panel -----------------------------------------------------------------------------------------------------------
        a1_gd = mbg[ADM0_CODE == admin_id,]
        
        a1_gd = merge(a1_gd, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
        a1_gd = merge(a1_gd, varnames)
        a1_gd[, ADM1_NAME := as.character(ADM1_NAME)]
        
        a1_gd[!is.na(coef), var_name := paste(var_name, round(coef, 2)) ]
        
        submodel_names <- unique(grep("[ ]", a1_gd$var_name, value=TRUE))
        a1_gd$var_name <- factor(a1_gd$var_name, levels = c(submodel_names, 'Stacking','Unraked','Raked'))
        
        childplot = ggplot(a1_gd, aes(x = year, y = value, color = ADM1_NAME)) + geom_line() + theme_bw() +
          ggtitle(paste(toupper(indicator), " Model Results for", adname, 'at admin 1 level')) + facet_wrap(~var_name) + theme(legend.position="bottom") + ylim(0, NA)
        
        if(nrow(ad1_dat[ADM0_CODE == admin_id,])>0){
          ad1_dat_breaks <-  seq(min(ad1_dat[ADM0_CODE == admin_id]$N, na.rm = T), max(ad1_dat[ADM0_CODE == admin_id]$N, na.rm = T), length.out = 5)
          childplot = childplot + geom_point(data = ad1_dat[ADM0_CODE == admin_id], alpha = 0.7, mapping = aes(x = year, y = outcome, size = N)) +
            scale_size_continuous(breaks = ad1_dat_breaks, name = "sample size", labels = signif(ad1_dat_breaks, 1))
        }
        
        plot(childplot)
        
        # loop over admin 1s
        for(a1 in unique(a1_gd$ADM1_CODE)){
          
          ad1name = unique(a1_gd[ADM1_CODE==a1, ADM1_NAME])
          
          # check for NA estimates for this admin in mbg, and write to a missing estimates list
          check <- mbg[ADM1_NAME == ad1name]
          
          if (all(is.na(mean(check$value)))) {
            
            write(paste0(ad1name, ' in ', adname, '; '), 
                  file = paste0(outputdir, 'diagnostic_plots/ad1_missing_estimates.txt'),
                  append = TRUE)
            
          } else {
            
            # Plot type (4): admin 1 level, 1 panel -----------------------------------------------------------------------------------------------------------------
            a1_plot = ggplot() + geom_line(data = a1_gd[ADM1_CODE==a1,], aes(x = year, y = value, color = var_name), size = 1.2) +
              geom_ribbon(data = a1_gd[ADM1_CODE==a1 & variable == 'Unraked',], aes(x = year, ymin = lower, ymax = upper), na.rm = TRUE, alpha = 0.2, fill = 'deepskyblue') + 
              theme_bw() + ggtitle(paste(toupper(indicator), ' Model Results for:',ad1name,'|',adname, '|', a1)) + ylim(0, NA) + xlim(1999, NA)
            
            if(nrow(ad1_dat[ADM1_CODE == a1,])>0){
              a1_plot = a1_plot + geom_label_repel(data = ad1_dat[ADM1_CODE == a1,], mapping = aes(x = year, y = outcome, label = Nrd))
            }
            
            plot(a1_plot)
            
            # Plot type (5): admin 1 level, 6 panel -----------------------------------------------------------------------------------------------------------------
            childplot = ggplot(a1_gd[ADM1_CODE==a1,], aes(x = year, y = value)) + geom_line(color = 'red') + theme_bw() +
              ggtitle(paste0(toupper(indicator), " Model Results for ", ad1name, ' (', adname, ') at admin 1 level')) +
              facet_wrap(~var_name) + theme(legend.position="bottom") + ylim(0, NA)
            
            if(nrow(ad1_dat[ADM1_CODE == a1,])>0){
              ad1_dat_breaks <-  seq(min(ad1_dat[ADM1_CODE == a1]$N, na.rm = T), max(ad1_dat[ADM1_CODE == a1]$N, na.rm = T), length.out = 5)
              if (length(unique(ad1_dat_breaks)) == 1) ad1_dat_breaks <- ad1_dat_breaks[1]
              childplot = childplot + geom_point(data = ad1_dat[ADM1_CODE == a1], alpha = 0.7, mapping = aes(x = year, y = outcome, size = N)) +
                scale_size_continuous(breaks = ad1_dat_breaks, name = "sample size", labels = signif(ad1_dat_breaks, 1))
              }
            
            plot(childplot)
          }
        }
      }
    }
    dev.off()
    message(paste0('Plots saved at ', outputdir, 'diagnostic_plots/admin_stacker_line_plots_', reg, '.pdf'))
  }
  message('stacker line plots complete!')
}


