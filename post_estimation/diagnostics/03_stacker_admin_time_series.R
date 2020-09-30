# ---------------------------------------------------------------------------------------------
#REDACTED
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
#                will aggregate data pulled from #REDACTED
#              if proivded, must be:
#                list of aggregated admin 0 and admin 1 data with data.table/data.frame 
#                elements named ad0 and ad1, respectively, and the following collumns:
#                ad0: svy_id, ADM0_NAME, ADM0_CODE, year, outcome, N
#                ad1: svy_id, ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_CODE, year, outcome, N
#                which can be created using mbg_central function input_aggregate_admin()
#
# Outputs:
# plots of data, stackers, mean unraked estimates with upper and lower credible intervals,
# and mean raked estimates (optionally) over time saved in #REDACTED
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
                                   draws = T,
                                   raked = T,
                                   credible_interval = 0.95,
                                   N_breaks = c(0, 10, 50, 100, 500, 1000, 2000, 4000)) {
  # -----------------------------------------------------------------------------------
  
  
  # --------------------------
  # Load packages
  pacman::p_load(data.table, ggplot2, ggthemes, ggrepel, magrittr, raster, rgeos, rgdal, sp)
  # --------------------------
  
  # ------------------------------------------------------------------------------------------------------------------------
  # Load and clean data and model information
  
  # set model output directories
  message('Loading and cleaning model information')
  outputdir = paste0('#REDACTED', indicator_group, '/', indicator, '/output/', run_date, '/')
  imagedir = paste0('#REDACTED', indicator_group, '/', indicator, '/model_image_history/')
  
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
    minn <- min(c(ad0_dat$N,ad1_dat$N))
    maxn <- max(c(ad0_dat$N,ad1_dat$N))
    
    # plot the sample size
    ad0_dat[, Ncut := cut(N, N_breaks)]
    ad1_dat[, Ncut := cut(N, N_breaks)]
    ad0_dat[, Nrd := round(N, 2)]
    ad1_dat[, Nrd := round(N,2)]

    # add NID label to sample size plots
    ad0_dat[, N_label := paste0(svy_id, ' [', Nrd, ']')]
    ad1_dat[, N_label := paste0(svy_id, ' [', Nrd, ']')]

    if (nrow(ad0_dat) != ad1_dat[,.(svy_id,ADM0_CODE, year)] %>% unique %>% nrow) warning('The number of years of admin 0 data does not match number of years of admin 1 data.')
    
    if (reg=='dia_wssa') browser()
    
    # for each country in the region
    for (admin_id in intersect(get_adm0_codes(reg), mbg[,ADM0_CODE] %>% unique)){

      adname <- unique(mbg[ADM0_CODE == admin_id, ADM0_NAME])
      message(adname)
      
      if (nrow(adm0[ADM0_CODE == admin_id & !is.na(value)]) == 0) message(paste0('Skipping ', adname, ' plot because there are not estimates here'))
      else {
        
        a0_gd <- adm0[ADM0_CODE == admin_id,] %>% 
          .[, ADM0_NAME := NULL]
        
        a0_gd <- merge(a0_gd, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
        a0_gd <- merge(a0_gd, varnames)
        
        a0_gd[!is.na(coef), var_name := paste(var_name, round(coef, 2)) ]
        
        # make an admin 0 plot of the stackers and whatnot
        maxval <- max(c(adm0[ADM0_CODE == admin_id,upper_rate], adm0[ADM0_CODE == admin_id,value], ad0_dat[ADM0_CODE == admin_id,outcome]), na.rm = TRUE)
        
        adplot <- ggplot(data = ad0_dat[ADM0_CODE == admin_id]) + 
          geom_line(data = a0_gd, aes(x = year, y = value, color = var_name), size = 1.2) +
          geom_ribbon(data = a0_gd[variable == 'Unraked'], aes(x = year, ymin = lower_rate, ymax = upper_rate), na.rm = TRUE, alpha = 0.2, fill = 'indianred3') + 
          ylim(0, maxval) +
          scale_color_brewer(palette = 'Set1') +
          ggtitle(paste(toupper(indicator), " Model Results for", adname, 'in', reg)) + 
          theme_bw()
        
        if (nrow(ad0_dat[ADM0_CODE == admin_id,])>0) {
          adplot <- adplot + 
            geom_point(aes(x = year, y = outcome, size = N)) +
            geom_text_repel(aes(x = year, y = outcome, label = N_label),
                            size=5, segment.colour = 'gray', segment.size = .2, direction = 'x', nudge_y=.05, force=5)
        }
        
        #print to dev
        plot(adplot)
        
        # make an admin 1 plot by child model
        a1_gd <- mbg[ADM0_CODE == admin_id,]
        
        a1_gd <- merge(a1_gd, coefs[region == reg, ], by.x = 'variable', by.y = 'child', all.x = T)
        a1_gd <- merge(a1_gd, varnames)
        a1_gd[, ADM1_NAME := as.character(ADM1_NAME)]
        
        a1_gd[!is.na(coef), var_name := paste(var_name, round(coef, 2)) ]
        maxval <- max(c(a1_gd[ADM0_CODE == admin_id,upper_rate], adm0[ADM0_CODE == admin_id,value], ad1_dat[ADM0_CODE == admin_id,outcome]), na.rm = TRUE)
        childplot <- ggplot(a1_gd, aes(x = year, y = value, color = ADM1_NAME)) + 
          geom_line() + 
          facet_wrap(~var_name) + 
          ylim(0, max(a1_gd[ADM0_CODE == admin_id,value])) +
          ggtitle(paste(toupper(indicator), " Model Results for", adname, 'at admin 1 level')) +
          theme(legend.position="bottom") +
          theme_bw()
        
        if (nrow(ad0_dat[ADM0_CODE == admin_id,])>0) {
          childplot <- childplot + 
            geom_point(data = ad1_dat[ADM0_CODE == admin_id,], mapping = aes(x = year, y = outcome, size = N)) +
            scale_size_continuous(breaks = N_breaks, limits = c(head(N_breaks, n = 1), tail(N_breaks, n = 1)))
        }
        
        plot(childplot)
        
        # make a plot per admin 1
        for(a1 in unique(a1_gd$ADM1_CODE)){
          
          ad1name <- a1_gd[ADM1_CODE==a1, ADM1_NAME] %>% unique
          
          # check for NA's
          check <- mbg[ADM1_NAME == ad1name]
          
          if (check$value %>% mean %>% is.na %>% all) {
            
            write(paste0(ad1name, ' in ', adname, '; '), 
                  file = paste0(outputdir, 'diagnostic_plots/ad1_missing_estimates.txt'),
                  append = TRUE)
            
          } else {

            # make an admin 1 plot of the stackers and whatnot
            maxval <- max(c(a1_gd[ADM1_CODE==a1,upper_rate], adm0[ADM0_CODE == admin_id,value], ad1_dat[ADM1_CODE == a1,outcome]), na.rm = TRUE)
            
            a1_plot <- ggplot(data = ad1_dat[ADM1_CODE == a1,]) + 
              geom_line(data = a1_gd[ADM1_CODE==a1,], aes(x = year, y = value, color = var_name), size = 1.2) +
              geom_ribbon(data = a1_gd[ADM1_CODE==a1 & variable == 'Unraked',], 
                          aes(x = year, ymin = lower_rate, ymax = upper_rate), na.rm = TRUE, alpha = 0.2, fill = 'indianred3') + 
              ylim(0, maxval) +
              scale_color_brewer(palette = 'Set1') +
              ggtitle(paste(toupper(indicator), ' Model Results for:',ad1name,'|',adname, '|', a1)) +
              theme_bw()
            
            if (nrow(ad1_dat[ADM1_CODE == a1,])>0) {
              a1_plot <- a1_plot + 
                geom_point(aes(x = year, y = outcome, size = N)) +
                geom_text_repel(aes(x = year, y = outcome, label = N_label),
                                size=5, segment.colour = 'gray', segment.size = .2, direction = 'x', nudge_y=.05, force=5)
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
