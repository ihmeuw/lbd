#custom function to interrogate covariates over time

covInterrogation <- function(scaled=centre_scale_covs, pixel_count=5) {

  #improve the covariate object naming
  covariate_list <- copy(all_cov_layers)
  tv_cov_names = c()
  for(lll in names(covariate_list)){
    if(class(covariate_list[[lll]]) == 'RasterBrick'){
      tv_cov_names = append(tv_cov_names, lll)
      if (nrow(period_map) > 1)  names(covariate_list[[lll]]) = paste(lll, year_list %>% as.character, sep='.') #use XXX as a place holder
    }
    
  }
  
  #custom function to process the cov rasters in data.tables
  processCovs <- function(i, centre_scale=F) {
    
    if (dim(i)[3] > 1) {
      
      dt <- as.data.frame(i, xy = TRUE) %>% 
        as.data.table
      
      #centre_scale the results
      if (centre_scale) {
        
        design_matrix = data.frame(dt[, -c('x', 'y'), with =F])
        cs_df <- getCentreScale(design_matrix)
        design_matrix <- centreScale(design_matrix, df = cs_df)
        
        #replace the df columns with the design matrix
        dt[, names(dt)[!grepl('x|y', names(dt))] := NULL]
        dt = cbind(dt, design_matrix)
        
      }
      
      melt(dt, id.vars = c('x','y')) %>% 
        .[, c("variable", "year") := tstrsplit(variable, ".", fixed=TRUE)] %>% 
        return
      
    } else return(NULL)
    
  }
  
  #process all covs
  dt <- lapply(covariate_list, processCovs, centre_scale=F) %>% rbindlist
  if(scaled) dt_scaled <- lapply(covariate_list, processCovs, centre_scale=T) %>% rbindlist
  
  #custom function to produce timeplots for n random pixels
  timeplotCovs <- function(dt, scaled=F, n=pixel_count) {
    
    title <- 'Covariate Interrogation Plots'
    title <- paste0(title, ifelse(!scaled, '', '(Centre-scaled)'))
    
    sample <- dt[!is.na(value), .(x,y)][sample(.N, n)] %>% 
      dt[., on=.(x,y)] %>%
      .[, id := paste0(x,y)]
    
    plot <- ggplot(sample, aes(x=year %>% as.numeric, y=value, color=paste0(x, y))) + 
      geom_point() + 
      geom_smooth() + 
      facet_wrap(~variable) + 
      labs(title=title, subtitle='5 Random Pixels') +
      theme_minimal()
    
    print(plot)
    
    plot <- ggplot(sample, aes(x=year %>% as.numeric, y=value, color=paste0(x, y))) + 
      geom_point() + 
      geom_smooth() + 
      facet_wrap(~variable, scales='free') + 
      labs(title=title, subtitle='5 Random Pixels') +
      theme_minimal()
    
    print(plot)
    
    #if worldpop included, make a plot using the top 5 pixels
    if('worldpop' %in% unique(dt$variable)) {
      
      urb_sample <- copy(dt[year==max(year) & variable=='worldpop' & !is.na(value)]) %>% 
        setkey(., value) %>% #sort by pop
        tail(n) %>% #take the top n
        .[, .(x,y)] %>% 
        dt[., on=.(x,y)] %>%
        .[, id := paste0(x,y)]
      
      plot <- ggplot(urb_sample, aes(x=year %>% as.numeric, y=value, color=paste0(x, y))) + 
        geom_point() + 
        geom_smooth(alpha=.2) + 
        facet_wrap(~variable) + 
        labs(title=title, subtitle='Top 5 Pixels by Population') +
        theme_minimal()
      
      print(plot)
      
      plot <- ggplot(urb_sample, aes(x=year %>% as.numeric, y=value, color=paste0(x, y))) + 
        geom_point() + 
        geom_smooth(alpha=.2) + 
        facet_wrap(~variable, scales='free') + 
        labs(title=title, subtitle='Top 5 Pixels by Population') +
        theme_minimal()
      
      print(plot)
      
    }
    
    return(NULL)
    
  }
  
  #make sure the folder has been made
  covdir <- file.path(outputdir, 'covariates') %T>% 
    dir.create
  
  pdf(paste0(covdir, '/', reg, '_', 'timeplots.pdf'), height=8, width=12)
  timeplotCovs(dt, scaled=F)
  if(scaled) timeplotCovs(dt_scaled, scaled=T)
  dev.off()
  
}