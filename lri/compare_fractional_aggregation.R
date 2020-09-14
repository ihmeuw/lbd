compare_frx_agg <- function(reg,
                            admin,
                            raked,
                            measure){
  
  if (raked) raked_ind <- paste0('raked_', measure) else raked_ind <- 'unraked'
  if (raked) raked_ind_leg <- 'raked' else raked_ind_leg <- 'unraked'
  
  #frx raked file
  load('<<<< FILEPATH REDACTED >>>>')
  frx_rake <- get(paste0('admin_', admin)) 
  
  #legacy raked file
  load('<<<< FILEPATH REDACTED >>>>')
  leg_rake <- get(paste0('admin_', admin))
  
  draw_cols <- grep('V', names(leg_rake), value = TRUE)
  frx_rake <- frx_rake[, mean := rowMeans(select(frx_rake, draw_cols))] %>%
    rename(frx_mean = mean)
  leg_rake <- leg_rake[, mean := rowMeans(select(leg_rake, draw_cols))] %>%
    rename(leg_mean = mean)
  
  compare <- merge(select(frx_rake, c('frx_mean', 'year', paste0('ADM', admin, '_CODE'))), select(leg_rake, c('leg_mean', 'year', paste0('ADM', admin, '_CODE'))), 
                   by = c( 'year', paste0('ADM', admin, '_CODE')))
  
  plot <- ggplot(compare, aes(x = frx_mean, y = leg_mean)) + geom_point(aes(color = year)) + ggtitle(paste0(reg, ' admin', admin, ' ', measure, ' ', raked_ind)) + 
    geom_abline(intercept = 0, slope = 1)
  plot(plot)
}

region_list <- c('cssa','essa','sssa', 'name', 'wssa')
admins <- c(0:2)
measures <- c('incidence', 'prevalence', 'mortality')
pdf('<<<< FILEPATH REDACTED >>>>')

#unraked estimates
raked <- FALSE
for (reg in region_list){
  for (admin in admins){
    compare_frx_agg(raked = raked,
                    reg = reg,
                    admin = admin,
                    measure = 'prevalence')
  }
}

#raked estimates
raked <- TRUE
for (reg in region_list){
  for (admin in admins){
    for (measure in measures){
      compare_frx_agg(raked = raked,
                      reg = reg,
                      admin = admin,
                      measure = measure)
    }
  }
}
dev.off()


