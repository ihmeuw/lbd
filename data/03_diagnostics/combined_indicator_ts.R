rm(list = ls())
# load libraries
libs <- c('ggplot2', 'data.table', 'dplyr', 'raster', 'tools')
lapply(libs, library, character.only = TRUE)

# define indicator groups
#indi_group <- 'water'
for(indi_group in c('water','sani')){
  if (indi_group == 'water') {
    indi <- c('w_network_cr','w_piped', 'w_imp_cr', 'w_unimp_cr')
    network_string <- 'piped on premises'
    piped_string <- 'piped to public'
  } else {
    indi <- c('s_network_cr','s_piped', 's_imp_cr', 's_unimp_cr')
    network_string <- 'sewer'
    piped_string <- 'septic'
  }
  
  inputdir <- file.path('<<<< FILEPATH REDACTED >>>>')
  setwd(inputdir)
  
  indicator <- indi[1]
  id <- list.files(pattern = paste0(indi[1], '.csv'))
  input_data <- fread(id)
  network <- input_data %>%
    group_by(nid, surv_year, country) %>%
    summarize(mean = weighted.mean(x = prop, w = sum_of_sample_weights*weight),
              ss = sum(weight*N)) %>%
    dplyr:::select(country, nid, surv_year, mean) %>%
    rename(network = mean)
  
  indicator <- indi[2]
  id <- list.files(pattern = paste0(indi[2], '.csv'))
  input_data <- fread(id)
  piped <- input_data %>%
    group_by(nid, surv_year, country) %>%
    summarize(mean = weighted.mean(x = prop, w = sum_of_sample_weights*weight),
              ss = sum(weight*N)) %>%
    dplyr:::select(country, nid, surv_year, mean) %>%
    rename(piped = mean)
  
  indicator <- indi[3]
  id <- list.files(pattern = paste0(indi[3], '.csv'))
  input_data <- fread(id)
  imp <- input_data %>%
    group_by(nid, surv_year, country) %>%
    summarize(mean = weighted.mean(x = prop, w = sum_of_sample_weights*weight),
              ss = sum(weight*N)) %>%
    dplyr:::select(country, nid, surv_year, mean) %>%
    rename(imp = mean)
  
  indicator <- indi[4]
  id <- list.files(pattern = paste0(indi[4], '.csv'))
  input_data <- fread(id)
  unimp <- input_data %>%
    group_by(nid, surv_year, country) %>%
    summarize(mean = weighted.mean(x = prop, w = sum_of_sample_weights*weight),
              ss = sum(weight*N)) %>%
    dplyr:::select(country, nid, surv_year, mean) %>%
    rename(unimp = mean)
  
  test <- subset(piped, !nid %in% network$nid) %>%
    rename(network = piped) %>%
    mutate(network = 0)
  network <- rbind(network, test)
  
  mydat <- left_join(network, piped)
  mydat <- left_join(mydat, imp)
  mydat <- left_join(mydat, unimp)
  mydat$imp <- ifelse(is.na(mydat$imp), 0, mydat$imp)
  mydat$unimp <- ifelse(is.na(mydat$unimp), 0, mydat$unimp)
  mydat$surface <- 0
  mydat <- mydat %>%
    mutate(piped_on_premises = piped*network) %>%
    mutate(piped_public = (1 - network)*piped) %>%
    mutate(imp = (1 - piped)*imp) %>%
    mutate(unimp = (1 - piped - imp)*unimp) %>%
    mutate(surface = 1 - piped - imp - unimp) %>%
    mutate(piped_on_premises_shape = ifelse(piped_on_premises == 0, 'triangle', 'dot'),
           piped_public_shape = ifelse(piped_public == 0, 'triangle', 'dot'),
           imp_shape = ifelse(imp == 0, 'triangle', 'dot'),
           unimp_shape = ifelse(unimp == 0, 'triangle', 'dot'),
           surface_shape = ifelse(surface == 0, 'triangle', 'dot'))
  

  regions <- fread('<<<< FILEPATH REDACTED >>>>',showProgress = FALSE)
  regions <- regions[,c('iso3','loc_id','mbg_reg')]
  colnames(regions) <- c('iso3','loc_id','region')

  mydat <- merge(mydat, regions, by.x = 'country', by.y = 'iso3', all.x = T)
  mydat <- subset(mydat, region %in% names(region_fix)) %>% arrange(region, country)
  
  library(ggplot2)
  library(ggrepel)
  date <- gsub("-", "_", Sys.Date())
  
  dir.create(paste0('<<<< FILEPATH REDACTED >>>>', date), showWarnings = FALSE)
  pdf(paste0('<<<< FILEPATH REDACTED >>>>',date,'/', indi_group,'.pdf'),
      8.5, 8.5)
  for (i in unique(mydat$country)) {
    print(i)
    print(
      ggplot(data = filter(mydat, country == i)) +
        scale_x_continuous(breaks = seq(2000, max(mydat$surv_year), 2), 
                           limits = c(2000, max(mydat$surv_year))) +
        geom_point(aes(x= surv_year, y = piped_on_premises, 
                       col = 'piped_on_premises', 
                       shape = as.character(piped_on_premises_shape)), 
                   size =3, alpha = .55) +
        geom_point(aes(x= surv_year, y = piped_public, 
                       col = 'piped_public', 
                       shape = as.character(piped_public_shape)), 
                   size =3, alpha = .55) +
        geom_point(aes(x= surv_year, y = imp, 
                       col = 'imp', 
                       shape = as.character(imp_shape)), 
                   size =3, alpha = .75) +
        geom_point(aes(x= surv_year, y = unimp, 
                       col = 'unimp', 
                       shape = as.character(unimp_shape)), 
                   size =3, alpha = .55) +
        geom_point(aes(x= surv_year, y = surface, 
                       col = 'surface', 
                       shape = as.character(surface_shape)), 
                   size =3, alpha = .55) +
        
        geom_line(aes(x= surv_year, y = piped_on_premises, col = 'piped_on_premises'), 
                  linetype = 'dashed', alpha = .65) +
        geom_line(aes(x= surv_year, y = piped_public, col = 'piped_public')
                  , linetype = 'dashed', alpha = .65) +
        geom_line(aes(x= surv_year, y = imp, col = 'imp')
                  , linetype = 'dashed', alpha = .85) +
        geom_line(aes(x= surv_year, y = unimp, col = 'unimp')
                  , linetype = 'dashed', alpha = .65) +
        geom_line(aes(x= surv_year, y = surface, col = 'surface')
                  , linetype = 'dashed', alpha = .65) +
        
        geom_text_repel(aes(x= surv_year, y = piped_on_premises, 
                            label = nid, col = 'piped_on_premises'), 
                        size =3, alpha = .65) +
        geom_text_repel(aes(x= surv_year, y = piped_public, 
                            label = nid, col = 'piped_public'), 
                        size =3, alpha = .65) +
        geom_text_repel(aes(x= surv_year, y = imp, 
                            label = nid, col = 'imp'), 
                        size =3, alpha = .85) +
        geom_text_repel(aes(x= surv_year, y = unimp, 
                            label = nid, col = 'unimp'), 
                        size =3, alpha = .65) +
        geom_text_repel(aes(x= surv_year, y = surface, 
                            label = nid, col = 'surface'), 
                        size =3, alpha = .65) +
        
        labs(color = 'Indicator', x = 'Year', y = 'Prevalence') +
        
        scale_color_manual(breaks = c('piped_on_premises','piped_public',
                                      'imp','unimp','surface'),
                           labels = c(network_string,piped_string,'imp',
                                      'unimp','surface'),
                           values = c('#fa9fb5','#e41a1c', '#377eb8', 
                                      '#4daf4a', '#984ea3')) +
        ylim(0, 1) +
        theme_gray() +
        theme(legend.position = "top") +
        guides(shape=FALSE) +
        ggtitle(paste0(toTitleCase(indi_group), 
                       ' Access in: ', i,'\n','Region: ', 
                       unique(filter(mydat, country == i)$region)))
    )
  }
  dev.off()
}

