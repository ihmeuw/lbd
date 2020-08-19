# Clear Environment
rm(list = ls())

# Load necessary libraries
package_list <- c('dplyr','ggrepel',
                  'openxlsx','data.table')

source('<<<< FILEPATH REDACTED >>>>')
mbg_setup(package_list = package_list, 
          repos='<<<< FILEPATH REDACTED >>>>')

regions <- read.csv('<<<< FILEPATH REDACTED >>>>')
regions <- regions[,c('iso3', 'region')]

codebook <- fread('<<<< FILEPATH REDACTED >>>>')

#SANI JMP NUMBERS from https://washdata.org/data/downloads#WLD
jmp_sani <- read.xlsx('<<<< FILEPATH REDACTED >>>>', sheet = 3)
iso3 <- jmp_sani$X2[4:3715]
year <- as.numeric(jmp_sani$X3[4:3715])
jmp_sani_basic <- as.numeric(jmp_sani$X6[4:3715])
jmp_sani_limited <- as.numeric(jmp_sani$X7[4:3715])
imp <- jmp_sani_basic + jmp_sani_limited
unimp <- as.numeric(jmp_sani$X8[4:3715])
jmp_sani <- data.table(iso3, year, imp, unimp)
jmp_sani <- melt(jmp_sani, id = c('iso3','year'), 
                 value.name = "jmp_prop", 
                 variable.name = 'indicator')

#WATER JMP NUMBERS from https://washdata.org/data/downloads#WLD
jmp_water <- read.xlsx('<<<< FILEPATH REDACTED >>>>', sheet = 2)
iso3 <- jmp_water$X2[4:3715]
year <- as.numeric(jmp_water$X3[4:3715])
jmp_water_basic <- as.numeric(jmp_water$X6[4:3715])
jmp_water_limited <- as.numeric(jmp_water$X7[4:3715])
imp <- jmp_water_basic + jmp_water_limited
unimp <- as.numeric(jmp_water$X8[4:3715])
jmp_water <- data.table(iso3, year, imp, unimp)
jmp_water <- melt(jmp_water, id = c('iso3','year'), 
                  value.name = "jmp_prop", 
                  variable.name = 'indicator')

#prepare data
for(indi_group in c('water','sani')){
  if(indi_group == 'water'){
    fam <- 'w'
  }else{
    fam <- 's'
  }
  for(indi in c('imp_cr','unimp_cr')){
    mydat <- fread(paste0('<<<< FILEPATH REDACTED >>>>',
                          fam, '_', indi, '.csv')) %>%
      group_by(nid, surv_year, country) %>%
      summarize(mean = weighted.mean(x = prop, 
                                     w = sum_of_sample_weights*weight),
                N = sum(weight*N)) %>%
      dplyr:::select(country, nid, surv_year, mean, N) %>%
      rename(prop = mean) %>% setDT()
    
    if (indi == 'imp_cr'){
      pipedat <- fread(paste0('<<<< FILEPATH REDACTED >>>>',
                              fam, '_', 'piped', '.csv')) %>%
        group_by(nid, surv_year, country) %>%
        summarize(piped_prop = weighted.mean(x = prop, 
                                             w = sum_of_sample_weights*weight),
                  piped_N = sum(weight*N)) %>%
        dplyr:::select(country, nid, surv_year, piped_prop, piped_N)
      
      mydat <- left_join(pipedat, mydat)
      mydat[is.na(mydat)] <- 0
      
      mydat <- as.data.table(mydat)
      mydat <- mydat[,prop := piped_prop + ((1 - piped_prop)*(prop))]
      mydat <- mydat[,N := N + piped_N]
      mydat <- mydat[,c('country', 'nid', 'surv_year', 'prop', 'N')]
    }
    
    mydat[,subnat := ifelse(nid %in% codebook$nid, 1, 0)]
    assign(paste0(indi_group, '_', indi), mydat)
  }
}

#make plots
today <- Sys.Date()
dir.create(file.path(paste0('<<<< FILEPATH REDACTED >>>>', today)), 
           showWarnings = FALSE)
for(indi_group in c('water','sani')){
  if (indi == 'water'){
    jmp <- copy(jmp_water)
  }else{
    jmp <- copy(jmp_sani)
  }
  
  imp <- get(paste0(indi_group, '_', 'imp_cr')) %>% 
    mutate(indicator = 'imp')
  unimp <- get(paste0(indi_group, '_', 'unimp_cr')) %>% 
    mutate(indicator = 'unimp')
  
  imp <- merge(imp, regions, by.x = 'country', 
               by.y = 'iso3', all.x = T)
  imp <- subset(imp, region %in% names(region_fix)) %>% 
    arrange(region, country, nid)
  unimp <- merge(unimp, regions, 
                 by.x = 'country', by.y = 'iso3', all.x = T)
  unimp <- subset(unimp, region %in% names(region_fix)) %>% 
    arrange(region, country, nid)
  unimp$prop <- (1 - imp$prop)*unimp$prop
  
  data <- rbind(imp, unimp)
  
  pdf(paste0('<<<< FILEPATH REDACTED >>>>', 
             today,'/',indi_group,'_national_ts.pdf'),15, 8.5)
  for(c in unique(imp$country)){
    jmp_temp <- jmp[iso3 == c,]
    wash_plot <-  ggplot(data = subset(data, country == c), 
                         aes(x = surv_year, y = prop)) +
      scale_x_continuous(breaks = seq(2000, 2018, 2), 
                         limits = c(2000, 2018)) +
      geom_point(aes(x = surv_year, y = prop, size = log(N), 
                     color = as.factor(subnat))) +
      geom_text_repel(aes(x = surv_year, y = prop, label = nid), 
                      size = 2) +
      facet_wrap(~indicator) +
      geom_line(aes(x = year, y = jmp_prop/100), 
                data = jmp_temp) +
      #xlim(2000, 2018) +
      ylim(0, 1) +
      ggtitle(paste0('Comparing ', indi_group, 
                     ' to JMP in ', c)) +
      theme_bw() +
      labs(color = 'Spatial Representation', 
           size = 'Sample Size (logged)') +
      #theme_gray() +
      theme(legend.position = "bottom")
    plot(wash_plot)
  }
  dev.off()
}
