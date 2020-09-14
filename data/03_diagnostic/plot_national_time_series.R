# Clear Environment
rm(list = ls())

# Load necessary libraries
package_list <- c('dplyr', 'ggrepel')

source('<<<< FILEPATH REDACTED >>>>')
mbg_setup(package_list = package_list, repos='<<<< FILEPATH REDACTED >>>>')

# Read in input data
mydat <- fread('<<<< FILEPATH REDACTED >>>>')

#shared funs
source("<<<< FILEPATH REDACTED >>>>/get_outputs.R")
source("<<<< FILEPATH REDACTED >>>>/get_location_metadata.R")

#location metadata
locs = get_location_metadata(location_set_id = 9, gbd_round_id = 5)

#gbd data for lrirrhea
gbd_dat <- get_outputs('cause',cause_id = 322, measure_id = 5, metric_id = 3, year_id = c(2000:2017), location_id = locs[,location_id], age_group_id = 1, sex_id = 3, gbd_round_id = 5, compare_version_id = 377)

# Calculate weighted means by country-year-nid-source as well as weighted sum of
# sample size post-processing
mydat3 <- mydat %>%
  group_by(nid, start_year, country, survey_series) %>%
  dplyr::summarize(prev = weighted.mean(x = has_lri/N, w = N, na.rm = T),
                   N = sum(N*weight, na.rm = T)) %>% subset(start_year > 1999) %>%
  mutate(survey_series = substr(survey_series, 1, 7))

#read in region list
regions <- fread('<<<< FILEPATH REDACTED >>>>',showProgress = FALSE)
regions <- regions[,c('iso3','loc_id','mbg_reg')]
regions <- regions[mbg_reg != '']
colnames(regions) <- c('iso3','loc_id','region')

region_fix <- list(
  'dia_afr_horn' = 'dji+eri+eth+sdn+som+ssd+yem',
  'dia_name' = 'dza+egy+esh+lby+mar+tun',
  'dia_sssa' = 'bwa+nam+zaf',
  'dia_mcaca' = 'blz+cri+cub+dma+dom+grd+gtm+hnd+hti+jam+lca+mex+nic+pan+slv+vct',
  'dia_central_asia' = 'kgz+tjk+tkm+uzb',
  'dia_se_asia' = 'khm+lao+mmr+mys+tha+vnm',
  'dia_malay' = 'idn+phl+png+tls',
  'dia_mid_east' = 'afg+irn+irq+jor+pse+syr',
  'dia_cssa' = 'ago+caf+cod+cog+gab+gnq+stp',
  'dia_wssa' = 'ben+bfa+civ+cmr+cpv+gha+gin+gmb+gnb+lbr+mli+mrt+ner+nga+sen+sle+tcd+tgo',
  'dia_s_america' = 'bol+bra+col+ecu+guf+guy+per+pry+sur+tto+ven',
  'dia_chn_mng' = 'chn+mng',
  'dia_south_asia' = 'bgd+btn+ind+lka+npl+pak',
  'dia_essa' = 'bdi+com+ken+lso+mdg+moz+mwi+rwa+swz+syc+tza+uga+zmb+zwe'
)

for (reg in 1:length(region_fix)){
  countries <- region_fix[[reg]]
  countries <- unlist(strsplit(countries, '+', fixed = T))
  countries <- toupper(countries)
  for(c in countries){
    regions[iso3 == c, region := names(region_fix)[[reg]]]
  }
}

regions <- subset(regions, region %like% 'dia')

#create plots by region
region_list <- unique(regions$region)
today <- Sys.Date()
today <- gsub("-", "_", today)
dir.create(paste0('<<<< FILEPATH REDACTED >>>>', today))
for (reg in region_list){
  temp_locs <- as.data.table(filter(regions, region == reg))
  max_lri <- subset(mydat3, country %in% temp_locs$iso3)
  if(nrow(max_lri) != 0){
    message('____________')
    message(reg)
    max_lri <- max(max_lri$prev, na.rm = TRUE)
    pdf('<<<< FILEPATH REDACTED >>>>')
    for (j in temp_locs$iso3) {
      message(j)
      test <- filter(mydat3, country == j)
      id <- unique(temp_locs[iso3 == j,]$loc_id)
      testgbd <- filter(gbd_dat, location_id == id)
      names(testgbd)[names(testgbd)=='year_id'] <- 'start_year'
      names(testgbd)[names(testgbd)=='val'] <- 'prev'
      max_lri <- max(max_lri, max(testgbd$upper))
      if (nrow(test) > 0) {
        print(
          p1 <-  ggplot(test, aes(x = start_year, y = prev)) +
            geom_ribbon(aes(ymin = lower, ymax = upper, fill = 'GBD Estimates 95% CI'), data = testgbd, alpha = 0.3) +
            scale_fill_manual("",values=alpha("grey12", .2)) +
            geom_line(aes(x = start_year, y = prev), data = testgbd, alpha = .4) +
            geom_point(aes(x = start_year, y = prev, size = log(N),
                           color = survey_series), alpha = 0.4) +
            geom_text_repel(aes(x = start_year, y = prev, label = nid), size = 2) +
            xlim(2000, 2018) + ylim(0,.03) +
            ggtitle(paste0(j, ' ')) +
            theme_bw() +
            labs(size = 'Sample Size (Logged)',color = 'Survey Series')
        )
      }
    }
    dev.off()
  }
}
