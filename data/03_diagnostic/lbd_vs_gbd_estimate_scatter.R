#####################################################################
####    Generate scatter plot of LBD LRI data vs GBD Estimates   ####
#####################################################################



# Clear Environment
rm(list = ls())

# Load necessary libraries
package_list <- c('dplyr', 'ggrepel', 'data.table')
source('<<<< FILEPATH REDACTED >>>>')
mbg_setup(package_list = package_list, repos='<<<< FILEPATH REDACTED >>>>')

# Read in LBD input data
mydat <- fread('<<<< FILEPATH REDACTED >>>>')

# Read in GBD Estimates

source("<<<< FILEPATH REDACTED >>>>/get_outputs.R")
source("<<<< FILEPATH REDACTED >>>>/get_location_metadata.R")

locs = get_location_metadata(location_set_id = 9, gbd_round_id = 5)

gbd_dat <- get_outputs('cause',cause_id = 322, measure_id = 5, metric_id = 3, 
                       year_id = c(2000:2016), location_id = locs[,location_id], 
                       age_group_id = 1, sex_id = 3, gbd_round_id = 5)

gbd_iso3 <- merge(gbd_dat, locs, all.x = TRUE,all.y = FALSE, by = 'location_name')
gbd_iso3 <- gbd_iso3[,c('ihme_loc_id', 'val', 'year_id')]
gbd_iso3 <- gbd_iso3[complete.cases(gbd_iso3),]

# Subset data to modeling period
mydat <- filter(mydat, year >= 2000)

# Merge data
all_dat <- as.data.table(merge(mydat, gbd_iso3, by.x = c('country', 'year'), by.y = c('ihme_loc_id', 'year_id'), allow.cartesian = TRUE, all.x = TRUE, all.y = FALSE))
# Set years into 5 year bins for easier to view plotting
all_dat <- all_dat[, year := ifelse(year < 2005, 2000, ifelse(year < 2010, 2005, ifelse(year < 2015, 2010, 2015)))]

mydat3 <- all_dat %>%
  group_by(nid, year, point, country, survey_series) %>%
  summarize(prev = weighted.mean(x = has_lri/N, w = N*weight, na.rm = T),
            N = sum(weight*N, na.rm = T), val = mean(val))




pdf('<<<< FILEPATH REDACTED >>>>')
print(
  ggplot(mydat3, aes(x = prev, y = val, color = as.factor(year))) +
    geom_point() + scale_colour_brewer(palette = "Set1") + ggtitle('LBD Data vs GBD Estimates') +
    geom_abline(intercept = 0, slope = 1, color="red", linetype="solid", size=0.5) +
    scale_x_continuous(limits=c(0,0.25)) +
    scale_y_continuous(limits=c(0,0.25)) +
    labs(x = 'LBD Data', y = 'GBD Estimates')
)
dev.off()


##########################################################################################
# Per country scatters

regions <- c('sssa_hi', 'cssa', 'name_hi', 'essa_hilo', 'wssa', 'latin_america', 's_asia', 'se_asia', 'middle_east')

for (region in regions){
  message(region)
  pdf('<<<< FILEPATH REDACTED >>>>')
  region_obj <- get(region)
  for (i in region_obj){
    test <- filter(mydat3, country == i)
    if (nrow(test) > 0){
      print(
        ggplot(test, aes(x = prev, y = val, color = as.factor(year))) +
          geom_point(aes(x = prev, y = val, size = N, shape = as.factor(point))) +
          scale_colour_brewer(palette = "Set1") + ggtitle(paste0('LBD Data vs GBD Estimates in ', i)) +
          geom_abline(intercept = 0, slope = 1, color="red", linetype="solid", size=1.5) +
          geom_text_repel(aes(x = prev, y = val, label = nid), color = "black", size = 3) +
          scale_x_continuous(limits=c(0,max(max(test$val),max(test$prev)) + .01)) +
          scale_y_continuous(limits=c(0,max(max(test$val),max(test$prev)) + .01)) +
          labs(x = 'LBD Data', y = 'GBD Estimates')
      )
    }
  }
  dev.off()
}

