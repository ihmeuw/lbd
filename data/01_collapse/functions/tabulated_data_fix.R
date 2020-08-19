#load packages
library(data.table)
library(dplyr)

# Read in data
mydat <- fread('<<<< FILEPATH REDACTED >>>>')

sb <- filter(mydat, survey_series %in% c('SWACHHTA_REPORT', 'NARSS')) %>%
  dplyr::select(survey_series, shapefile, location_code, imp_cw, iso3, 
                nid, lat, long, total_hh, 
                year = year_start)

# Read in definitions
setwd('<<<< FILEPATH REDACTED >>>>')

sani_cw <- fread('<<<< FILEPATH REDACTED >>>>')
ind_cw <- filter(sani_cw, iso3 == 'IND')

cw_ratios <- ind_cw %>%
  mutate(total = s_piped + imp + unimp + od + latrine_imp +
           latrine_unimp + latrine_cw + flush_unimp + 
           flush_cw*(flush_unimp/(flush_unimp+flush_imp))) %>%
  mutate(
    s_network = network,
    s_piped = s_piped,
    s_imp_other = imp + latrine_imp + latrine_cw*
      (latrine_imp/(latrine_imp+latrine_unimp)),
    s_unimp = unimp + latrine_cw*(latrine_unimp/(latrine_unimp+latrine_unimp)) + 
      latrine_unimp +flush_cw*(flush_unimp/(flush_unimp+flush_imp)) + flush_unimp,
    s_od = od
  ) %>%
  select(total, s_network, s_piped, s_imp_other, s_unimp, s_od) %>%
  mutate(
    s_network = s_network/s_piped,
    s_piped = s_piped/total,
    s_imp_other = s_imp_other/total,
    s_unimp = s_unimp/total,
    s_od = s_od/total,
    total = total/total) 
sb2 <- mutate(sb,
              s_piped = imp_cw *(cw_ratios$s_piped/(cw_ratios$s_piped+cw_ratios$s_imp_other)),
              s_network = s_piped * cw_ratios$s_network,
              s_imp_other = imp_cw *(cw_ratios$s_imp_other/(cw_ratios$s_piped+cw_ratios$s_imp_other)),
              s_unimp = (1-imp_cw)*(cw_ratios$s_unimp/(cw_ratios$s_unimp + cw_ratios$s_od)),
              s_od = (1-imp_cw)*(cw_ratios$s_od/(cw_ratios$s_unimp + cw_ratios$s_od))) 

sb2$total_hh <- ifelse(is.na(as.numeric(sb2$total_hh)), 0, as.numeric(sb2$total_hh))

sb2 %>%
  summarize(
    s_network = weighted.mean(x = s_network, w = total_hh),
    s_piped = weighted.mean(x = s_piped, w = total_hh),
    s_imp_other = weighted.mean(x = s_imp_other, w = total_hh),
    s_unimp = weighted.mean(x = s_unimp, w = total_hh),
    s_od = weighted.mean(x = s_od, w = total_hh)
  )

sb3 <- sb2 %>%
  mutate(
    s_network_cr = s_network,
    s_imp_cr = s_imp_other/(1 - s_piped),
    s_unimp_cr = s_unimp/(1 - s_piped - s_imp_other)
  ) %>%
  select(shapefile, location_code,
         iso3, nid, total_hh, year,
         s_piped, s_imp_cr, s_unimp_cr, s_network_cr)