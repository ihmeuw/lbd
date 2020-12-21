# --------------------------------------------------------------------------------
# Script to run full analysis of ORT changes in Sierra Leone, Mali, and Senegal
# --------------------------------------------------------------------------------


# ------------------------------------------------
# Set up

# clear space
rm(list=ls())

# set repo
repo <- '<<<< FILEPATH REDACTED >>>>'

# set model arguments
run_date          <- '2019_12_13_12_00_46'
shapefile_version <- '2019_09_10'

# set countries
countries <- c('Sierra Leone', 'Mali', 'Senegal')
# ------------------------------------------------


# --------------------------------------------------------------------------------
# Packages and functions

# load required packages
libs <- list('raster', 'sf', 'sp', 'dplyr', 'fasterize', 'ggplot2', 'scales',
             'RColorBrewer', 'rgeos', 'rgdal', 'data.table', 'viridis', 'ggsn')
lapply(libs, library, character.only = TRUE)
library(scico, lib.loc = '<<<< FILEPATH REDACTED >>>>')

# load custom functions
funcs <- paste0(repo, 'post_estimation/',
                list('squeeze_outputs', 'area_plots', 'anchor_maps',
                     'treatment_change_analysis', 'efficacy_analysis'),
                '.R')
lapply(funcs, source)
# --------------------------------------------------------------------------------


# ------------------------------------------------------
# Squeeze outputs

squeeze_outputs(run_date = run_date,
                shapefile_version = shapefile_version)
# ------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------
# Load analysis inputs

# load config
config <- fread(paste0(repo, 'post_estimation/key_dates.csv'))

# load admin 2 shapefile
ad2_shp <- shapefile(paste0('<<<< FILEPATH REDACTED >>>>', 
                            shapefile_version, '/lbd_standard_admin_2.shp'))

# load any ORS coverage estimates
ors <- fread(paste0('<<<< FILEPATH REDACTED >>>>/any_ors_admin_2_squeezed_summary.csv'))

# load only RHF coverage estimates
rhf <- fread(paste0('<<<< FILEPATH REDACTED >>>>/rhf_only_admin_2_squeezed_summary.csv'))

# load no ORT coverage estimates
no_ort <- fread(paste0('<<<< FILEPATH REDACTED >>>>/no_ort_admin_2_squeezed_summary.csv'))

# create subnational groups
groups_master <- list(
  list(
    'N' = c('Tonkolili', 'Port Loko', 'Koinadugu', 'Kambia', 'Bombali'),
    'C' = c('Western Rural', 'Western Urban'),
    'S' = c('Bo', 'Bonthe', 'Moyamba', 'Pujehun', 'Kailahun', 'Kenema', 'Kono')
  ),
  groups <- list(
    'N' = unique(ors[ADM1_NAME %in% c('Timbuktu', 'Kidal', 'Gao', 'Mopti'), ADM2_NAME]),
    'S' = unique(ors[ADM1_NAME %in% c('Kayes', 'Koulikoro', 'Sikasso', 'Ségou'), ADM2_NAME]),
    'C' = unique(ors[ADM2_NAME %in% c('Bamako'), ADM2_NAME])
  ),
  groups <- list(
    'N' = unique(ors[ADM1_NAME %in% c('Saint-Louis', 'Louga', 'Matam'), ADM2_NAME]),
    'C' = unique(ors[ADM1_NAME %in% c('Thiès', 'Dakar', 'Kaolack', 'Fatick', 'Kaffrine', 'Thiès', 'Dakar', 'Diourbel', 'Tambacounda', 'Kédougou', 'Kolda'), ADM2_NAME]),
    'S' = unique(ors[ADM1_NAME %in% c('Sédhiou', 'Ziguinchor'), ADM2_NAME])
  )
)
names(groups_master) <- countries
# ------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------
# Make summary csv for Supplementary File 2

sf <- rbindlist(
  list(ors[, treatment := 'any_ORS'],
       rhf[, treatment := 'only_RHF'],
       no_ort[, treatment := 'no_ORT'])
)
sf <- sf[, c('ADM0_NAME', 'ADM1_NAME', 'ADM2_NAME', 'year', 'treatment', 'mean', 'upper', 'lower')]
setnames(sf, c('ADM0_NAME', 'ADM1_NAME', 'ADM2_NAME'), 
         c('country', 'first_admin_unit', 'second_admin_unit'))
write.csv(sf, '<<<< FILEPATH REDACTED >>>>/BMED_RTR_Supplementary_File_2.csv')
# ------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------
# Make line graphs for Sierra Leone for figure 2

pdf('<<<< FILEPATH REDACTED >>>>/sle_time_plots_rtr.pdf', width = 8.25, height = 6.375)

# plotting function
lineplot <- function(dat, group_name = '', ind_name,
                     group_colors = c('#58A6A6', '#EFA355', '#421E22'),
                     xlabels) {
  print(
    ggplot(dat) +
      geom_ribbon(aes(x = as.numeric(year), ymin = min_ad_S, ymax = max_ad_S), fill = group_colors[[3]], alpha = 0.55) +
      geom_line(aes(x = as.numeric(year), y = min_ad_S), size = 1, color = group_colors[[3]], alpha = 0.9) +
      geom_line(aes(x = as.numeric(year), y = max_ad_S), size = 1, color = group_colors[[3]], alpha = 0.9) +
      geom_ribbon(aes(x = as.numeric(year), ymin = min_ad_N, ymax = max_ad_N), fill = group_colors[[2]], alpha = 0.6) +
      geom_line(aes(x = as.numeric(year), y = min_ad_N), size = 1, color = group_colors[[2]], alpha = 0.8) +
      geom_line(aes(x = as.numeric(year), y = max_ad_N), size = 1, color = group_colors[[2]], alpha = 0.8) +
      geom_ribbon(aes(x = as.numeric(year), ymin = min_ad_C, ymax = max_ad_C), fill = group_colors[[1]], alpha = 0.45) +
      geom_line(aes(x = as.numeric(year), y = min_ad_C), size = 1, color = group_colors[[1]], alpha = 0.8) +
      geom_line(aes(x = as.numeric(year), y = max_ad_C), size = 1, color = group_colors[[1]], alpha = 0.8) +
      scale_x_continuous(name = '', breaks = c(1:4), labels = xlabels, expand = c(0.1, 0.1)) +
      scale_y_continuous(name = 'Coverage') +
      geom_vline(xintercept = c(1:4), linetype = 'dashed') +
      ggtitle(paste0('Changes in treatment with ', ind_name)) +
      theme_minimal(base_size = 22) +
      theme(axis.text.x=element_text(size=18, color = 'black'),
            axis.title.y=element_text(size=22, color = 'black'),
            axis.text.y=element_text(size=18, color = 'black'),
            title=element_text(size=22, color = 'black'),
            panel.grid.major = element_line(size = 0.6, colour = 'grey85'),
            panel.grid.minor = element_line(size = 0.3, colour = 'grey85')
      )
  )
}

# data prep function
prep_data <- function(dat) {
  dat <- dat[ADM0_NAME == c & year %in% y]
  dat[, year := as.factor(year)]
  dat[, mean := mean*100]
  levels(year) <-  y
  for (i in names(groups)) dat[ADM2_NAME %in% groups[[i]], group := i]
  dat[, min_ad := min(mean), by = c('group', 'year')]
  dat[, max_ad := max(mean), by = c('group', 'year')]
  setorderv(dat, 'group')
  dat <- unique(dat[, c('group', 'year', 'min_ad', 'max_ad')])
  dat <- dcast(data = dat, formula = year ~ group, value.var = c('min_ad', 'max_ad'))
  dt <- copy(dat)
  return(dat)
}

# map plotting function
plot_map <- function(map_colors = c('#421E22', '#EFA355', '#58A6A6')) {
    # subset shapefile
    shp_cty <- ad2_shp[which(ad2_shp@data$ADM0_NAME == c), ]
    shp_cty@data$ADM2_CODE <- as.numeric(shp_cty@data$ADM2_CODE)
    shp_cty <- st_as_sf(shp_cty)
    
    # add to data frame
    cty_groups <- copy(ors)
    for (i in names(groups)) cty_groups[ADM2_NAME %in% groups[[i]], group := i]
    map <- left_join(shp_cty, unique(cty_groups[, c('ADM2_CODE', 'group')]), by = 'ADM2_CODE')
    
    # plot
    print (
      ggplot(map) +
        geom_sf(aes(fill = group), color = 'black', show.legend = FALSE) +
        scale_fill_manual(values = map_colors) +
        theme_minimal() +
        coord_sf(datum = NA)
    )
}

# get year function
get_years <- function(c) {
  y <- config[name == c, grep('date', names(config)), with = FALSE]
  y <- c(y[,1][[1]], y[,2][[1]], y[,3][[1]], y[,4][[1]])
  return(y)
}

# get years
c <- 'Sierra Leone'
message(c)
y <- get_years(c)

# get groups
groups <- groups_master[[c]]

# plot ORS
dt <- prep_data(ors)
lineplot(dt, '', 'any ORS', group_colors = c('#421E22', '#EFA355', '#58A6A6'),
         xlabels = c('Start of\nstudy\n(2000)', 'Before \npolicy change\n(2009)', '3 years after \npolicy change\n(2013)', '2 years after\nEbola epidemic\n(2017)'))

# plot RHF
dt <- prep_data(rhf)
lineplot(dt, '', 'only RHF', group_colors = c('#421E22', '#EFA355', '#58A6A6'),
         xlabels = c('Start of\nstudy\n(2000)', 'Before \npolicy change\n(2009)', '3 years after \npolicy change\n(2013)', '2 years after\nEbola epidemic\n(2017)'))

# plot map
plot_map()

dev.off()
# ------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------
# Make line graphs for Mali for figure 3

pdf('<<<< FILEPATH REDACTED >>>>/mli_time_plots_rtr.pdf', width = 8.25, height = 6.375)

# get dates
c <- 'Mali'
message(c)
y <- get_years(c)

# get groups
groups <- groups_master[[c]]

# plot ORS
dt <- prep_data(ors)
lineplot(dt, '', 'any ORS', group_colors = c('#421E22', '#EFA355', '#58A6A6'),
         xlabels = c('Before\ninterventions\n(2001)', '1 year after\ninterventions\n(2004)', '8 years after\ninterventions\n(2011)', '6 years after\nwar in North\n(2018)'))

# plot RHF
dt <- prep_data(rhf)
lineplot(dt, '', 'only RHF', group_colors = c('#421E22', '#EFA355', '#58A6A6'),
         xlabels = c('Before\ninterventions\n(2001)', '1 year after\ninterventions\n(2004)', '8 years after\ninterventions\n(2011)', '6 years after\nwar in North\n(2018)'))

# plot map
plot_map()

dev.off()
# ------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------
# Make line graphs for Senegal for figure 4

pdf('<<<< FILEPATH REDACTED >>>>/sen_time_plots_rtr.pdf', width = 8.25, height = 6.375)

# plotting function
lineplot <- function(dat, group_name = '', ind_name,
                     group_colors = c('#58A6A6', '#EFA355', '#421E22'),
                     xlabels) {
  print(
    ggplot(dat) +
      geom_ribbon(aes(x = as.numeric(year), ymin = min_ad_C, ymax = max_ad_C), fill = group_colors[[1]], alpha = 0.55) +
      geom_line(aes(x = as.numeric(year), y = min_ad_C), size = 1, color = group_colors[[1]], alpha = 0.9) +
      geom_line(aes(x = as.numeric(year), y = max_ad_C), size = 1, color = group_colors[[1]], alpha = 0.9) +
      geom_ribbon(aes(x = as.numeric(year), ymin = min_ad_N, ymax = max_ad_N), fill = group_colors[[2]], alpha = 0.45) +
      geom_line(aes(x = as.numeric(year), y = min_ad_N), size = 1, color = group_colors[[2]], alpha = 0.8) +
      geom_line(aes(x = as.numeric(year), y = max_ad_N), size = 1, color = group_colors[[2]], alpha = 0.8) +
      geom_ribbon(aes(x = as.numeric(year), ymin = min_ad_S, ymax = max_ad_S), fill = group_colors[[3]], alpha = 0.6) +
      geom_line(aes(x = as.numeric(year), y = min_ad_S), size = 1, color = group_colors[[3]], alpha = 0.8) +
      geom_line(aes(x = as.numeric(year), y = max_ad_S), size = 1, color = group_colors[[3]], alpha = 0.8) +
      scale_x_continuous(name = '', breaks = c(1:4), labels = xlabels, expand = c(0.1, 0.1)) +
      scale_y_continuous(name = 'Coverage') +
      geom_vline(xintercept = c(1:4), linetype = 'dashed') +
      ggtitle(paste0('Changes in treatment with ', ind_name)) +
      theme_minimal(base_size = 22) +
      theme(axis.text.x=element_text(size=18, color = 'black'),
            axis.title.y=element_text(size=22, color = 'black'),
            axis.text.y=element_text(size=18, color = 'black'),
            title=element_text(size=22, color = 'black'),
            panel.grid.major = element_line(size = 0.6, colour = 'grey85'),
            panel.grid.minor = element_line(size = 0.3, colour = 'grey85')
            )
  )
}

# get dates
c <- 'Senegal'
message(c)
y <- get_years(c)

# get groups
groups <- groups_master[[c]]

# plot ORS
dt <- prep_data(ors)
lineplot(dt, '', ind_name = 'any ORS', group_colors = c('#EFA355', '#421E22', '#58A6A6'),
         xlabels = c('Start of\nstudy\n(2000)', 'Before\npolicy change\n(2006)', '4 years after \npolicy change\n(2012)', '2 years after\ninterventions\n(2017)'))

# plot RHF
dt <- prep_data(rhf)
lineplot(dt, '', 'only RHF', group_colors = c('#EFA355', '#421E22', '#58A6A6'),
         xlabels = c('Start of\nstudy\n(2000)', 'Before\npolicy change\n(2006)', '4 years after \npolicy change\n(2012)', '2 years after\ninterventions\n(2017)'))

# plot map
plot_map(map_colors = c('#EFA355', '#421E22', '#58A6A6'))

dev.off()
# ------------------------------------------------------------------------------------------


# -------------------------------------
# Make area plots for figure 1

make_area_plots(run_date = run_date)
# -------------------------------------


# -----------------------------------------------------------------------
# Treatment change analysis

# run for full study period
for (i in unique(config$iso3)) {
  message(i)
  analyze_treatment_changes(run_date = run_date,
                            year_start = config[iso3 == i, date_start],
                            year_end = config[iso3 == i, date_end],
                            country = i)
}

# run for up to pre intervention
for (i in unique(config$iso3)) {
  message(i)
  analyze_treatment_changes(run_date = run_date,
                            year_start = config[iso3 == i, date_start],
                            year_end = config[iso3 == i, date_pre],
                            country = i)
}

# run for up to post intervention
for (i in unique(config$iso3)) {
  message(i)
  analyze_treatment_changes(run_date = run_date,
                            year_start = config[iso3 == i, date_start],
                            year_end = config[iso3 == i, date_post],
                            country = i)
}

# run for during intervention
for (i in unique(config$iso3)) {
  message(i)
  analyze_treatment_changes(run_date = run_date,
                            year_start = config[iso3 == i, date_pre],
                            year_end = config[iso3 == i, date_post],
                            country = i)
}

# run for post intervention
for (i in unique(config$iso3)) {
  message(i)
  analyze_treatment_changes(run_date = run_date,
                            year_start = config[iso3 == i, date_post],
                            year_end = config[iso3 == i, date_end],
                            country = i)
}
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
# ORS:RHF efficacy ratio needed analysis

# run for full study period
for (i in unique(config$iso3)) {
  message(i)
  analyze_efficacy(run_date = run_date,
                   year_start = config[iso3 == i, date_start],
                   year_end = config[iso3 == i, date_end],
                   country = i)
}

# run up to pre intervention
for (i in unique(config$iso3)) {
  message(i)
  analyze_efficacy(run_date = run_date,
                   year_start = config[iso3 == i, date_start],
                   year_end = config[iso3 == i, date_pre],
                   country = i)
}

# run up to post intervention
for (i in unique(config$iso3)) {
  message(i)
  analyze_efficacy(run_date = run_date,
                   year_start = config[iso3 == i, date_start],
                   year_end = config[iso3 == i, date_post],
                   country = i)
}

# run for during intervention
for (i in unique(config$iso3)) {
  message(i)
  analyze_efficacy(run_date = run_date,
                   year_start = config[iso3 == i, date_pre],
                   year_end = config[iso3 == i, date_post],
                   country = i)
}

# run post intervention
for (i in unique(config$iso3)) {
  message(i)
  analyze_efficacy(run_date = run_date,
                   year_start = config[iso3 == i, date_post],
                   year_end = config[iso3 == i, date_end],
                   country = i)
}
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
# Make treatment change and efficacy maps

# set limits by country
ind_limits <- list(c(-0.04, 0.04), c(-0.04, 0.04), c(-0.04, 0.04))
names(ind_limits) <- unique(config$iso3)

# run for full study period
for (i in unique(config$iso3)) {
  message(i)
  make_anchor_maps(run_date = run_date,
                   year_start = config[iso3 == i, date_start],
                   year_end = config[iso3 == i, date_end],
                   country = i,
                   country_name = config[iso3 == i, name],
                   shp = ad2_shp,
                   change_limits = ind_limits[[i]])
}

# run up to pre intervention
for (i in unique(config$iso3)) {
  message(i)
  make_anchor_maps(run_date = run_date,
                   year_start = config[iso3 == i, date_start],
                   year_end = config[iso3 == i, date_pre],
                   country = i,
                   country_name = config[iso3 == i, name],
                   shp = ad2_shp,
                   change_limits = ind_limits[[i]])
}

# run up to post intervention
for (i in unique(config$iso3)) {
  message(i)
  make_anchor_maps(run_date = run_date,
                   year_start = config[iso3 == i, date_start],
                   year_end = config[iso3 == i, date_post],
                   country = i,
                   country_name = config[iso3 == i, name],
                   shp = ad2_shp,
                   change_limits = ind_limits[[i]])
}

# run for during intervention
for (i in unique(config$iso3)) {
  message(i)
  make_anchor_maps(run_date = run_date,
                   year_start = config[iso3 == i, date_pre],
                   year_end = config[iso3 == i, date_post],
                   country = i,
                   country_name = config[iso3 == i, name],
                   shp = ad2_shp,
                   change_limits = ind_limits[[i]])
}

# run post intervention
for (i in unique(config$iso3)) {
  message(i)
  make_anchor_maps(run_date = run_date,
                   year_start = config[iso3 == i, date_post],
                   year_end = config[iso3 == i, date_end],
                   country = i,
                   country_name = config[iso3 == i, name],
                   shp = ad2_shp,
                   change_limits = ind_limits[[i]])
}
# -----------------------------------------------------------------------


# --------------------------------------------------------------------------------------------
# Make summary tables and figures for admin 2 results of treatment change and RHF replacement

# directories
indis <- c('no_ort', 'rhf_only', 'any_ors')
dir1 <- paste0('<<<< FILEPATH REDACTED >>>>')
names(dir1) <- indis
dir2 <- paste0('<<<< FILEPATH REDACTED >>>>')

# changes in treatment data
treat <- list()
for (i in indis) {
  files <- list.files(dir1[[i]], pattern = paste0('change_summary_ad2'))
  # make sure we're just pulling the dates we want for Senegal
  dt <- lapply(paste0(dir1[[i]], files), fread)
  matches <- regmatches(files, gregexpr("[[:digit:]]+", files))
  for (j in 1:length(dt)) dt[[j]][, year_range := paste0(matches[[j]][1], '_', matches[[j]][2])]
  dt <- rbindlist(dt)
  treat[[i]] <- dt
}
names(treat) <- indis

# clean treatment change data
for (i in indis) treat[[i]][, indicator := i]
treat <- rbindlist(treat)
dt <- dcast(treat, ADM0_NAME + ADM2_CODE + ADM2_NAME + year_range ~ indicator, value.var = c('mean', 'upper', 'lower'))

# get mean, upper, and lower for replacement
for (i in c('mean', 'upper', 'lower')) {
  dt[, (paste0(i, '_rhf_replaced')) := get(paste0(i, '_no_ort')) < 0]
}

# add in draw summaries to create master table
files <- list.files(dir2, pattern = 'draw_summary')
sumz <- lapply(paste0(dir2, files), fread)
matches <- regmatches(files, gregexpr("[[:digit:]]+", files))
for (j in 1:length(sumz)) sumz[[j]][, year_range := paste0(matches[[j]][1], '_', matches[[j]][2])]
sumz <- rbindlist(sumz)
sumz <- sumz[agg_level == 'ADM2']
sumz[, c('V1', 'name', 'agg_level') := NULL]
dt <- merge(dt, sumz, by.x = c('ADM2_CODE', 'year_range'), by.y = c('code', 'year_range'))

# set study period by year range
setorderv(dt, c('ADM2_CODE', 'year_range'))
unique(dt[, c('ADM0_NAME', 'year_range')])
dt$period <- rep(c('pre', 'post', 'full', 'during', 'after'), times = nrow(dt)/5)

# clean up and save treatment change table
dt1 <- dt[, grep('ADM|period|mean|upper|lower', names(dt)), with = FALSE]
dt1 <- dt1[, !grep('replaced', names(dt1)), with = FALSE]
dt1 <- dt1[period != 'post' & period != 'full']
write.csv(dt1, paste0(dir1[['no_ort']], 'cleaned_treatment_change_summary_table.csv'))
  
# clean up and save RHF replacement table
dt2 <- dt[, grep('ADM|period|replaced', names(dt)), with = FALSE]
dt2 <- dt2[, !grep('num', names(dt2)), with = FALSE]
names(dt2) <- c('ADM2_CODE', 'ADM0_NAME', 'ADM2_NAME', 'mean', 'upper', 'lower', 'prop', 'period')
dt2 <- dt2[period != 'after' & period != 'during']
write.csv(dt2, paste0(dir1[['no_ort']], 'cleaned_rhf_replacement_summary_table_v2.csv'))

# clean up and prep data for figures 2-4
dt3 <- dt[(period %in% c('pre', 'post', 'full') & ADM0_NAME == 'Sierra Leone') |
            (period %in% c('pre', 'full') & (ADM0_NAME == 'Mali' | ADM0_NAME == 'Senegal'))]
dt3[, prev_state := data.table::shift(rhf_replaced_prop, type = 'lag')]
dt3[rhf_replaced_prop > 0.95, map_var := '2. RHF replaced']
dt3[period != 'pre' & prev_state > 0.95, map_var := '1. RHF previously replaced']
dt3[rhf_replaced_prop < 0.05, map_var := '3. RHF not replaced']
dt3[rhf_replaced_prop <= 0.95 & rhf_replaced_prop >= 0.05, map_var := '4. Uncertain']
write.csv(dt3, paste0(dir1[['no_ort']], 'summary_data_for_figs.csv'))
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Make map panels for manuscript figures 2-4

# set up color palette
map_colors = list('1. RHF previously replaced' = '#1b7837', # dark green 
                  '2. RHF replaced' = '#5aae61', # light green
                  '3. RHF not replaced' = '#9970ab', # dark purple
                  '4. Uncertain' = '#e7d4e8') #light purple

# add colors to data table
for (i in names(map_colors)) dt3[map_var == i, map_color := map_colors[[i]]]
setorderv(dt3, 'map_var')
dt3[, map_color := as.factor(map_color)]

# plot map
plot_map <- function(scale_dist = 50, 
                     scale_pos = 'bottomleft') {
    # subset shapefile
    shp_cty <- ad2_shp[which(ad2_shp@data$ADM0_NAME == c), ]
    shp_cty@data$ADM2_CODE <- as.numeric(shp_cty@data$ADM2_CODE)
    shp_cty <- st_as_sf(shp_cty)
    
    # add to data frame
    dt_map <- dt3[ADM0_NAME == c & period == p]
    map <- left_join(shp_cty, dt_map[, c('ADM2_CODE', 'map_var')], by = 'ADM2_CODE')
    
    # get mapping colors
    mapping_colors <- as.character(unique(dt_map$map_color))
    
    # plot
    print(
      ggplot(map) +
        geom_sf(aes(fill = map_var), color = 'black', show.legend = FALSE) +
        scale_fill_manual(values = mapping_colors) +
        theme_minimal() + labs(x = '', y = '') +
        coord_sf(datum = NA) +
        ggsn::scalebar(map, dist = scale_dist, st.size=4, height=0.03, st.dist = 0.03,
                       transform = TRUE, dist_unit = 'km', model = 'WGS84',
                       location = scale_pos)
    )
}

# save
pdf('<<<< FILEPATH REDACTED >>>>/manuscript_draft_figs_rtr.pdf')

c <- 'Sierra Leone'
for (p in c('pre', 'post', 'full')) plot_map()

c <- 'Mali'
for (p in c('pre', 'full')) {
  plot_map(250, 'bottomright')
}

c <- 'Senegal'
for (p in c('pre', 'full')) {
  plot_map(75, 'topright')
}
dev.off()
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Clean up, prep, and plot data for figure 5

# load and prep data
dia <- fread('<<<< FILEPATH REDACTED >>>>/had_diarrhea_admin_2_raked_prevalence_summary.csv')
dia <- dia[year == 2017 & (ADM0_NAME == 'Senegal' | ADM0_NAME == 'Mali')]
dia[, num_kids := mean*pop]
dt4 <- merge(ors, dia[, c('ADM2_CODE', 'year', 'num_kids')], by = 'ADM2_CODE', allow.cartesian = T)

no_ort_sub <- no_ort[((year == 2000 | year == 2017) & ADM0_NAME == 'Senegal') | 
                       (ADM0_NAME == 'Mali' & (year == 2001 | year == 2018))]
no_ort_sub[, year := ifelse(year < 2005, 'start', 'end')]
no_ort_sub <- dcast(no_ort_sub, formula = ADM2_CODE ~ year, value.var = 'mean')

dt4 <- merge(dt4, no_ort_sub[, c('ADM2_CODE', 'start')], by = 'ADM2_CODE')

# number of kids not treated with ORS
dt4[, no_ors := (1 - mean) * num_kids, by = .I]
tmp <- dt3[period == 'full']
dt4 <- merge(dt4, tmp[, c('ADM2_CODE', 'map_var')], by = 'ADM2_CODE')
dt4[map_var == '3. RHF not replaced', would_have_ort := no_ors*(1-start)]
dt4$map_var <- NULL

# plot map
plot_map <- function(map_var, 
                     scale_dist = 50,
                     scale_pos = 'bottomleft',
                     bins = c(0, 10, 250, 500, 1000, 1500, 2000, 2500, 5000, 10000),
                     map_colors = brewer.pal('Reds', n = 9)) {
    # subset shapefile
    shp_cty <- ad2_shp[which(ad2_shp@data$ADM0_NAME == c), ]
    shp_cty@data$ADM2_CODE <- as.numeric(shp_cty@data$ADM2_CODE)
    shp_cty <- st_as_sf(shp_cty)
    
    # add to data frame
    dt_map <- dt4[ADM0_NAME == c]
    dt_map[, variable := findInterval(get(map_var), bins)]
    map <- left_join(shp_cty, dt_map[, c('ADM2_CODE', 'variable', map_var), with = FALSE], by = 'ADM2_CODE')

    # get colors
    rr <- sort(unique(dt_map$variable))
    map_colors <- map_colors[rr]
    
    # plot
    print(
      ggplot(map) +
        geom_sf(aes(fill = as.factor(variable)), color = 'black', show.legend = FALSE) +
        scale_fill_manual(values = map_colors, name = '', na.value = 'grey') +
        theme_minimal() + labs(x = '', y = '') +
        coord_sf(datum = NA)
    )
}

pdf('<<<< FILEPATH REDACTED >>>>/manuscript_draft_fig5_rtr.pdf')

# senegal
c <- 'Senegal'
plot_map('no_ors', 75, 'topright')
plot_map('would_have_ort', 75, 'topright')

# mali
c <- 'Mali'
plot_map('no_ors', 225, 'bottomright')
plot_map('would_have_ort', 250, 'bottomright')

dev.off()
# --------------------------------------------------------------------------------
