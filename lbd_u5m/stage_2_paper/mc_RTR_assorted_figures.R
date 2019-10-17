
##LOAD IN DQ_META FROM EXPLORATORY_SHOW_Q_BY_D.R BEFORE RUNNING THIS
library(scales)
out_dir <- "<<<< FILEPATH REDACTED >>>>"

source("<<<< FILEPATH REDACTED >>>>")

S2_BOUNDS <- list(
  'lat_min'  = -57,
  'lat_max'  = 63,
  'long_min' = -129,
  'long_max' = 165
)

pull_format_admin_shp <- function(
  adm_level,
  shp_version = 'current',
  simp_prop   = NULL,
  stage2_only = FALSE,
  include_non_gbd = TRUE
){
  # Load the shapefile, simplifying Admin 0 shapefiles a bit by default
  if(is.null(simp_prop) & (adm_level==0)) simp_prop <- 0.05
  shp <- fast_load_shapefile(
    shp_path = get_admin_shapefile(admin_level=adm_level, version=shp_version),
    simplify_tol=simp_prop
  )
  # Pull the admin0 lookup table 
  lookup_table <- load_adm0_lookup_table()
  adm_field <- ifelse(
    detect_adm_shapefile_date_type(
      shpfile_path=get_admin_shapefile(admin_level=0, version=shp_version)
    )$shpfile_type=='gaul',
    'GAUL_CODE',
    'gadm_geoid'
  )
  setnames(lookup_table, adm_field, 'ADM0_CODE')
  if(include_non_gbd==FALSE){
    # Drop all non-GBD locations, which will not have a valid location ID
    lookup_table <- lookup_table[ loc_id > 0, ]
  }
  if(stage2_only){
    lookup_table <- lookup_table[ Stage %in% c('1','2a','2b'), ]
  }
  # Fortify for ggplotting, but don't add fields
  shp$admin_code <- shp@data[, paste0('ADM',adm_level,'_CODE')]
  shp@data       <- shp@data[, c('ADM0_CODE','admin_code')]
  shp_fort <- prep_shp_data_for_mapping(
    shp       = shp,
    dataset   = lookup_table[, .(ADM0_CODE, location_name)],
    merge_var ='ADM0_CODE'
  )
  # Subset to the Stage 2 boundaries
  shp_fort[lat  < S2_BOUNDS$lat_min,  lat  := S2_BOUNDS$lat_min  ]
  shp_fort[lat  > S2_BOUNDS$lat_max,  lat  := S2_BOUNDS$lat_max  ]
  shp_fort[long < S2_BOUNDS$long_min, long := S2_BOUNDS$long_min ]
  shp_fort[long > S2_BOUNDS$long_max, long := S2_BOUNDS$long_max ]
  message(sprintf("Done prepping admin%s data.\n",adm_level))
  return(shp_fort)
}

modeling_shapefile_version <- "<<<< FILEPATH REDACTED >>>>"

bg_for_mapping_adm0 <- pull_format_admin_shp(
  adm_level   = 0,
  shp_version = modeling_shapefile_version,
  simp_prop   = NULL,
  stage2_only = FALSE,
  include_non_gbd = TRUE
)

bg_shp <- pull_format_admin_shp(
  adm_level   = 2,
  shp_version = modeling_shapefile_version,
  simp_prop   = NULL,
  stage2_only = TRUE,
  include_non_gbd = TRUE
)

drop_countries <- c("French Guiana", "Western Sahara", "Malaysia")
bg_shp <- bg_shp[!(location_name %in% drop_countries),]

mapping_2017 <- dq_meta[year == 2017,]
mapping_2000 <- dq_meta[year == 2000,]

mapping_2017 <- merge(mapping_2017, bg_shp, by.x = "ADM_CODE", by.y = "admin_code")
mapping_2017[, over80 := "no"]
mapping_2017[q > 80, over80 := "yes"]
mapping_2017$over80 <- as.factor(mapping_2017$over80)

mapping_2000 <- merge(mapping_2000, bg_shp, by.x = "ADM_CODE", by.y = "admin_code")
mapping_2000[, over80 := "no"]
mapping_2000[q > 80, over80 := "yes"]
mapping_2000$over80 <- as.factor(mapping_2000$over80)


admin2030 <- ggplot() +
  geom_polygon(data = mapping_2000,
               aes(x=long, y=lat, group=group, fill=as.factor(over80)), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_manual(values = c("grey", "#009999"), name = "", labels = c("5q0 Under 80", "5q0 Over 80"))+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  )

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()


admin2030 <- ggplot() +
  geom_polygon(data = mapping_2017,
               aes(x=long, y=lat, group=group, fill=as.factor(over80)), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_manual(values = c("grey", "#009999"), name = "", labels = c("5q0 Under 80", "5q0 Over 80"))+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  )

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()



#######################################################################################################################
# population disaggregation

pop <- fread("<<<< FILEPATH REDACTED >>>>", drop = "V1")
dq_pop <- merge(dq_meta, pop, by.x = c("ADM_CODE", "year"), by.y = c("ADM2_CODE", "year"))

dq_agg <- dq_pop[, .(pop = sum(pop, na.rm = T)), by = c("year", "q_bin")]

for(i in seq(200,315, 5)){
  dq_agg <- rbind(dq_agg, list(2017, i, 1))
}
dq_agg$year <- as.factor(dq_agg$year)
dq_agg$pop <- dq_agg$pop / 1000000

line_colors <- c("orangered4", "dodgerblue4")
plot <- ggplot(dq_agg, aes(q_bin, pop, fill=year)) + 
  geom_bar(stat = "identity", position = 'dodge') +
  theme_minimal() +
  theme(legend.position='right',
        legend.title = element_blank()) +
  scale_fill_manual(values = line_colors) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_x_continuous(breaks = c(0,27.5,102.5,202.5,302.5), expand=c(0,0), labels = c(0,25,100,200,300))+
  coord_cartesian(xlim = c(0, 315)) +
  xlab("Under-5 Mortality Rate (per 1,000 LB)") +
  ylab("Under-5 Population (millions)")
  

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(plot)
dev.off()

#########################################################################################################################
# absolute difference maps
pops <- fread("<<<< FILEPATH REDACTED >>>>", drop = "V1")

dq_meta[, q_bin := round(q/10,0)*10]

total <- dq_meta_pop[, .(total = sum(pop, na.rm = T)), by = "year"]
dq_meta_pop <- merge(dq_meta, pops, by.x = c("ADM_CODE", "year"), by.y = c("ADM2_CODE", "year"), all.x =T)

dq_meta_pop[, q_bin := NULL]
dq_meta_pop <- merge(dq_meta_pop, adm2_bin, by = "ADM_CODE")

#reshape functions are for losers
dq_meta_pop_2000 <- dq_meta_pop[year == 2000,]
dq_meta_pop_2017 <- dq_meta_pop[year == 2017,]
setnames(dq_meta_pop_2000, c("d", "pop"), c("d_2000", "pop_2000"))
setnames(dq_meta_pop_2017, c("d", "pop"), c("d_2017", "pop_2017"))
dq_meta_pop_2000 <- dq_meta_pop_2000[, !c("year")]
dq_meta_pop_2017 <- dq_meta_pop_2017[, !c("year")]
dq_meta_pop_2000 <- dq_meta_pop_2000[, c("ADM_CODE", "d_2000", "pop_2000")]
dq_meta_pop_2017 <- dq_meta_pop_2017[, c("ADM_CODE", "d_2017", "pop_2017")]
dq_meta_pop_spread <- merge(dq_meta_pop_2000, dq_meta_pop_2017, by = "ADM_CODE", all = T)

dq_meta_pop_spread[, d_dis := d_2000 * pop_2017 / pop_2000]

mapping <- merge(dq_meta_pop_spread, bg_shp, by.x = "ADM_CODE", by.y = "admin_code")
mapping[, pop_diff := pop_2017 - pop_2000]
mapping[, d_dis_d_2000 := d_dis - d_2000]
mapping[, d_2017_d_dis := d_2017 - d_dis]
mapping[, d_2017_d_2000 := d_2017 - d_2000]

mapping[, pop_diff_rel := pop_2017 / pop_2000]
mapping[, d_dis_d_2000_rel := d_dis / d_2000]
mapping[, d_2017_d_dis_rel := d_2017 / d_dis]
mapping[, d_2017_d_2000_rel := d_2017 / d_2000]

total <- dq_meta_pop_spread[, .(d_2000 = sum(d_2000, na.rm=T), d_dis = sum(d_dis, na.rm=T), d_2017 = sum(d_2017, na.rm = T))]

#divide everything by 10 because deaths are already / 100,000
line_86 <- sprintf("The decline in child deaths has been a function of mortality rate reduction outpacing population growth. For example, if the 2017 under-5 population was exposed to the mortality rates observed in 2000, there would have been %.1f million deaths in 2017, which is %.1f million more deaths than in 2000 (%.1f million), but rather we observe %.1f million deaths in 2017.", total$d_dis[[1]] / 10, (total$d_dis[[1]] / 10) - (total$d_2000[[1]]/10), total$d_2000[[1]]/10, total$d_2017[[1]]/10)

sink("<<<< FILEPATH REDACTED >>>>", append = T)
cat("\n\nline 86 pop disaggregation\n")
cat(line_86)
sink()

admin2030 <- ggplot() +
  geom_polygon(data = mapping,
               aes(x=long, y=lat, group=group, fill=pop_diff), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_gradient2(midpoint = 0, mid = "lightyellow", high = muted("red"), low = muted("blue"), na.value = "grey50")+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  ) +
  labs(fill = "2017 pop -\n2000 pop")

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()


admin2030 <- ggplot() +
  geom_polygon(data = mapping,
               aes(x=long, y=lat, group=group, fill=d_dis_d_2000), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_gradient2(midpoint = 0, mid = "lightyellow", high = muted("red"), low = muted("blue"), na.value = "grey50")+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  ) +
  labs(fill = "Expected 2017 deaths with \npopulation growth -\n2000 deaths")

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()


admin2030 <- ggplot() +
  geom_polygon(data = mapping,
               aes(x=long, y=lat, group=group, fill=d_2017_d_dis), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_gradient2(midpoint = 0, mid = "lightyellow", high = muted("red"), low = muted("blue"), na.value = "grey50")+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  ) +
  labs(fill = "2017 deaths - \nExpected 2017 deaths with \npopulation growth")

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()

admin2030 <- ggplot() +
  geom_polygon(data = mapping,
               aes(x=long, y=lat, group=group, fill=d_2017_d_2000), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_gradient2(midpoint = 0, mid = "lightyellow", high = muted("red"), low = muted("blue"), na.value = "grey50")+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  ) +
  labs(fill = "2017 deaths - \n2000 deaths")

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()


############## Relative ########################
admin2030 <- ggplot() +
  geom_polygon(data = mapping,
               aes(x=long, y=lat, group=group, fill=pop_diff_rel), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_gradient2(midpoint = 1, mid = "lightyellow", high = muted("red"), low = muted("blue"), na.value = "grey50")+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  ) +
  labs(fill = "2017 pop /\n2000 pop")

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()


admin2030 <- ggplot() +
  geom_polygon(data = mapping,
               aes(x=long, y=lat, group=group, fill=d_dis_d_2000_rel), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_gradient2(midpoint = 1, mid = "lightyellow", high = muted("red"), low = muted("blue"), na.value = "grey50")+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  ) +
  labs(fill = "Expected 2017 deaths with \npopulation growth /\n2000 deaths")

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()


admin2030 <- ggplot() +
  geom_polygon(data = mapping,
               aes(x=long, y=lat, group=group, fill=d_2017_d_dis_rel), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_gradient2(midpoint = 1, mid = "lightyellow", high = muted("red"), low = muted("blue"), na.value = "grey50")+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  ) +
  labs(fill = "2017 deaths / \nExpected 2017 deaths with \npopulation growth")

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()

admin2030 <- ggplot() +
  geom_polygon(data = mapping,
               aes(x=long, y=lat, group=group, fill=d_2017_d_2000_rel), show.legend = T) + 
  geom_path(data = bg_for_mapping_adm0,
            aes(x=long, y=lat, group=group),
            color='#000000',
            lwd=0.1) +
  scale_fill_gradient2(midpoint = 1, mid = "lightyellow", high = muted("red"), low = muted("blue"), na.value = "grey50")+
  theme_map() +
  theme(legend.position = c(0.95, 0.80),
        legend.title = element_text(size=9, color='#222222'),
        legend.text  = element_text(size=6, color='#222222')) +
  coord_map(
    projection = "mollweide",
    xlim = c(-108, 157.5), 
    ylim = c(-35.5,  53)
  ) +
  labs(fill = "2017 deaths / \n2000 deaths")

# Plot
png(
  sprintf("<<<< FILEPATH REDACTED >>>>"), 
  height=5.5, width=12, units='in', res=800
)
print(admin2030)
dev.off()

