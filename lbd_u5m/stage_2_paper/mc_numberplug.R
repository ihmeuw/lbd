## #############################################################################
##
## STAGE 2 PAPER: MANUSCRIPT NUMBER PLUGGING
## Date: January 18, 2019
## Purpose: Number-plug the U5M Stage 2 manuscript. This task was split between
##   "<<<< NAME REDACTED >>>>"
##
## #############################################################################

core_repo <- "<<<< FILEPATH REDACTED >>>>"
ig_repo <- "<<<< FILEPATH REDACTED >>>>"
run_date <- "<<<< FILEPATH REDACTED >>>>"
shapefile_version <- "<<<< FILEPATH REDACTED >>>>"

output_path <- "<<<< FILEPATH REDACTED >>>>"

regions <- c("ansa+trsa-bra", "caca-mex", "cssa", "essa-yem", "mide+yem", "noaf", "ocea+seas-mys", "soas", "sssa", "stan+mng", "wssa")

remove_china <- T

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              SETUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("<<<< FILEPATH REDACTED >>>>")
pl <- c("rgeos",      "data.table", "raster",     "rgdal",      "INLA",
        "seegSDM",    "seegMBG",    "dismo",      "gbm",        "foreign",
        "parallel",   "doParallel", "grid",       "gridExtra",  "pacman",
        "gtools",     "glmnet",     "ggplot2",    "RMySQL",     "plyr",
        "tictoc",     "dplyr",      "magrittr",   "tidyr",      "sp",
        "sf",         "matrixStats", "fasterize")
mbg_setup(package_list = pl, repos = core_repo)
source("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")

adm_level = 2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Get those numbers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stage_list <- fread("<<<< FILEPATH REDACTED >>>>")
region_keeper <- stage_list[gadm_geoid %in% get_adm0_codes(regions),]

gbd <- get_gbd_q(age = "under5", year_list = 2017)

gbd <- merge(gbd, stage_list, by.x = "name", by.y = "loc_id", all.x = T)
gbd <- gbd[Stage %in% c("1", "2a", "2b")]
maxrow <- gbd[mean == max(mean),] 
minrow <- gbd[mean == min(mean),] 

#Intro, paragraph 1, last sentence
line_44 <- sprintf("In 2017, national mortality rates for children in LMICs ranged from %d deaths per 1000 livebirths in %s to %d deaths per 1000 livebirths in the %s.", round(maxrow$mean * 1000), maxrow$location_name, round(minrow$mean * 1000), minrow$location_name)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop <- get_population(age_group_id = 1, 
                      location_id = c(1, gbd[Stage == 1,]$name), 
                      year_id = 2017, 
                      sex_id = 3, 
                      gbd_round_id = 5)

africa_pop <- pop[location_id != 1,]
africa_pop <- sum(africa_pop$population)
nonafrica_per <- round(((pop[location_id == 1, population] - africa_pop) / pop[location_id == 1, population]) * 100)

lt <- data.table(get_envelope(age_group_id = 1, 
                              year_id = 2017, 
                              sex_id = 3, 
                              location_id = c(1, gbd[Stage == 1,]$name), 
                              gbd_round_id = 5,
                              with_shock = 1,
                              with_hiv = 1))

africa_deaths <- lt[location_id != 1,]
africa_deaths <- sum(africa_deaths$mean)
nonafrica_deaths_per <- round(((lt[location_id == 1, mean] - africa_deaths) / lt[location_id == 1, mean]) * 100)

#Precision public health and child mortality, Paragraph 2, end
line_61 <- sprintf("However, heterogeneity at subnational scales remains undescribed for the areas outside of Africa which encompass %d%% of the worlds' children and %d%% of total global child deaths.", nonafrica_per, nonafrica_deaths_per)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#data that feeds directly into the model
train_list <- list()
for(reg in regions){
  train <-fread(sprintf("<<<< FILEPATH REDACTED >>>>")
  subsetter <- region_keeper[gadm_geoid %in% get_adm0_codes(reg),]$iso3
  train2 <- train[country %in% subsetter,]
  train_list[[reg]] <- train2
}
com <- do.call("rbind", c(train_list, list(fill = T)))

com <- com[year >= 2000,]
com_ceb <- com[age == 1, sum(N)]
com_ced <- com[, sum(died)]

com_ceb <- com[, sum(N)]

#Precision public health and child mortality, last paragraph, 2nd to last line
line_74 <- sprintf("containing records of approximately %.1f million births and %.1f million child deaths from 2000 to 2017.", round(com_ceb / 1000000,1), round(com_ced / 1000000,1))

#take unprepared and do histogram by year of survey of ceb and cbh year of birth

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Finding locations in stage 1 and stage 2 where the proportion of adm2s with q > 0.1 is highest
adm2_raked <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
admin_shp <- sf::st_read(get_admin_shapefile(admin_level = 2, version = shapefile_version), quiet=T)
admin_shp_data <- as.data.table(admin_shp)
adm2 <- merge(adm2_raked, admin_shp_data, by = "ADM2_CODE", all.x = T)

adm2 <- adm2[, adm2s := .N, by = c("ADM0_NAME", "ADM0_CODE")]
adm2_sum <- adm2[mean > 0.1, .(count = .N), by = c("ADM0_NAME", "ADM0_CODE", "adm2s")]
adm2_sum <- merge(adm2_sum, stage_list, by.x = "ADM0_CODE", by.y = "gadm_geoid", all.x = T)
adm2_sum[, per_100 := count / adm2s]

ssa <- adm2_sum[mbg_reg %in% c("cssa", "wssa", "sssa", "essa"),]
setorder(ssa, -"per_100")

non_ssa <- adm2_sum[!(mbg_reg %in% c("cssa", "wssa", "sssa", "essa", "noaf")),]
setorder(non_ssa, -"per_100")

ssa_countries <- paste(as.character(ssa$ADM0_NAME[1:4]), collapse=", ")
non_ssa_countries <- paste0(paste(sort(as.character(non_ssa$ADM0_NAME[1:4])), collapse=", "), ", and ", non_ssa$ADM0_NAME[5])

#Unequal rates of under-5 mortality, 1st paragraph, last sentence
line_118 <- sprintf("Despite this broad success in mortality rate reduction, the highest rates of under-5 mortality in 2017 were still largely concentrated where rates were highest in 2000 (Fig. 1a and 1b). We estimated that there were mortality rates in excess of 100 deaths per 1,000 live births across large geographic regions within western and central sub-Saharan Africa (SSA), and within %s (Fig. 1b.)", non_ssa_countries)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
adm2_2017 <- adm2[year == 2017]
max_adm2_2017 <- adm2[mean == max(adm2_2017$mean, na.rm=T),]

#Unequal rates of under-5 mortality, 2nd_paragraph, middle
line_101 <- sprintf("By 2017, the second administration sub-division with the highest mortality across LMICs was %s in %s state, %s, with %.1f deaths (%.1f-%.1f) per 1,000 live births.",max_adm2_2017$ADM2_NAME, max_adm2_2017$ADM1_NAME, max_adm2_2017$ADM0_NAME, max_adm2_2017$mean * 1000, max_adm2_2017$lower * 1000, max_adm2_2017$upper * 1000)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Unequal Burden of child deaths, last paragraph, end
adm2_raked <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
admin_shp <- sf::st_read(get_admin_shapefile(admin_level = 2, version = shapefile_version), quiet=T)
admin_shp_data <- as.data.table(admin_shp)
adm2 <- merge(adm2_raked, admin_shp_data, by = "ADM2_CODE", all.x = T)

adm2_sum <- adm2[year == 2017,.(mean = mean(mean, na.rm= T), lower = mean(lower, na.rm = T), upper = mean(upper, na.rm = T)), by = c("ADM0_NAME", "ADM0_CODE", "year")]
adm2_sum[, CI := upper - lower]

adm2_mer <- merge(adm2_sum, stage_list, by.x = "ADM0_CODE", by.y = "gadm_geoid")

ssa_high_u_high_q <- adm2_mer[mbg_reg %in% c("cssa", "wssa", "sssa", "essa"),]
non_ssa_high_u_high_q <- adm2_mer[!(mbg_reg %in% c("cssa", "wssa", "sssa", "essa", "noaf")),]

setorder(ssa_high_u_high_q , -CI, -mean)
setorder(non_ssa_high_u_high_q, -CI, -mean)

low_u_low_q <- adm2_mer
setorder(low_u_low_q , CI, mean)

low_u_low_q[1:10, spr_reg_nm]

low_u_high_q <- adm2_mer[mbg_reg %in% c("caca", "ansa", "trsa"),]
setorder(low_u_high_q, -mean, CI)

line_244 <- sprintf("High mortality rates and high uncertainty were estimated across much of %s and %s in sub-Saharan Africa, and in %s, %s, and %s outside of sub-Saharan Africa (Figure X). In contrast, we found that under-5 mortality rates and uncertainty associated with those estimates were relatively low across large areas of Latin and South America. Within those broad geographic regions, large areas with contrary arrangements can be found: either high rates of mortality and low uncertainty (such as within country, country, and country) or low mortality rates but with high uncertainty around those estimates (such as within country, country, and country).", as.character(ssa_high_u_high_q$ADM0_NAME)[[1]], as.character(ssa_high_u_high_q$ADM0_NAME)[[2]], as.character(non_ssa_high_u_high_q$ADM0_NAME)[[1]], as.character(non_ssa_high_u_high_q$ADM0_NAME)[[2]], as.character(non_ssa_high_u_high_q$ADM0_NAME)[[3]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(remove_china){
  adm2 <- adm2[ADM0_NAME != "China"]
}
adm2_2000 <- adm2[year == 2000,]
adm2_2017 <- adm2[year == 2017,]
total_adm2_units <- nrow(adm2_2000)
adm2_2000_greater_than_100q <- nrow(adm2_2000[mean > .08,])
adm2_2017_greater_than_100q <- nrow(adm2_2017[mean > .08,])
adm2_2000_greater_than_100_percent <- round((adm2_2000_greater_than_100q / total_adm2_units) * 100, 1)
adm2_2017_greater_than_100_percent <- round((adm2_2017_greater_than_100q / total_adm2_units) * 100, 1)


line_103 <- sprintf("Overall, the total number of second administrative divisions with an under-5 mortality rate in excess of 100 deaths per 1,000 live births decreased from %.1f%% (%i out of %i) of second administrative divisions units in 2000 to %.1f%% (%i out of %i) of these areas units in 2017.", adm2_2000_greater_than_100_percent, adm2_2000_greater_than_100q, total_adm2_units, adm2_2017_greater_than_100_percent, adm2_2017_greater_than_100q, total_adm2_units)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lt_deaths <- data.table(get_envelope(age_group_id = 1, 
                              year_id = c(2000, 2017),
                              sex_id = 3, 
                              location_id = region_keeper$loc_id, 
                              gbd_round_id = 5,
                              with_shock = 1,
                              with_hiv = 1))
deaths_2000 <- sum(lt_deaths[year_id == 2000,]$mean, na.rm = T)
deaths_2017 <- sum(lt_deaths[year_id == 2017,]$mean, na.rm = T)
fewer_deaths <- round((deaths_2000 - deaths_2017) / 100000) /10
per_fewer_deaths <- round(((deaths_2000 - deaths_2017) / deaths_2000) * 100)

million_2000 <- round(deaths_2000 / 1000000, 1)
million_2017 <- round(deaths_2017 / 1000000, 1)

lt_2017 <- data.table(get_envelope(age_group_id = 1, 
                              year_id = 2017, 
                              sex_id = 3, 
                              location_id = region_keeper$loc_id, 
                              gbd_round_id = 5,
                              with_shock = 1,
                              with_hiv = 1))

lt_2017 <- merge(lt_2017, stage_list, by.x = "location_id", by.y = "loc_id")
setorder(lt_2017, -mean)
first_2017 <- lt_2017[1, c("mean", "upper", "lower", "location_name")]
second_2017 <- lt_2017[2, c("mean", "upper", "lower", "location_name")]
third_2017 <- lt_2017[3, c("mean", "upper", "lower", "location_name")]
fourth_2017 <- lt_2017[4, c("mean", "upper", "lower", "location_name")]

#Prospects for future reductions in child mortality, 1st paragraph, start
line_135 <- sprintf("The goal of mortality reduction efforts is ultimately to prevent premature deaths, and not just to reduce mortality rates. Across the countries studied here, there were %.1f million (%i%%) fewer under-5 deaths in 2017 than in 2000 (Fig. 4) (%.1f million compared to %.1f million). At the national level, the largest number of child deaths in 2017 occurred in %s (%s deaths [%s-%s]), %s (%s deaths [%s-%s]), %s (%s deaths [%s-%s]), and %s (%s deaths [%s-%s])", fewer_deaths, per_fewer_deaths, million_2017, million_2000, first_2017$location_name[[1]], format(signif(first_2017$mean[[1]],3), big.mark=",", scientific=FALSE), format(signif(first_2017$lower[[1]],3), big.mark=",", scientific=FALSE), format(signif(first_2017$upper[[1]],3), big.mark=",", scientific=FALSE), second_2017$location_name[[1]], format(signif(second_2017$mean[[1]],3), big.mark=",", scientific=FALSE), format(signif(second_2017$lower[[1]],3), big.mark=",", scientific=FALSE), format(signif(second_2017$upper[[1]],3), big.mark=",", scientific=FALSE), third_2017$location_name[[1]], format(signif(third_2017$mean[[1]],3), big.mark=",", scientific=FALSE), format(signif(third_2017$lower[[1]],3), big.mark=",", scientific=FALSE), format(signif(third_2017$upper[[1]],3), big.mark=",", scientific=FALSE), fourth_2017$location_name[[1]], format(signif(fourth_2017$mean[[1]],3), big.mark=",", scientific=FALSE), format(signif(fourth_2017$lower[[1]],3), big.mark=",", scientific=FALSE), format(signif(fourth_2017$upper[[1]],3), big.mark=",", scientific=FALSE))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prospects for future reductions in child mortality, second paragraph, all

assemble_q_d <- function(run_date, adm_level){
  # Define q and d input files
  in_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  fps <- list(
    q = paste0("<<<< FILEPATH REDACTED >>>>"),
    d = paste0("<<<< FILEPATH REDACTED >>>>"),
    q_lower = paste0("<<<< FILEPATH REDACTED >>>>"),
    d_lower = paste0("<<<< FILEPATH REDACTED >>>>"),
    q_upper = paste0("<<<< FILEPATH REDACTED >>>>"),
    d_upper = paste0("<<<< FILEPATH REDACTED >>>>")
  )
  dts <- list(q=NA, d=NA)
  # Read in and prep data.tables for each entity (q and d)
  for(val in names(dts)){
    dts[[val]] <- fread(fps[[val]])
    lower <- fread(fps[[paste0(val, "_lower")]])
    upper <- fread(fps[[paste0(val, "_upper")]])
    setnames(lower, "value", paste0("lower_", val))
    setnames(upper, "value", paste0("upper_", val))
    dts[[val]] <- merge(dts[[val]], lower, by = c("ADM2_CODE", "year"))
    dts[[val]] <- merge(dts[[val]], upper, by = c("ADM2_CODE", "year"))
    
    setnames(dts[[val]], paste0('ADM',adm_level,"_CODE"), "ADM_CODE")
    dts[[val]] <- dts[[val]][, .(ADM_CODE, year, value,  get(paste0("lower_", val)),  get(paste0("upper_", val)))]
    setnames(dts[[val]], c('value', "V4", "V5"), c(val, paste0("lower_", val), paste0("upper_", val)))
    dts[[val]][, level := adm_level]
    
  }
  # Merge q and d by admin code and year
  q_d_merged <- merge(
    x   = dts[['q']],
    y   = dts[['d']],
    by  = c("ADM_CODE","year","level"),
    all = TRUE
  )
  # Return prepped data.table
  return(q_d_merged)
}


## Read and prep file
dq <- assemble_q_d(
  run_date = run_date,
  adm_level = 2
)

# Clean and format columns
dq[, year_factor := as.factor(year)]
dq[, q := q * 1000][,lower_q := lower_q * 1000][,upper_q := upper_q * 1000]
dq[, d := d / 1E5][,lower_d := lower_d / 1E5][,upper_d := upper_d / 1E5]
dq <- dq[order(year,q)]
dq[, q_bin := round(q/5,0)*5]

# For first and last years only, merge on region
# Load tables connecting admin2 units to countries and countries to regions
dq_firstlast <- dq[year %in% c(2000, 2017)]

hierarchy <- admin_shp_data
setnames(hierarchy, paste0('ADM',2,'_CODE'), "ADM_CODE")
# Subset to only subnational admin units and national admin codes
hierarchy <- hierarchy[, .(ADM_CODE, ADM0_CODE)]
hierarchy$ADM_CODE <- as.integer(as.character(hierarchy$ADM_CODE))
hierarchy$ADM0_CODE <- as.integer(as.character(hierarchy$ADM0_CODE))

# - Admin0 to country name and region
by_stage <- region_keeper
setnames(by_stage, 'gadm_geoid','ADM0_CODE')

# - Combine to get merge table from adm2 to all metadata
all_meta <- merge(
  x   = hierarchy,
  y   = by_stage,
  by  = c('ADM0_CODE')
)

all_meta$ADM_CODE <- as.integer(as.character(all_meta$ADM_CODE))
# Merge onto the admin2 table
dq_meta <- merge(
  x     = dq_firstlast,
  y     = all_meta,
  by    = c('ADM_CODE'),
  all.x = TRUE
)

if(remove_china){
  dq_meta <- dq_meta[location_name != "China"]
}

total_deaths <- dq_meta[,.(d = sum(d * 1E5, na.rm = T)), by = "year"]

less_than_80_2000 <- dq_meta[(q < 80) & (year == 2000), .(d80 = sum(d * 1E5, na.rm = T), d80_lower = sum(lower_d * 1E5, na.rm = T), d80_upper = sum(upper_d * 1E5, na.rm = T)), by = c("year")]

percent_deaths_2000_under_80 <- (less_than_80_2000$d80[[1]]/ total_deaths[year == 2000, d]) * 100
less_than_80_2000_d <- format(signif(less_than_80_2000$d80[[1]],3), big.mark=",", scientific=FALSE)
less_than_80_2000_d_lower <- format(signif(less_than_80_2000$d80_lower[[1]],3), big.mark=",", scientific=FALSE)
less_than_80_2000_d_upper <- format(signif(less_than_80_2000$d80_upper[[1]],3), big.mark=",", scientific=FALSE)

less_than_80_2017 <- dq_meta[(q < 80) & (year == 2017), .(d80 = sum(d * 1E5, na.rm = T), d80_lower = sum(lower_d * 1E5, na.rm = T), d80_upper = sum(upper_d * 1E5, na.rm = T)), by = c("year")]

percent_deaths_2017_under_80 <- (less_than_80_2017$d80[[1]]/ total_deaths[year == 2017, d]) * 100

less_than_75_2017 <- dq_meta[(q < 75) & (year == 2017), .(d75 = sum(d * 1E5, na.rm = T), d75_lower = sum(lower_d * 1E5, na.rm = T), d75_upper = sum(upper_d * 1E5, na.rm = T)), by = c("year")]
less_than_75_2017_d <- format(signif(less_than_75_2017$d75[[1]],3), big.mark=",", scientific=FALSE)
less_than_75_2017_d_lower <- format(signif(less_than_75_2017$d75_lower[[1]],3), big.mark=",", scientific=FALSE)
less_than_75_2017_d_upper <- format(signif(less_than_75_2017$d75_upper[[1]],3), big.mark=",", scientific=FALSE)

max_stage2_q <- ceiling(max(dq_meta[Stage != 1 & year == 2017, q], na.rm = T))
more_than_max_stage2_2017 <- dq_meta[(q > max_stage2_q) & (year == 2017), .(d = sum(d * 1E5, na.rm = T), d_lower = sum(lower_d * 1E5, na.rm = T), d_upper = sum(upper_d * 1E5, na.rm = T)), by = c("year")]
more_than_max_stage2_2017_d <- format(signif(more_than_max_stage2_2017$d[[1]],3), big.mark=",", scientific=FALSE)
more_than_max_stage2_2017_d_lower <- format(signif(more_than_max_stage2_2017$d_lower[[1]],3), big.mark=",", scientific=FALSE)
more_than_max_stage2_2017_d_upper <- format(signif(more_than_max_stage2_2017$d_upper[[1]],3), big.mark=",", scientific=FALSE)

ssa_2017 <- dq_meta[mbg_reg %in% c("sssa", "cssa", "essa-yem", "wssa") & year == 2017]
ssa_d_2017 <- sum(ssa_2017$d * 1E5, na.rm = T)
percent_ssa_d <- (more_than_max_stage2_2017$d[[1]] / ssa_d_2017) * 100

overs <- paste0('V',1:1000)

#find cutoff of q where the same proportion of deaths occured as in Africa where q is greater than occurs in stage 2 
#sort soas by q and cumsum deaths
soas_2017 <- dq_meta[mbg_reg %in% c("soas") & year == 2017]
setorder(soas_2017, q)
soas_2017_agg <- soas_2017[, d_cum := cumsum(d)]
soas_2017_agg[, d_cum := d_cum / max(d_cum)]
cutoff_soas <- soas_2017_agg[d_cum > ((100 - percent_ssa_d) / 100), q][[1]]

draws_soas <- readRDS("<<<< FILEPATH REDACTED >>>>")
soas_names <- soas_2017[q > cutoff_soas[[1]], ADM_CODE]
draws_soas_agg_over <- draws_soas[ADM2_CODE %in% soas_names & year == 2017,]
soas_d_2017 <- sum(soas_2017$d * 1E5, na.rm = T)

draws_soas_agg <- draws_soas_agg_over[, lapply(overs, function(x) sum(get(x), na.rm = T))]
draws_soas_per <- draws_soas_agg / soas_d_2017
soas_per_upper <- quantile(draws_soas_per[1,], .975)[[1]] * 100
soas_per_lower <- quantile(draws_soas_per[1,], .025)[[1]] * 100

less_than_25_2017 <- dq_meta[(q < 25) & (year == 2017), .(d25 = sum(d * 1E5, na.rm = T), d25_lower = sum(lower_d * 1E5, na.rm = T), d25_upper = sum(upper_d * 1E5, na.rm = T)), by = c("year")]
less_than_25_2017_d <- format(signif(less_than_25_2017$d25[[1]],3), big.mark=",", scientific=FALSE)
less_than_25_2017_d_lower <- format(signif(less_than_25_2017$d25_lower[[1]],3), big.mark=",", scientific=FALSE)
less_than_25_2017_d_upper <- format(signif(less_than_25_2017$d25_upper[[1]],3), big.mark=",", scientific=FALSE)

percent_deaths_2017_under25 <- (less_than_25_2017$d25[[1]]/ total_deaths[year == 2017, d]) * 100
percent_deaths_2017_under25_upper <- (less_than_25_2017$d25_upper[[1]]/ total_deaths[year == 2017, d]) * 100
percent_deaths_2017_under25_lower <- (less_than_25_2017$d25_lower[[1]]/ total_deaths[year == 2017, d]) * 100

less_than_25_2000 <- dq_meta[(q < 25) & (year == 2000), .(d25 = sum(d * 1E5, na.rm = T), d25_lower = sum(lower_d * 1E5, na.rm = T), d25_upper = sum(upper_d * 1E5, na.rm = T)), by = c("year")]
less_than_25_2000_d <- format(signif(less_than_25_2000$d25[[1]],3), big.mark=",", scientific=FALSE)
less_than_25_2000_d_lower <- format(signif(less_than_25_2000$d25_lower[[1]],3), big.mark=",", scientific=FALSE)
less_than_25_2000_d_upper <- format(signif(less_than_25_2000$d25_upper[[1]],3), big.mark=",", scientific=FALSE)

percent_deaths_2000_under25 <- (less_than_25_2000$d25[[1]]/ total_deaths[year == 2000, d]) * 100
percent_deaths_2000_under25_upper <- (less_than_25_2000$d25_upper[[1]]/ total_deaths[year == 2000, d]) * 100
percent_deaths_2000_under25_lower <- (less_than_25_2000$d25_lower[[1]]/ total_deaths[year == 2000, d]) * 100

less_than_25_2017_region <- dq_meta[(q < 25) & (year == 2017), .(d25 = sum(d * 1E5, na.rm = T), d25_lower = sum(lower_d * 1E5, na.rm = T), d25_upper = sum(upper_d * 1E5, na.rm = T)), by = c("year", "spr_reg_nm")]
setorder(less_than_25_2017_region, -d25)
less_than_25_2017_first <- format(signif(less_than_25_2017_region[1, d25],3), big.mark=",", scientific=FALSE)
less_than_25_2017_second <- format(signif(less_than_25_2017_region[2, d25],3), big.mark=",", scientific=FALSE)
less_than_25_2017_first_name <- less_than_25_2017_region[1, spr_reg_nm]
less_than_25_2017_second_name <- less_than_25_2017_region[2, spr_reg_nm]
less_than_25_2017_ssa <- less_than_25_2017_region[spr_reg_nm == "Sub-Saharan Africa", d25]

per_less_than_25_ssa <- (less_than_25_2017_ssa / ssa_d_2017) * 100
line_153 <- sprintf("In 2000, %.1f%% of child deaths in LMIC - representing %s deaths (%s-%s) - occurred where mortality rates were less than 80 deaths per 1,000 live births in that year(Fig. 4). By comparison, in 2017, %.1f%% of child deaths occurred in areas where the mortality rate was below 80 deaths per 1,000 live births. A growing proportion of under-5 deaths are occurring in 'low' mortality areas; %.1f%% (%.1f%%-%.1f%%) of all under-5 deaths in 2017 occurred in locations where the mortality rate was below the SDG target rate of 25 deaths per 1,000 live births, compared to %.1f%% (%.1f%%-%.1f%%) in 2000."
, percent_deaths_2000_under_80, less_than_80_2000_d, less_than_80_2000_d_lower, less_than_80_2000_d_upper, percent_deaths_2017_under_80, percent_deaths_2017_under25, percent_deaths_2017_under25_lower, percent_deaths_2017_under25_upper, percent_deaths_2000_under25, percent_deaths_2000_under25_lower, percent_deaths_2000_under25_upper)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pull_age_binned_deaths <- function(run_date, adm_level) {
  infant <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  neonatal <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  under5 <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  
  setnames(neonatal, "value", "neonatal")
  setnames(infant, "value", "infant")
  setnames(under5, "value", "under5")
  
  under1 <- merge(neonatal, infant, by = c("ADM2_CODE", "year"))
  under1[, value := infant - neonatal]
  under1[, age := "Post-Neonatal"]
  child <- merge(infant, under5, by = c("ADM2_CODE", "year"))
  child[, value := under5 - infant]
  child[, age := "Child"]
  
  setnames(neonatal, "neonatal", "value")
  neonatal[,age := "Neonatal"]
  neonatal[,V1 := NULL]
  
  under1 <- under1[, c("ADM2_CODE", "year", "value", "age")]
  child <- child[, c("ADM2_CODE", "year", "value", "age")]
  
  all_ages <- rbind(neonatal, under1, child)
  return(all_ages)
}

d <- pull_age_binned_deaths(run_date, adm_level)
dq <- assemble_q_d(
  run_date = run_date,
  adm_level = adm_level
)
dq[,d:=NULL]
dq <- merge(d, dq, by.x = c("ADM2_CODE", "year"), by.y = c("ADM_CODE", "year"))
setnames(dq, "value", "d")

# Clean and format columns
dq[, year_factor := as.factor(year)]
dq[, q := q * 1000][,lower_q := lower_q * 1000][,upper_q := upper_q * 1000]
dq[, d := d / 1E5][,lower_d := lower_d / 1E5][,upper_d := upper_d / 1E5]
dq[, q_bin := round(q/5,0)*5]

# For first and last years only, merge on region
# Load tables connecting admin2 units to countries and countries to regions
dq_firstlast <- dq[year %in% c(2000, 2017)]

hierarchy <- admin_shp_data

# Subset to only subnational admin units and national admin codes
hierarchy <- hierarchy[, .(ADM_CODE, ADM0_CODE)]
hierarchy$ADM_CODE <- as.integer(as.character(hierarchy$ADM_CODE))
hierarchy$ADM0_CODE <- as.integer(as.character(hierarchy$ADM0_CODE))

# - Combine to get merge table from adm2 to all metadata
all_meta <- merge(
  x   = hierarchy,
  y   = by_stage,
  by  = c('ADM0_CODE')
)

setnames(dq_firstlast, "ADM2_CODE", "ADM_CODE")
all_meta$ADM_CODE <- as.integer(as.character(all_meta$ADM_CODE))
# Merge onto the admin2 table
dq_meta <- merge(
  x     = dq_firstlast,
  y     = all_meta,
  by    = c('ADM_CODE'),
  all.x = TRUE
)

if(remove_china){
  dq_meta <- dq_meta[location_name != "China"]
}

dq_meta_2000 <- dq_meta[year==2000,]
dq_meta_2000 <- dq_meta_2000[Stage != 3,]
dq_meta_2017  <- dq_meta[year==2017,]
dq_meta_2017 <- dq_meta_2017[Stage != 3,]

dq_meta_2000_age <- dq_meta_2000[, .(d=sum(d, na.rm = T)), by=c("age")]
dq_meta_2017_age <- dq_meta_2017[, .(d=sum(d, na.rm = T)), by=c("age")]
nn_pro_2000 <- round((sum(dq_meta_2000_age[1, "d"], na.rm = T) / sum(dq_meta_2000_age[, "d"])) * 100, 1)
nn_pro_2017 <- round((sum(dq_meta_2017_age[1, "d"]) / sum(dq_meta_2017_age[, "d"])) * 100, 1)

dq_meta_2017_age_under50 <- dq_meta_2017[q_bin < 50, .(d=sum(d, na.rm = T)), by=c("age")]
dq_meta_2017_age_over50 <- dq_meta_2017[q_bin > 50, .(d=sum(d, na.rm = T)), by=c("age")]
nn_pro_2017_under50 <- round((sum(dq_meta_2017_age_under50[1, "d"]) / sum(dq_meta_2017_age_under50[, "d"])) * 100, 1)
nn_pro_2017_over50 <- round((sum(dq_meta_2017_age_over50[1, "d"]) / sum(dq_meta_2017_age_over50[, "d"])) * 100, 1)

line_224 <- sprintf("Neonatal mortality rates have declined but failed to keep pace with reductions in rates of older children, leading to a higher proportion of under-5 deaths occurring within the first four weeks of life: from %.1f%% in 2000 to %.1f%% in 2017", nn_pro_2000, nn_pro_2017)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Write Text to file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("<<<< FILEPATH REDACTED >>>>"), append = T)
cat("#Intro, paragraph 1, last sentence\n")
cat(line_44)
cat("\n\n (line 61) Precision public health and child mortality, Paragraph 2, middle\n")
cat(line_61)
cat("\n\n (line 74) Precision public health and child mortality, last paragraph, start\n")
cat(line_74)
cat("\n\n (line 101) Unequal rates of under-5 mortality, 2nd_paragraph, middle\n")
cat(line_101)
cat("\n\n (line 103) Unequal rates of under-5 mortality, 2nd_paragraph, middle\n")
cat(line_103)
cat("\n\n (line 118) Unequal rates of under-5 mortality, 2nd to last paragraph\n")
cat(line_118)
cat("\n\n (line 135) Distribution of under-5 deaths may not follow rates, 1st paragraph, start\n")
cat(line_135)
cat("\n\n (line 153) Distribution of under-5 deaths may not follow rates, second paragraph\n")
cat(line_153)
cat("\n\nDiscussion, limitations, and future work, 3rd paragraph, middle\n")
cat(line_224)
cat("\n\n# (line 244) Discussion- do manually/ double check\n")
cat(line_244)
sink()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         SI Raking factors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("<<<< FILEPATH REDACTED >>>>")

subnats <- c(4841, 43908, 43872, 4842, 43909, 43873, 4843, 43910, 43874, 4844, 43911, 43875, 4846, 43913, 43877, 4849, 43916, 43880, 4850, 43917, 43881, 4851, 43918, 43882, 4852, 43919, 43883, 4853, 43920, 43884, 4854, 43921, 43885, 4855, 43922, 43886, 4856, 43923, 43887, 4857, 43924, 43888, 4859, 43926, 43890, 4860, 43927, 43891, 4861, 43928, 43892, 4862, 43929, 43893, 4863, 43930, 43894, 4864, 43931, 43895, 4865, 43932, 43896, 4867, 43934, 43898, 4868, 43935, 43899, 4869, 43936, 43900, 4870, 43937, 43901, 4871, 43938, 43902, 4872, 43939, 43903, 44538, 44539, 44540, 4873, 43940, 43904, 4874, 43941, 43905, 4875, 43942, 43906, 4709, 4726, 4717, 4725, 4715, 4737, 4720, 4713, 4721, 4722, 4724, 4729, 4731, 4730, 4732, 4719, 4718, 4716, 4739, 4740, 4727, 4728, 4742, 4741, 4712, 4738, 4735, 4734, 4736, 4733, 4711, 4714, 4710, 4723
)

file.list <- list.files(, pattern = '*rf.csv$')

## Loop through and bind everything together
for(i in 1:length(file.list)) {
  file.name <- file.list[i]
  mydata <- fread(file.name)
  message(paste("processing file ", i, " of ", length(file.list) ))
  if(file.name == file.list[1]){
    combined <- mydata 
  } else{
    combined <- rbind(combined, mydata, fill=TRUE)
  }
}

sink("<<<< FILEPATH REDACTED >>>>"), append = T)
print("Under 5 raking factors, national and subnational")
print(median(combined$raking_factor, na.rm=T))
print(summary(combined$raking_factor,na.rm=T))

combined <- combined[loc %in% subnats,]
print(median(combined$raking_factor, na.rm=T))
      print(summary(combined$raking_factor,na.rm=T))
sink()

setwd("<<<< FILEPATH REDACTED >>>>")

file.list <- list.files(, pattern = '*rf.csv$')

## Loop through and bind everything together
for(i in 1:length(file.list)) {
  file.name <- file.list[i]
  mydata <- fread(file.name)
  message(paste("processing file ", i, " of ", length(file.list) ))
  if(file.name == file.list[1]){
    combined <- mydata 
  } else{
    combined <- rbind(combined, mydata, fill=TRUE)
  }
}

sink("<<<< FILEPATH REDACTED >>>>"), append = T)
print("neonatal raking factors, national and subnational")
print(median(combined$raking_factor, na.rm=T))
print(summary(combined$raking_factor,na.rm=T))

combined <- combined[loc %in% subnats,]
print(median(combined$raking_factor, na.rm=T))
print(summary(combined$raking_factor,na.rm=T))
sink()

setwd("<<<< FILEPATH REDACTED >>>>")

file.list <- list.files(, pattern = '*rf.csv$')

## Loop through and bind everything together
## Loop through and bind everything together
for(i in 1:length(file.list)) {
  file.name <- file.list[i]
  mydata <- fread(file.name)
  message(paste("processing file ", i, " of ", length(file.list) ))
  if(file.name == file.list[1]){
    combined <- mydata 
  } else{
    combined <- rbind(combined, mydata, fill=TRUE)
  }
}

sink("<<<< FILEPATH REDACTED >>>>"), append = T)
print("infant raking factors, national and subnational")
print(median(combined$raking_factor, na.rm=T))
print(summary(combined$raking_factor,na.rm=T))

combined <- combined[loc %in% subnats,]
print(median(combined$raking_factor, na.rm=T))
print(summary(combined$raking_factor,na.rm=T))
sink()