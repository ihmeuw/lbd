# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: #REDACTED
# Date: 04/08/2020
# Purpose: Produce HAP paper results
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

#REDACTED

#use cairo to render instead of quartz (quartz causes big slowdowns with geom_sf)
if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}

#REDACTED
commondir       <- paste(core_repo, '#REDACTED', sep = '/')

#load packages
package_lib    <- sprintf('%s_code/_lib/pkg',h_root)
## Load libraries and  MBG project functions.
.libPaths(package_lib)
pacman::p_load(data.table, fst, scales, ggplot2, RColorBrewer, sf, viridis, farver, reldist) 
package_list    <- package_list <- fread('/#REDACTED/package_list.csv') %>% t %>% c

# Use setup.R functions to load common LBD packages and mbg_central "function" scripts
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#capture date
today <- Sys.Date() %>% gsub("-", "_", .)

#options
var_types <- c('lower', 'mean', 'upper')
new_gbd_estimates <- F

indicator_group <- 'cooking'
indicator <- 'cooking_fuel_solid'

config_par <- 'hap_sp_fine'
cov_par <- 'cooking_VNM'
type <- 'mean'
raked <- F
start_year <- 2000
end_year <- 2019
analysis_year <- 2018
cores <- 10
modeling_shapefile_version <- "2019_09_10"
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
data.dir <- file.path('/#REDACTED/cooking/post', run_date)
global_link_dir <- file.path('/#REDACTED/admin_shapefiles', modeling_shapefile_version) 
hap.paths <- data.table(admin2=file.path(data.dir, 'new_admin_2_summary.csv'))
hap.paths.d <- data.table(admin2=file.path(data.dir, 'admin_2_delta_summary.csv'))

###Output###
out.dir  <- file.path('/#REDACTED/maps', run_date) %T>% dir.create(recursive = T)
share.model.dir  <- file.path('/#REDACTED/input_data/')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
#PE functions#
file.path(my_repo, '_lib', 'post', 'map_fx.R') %>% source

#gbd fx
gbd.shared.function.dir <- '/#REDACTED'
file.path(gbd.shared.function.dir, 'get_location_metadata.R') %>% source
file.path(gbd.shared.function.dir, 'get_covariate_estimates.R') %>% source
file.path(gbd.shared.function.dir, 'get_age_metadata.R') %>% source
file.path(gbd.shared.function.dir, 'get_draws.R') %>% source
file.path(gbd.shared.function.dir, 'get_outputs.R') %>% source
file.path(gbd.shared.function.dir, 'get_population.R') %>% source

#helper fx to pull/prep the appropriate files from our list of SDG projection objects
prepCasts <- function(id, type, list=sdg_files, id_dt=NA, id_var=NA) {
  
  #format ID var if necessary
  if(nchar(id)==4) id <- as.character(id) #if the ID is a year, format as character
  
  #helper function to extract the correct object
  extractObj <- ifelse(type!='aroc',
                       function(x) list[[x]][[type]][[id]] %>% as.data.table,
                       function(x) list[[x]][[type]] %>% as.data.table ) #aroc only has one object
  
  #do the formatting and extractions
  lapply(1:length(list), extractObj) %>% 
    rbindlist %>% 
    { if(id_var %>% is.na) cbind(., id_dt[id,]) else .[, (id_var) := id] } %>% 
    return
  
}

#spatial functions
#return the max/min districts in a country and the distance between them
exploreRange <- function(explore_country, shp, dt, var, types=var_types, stub=NULL) {

 if (stub %>% is.null) new_vars <- paste(var, types, sep='_')
 else new_vars <- paste(var, types, stub, sep='_')
  
  out <- dt %>% 
    copy %>% 
    setnames(., new_vars, c('var_lower', 'var_mean', 'var_upper')) %>% 
    .[iso3==explore_country & !is.na(var_mean), .(var_lower, var_mean, var_upper, ADM0_NAME, ADM1_NAME,  ADM2_NAME)] %>% 
    .[order(var_mean)] %>% 
    .[c(1, nrow(.)), .(var_lower, var_mean, var_upper, ADM0_NAME, ADM1_NAME, ADM2_NAME)]

  st_distance(filter(shp, NAME_2==out[1, ADM2_NAME]), 
              filter(shp, NAME_2==out[2, ADM2_NAME])) %>% 
    as.numeric %>% 
    round %>% 
    {if (length(.)>1) message('warning: multipolygon, returning min distance'); min(.)} %>% 
    message('\nDistance is ', ./1e3, ' km')
  
  return(out)
  
}

#function to find polygons that share a border (rook neighbors)
st_rook <- function(a, b = a) st_relate(a, b, pattern = "F***1****") 

tabulateR <- function(results_dt=results, 
                      lvl='ad2', years=analysis_year, types='HAP', terms='lvl',
                      ind='prev', metrics=var_types, stub='',
                      filter=NULL, #should be provided as list(var=,vals=)
                      cleanup=T,
                      sorted='down') {
  
  #subset DT to requested results
  dt <- results_dt[dimension==lvl & year%in%years & type%in%types & term%in%terms]
  
  #define ind vars in order to remove irrelevant ones
  ind_vars <- names(dt) %>% .[. %like% paste(metrics, collapse='|')]
  irrel_vars <- ind_vars %>% .[!(. %like% ind)]
  rel_vars <- paste(ind, metrics, sep="_")
  mean_var <- rel_vars %>% .[. %like% 'mean']
  name_vars <- names(dt) %>% .[. %like% 'name|NAME']
  
  #remove irrelevant vars and missing rows
  dt[is.na(get(mean_var)), .N] %>% message('missing #', ., ' rows of ', mean_var)
  dt <- dt[!is.na(get(mean_var)), -c(irrel_vars), with=F] 
  setcolorder(dt, neworder=c('year', 'type', rel_vars, name_vars)) #reorder for legibility
  
  #filter if requested
  if(!is.null(filter)) { 
    if(filter$type) dt <- dt[get(filter$var) %in% filter$vals] #inclusive
    else dt <- dt[!(get(filter$var) %in% filter$vals)] #exclusive
  }
  
  #cleanup irrelevant vars if requested
  if(cleanup) dt <- Filter(function(x) !all(is.na(x)), dt)

  #return table, sorted if requested
  if(sorted=='down') setorderv(dt, cols=c('term', mean_var), order = -1) #descending order
  else if(sorted=='up') setorderv(dt, cols=c('term', mean_var)) #ascending order
  
  
  return(dt)
  
}

#helper functions
absChange <- function(x) x-data.table::shift(x,n=1)
relChange <- function(x) (x-data.table::shift(x, n=1))/data.table::shift(x,n=1)
#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
##read in and prep datasets for analysis##
## Read config file and save all parameters in memory
config <- set_up_config(repo            = my_repo,
                        indicator_group = indicator_group,
                        indicator       = indicator,
                        config_name     = paste0('/model/configs/config_', config_par),
                        covs_name       = paste0('/model/configs/covs_', cov_par),
                        run_tests       = F,
                        post_est_only   = T
)

#read in the proper annotations (borders, lakes, mask)
annotations_path <- file.path(out.dir, 'annotations.RDs')
check <- file.exists(annotations_path)
annotations <- ifelse(
  check,
  readRDS(annotations_path),
  load_map_annotations()
)
if(!check) saveRDS(annotations, file=annotations_path)

#merge sr region names/IDs/SDI quintiles
locs <- get_location_metadata(location_set_id = 35, gbd_round_id = 6) %>% 
  .[, .(iso3=ihme_loc_id, location_id, location_name, super_region_id, super_region_name, region_id, region_name)] #subset to relevant columns

locs_sdi <- get_location_metadata(location_set_id = 40, gbd_round_id = 6) %>% 
  .[, .(location_id, location_name, parent_id, level)] #subset to relevant columns

locs_sdi <- locs_sdi[level!=0] %>% 
  merge(., 
        locs_sdi[level==0, .(location_id, sdi_quintile=location_name)], 
        by.x='parent_id',
        by.y='location_id') %>% 
  .[, `:=` (parent_id=NULL, location_name=NULL, level=NULL)]

#read in shps
#REDACTED
adm2 <- rbind(stage1, stage2)

#read in mbg region info
stages <- file.path(j_root, '#REDACTED/stage_master_list.csv') %>% fread #read info about stages

#read in link_table
global_link_table <- file.path(global_link_dir, "lbd_full_link.rds") %>% readRDS %>% as.data.table
adm_links <- global_link_table[, .(ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_CODE)] %>% unique

#create file to crosswalk AD0 to iso3
iso3_map <- dplyr::select(adm2, iso3, ADM0_CODE=gadm_geoid) 
iso3_map$geometry <- NULL
iso3_map <- as.data.table(iso3_map) %>% unique
locs <- merge(locs, iso3_map, by='iso3')
adm_links <- merge(adm_links, locs, by=c('ADM0_CODE'), all.x=T)

#read in the input data
input_dt <- file.path(share.model.dir, "cooking_fuel_solid.csv") %>% fread
ker_dt <- file.path(share.model.dir, "cooking_fuel_kerosene.csv") %>% fread

#read in results
results <- file.path(data.dir, 'all_summary.fst') %>% read_fst(as.data.table=T)
dt <- results[dimension=='ad2'] %>% Filter(function(x) !all(is.na(x)), .)

# #merge sr region names/IDs
dt <- merge(dt, adm_links, by=c('ADM0_CODE', 'ADM2_CODE'), all.x=T)

#read in input data and prepare it for mapping
data <- load_map_results(indicator, indicator_group, run_date, raked, 
                   year_list=c(2000:2018),
                   custom_path = list('admin2'=dt),
                   geo_levels=c('admin2'),
                   cores=cores)
data_d <-
  load_map_results(indicator, indicator_group, run_date, raked, 
                   year_list=2018,
                   custom_path = list('admin2'=dt_d),
                   geo_levels=c('admin2'),
                   cores=cores)
#define extent of map
zoom.afr <- data.table(x1=-10, x2=50, y1=-20, y2=40)
zoom.global <- data.table(x1=-120, x2=150, y1=-40, y2=55)
#***********************************************************************************************************************

# ---LOAD GBD RESULTS---------------------------------------------------------------------------------------------------
#also pull the GBD2019 attributable results for comparison
#build the connector in order to define which locs to pull
# add connector object between ADM codes and GBD loc ids
connector <-
  data.table(get_gbd_locs(rake_subnational = T,
                          reg = "all",
                          shapefile_version = raking_shapefile_version))

ad0_dt <- results[dimension=='ad0' & term=='lvl' & type=='HAP'] %>% 
  Filter(function(x) !all(is.na(x)), .)

ad1_dt <- results[dimension=='ad1' & term=='lvl' & type=='HAP'] %>% 
  Filter(function(x) !all(is.na(x)), .)

#merge the connector, using the rake levels for most detailed identifier
#then collapse results to most detailed level
gbd_dt <-  
  list(
    merge(ad1_dt, 
          connector[rak_level==1, .(location_id, ADM0_CODE, ADM1_CODE, rak_level)], 
          by=c('ADM0_CODE','ADM1_CODE')),
    merge(ad0_dt, 
          connector[rak_level==0, .(location_id, ADM0_CODE, rak_level)], 
          by='ADM0_CODE')
  ) %>% rbindlist(use.names=T, fill=T)

#reload estimates from the central db if they have changed
if (new_gbd_estimates) {
  
  #REDACTED
  
} else {

  gbd <- file.path('#REDACTED', indicator_group, 'data', 
                   paste0('gbd_2019_best_', indicator, '_attrib_mortality.csv')) %T>% 
    message('reading GBD best estimates from this path\n', .) %>% 
    fread 
  
}

#merge ad2 results and GBD results using connector
gbd_dt <- merge(gbd_dt, gbd, by=c('location_id', 'year'))
#***********************************************************************************************************************

# ---GENERAL METRICS----------------------------------------------------------------------------------------------------
#number of LMICs
stages[Stage %in% c('1', '2a', '2b'), uniqueN(iso3)]

#number of countries with data
uniqueN(input_dt$ihme_loc_id)

#number of surveys
uniqueN(input_dt$nid)

#number of people
sum(input_dt$N, na.rm=T)

#number of points/polys
table(input_dt$point)

#explore kerosene data
ker_dt <- merge(ker_dt, locs, by.x='ihme_loc_id', by.y='iso3')
ker_dt[, ad0_kerosene := weighted.mean(cooking_fuel_kerosene, w=N), by=.(ihme_loc_id, year)]
ggplot(ker_dt, aes(x=year, y=ad0_kerosene, color=super_region_name)) +
  geom_hline(yintercept=.1, color='grey10', linetype='dashed') +
  geom_point() +
  geom_smooth() +
  scale_color_brewer(type='qual', palette = 'Paired') +
  facet_wrap(~ihme_loc_id) +
  theme_minimal()
#***********************************************************************************************************************

# ---SPATIAL VARIATION--------------------------------------------------------------------------------------------------
#country level patterns
#abs/rel change in prevalence
tabulateR(lvl='global', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018))
tabulateR(lvl='super_region', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018))
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018))

#district level
#counts based on threshold
threshold <- .95

#how many districts were above .95 in 2018
dt[year==analysis_year & type=='HAP' & term=='lvl', lapply(.SD, function(x) sum(x>threshold, na.rm=T)/.N), 
   .SDcols=paste0('prev_', var_types)]

#worst country
dt[year==analysis_year & type=='HAP' & term=='lvl', 
        .(count_lower=sum(prev_lower>threshold, na.rm=T), 
        count_mean=sum(prev_mean>threshold, na.rm=T),
        count_upper=sum(prev_upper>threshold, na.rm=T), N=.N), by=ADM0_NAME] %>% 
  .[, .(count_lower, count_mean, count_upper, N, pct=count_mean/N), by=ADM0_NAME] %>% 
  .[order(-pct)]

#worst country outside of SSA
dt[year==analysis_year & type=='HAP' & term=='lvl' & super_region_id!=166, 
   .(count_lower=sum(prev_lower>threshold, na.rm=T), 
     count_mean=sum(prev_mean>threshold, na.rm=T),
     count_upper=sum(prev_upper>threshold, na.rm=T), N=.N), by=ADM0_NAME] %>% 
  .[, .(count_lower, count_mean, count_upper, N, pct=count_mean/N), by=ADM0_NAME] %>% 
  .[order(-pct)]

#worst country in Americas
dt[year==analysis_year & type=='HAP' & term=='lvl' & super_region_id==103, 
   .(count_lower=sum(prev_lower>threshold, na.rm=T), 
     count_mean=sum(prev_mean>threshold, na.rm=T),
     count_upper=sum(prev_upper>threshold, na.rm=T), N=.N), by=ADM0_NAME] %>% 
  .[, .(count_lower, count_mean, count_upper, N, pct=count_mean/N), by=ADM0_NAME] %>% 
  .[order(-pct)]

#counts of exposed
dt[term=='lvl' & year==analysis_year & type=='HAP', .(count_lower=sum(prev_lower*pop_total, na.rm=T),
                                                      count=sum(prev_mean*pop_total, na.rm=T),
                                                      count_upper=sum(prev_upper*pop_total, na.rm=T)
                                                      ), by=.(ADM2_CODE, ADM2_NAME, ADM0_NAME, ADM1_NAME)] %>% 
  .[order(count)]
#***********************************************************************************************************************

# ---INEQUALITY---------------------------------------------------------------------------------------------------------
##produce metrics of inequality for 2017##
#calculate GINI/MAD at country level
dt_ineq <- dt[term=='lvl' & year==analysis_year & type=='HAP',
              .(iso3, year, ADM0_CODE, ADM2_CODE, ADM2_NAME, prev_mean, 
                super_region_id, super_region_name, region_id, region_name, pop_total)]
dt_ineq[, gini := gini(prev_mean, weights = pop_total), by=.(iso3, year)]
dt_ineq[, mean := weighted.mean(prev_mean, weights = pop_total), by=.(iso3, year)]
dt_ineq[, mad := mad(prev_mean, center = mean(prev_mean)), by=.(iso3, year)]
dt_ineq[, mean := mean(prev_mean, na.rm=T), by=.(iso3, year)]
dt_ineq[, max := max(prev_mean, na.rm=T), by=.(iso3, year)]
dt_ineq[, min := min(prev_mean, na.rm=T), by=.(iso3, year)]
dt_ineq[, range := max-min]

#range results
dt_ineq[year==analysis_year, .(range=max-min), by=iso3] %>% 
  unique %>% 
  .[order(range)]

exploreRange('MRT', shp=adm2, dt=dt[term=='lvl' & year==analysis_year & type=='HAP'], var='prev')

#AID results
#note that AID = gini * 2 * mean
aid.dt <- dt[cause=='lri' & grouping=='under5' & type=='HAP' & term=='lvl', 
             .(iso3, year, ADM0_CODE, ADM0_NAME, ADM2_CODE, ADM2_NAME, prev_mean, 
               super_region_id, super_region_name, region_id, region_name, pop_total)] %>% 
    na.omit(., cols='prev_mean') %>% 
  .[, mean := weighted.mean(prev_mean, weights = pop_total, na.rm=T), by=.(iso3, year)] %>% 
  .[, .(mean,
        aid=gini(prev_mean, weights=pop_total) * 2 * mean,
        pop=sum(pop_total, na.rm=T)), 
    by=.(super_region_id, region_name, iso3, year)] %>% 
  unique(by=c('iso3', 'year')) %>% 
  .[order(aid),] %>% 
  .[year==min(year), aid_start := aid] %>% 
  .[, aid_start := mean(aid_start, na.rm=T), by=.(iso3)] %>%
  .[, aid_d := (aid-aid_start)] %>% 
  .[, aid_dr := aid_d/aid_start]

#average change in AID
aid.dt[year==max(year)] %>% .[,weighted.mean(aid_dr, weights=pop)]

#regional changes
aid.dt[year==max(year) & region_name=='South Asia'] #south asia
aid.dt[year==max(year) & region_id==167] #central sub-saharan africa
#***********************************************************************************************************************

# ---TEMPORAL VARIATION-------------------------------------------------------------------------------------------------
#global results
#abs/rel change in prevalence
tabulateR(lvl='global', types='HAP', ind='prev', terms=c('change', 'change_rate', 'lvl'), years=c(2000, 2018))

#absolute decrease in exposure
tabulateR(lvl='global', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018)) %>% 
  .[, lapply(.SD, function(x) x * pop_total), .SDcols=paste('prev', var_types, sep='_'), by=year] %T>% 
  print %>% 
  .[, lapply(.SD, absChange), .SDcols=paste('prev', var_types, sep='_')] %>% .[2]

#absolute decrease in proportion at global level
tabulateR(lvl='global', types='HAP', ind='share', terms=c('change', 'change_rate', 'lvl'), years=c(2000, 2018))

#regional/ad0 results
tabulateR(lvl='super_region', types='HAP', ind='prev', terms=c('change_rate'))
tabulateR(lvl='region', types='HAP', ind='prev', terms=c('change', 'change_rate'))
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('change', 'change_rate'))

#best in SSA
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('change_rate'),
          filter=list(var='ADM0_CODE', vals=unique(locs[super_region_id==166, ADM0_CODE]), type=T),
          sorted='up')

#ad2 results
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('change_rate'))

#what percent of districts stagnated?
#counts based on threshold of improving less than 1%
threshold <- -.01
dt[year==analysis_year & term=='change' & type=='HAP', lapply(.SD, function(x) sum(x >= threshold, na.rm=T)), 
   .SDcols=paste('prev', var_types, sep='_')]
dt[year==analysis_year & term=='change' & type=='HAP', lapply(.SD, function(x) sum(x >= threshold, na.rm=T)/.N), 
   .SDcols=paste('prev', var_types, sep='_')]


#subset to upper quartile
q75_sfu <- results[dimension=='ad0' & term=='lvl' & year==min(year) & type=='HAP',  quantile(prev_mean, p=.75, na.rm=T)]
q75_countries <-  results[dimension=='ad0' & year==min(year) & type=='HAP' & prev_mean>=q75_sfu, unique(ADM0_NAME)]

dt[year==analysis_year & term=='change' & ADM0_NAME %in% q75_countries & type=='HAP', lapply(.SD, function(x) sum(x >= threshold)), 
     .SDcols=paste('prev', var_types, sep='_')]
dt[year==analysis_year & term=='change' & ADM0_NAME %in% q75_countries & type=='HAP', lapply(.SD, function(x) sum(x >= threshold)/.N), 
     .SDcols=paste('prev', var_types, sep='_')]

dt[year==analysis_year & term=='change' & type=='HAP'] %>% 
  .[, .(count=sum(prev_mean>threshold, na.rm=T), N=.N), by=ADM0_NAME] %>% 
  .[, .(pct=count/N), by=ADM0_NAME] %>% 
  .[order(pct)]

#with C.I
dt[year==analysis_year & term=='change' & type=='HAP'] %>% 
  .[, .(prev_lower=sum(prev_lower>threshold, na.rm=T),
        prev_mean=sum(prev_mean>threshold, na.rm=T),
        prev_upper=sum(prev_upper>threshold, na.rm=T),
        N=.N), by=ADM0_NAME] %>% 
  .[, lapply(.SD, function(x) x/N), .SDcols=paste('prev', var_types, sep='_'), by=ADM0_NAME] %>% 
  .[order(prev_mean)]

#range results
dt[year==analysis_year & term=='change' & type=='HAP', .(range=max(prev_mean, na.rm=T)-min(prev_mean, na.rm=T)), by=iso3] %>%
  unique %>%
  .[order(range)]

#explore the range to find inequality
exploreRange('BTN', shp=adm2, 
             dt=dt[term=='change' & year==analysis_year & type=='HAP'],
             var='prev')
#***********************************************************************************************************************

# ---SDG PROJECTIONS----------------------------------------------------------------------------------------------------
#SDG projection probabilities
#append the ADM2 files
sdg_files <-
  file.path(data.dir, 'sdg_projections') %>% list.files(pattern='admin_2', full.names = T) %>% 
  lapply(., readRDS)

#extract goal obj to index over
goals <- lapply(1:length(sdg_files), function(x) sdg_files[[x]]$goals) %>% rbindlist %>% unique

#create a dt with all probabilities
probs <- lapply(1:nrow(goals), prepCasts, type='probs', id_dt=goals) %>% 
  rbindlist %>% 
  setnames(c('target_year', 'spatial_idx'), c('year', 'ADM2_CODE')) %>% 
  merge(., adm_links[, .(ADM2_CODE, ADM0_CODE)], by='ADM2_CODE')  %>% 
  merge(., locs, by='ADM0_CODE') %>% 
  merge(., results[dimension=='ad2' & year==analysis_year & type=='HAP' & term=='lvl', pop_total, by=ADM2_CODE],
        by='ADM2_CODE')

#create a dt with all projections
projs <- 
  lapply(c(2018, seq(2020, 2030, 5)), prepCasts, type='proj', id_var='year') %>% 
  rbindlist %>% 
  setnames('spatial_idx', 'ADM2_CODE') %>% 
  melt(measure = patterns("V"), variable.name = "draw", value.name='sev')

#create a dt with aroc and combine
projs <- prepCasts(2018, type='aroc', id_var='year') %>% 
  melt(measure = patterns("V"), variable.name = "draw", value.name='aroc') %>% 
  merge(., projs, by=c('ADM2_CODE', 'year', 'draw'), all.y=T) %>% 
  merge(., adm_links[, .(ADM2_CODE, ADM0_CODE)], by='ADM2_CODE') 

#generate mean/ci
cols <- c('aroc', 'sev')
setkey(projs, year, ADM2_CODE)
projs[, paste0(cols, '_mean') := lapply(.SD, mean, na.rm=T), .SDcols=cols, by=key(projs)]
projs[, paste0(cols, '_lower') := lapply(.SD, quantile, probs=.025, na.rm=T), .SDcols=cols, by=key(projs)]
projs[, paste0(cols, '_upper') := lapply(.SD, quantile, probs=.975, na.rm=T), .SDcols=cols, by=key(projs)]
projs <- unique(projs, by=key(projs)) %>% 
  .[, c(cols, 'draw') := NULL]

#calculate relative uncertainty of SEV
projs[, sev_rel_uncertainty := (sev_upper-sev_lower)/sev_mean]
projs[sev_rel_uncertainty>2, sev_rel_uncertainty := 2] #cap at 2

#district level results
target_threshold <- .05
prob_threshold <- .95

#how many countries will succeed
ad0_probs <-
probs[target==target_threshold & year==2030, 
      .(prob=weighted.mean(absolute_goal_prob, w=pop_total), 
        pop=sum(pop_total, na.rm=T)), by=.(super_region_id, region_id, iso3)] %>% 
  .[, global_pop := sum(pop, na.rm=T)] %>% 
  merge(., iso3_map, by='iso3')

ad0_probs[, .(success=sum(prob>=prob_threshold, na.rm=T), fail=sum(prob<=(1-prob_threshold), na.rm=T))]
ad0_probs[prob>=prob_threshold, .(regs=uniqueN(region_id), pop=sum(pop, na.rm=T), global_pop=mean(global_pop))] %>% 
  .[, pop_share := pop/global_pop] %>% print
ad0_probs[prob<=(1-prob_threshold), .(regs=uniqueN(region_id), pop=sum(pop, na.rm=T), global_pop=mean(global_pop))] %>% 
  .[, pop_share := pop/global_pop] %>% print

#what is the country with the highest SFU in 2000 that is on track to succeed
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('lvl'), years=c(2000)) %>% 
  .[ADM0_CODE %in% ad0_probs[prob>=prob_threshold, unique(ADM0_CODE)]]

#how many districts have met in 2018
probs[target==target_threshold & year==2018] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), N=.N)] %>% 
  .[, .(count, N, pct=count/N)]
#what share of pop
probs[target==target_threshold & year==2018 & absolute_goal_prob>=prob_threshold, 
      sum(pop_total, na.rm=T)] / probs[target==target_threshold & year==2018, sum(pop_total, na.rm=T)]

#what share of unmet districts will meet between 2018-2030
unmet_districts <- probs[target==target_threshold & year==2018 & absolute_goal_prob<prob_threshold, unique(ADM2_CODE)]
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), N=.N)] %>% 
  .[, .(count, N, pct=count/N)]
#what share of pop
probs[target==target_threshold & year==2030 & absolute_goal_prob>=prob_threshold & ADM2_CODE%in%unmet_districts, 
      sum(pop_total, na.rm=T)] / probs[target==target_threshold & year==2018, sum(pop_total, na.rm=T)]

#superregional breakdown
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), N=.N), by=super_region_name] %>% 
  .[, .(count, N, pct=count/N), by=super_region_name] %>% 
  .[order(pct)]

#regional breakdown
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), N=.N), by=region_name] %>% 
  .[, .(count, N, pct=count/N), by=region_name] %>% 
  .[order(pct)]

#country breakdown
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), N=.N), by=iso3] %>% 
  .[, .(count, N, pct=count/N), by=iso3] %>% 
  .[order(pct)]

#country breakdown for SSA
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts & region_name %like% 'Latin'] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), N=.N), by=iso3] %>% 
  .[, .(count, N, pct=count/N), by=iso3] %>% 
  .[order(pct)]

#country failure breakdown for 2030
probs[target==target_threshold & year==2030] %>% 
  .[, .(count=sum(absolute_goal_prob<=(1-prob_threshold), na.rm=T), N=.N), by=iso3] %>% 
  .[, .(count, N, pct=count/N), by=iso3] %>% 
  .[order(pct)]

#country failure breakdown for 2030 (ESSA)
probs[target==target_threshold & year==2030 & region_id %in% c(174)] %>% 
  .[, .(count=sum(absolute_goal_prob<=(1-prob_threshold), na.rm=T), N=.N), by=iso3] %>% 
  .[, .(count, N, pct=count/N, resid=N-count), by=iso3] %>% 
  .[order(pct)]

#country with the largest divide
probs[target==target_threshold & year==2030] %>% 
  .[, .(count_fail=sum(absolute_goal_prob<=(1-prob_threshold), na.rm=T), 
        count_success=sum(absolute_goal_prob>=prob_threshold, na.rm=T),
        N=.N), by=iso3] %>% 
  .[, ratio := count_fail/count_success] %>% 
  .[count_fail!=0&count_success!=0] %>% 
  .[order(ratio)]
#***********************************************************************************************************************

# ---AIR POLLUTION------------------------------------------------------------------------------------------------------
##analyze relationship to AAP; TAP; HAP_SHARE##

##global levels/change
#current tap_pc/tap_paf global avg
tabulateR(lvl='global', types='TAP', ind='pm_pc', terms=c('lvl', 'change', 'change_rate'))  #highest TAP
tabulateR(lvl='global', types='TAP', ind='paf', terms=c('lvl', 'change', 'change_rate')) #highest TAP PAF
tabulateR(lvl='global', types='HAP', ind='share', terms=c('lvl', 'change', 'change_rate')) #highest HAP share

tabulateR(lvl='super_region', types='TAP', ind='pm_pc') #highest TAP
tabulateR(lvl='super_region', types='TAP', ind='paf', years=c(start_year, end_year)) #highest TAP PAF
tabulateR(lvl='super_region', types='HAP', ind='share', terms=c('lvl', 'change', 'change_rate')) #highest HAP share

tabulateR(lvl='region', types='TAP', ind='pm_pc') #highest TAP
tabulateR(lvl='region', types='TAP', ind='paf', years=c(start_year, end_year)) #highest TAP PAF
tabulateR(lvl='region', types='HAP', ind='share', terms=c('lvl', 'change')) #highest HAP share

#examine south asia
tabulateR(lvl='region', types='HAP', ind='share', terms=c('lvl', 'change'), years=c(start_year, end_year)) %>% 
  .[region_id==159]#change in HAP share in South Asia
tabulateR(lvl='region', types='HAP', ind='prev', terms=c('change_rate')) %>% 
  .[region_id==159]#change in SFU for South Asia
tabulateR(lvl='region', types='AAP', ind='pm_pc', terms=c('lvl', 'change'), years=c(start_year, end_year)) %>% 
  .[region_id==159]#change in AAP dose in South Asia

#deeper dive to countries in south asia
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('change_rate'), 
          filter=list(var='ADM0_CODE', vals=unique(locs[super_region_name=='South Asia', ADM0_CODE]), type=T))

##country levels
tabulateR(lvl='ad0', types='TAP', ind='pm_pc') #highest TAP
tabulateR(lvl='ad0', types='TAP', ind='pm_pc', 
          filter=list(var='super_region_id', vals=166, type=F)) #highest TAP excluding SSA

tabulateR(lvl='ad0', types='HAP', ind='share', 
          filter=list(var='super_region_id', vals=166, type=F)) #highest HAP share excluding SSA

tabulateR(lvl='ad0', types='TAP', ind='pm_pc', sorted='up') #lowest TAP
tabulateR(lvl='ad0', types='AAP', ind='share', sorted='up') #lowest AAP share
tabulateR(lvl='ad0', types='AAP', ind='share') #highest AAP share

#country level HAP:AAP ratio, what percent of countries is HAP the main contributor vs AAP
tabulateR(lvl='ad0', types='HAP', ind='share') %>% .[share_mean>.5, .N] #2018
tabulateR(lvl='ad0', types='HAP', ind='share', years=start_year) %>% .[share_mean>.5, .N] #2000

#pct of population that lives below WHO threshold 
threshold <- 10
dt[year==analysis_year & term=='lvl' & type=='TAP', 
   lapply(.SD, function(x) sum((x < threshold)*pop, na.rm=T)/sum(pop, na.rm=T)),
   .SDcols=paste0('pm_pc_', var_types)]

#pct of population that lives below WHO threshold interim-1
threshold <- 35
dt[year==analysis_year & term=='lvl' & type=='TAP', 
   lapply(.SD, function(x) sum((x < threshold)*pop, na.rm=T)/sum(pop, na.rm=T)),
          .SDcols=paste0('pm_pc_', var_types)]
dt[year==min(year) & term=='lvl' & type=='TAP', 
   lapply(.SD, function(x) sum((x < threshold)*pop, na.rm=T)/sum(pop, na.rm=T)),
   .SDcols=paste0('pm_pc_', var_types)]

#***********************************************************************************************************************

# ---ATTRIBUTABLE LRI---------------------------------------------------------------------------------------------------
##analyses of attributable LRI##

#what was the u5 mortality rate of LRI in regions where more than half of LRI was attributed to TAP

#how many children died from LRI attributable to TAP
tabulateR(lvl='global', types='TAP', ind='atr_count', years=c(analysis_year, start_year),
          terms=c('lvl', 'change', 'change_rate')) 
tabulateR(lvl='global', types=c('HAP', 'AAP'), ind='atr_count', years=c(analysis_year, start_year),
          terms=c('lvl', 'change', 'change_rate')) 

tabulateR(lvl='super_region', types='TAP', ind='atr_count') #highest TAP count
tabulateR(lvl='region', types='TAP', ind='atr_count') #highest TAP count
tabulateR(lvl='ad0', types='TAP', ind='atr_count') #highest TAP count

#how did the HAP share change
tabulateR(lvl='global', types='HAP', ind='share', terms=c('lvl', 'change', 'change_rate')) 

#how did the PAFs change
tabulateR(lvl='global', types='TAP', ind='paf', terms=c('lvl', 'change', 'change_rate'), years=c(analysis_year, start_year))

#PAFs by reg/country/district
tabulateR(lvl='super_region', types='TAP', ind='paf') 
tabulateR(lvl='region', types='TAP', ind='paf')
tabulateR(lvl='ad0', types='TAP', ind='paf', sorted='up') 
tabulateR(lvl='ad2', types='TAP', ind='paf')

#How many districts are majority TAP attributable
dt[year==analysis_year & type=='TAP' & term=='lvl', lapply(.SD, function(x) sum(x>.5, na.rm=T)/.N), 
   .SDcols=paste0('paf_', var_types)]

#proportion of total deaths from HAP
tabulateR(lvl='global', types='HAP', ind='atr_count', years=c(analysis_year, start_year)) %>%
  .[, c(3:5)] %>%
  .[] / (tabulateR(lvl='global', types='TAP', ind='atr_count', years=c(analysis_year, start_year)) %>% .[, c(3:5)])

#investigate countries where LRI rate fell but SFU was stable
tabulateR(lvl='ad0', types='TAP', ind='rate', terms='change_rate', years=c(analysis_year, start_year), sorted='up')
tabulateR(lvl='ad0', types='HAP', ind='prev', terms='change_rate', years=c(analysis_year, start_year))

#what is the lowest HAP paf in a country where LRI rate remains above 4/1000
results[year==analysis_year & dimension=='ad0' & type=='HAP' & term=='lvl' & rate_mean>2/1e3] %>% 
  .[order(paf_mean), .(ADM0_NAME, 
                        rate_lower=rate_lower*1e3, rate_mean=rate_mean*1e3, rate_upper=rate_upper*1e3, 
                        prev_lower, prev_mean, prev_upper,
                        paf_lower, paf_mean, paf_upper)]
#***********************************************************************************************************************