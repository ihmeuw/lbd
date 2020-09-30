# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: #REDACTED
# Date: 11/01/2019
# Purpose: Calculate TAP PAFs
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

#REDACTED

#load external packages
pacman::p_load(assertthat, ccaPP, data.table, dplyr, fst, mgsub, sf, stringr, magrittr)

#detect if running interactively
interactive <- F  %>% #manual override
  ifelse(., T, !length(commandArgs())>2) %>%  #check length of arguments being passed in
  ifelse(., T, !(is.na(Sys.getenv("RSTUDIO", unset = NA)))) #check if IDE

if (interactive) {
  
  #REDACTED
  
} else {
  
  ## Set repo location, indicator group, and some arguments
  user            <- commandArgs()[4]
  core_repo       <- commandArgs()[5]
  indicator_group <- commandArgs()[6]
  indicator       <- commandArgs()[7]
  config_par      <- commandArgs()[8]
  cov_par         <- commandArgs()[9]
  region          <- commandArgs()[10]
  run_date        <- commandArgs()[11]
  measure         <- commandArgs()[12]
  holdout         <- as.numeric(commandArgs()[13])
  my_repo         <- commandArgs()[14]
  age             <- 0
  
}

# collect date
today <- Sys.Date()
#***********************************************************************************************************************

# ---OPTIONS------------------------------------------------------------------------------------------------------------
## set arguments
reformat <- F #set T if needing to reformat the cellpreds to long data.table
prep_gbd_files <- F #set T if change on GBD side requires reprep of GBD info

#countries we want to process by state 
big_countries <- c(
  143, #Mexico
  44, #China
  47, #DRC
  65, #Algeria
  138, #Libya
  193, #Sudan
  33, #Brazil
  108, #Iran
  150, #Mongolia
  105 #India
)

super_big_countries <- c(44, 33) #run these states singlecore

#mbg options
rk                       <- (indicator=='cooking_fuel_solid')
suffix                   <- '_eb_bin0_0'
draw_number              <- 50

#covariates to pull
covs <- data.table(
  covariate = c("ihmepm25"),
  measure = c("median"),
  gbd = c(FALSE),
  include = c(TRUE),
  release = c("2019_11_05")
)
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
#load gbd fxs
#REDACTED

# load MBG packages
#REDACTED
package_list <- c(t(read.csv('#REDACTED/package_list.csv',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# load post functions
file.path(my_repo, '_lib', 'post', 'format_cell_preds.R') %>% source
file.path(my_repo, '_lib', 'post', 'risk_calculations.R') %>% source
file.path(my_repo, '_lib', 'post', 'make_projections.R') %>% source

#load erf custom fx
file.path(h_root, '_code/risks/erf/air/_lib/misc.R') %>% source
file.path(h_root, '_code/risks/erf/air/paf/_lib/paf_helpers.R') %>% source

#use your own diacritics fx, due to inscrutable error
#note: requires mgsub pkg
fix_diacritics <<- function(x) {
  
  #first define replacement patterns as a named list
  defs <-
    list('Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 
         'Ç'='C', 'È'='E', 'É'='E','Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 
         'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U','Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Þ'='B', 
         'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'ç'='c','è'='e', 'é'='e', 'ê'='e', 
         'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
         'ö'='o', 'ø'='o', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y', 'ß'='Ss')
  
  #then force conversion to UTF-8 and replace with non-diacritic character
  enc2utf8(x) %>% 
    mgsub(., pattern=enc2utf8(names(defs)), replacement = defs) %>% 
    return
  
}

#custom function for saving country level fst file
saveCountry <- function(country, type, dt, debug=F) {
  
  if (debug) browser()
  
  out.path <- sprintf('%s/%s_%s.fst', paste0(outputdir, '/tmp'), country, type)
  
  #write fst file for country if data available for this indicator
  if (dt[ADM0_CODE==country, .N]>0) write.fst(dt[ADM0_CODE==country], path = out.path) #write file in fst format
  else message('no ', type, ' data present for #', country)
  
  return(NULL)
  
}

#custom fucntion for verifying files all exist
checkTmpFiles <- function(x){ sapply(x, function(i) paste0(outputdir, '/tmp/', i, '_sfu.fst') %>% file.exists) %>% all}
#***********************************************************************************************************************

# ---MBG PREP-----------------------------------------------------------------------------------------------------------
#Prep MBG inputs/Load Data
PID <- Sys.getpid()
tic("Entire script") # Start master timer

## Set seed for reproducibility
message('Setting seed 98118 for reproducibility')
set.seed(98118)

## Read config file and save all parameters in memory
config_filepath <- 'cooking/model/configs/'
config <- set_up_config(repo            = my_repo,
                        indicator_group = indicator_group,
                        indicator       = indicator,
                        config_name     = paste0('/model/configs/config_', config_par),
                        covs_name       = paste0('/model/configs/covs_', cov_par),
                        run_tests       = F,
                        post_est_only   = T
                        )

# Get the necessary variables out from the config object into global env
rake_countries <- eval(parse(text = config[V1 == 'rake_countries', V2]))
rake_subnational <- eval(parse(text = config[V1 == 'subnational_raking', V2]))
modeling_shapefile_version <- config[V1 == 'modeling_shapefile_version', V2]
raking_shapefile_version <- config[V1 == 'raking_shapefile_version', V2]
countries_not_to_rake <- config[V1 == 'countries_not_to_rake', V2]
countries_not_to_subnat_rake <- config[V1 == 'countries_not_to_subnat_rake', V2]
year_list <- eval(parse(text = config[V1 == 'year_list', V2]))
metric_space <- config[V1 == 'metric_space', V2]
summstats <- eval(parse(text = config[V1 == 'summstats', V2]))

#***********************************************************************************************************************

# ---IN/OUT-------------------------------------------------------------------------------------------------------------
#input dirs
#REDACTED
global_link_dir <- file.path('#REDACTED', modeling_shapefile_version)

#intermediate data
data.dir <- file.path('#REDACTED', indicator_group, 'data')

# create outputdir
outputdir <- paste0('#REDACTED', indicator_group, '/post/', run_date)
  dir.create(paste0(outputdir, '/tmp/'), recursive = T)
  
# prep GBD files and values ------------------------------------------------------------
##define share directory
share_dir <- paste0('#REDACTED', indicator_group, '/', indicator, '/output/', run_date, '/')

## start master timer
tic('Master timer')

#***********************************************************************************************************************

# ---DATA PREP----------------------------------------------------------------------------------------------------------
if(prep_gbd_files) {

  tic('Processing GBD datasets')
  # load the gbd location hierarchy
  # note that these IDs are the reporting hierarchy for GBD2019
  locs <- get_location_metadata(location_set_id = 35, gbd_round_id = 6) %>% 
    .[, .(location_id, iso3=ihme_loc_id)] #subset to relevant columns
  
  # load and format the GBD results
  gbd.hap.version <- '092419'
  gbd.hap.pm <- file.path('/#REDACTED', gbd.hap.version, 
                          paste0('lm_map_', gbd.hap.version, '.csv')) %>% 
    fread %>% 
    .[grouping!='indoor', .(location_id, year=year_id, hap_excess_pm25=mean, grouping)] %>% 
    merge(., locs, by='location_id', all.x=T) %>% 
    .[nchar(iso3)<4] %>% 
    .[, ADM0_CODE := get_adm0_codes(iso3), by=iso3] #merge on ad0 code
  
  #duplicate children to make specific u5 category
  gbd.hap.pm <- list(
    gbd.hap.pm,
    copy(gbd.hap.pm[grouping=='child']) %>% .[,grouping := 'under5']
  ) %>% rbindlist
  
  #choose ages to merge for RRs
  gbd.hap.pm[grouping=='under5', `:=` (age_group_id=2, sex_id=1)]
  gbd.hap.pm[grouping=='child', `:=` (age_group_id=2, sex_id=1)]
  gbd.hap.pm[grouping=='female', `:=` (age_group_id=10, sex_id=2)]
  gbd.hap.pm[grouping=='male', `:=` (age_group_id=10, sex_id=1)]
  
  #load the rr max for SEV
  rr.max <- fread('/#REDACTED/air_hap.csv') %>% 
    setnames(., 'rr', 'rr_max') %>% 
    .[, draw := paste0("V", draw)] %>% 
    .[, rei := NULL] %>% #to save space
    setnames('cause_id', 'cause')
  
  #replace the cause ids with their code
  cause.codes <- c('cvd_ihd',
                   "cvd_stroke_isch",
                   "cvd_stroke_intracerebral",
                   "cvd_stroke_subarachnoid",
                   "lri",
                   'neo_lung',
                   'resp_copd',
                   't2_dm')
  
  cause.ids <- c(493,
                 495,
                 496,
                 497,
                 322,
                 426,
                 509,
                 976)
  
  # then pass to your custom function
  rr.max <- findAndReplace(rr.max,
                           cause.ids,
                           cause.codes,
                           "cause", 
                           'cause')
  
  #keep only one type of stroke so as not to triple weight this in the SEV collapse
  rr.max <- rr.max[!(cause %in% c('cvd_stroke_intracerebral', 'cvd_stroke_subarachnoid'))]
  rr.max[cause %like% 'stroke', cause := 'cvd_stroke']
  
  rr.dt <- 
    list.files(rr.dir) %>% 
    lapply(formatRR, dir=rr.dir) %>% 
    rbindlist(use.names=T, fill=T)
  
  #remove the RRs that are mediated through shifts
  #these include low birthweight and gestational age
  #we will use the GBD pafs for these instead
  rr.dt <- rr.dt[!(cause %in%  c('bw', 'ga', 'lbw', 'ptb'))]
  
  #expand to the appropriate age groups using a custom function
  rr.dt <- unique(rr.dt$cause) %>% 
    lapply(causeExpansion, input.table=rr.dt) %>% 
    rbindlist(use.names=T, fill=T)
  
  #now subset to the age/sex combos that you have specific hap pm data for so as not to create a massive join
  #these are children, women, and men
  rr.dt <- rr.dt[(age_group_id==2 & sex_id==1) | age_group_id==10]
  
  #prep for merge with hap data
  setnames(rr.dt, 'exposure_spline', 'tap_pc')
  
  #combine necessary GBD info
  gbd_data <- list(
    'pm'=gbd.hap.pm,
    'rr'=rr.dt,
    'rr_max'=rr.max
  )

  #save 
  saveRDS(gbd_data,
          file=file.path(data.dir, 'gbd_data.RDS'))
  
  #cleanup
  rm(gbd.hap.pm, 
     rr.dt, 
     rr.max)

} else gbd_data <- file.path(data.dir, 'gbd_data.RDS') %>% readRDS

#read in link_table
global_link_table <- file.path(global_link_dir, "lbd_full_link.rds") %>% readRDS %>% as.data.table
adm_links <- global_link_table[, .(ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_NAME, ADM2_CODE)] %>% unique

## find the adm0s/iso3s of countries in this region
adm0s <- get_adm0_codes(region, shapefile_version = modeling_shapefile_version)
countries <- pull_custom_modeling_regions(region) %>% unlist %>% str_split(., pattern='\\+', simplify = T)

#***********************************************************************************************************************

# ---FORMAT PREDS-------------------------------------------------------------------------------------------------------
#format cell preds into long DTs
if(!reformat & adm0s %>% checkTmpFiles) {
  
  message('skipping format stage as results are preformatted')
    
} else {

  ## create the long data.tables and format them for the TAP calculations
  tic('Make table(s) for HAP')
  
  hap.dt <- link_cell_pred(ind_gp = indicator_group,
                           ind = indicator,
                           rd = run_date,
                           reg = region,
                           measure = measure,
                           pop_measure = 'hap', #men/women/children/under5
                           covs = covs,
                           n_draws = draw_number, #reduce the # of draws to speed up processing
                           year_start = year_list %>% min,
                           year_end = year_list %>% max,
                           rk = rk,
                           shapefile_version = modeling_shapefile_version,
                           coastal_fix = T,
                           debug=F)
  toc(log = TRUE)
  
  tic('Saving results at country level with .fst')
  lapply(adm0s, saveCountry, type='sfu', dt=hap.dt)
  rm(hap.dt) #save memory
  toc(log=TRUE)
  
}
#***********************************************************************************************************************

# ---TAP/RISK CALC------------------------------------------------------------------------------------------------------
#custom function to calculate the ad0/2 level TAP results
#produces the TAP PAF for LRI, the LRI rate attributable to TAP, and the LRI counts attributable to TAP
calcTAP <- function(country, dir, rr_data=rr.dt,
                    adm_info=adm_links,
                    hap_data=gbd_data,
                    debug=F) {
  
  if(debug) browser()
  
  #makes sure file exists
  country_file <- sprintf('%s/tmp/%s_sfu.fst', dir, country)
           
  if (country_file %>% file.exists) {
    
    #read in the long data.tables for the appropriate country
    message('\n********\nstarting process for ', 
            adm_info[ADM0_CODE==country, ADM0_NAME %>% unique], '\n********\n')
    message('reading data from\n')
    message(country_file)
    country_dt <- read_fst(country_file, as.data.table = T) %>% .[, ID := .I]
    states <- unique(country_dt$ADM1_CODE) %>% sort

    #for larger countries, we will process them by state
    stateLoop <- function(i, state_list, debug=F) {
      
      if(debug) browser()
      
      #which state
      state <- state_list[i]
      
      #display progress based on quantiles
      prog_list <- quantile(state_list, p=c(.2, .4, .6, .8)) %>% as.list
      status <- sapply(prog_list, function(x) state_list[which.min(abs(state_list-x))]) %>% 
        sapply(., function(x) abs(state-x) < 1e-5) %>% 
        .[.] %>% 
        names
      
      #if at quantile, give messages
      msg <- length(status)==1
      
      if (msg) message(status, ' complete, starting on state=', state)

      #subset to working adm1
      dt <- country_dt[ADM1_CODE %in% state]
      
      #reshape long draws and key on the pixel ID/draw/year
      if (msg) message('-melting long')
      dt <- melt(dt,
                 measure = patterns("V"),
                 variable.name = "draw",
                 value.name = 'sfu') %>% 
        setkey(., pixel_id, draw, year)

      #merge on the HAP excess PM2.5 values for each ad0
      #note that this expands the dt to have man/woman/child/under5 values for each pixel
      if (msg) message('-merging PM2.5 data')
      dt <- merge(dt, hap_data$pm, by=c('ADM0_CODE', 'year', 'grouping'), all.x=T, allow.cartesian=T)
      
      #calculate the TAP PM2.5/capita values
      if (msg) message('-calculating TAP PM2.5/capita')
      dt[, aap_pm := pop * ihmepm25]
      dt[, hap_pm := pop * sfu * hap_excess_pm25]
      dt[, tap_pm := (aap_pm + hap_pm)]
      dt[, tap_pc := tap_pm / pop]
      dt[pop==0, tap_pc := 0] #must assume no exposure in zero population cells
      
      #calculate the hap proportion and other key hap indicators
      dt[, hap_pct := hap_pm/(aap_pm+hap_pm)]

      #calculate the TAP RR for LRI
      if (msg) message('-merging risk files and expanding dt')
      #first merge on the RR.max to expand the causes for our cellpred dt
      dt <- merge(dt, hap_data$rr_max, by=c('draw', 'age_group_id', 'sex_id'), allow.cartesian=T)
      
      #merge on the mrBRT RR predictions using the nearest spline cutpoint
      if (msg) message('-merging mrBRT predictions in order to estimate RR')
      dt <- hap_data$rr[dt, on=.(draw, age_group_id, sex_id, cause, tap_pc), roll='nearest'] 
      
      #calculate the TAP PAF for each cause
      if (msg) message('-making PAF calculations')
      dt[, tap_paf := (rr-1)/rr]
  
      #calculate the SEV by cause using the HAP PAF and then average across causes
      if (msg) message('-making SEV calculations')
      dt[, hap_sev := (tap_paf*hap_pct/(1-tap_paf*hap_pct))/(rr_max-1)]
      dt[, hap_sev := mean(hap_sev), by=.(pixel_id, year, draw, grouping)]
      
      #population weight using the age/sex distribution per pixel to collapse over causes
      dt[, hap_sev := sum(hap_sev*pop)/sum(pop), by=.(pixel_id, year, draw)]

      #we will output SEVs at the cell_pred/admin_2 level in order to make forecasts
      if (msg) message('-producing SEV cell_pred')
      #child grouping is arbitrary since SEV is now collapsed, we just choose this one to have a single SEV per cellid
      sev_cell_pred <- dt[grouping=='child', .(cell_pred_id, year, draw, hap_sev, area_fraction)] %>% 
        data.table::dcast(area_fraction+cell_pred_id+year~draw, value.var='hap_sev')
      
      if (msg) message('-producing admin_2 SEV')
      
      #fractional aggregation
      sev_admin_2 <- dt[grouping=='child', .(ADM2_CODE, year, draw, hap_sev, pop_total, area_fraction)] %>% 
        setkeyv(., cols=c('ADM2_CODE', 'year', 'draw')) %>% 
        .[, hap_sev := weighted.mean(hap_sev, w=pop_total*area_fraction, na.rm=T), by=key(.)] %>% 
        unique(., by=key(.)) %>% #keep one row per adm2/year/draw now that they are collapsed
        .[, `:=` (area_fraction=NULL, pop_total=NULL)] %>% #no longer needed
        data.table::dcast(ADM2_CODE+year~draw, value.var='hap_sev') #cast back to wide
      
      #add on the identifying info
      sev_admin_2 <- merge(sev_admin_2, adm_links, by='ADM2_CODE')
      
      #drop all cause/groupings except child LRIs, we only use the others to calculate the SEV
      dt <- dt[cause=='lri' & grouping=='under5']

      #first remove any unnecessary cols
      null_cols <- c('ihmepm25', 'hap_excess_pm25', 'aap_pm', 'hap_pm', 'tap_pm',
                     'rr', 'rr_max', 'ID', 'age_group_id', 'sex_id')
      dt <- dt[, (null_cols) := NULL]

      #helper function for aggregating columns from pixel-draws to adms
      aggResults <- function(dt, by_cols, agg_cols, mbg_var, make_ci=T, make_sums=T) {
        
        agg <- copy(dt)
        
        # aggregate to ad2
        if (msg) message('-Aggregating at the level of ', paste(by_cols, collapse=' / '))
        
        #distinguish that count columns should be a weighted sum instead of a weighted mean
        sum_cols <- agg_cols %>% .[. %like% '_c$'] 
        mean_cols <- agg_cols %>% .[!(. %in% sum_cols)] 
        
        #which columns will no longer be relevant after this collapse?
        null_cols <- c('pixel_id', 'area_fraction', 
                       names(agg) %>% .[(. %like% 'ADM')] %>% .[!(. %in% by_cols)])
        
        #check which pop columns were returned (pop_total is produced if using a non-total pop to aggregate)
        #append them to sum_cols, they will be likewise summed in the aggregation step
        pop_cols <- names(agg) %>% .[(. %like% 'pop')]
        sum_cols <- c(sum_cols, pop_cols)
        
        #find out prev of missingness in primary variable
        if (msg) message('-checking missingness in ', mbg_var)
        setkeyv(agg, by_cols)
        agg[, miss := mbg_var %>% get %>% is.na %>% sum %>% as.numeric, by=key(agg)]
        agg[, miss := miss/.N, by=by_cols]
        
        #apply area fraction to weight the summed variables (pops)
        agg[, (sum_cols) := lapply(.SD, function(x) x * area_fraction), .SDcols=sum_cols]

        #aggregate means w/ fractional aggregation
        if (msg) message('-producing means')
        agg[miss<1, paste0(mean_cols) := lapply(.SD,  Hmisc::wtd.mean, weights=pop), .SDcols=mean_cols, by=key(agg)] 

        if (make_sums) agg[, (sum_cols) := lapply(.SD, sum, na.rm=T), .SDcols=sum_cols, by=key(agg)] 
        
        #select which cols to keep
        keep_cols <- c(by_cols, mean_cols)
        if (make_sums) keep_cols <- c(keep_cols, sum_cols)
        
        #collapse and return
        unique(agg, by=key(agg)) %>% 
          .[, keep_cols, with=F] %>% 
          return
        
      }
      
      #agg ad2s by grouping
      agg <- aggResults(dt, 
                        by_cols=c('ADM0_CODE', 'ADM2_CODE', 'year','draw'), 
                        agg_cols=c('hap_pct', 'tap_pc', 'sfu', 'tap_paf'),
                        mbg_var='sfu', 
                        make_sums=T)

      rm(dt) #cleanup to save space
      
      #reformatting results
      #first distinguish between pollution types in long format
      agg <- tidyr::crossing(agg, data.table(type=c('TAP', 'AAP', 'HAP'))) %>% as.data.table
      
      #redefine variables based on particle type
      setnames(agg, c('hap_pct', 'tap_pc', 'tap_paf', 'sfu'), c('share', 'pm_pc', 'paf', 'prev'))
      
      #define share of particulates
      agg[type=='TAP', share := 1]
      agg[type=='AAP', share := 1-share]
      
      #define mass of particulates
      agg[, pm_pc := pm_pc * share]
      
      #define pafs
      agg[, paf := paf * share]
      
      #assume all are exposed to ambient
      agg[type=='AAP', prev := 1]


      list('sev_cell_pred'=sev_cell_pred,
           'sev_admin_2'=sev_admin_2,
           'admin_2'=agg) %>% 
        return
      
    }  

    #define country parallelization
    if(country%in%big_countries) country_cores <- ifelse(country%in%super_big_countries, 1, 3) 
    else country_cores <- 7

    #run all states in country
    out_states <- mclapply(1:length(states), stateLoop, state_list=states, mc.cores=country_cores)

    #combine results
    admin_2 <- lapply(out_states, function(x) x[['admin_2']]) %>% rbindlist
    sev_cell_pred <- lapply(out_states, function(x) x[['sev_cell_pred']]) %>% rbindlist
    sev_admin_2 <- lapply(out_states, function(x) x[['sev_admin_2']]) %>% rbindlist

    #output
    list('admin_2'=admin_2,
         'sev_admin_2'=sev_admin_2,
         'sev_cell_pred'=sev_cell_pred) %>% 
      return
      
  } else {message(country_file[!existance], ' does not exist..skipping!'); return(NULL)}
  
}  

#loop over countries and produce ad0/2 results for TAP
tic('Calculating TAP for each country')
out <- lapply(adm0s, calcTAP, dir=outputdir, debug=F)
toc(log=TRUE)

#bind results
tic('Extracting results from lists')
out_ad2 <- lapply(out, function(x) x[['admin_2']]) %>% 
  .[!sapply(., is.null)] %>% #remove the null tables (missing raster values)
  rbindlist(use.names=T, fill=T)

sev_admin_2 <- lapply(out, function(x) x[['sev_admin_2']]) %>% 
  .[!sapply(., is.null)] %>% #remove the null tables (missing raster values)
  rbindlist(use.names=T, fill=T)

out_sev <- lapply(out, function(x) x[['sev_cell_pred']] %>% as.data.table) %>%
  .[!sapply(., is.null)] %>% #remove the null tables (missing raster values)
  rbindlist(use.names=T, fill=T) 

toc(log=TRUE)

#recombine SEVs into a deduped cell_pred (linking process generates duplicates)
sev_cell_pred <- names(out_sev) %>% .[. %like% 'V'] %>%
  unlink_cell_pred(out_sev, cols=.)

# finish up and save
tic('Saving ad2')
out_path <- file.path(outputdir, paste0(region, '_ad2_draws.fst'))
message('-> finished calculating TAP for ad2 level, now saving as \n...', out_path)
write.fst(out_ad2, path = out_path)
toc(log = TRUE)

tic('Saving SEV')
sev_file <- '<<< FILEPATH REDACTED >>>'
save(sev_cell_pred, file=sev_file)

#***********************************************************************************************************************

# ---SDG PROJECTIONS----------------------------------------------------------------------------------------------------
# Define goals: start by initializing goal object
goals <- add_goal(target_year = 2030, 
                  target = 0.01,
                  target_type = "less",
                  abs_rel = "absolute",
                  pred_type = c("admin_2"))

goals <- add_goal(goal_obj=goals,
                  target_year = 2030,
                  target = 0.05,
                  target_type = "less",
                  abs_rel = "absolute",
                  pred_type = c("admin_2"))

goals <- add_goal(goal_obj=goals,
                  target_year = 2018,
                  target = 0.01,
                  target_type = "less",
                  abs_rel = "absolute",
                  pred_type = c("admin_2"))

goals <- add_goal(goal_obj=goals,
                  target_year = 2018,
                  target = 0.05,
                  target_type = "less",
                  abs_rel = "absolute",
                  pred_type = c("admin_2"))

#make AROC predictions
aroc <- 
  make_aroc(          
    region = region,
    cell_pred = sev_admin_2,
    type = "admin",
    admin_lvl = 2, 
    year_list = year_list,
    uselogit = TRUE,
    weighting_res = 'domain',
    weighting_type = 'exponential',
    pow = 1,
    debug=F
  )

proj <-
  make_proj(
    aroc_draws=aroc,
    cell_pred=sev_admin_2,
    region = region,
    type = "admin",
    admin_lvl=2,
    proj_years = c(2018, seq(2020, 2030, 5)),
    year_list = year_list,
    uselogit = TRUE,
    debug=F
  ) 

probs <- 
  compare_to_target(  
    obj=proj,
    region = region,
    goal_obj = goals,
    year_list = year_list,
    uselogit = TRUE,
    debug=F
  )

#save outputs
sdg_dir <- file.path(outputdir, 'sdg_projections') %T>% dir.create(recursive = T) 
list('goals'=goals,
     'aroc'=aroc,
     'proj'=proj,
     'probs'=probs) %>% 
  saveRDS(., file=file.path(sdg_dir, paste0(region, '_admin_2_outputs.RDS')))

## Attempt to calculate at cell level as well
# Define goals: start by initializing goal object
goals <- add_goal(target_year = 2030, 
                  target = 0.01,
                  target_type = "less",
                  abs_rel = "absolute",
                  pred_type = c("cell"))

#run comparisons
cell_probs <- 
  make_aroc(          
    region = region,
    cell_pred = sev_cell_pred,
    type = "cell",
    admin_lvl = 2, 
    year_list = year_list,
    uselogit = TRUE,
    weighting_res = 'domain',
    weighting_type = 'exponential',
    pow = 1,
    debug=F
  ) %>% 
  make_proj(
    aroc_draws=.,
    cell_pred=sev_cell_pred,
    region = region,
    type = "cell",
    admin_lvl=2,
    proj_years = seq(2020, 2030, 5),
    year_list = year_list,
    uselogit = TRUE,
    debug=F
  ) %>% 
  compare_to_target(  
    obj=.,
    region = region,
    goal_obj = goals,
    year_list = year_list,
    uselogit = TRUE,
    debug=F
  )

#save outputs
saveRDS(cell_probs, file=file.path(sdg_dir, paste0(region, '_sev_cell_probs.RDS')))


toc() # End master timer
#*********************************************************************************************************************** 
