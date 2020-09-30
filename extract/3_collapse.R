# ----HEADER------------------------------------------------------------------------------------------------------------
# Author:#Redacted
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

#Redacted

#load packages
#Redacted
pacman::p_load(data.table, dplyr, magrittr, feather, fst, googledrive, readxl, gargle, ellipsis, sf) 
package_list <- paste0(core_repo, '#Redacted/package_list.csv') %>% fread %>% t %>% c

#capture date
today <- Sys.Date() %>% gsub("-", "_", .)
#options
cores <- 10
this.family='cooking'
modeling_shapefile_version <- "2019_09_10"
#manual_date <- "2019_12_23" #set this value to use a manually specified extract date
latest_date <- T #set to TRUE in order to disregard manual date and automatically pull the latest value
save_intermediate <- F
run_collapse <- T #set to TRUE if you have new data and want to recollapse
run_resample <- T #set to TRUE if you have new data and want to rerun polygon resampling
save_diagnostic <- F #set to TRUE to save the problematic survey diagnostic
new_vetting <- F #set to TRUE to refresh the vetting diagnostic
redownload <- T
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
#Redacted

###Output###
#Redacted
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
#general functions#
central.function.dir <- file.path(h_root, "_code/_lib/functions")
# this pulls the general misc helper functions
file.path(central.function.dir, "misc.R") %>% source
#hap functions#
hap.function.dir <- file.path(h_root, '_code/lbd/hap/extract/functions')
#this pulls hap collapse helper functions
file.path(hap.function.dir, '/collapse_fx.R') %>% source
#shared functions#
gbd.shared.function.dir <- file.path(j_root,  "#Redacted")
file.path(gbd.shared.function.dir, 'get_location_metadata.R') %>% source
file.path(gbd.shared.function.dir, 'get_ids.R') %>% source
file.path(gbd.shared.function.dir, 'get_covariate_estimates.R') %>% source

lbd.shared.function.dir <- file.path(h_root, "_code/lbd/hap/mbg_central")
file.path(lbd.shared.function.dir, 'setup.R') %>% source
mbg_setup(repo=lbd.shared.function.dir, package_list=package_list) #load mbg functions
#***********************************************************************************************************************

# ---COLLAPSE-----------------------------------------------------------------------------------------------------------
#read in codebook
#Redacted
codebook <- file.path(raw.dir, 'hap.xlsx') %>% read_xlsx(., sheet='codebook') %>% as.data.table
stages <- file.path(j_root, '#Redacted/stage_master_list.csv') %>% fread #read info about stages

#automatically pull latest date if manual date not provided
# get input version from most recently modified data file
if (latest_date) { 

  #pull latest filedate
  file_date <- list.files(data.dir, pattern = '.fst') %>% 
    sort(., decreasing=T) %>% 
    .[1] %>% #pull only the latest file to save
    str_sub(., start=-14, end=-5) #extract from the location of the filedate (end of file -'.fst')
  
} else file_date <- manual_date
  
#loop over points and polygons to collapse
collapseData <- function(this.family,
                         census.file=NULL,
                         out.temp=NULL,
                         subcountry=NULL,
                         debug=F) {
  
  message("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

  #ipums files are done in parallel due to size
  if (!is.null(census.file)) {
    message("CENSUS FILE=", census.file)
    census <- T
    # Load data
    raw <- read.fst(census.file, as.data.table=T)
    dt <- initialClean(raw, var.fam=this.family)
    
  } else {
    census <- F
    # Load data
    raw <- list.files(data.dir, pattern=paste0(file_date, '.fst'), full.names = T) %>% 
      lapply(., read.fst, as.data.table=T) %>% 
      rbindlist(use.names = T, fill=T)

    dt <- list(
      initialClean(raw[!is.na(lat) & !is.na(long)], var.fam=this.family),
      initialClean(raw[is.na(lat) & is.na(long)], var.fam=this.family)
    ) %>% rbindlist
    
  } 
  
  message("Loading data...[census=", census, "]")

  #loop over various families of indicators
  message(paste('->Processing:', this.family))
  #### Subset & Shape Data ####
  #subset to countries, generally used interactively in order to see why surveys are being dropped
  if (!(is.null(subcountry))) {
    dt <- dt[ihme_loc_id==subcountry]
    raw <- raw[iso3==subcountry]
  }
  
  #launch browser to debug interactively
  if (debug) browser()

    
  #output an intermediate file prior to collapse/indicator definition for preliminary analysis
  if (!is.null(out.temp)) {
    message('----->Save raw data to temp folder')
    
    
    out.dir <- file.path(out.temp, this.family) 
    if (!out.dir %>% dir.exists) dir.create(out.dir, recursive=T) #create dir if missing
    
    #build filename and then save a feather
    paste0(out.dir, '/uncollapsed_',
           ifelse(census, tools::file_path_sans_ext(basename(census.file)),
                  ifelse(point, 'points', 'poly')),
           '.feather') %>% write_feather(raw, path=.)
    
    paste0("saved intermediate files to:", out.temp) %>% return #end process here if saving int files
    
  } else {  

    #define the indicators based on the intermediate variables youve extracted  
    dt <- defIndicator(dt, var.fam=this.family, definitions=def.file, debug=F)
    
    #### Address Missingness ####
    message("\nBegin Addressing Missingness...")
    
    # ID clusters with more than 20% weighted missingness
    missing.vars <- idMissing(dt, this.var="cooking_fuel_dirty", criteria=.2, wt.var='hh_size')
    
    #ID cluster_ids with missing hhweight
    #decided to use 20% unweighted criteria instead of 0 tolerance
    missing.wts <- idMissing(dt, this.var="hhweight", criteria=.2, wt.var=NA)
    
    #ID points with hhweight|hh_size<=0 (invalid)
    #drop clusters with more than 20% invalid, then drop invalid rows 
    invalid.wts <- idMissing(dt, this.var="hhweight", criteria=.2, wt.var=NA, check.threshold = T, threshold=0)
    invalid.sizes <- idMissing(dt, this.var="hh_size", criteria=.2, wt.var=NA, check.threshold = T, threshold=0)
    
    #ID missing hh sizes, then crosswalk values
    missing.sizes <- idMissing(dt, this.var="hh_size", criteria=.2, wt.var=NA)
    
    #also print the # of hh sizes that are missing (rowwise):
    message('There are #', nrow(dt[is.na(hh_size)]), '(',
            round(nrow(dt[is.na(hh_size)])/nrow(dt)*100), '%) rows missing hh_size') 
    
    #output diagnostics regarding invalid clusters
    remove.clusters <- c(missing.vars, 
                         missing.wts,
                         invalid.wts,
                         invalid.sizes) %>% unique
    
    #remove these clusters and proceed
    message('dropping ', length(remove.clusters), 
            ' clusters based on variable missingness/invalidity above cluster-level criteria thresholds')
    dt <- dt[!(cluster_id %in% remove.clusters)]
    
    #remove invalid rows that were insufficient in number to trigger criteria thresholds
    message('dropping additional ', dt[(hhweight<=0)] %>% nrow, 
            ' rows based on hhweight missingness/invalidity below cluster-level criteria thresholds')
    dt <- dt[hhweight>0] #drop invalid rows as well
    message('dropping additional ', dt[(hh_size<=0)] %>% nrow, 
            ' rows based on hhsize invalidity below cluster-level criteria thresholds')
    dt <- dt[hh_size>0] #drop invalid rows as well
    
    #subset to years that are >= 2000 as we dont model before this time period
    message('\nCreate column with median year of each cluster. Subset to >2000')
    dt[, year_median := median(int_year, na.rm=T) %>% floor, by=cluster_id]
    dt[, year_median := weighted.mean(year_median, w=hhweight*hh_size) %>% floor, by=cluster_id]
    dt <- dt[year_median>=2000]
  
    #### Aggregate Data ####
    # Aggregate indicator to cluster level
    message("\nBegin Collapsing Variables")
    agg.dt <- aggIndicator(dt, var.fam=this.family, debug=F) #list of variables to aggregate
    message("->Complete!")

    # Standardize the year variable for each survey using the weighted mean of the NID
    # Weight by sum of sample weights
    message("\nStandardizing Year Variable")
    agg.dt[, year := weighted.mean(year_median, w=sum_of_sample_weights) %>% floor, by=nid]
  
    # Skip the rest of the process if no rows of data are left
    if (nrow(dt) == 0) {message('no data left to return!'); return(NULL)}
    else return(agg.dt)
    
  }
  
}
  
#populate vector of IPUMS filepaths
ipums.files = list.files(census.dir, pattern='*.fst', full.names = T)

if (run_collapse) {
  
  #Run fx for each point/poly
  cooking <- collapseData('cooking')
  
  #Run fx for each census file
  cooking.census <- mcmapply(collapseData, census.file=ipums.files, this.family='cooking', SIMPLIFY=F, mc.cores=1) %>% 
    rbindlist
    
  #combine all and redefine the row ID
  cooking <- list(cooking, cooking.census) %>% 
    rbindlist %>% 
    .[, row_id := .I] %>% 
    setkey(., row_id)
  
  #cleanup
  rm(cooking.census)

  #save poly and point collapses
  paste0(out.dir, "/", "collapsed_data_", this.family, ".fst") %>%
    write.fst(cooking, path=.)

  #combine diagnostics and cleanup intermediate files
  col.diagnostic <- 
    collapseCleanup(this.family, codebook=codebook, test.vars=c('cooking_fuel_dirty', 'hh_size'), cleanup=T, debug=F)
  
} else cooking <- paste0(out.dir, "/", "collapsed_data_", this.family, ".fst") %>% read.fst(as.data.table=T)
#***********************************************************************************************************************

# ---DATA EXCLUSION-----------------------------------------------------------------------------------------------------
#exclude datapoints based on HAP vetting
setwd(doc.dir)
#Redacted
vetting <- file.path(doc.dir, 'HAP Tracking Sheet.xlsx') %>% read_xlsx(sheet='1. Vetting', skip=1) %>% as.data.table
excluded_nids <- vetting[`HAP Vetting Status`=='Excluded', nid] %>% unique #get list of excluded points
message('Excluding the following NIDs:')
print(excluded_nids)
cooking <- cooking[!(nid %in% excluded_nids)] #remove any excluded datapoints
#***********************************************************************************************************************

# ---RESAMPLE-----------------------------------------------------------------------------------------------------------
#prep for resampling
vars <- names(cooking) %>% .[. %like% 'cooking']
vars <- c('cooking_fuel_solid', 'cooking_fuel_dirty', 'cooking_fuel_kerosene')

#convert to count space
cooking[, (vars) := lapply(.SD, function(x, count.var) {x*count.var}, count.var=N), .SDcols=vars]

dt <- cooking[iso3 %in% unique(stages[Stage %in% c('1', '2a', '2b'), iso3])] %>% 
  setnames(.,  c('lat', 'long'), c('latitude', 'longitude')) %>% #mbg formatting requirement
  .[!(shapefile %in% shapefile_issues)] 

#resample the polygons using internal fx
if (run_resample) {
  
  pt <- dt[polygon==T] %>% #only pass the poly rows to the fx, pts are dropping
    resample_polygons(data = .,
                      cores = 20,
                      indic = vars,
                      density = 0.001,
                      gaul_list = lapply(unique(dt$iso3) %>% tolower, get_adm0_codes) %>% unlist %>% unique) %>% 
    .[, pseudocluster := NULL] #redundant

} else stop('run_resample=FALSE, cannot save model data')

#combine all points
dt <- list(dt[polygon==F], pt) %>% 
  rbindlist(use.names=T, fill=T) %>% 
  .[polygon==F, weight := 1] #weights are produced by the resample polygons fx

for (iso in unique(dt$ihme_loc_id) %>% sort) {
  
  dropped.nids <- unique(cooking[ihme_loc_id==iso, nid]) %>% 
    .[!(. %in% unique(dt[ihme_loc_id==iso, nid]))]
  
  if (length(dropped.nids)>0) { 
    
    message('\n', iso, '...NIDs dropped by resample:\n') 
    cat(dropped.nids, sep=' | ')
    #note that WHO WHS is a known issue that currently cannot be geolocated
    codebook[nid %in% dropped.nids & survey_name!='WHO_WHS', .(nid, year_start, survey_name)] %>% print
    
  }
    
}

#redefine row ID after resampling
dt[, row_id := .I]
setkey(dt, row_id)

#save resampled data
paste0(model.dir, "/", "resampled_data_", this.family, ".fst") %>%
  write.fst(dt, path=.)

#prep for MDG
setnames(dt,
         c('iso3'),
         c('country'))

#Redacted

#save into MDG dir
#save each one for modelling in binary/ordinal space
#Redacted
#***********************************************************************************************************************************