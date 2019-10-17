
rm(list=ls())

# Imports
library(data.table)
library(raster)
library(ggplot2)
library(ggrepel)
library(rgdal)
library(rgeos)

# Source LBD core
source("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")
# Source U5M-specific functions
source("<<<< FILEPATH REDACTED >>>>")
# Source GBD shared functions
source('<<<< FILEPATH REDACTED >>>>')
source("<<<< FILEPATH REDACTED >>>>")

## SET ARGS
yl    <- 2000:2016 # Analysis years
rd    <- '<<<< REDACTED >>>>'
group <- 'under5'
model_start_year <- 2000 # MODEL start year
model_end_year   <- 2017 # MODEL end year
interval_mo <- 12
resume <- TRUE


## SET SAVE LOCATION
save_folder <- '<<<< FILEPATH REDACTED >>>>'
sub_folder <- NULL # Default is a timestamp
prepped_d_file <- paste0(save_folder, '<<<< FILEPATH REDACTED >>>>')


## ALL CODE BELOW THIS POINT WILL BE CONVERTED TO FUNCTIONS ~~~~~~~~~~~~~~~~~~

## Set up save folder
# Generate default folder if save folder is NULL
if( is.null(sub_folder) ) sub_folder <- gsub('-','_',Sys.Date())
# Create folder
save_dir_full <- paste0(save_folder, sub_folder,'/')
dir.create(save_dir_full, recursive = TRUE)


## FORMAT AND LOAD DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## FUNCTION TO PULL ADMIN UNITS
## Function for pulling adm(X) rasters for Stage 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pull_adm_units <- function(save_folder, resume, gauls, compare_raster){
  message("Assembling administrative rasters...")
  out_list <- list('adm0'=NA, 'adm1'=NA, 'adm2'=NA)
  for(lev in 0:2){
    backup_file <- paste0(save_folder,'inputs/adm',lev,'_boundaries.Rdata')
    if((resume==TRUE) & (file.exists(backup_file))){
      # Get admin file from backup
      this_adm_ras <- get(load(backup_file))
    } else {
      # Create adm(X) rasters
      this_adm_ras <- GetAdmin(0,compare_raster,gauls)[['rast']]
    }
    # Ensure that the comparison raster has the same dimensions and extent
    if( !all.equal(dim(this_adm_ras),dim(compare_raster)) ){
      stop(paste0("Admin",lev," raster does not have the correct dimensions."))
    }
    # Add to list
    out_list[[paste0('adm',lev)]] <- this_adm_ras
  }
  # Return all levels in a list
  return(out_list)
}


# Load country lists of interest
lookup_table <- load_adm0_lookup_table()
stage_2_gauls <- lookup_table[Stage !='3', GAUL_CODE]

## IF the prepped data file already exists, load it
## Otherwise, create the data file
if((resume==TRUE) & file.exists(prepped_d_file)){
  # Load the prepped data file  
  d <- get(load(prepped_d_file))
  
} else {
  # Create the prepped data file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load in D and Q rasters
  in_dir <- sprintf('<<<< FILEPATH REDACTED >>>>',group,rd)
  D <- brick(sprintf('<<<< FILEPATH REDACTED >>>>',
                     in_dir, group, model_start_year, model_end_year))
  Q <- brick(sprintf('<<<< FILEPATH REDACTED >>>>',
                     in_dir, group, model_start_year, model_end_year))
  # Ensure that these bricks have the same extent
  if(extent(D)!=extent(Q)){
    stop("The input raster bricks do not have the same dimensions - please fix.")
  }
  # Load population
  pop <- load_and_crop_covariates_annual(
    covs           = 'worldpop',
    measures       = 'a0004t',
    simple_polygon = D,
    start_year     = model_start_year,
    end_year       = model_end_year,
    interval_mo    = as.numeric(12),
    agebin=1
  )[[1]]
  
  # Load in country data (admin0)
  # This function lives in team repo: 
  adm_units <- pull_adm_units(
    save_folder = save_folder,
    resume = TRUE,
    gauls = stage_2_gauls,
    compare_raster = D
  )
  countries <- adm_units[['adm0']]
  
  # Subset only to years in the year list
  year_idx <- which(model_start_year:model_end_year %in% yl)
  D_sub         <- D[[year_idx]]
  Q_sub         <- Q[[year_idx]]
  pop_sub       <- pop[[year_idx]]
  countries_sub <- countries[[year_idx]]
  # Remove the larger raster objects
  rm(list=c('D','Q','pop','countries'))
  
  # Vectorize all data into a data.table
  d  <- data.table(d         = as.vector(D_sub),
                   q         = as.vector(Q_sub),
                   pop       = as.vector(pop_sub),
                   year      = rep(yl, each = dim(D_sub)[1] * dim(D_sub)[2]),
                   GAUL_CODE = as.vector(countries_sub)
  )
  d  <- na.omit(d)
  
  # Pull in mapping for country names
  loc <- load_adm0_lookup_table()
  loc <- loc[, .(GAUL_CODE, location_name, spr_reg_nm, loc_id)]
  setnames(loc,'location_name','country_name')
  # Merge on country names
  d <- merge(
    x     = d,
    y     = loc,
    by    = c('GAUL_CODE')
  )
  # Set under-5 mortality rate as deaths/population
  d[, u5mr := d/pop]
  
  # Save prepped data file
  if(!file.exists(prepped_d_file)) save(d, file=prepped_d_file)
}



## DETERMINE GINI COEFFICIENT GLOBALLY AND FOR EACH COUNTRY ~~~~~~~~~~~~~~~~~~~~

## calc_gini () ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' Given cumulative proportions of the population and some indicator of
#'  interest, calculate the Gini coefficient for that group. This calculation
#'  comes from the World Bank Poverty Manual, p. 97:
#'  see http://siteresources.worldbank.org/PGLP/Resources/PMch6.pdf
#'
#' @param sorted_pop A vector of all population, sorted from lowest to highest
#'   on the amount of <var> per person within each group
#' @param sorted_var A vector of the amount of <var> attributable to each group,
#'   sorted in the same order as sorted_pop
calc_gini <- function(sorted_pop, sorted_var){
  # Ensure that the vectors are of the same length
  if(length(sorted_pop) != length(sorted_var)){
    stop(paste("When calculating the gini coefficient, the sorted_pop and",
               " sorted_var vectors must be of the same length."))
  }
  # Get the total amount of population and the given variable
  total_var <- sum(sorted_var)
  total_pop <- sum(sorted_pop)
  # Normalize population and the variable of interest, then get the cumulative
  #  sum of both within the sorted vectors
  cumul_var <- cumsum(sorted_var) / total_var
  cumul_pop <- cumsum(sorted_pop) / total_pop
  # Start at (0,0)
  cumul_var <- c(0, cumul_var)
  cumul_pop <- c(0, cumul_pop)
  # Get the contribution of the Gini coefficient from each value 
  # Calculate the Gini formula recursively
  g <- 1
  for(ii in 2:length(cumul_var)){
    g <- g - (cumul_pop[ii]-cumul_pop[ii-1])*(cumul_var[ii]+cumul_var[ii-1])
  }
  # Return the calculated Gini coefficient
  return(round(g,6))
}

## MERGE A COVARIATE ONTO THE DATA BY COUNTRY NAME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merge_gbd_covariate <- function(in_data, covariate_id, loc_field='country_name'){
  # Get covariate estimates
  cov_est <- get_covariate_estimates(
    covariate_id = covariate_id,
    sex_id = 3
  )
  # Get covariate name
  cov_name <- cov_est[1,covariate_name_short]
  # Set column names to match input data
  setnames(cov_est, 'location_name', loc_field)
  setnames(cov_est, 'year_id', 'year')
  ## MERGE THE COVARIATE ON  
  if(length(unique(cov_est[,year])) > 1){
    # If there are multiple years, merge onto the dataset by country name and year
    cov_est <- cov_est[, .(country_name, year, mean_value)]
    merged_data <- merge(
      x     = in_data,
      y     = cov_est,
      by    = c('country_name','year'),
      all.x = TRUE
    )
  } else {
    # If there is only one year, merge onto the dataset by country name
    cov_est <- cov_est[, .(country_name, mean_value)]
    merged_data <- merge(
      x     = in_data,
      y     = cov_est,
      by    = c('country_name'),
      all.x = TRUE
    )
  }
  # Set the covariate name to something descriptive
  setnames(merged_data, 'mean_value', cov_name)
  # Return dataset
  return(merged_data)
}


## Determine Gini coefficient of mortality globally for each year ~~~~~~~~~~~~~~
gini_global <- copy(d)
gini_global <- gini_global[order(year, u5mr)]
gini_global <- gini_global[, .(gini = calc_gini(pop, d)), by=year]
# Create field normalized to 2000
gini_global[, merge_var := 1]
gini_global <- merge(
  x        = gini_global,
  y        = gini_global[year==min(yl),.(gini, merge_var)],
  by       = c('merge_var'),
  suffixes = c('','_fy'),
  all.x    = TRUE
)
gini_global[, gini_normalized := gini/gini_fy]
gini_global[, c('merge_var','gini_fy') := NULL]
  
## Determine Gini coefficient of mortality by super region and year ~~~~~~~~~~~~
gini_sr <- copy(d)
gini_sr <- gini_sr[order(spr_reg_nm, year, u5mr)]
gini_sr <- gini_sr[, .(gini = calc_gini(pop,d)), by=c('spr_reg_nm','year')]
# Create field normalized to 2000
gini_sr <- merge(
  x        = gini_sr,
  y        = gini_sr[year==min(yl), .(gini, spr_reg_nm)],
  by       = c('spr_reg_nm'),
  suffixes = c('','_fy'),
  all.x    = TRUE
)
gini_sr[, gini_normalized := gini/gini_fy]
gini_sr[, gini_fy := NULL]

## Determine Gini coefficient of mortality by country and year ~~~~~~~~~~~~~~~~~
gini_country <- copy(d)
gini_country <- gini_country[order(country_name, year, u5mr)]
gini_country <- gini_country[, .(gini = calc_gini(pop, d)),
                             by=c('country_name','year','spr_reg_nm')]
# Create field normalized to 2000
gini_country <- merge(
  x        = gini_country,
  y        = gini_country[year==min(yl), .(gini, country_name)],
  by       = c('country_name'),
  suffixes = c('','_fy'),
  all.x    = TRUE
)
gini_country[, gini_normalized := gini/gini_fy]
gini_country[, gini_fy := NULL]

gini_country <- merge_gbd_covariate(gini_country, covariate_id=881) # SDI
gini_country <- merge_gbd_covariate(gini_country, covariate_id=57) # LDI
gini_country <- merge_gbd_covariate(gini_country, covariate_id=463) # Maternal edu
# Shorten a long variable name
setnames(gini_country, 'maternal_educ_yrs_pc','mat_edu')

# Check how countries have changed over time
measure_vars <- c('gini','sdi','LDI_pc','mat_edu')
gini_country_wide <- dcast.data.table(
  data      = gini_country,
  formula   = spr_reg_nm + country_name ~ year,
  value.var = measure_vars
)
# Get ARC for variables, 2000-2016
for(nm in measure_vars){
  var_start <- paste0(nm,'_',yl[1])
  var_end   <- paste0(nm,'_',yl[length(yl)])
  newvar <- paste0('arc_',nm)
  gini_country_wide[, (newvar) := log(get(var_end)/get(var_start))/(yl[length(yl)]-yl[1])]
}



## CREATE GLOBAL GINI PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d_global <- copy(d)
# sort it and get cumulative percents by country year
d_global <- d_global[order(year, u5mr)]
d_global[, global_pop := sum(pop), by=year]
d_global[, cumul_prop_pops := cumsum(pop)/global_pop, by=year]
# Bin into %1 proportions of total population for easier graphing
d_global[, prop_pop_rounded :=round(cumul_prop_pops,2)]
d_global_agg <- d_global[, .(d = sum(d)),
                             by=c('year','prop_pop_rounded')]
d_global_agg[, global_deaths     := sum(d), by=year]
d_global_agg[, cumul_prop_deaths := cumsum(d)/global_deaths, by=year]
d_global_agg[, year_factor       := as.factor(year)]

## PLOT GLOBAL GINI COEFFICIENT
pdf(paste0(save_dir_full,"<<<< FILEPATH REDACTED >>>>"),
    width=6,
    height=6,
    onefile=TRUE)
fig <- ggplot(d_global_agg, aes(x=prop_pop_rounded, y=cumul_prop_deaths)) +
  geom_line(aes(color=year_factor)) +
  geom_abline(slope=1, intercept=0, color='#AAAAAA') +
  xlim(c(0,1)) + ylim(c(0,1)) +
  labs(title='Plots of mortality inequality across Stage 2',
       subtitle=paste0("Stage 2 Gini: ",round(gini_global[year==min(yl),gini],3),
                      " (",min(yl),") , ",round(gini_global[year==max(yl),gini],3),
                      " (",max(yl),")"),
       x="Proportion of Total Population in Stage 2",
       y="Proportion of Total Deaths in Stage 2",
       color="Year") +
  theme_minimal()
print(fig)
dev.off()




## CREATE GINI PLOTS BY COUNTRY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sort it and get cumulative percents by country year
d <- d[order(GAUL_CODE, year, u5mr)]
d[, country_pop := sum(pop), by=c('GAUL_CODE','year')]
d[, cumul_prop_pops := cumsum(pop)/country_pop, by=c('GAUL_CODE','year')]
# Bin into %1 proportions of total population for easier graphing
d[, prop_pop_rounded :=round(cumul_prop_pops,2)]
d_agg <- d[, .(d = sum(d)),
           by=c('GAUL_CODE','country_name','year','prop_pop_rounded')]
d_agg[, country_deaths    := sum(d), by=c('GAUL_CODE','year')]
d_agg[, cumul_prop_deaths := cumsum(d)/country_deaths, by=c('GAUL_CODE','year')]
d_agg[, year_factor       := as.factor(year)]
# Merge country-level Gini data onto dataset
gini_for_merge <- gini_country_wide[, .(country_name, gini_2016)]
d_agg <- merge(
  x     = d_agg,
  y     = gini_for_merge,
  by    = c('country_name'),
  all.x = TRUE
)
d_agg[, country_title := paste0(country_name,"\n(2016 Gini: ",round(gini_2016,3),")")]

## CREATE THE FIGURE BY COUNTRY CHUNKS
all_countries <- unique(d_agg[,country_name])
# Split into chunks of 16
country_chunks <- split(all_countries, ceiling(seq_along(all_countries)/16))

pdf(paste0(save_dir_full,"<<<< FILEPATH REDACTED >>>>"),
    width=10,
    height=10,
    onefile=TRUE)

d_agg_2 <- d_agg[year %in% c(min(yl), max(yl)),]

for(chunk in country_chunks){
  d_sub <- d_agg_2[country_name %in% chunk,]
  # PLOT SUBSET
  fig <- ggplot(d_sub, aes(x=prop_pop_rounded, y=cumul_prop_deaths)) +
    geom_line(aes(color=year_factor)) +
    facet_wrap('country_title', ncol=4) +
    geom_abline(slope=1, intercept=0, color='#AAAAAA') +
    xlim(c(0,1)) + ylim(c(0,1)) +
    labs(title='Plots of mortality inequality within countries',
         x="Proportion of Total Population in Country",
         y="Proportion of Total Deaths in Country",
         color="Year") +
    theme_minimal()
  print(fig)
}

dev.off()


## SHOW NATIONAL INEQUALITY AND ARC BY NATIONAL COVARIATES ~~~~~~~~~~~~~~~~~~~~~

# (1) Show national inequality by 2016 SDI
fig1 <- ggplot(gini_country_wide,
               aes(x=sdi_2016, y=gini_2016, color=spr_reg_nm)) +
  geom_point() + 
  labs(
    title='Mortality Gini coefficient in 2016 against SDI in 2016',
    x='Socio-demographic index (2016)',
    y='Mortality Gini coefficient (2016)',
    color='GBD Super Region') + 
  theme_minimal() + 
  theme(legend.position='bottom')

# (2) Show ARC in inequality by ARC in SDI
fig2 <- ggplot(gini_country_wide,
               aes(x=arc_sdi, y=arc_gini, color=spr_reg_nm)) +
  geom_point() + 
  labs(
    title='Annualized rate of change in Mortality Gini against ARC in SDI',
    subtitle='Period considered: 2000-2016',
    x='ARC of Socio-demographic index (2000-2016)',
    y='ARC of mortality Gini coefficient (2000-2016)',
    color='GBD Super Region') + 
  theme_minimal() + 
  theme(legend.position='bottom')

# (3) Show ARC in inequality by ARC in LDI
fig3 <- ggplot(gini_country_wide,
               aes(x=arc_LDI_pi, y=arc_gini, color=spr_reg_nm)) +
  geom_point() + 
  labs(
    title='Annualized rate of change in Mortality Gini against ARC in LDI per capita',
    subtitle='Period considered: 2000-2016',
    x='ARC of Lag-distributed income per capita (2000-2016)',
    y='ARC of mortality Gini coefficient (2000-2016)',
    color='GBD Super Region') + 
  theme_minimal() + 
  theme(legend.position='bottom')

# (4) Show ARC in inequality by ARC in maternal education
fig4 <- ggplot(gini_country_wide,
               aes(x=mat_edu_2000, y=arc_gini, color=spr_reg_nm)) +
  geom_point() + 
  labs(
    title='ARC in Mortality Gini against ARC in Maternal Education (2000)',
    subtitle='Period considered: 2000-2016',
    x='ARC of Maternal Education per capita (2000-2016)',
    y='ARC of mortality Gini coefficient (2000-2016)',
    color='GBD Super Region') + 
  theme_minimal() + 
  theme(legend.position='bottom')



## SHOW CHANGE IN GINI COEFFICIENT OVER TIME, BY COUNTRY, SR, AND WORLDWIDE ~~~

fig_change_sr <- ggplot() + 
  geom_line(data=gini_sr, aes(x=year, y=gini, color=spr_reg_nm), size=1) + 
  geom_line(data=gini_global, aes(x=year, y=gini), size=1.5, color='#000000') + 
  labs(
    title    = paste('Changes in mortality inequality across country groups,',
                     min(yl),'to',max(yl)),
    subtitle = paste('Overall mortality inequality shown for all LMICs (in black)',
                     'and by GBD Super Region'),
    x        = 'Year',
    y        = 'Gini coefficient for mortality',
    color    = 'GBD Super Region'
  ) +
  theme_minimal() + 
  theme(legend.position='bottom')

fig_change_sr_norm <- ggplot() + 
  geom_line(data=gini_sr, aes(x=year, y=gini_normalized, color=spr_reg_nm), size=1) + 
  geom_line(data=gini_global, aes(x=year, y=gini), size=1.5, color='#000000') + 
  labs(
    title    = paste('Changes in mortality inequality across country groups,',
                     min(yl),'to',max(yl)),
    subtitle = paste('Overall mortality inequality shown for all LMICs (in black)',
                     'and by GBD Super Region'),
    x        = 'Year',
    y        = 'Gini coefficient for mortality',
    color    = 'GBD Super Region'
  ) +
  theme_minimal() + 
  theme(legend.position='bottom')

  
fig_change_sr


