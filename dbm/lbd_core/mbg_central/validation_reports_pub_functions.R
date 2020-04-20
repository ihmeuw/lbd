plot_violins <- function(indicator, indicator_proper_name, indicator_group, run_date, this_gaul, shapefile_verison = 'current') {
  
  print(this_gaul)
  
  in_dir  <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date)
  default_rf_path <- paste0(in_dir, '/', indicator, '_rf.csv')
  all_rfs <- fread(default_rf_path)
  
  ## Define path to .tif of results and raked results, load in raster bricks.
  default_raked_results_path <- brick(paste0(in_dir, '/', indicator, '_mean_raked_raster.tif'))
  results_raked <- brick(default_raked_results_path)
  default_results_path <- brick(paste0(in_dir, '/', indicator, '_mean_raster.tif'))
  results <- brick(default_results_path)
  total_periods <- length(names(results))
  regions <- get_output_regions(in_dir)
  for(reg in regions) {
    if(this_gaul %in% get_adm0_codes(reg, shapefile_version = shapefile_version)) this_region <- reg
  }
  input_data <- list.files(in_dir, pattern = paste0("input_data[a-zA-Z0-9_]*", this_region), full.names = T) %>% fread
  
  gaul_to_loc_id <- get_location_code_mapping()
  gauls <- this_gaul
  this_country <- gaul_to_loc_id[GAUL_CODE %in% gauls, ihme_lc_id]
  
  if(indicator != 'edu_mean') input_data <- input_data[, rate := get(indicator) / N]
  if(indicator == 'edu_mean') input_data <- input_data[, rate := get(indicator)]
  country_input <- input_data[country %in% this_country, ]
  country_input <- subset(country_input, year >= 1998)
  names(country_input)[names(country_input) == "year"] = "original_year"
  country_input <- country_input[original_year >= 1998 & original_year < 2003, year := 2000]
  country_input <- country_input[original_year >= 2003 & original_year < 2008, year := 2005]
  country_input <- country_input[original_year >= 2008 & original_year < 2013, year := 2010]
  country_input <- country_input[original_year >= 2013 & original_year < 2018, year := 2015]
  country_input_dt <- country_input
  
  if(nrow(country_input_dt) == 0){
    return(NULL)
  }
  
  coordinates(country_input) = ~longitude+latitude
  
  ## Calculate variogram for predictions
  ## Extract preds for each year in country data
  default_raked_results_path <- paste0(in_dir, '/', indicator, '_mean_raked_raster.tif')
  results_raked <- brick(default_raked_results_path)
  default_results_path <- paste0(in_dir, '/', indicator, '_mean_raster.tif')
  results <- brick(default_results_path)
  extract_year_preds <- function(country_year) {
    period <- country_year - 2000 + 1
    input_data_year <- country_input_dt[year == country_year, ]
    preds_year <- results[[period]]
    preds_at_points <- raster::extract(preds_year, input_data_year[, c('longitude','latitude'), with = F])
    input_data_year <- input_data_year[, pred := preds_at_points]
    return(input_data_year)
  }
  country_input_wpreds <- rbindlist(lapply(unique(country_input_dt[,year]), extract_year_preds))
  country_input_wpreds <- country_input_wpreds[!is.na(pred),]
  country_input_wpreds <- country_input_wpreds[, data := rate]
  violin_vars <- c('data', 'pred')
  country_input_wpreds <- country_input_wpreds[, c('latitude','longitude','weight', violin_vars ,'year'), with=F]
  country_input_wpreds <- melt(country_input_wpreds, id.vars = c('longitude','latitude','weight','year'), measure.vars = violin_vars)
  country_input_wpreds <- country_input_wpreds[, country := this_country]
  if(file.exists(paste0(results_dir,'/', indicator, '_rf.csv'))) {
    rfs <- fread(paste0(results_dir,'/', indicator, '_rf.csv'))
    country_input_wpreds <- country_input_wpreds[, name := as.integer(this_gaul)]
    country_input_wpreds <- merge(country_input_wpreds, rfs, by=c('name','year'))
  }
  violin_countries <- this_country
  loc_names <- get_location_code_mapping()
  convert_to_iso3 <- function(x) {
    if(nchar(x) > 3) x <- loc_names[loc_name==x, ihme_lc_id]
    return(x)
  }
  
  
  #make seperate graphs per year rather than facet wrapping
  graphs = list()
  iter = 1
  country_input_wpreds[variable == 'pred', variable := 'prediction']
  for(yyy in unique(country_input_wpreds[,year])){
    graph_dat = country_input_wpreds[year == yyy,]
    violin <- ggplot(data=graph_dat) +
      geom_violin(aes(x = variable, y = value, fill = country, weight = weight))
    if(file.exists(paste0(results_dir,'/', indicator, '_rf.csv'))) {
      violin <- violin + geom_hline(aes(yintercept=rake_to_mean))
    }
    violin <- violin +
      theme(axis.text.x = element_text(face="bold", size=10, angle=90)) +
      theme(legend.position="none") +
      xlab("") +
      ylab(paste0('Mean: ', indicator_proper_name)) +
      facet_wrap(~country + year) +theme_bw() + theme(legend.position="none")
    
    graphs[[iter]] = violin
    iter = iter+1
    
  }
  ## Return ggs for violin, vario, and pacf model object to plot() later.
  return(graphs)
  
}
