#Extract Covariates
#About: Reads in a data table (with latitude and longitude) and a list of rasters/bricks/stacks and extracts the values at the points
#df: a data table. Should have columns named latitude and longitude at a miniumum
#covariate list: a list object of rasters-like formats (e.g. layer, brick, stack) or just one of that raster-like object
#unlike some of the other MBG covariate functions, this extract_covariates function should be able to treat
#raster bricks as a de facto list-- only if 1 brick is passed the function doesn't automatically assume time varying)
#instead, that is implied from the layer names within the brick.
#For best use, provide a list OR a single brick when trying to analyze multiple sets of covariates
#reconcile_timevarying: If true, the functio n enforces period matching. For example, if evi is in 4 periods, there will be one EVI column
#where the values of the column are period/time/year specific. Period in this case is usually identified by
# a .#### at the end of the layer name (usually within a brick). As of 11/2, this is .1, .2, .3 ,.4 but it should
# work for something like .1990, .2000, etc. more testing is needed
#period_var: is there a column in the passed data frame designating the time period that matches the rasters
#return_only_results: return only the covariate columns? (better for cbinding)
extract_covariates = function(df, covariate_list, id_col = NULL, reconcile_timevarying = T, period_var = NULL, return_only_results = F, centre_scale = F, period_map = NULL){
  #Load libraryd packages. Should be redundent
  
  library('raster')
  library('data.table')
  df = copy(df) #data table does weird scoping things. Do this just in case
  
  #create an id column
  if(is.null(id_col)){
    df[,rrrr := 1:nrow(df)]
    id_col = 'rrrr'
  }
  
  # fill in standard 5-year period map if it's missing
  if(is.null(period_map)) period_map <- make_period_map(modeling_periods = c(2000,2005,2010,2015))
  
  #the periods thing is kind of dumb. for now, confirm that the year variable is in four periods to match the raster bricks
  # if(length(unique(df[,year]))!=4){
  #   stop("Your input data frame's year column is not in the 4 period mode to match the 4 periods of likely raster bricks \n Please fix and rerun")
  # }
  
  
  #if covariate list is not actually a list, convert it to a list
  if(class(covariate_list) != 'list'){
    
    #if a brick or a stack, extract the basename
    if(class(covariate_list) %in% c('RasterBrick', 'RasterStack', "RasterLayer")){
      #get the list of brick names
      #find the prefixes, remove the .### at the end of the brick
      brick_names = unique(gsub("(\\.[0-9]+)$","",names(covariate_list)))
      
      #convert the brick into a list the rest of the function expects (e.g. time varying vs. single)
      #assuming I've done this properly, this will convert the raster brick into a list where items in the list are seperated by the prefix
      #(assumes a .### at the end of the name in the raster brick object refers to the period)
      
      covariate_list = setNames(lapply(brick_names, function(x) covariate_list[[grep(x, names(covariate_list), value =T )]]),brick_names)
      
      
    }
    
  }
  
  #rename the raster bricks
  #borrowed from some of the covariate functions
  #converts the time varying covariate names into something that places nicer with the extract function
  tv_cov_names = c()
  for(lll in names(covariate_list)){
    if(class(covariate_list[[lll]]) %in% c('RasterBrick', "RasterStack")){
      tv_cov_names = append(tv_cov_names, lll)
      if (nrow(period_map) > 1)  names(covariate_list[[lll]]) = paste0('XXX',1:length(names(covariate_list[[lll]]))) #use XXX as a place holder
    }
    
  }
  
  #extract the rasters by points
  locations = SpatialPoints(coords = as.matrix(df[,.(longitude,latitude)]), proj4string = CRS(proj4string(covariate_list[[1]])))
  
  cov_values = as.data.frame(lapply(covariate_list, function(x) raster::extract(x, locations)))
  
  #fix the names of the time varying covariates
  if (nrow(period_map) > 1)  names(cov_values) = sub('XXX', '', names(cov_values))
  
  #right now 4 periods are assumed-- make more flexible later
  
  #reconcile time varying covariates
  #start by converting the year variable into periods
  #if period is empty -- infer
  if(is.null(period_var)){
    
    year_map = as.data.frame(sort(unique(df[,year])))
    year_map$period_hold = 1:nrow(year_map)
    names(year_map) = c('year','period_hold')
    
    df = merge(df, year_map, by = 'year', sort =F) #sorting screws up things relative to cov values
    period_var = 'period_hold'
  }
  
  #sort the dataset by rid
  
  setorderv(df, cols = c(id_col))
  
  
  if(return_only_results){
    df = df[, c(id_col, 'longitude', 'latitude', period_var), with = F]
  }
  
  #combine the dataset with the covariate extractions
  df = cbind(df, cov_values)
  
  #reshape to fill out the columns-- these just get the tv cov names in a way that is needed for two different efforts
  #make this less likey to cause a hiccup
  tv_cov_colnames = grep(paste(tv_cov_names, collapse = "|"), names(df), value = T) #unlisted
  tv_cov_colist = lapply(tv_cov_names, function(x) grep(x, names(df), value = T))
  
  #Keep only where period of the data matches the period of the covariate for time varying covariates
  if(reconcile_timevarying & !is.null(tv_cov_names)){
    
    #reshapes the time varying covariates long
    df = melt(df, id.vars = names(df)[!names(df) %in% tv_cov_colnames], measure = tv_cov_colist, value.name = tv_cov_names, variable.factor =F)
    
    #melt returns different values of variable based on if its reshaping 1 or 2+ columns.
    #enforce that it must end with the numeric after the period
    df <- df[,variable:= as.numeric(substring(variable,regexpr("(\\.[0-9]+)$", variable)[1]+1))]
    
    #if(length(tv_cov_colist)==1 & (packageVersion("data.table") < package_version("1.10.0"))){
    #  setnames(df,paste0(tv_cov_names[[1]],'1'), tv_cov_names)
    #}
    
    #keep only where the periods match
    df <- merge(df, period_map, by.x = period_var, by.y = 'data_period')
    df = df[period_id == variable, ]
    
    #clean up
    df = df[,variable := NULL]
    
  }
  
  #clean up
  if(period_var=='period_hold'){
    df = df[,period_hold := NULL]
  }
  
  #centre_scale the results
  if(centre_scale){
    design_matrix = data.frame(df[,names(covariate_list)[!grepl('gaul', names(covariate_list))], with =F])
    cs_df <- getCentreScale(design_matrix)
    design_matrix <- centreScale(design_matrix, df = cs_df)
    
    #replace the df columns with the design matrix
    df[, names(covariate_list)[!grepl('gaul', names(covariate_list))] := NULL]
    df = cbind(df, design_matrix)
    
  }
  
  df <- df[, year := NULL]
  
  #return the data frame with the year-location specific covariates
  if(return_only_results & reconcile_timevarying){
    df = df[, names(df)[!names(df) %in% c('latitude','longitude')], with = F]
  }
  
  #sort df just to be sure
  setorderv(df, cols = c(id_col))
  if(id_col == 'rrrr') df[,rrrr := NULL]
  
  if(centre_scale){
    ## add on rows to covs_cs_df for the gaul_codes
    to.add <- colnames(df)[grepl("gaul_code_", colnames(df))]
    m.add  <- rep(0, length(to.add))
    s.add  <- rep(1, length(to.add))
    to.add <- data.frame(name = to.add, mean = m.add, sd = s.add)
    cs_df <- rbind(cs_df, to.add)
    return(list(covs = df, cs_df = cs_df))
  } else{
    return(df)
  }
}
