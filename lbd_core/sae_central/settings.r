####################################################################################################
## Description: Define two functions related to model settings:
##
##              'get_settings':
##                load settings from the 'settings' CSV into the R global environment and evaluate
##                to the appropriate data type.
##
##              'check_settings':
##                run checks on the 'settings' CSV file. That is, check that everything that needs
##                to be specified is there and in the right format. Similarly, check that all files
##                listed in settings exist and are in the right format and that all files and
##                settings are self-consistent. Note that 'check_settings' runs 'get_settings' so
##                all settings are loaded into the global environment as a side-effect.
##
## Inputs:      main_dir [character] -- the directory where 'settings.csv' is located.
##
## NOTES:       you must be in the sae_central directory for check_settings() to work (specifically
##                because check_settings() verifies that the appropriate model code is present).
####################################################################################################

require(data.table)
require(Matrix)
require(sp)

get_settings <- function(main_dir) {

  # load settings
  settings <- read.csv(paste0(main_dir, "/settings.csv"), stringsAsFactors=F, header=F)

  # loop over settings
  for (var in 1:nrow(settings)) {
    arg_name <- settings[var, 1]
    arg_value <- settings[var, 2]

    # assign settings in the global environment
    exp <- try(assign(arg_name, eval(parse(text = arg_value)), envir = .GlobalEnv), silent = T)
    if(class(exp) == "try-error") assign(arg_name, arg_value, envir = .GlobalEnv)
  }

  # create (theoretically) optional variables if not provided
  optional <- c("covars", "covars_as", "covars_trans", "covar_file", "covar_as_file", "geoagg_files")
  for (var in setdiff(optional, ls(envir=.GlobalEnv))) assign(var, NULL, envir=.GlobalEnv)

  return("Settings loaded")
}

check_settings <- function(main_dir) {

  # check that the settings file exists
  if(!file.exists(paste0(main_dir, "/settings.csv"))) stop(paste0(main_dir, "/settings.csv does not exist"))

  # load settings
  get_settings(main_dir)

  # check that all necessary settings are present
  all_settings <- c("model", "area_var", "years", "ages", "sexes",
                    "covars", "covars_as", "covars_trans", "n.sims",
                    "adjmat_file", "events_file", "pop_file", "covar_file", "covar_as_file",
                    "geoagg_files", "shape_file", "age_std_file", "temp_dir")

  miss <- all_settings[!all_settings %in% ls(envir=.GlobalEnv)]
  if(length(miss)) stop(paste("Settings missing:", paste(miss, collapse="; ")))

  # check that all settings are the right R type
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    if(class(x) == "integer") return(TRUE)
    else if(class(x) == "numeric" & max(abs(x - round(x))) < tol) return(TRUE)
    else return(FALSE)
  }

  if(class(area_var) != "character") stop("'area_var' should be type character")
  if(!is.wholenumber(years)) stop("'years' should be a whole number of type integer or numeric")
  if(!is.wholenumber(ages)) stop("'ages' should be a whole number of type integer or numeric")
  if(!is.wholenumber(sexes)) stop("'sexes' should be a whole number of type integer or numeric")
  if(!is.null(covars)) if(class(covars) != "character") stop("'covars' should be type character or NULL")
  if(!is.null(covars_as)) if(class(covars_as) != "character") stop("'covars_as' should be type character or NULL")
  if(!is.null(covars_trans)) if(class(covars_trans) != "character") stop("'covars_trans' should be type character or NULL")
  if(!is.wholenumber(n.sims)) stop("'n.sims' should be a whole number of type integer or numeric")
  if(class(adjmat_file) != "character") stop("'adjmat_file' should be type character")
  if(class(events_file) != "character") stop("'events_file' should be type character")
  if(class(pop_file) != "character") stop("'pop_file' should be type character")
  if(!is.null(covar_file)) if(class(covar_file) != "character") stop("'covar_file' should be type character or NULL")
  if(!is.null(covar_as_file)) if(class(covar_as_file) != "character") stop("'covar_as_file' should be type character or NULL")
  if(!is.null(geoagg_files)) if(class(geoagg_files) != "character") stop("'geoagg_files' should be type character or NULL")
  if(class(age_std_file) != "character") stop("'age_std_file' should be type character")
  if(class(temp_dir) != "character") stop("'temp_dir' should be type character")

  # check that all settings are the right length
  if(length(model) != 1) stop("'model' should be length 1")
  if(length(area_var) != 1) stop("'area_var' should be length 1")
  if(length(years) <= 1) stop("'years' should be length >2")
  if(length(ages) < 1) stop("'ages' should be length >=1")
  if(length(sexes) < 1 | length(sexes) > 2) stop("'sexes' should be length 1 or 2")
  if(!is.null(covars)) if(length(covars) < 1) stop("'covars' should be length >=1 or NULL")
  if(!is.null(covars_as)) if(length(covars_as) < 1) stop("'covars_as' should be length >=1 or NULL")
  if(!is.null(covars_trans)) if(length(covars_trans) < 1) stop("'covars_trans' should be length >=1 or NULL")
  if(length(n.sims) != 1) stop("'n.sims' should be length 1")
  if(length(adjmat_file) != 1) stop("'adjmat_file' should be length 1")
  if(length(events_file) != 1) stop("'events_file' should be length 1")
  if(length(pop_file) != 1) stop("'pop_file' should be length 1")
  if(!is.null(covar_file)) if(length(covar_file) != 1) stop("'covar_file' should be length 1 or NULL")
  if(!is.null(covar_as_file)) if(length(covar_as_file) != 1) stop("'covar_as_file' should be length 1 or NULL")
  if(!is.null(geoagg_files)) if(length(geoagg_files) < 1) stop("'geoagg_files' should be length >= 1 or NULL")
  if(length(shape_file) != 1) stop("'shape_file' should be length 1")
  if(length(age_std_file) != 1) stop("'age_std_file' should be length 1")
  if(length(temp_dir) != 1) stop("'temp_dir' should be length 1")

  # check that n.sims > 0 but not crazy big (which will cause memory problems)
  if(n.sims <= 0) stop("'n.sims' must be >= 1")
  if(n.sims > 1e4) stop("setting 'n.sims' too high which may cause memory problems")

  # check that all necessary model code exists for the selected model
  all_mod_files <- c(paste0("mod_", model, ".cpp"), paste0("fit_mod_", model, ".r"))
  miss <- setdiff(all_mod_files, dir("models"))
  if(length(miss)) stop(paste("Model code missing:", paste(miss, collapse="; ")))

  # check that the model selected is appropriate
  if (length(ages) > 1 & model %in% c("1b", "2b")) stop(paste("Model", model, "is only for use with a single age group (i.e., 'ages' is length 1)"))
  if (length(ages) == 1 & !model %in% c("1b", "2b")) stop(paste("Model", model, "is only for use with multiple age groups (i.e., 'ages' is length 2+)"))

  # check that all input files exist
  all_files <- list("adjmat_file" = adjmat_file, "events_file" = events_file, "pop_file" = pop_file,
                    "shape_file" = shape_file, "age_std_file" = age_std_file)
  if(!is.null(geoagg_files)) all_files[["geoagg_files"]] <- geoagg_files
  if(!is.null(covar_file)) all_files[["covar_file"]] <- covar_file
  if(!is.null(covar_as_file)) all_files[["covar_as_file"]] <- covar_as_file
  miss <- sapply(all_files, function(x) sum(!file.exists(x)))
  miss <- miss[miss > 0]
  if(length(miss)) stop(paste("Input files missing:", paste(names(miss), collapse="; ")))

  # check that the adjacency matrix is correctly formatted
  if(!grepl(".rdata$", tolower(adjmat_file))) stop("adjmat_file must be an .rdata")
  var <- load(adjmat_file)
  if(var != "adjmat") stop("object in 'adjmat_file' should be named 'adjmat'")
  if(!class(adjmat)[1] == "dgTMatrix") stop("object in 'adjmat_file' should be a dgTMatrix")
  if(!isSymmetric(adjmat)) stop("adjacency matrix in 'adjmat_file' should be symmetric")
  areas <- 1:nrow(adjmat) - 1
  rm(adjmat)

  # check that the events file is formatted correctly and consistent with the settings for area_var and years and with the adjacency matrix
  if(!grepl(".rdata$", tolower(events_file))) stop("events_file must be an .rdata")
  var <- load(events_file)
  if(var != "events") stop("object in 'events_file' should be named 'events'")
  if(!class(events)[1] == "data.table") stop("object in 'events_file' should be a data.table")
  all_vars <- c(area_var, "year", "sex", "age", "events")
  miss <- setdiff(all_vars, names(events))
  if(length(miss)) stop(paste("Variables missing from 'events_file':", paste(miss, collapse="; ")))
  if(sum(is.na(events[, all_vars, with=F]))) stop("NAs in 'events_file'")
  setkeyv(events, c(area_var, "year", "sex", "age"))
  if(anyDuplicated(events, by=key(events))) stop("there are duplicated area-year-sex-ages in 'events_file'")
  if(sum(!events[[area_var]] %in% areas)) stop("not all areas in 'events_file' are in 'adjmat_file'")
  if(sum(!years %in% events$year)) stop("not all years in 'years' are present in 'events_file'")
  if(sum(!sexes %in% events$sex)) stop("not all sexes in 'sexes' are present in 'events_file'")
  if(sum(!ages %in% events$age)) stop("not all ages in 'ages' are present in 'events_file'")
  if(events[, min(events)] < 0) stop("there are negative event counts in 'events_file'")
  rm(events)

  # check that the population file is formatted correctly and consistent with the settings for area_var, years, sexes, and ages and with the adjacency matrix
  if(!grepl(".rdata$", tolower(pop_file))) stop("pop_file must be an .rdata")
  var <- load(pop_file)
  if(var != "pop") stop("object in 'pop_file' should be named 'pop'")
  if(!class(pop)[1] == "data.table") stop("object in 'pop_file' should be a data.table")
  all_vars <- c(area_var, "year", "sex", "age", "pop")
  miss <- setdiff(all_vars, names(pop))
  if(length(miss)) stop(paste("Variables missing from 'pop_file':", paste(miss, collapse="; ")))
  if(sum(is.na(pop[, all_vars, with=F]))) stop("there are NAs in 'pop_file'")
  setkeyv(pop, c(area_var, "year", "sex", "age"))
  if(anyDuplicated(pop, by=key(pop))) stop("there are duplicated area-year-sex-ages in 'pop_file'")
  if(pop[CJ(areas, years, sexes, ages), sum(is.na(pop))]) stop("there are missing area-year-sex-age combinations in 'pop_file'")
  if(!setequal(areas, pop[[area_var]])) stop("areas in 'pop_file' do not match areas in 'adjmat_file'")
  if(pop[, min(pop)] < 0) stop("there are negative populations in 'pop_file'")

  # check that covariates file is formatted correctly and consistent with the settings for area_var, years, covars, and with the adjacency matrix
  if(!is.null(covars) & is.null(covar_file)) stop("covariates specified but there is not a covariates file specified")
  if(!is.null(covars)) {
    if(!grepl(".rdata$", tolower(covar_file))) stop("covar_file must be an .rdata")
    var <- load(covar_file)
    if(var != "covar") stop("object in 'covar_file' should be named 'covar'")
    if(!class(covar)[1] == "data.table") stop("object in 'covar_file' should be a data.table")
    all_vars <- c(area_var, "year", covars)
    miss <- setdiff(all_vars, names(covar))
    if(length(miss)) stop(paste("Variables missing from 'covar_file':", paste(miss, collapse="; ")))
    if(sum(is.na(covar[, all_vars, with=F]))) stop("there are NAs in 'covar_file'")
    setkeyv(covar, c(area_var, "year"))
    if(anyDuplicated(covar, by=key(covar))) stop("there are duplicated area-years in 'covar_file'")
    if(covar[CJ(areas, years), sum(is.na(get(covars[1])))]) stop("there are missing area-year combinations in 'covar_file'")
    if(!setequal(areas, covar[[area_var]])) stop("areas in 'covar_file' do not match areas in 'adjmat_file'")
    rm(covar)
  }

  # check that age-sex covariates file is formatted correctly and consistent with the settings for area_var, years, covars_as, and with the adjacency matrix
  if(!is.null(covars_as) & is.null(covar_as_file)) stop("age-sex covariates specified but there is not a age-sex covariates file specified")
  if(!is.null(covars_as)) {
    if(!grepl(".rdata$", tolower(covar_as_file))) stop("covar_as_file must be an .rdata")
    var <- load(covar_as_file)
    if(var != "covar") stop("object in 'covar_as_file' should be named 'covar'")
    if(!class(covar)[1] == "data.table") stop("object in 'covar_as_file' should be a data.table")
    all_vars <- c(area_var, "year", "sex", "age", covars_as)
    miss <- setdiff(all_vars, names(covar))
    if(length(miss)) stop(paste("Variables missing from 'covar_as_file':", paste(miss, collapse="; ")))
    if(sum(is.na(covar[, all_vars, with=F]))) stop("there are NAs in 'covar_as_file'")
    setkeyv(covar, c(area_var, "year", "sex", "age"))
    if(anyDuplicated(covar, by=key(covar))) stop("there are duplicated area-year-sex-ages in 'covar_as_file'")
    if(covar[CJ(areas, years, sexes, ages), sum(is.na(get(covars[1])))]) stop("there are missing area-year-sex-age combinations in 'covar_as_file'")
    if(!setequal(areas, covar[[area_var]])) stop("areas in 'covar_file' do not match areas in 'adjmat_file'")
    rm(covar)
  }

  # check that if both covars and covars_as are specified such that there is no overlap
  if(length(intersect(covars, covars_as))) stop("'covars' and 'covars_as' cannot contain any of the same items")

  # check that covars_trans is appropriately specified and consistent with covars, covars_as, covar_file, and covar_as_file
  if(!is.null(covars_trans)) {
    if(is.null(names(covars_trans))) stop("'covars_trans' must be a named character vector")
    if(sum(!names(covars_trans) %in% c(covars, covars_as))) stop("the name of each item in 'covars_trans' must also be in 'covars' or 'covars_as'")

    covar <- lapply(c(covar_file, covar_as_file), function(x) { load(x); return(covar) })
    if(length(covar) > 1) covar <- merge(covar[[1]], covar[[2]], by=c(area_var, "year"))
    else covar <- covar[[1]]

    for (var in names(covars_trans)) {
      fun <- try(eval(parse(text=paste("function (x)", covars_trans[var]))), silent=T)
      if(class(fun) == "try-error") stop("each item in 'covars_trans' should be a string specifying a function of x")
      temp <- try(fun(covar[[var]]))
      if(class(temp) == "try-error") stop("each item in 'covars_trans' should be a string specifying a function of x")
      if(sum(is.na(temp))) stop("'covars_trans' settings introduce NAs in 'covar_file' or 'covar_as_file'")
    }
    rm(covar)
  }

  # check that the geoagg files are formatted correctly, consistent with the settings for area_var, sexes, and ages, and consistent with pop_file
  if(!is.null(geoagg_files)) {
    if(is.null(names(geoagg_files))) stop("geoagg_files must be a named vector")
    for (level in names(geoagg_files)) {
      if(!grepl(".rdata$", tolower(geoagg_files[level]))) stop("each element in geoagg_files must be an .rdata")
      var <- load(geoagg_files[level])
      if(var != "weights") stop("object in 'geoagg_files' should be named 'weights'")
      if(!class(weights)[1] == "data.table") stop("object in 'geoagg_files' should be a data.table")
      all_vars <- c(level, area_var, "year", "sex", "age", "pop")
      miss <- setdiff(all_vars, names(weights))
      if(length(miss)) stop(paste("'geoagg_files' is missing variables:", paste(miss, collapse=", ")))
      if(sum(is.na(weights[, all_vars, with=F]))) stop("NAs in 'geoagg_files'")
      setkeyv(weights, c(level, area_var, "year", "sex", "age"))
      if(anyDuplicated(weights, by=key(weights))) stop("there are duplicated level-area-year-sex-ages in 'geoagg_files'")
      setkeyv(weights, c(level, "year", "sex", "age"))
      if(weights[CJ(unique(get(level)), years, sexes, ages), sum(is.na(pop))]) stop("there are missing level-year-sex-age combinations in 'geoagg_files'")
      if(sum(!ages %in% weights$age)) stop("not all ages in 'ages' are present in 'geoagg_files'")
      if(weights[, min(pop)] < 0) stop("there are negative populations in 'geoagg_files'")
      weights <- merge(weights[, list(geo_pop = sum(pop)), keyby=c(area_var, "year", "sex", "age")],
                       pop[, list(pop = sum(pop)), keyby=c(area_var, "year", "sex", "age")], all.x=T)
      if(weights[, max(abs(geo_pop - pop)) > 1e-10]) stop("populations at the 'area_var' level should be the same in 'geoagg_files' as in 'pop_file'")
      rm(weights)
    }
  }

  # check that the shapefiles are formatted correctly and consistent with the adjacency matrix
  if(!grepl(".rdata$", tolower(shape_file))) stop("shape_file must be an .rdata")
  map <- get(load(shape_file))
  if(class(map)[1] != "SpatialPolygonsDataFrame") stop("Object in 'shape_file' should be a SpatialPolygonsDataFrame")
  if(length(map) != length(areas)) stop("object in 'shape_file' does not have the same number of areas as 'adjmat' in 'adjmat_file'")
  if(sum(sapply(map@polygons, function(x) x@ID) != areas)) stop("polygon IDs in 'shape_file' do not have the correct area variable coding (i.e., sequentially from 0) or are not sorted properly")

  # check that the age standard file is correctly formatted and consistent with the 'ages' setting
  if(!grepl(".csv$", tolower(age_std_file))) stop("age_std_file must be an .csv")
  age_std <- fread(age_std_file)
  all_vars <- c("age", "wt")
  miss <- setdiff(all_vars, names(age_std))
  if(length(miss)) stop(paste("Variables missing from 'age_std_file':", paste(miss, collapse="; ")))
  if(sum(is.na(age_std[, all_vars, with=F]))) stop("NAs in 'age_std_file'")
  if(sum(!ages %in% age_std$age)) stop("There are ages in 'ages' missing from 'age_std_file'")
  if(nrow(age_std) != uniqueN(age_std$age)) stop("there are duplicated ages in 'age_std_file'")
  if(sum(age_std$wt) - 1 > 10e-10) stop("weights in 'age_std_file' do not sum to 1")
  rm(age_std)

  return("checks on settings passed")
}
