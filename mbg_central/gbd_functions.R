
###################
### Setting up ####
###################
#pacman::p_load(data.table, dplyr, parallel, plyr, RMySQL, stringr)

####################################################################################################################################################
# 															   Table of Contents
####################################################################################################################################################

## Base
## get_con :
## run_query :

## Pulls
## get_ages
## get_sexes
## get_location_hierarchy
## get_demographics
## get_populations
## get_covariate_metadata
## get_covariates
## get_citations

## Utilities
## detect_demographics
## age_sex_spec
## age_sex_merge
## agg_vars
## get_names
## make_square


####################################################################################################################################################
# 																	 Base
####################################################################################################################################################


get_con <- function(dbname, host) {
  con <- suppressWarnings(src_mysql(dbname = dbname, host = host, user = "<<<<REDACTED>>>>>", password = "<<<<REDACTED>>>>>"))
  return(con)
}

run_query <- function(dbname, host, query) {
  con <- get_con(dbname, host)
  con %>% tbl(sql(query)) %>% collect(n=Inf) %>% data.table
}


####################################################################################################################################################
# 																	 Pulls
####################################################################################################################################################

get_ages <- function() {
  dbname <- "<<<<REDACTED>>>>>"
  host   <- "<<<<REDACTED>>>>>"
  query  <- "SELECT * FROM shared.age_group"
  run_query(dbname, host, query)
}

#####################################################################################################################################################

get_sexes <- function() {
  dbname <- "<<<<REDACTED>>>>"
  host   <- "<<<<REDACTED>>>>"
  query  <- "SELECT * FROM shared.sex"
  run_query(dbname, host, query)
}

#####################################################################################################################################################

get_location_hierarchy <- function(location_set_version_id, china.fix=FALSE) {
  dbname <- "<<<<REDACTED>>>>"
  host   <- "<<<<REDACTED>>>>"
  ## Find active version id if not specified
  query  <- paste0("SELECT * FROM shared.location_hierarchy_history WHERE location_set_version_id=", location_set_version_id)
  df <- suppressWarnings(run_query(dbname, host, query))
  
  ## China hierarchy fix
  if (china.fix) df <- china_hierarchy_fix(df)
  
  ## Create hierarchy
  hierarchy <- str_split_fixed(df$path_to_top_parent, ",", max(df$level) + 1) %>% data.table
  hierarchy <- hierarchy[, lapply(.SD, as.numeric)]
  setnames(hierarchy, names(hierarchy), paste0("level_", seq(0, max(df$level))))
  df <- cbind(df, hierarchy)
  return(df)
}

#####################################################################################################################################################

get_demographics <- function(location_set_version_id, 
                             year_start, year_end, 
                             by_sex=1, by_age=1, 
                             custom_sex_id=NULL, custom_age_group_id=NULL) {
  ## Locations
  locs <- get_location_hierarchy(location_set_version_id)[level >= 3 & level <6,]$location_id
  ## Years
  years <- seq(year_start, year_end,1)
  ## Sexes
  if (is.blank(custom_sex_id)) { 	
    if (by_sex == 1) {
      sexes <- c(1,2)
    } else if (by_sex == 0) {
      sexes <- 3
    }
  } else {
    sexes <- custom_sex_id
  }
  ## Ages
  if (is.blank(custom_age_group_id)) {
    if (by_age==1) {
      dbname <- "<<<<REDACTED>>>>"
      host   <- "<<<<REDACTED>>>>"
      query  <- "SELECT age_group_id FROM shared.age_group_set_list WHERE age_group_set_id=1"
      ages <- run_query(dbname, host, query)$age_group_id
    } else if (by_age==0) {
      ages <- 22
    }
  } else {
    ages <- custom_age_group_id
  }
  
  ## Expand
  df <- data.table(expand.grid(location_id = locs, year_id=years, sex_id=sexes, age_group_id=ages))
  ## Force integer
  df <- df[, lapply(.SD, as.character)]
  df <- df[, lapply(.SD, as.integer)]
  
  return(df)
}

#####################################################################################################################################################

get_populations <- function(location_set_version_id,
                            year_start, year_end,
                            by_sex=1, by_age=1,
                            custom_sex_id=NULL, custom_age_group_id=NULL) {
  
  ## Make Frame
  df <- get_demographics(location_set_version_id, year_start, year_end, by_sex, by_age, custom_sex_id, custom_age_group_id)
  for (ids in c("location_id", "year_id", "age_group_id", "sex_id")) {
    assign(ids, paste0(unique(df[[ids]]), collapse=","))
  }
  
  ## Pull
  dbname <- "<<<<REDACTED>>>>"
  host   <- "<<<<REDACTED>>>>"
  query <- paste0("SELECT 
                  o.age_group_id,
                  year_id,
                  o.location_id,
                  o.sex_id,
                  pop_scaled 
                  FROM
                  mortality.output o
                  LEFT JOIN
                  mortality.output_version ov using (output_version_id)
                  LEFT JOIN
                  shared.age_group a using (age_group_id)
                  LEFT JOIN
                  shared.location l using (location_id)
                  LEFT JOIN
                  shared.sex s using (sex_id)
                  WHERE
                  ov.is_best = 1
                  and year_id in (", year_id, ")
                  and o.location_id in (", location_id, ")
                  and o.sex_id in (", sex_id, ")
                  and o.age_group_id in (", age_group_id, ")")
  run_query(dbname, host, query)
}

#####################################################################################################################################################

get_covariate_metadata <- function(list=NULL) {
  ## Where clause
  if (is.blank(list)) {
    where <- ""
  } else {
    where <- paste0("WHERE LOWER(covariate_name_short) in (", tolower(toString(shQuote(list))), ")")
  }
  
  ## Pull
  dbname <- "<<<<REDACTED>>>>"
  host <- "<<<<REDACTED>>>>"
  query <- paste0("SELECT * FROM shared.covariate ", where)
  df <- run_query(dbname, host, query)
  
  ## Convert names to lowercase
  df$covariate_name_short <- tolower(df$covariate_name_short)
  
  ## Throw an error check if not int output
  if (!is.blank(list)) {
    return <- df$covariate_name_short
    if (length(list[which(!list %in% return)]) > 0) {
      stop(paste0("The following covariates do not exist in the db: ", toString(list[which(!list %in% return)])))
    }
  }
  return(df)
}

#####################################################################################################################################################

pull_covariate <- function(cov, ci=FALSE) {
  
  dbname <- "<<<<REDACTED>>>>"
  host <- "<<<<REDACTED>>>>"
  
  if (ci) ci_query  <- paste0(", model.lower_value AS '", tolower(paste0(cov, "_lower")), "',
                              model.upper_value AS '", tolower(paste0(cov, "_upper")), "'")
  if (!ci) ci_query <- ""
  
  query <- paste0("SELECT
                  model.location_id,
                  model.year_id,
                  model.age_group_id,
                  model.sex_id,
                  model.mean_value AS '", tolower(cov), "'",
                  ci_query,               
                  "FROM covariate.model
                  JOIN covariate.model_version ON model.model_version_id=model_version.model_version_id
                  JOIN covariate.data_version ON model_version.data_version_id=data_version.data_version_id
                  JOIN shared.covariate ON data_version.covariate_id=covariate.covariate_id
                  JOIN shared.location ON model.location_id=location.location_id
                  JOIN shared.age_group ON model.age_group_id=age_group.age_group_id
                  WHERE covariate.last_updated_action!='DELETE' 
                  AND is_best=1 AND covariate.covariate_name_short = '", cov, "'")
  data <- run_query(dbname, host, query)
  
  if (!(4749 %in% unique(data$location_id))|!(44533 %in% unique(data$location_id))) data <- janky_covariate_fix(data, cov)
  
  return(data)
  
}

#####################################################################################################################################################

get_covariates <- function(list, ci=FALSE) {
  
  ## Get metadata
  meta <- get_covariate_metadata(list)
  
  ## Parse into age/sex specific, age specific, sex specific
  meta <- data.table(meta)
  age_sex <- meta[by_sex == 1 & by_age == 1, .(covariate_name_short)]
  age <- meta[by_sex == 0 & by_age==1, .(covariate_name_short)]
  sex <- meta[by_sex == 1 & by_age==0, .(covariate_name_short)]
  all <- meta[by_sex == 0 & by_age==0, .(covariate_name_short)]
  
  ## Reorder from the most specific to most general
  list <- unique(rbind(age_sex, age, sex, all)$covariate_name_short)
  
  # Pull and Merge Covariates
  flag <- 0
  output <- NULL
  for (cov in list) {
    
    data <- pull_covariate(cov, ci)
    
    if (flag == 0) {
      output <- data
      flag <- 1
    } else {
      output <- age_sex_merge(output, data)
    }
  }
  
  return(output)
}

#####################################################################################################################################################

get_covariate_version <- function(list) {
  dbname <- "<<<<REDACTED>>>>"
  host <- "<<<<REDACTED>>>>"
  
  query <- paste0("SELECT
                  data_version.covariate_id,
                  covariate_name_short,
                  model_version.data_version_id,
                  model_version.model_version_id,
                  model_version.last_updated,
                  model_version.description 			
                  FROM covariate.model_version 
                  JOIN covariate.data_version ON model_version.data_version_id=data_version.data_version_id
                  JOIN shared.covariate ON data_version.covariate_id=covariate.covariate_id
                  WHERE covariate.last_updated_action!='DELETE' 
                  AND is_best=1 AND LOWER(covariate_name_short) in (", tolower(toString(shQuote(list))), ")")
  
  run_query(dbname, host, query)
}


#####################################################################################################################################################

get_citations <- function(nids) {
  dbname <- "<<<<REDACTED>>>>"
  host <- "<<<<REDACTED>>>>"
  query <- paste0("SELECT 
                  nid, field_citation_value, series_title
                  FROM
                  shared.mv_citation
                  WHERE nid in (", nids, ")")
  run_query(dbname, host, query)
}


####################################################################################################################################################
# 																	 Utility
####################################################################################################################################################

is.blank <- function(x) {
  any(is.null(x))  || any(is.na(x))  || any(is.nan(x)) 
}

#####################################################################################################################################################


detect_demographics <- function(df) {
  vars <- c("location_id", "year_id", "age_group_id", "sex_id")
  ## Check if data frame has id's and only return on those that exist
  vars <- names(df)[names(df) %in% vars]
  demos <- lapply(vars, function(x) unique(df[[x]]))
  names(demos) <- vars
  return(demos)
}

#####################################################################################################################################################


age_sex_spec <- function(df) {	
  demos <- detect_demographics(df)
  by_age <- ifelse(22 %in% demos$age_group_id, 0, 1)
  by_sex <- ifelse(3 %in% demos$sex_id, 0 , 1)
  spec <- cbind(by_age, by_sex) %>% data.table
  return(spec)
}

#####################################################################################################################################################

age_sex_merge <- function(df1, df2) {
  
  ## Find specificity
  spec1 <- age_sex_spec(df1)
  spec2 <- age_sex_spec(df2)
  test <- data.table(spec1 == spec2)
  
  ## Merge
  cols <- c("location_id", "year_id", "age_group_id", "sex_id")
  drop_cols <- NULL
  merge_cols <- col
  
  ## If age and sex match
  if (test$by_age & test$by_sex) {
    df <- merge(df1, df2, by=cols)
  } else {
    ## If age matches but sex doesn't match
    if (test$by_age & !test$by_sex) {
      drop_cols <- "sex_id"
    }	
    ## If age doesnt match and sex matches
    else if (!test$by_age & test$by_sex) {
      drop_cols <- "age_group_id"
    }
    ## If neither match
    else {
      drop_cols <- c("sex_id", "age_group_id")
    }
    ## Merge
    merge_cols <- cols[!cols %in% drop_cols]
    df <- merge(df1, df2[, drop_cols := NULL, with=F], by=merge_cols)
  }
  return(df)
}

#####################################################################################################################################################

agg_vars <- function(df, location_set_version_id, vars, parent_ids) {
  
  ## Merge on populations if need be
  if (!("pop_scaled" %in% names(df))) {
    demos <- detect_demographics(df)
    year_start <- min(demos$year_id)
    year_end <- max(demos$year_id)
    custom_sex_id <- demos$sex_id
    custom_age_group_id <- demos$age_group_id
    pops <- get_populations(location_set_version_id, year_start, year_end, custom_sex_id=custom_sex_id, custom_age_group_id=custom_age_group_id)
    df <- merge(df, pops, by=c("location_id", "year_id", "age_group_id", "sex_id")) 
  }	
  
  ## Merge on parent_ids if need be
  if (!("parent_id" %in% names(df))) {
    locs <- get_location_hierarchy(location_set_version_id)[, .(location_id, parent_id)]
    df <- merge(df, locs, by="location_id")
  }
  
  ## Subset to requested parent_ids
  df <- df[parent_id %in% parent_ids]
  
  ## Aggregate estimates to the parent_id [ sum(var * pop) /sum(pop) ]
  df[, (vars) := lapply(.SD, function(x) sum(x * df[['pop_scaled']])/sum(df[['pop_scaled']])), .SDcols=vars, by=c("parent_id", "year_id", "age_group_id", "sex_id")]
  ## De-duplicate so get one set of estimates
  df <- unique(df[, c("parent_id", "year_id", "age_group_id", "sex_id", vars), with=F])
  ## Rename parent_id -> location_id
  setnames(df, "parent_id", "location_id")
  
  return(df)        
}

#####################################################################################################################################################

get_names <- function(df) {
  
  ## Given a data.table with location, age_group, sex, merge on namescu
  required <- c("location_id", "age_group_id", "sex_id") 
  missing <- required[!required %in% names(df)]
  if (length(missing) > 0) stop(paste0("Missing required columns: ", toString(required)))
  
  ## Detect what needs names
  cols <- c("ihme_loc_id", "location_name", "age_group_name", "sex")
  need <- cols[!cols %in% names(df)]
  
  ## Names
  if ("ihme_loc_id" %in% need) df <- merge(df, get_location_hierarchy(41)[,.(location_id, ihme_loc_id)], by="location_id", all.x=T)
  if ("location_name" %in% need) df <- merge(df, get_location_hierarchy(41)[,.(location_id, location_name)], by="location_id", all.x=T)
  if ("age_group_name" %in% need) df <- merge(df, get_ages()[,.(age_group_id, age_group_name)], by="age_group_id", all.x=T)
  if ("sex" %in% need) df <- merge(df, get_sexes()[,.(sex_id, sex)], by="sex_id", all.x=T)
  
  return(df)
}

#####################################################################################################################################################

make_square <- function(location_set_version_id,
                        year_start, year_end,
                        by_sex=1, by_age=1,
                        custom_sex_id=NULL, custom_age_group_id=NULL,
                        covariates=NULL, population=FALSE) {	
  ## Skeleton
  df <- get_demographics(location_set_version_id, 
                         year_start, year_end, 
                         by_sex, by_age, 
                         custom_sex_id, custom_age_group_id)
  
  ## Covariates
  if (!is.null(covariates)) {
    covs <- get_covariates(covariates)
    df <- age_sex_merge(df, covs)
  }
  
  ## Population
  if (population) {
    pops <- get_populations(location_set_version_id, year_start, year_end, by_sex, by_age)
    df <- age_sex_merge(df, pops)
  }
  
  return(df)
  
}

################################################
## Download covariate estimates
################################################

# Sample usage
# x <- get_cov_estimates('LDI_pc')
# y <- get_cov_estimates('LDI_pc', filters=list(year_id=c(1990,2013), location_id=c(6,144)))
# z <- get_cov_estimates('LDI_pc', filters=list(year_id=c(1990,2013), location_id=6))

get_cov_estimates <- function(covariate_name_short, filters=list()) {
  
  # Setup base query
  base_query <- sprintf("
                        SELECT
                        model.model_version_id,
                        covariate.covariate_id,
                        covariate.covariate_name_short,
                        model.location_id,
                        location.location_name,
                        model.year_id,
                        model.age_group_id,
                        age_group.age_group_name,
                        model.sex_id,
                        model.mean_value,
                        model.lower_value,
                        model.upper_value
                        FROM covariate.model
                        JOIN covariate.model_version ON model.model_version_id=model_version.model_version_id
                        JOIN covariate.data_version ON model_version.data_version_id=data_version.data_version_id
                        JOIN shared.covariate ON data_version.covariate_id=covariate.covariate_id
                        JOIN shared.location ON model.location_id=location.location_id
                        JOIN shared.age_group ON model.age_group_id=age_group.age_group_id
                        WHERE covariate_name_short='%s'", covariate_name_short)
  
  query_params <- c()
  # Default to best model if not specified
  if (!("model_version_id" %in% names(filters))) {
    query_params <- c(query_params, "AND is_best=1")
  }
  
  # Implement additional filters
  for (f in names(filters)){
    
    fvalues <- paste(unlist(filters[f]),collapse=",")
    
    # Fix table ambiguity for known join columns
    if (f %in% c("age_group_id","model_version_id","data_version_id","location_id")) {
      f <- paste0("model.",f)
    }
    
    ext_param <- sprintf("AND %s IN (%s)", f, fvalues)
    query_params <- c(query_params, ext_param)
  }
  
  # Construct filtered query
  query <- c(base_query, query_params)
  query <- paste(query, collapse=" ")
  
  # Get and return estimates
  conn <- dbConnect(RMySQL::MySQL(), host="<<<<REDACTED>>>>", username="<<<<REDACTED>>>>", password="<<<<REDACTED>>>>")
  dbSendQuery(conn, "SET NAMES utf8")
  estimates <- dbGetQuery(conn, query)
  dbDisconnect(conn)
  
  Encoding(estimates$location_name) <- "UTF-8"
  
  return(estimates)
}


####################################################################################################################################################
# 															     Janky Shit
####################################################################################################################################################

# ## 4749, 44533 fix to create if missing
janky_covariate_fix <- function(df, var) {
  
  ## Merge on populations
  spec <- age_sex_spec(df)
  year_start <- min(df$year_id)
  year_end <- max(df$year_id)
  pops <- get_populations(41, year_start, year_end, spec$by_sex, spec$by_age)
  df <- merge(df, pops, by=c("location_id", "year_id", "sex_id", "age_group_id"))
  ## Aggregate
  if (!(4749 %in% unique(df$location_id))) {
    gbr_4749 <- agg_vars(df, location_set_version_id=41, vars=var, parent_id=4749)
    df <- rbind(df, gbr_4749, fill=T)
  }
  if (!(44533 %in% unique(df$location_id))) {
    chn_44533 <- agg_vars(df, location_set_version_id=41, vars=var, parent_id=44533)
    df <- rbind(df, chn_44533, fill=T)
  }
  # Drop the pop_scaled column
  df[, pop_scaled := NULL]
  
  return(df)
}

china_hierarchy_fix <- function(df) {
  
  ## Reroute China subnationals to level 4, make CHN w/o HKG and MAC, HKG, MAC level 3
  ## Remove China
  df <- df[location_id != 6,]
  ## Remove from path to top parent
  df <- df[, path_to_top_parent := gsub(",6,", ",", path_to_top_parent)]
  ## Reroute 44533, 354, 361 to level 3
  df <- df[location_id %in% c(44533, 354, 361), level := 3]
  ## Reroute the other subnationals to level 4
  df <- df[grepl("CHN", df$ihme_loc_id) & level == 5, level := 4]
  
  return(df)
  
}
