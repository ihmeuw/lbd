# ##############################################################################
# PLOT DATA GAPS BY TEAM
# Created: April 2, 2018
# Purpose: Determine the total number of sources (NIDs) by country and team
#  for the purposes of prioritizing data seeking
# ##############################################################################

library(data.table)

# Helper function to open a csv path and only keep the unique NIDs (or svy_ids)
get_file_sources <- function(filepath, use_datayears = FALSE){
  df <- fread(filepath)
  names(df) <- tolower(names(df))
  # Use a numeric NID field if it exists
  if ("nid" %in% names(df)){
    if (("numeric" %in% class(df[,nid]))|("integer" %in% class(df[,nid]))){
      df[,svy_id:=nid]
    }
  }
  # If use_datayears is true, then svy_id is unreliable
  # Use years as the "unique" identifier instead
  if (use_datayears){
    if ("year_var" %in% names(df)){
      df[,svy_id:=year_var]
    } else {
      return(data.table())
    }
  }
  # If no svy_id field exists at this point, return an empty data.table
  if (!("svy_id" %in% names(df))){
    return(data.table())
  }
  # Coerce svy_id field to numeric
  df[,svy_id:=as.numeric(svy_id)]
  # Coerce IHME codes to country codes (ISO3)
  df[,country:=substr(as.character(country),1,3)]
  # Keep only unique countries and survey IDs
  df <- unique(df[,c('country','svy_id'),with=F])
  df <- df[!is.na(svy_id),]
  return(df)
}

# Helper function to plot data gaps for each individual indicator
get_indicator_sources <- function(
    ind_name,
    regex_matches,
    dcp_folder="<<<< FILEPATH REDACTED >>>>",
    include_vr=TRUE
  ){
  # Get the list of all indicator sources
  all_matches <- unlist(
    lapply(
      regex_matches,
      function(r) file.path(dcp_folder,grep(r,list.files(dcp_folder),value=T,ignore.case=TRUE))
    )
  )
  if(include_vr==FALSE){
    all_matches <- grep('vr', all_matches, ignore.case=TRUE, value=TRUE, invert=TRUE)
  }
  if (length(all_matches)==0){
    return(data.table())
  }
  # Get data.tables with full list of countries and sources for each filepath
  if (tolower(ind_name) %in% c('focal3')){
    sources_list <- lapply(all_matches,
                           function(x) get_file_sources(x, use_datayears=T))
  } else {
    sources_list <- lapply(all_matches,get_file_sources)
  }
  sources <- rbindlist(sources_list, fill=T)
  sources <- unique(sources)
  sources[,num_sources := 1]
  sources <- sources[,lapply(.SD,sum), .SDcols = c('num_sources'),
                     by=c('country')]
  return(sources)
}


# THIS IS THE MAIN FUNCTION
plot_data_gaps <- function(include_vr=TRUE){
  # You should only be using the cluster for MBG central functions!
  if ("windows" %in% tolower(Sys.info()[1])){
    stop("Please don't run this using Windows.")
  }
  # Load packages
  library(data.table)
  # Source GBD shared functions
  source("<<<< FILEPATH REDACTED >>>>")

  # Define input and output filepaths:
  # - Template file with all ISO3s by stage
  template_file <- "<<<< FILEPATH REDACTED >>>>"
  # - Output file
  out_file <- "<<<< FILEPATH REDACTED >>>>"

  # Read in the template file and format
  templ <- fread(template_file)
  templ <- templ[,c('iso3','Stage','location_name','loc_id'),with=F]
  setnames(templ, 'iso3', 'country')
  templ[, Stage := as.numeric(substr(as.character(Stage),1,1))]
  templ <- templ[Stage < 3,]
  # Get populations and merge on
  templ[,loc_id:=as.integer(loc_id)]
  locs <- templ[,loc_id]
  locs <- locs[!is.na(locs)]
  pops <- get_population(age_group_id = 22,
                         location_id  = locs,
                         sex_id       = 3,
                         year_id      = 2016)
  pops <- pops[,c('location_id','population'),with=F]
  pops[,SmallPop := 1 * (population <= 2.5*1E5)]
  pops[,population:=as.integer(round(population/1E6))]
  names(pops) <- c('loc_id','Pop(Millions)','SmallPop')
  templ <- merge(x=templ,
                 y=pops,
                 by=c('loc_id'),
                 all.x=TRUE)
  templ[country=="GUF",c('Pop(Millions)','SmallPop'):=list(0,1)]
  templ[country=="ESH",c('Pop(Millions)','SmallPop'):=list(0,0)]

  # Define list of projects and relevant regex statements that can be used
  #  to find all data coverage plot files from that project
  lbd_project <- list(
      "Edu"      = list("^edu"),
      "CGF"      = list("^stunting_","^wasting_","^underweight_","^overweight",
                     "^obese","^undernutrition","^low_bw"),
      "Vaccines" = list("^vax_","^dpt","^mcv","^hib3","^pcv","^polio"),
      "U5M"      = list("^died_","^dead_","cbh","sbh"),
      # "CGF:EBF"  = list("^bf","bf_"),
      "CGF:EBF"  = list("^ebf"),
      "Rota/PCV" = list("^rota","^pcv"),
      "Diarrhea" = list("had_diarrhea"),
      "LRI"      = list("^lri","^has_lri"),
      "Focal3"   = list("^lf_","^Oncho","^oncho","^Lymphatic Filariasis"),
      "HIV"      = list("^hiv_test"),
      "HIV_ANC"  = list("^anc_hiv_"),
      "HIVcovs"  = list("^male_circumcision", "^sti_symptoms", "^had_sti_symptoms",
                     "^in_union", "^partner_away", "^condom_last_time",
                     "^multiple_partners", "^num_partners_lifetime",
                     "^debut_age", "^client_sex_worker"),
      "WaSH"     = list("^water_","^sani_data","_water_","^imp"),
      "TB"       = list("^tb_")
  )

  # Get a dataframe containing the number of unique sources (NIDs) for each
  #  project by country
  get_indicator_wrapper <- function(n){
   message(n)
   data <- lbd_project[[n]]
   num_sources <- get_indicator_sources(n,data,include_vr=include_vr)
   num_sources[,project:=n]
   num_sources <- num_sources[,c('project','country','num_sources'),with=F]
   return(num_sources)
  }
  sources_list <- lapply(names(lbd_project),get_indicator_wrapper)
  all_sources <- rbindlist(sources_list)

  # Cast by project name
  project_names <- unique(all_sources[,project])
  project_names_mean_stage1 <- project_names[!(project_names %in% c("Focal3","HIV_ANC"))]
  project_names_mean_stage2 <- project_names[!(project_names %in% c("Focal3","CGF:EBF","HIV_ANC",
                                                                    "HIVcovs","HIV","TB"))]
  casted <- dcast(all_sources, country ~ project, value.var='num_sources',
                  drop=FALSE, fill=0, fun=sum)

  # Merge with template data.table
  full_data <- merge(x=templ,
                     y=casted,
                     by=c('country'),
                     all.x=TRUE)
  setnames(full_data,'country','ISO')
  setnames(full_data,'location_name','Country Name')
  for (source in names(lbd_project)){
    full_data[is.na((source)),(source):=0]
  }
  # Split into Stage 1 and Stage2 data
  stage1 <- full_data[Stage == 1, ]
  stage2 <- full_data[Stage == 2, ]
  # Get measures of data sparsity for each stage and bind together
  stage1$Mean.All <- round(rowMeans(stage1[,project_names_mean_stage1,with=F]),2)
  stage2$Mean.All <- round(rowMeans(stage2[,project_names_mean_stage2,with=F]),2)
  full_data <- rbindlist(list(stage1, stage2))
  full_data <- full_data[order(Stage,Mean.All)]
  full_data <- full_data[,c('Stage','Country Name','ISO','Mean.All','Pop(Millions)',
                            'SmallPop',project_names),with=F]

  # Save output data
  write.csv(full_data, file=out_file, row.names=F)
}

