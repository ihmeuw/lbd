collapse_single_year <- function(cores, in_dir, out_dir='') {

  #load data
  load_edu_dta <- function(file) {
    data <- data.table(read_dta(file))
    print(file)
    if(!"sex" %in% colnames(data)) data[, sex := sex_id]
    if(!"year" %in% colnames(data)) data[, year := year_n]
    if(!"survey_name" %in% colnames(data)) data[,survey_name := type]
    if('age_start' %in% colnames(data)) data[, age := age_start]
    if('age_year' %in% colnames(data)) data[, age := age_year]
    data[,count:=1*pweight]
    data[,num_persons:=1]

    
    # make sure that any unreasonable values of edu_min are set to missing 
    # since they likely represent missing codes that weren't properly handled
    data[edu_min > 50 | edu_min < 0, edu_min := NA]
    
    data <- data[,.(count=sum(count),
                    num_persons=sum(num_persons), 
                    total_weight_squared=sum(count, na.rm = T)^2, 
                    sum_weight_squared=sum(count^2, na.rm = T), 
                    iso3 = unique(iso3)),
                 by=.(age,sex,year,edu_min,edu_max,survey_name,nid,geospatial_id)]
    data[, kish_N := total_weight_squared/sum_weight_squared]
    return(data)
  }
  #load data
  files <- list.files(in_dir,full.names = T,pattern='.dta')
  data.list <- mclapply(files, load_edu_dta ,mc.cores=ifelse(Sys.info()[1]=="Windows", 1, cores))
  edu.data <- rbindlist(data.list,fill=T)

  data <- copy(edu.data)
  #subset to only keep geo_ids (clusters) that only have single year education data
  data <- data[age>14 & geospatial_id != "" & sex %in% c(1,2) & !is.na(age) & !is.na(edu_min) & edu_min >= 0 & age < 140]
  data[,bins := ifelse(edu_min == edu_max, 0, 1)]
  bins <- data[,.(bins = max(bins)), by = .(year,geospatial_id,survey_name,nid)]
  data[,bins := NULL]
  data <- data.table(merge(data, bins, by = c("year","geospatial_id", "survey_name","nid")))
  data <- data[bins == 0]
  
  # make sure edu_min is capped at 18 and is an integer
  data[edu_min>18,edu_min:=18]
  data[,edu_min:=floor(edu_min)]

  return(data)
  
}

