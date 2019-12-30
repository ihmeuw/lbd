################ parallel_collapse_binned.R ###########################################
### Revised 5/2018 by Kate Wilson (kwilson7@uw.edu) 
### Crosswalks binned education (e.g., 10-12 years of education) into
### single years. Crosswalks each cluster in parallel (within each 
### cluster, each bin is crosswalked in serial). 
########### Usually launched by submit_collapse_binned for each binned survey ##########
########################################################################################

######################### PREP ############################################
package_list <- c("data.table", "gtools","parallel","SDMTools", "readstata13")

# Use the package imports function in mbg_central that figures out what setup you're on for the library location
source(paste0("<<<< FILEPATH REDACTED >>>>", '/lbd_core/mbg_central/setup.R'))
load_R_packages(package_list)

cores <- 3

## Get parameters for survey from command line
print(commandArgs())
c.year <- as.numeric(commandArgs()[4])
c.file <- commandArgs()[5]
## Set parameters for training set
space_weight <- .6
time_weight <- (1 - space_weight)
training_size <- 12

######################## LOAD AND PREP DATA ##################################

#load data
if (regexpr("*tabulations*", c.file)==-1) {
  df <- data.table(read.dta13(c.file))
  df[,count := 1*pweight]
  df[, num_persons := 1]
} else {
  df <- data.table(read.dta13(c.file))
  df[, num_persons := count]
  df[, pweight := count]
  df[, level_3 := iso3]
}
#rename columns to standardize
if('age_start' %in% colnames(df)){
  df <- df[, age := age_start]
}
if('age_year' %in% colnames(df)){
  df <- df[, age := age_year]
}
if(!'sex' %in% colnames(df)){
  df <- df[, sex := sex_id]
}
if(!'year' %in% colnames(df)){
  df <- df[, year := year_n]
}
if(!'survey_name' %in% colnames(df)){
  df <- df[, survey_name := type]
}

# make sure that any unreasonable values of edu_min are set to missing 
# since they likely represent missing codes that weren't properly handled
df[edu_min > 50 | edu_min < 0, edu_min := NA]

# drop data missing ages
df <- df[!is.na(age),]

# calculate nearest 5 year age bin for training
df[,age_group := 5*floor(age/5)]
df[age_group>95,age_group:=95]

if (!'age_end' %in% colnames(df)) {
  df[,age_group_end := age_group]
} else {
  df[,age_group_end := 5*floor(age_end/5)]
  df[age_group_end>95,age_group_end:=95]
}


######################## LOAD AND IDENTIFY TRAINING SET FOR SURVEY ##################################
# do this before collapsing data because have to use raw age_group and age_group_end in data to collapse training

####### Read in training data
train <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
setnames(train, "age_start", "age_group")
####### Read in location hierarchy
locs <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

####### Identify training set for survey
# The way we're extracting and saving surveys ensures that each file has at most 1 ISO3 in it
cc.iso <- unique(df$level_3)[1]

# if our starting and ending age groups in data aren't the same, collapse training set to match these groups
if (sum(df$age_group != df$age_group_end, na.rm=T) > 0) {
  unique_age_groups <- unique(df[,c("age_group", "age_group_end")])
  train$age_group_corrected <- apply(train[,c("age_group")], 1, function(d) unique_age_groups[age_group <= d & age_group_end >= d, age_group])
  train <- train[, .(count=sum(count), sample_size=sum(sample_size), proportion=sum(proportion)), by=.(age_group_corrected, sex, year, nid, ihme_loc_id, type, edu_yrs, total, total_ss, mean, sd, mean_se)]
  setnames(train, "age_group_corrected", "age_group")
}

# get region corresponding to country (cc.iso)
target.r <- unique(locs[ihme_loc_id==cc.iso,region_name])

# get super region corresponding to country (cc.iso)
target.sr <- unique(locs[ihme_loc_id==cc.iso,super_region_name])

# spatial distance is 1 for locs outside the super region our survey is in
train[,spatial_distance:=1]

# spatial distance is .66 for locs in same super region as our survey
train[ihme_loc_id %in% unique(locs[super_region_name==target.sr,ihme_loc_id]),spatial_distance:=.66]

# spatial distance is .33 for locs in same region as our survey
train[ihme_loc_id %in% unique(locs[region_name==target.r,ihme_loc_id]),spatial_distance:=.33]

# spatial distance is 0 for locs in same country as our survey
train[ihme_loc_id == cc.iso,spatial_distance:=0]

# calculate temporal distance relative to survey year
train[, temporal_distance := abs(year - c.year) / (max(year) - min(year))]

# calculate the total distance by weighting the calculated spatial and temporal distances by their corresponding passed weights
train[, distance := (spatial_distance * space_weight) + (temporal_distance * time_weight)]

# collapse training by location, survey series (type), and year
rank <- train[,.(distance=mean(distance)),by=.(ihme_loc_id,type,year)]

# create a ranking over collapsed training set, break ties at random
rank[,rank:=rank(distance,ties.method = "random")]

# apply ranking to training set
train <- data.table(merge(train,rank,by=c("ihme_loc_id","type","year","distance")))

# only keep items in training set with rank less than passed training size 
train <- train[rank<=training_size]



##################### COLLAPSE DATA to desired geography, age groups, sex, bins #######


df <- df[,.(count=sum(count),
            num_persons=sum(num_persons),
            total_weight_squared=sum(count, na.rm = T)^2, 
            sum_weight_squared=sum(count^2, na.rm = T), 
            level_3 = unique(level_3)),
         by=.(age,age_group,sex,year,edu_min,edu_max,survey_name,nid,geospatial_id)]
df[, kish_N := total_weight_squared/sum_weight_squared]

df <- df[!is.na(edu_min),]

## drop any clusters that aren't binned (corresponds to edu_collapse_functions where only geo_ids with no bins are kept) ##
df[,bins := ifelse(edu_min == edu_max, 0, 1)]
bins <- df[,.(bins = max(bins)), by = .(year,geospatial_id,survey_name,nid)]
df[,bins := NULL]
df <- data.table(merge(df, bins, by = c("year","geospatial_id", "survey_name","nid")))
df <- df[bins == 1]
df[,bins := NULL]


######################## CROSSWALKING FUNCTIONS ##################################

###### Function to crosswalk data of a single bin (s.data) using the passed training set train
crosswalk_bin <- function(s.data, train) {
  # get the matching training set for the bin (all education years within the bin, inclusive)
  s.train <- train[edu_yrs >= min(s.data[,edu_min]) & edu_yrs <= max(s.data[,edu_max])]
  s.train <- s.train[,.SD,.SDcols=c("rank","edu_yrs","count","age_group","sex")]
  s.train.tot <- s.train[,.(total=sum(count)),by=.(rank,age_group,sex)]  
  s.train <- data.table(merge(s.train,s.train.tot,by=c("rank","age_group","sex")))  
  s.train[,prop:=count/total]
  s.train <- s.train[,.(prop=mean(prop)),by=.(edu_yrs,age_group,sex)]
  s.data <- data.table(merge(s.data,s.train,by=c("age_group","sex"),allow.cartesian=T))
  s.rake <- s.data[,.(rake=sum(prop)),by=.(age,age_group,sex)]
  s.data <- data.table(merge(s.data,s.rake,by=c("age","age_group","sex")))
  s.data[,prop := prop/rake]
  s.data[,count := count*prop]
  s.data[,total_weight_squared := total_weight_squared*prop]
  s.data[,sum_weight_squared := sum_weight_squared*prop]
  s.data[,kish_N := kish_N*prop]
  s.data[,sample_size := sample_size*prop]
  return(s.data)
}

###### Function to crosswalk each cluster (geospatial_id) 
## Params: df = data for a single cluster (geospatial_id); train = training set for the survey
# train = train
train_cluster <- function(df, train) {
  c.cluster <- unique(df$geospatial_id)[1]
  print(c.cluster)
  
  data <- df[age > 14 & geospatial_id != "" & sex %in% c(1,2) & !is.na(age) & !is.na(edu_min) & edu_min >= 0 & age < 140]
  data[edu_min > 18,edu_min:=18]
  data[,edu_min:=floor(edu_min)]
  data[, sample_size := num_persons]
  
  if(nrow(data) != 0){
    ###### Crosswalk each bin (in serial)
    data[,bin:=paste0(edu_min,"_",edu_max)]
    bin_data <- split(data, f=data$bin)
    split.data <- lapply(bin_data, crosswalk_bin, train)
    data <- rbindlist(split.data)
    
    #get total sample size and totals for calculations
    total <- data[,.(total=sum(count),
                     total_ss=sum(sample_size)),
                  by=.(age,age_group,sex,year,nid,geospatial_id,survey_name)]
    data <- data.table(merge(data,total,by=c("age","age_group","sex","year","nid","geospatial_id","survey_name")))
    data[,proportion:=count/total]
    #standard error of proportions
    data[,orig_proportion_se:=sqrt(proportion * (1 - proportion)/total_ss)]
    data[,prop_se:=sqrt(prop * (1 - prop)/training_size)]
    data[,proportion_se:= sqrt((prop_se^2) + (orig_proportion_se)^2)]
    
    #mean and standard error of mean
    means <- data[,.(mean=weighted.mean(x=edu_yrs,w=proportion),
                     sd=wt.sd(x=edu_yrs,w=proportion)),
                  by=.(age,age_group,sex,year,nid,survey_name,geospatial_id)]
    data <- data.table(merge(data,means,by=c("age","age_group","sex","year","nid","geospatial_id","survey_name")))
    data[,mean_se:=sd / sqrt(total_ss)]
    return(data)
  }
}



######################## CROSSWALK EACH CLUSTER IN SURVEY ##################################

###### Split our survey into a list of subsets for each geospatial id
# list where each index is the data at a unique geospatial_id 
df_locs <- split(df, f=df$geospatial_id)

full_df <- mclapply(df_locs, train_cluster, train = train, mc.cores=cores)


######################## BIND AND SAVE DATA ##################################

save <- rbindlist(full_df)
slashes <- gregexpr("/",c.file)[[1]]
c.name <- substr(c.file, slashes[length(slashes)]+1, nchar(c.file)-4)
write.csv(save,paste0("<<<< FILEPATH REDACTED >>>>",c.name,".csv"), row.names=F)
