###### Perform age crosswalk
source(paste0(repo, "/collapse_functions/crosswalk/age_crosswalk.R"))

min_age_inclusive <- 0 ### Set the minimum age (inclusive) to which the data should be crosswalked
max_age_exclusive <- 5 ### Set the maximum age (exclusive) to which the data should be crosswalked; e.g., under-5 or [0-5) would have min_age_inclusive = 0 and max_age_exclusive = 5

#### Initial setup
temp_coverage_data <- coverage_data
temp_coverage_data$had_indicator <- temp_coverage_data$has_lri
age_crosswalk_setup()
tracking <- get_tracking_sheet(sheet_name="LRI Vetting Sheet")

#### Subset tracking to surveys with non-standard age ranges
surveys_non_standard <- unique(tracking[!(((Youngest_age) == as.character(min_age_inclusive)) & ((Oldest_age) == as.character(max_age_exclusive))),])

#### Subset to data with non-standard ages
non_standard_data <- temp_coverage_data[svy_id %in% surveys_non_standard$nid]
nids_to_adjust <- unique(non_standard_data$svy_id)

#### Check that non-standard surveys have numeric values in both fields
#View(surveys_non_standard[nid %in% non_standard_data$svy_id])
nids_to_drop <- list()
for (a in 1:nrow((surveys_non_standard[nid %in% non_standard_data$svy_id]))) {
  if (is.na(as.numeric(surveys_non_standard[nid %in% non_standard_data$svy_id][a]$Youngest_age)) | is.na(as.numeric(surveys_non_standard[nid %in% non_standard_data$svy_id][a]$Oldest_age)))  {
    print(paste0("nid ", surveys_non_standard[nid %in% non_standard_data$svy_id][a]$nid, " has non-numeric start and/or end age(s). Skipping for now, but you may want to investigate this and rerun the collapse code."))
    nids_to_drop <- c(nids_to_drop, surveys_non_standard[nid %in% non_standard_data$svy_id][a]$nid)
  }
}
surveys_non_standard <- surveys_non_standard[!(nid %in% nids_to_drop)]
non_standard_data <- temp_coverage_data[svy_id %in% surveys_non_standard$nid]
surveys_non_standard <- surveys_non_standard[(nid %in% non_standard_data$svy_id)]

#### Retrieve population and prevalence data from GBD
modelable_entity <- 1258 # lower respritory infections
GBD_pop_prev <- retrieve_GBD_pop_prev(non_standard_data, min_age_inclusive, max_age_exclusive, gbd_round_id=5, modelable_entity=modelable_entity)
non_standard_data <- GBD_pop_prev[["non_standard_data"]] 
sex_collapsed <- GBD_pop_prev[["sex_collapsed"]]
popNumbers <- GBD_pop_prev[["popNumbers"]]
age_metadata <- GBD_pop_prev[["age_metadata"]]
location_hierarchy <- GBD_pop_prev[["location_hierarchy"]]
age_ids <- GBD_pop_prev[["age_ids"]]

#### Summarize non-standard data by nid and year (don't assume that all data rows for a given nid have the same year)
nid_year_unique <- unique(non_standard_data[, c("svy_id", "int_year")])
nid_year_unique <- nid_year_unique[int_year > 1998]


#testing in parallel
run_cw_on_nid_year <- function(a, dt = temp_coverage_data){
  data_unadjusted <- copy(temp_coverage_data[svy_id == nid_year_unique[a, svy_id] & floor(int_year) == nid_year_unique[a, int_year]])
  adj_prevs <- age_crosswalk_by_nid_year(data_unadjusted=data_unadjusted, min_age_inclusive=min_age_inclusive, 
                            max_age_exclusive=max_age_exclusive, fractional_prev_input=TRUE, 
                            adjust_logit_space=FALSE)$adjusted_had_indicator
  data_unadjusted <- copy(temp_coverage_data[svy_id == nid_year_unique[a, svy_id] & floor(int_year) == nid_year_unique[a, int_year]])
  data_unadjusted$had_indicator <- adj_prevs
  return(data_unadjusted)
}
start_time <- Sys.time()
adj_coverage_data <- mclapply(1:nrow(nid_year_unique), run_cw_on_nid_year, dt = temp_coverage_data, mc.cores = cores)
end_time <- Sys.time()
end_time - start_time
adj_coverage_data <- rbindlist(adj_coverage_data)

adj_coverage_data$has_lri <- adj_coverage_data$had_indicator
coverage_data <- temp_coverage_data[(!svy_id %in% nid_year_unique$svy_id)]
coverage_data <- rbind(adj_coverage_data, coverage_data)
