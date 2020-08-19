##########################################################################################
####################  LF age and diagnostic crosswalk data preparation  ##################
##########################################################################################

##################################################################################
###### I. Setup
##################################################################################

user <- Sys.info()[["user"]] ## Get current user name
core_repo <- paste0(<<<< FILEPATH REDACTED >>>>)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- sprintf(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))
package_list <- c(package_list, "sf")
path <- paste0(<<<< FILEPATH REDACTED >>>>)

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

library(jsonlite)
mbg_setup(package_list = package_list, repos = core_repo)

library(ggthemes)
library(fasterize)
library(tidyverse)
library(fda)
library(subplex)
library(faraway)
library(BayesianTools)
library(spdep)
library(boot)

#### Load shared db functions
source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)

##################################################################################
###### II. Derive sex-collapsed raw data set
##################################################################################

### Prepare original sex-collapsed data file
data1 <- fread(<<<< FILEPATH REDACTED >>>>)
data1 <- data1[!is.na(uniqueID)]
data1$cases <- round(data1$cases, 0)
data1$sample_size <- round(data1$sample_size, 0)
setnames(data1, "note_modler", "note_modeler")

### Correcting extraction errors
data1[nid == 136597 & age_start == 10 & group == 13 & diagnostic == "mf", cases := 3]
data1[nid == 136449 & age_start == 6 & group == 407 & diagnostic == "mf", cases := 2]
data1[nid == 136449 & age_start == 16 & group == 407 & diagnostic == "mf", cases := 1]

### Compute missing cases or sample size
data1[is.na(cases) & !is.na(mean) & !is.na(sample_size), cases := round(mean * sample_size, 0)]
data1[!is.na(cases) & !is.na(mean) & is.na(sample_size), sample_size := round(cases / mean, 0)]

data1_subset <- data1[, c("nid", "location_name", "location_id", "ihme_loc_id", "latitude", "longitude", "site_memo", "case_diagnostics", "case_definition", "case_name", "sex", "year_start", "year_end", "age_start", "age_end", "diagnostic", "group", "cases", "sample_size", "note_modeler")]

### Prepare new extractions (2017-2018)
extraction <- fread(<<<< FILEPATH REDACTED >>>>)
extraction3 <- fread(<<<< FILEPATH REDACTED >>>>)
extraction4 <- fread(<<<< FILEPATH REDACTED >>>>)
extraction5 <- fread(<<<< FILEPATH REDACTED >>>>)
extraction6 <- fread(<<<< FILEPATH REDACTED >>>>)
extraction7 <- fread(<<<< FILEPATH REDACTED >>>>)

setnames(extraction, c("lat", "long"), c("latitude", "longitude"))
setnames(extraction3, c("lat", "long"), c("latitude", "longitude"))
extraction$row_num <- paste0("ext1_", 1:nrow(extraction))
setnames(extraction3, c("mean prevalence (decimal)"), c("mean prevalence (%)"))
extraction3 <- extraction3[c(2:nrow(extraction3)), ]
extraction3$latitude <- as.numeric(extraction3$latitude)
extraction3$longitude <- as.numeric(extraction3$longitude)
extraction3$row_num <- paste0("ext3_", 1:nrow(extraction3))
setnames(extraction4, c("mean prevalence (decimal)", "lat", "long"), c("mean prevalence (%)", "latitude", "longitude"))
extraction4 <- extraction4[c(2:nrow(extraction4)), ]
extraction4$latitude <- as.numeric(extraction4$latitude)
extraction4$longitude <- as.numeric(extraction4$longitude)
extraction4$row_num <- paste0("ext4_", 1:nrow(extraction4))
extraction <- rbind(extraction, extraction3, use.names = T, fill = T)
extraction[, diagnostic := NA]
extraction <- rbind(extraction, extraction4, use.names = T, fill = T)
setnames(extraction5, c("mean prevalence (decimal)", "lat", "long"), c("mean prevalence (%)", "latitude", "longitude"))
extraction5$latitude <- as.numeric(extraction5$latitude)
extraction5$longitude <- as.numeric(extraction5$longitude)
extraction5$row_num <- paste0("ext5_", 1:nrow(extraction5))
setnames(extraction6, c("country"), c("location_id"))
extraction6$latitude <- as.numeric(extraction6$latitude)
extraction6$longitude <- as.numeric(extraction6$longitude)
extraction6$row_num <- paste0("ext6_", 1:nrow(extraction6))
setnames(extraction7, c("mean prevalence (decimal)", "lat", "long"), c("mean prevalence (%)", "latitude", "longitude"))
extraction7$latitude <- as.numeric(extraction7$latitude)
extraction7$longitude <- as.numeric(extraction7$longitude)
extraction7$row_num <- paste0("ext7_", 1:nrow(extraction7))
extraction5 <- extraction5[!is.na(nid)]
extraction7 <- extraction7[!(title == "")]
extraction <- rbind(extraction, extraction5, use.names = T, fill = T)
extraction <- rbind(extraction, extraction6, use.names = T, fill = T)
extraction <- rbind(extraction, extraction7, use.names = T, fill = T)
extraction <- extraction[c(2:nrow(extraction))]
extraction <- extraction[!is.na(row_num)]

data_new <- copy(extraction)
setnames(data_new, "mean prevalence (%)", "mean_prev")
data_new$cases <- round(as.numeric(data_new$cases), 0)
data_new$sample_size <- round(as.numeric(data_new$sample_size), 0)
data_new$mean_prev <- as.numeric(data_new$mean_prev) / 100
data_new[is.na(cases) & !is.na(mean_prev), cases := round(mean_prev * sample_size, 0)]
data_new$diagnostic <- NA
data_new$note_modeler <- NA
data_new$sex <- tolower(data_new$sex)
data_new$had_lf_poly <- as.numeric(data_new$had_lf_poly)
data_new$N <- as.numeric(data_new$N)
data_new[is.na(cases) & !is.na(had_lf_poly), cases := had_lf_poly]
data_new[is.na(sample_size) & !is.na(N), sample_size := N]

data_new <- data_new[!is.na(sample_size) & !is.na(cases) & !(cases > sample_size)]

data_new_subset <- data_new[, c("nid", "location_name", "location_id", "ihme_loc_id", "latitude", "longitude", "site_memo", "case_diagnostics", "case_definition", "case_name", "sex", "year_start", "year_end", "age_start", "age_end", "diagnostic", "group", "cases", "sample_size", "note_modeler", "shape_type", "poly_reference", "poly_id_field_name")]
data_new_subset <- data_new_subset[!(case_name %in% c("hydrocele", "lymphoedema"))] # Drop morbidity data

### Sex collapse data2_subset
sex_counts <- aggregate(sex ~ nid, data = unique(data_new_subset[sex %in% c("male", "female"), c("sex", "nid")]), length)
sex_agg <- merge(
  aggregate(cases ~ nid + site_memo + case_diagnostics + case_definition + case_name + year_start + year_end + age_start + age_end, data_new_subset[sex %in% c("male", "female")], sum),
  aggregate(sample_size ~ nid + site_memo + case_diagnostics + case_definition + case_name + year_start + year_end + age_start + age_end, data_new_subset[sex %in% c("male", "female")], sum)
)
sex_agg_merged <- merge(data_new_subset, sex_agg, by = c("nid", "site_memo", "case_diagnostics", "case_definition", "case_name", "year_start", "year_end", "age_start", "age_end"), all.x = TRUE)
sex_male <- sex_agg_merged[sex %in% c("male")] # Select one row of each pair
sex_male[, cases.x := NULL]
sex_male[, sample_size.x := NULL]
sex_male$sex <- "both"
setnames(sex_male, "cases.y", "cases")
setnames(sex_male, "sample_size.y", "sample_size")
data_new_sex_collapsed <- rbind(data_new_subset[sex == "both"], sex_male[, c("nid", "location_name", "location_id", "ihme_loc_id", "latitude", "longitude", "site_memo", "case_diagnostics", "case_definition", "case_name", "sex", "year_start", "year_end", "age_start", "age_end", "diagnostic", "group", "cases", "sample_size", "note_modeler", "shape_type", "poly_reference", "poly_id_field_name")])

### Append data files
data1_subset[, c("shape_type", "poly_reference", "poly_id_field_name") := list(NA, NA, NA)]
sex_collapsed_all <- rbind(data1_subset, data_new_sex_collapsed)

##################################################################################
###### III. Assemble training data set for crosswalk model
##################################################################################
### Clean up ihme_loc_id
sex_collapsed_all[length(ihme_loc_id) > 3, ihme_loc_id := substr(ihme_loc_id, 1, 3) ]

sex_collapsed_all[grep("Sri Lanka", location_name), ihme_loc_id := "LKA"]
sex_collapsed_all[grep("Papua New Guinea", location_name), ihme_loc_id := "PNG"]
sex_collapsed_all[grep("Haiti", location_name), ihme_loc_id := "HTI"]
sex_collapsed_all[grep("Brazil", location_name), ihme_loc_id := "BRA"]
sex_collapsed_all[grep("Kenya", location_name), ihme_loc_id := "KEN"]
sex_collapsed_all[grep("Malawi", location_name), ihme_loc_id := "MWI"]
sex_collapsed_all[grep("Nigeria", location_name), ihme_loc_id := "NGA"]
sex_collapsed_all[grep("Thailand", location_name), ihme_loc_id := "THA"]
sex_collapsed_all[grep("Vanuatu", location_name), ihme_loc_id := "VUT"]
sex_collapsed_all[grep("Republic of the Congo", location_name), ihme_loc_id := "COD"]
sex_collapsed_all[grep("Ghana", location_name), ihme_loc_id := "GHA"]
sex_collapsed_all[grep("Indonesia", location_name), ihme_loc_id := "IDN"]
sex_collapsed_all[grep("Sierra Leone", location_name), ihme_loc_id := "SLE"]
sex_collapsed_all[grep("Madagascar", location_name), ihme_loc_id := "MDG"]
sex_collapsed_all[grep("India", location_name), ihme_loc_id := "IND"]
sex_collapsed_all[grep("American Samoa", location_name), ihme_loc_id := "ASM"]
sex_collapsed_all[grep("Cook Islands", location_name), ihme_loc_id := "COK"]
sex_collapsed_all[grep("Togo", location_name), ihme_loc_id := "TGO"]
sex_collapsed_all[grep("Mali", location_name), ihme_loc_id := "MLI"]
sex_collapsed_all[grep("Uganda", location_name), ihme_loc_id := "UGA"]
sex_collapsed_all[grep("Tanzania", location_name), ihme_loc_id := "TZA"]
sex_collapsed_all[grep("Equatorial Guinea", location_name), ihme_loc_id := "GNQ"]
sex_collapsed_all[grep("Myanmar", location_name), ihme_loc_id := "MMR"]
sex_collapsed_all[grep("The Gambia", location_name), ihme_loc_id := "GMB"]
sex_collapsed_all[grep("Egypt", location_name), ihme_loc_id := "EGY"]
sex_collapsed_all[location_name == "Federal States of Micronesia", location_name := "Federated States of Micronesia"]
sex_collapsed_all[grep("Federated States of Micronesia", location_name), ihme_loc_id := "FSM"]
sex_collapsed_all[grep("Cameroon", location_name), ihme_loc_id := "CMR"]
sex_collapsed_all[grep("Wallis and Futuna", location_name), ihme_loc_id := "WLF"]
sex_collapsed_all[grep("Niue", location_name), ihme_loc_id := "NIU"]
sex_collapsed_all[grep("New Caledonia", location_name), ihme_loc_id := "NCL"]
sex_collapsed_all[grep("Dominican Republic", location_name), ihme_loc_id := "DOM"]
sex_collapsed_all[grep("Benin", location_name), ihme_loc_id := "BEN"]
sex_collapsed_all[grep("Yemen", location_name), ihme_loc_id := "YEM"]
sex_collapsed_all[ihme_loc_id == "Fij", ihme_loc_id := "FJI"]
sex_collapsed_all[ihme_loc_id == "Hai", ihme_loc_id := "HTI"]
sex_collapsed_all[location_id == "IND", ihme_loc_id := "IND"]
sex_collapsed_all[location_id == "IDN", ihme_loc_id := "IDN"]
sex_collapsed_all[ihme_loc_id == "Pap", ihme_loc_id := "PNG"]

unique(sex_collapsed_all[, c("location_name", "ihme_loc_id")])
sex_collapsed_all$ihme_loc_id <- toupper(sex_collapsed_all$ihme_loc_id)

sex_collapsed_all[is.na(site_memo), site_memo := as.character(nid)]

### Clean up diagnostic type
sex_collapsed_all[nid == 288836, diagnostic := "ict"]
sex_collapsed_all[nid == 288809, diagnostic := "mf"]
sex_collapsed_all[nid == 159048, diagnostic := "mf"]
sex_collapsed_all[nid == 288891, diagnostic := "mf"]
sex_collapsed_all[nid == 136493, diagnostic := "ict"]
sex_collapsed_all[nid == 288817, diagnostic := "mf"]
sex_collapsed_all[nid == 136521 & case_definition == "ICT test", diagnostic := "ict"]
sex_collapsed_all[nid == 136521 & case_definition == "Mf, thick blood film 60ul", diagnostic := "mf"]
sex_collapsed_all[nid == 147752, diagnostic := case_name]
sex_collapsed_all[nid == 136548 & cases %in% c(41, 103, 71, 112), diagnostic := "ict"]
sex_collapsed_all[nid == 136548 & cases %in% c(9, 40, 18, 49), diagnostic := "mf"]
sex_collapsed_all[nid == 136449, diagnostic := case_name]
sex_collapsed_all[nid == 136449 & is.na(diagnostic) & age_start == 2 & sample_size == 7, diagnostic := "ict"]
sex_collapsed_all[nid == 136449 & is.na(diagnostic) & note_modeler == 100008, diagnostic := "ict"]
sex_collapsed_all[nid == 136449 & is.na(diagnostic) & note_modeler == 100022, diagnostic := "mf"]
sex_collapsed_all[nid == 147658 & note_modeler == 100012, diagnostic := "ict"]
sex_collapsed_all[nid == 147658 & note_modeler == 100026, diagnostic := "mf"]
sex_collapsed_all[nid == 143009, diagnostic := "ict"]
sex_collapsed_all[nid == 222125, diagnostic := case_name]
sex_collapsed_all[nid == 288878 & case_diagnostics == "ICT test postive", diagnostic := "ict"]
sex_collapsed_all[nid == 288878 & case_diagnostics == "FTS test", diagnostic := "fts"]
sex_collapsed_all[nid == 293289, diagnostic := "ict"]
sex_collapsed_all[nid == 136588, diagnostic := "mf"]
sex_collapsed_all[nid == 288844, diagnostic := "mf"]
sex_collapsed_all[nid == 222127, diagnostic := "mf"]
sex_collapsed_all[nid == 222538, diagnostic := "ict"]
sex_collapsed_all[case_diagnostics == "ICT", diagnostic := "ict"]
sex_collapsed_all[case_diagnostics == "FTS", diagnostic := "fts"]
sex_collapsed_all[case_diagnostics == "blood smear", diagnostic := "mf"]
sex_collapsed_all[case_definition == "ICT positive for W bancroft CFA", diagnostic := "ict"]
sex_collapsed_all[case_definition == "Ag prevalence by ICT", diagnostic := "ict"]
sex_collapsed_all[case_diagnostics == "LF IgG4 RDT", diagnostic := "ab"]
sex_collapsed_all[nid == 389542, diagnostic := "ict"]
sex_collapsed_all[nid == 136427, diagnostic := "mf"]
sex_collapsed_all[nid == 389545, diagnostic := "ab"]
sex_collapsed_all[nid == 389547 & case_name == "mf night blood smear", diagnostic := "mf"]
sex_collapsed_all[nid == 389547 & case_definition == "positive brugia rapid test", diagnostic := "br"]
sex_collapsed_all[nid == 389547 & case_name == "IgG4 antibodies for B. malayi", diagnostic := "br"]
sex_collapsed_all[case_name == "mf", diagnostic := "mf"]
sex_collapsed_all[case_name %in% c("ict", "ICT"), diagnostic := "ict"]
sex_collapsed_all[case_diagnostics == "rapid format card test (BinazNOW Filariasis, Alere. Inc)", diagnostic := "ict"]
sex_collapsed_all[case_diagnostics == "ELISA", diagnostic := "ab"]
sex_collapsed_all[case_diagnostics == "ICT test (100 microliters)", diagnostic := "ict"]
sex_collapsed_all[case_name == "Og4C3 Ag > 32 units", diagnostic := "ab"]
sex_collapsed_all[case_diagnostics == "InBios Wb123 ELISA", diagnostic := "ab"]
sex_collapsed_all[case_diagnostics == "ELISA tests (CDC in house version)", diagnostic := "ab"]
sex_collapsed_all[case_definition == "elephantiasis", diagnostic := "elephantiasis"]
sex_collapsed_all[case_diagnostics == "Binax now filariasis antigen test", diagnostic := "ict"]
sex_collapsed_all[diagnostic == "immunochromatographic test (ict) results", diagnostic := "ict"]
sex_collapsed_all[diagnostic == "microfilaria (mf) counts", diagnostic := "mf"]

### Revisions based on closer examination of studies in final crosswalk data set
sex_collapsed_all[nid == 271326, diagnostic := "CAg (LIPS)"]
sex_collapsed_all[nid == 136427, diagnostic := "mf"]
sex_collapsed_all <- sex_collapsed_all[!(nid == 293113)]

### Revisions based on inclusion of new extractions
# Drop morbidity data
sex_collapsed_all <- sex_collapsed_all[!(diagnostic %in% c("elephantiasis", "hydrocele", "lymphoedema", "hydrocoele", "adl", "filarial_fever", "other", "scrotal elephantiasis"))]
sex_collapsed_all <- sex_collapsed_all[!(case_diagnostics %in% c("Filarial disease", "Clinical filarias"))]

sex_collapsed_all[case_diagnostics %in% c("night blood smear", "thick blood smear"), diagnostic := "mf"]
sex_collapsed_all[diagnostic == "polymerase chain reaction (pcr) assay", diagnostic := "pcr"]
sex_collapsed_all[case_diagnostics %in% c("Brugia Rapid test"), diagnostic := "br"]
sex_collapsed_all[case_name %in% c("Wb123", "wb123 antibody"), diagnostic := "ab"]
sex_collapsed_all[case_name %in% c("mf positive"), diagnostic := "mf"]
sex_collapsed_all[case_name %in% c("FTS positive for cfa", "FTS positive"), diagnostic := "fts"]
sex_collapsed_all[case_name %in% c("ICT positive"), diagnostic := "ict"]
sex_collapsed_all[case_name %in% c("mf positive", "mf+"), diagnostic := "mf"]
sex_collapsed_all[case_definition %in% c("w bancrofti antigen presence") &  nid == 389545, diagnostic := "elisa"]
sex_collapsed_all[case_diagnostics %in% c("ELISA"), diagnostic := "elisa"]
sex_collapsed_all[case_name %in% c("Bm 14 positive", "Bm14 positive"), diagnostic := "elisa"]
sex_collapsed_all[case_diagnostics %in% c("ICT, 100 microliters of finger prick blood", "ICT, BinaxNow Filariasis", "ICT, BinaxNow Filariasis "), diagnostic := "ict"]
sex_collapsed_all[case_diagnostics %in% c("Brugia Rapid"), diagnostic := "br"]
sex_collapsed_all[case_definition %in% c("mf positive", "mf prevalence rate"), diagnostic := "mf"]
sex_collapsed_all[case_diagnostics %in% c("thick night blood smear"), diagnostic := "mf"]
sex_collapsed_all[case_name %in% c("antibody response to the W bancrofti antigen Wb123"), diagnostic := "ab"]
sex_collapsed_all[case_name %in% c("microfilaraemia prevalence"), diagnostic := "mf"]

### Adjust location for Cook Islands observations
sex_collapsed_all[nid == 271326, "location_name" := "Cook Islands"]
sex_collapsed_all[nid == 271326, "ihme_loc_id" := "COK"]

sex_collapsed_all$location_id <- as.integer(sex_collapsed_all$location_id)

### Standardize on country-level location ids
sex_collapsed_all[ihme_loc_id == "BEN", location_id := 200]
sex_collapsed_all[ihme_loc_id == "BGD", location_id := 161]
sex_collapsed_all[ihme_loc_id == "BRA", location_id := 135]
sex_collapsed_all[ihme_loc_id == "CMR", location_id := 202]
sex_collapsed_all[ihme_loc_id == "COD", location_id := 170]
sex_collapsed_all[ihme_loc_id == "COK", location_id := 21] ### Use population distribution for Oceania, as GBD doesn't break out this country
sex_collapsed_all[ihme_loc_id == "DOM", location_id := 111]
sex_collapsed_all[ihme_loc_id == "EGY", location_id := 141]
sex_collapsed_all[ihme_loc_id == "GHA", location_id := 207]
sex_collapsed_all[ihme_loc_id == "HTI", location_id := 114]
sex_collapsed_all[ihme_loc_id == "IDN", location_id := 11]
sex_collapsed_all[ihme_loc_id == "IND", location_id := 163]
sex_collapsed_all[ihme_loc_id == "KEN", location_id := 180]
sex_collapsed_all[ihme_loc_id == "LKA", location_id := 17]
sex_collapsed_all[ihme_loc_id == "MLI", location_id := 211]
sex_collapsed_all[ihme_loc_id == "MWI", location_id := 182]
sex_collapsed_all[ihme_loc_id == "NGA", location_id := 214]
sex_collapsed_all[ihme_loc_id == "PHL", location_id := 16]
sex_collapsed_all[ihme_loc_id == "PNG", location_id := 26]
sex_collapsed_all[ihme_loc_id == "SLE", location_id := 217]
sex_collapsed_all[ihme_loc_id == "THA", location_id := 18]
sex_collapsed_all[ihme_loc_id == "TZA", location_id := 189]
sex_collapsed_all[ihme_loc_id == "UGA", location_id := 190]
sex_collapsed_all[ihme_loc_id == "VUT", location_id := 30]
sex_collapsed_all[ihme_loc_id == "WSM", location_id := 425]
sex_collapsed_all[ihme_loc_id == "PHL", location_id := 16]
sex_collapsed_all[ihme_loc_id == "ASM", location_id := 298]
sex_collapsed_all[ihme_loc_id == "GNQ", location_id := 172]
sex_collapsed_all[ihme_loc_id == "TGO", location_id := 218]
sex_collapsed_all[ihme_loc_id == "MMR", location_id := 15]
sex_collapsed_all[ihme_loc_id == "GMB", location_id := 206]
sex_collapsed_all[ihme_loc_id == "FSM", location_id := 25]
sex_collapsed_all[ihme_loc_id == "MDG", location_id := 181]
sex_collapsed_all[ihme_loc_id == "WLF", location_id := 21] ### Use population distribution for Oceania, as GBD doesn't break out this country
sex_collapsed_all[ihme_loc_id == "NIU", location_id := 21] ### Use population distribution for Oceania, as GBD doesn't break out this country
sex_collapsed_all[ihme_loc_id == "NCL", location_id := 21] ### Use population distribution for Oceania, as GBD doesn't break out this country
sex_collapsed_all[ihme_loc_id == "FJI", location_id := 22]
sex_collapsed_all[ihme_loc_id == "YEM", location_id := 157]
sex_collapsed_all[ihme_loc_id == "SEN", location_id := 216]
sex_collapsed_all[ihme_loc_id == "TON", location_id := 29]

sex_collapsed_all <- sex_collapsed_all[diagnostic %in% c("ict", "mf", "fts", "br"), ]

sex_collapsed_all <- sex_collapsed_all[!is.na(cases) & !is.na(sample_size)]

### Year adjustment - under assumption
sex_collapsed_all$year_end <- as.numeric(as.character(sex_collapsed_all$year_end))
sex_collapsed_all$year_start <- as.numeric(as.character(sex_collapsed_all$year_start))

sex_collapsed_all[, year := year_end - year_start]

sex_collapsed_all[year == 0 | year == 1, year := year_start]
sex_collapsed_all[year == 2 | year == 3, year := year_start + 1]
sex_collapsed_all[year == 4, year := year_start + 2]

sex_collapsed_all[, uid_crosswalk := 1:nrow(sex_collapsed_all)]
sex_collapsed_all[, cohort := .GRP, by = list(nid, site_memo, year_start, year_end, group)]

### Pull population totals
loc_ids <- unique(data.table(
  country = sex_collapsed_all$ihme_loc_id,
  id = sex_collapsed_all$location_id,
  year = sex_collapsed_all$year
))

age_ids <- c(49:142) # Age 1 = ID 49, Age 94 = ID 142; get_population does not return single_year_age for ages 95-99

popNumbers <- data.table(get_population(age_group_id = age_ids, location_id = unique(loc_ids$id), year_id = unique(loc_ids$year), sex_id = 3, single_year_age = T, gbd_round = 6, decomp_step = "iterative"))
popNumbers_under_1 <- data.table(get_population(age_group_id = 28, location_id = unique(loc_ids$id), year_id = unique(loc_ids$year), sex_id = 3, single_year_age = F, gbd_round = 6, decomp_step = "iterative"))

popNumbers_concat <- rbind(popNumbers, popNumbers_under_1)

### Rounding numbers
sex_collapsed_all$sample_size <- ceiling(sex_collapsed_all$sample_size)
sex_collapsed_all$cases <- ceiling(sex_collapsed_all$cases)
popNumbers_concat$population <- ceiling(popNumbers_concat$population)

## Get and combine with age group names
ages <- get_ids("age_group")
pop_age_merged <- merge(popNumbers_concat, ages, by = "age_group_id")

## Calculate total population (0-94 years of age)
pop_agg <- aggregate(population ~ location_id + year_id, data = pop_age_merged, sum)
colnames(pop_agg)[3] <- "total_population"
pop_age_merged <- merge(pop_age_merged, pop_agg, by = c("location_id", "year_id"))

## Calculate % of pop. in each age year
pop_age_merged$age_prop <- pop_age_merged$population / pop_age_merged$total_population

## Convert to wide format
pop_age_merged_wide <- reshape(subset(pop_age_merged, select = -c(age_group_id, sex_id, run_id, population, total_population)), idvar = c("location_id", "year_id"), timevar = "age_group_name", direction = "wide")
colnames(pop_age_merged_wide)[2] <- "year"

## Merge age dist. with prev data
data_final <- merge(sex_collapsed_all, pop_age_merged_wide, by = c("location_id", "year"))

## Final clean-ups
setnames(data_final, "age_prop.<1 year", "age_prop.0")
data_final$age_start <- as.integer(data_final$age_start)
data_final$age_end <- as.integer(data_final$age_end)
data_final[age_start > 94, age_start := 94]
data_final[age_end > 94, age_end := 94]
data_final <- data_final[!is.na(sample_size)]
data_final$lat <- as.numeric(data_final$lat)

## Post-hoc corrections identified by review of preliminary crosswalk models
data_final <- data_final[!((nid == 135493) & (age_start == 6) & age_end == 15)]
data_final <- data_final[!((nid == 135493) & (age_start == 15) & age_end == 94)]
data_final <- data_final[!((nid == 136548) & (site_memo == "Nchacha") & (diagnostic == "ict") & (age_start == 1) & age_end == 14)]
data_final <- data_final[!((nid == 136548) & (site_memo == "Nchacha") & (diagnostic == "ict") & (age_start == 15) & age_end == 94)]
data_final <- data_final[!((nid == 136548) & (site_memo == "Belo") & (diagnostic == "ict") & (age_start == 1) & age_end == 14)]
data_final <- data_final[!((nid == 136548) & (site_memo == "Belo") & (diagnostic == "ict") & (age_start == 15) & age_end == 94)]
data_final <- data_final[!((nid == 136548) & (site_memo == "136548"))]
data_final[((nid == 293289) & (note_modeler == 100016)), diagnostic := "mf"]
data_final <- data_final[!((nid == 389386) & (sample_size %in% c(1143, 2496)))]
data_final <- data_final[!((nid == 147752) & (location_name == "Kilifi|KEN_35630"))]
data_final <- data_final[!((nid == 136493) & (age_start == 15) & (age_end == 94))]
data_final <- data_final[!((nid == 136493) & (age_start == 1) & (age_end == 14))]
data_final <- data_final[!((nid == 389700) & (age_start == 2) & age_end == 94)]
data_final <- data_final[!((nid == 136493) & (age_start == 6) & age_end == 15)]

## Drop all cohorts with 0 total cases
data_final$cases <- ceiling(data_final$cases)
data_final$sample_size <- ceiling(data_final$sample_size)
total_cases_by_cohort <- as.data.table(aggregate(cases ~ cohort, data_final, sum))
data_final <- data_final[cohort %in% total_cases_by_cohort[cases != 0, cohort]]

## Save crosswalk training data set to file
write.csv(data_final, file=<<<< FILEPATH REDACTED >>>>, row.names=F)
