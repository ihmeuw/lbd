###
### ALL ESPEN API CALLS TOGETHER 
###
### Written to easily use previously written API calls 


### All necessary packages and also the espen_api.R file location ### ----
# library(httr)
# library(jsonlite)
# library(lubridate)
# library(sf)
# library(geojsonsf)
# library(geojson)
# library(data.table)
# source(<<<< FILEPATH REDACTED >>>>)

# This is needed to read what the api returns p sure
options(stringsAsFactors = F)

### ALL CALLS TOGETHER ### ----

diseases <- c("oncho", "lf", "loa", "sch", "sth")

for (dis in diseases) {
  espenAPIData(dis, "sitelevel")
  espenAPIData(dis, "iu")
  espenAPIMap(dis, "iu", "endemicity", 2017)
}



### ALL DATA CALLS ### ----
# espenAPIData("oncho", "sitelevel")
# espenAPIData("oncho", "iu")
# espenAPIData("lf", "sitelevel")
# espenAPIData("lf", "iu")
# espenAPIData("loa", "sitelevel")
# espenAPIData("loa", "iu")
# espenAPIData("sch", "sitelevel")
# espenAPIData("sch", "iu")
# espenAPIData("sth", "sitelevel")
# espenAPIData("sth", "iu")

# No trachoma info
# espenAPIData("trachoma", "sitelevel")
# espenAPIData("trachoma", "iu")

### ALL MAP CALLS ### ----

# espenAPIMap("oncho", "iu", "endemicity", 2017)
# espenAPIMap("loa", "iu", "endemicity", 2017)
# espenAPIMap("sch", "iu", "endemicity", 2017)
# espenAPIMap("sth", "iu", "endemicity", 2017)
# espenAPIMap("lf", "iu", "endemicity", 2017)


### SHOULD NOT NEED THESE ###

# espenAPIMap("lf", "sitelevel", "tas")
# espenAPIMap("lf", "sitelevel", "sentinel_sites")
# espenAPIMap("lf", "sitelevel", "mapping_surveys")
# espenAPIMap("lf", "iu", "mda_pc_rounds", "projections")
# espenAPIMap("lf", "iu", "mda_pc_rounds", "therapeutic")
# espenAPIMap("lf", "iu", "mda_pc_rounds", "geographic")

# espenAPIMap("oncho", "sitelevel", "impact_assessment", "skin_biopsy")
# espenAPIMap("oncho", "sitelevel", "impact_assessment", "anti_ov16_test")
# espenAPIMap("oncho", "sitelevel", "impact_assessment", "nodule_palpation")
# espenAPIMap("oncho", "sitelevel", "mapping_surveys", "skin_biopsy")
# espenAPIMap("oncho", "sitelevel", "mapping_surveys", "anti_ov16_test")
# espenAPIMap("oncho", "sitelevel", "mapping_surveys", "nodule_palpation")
# espenAPIMap("oncho", "iu", "mda_pc_coverage", "geographic")
# espenAPIMap("oncho", "iu", "mda_pc_rounds", "geographic")
# espenAPIMap("oncho", "iu", "mda_pc_coverage", "therapeutic")
# espenAPIMap("oncho", "iu", "mda_pc_rounds", "therapeutic")

# espenAPIMap("loa", "sitelevel", "mapping_surveys", "ewh_questionnaire")
# espenAPIMap("loa", "sitelevel", "mapping_surveys", "blood_smear")

# espenAPIMap("sch", "sitelevel", "mapping_surveys", "all_species")
# espenAPIMap("sch", "sitelevel", "mapping_surveys", "s_haematobium")
# espenAPIMap("sch", "sitelevel", "mapping_surveys", "s_mansoni")
# espenAPIMap("sch", "iu", "mda_pc_coverage", "geographic_sac")
# espenAPIMap("sch", "iu", "mda_pc_rounds", "geographic_sac")
# espenAPIMap("sch", "iu", "mda_pc_coverage", "therapeutic_sac")
# espenAPIMap("sch", "iu", "mda_pc_rounds", "therapeutic_sac")
# espenAPIMap("sch", "iu", "mda_pc_coverage", "geographic_total")
# espenAPIMap("sch", "iu", "mda_pc_rounds", "geographic_total")
# espenAPIMap("sch", "iu", "mda_pc_coverage", "therapeutic_total")
# espenAPIMap("sch", "iu", "mda_pc_rounds", "therapeutic_total")

# espenAPIMap("sth", "sitelevel", "mapping_surveys", "all_species")
# espenAPIMap("sth", "sitelevel", "mapping_surveys", "ascaris")
# espenAPIMap("sth", "sitelevel", "mapping_surveys", "hookworms")
# espenAPIMap("sth", "sitelevel", "mapping_surveys", "trichuris")
# espenAPIMap("sth", "iu", "mda_pc_coverage", "geographic_sac")
# espenAPIMap("sth", "iu", "mda_pc_rounds", "geographic_sac")
# espenAPIMap("sth", "iu", "mda_pc_coverage", "therapeutic_sac")
# espenAPIMap("sth", "iu", "mda_pc_rounds", "therapeutic_sac")
# espenAPIMap("sth", "iu", "mda_pc_coverage", "geographic_total")
# espenAPIMap("sth", "iu", "mda_pc_rounds", "geographic_total")
# espenAPIMap("sth", "iu", "mda_pc_coverage", "therapeutic_total")
# espenAPIMap("sth", "iu", "mda_pc_rounds", "therapeutic_total")

