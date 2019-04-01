# HEADER ------------------------------------------------------------------
# 03_graph_coverage.R
# Purpose: Create data coverage plots
#**************************************************************************

# Note: run from 00_master.R, which sets up all directories & files
# Outputs: vaccine-specific plots (maps and scatter plots)

vaccines_to_graph <- c("dpt3_cov")
vaccine_titles <- c("DPT3")

for (i in length(vaccines_to_graph)) {
  
  # Define objects 
  vax <- vaccines_to_graph[i]
  vax_title <- vaccine_titles[i]
  vax_prefix <- str_match(vax, "^([a-z]*)[0-9]?_cov")[2]
  vax_dose <- str_match(vax, "^[a-z]*([0-9]?)_cov")[2] %>% as.numeric

  # Load data
  coverage_data <- readRDS(paste0(vaccine_cleaned_file_prefix, vax_prefix, ".rds"))

  # Prep data for graphing
  coverage_data <- coverage_data[, latitude := as.numeric(latitude)]
  coverage_data <- coverage_data[, longitude := as.numeric(longitude)]

  coverage_data[grep("/", source), source := "COUNTRY-SPECIFIC"]

  # First, drop year
  coverage_data[,year := NULL]

  ### Create vaccine variable
  
  # Figure out which doses are present in the data & grab max
  dose_vars <- names(coverage_data)[grepl(paste0(vax_prefix, "_dose"), names(coverage_data))]
  max_doses <- str_match(dose_vars, paste0(vax_prefix, "_dose_([0-9])"))[,2] %>%
                  as.numeric %>% max

  # Want to sum up this vaccine and all vaccines above it (e.g. DPT1 = those who rec'd 1,2, or 3 doses)
  sum_vars <- paste0(vax_prefix, "_dose_", vax_dose:max_doses)
  coverage_data[, tmp := rowSums(.SD), .SDcols = sum_vars]

  # Remove other vaccine dose variables
  coverage_data <- subset(coverage_data, select = !(names(coverage_data) %in% dose_vars))

  # Rename
  setnames(coverage_data, "tmp", vax)

  # Next, sum up over all variables but outcome & N
  coverage_data <- coverage_data[, list(N = sum(N), outcome = sum(get(vax))), 
                                  by = setdiff(names(coverage_data), c("N", vax))]

  # Restore names
  setnames(coverage_data, "outcome", vax)

  # Finally, convert coverage to percents
  coverage_data[, eval(vax) := get(vax) / N]

  # Truncate long names
  coverage_data[source == "GLOBAL_FUND_HOUSEHOLD_SURVEY", source := "GLOBAL_FUND"]
  coverage_data[source == "ARAB_LEAGUE_PAPFAM", source := "ARAB_LG_PAPFAM"]
  coverage_data[source == "MACRO_DHS_SP", source := "MACRO_DHS"]
  coverage_data[source == "WB_LSMS_ISA", source := "WB_LSMS"]

  if("original_year" %in% names(coverage_data)) {
     (setnames(coverage_data, "original_year", "svy_year") )
  }

 # for (reg in c("africa", "middle_east", "south_asia", "south_america", "central_america", "se_asia")) {
  for (reg in c("africa","south_asia","south_america","central_america","se_asia")) {

    # Make graphs
    # Note - using raw data (not point-poly processed)
    coverage_maps <- graph_data_coverage_values(df = coverage_data,
                                                var = vax,
                                                title = vax_title,
                                                year_min = '1997',
                                                year_max = '2016',
                                                year_var = 'svy_year',
                                                region = reg,
                                                cores = cores,
                                                indicator = vax,
                                                since_date = '2017-06-27',
                                                high_is_bad = FALSE,
                                                return_maps = TRUE,
                                                legend_title = "Vaccine \nCoverage",
                                                endemic_gauls = NULL,
                                                map_point_size = 0.8,
                                                fast_shapefiles = T,
                                                simplify_polys = T,
                                                tolerance = 0.01,
                                                save_on_share = F,
                                                base_font_size = 18,
                                                out_dir = "<<<< FILEPATH REDACTED >>>>",
                                                prep_shiny = F,
                                                color_scheme = "classic")
  }
}
