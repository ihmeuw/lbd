#########################################
# Scatter Plots for Data Quality Checks (comparisons with STATCompiler reports)
# Created by Krista Steuben 10/24/2017
# Updated by Kate Wilson 12/21/2017 for male circumcision
# updated by Michael Cork 05/042019
#########################################
rm(list=ls())

# Load required packages
libs <- c("dplyr", "gdata", "openxlsx", "readxl", "ggrepel")
sapply(libs, require, character.only = T)

topic <- "male_circumcision"
core_repo  <- "/lbd_core/"
indic_repo <- "/lbd_hiv/"
source(paste0(indic_repo, "mbg/functions/collapse_functions.r"))
geomatched_version <- most_recent_date("<<<< FILEPATH REDACTED >>>>")

micro_data <- readRDS("<<<< FILEPATH REDACTED >>>>")


# micro_zaf <- 
#   micro_data %>% filter(nid == 313076)
# 
# micro_zaf %>%
#   filter(pweight > 0) %>%
#   mutate(age_group = case_when(between(age_year, 15, 19) ~ 1,
#                                between(age_year, 20, 24) ~ 2,
#                                between(age_year, 25, 49) ~ 3,
#                                age_year > 49 ~ 4,
#                                TRUE ~ NA_real_)) %>%
#   group_by(age_group) %>%
#   summarize(n = n(),
#             pct = 100*weighted.mean(male_circumcision, pweight, na.rm = T))
# 
# 
# micro_zaf %>%
#   filter(pweight > 0) %>%
#   mutate(age_group = case_when(between(age_year, 15, 19) ~ 1,
#                                between(age_year, 20, 24) ~ 2,
#                                between(age_year, 25, 49) ~ 3,
#                                age_year > 49 ~ 4,
#                                TRUE ~ NA_real_)) %>%
#   count(age_group, male_circumcision) %>% 
#   group_by(age_group) %>% 
#   mutate(total = n / sum(n))
# 
# micro_MWI <- 
#   micro_data %>% 
#   filter(nid == 287629)
# 
# micro_MWI %>% 
#   filter(pweight > 0) %>% 
#   count(male_circumcision) %>% 
#   mutate(n = n / sum(n))
# 
# micro_MWI %>% 
#   filter(pweight > 0) %>% 
#   summarize(mc = mean(male_circumcision, na.rm = T),
#             n = n())
# 
# 
# micro_zaf %>% 
#   filter(between(age_year, 20, 24)) %>%
#   filter(!is.na(male_circumcision)) %>%
#   summarize(mc = mean(male_circumcision))

# Pull in table of report data
stat_compile_file <- '/lbd_hiv/data/covariate_surveys/microdata/male_circumcision/'
stat_report <- 
  read.xlsx(paste0(stat_compile_file, "4_statcompiler_mc_indicators.xlsx")) %>% 
  mutate(start_age = 15, end_age = 49)

# transform statcompiler download in db form to get sample sizes (weighted) and bind to table with nids
sc_db <- read.xlsx(paste0(stat_compile_file, "4_statcompiler_mc_db.xlsx"))
  
# drop all rows except those with weighted sample sizes
sc_db <- sc_db[grepl("Number*", sc_db$Indicator) & !grepl("*unweighted*", sc_db$Indicator), ]

# merge sc_db and sample
stat_report <- merge(stat_report, sc_db, by = c("Country.Name", "Survey.Name"))

# read in map of country names to IHME country codes and bind codes to sc_db
stat_report <- 
  read_xlsx("<<<< FILEPATH REDACTED >>>>") %>% 
  right_join(stat_report, by = "Country.Name")

# replace Survey Name with survey series to match extracted data
stat_report$survey_name <- ifelse(grepl("*DHS*", stat_report$Survey.Name), "MACRO_DHS", "MACRO_AIS")

# drop unnecessary rows & rename remaining 
stat_report <- 
  stat_report %>% 
  dplyr::select(survey_name, year = Survey.Year,
         sample_size = Value, country, mc, start_age, end_age) %>% 
  mutate(nid = NA)

# Read in non statcompiler report data and select the rows corresponding to sample in the same order
non_stat_report <- 
  read.xlsx(paste0(stat_compile_file, "4_non_statcompiler_report_data_mc.xlsx")) %>% 
  dplyr::select(survey_name, year, sample_size, country, mc, start_age, end_age, nid)


# bind the non statcompiler data to the statcompiler data and drop any rows that don't have a sample size (i.e., surveys that didn't report mc data)
report_data <- rbind(stat_report, non_stat_report) 
report_data <- filter(report_data, !is.na(mc))

# Only keep surveys that were extracted
extracted_surveys <-
  micro_data %>%
  dplyr::count(survey_name, year, country, nid) %>%
  dplyr::mutate(year = as.numeric(year)) %>%
  dplyr::select(-n)

# Make sure none are being dropped becuase of faulty match
excluded <- 
  report_data %>% 
  anti_join(extracted_surveys, by = c("year", "country"))

if (excluded %>% semi_join(extracted_surveys, by = c("year", "country")) %>% nrow() != 0){
  warning("Make sure that survey names are correct")
}

# Only reports that were extracted and appear in the microdata
report_data <- 
  report_data %>% 
  semi_join(extracted_surveys, by = c("survey_name", "year", "country")) %>% 
  dplyr::select(-nid) %>% 
  left_join(extracted_surveys, by = c("survey_name", "year", "country")) 

# Finding surveys by survey_series, country, year, & survey_module (should be enough to uniquely id) and means I don't have to collapse all by nid and mod only to merge nid onto stat_report
# Calculate the sample size and prevalence rates for each var of interest for each extracted survey
extracted_data <-
  rbindlist(apply(report_data, 1, function(report) {
    # Subset microdata to survey
    survey <-
      micro_data %>%
      filter(survey_name == report[1],
             country == report[4],
             year == report[2] | end_year == report[2],
             pweight > 0)
    # Subset to correct ages if specified
    if (!is.na(report[6])) survey <- filter(survey, age_year >= report[6])
    if (!is.na(report[7])) survey <- filter(survey, age_year <= report[7])
    if (unique(survey$nid) == 27987) survey <- filter(survey, psu <= 105)
    survey <- data.table(survey)
    survey <-
      survey[, .(mc_micro = 100*weighted.mean(male_circumcision, pweight, na.rm = T),
                 N_micro = .N,
                 year = as.numeric(min(year, na.rm = T))),
      by = .(nid, country)]
  }))

# join extracted data with reports
compare_dataset <- 
  report_data %>% 
  dplyr::rename(mc_report = mc, N_report = sample_size) %>% 
  left_join(extracted_data, by = c("year", "nid", "country"))


nids <- c(218565, 287630, 287629, 415531, 327591)
compare_dataset %>% 
  filter(nid %in% c(218565, 287630, 287629, 415531, 327591))

# micro_data %>% 
#   filter(nid == 313076, between(ages, 15, 19)) %>% 

# Make plots comparing reported and extracted circumcision prevalence and sample size
gg_mc <- 
  compare_dataset %>% 
  ggplot(aes(x = mc_report, y = mc_micro)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_text_repel(data = filter(compare_dataset, (abs(mc_report - mc_micro) / mc_report) > 0.05),
                  aes(label = nid)) + 
  labs(title = "Comparing prevalence between report data and extracted microdata", 
       x = "Report data", y = "Extracted microdata") + 
  coord_equal() + 
  theme_bw()

gg_N <- 
  compare_dataset %>% 
  ggplot(aes(x = N_report, y = N_micro)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_text_repel(data = filter(compare_dataset, (abs(N_report - N_micro) / N_report) > 0.05),
                  aes(label = nid)) + 
  labs(title = "Comparing sample size of report data to extracted microdata", 
       x = "Report data", y = "Extracted microdata") + 
  coord_equal() + 
  theme_bw()

# plot results 
pdf("<<<< FILEPATH REDACTED >>>>")
gg_mc
gg_N
dev.off()
