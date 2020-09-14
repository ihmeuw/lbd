## ######################################################################################################################
##
## StatCompiler Comparison Scatterplots
##
## Purpose: Creates a pdf of scatterplots comparing prevalence and sample size
##          in extracted microdata vs. StatCompiler and/or statcompilers
##
##  Prior to running this code, dowload the StatCompiler table for your indicator
##  (see "StatCompiler Comparison Checks" page on the HUB)
##
## #######################################################################################################################

rm(list=ls())

# Load required packages
libs <- c("plyr", "dplyr", "gdata", "openxlsx", "tidyr")
sapply(libs, require, character.only = T)
j <- "<<<< FILEPATH REDACTED >>>>"

microdata_filepath     <- "..."  # Fill in with filepath to latest microdata extraction (including file name)
statcompiler_filepath  <- "..."  # Fill in with filepath to downloaded Statcompiler database for indicator (including file name)
pdf_filepath           <- "..."  # Fill in with filepath to put final scatterplot pdf (including file name)
indicator_microdata    <- "..."  # Fill in name of indicator from microdata
indicator_statcompiler <- "..."  # Fill in variable name of indicator from StatCompiler
weighted               <- F      # Change to T is using weighted denominators
seperate_sex           <- F      # Change to T if indicator seperates men and women
indicator_sex          <- "..."  # 'WN', 'MN', or 'BOTH' (leave as "" if seperate_sex = T)
lower_age_range        <- 15     # Fill in with number if indicator is subset by age in Statcompiler, NA if there is no age range
upper_age_range        <- 49     # Fill in with number if indicator is subset by age in Statcompiler, NA if there is no age range


## Load extracted microdata-----------------------------------------------------------------------------------------------
# Pull in extracted data
all <- readRDS(microdata_filepath)

# If applicable, add indicator we are interested in to the microdata
all$indicator <- ... # Fill in. See examples below.

# Example:
## Indicator already in microdata
# all$indicator <- all$hiv_test
## Create indicator: in_union - true if currently married/living with partner, false otherwise (missing if marital_status is missing)
# all$indicator <- ifelse(all$marital_status < 3, 1, 0)


## Format StatCompiler data to compare with aggregated microdata---------------------------------------------------------
# Pull in table of statcompiler data
statcompiler <- read.xlsx(statcompiler_filepath) %>%
                  subset(!is.na(Value) & Characteristic.Category == "Total 15-49") # Change Characteristic.Category to reflect your indicator

statcompiler <- statcompiler %>%
                  subset(Survey.Year >= 1998) %>%
                    dplyr::select(Country.Name, Survey.Year, Survey.Name, Indicator, Value) # Change year if neccessary

# Keep only weighted (or unweighted) sample size
if (weighted == F){
  statcompiler <- statcompiler[grepl("Number*", statcompiler$Indicator) & grepl("*unweighted*", statcompiler$Indicator) | !grepl("Number*", statcompiler$Indicator), ]
} else{
  statcompiler <- statcompiler[grepl("Number*", statcompiler$Indicator) & !grepl("*unweighted*", statcompiler$Indicator) | !grepl("Number*", statcompiler$Indicator), ]
}

# Indicate if row represents men, women, or both
if (seperate_sex == T){
  statcompiler$sex <- ifelse(grepl("*Women*", statcompiler$Indicator, ignore.case=T), "WN", "MN")
} else {
  statcompiler$sex <- indicator_sex
}

# Drop rows that are the % of 'missing' responses (sometimes these are NA, sometimes they are 0)
statcompiler <- statcompiler[!grepl("*Missing*", statcompiler$Indicator, ignore.case=T), ]

# Drop rows that have 'total' in the indicator column (assuming that value is 100 for all of these. check to be sue this is true)
statcompiler <- statcompiler[!grepl("*Total*", statcompiler$Indicator, ignore.case=T), ]

# Rename indicator
statcompiler$Indicator <- ifelse(grepl(paste("*",indicator_statcompiler, sep = ""), statcompiler$Indicator, ignore.case=T),
                                 "indicator",
                                 "sample_size")

# Reshape long to wide
statcompiler <- reshape(statcompiler,
                        idvar=c("Country.Name", "Survey.Year", "Survey.Name", "sex"),
                        v.names="Value", timevar="Indicator", direction="wide")

# replace Survey Name with survey series to match extracted data
statcompiler$survey_name <- ifelse(grepl("*DHS*", statcompiler$Survey.Name), "MACRO_DHS", "MACRO_AIS")

# rename columns
statcompiler <- statcompiler %>%
                dplyr::select(Country.Name,
                              year=Survey.Year,
                              survey_name,
                              sex,
                              statcompiler_indicator = Value.indicator,
                              statcompiler_sample_size=Value.sample_size)

statcompiler$start_age <- lower_age_range
statcompiler$end_age <- upper_age_range

# read in map of country names to IHME country codes and bind codes to sc_db
ctry_codes <- read.csv("<<<< FILEPATH REDACTED >>>>") #Needs to be updated for Stage 2 and 3 countries
ctry_codes <- ctry_codes[,c("iso3","location_name")]
ctry_codes$location_name <- as.character(ctry_codes$location_name)
ctry_codes$location_name[ctry_codes$location_name == "Republic of the Congo"] <- "Congo"
ctry_codes$location_name[ctry_codes$location_name == "Democratic Republic of the Congo"] <- "Congo Democratic Republic"
ctry_codes$location_name[ctry_codes$location_name == "The Gambia"] <- "Gambia"
colnames(ctry_codes) <- c("country", "Country.Name")
statcompiler <- merge(statcompiler, ctry_codes, by = "Country.Name", all.x = TRUE)


statcompiler$nid <- NA


## Aggregate the microdata and merge to StatCompiler data for comparison ------------------------------------

# Finding surveys by survey_series, country, year, & sex
# Calculate the statcompiler size and prevalence rates for each var of interest for each extracted survey
for (i in 1:nrow(statcompiler)){
  # Identify the survey matching current row in statcompiler
  if (seperate_sex == T | indicator_sex != "BOTH"){
    statcompiler_sex_id = ifelse(statcompiler$sex[i] == "MN", 1, 2) # AIS have HHM mods but statcomp splits into mn & wn
    survey <- subset(all, survey_name == statcompiler$survey_name[i] & country == statcompiler$country[i] &
                       (year == statcompiler$year[i] | end_year == statcompiler$year[i]) & sex_id == statcompiler_sex_id)
  } else {
    survey <- subset(all, survey_name == statcompiler$survey_name[i] & country == statcompiler$country[i] &
                       (year == statcompiler$year[i] | end_year == statcompiler$year[i]))
  }

  # Subset to only include ages included in statcompiler
  if (!is.na(statcompiler$start_age[i]))
    survey <- subset(survey, age_year >= statcompiler$start_age[i])

  if (!is.na(statcompiler$end_age[i]))
    survey <- subset(survey, age_year <= statcompiler$end_age[i])

  # Add on the nid from the microdata
  if (is.na(statcompiler$nid[i]))
    statcompiler$nid[i] = survey$nid[1]

  statcompiler$microdata_sample_size[i] = nrow(survey)
  statcompiler$microdata_indicator[i] = 100 *sum(survey$pweight * survey$indicator, na.rm=T)/sum(survey$pweight, na.rm=T)
}

# Only look at the surveys with extracted microdata
statcompiler <- subset(statcompiler, !is.na(nid))

# Set bounds for plot
boundary_prev = 2
boundary_size_percent = 5
boundary_size = 1 + (boundary_size_percent) / 100

#Saves plots in pdf
pdf(pdf_filepath)

# subset to just the surveys that actually asked the variable we're plotting
#Create prevalence plot
plot(statcompiler$statcompiler_indicator, statcompiler$microdata_indicator,
     main=paste(indicator_microdata, "Prevalences"),
     xlab = "StatCompiler", ylab= "Microdata")
rep <- subset(statcompiler, abs(statcompiler_indicator - microdata_indicator) > boundary_prev)
text(rep$statcompiler_indicator, rep$microdata_indicator,
     paste(rep$nid,rep$sex, sep="_"), cex=0.6, pos=3, col="red") # label surveys outside boundary
abline(0, 1, col = "blue")
abline(boundary_prev, 1)
abline(-boundary_prev, 1)

# Create statcompiler size plot
plot(statcompiler$statcompiler_sample_size, statcompiler$microdata_sample_size,
     main=paste(indicator_microdata, "Sample Sizes"),
     xlab = "StatCompiler", ylab = "Microdata")
rep <- subset(statcompiler, abs(statcompiler_sample_size - microdata_sample_size)/statcompiler_sample_size > boundary_size - 1)
text(rep$statcompiler_sample_size, rep$microdata_sample_size,
     paste(rep$nid, rep$sex, sep="_"), cex=0.6, pos=4, col="red") # label surveys outside boundary
abline(0, 1, col = "blue")
abline(0, boundary_size)
abline(0, 1 / boundary_size)

dev.off()
