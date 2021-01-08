# This is launch script for SAE models, where it creates a model directory under sae/stage2/country/outputs


library(Matrix)
library(data.table)
library(dplyr)
library(readr)

rm(list=ls())

## Print message for clarity
message(commandArgs()[4])
message(commandArgs()[5])
message(commandArgs()[6])
message(commandArgs()[7])

## Set indicator
country       <- commandArgs()[4]
model_number  <- commandArgs()[5]
run_date      <- commandArgs()[6]
resub         <- commandArgs()[7]

# Define model directory 
model_dir <- paste0("<<<< FILEPATH REDACTED >>>>")

## Get settings ------------------------------------------------------------------------------------
# Make sure everything is in order
setwd("<<<< FILEPATH REDACTED >>>>/lbd_core/sae_central/")

most_recent_date <- function(dir, date_format = "%Y_%m_%d", out_format = "%Y_%m_%d", file_pattern = NULL) {
  date_pattern <- gsub("%y|%m|%d", "[[:digit:]]{2}", gsub("%Y", "[[:digit:]]{4}", date_format))
  dates <- dir(dir, pattern = date_pattern)
  if (!is.null(file_pattern)) dates <- grep(file_pattern, dates, value = T)
  dates <- gsub(paste0("(.*)(", date_pattern, ")(.*)"), "\\2", dates)
  dates <- as.Date(dates, date_format)
  format(max(dates), out_format)
}

# Grab settings folder 
data_run_date <- most_recent_date(paste0("<<<< FILEPATH REDACTED >>>>"))
settings_folder <- paste0("<<<< FILEPATH REDACTED >>>>")

# Check settings 
source("<<<< FILEPATH REDACTED >>>>")
get_settings(settings_folder)
check_settings(settings_folder)

# Temp directory is assigned to model output directory, save settings there 
settings <- read.csv(paste0(settings_folder, "settings.csv"), header = F, stringsAsFactors = F) %>% data.table()
settings[V1 == "temp_dir", V2 := paste0(model_dir, "/temp_dir/")]

# Change model number, save settings to model directory 
settings[V1 == "model", V2 := as.character(model_number)]
write_csv(settings, path = paste0(model_dir, "/settings.csv"))

# Specify model arguments 
type <- "all"
proj <- "proj_geo_nodes"
core_dir <- "<<<< FILEPATH REDACTED >>>>"

error_dir <- paste0(model_dir, "errors/")
dir.create(error_dir)
output_dir <- paste0(model_dir, "output/")
dir.create(output_dir)
shell <- paste0('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/share_scripts/shell_sing.sh')

code <- "<<<< FILEPATH REDACTED >>>>/lbd_core/sae_central/submit.r"


qsub <- paste('qsub -e', error_dir, '-o', output_dir, '-l m_mem_free=90G -l fthread=5', 
              paste0("-N sae_", country, "_model_", as.numeric(settings[V1 == "model", V2])),
              '-P proj_geo_nodes -q geospatial.q -l h_rt=08:00:00:00 -l archive=TRUE',
              '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
              shell, code,
              model_dir, type, resub, proj)



system(qsub)