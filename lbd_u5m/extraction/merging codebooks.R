## ##############################################################
## Inputs: geography codebooks
## Outputs: Stata dataset of rbinded and formatted geographies
## Description: This script combines geography codebooks 
##    for use in geomatching CBH and BH extractions in stata
## Dependencies: haven
## ##############################################################

################################################################
######################### Setup ################################
################################################################
rm(list = ls())

library(haven)

setwd("<<<< FILEPATH REDACTED >>>>")
outpath <- "<<<< FILEPATH REDACTED >>>>"

################################################################
#################### Combine Geography Codebooks ###############
################################################################
# List csv's
fileslist <- list.files(, pattern = '*csv$')

for (i in 1:length(fileslist)) {
  file_in <- fileslist[i]
  
  mydata <- read.csv(file_in, stringsAsFactors = FALSE)
  
  #rename lat long columns as stata doesnt like long as a var name
  colnames(mydata)[which(colnames(mydata)=="lat")] <- "latnum"
  colnames(mydata)[which(colnames(mydata)=="long")] <- "longnum"
  colnames(mydata)[which(colnames(mydata)=="geospatial_id")] <- "cluster_number"
  mydata$cluster_number <- as.character(mydata$cluster_number)
  
  ##select required variables & remove duplicates:
  
  mydata <- unique(mydata[c("nid", "iso3", "cluster_number", "point", "latnum", "longnum",  
                            "old_uncertainty", "old_buffer","location_name", "location_code",
                            "admin_level", "shapefile" )])
  
  if (file_in == fileslist[1]){
    combineddata <- mydata
  }  else{
    combineddata <- rbind(combineddata, mydata)
  }
  
}

#####################################################################
########################## Cleaning #################################
#####################################################################

##Remove the NA's and change numeric variables to numeric
combineddata$latnum[combineddata$latnum == "NA"] <- ""
combineddata$longnum[combineddata$longnum == "NA"] <- ""
combineddata$nid[combineddata$nid == "NA"] <- ""
combineddata$point[combineddata$point == "NA"] <- ""
combineddata$latnum <- as.numeric(combineddata$latnum)
combineddata$longnum <- as.numeric(combineddata$longnum)
combineddata$nid <- as.numeric(combineddata$nid)
combineddata$point <- as.numeric(combineddata$point)

##Change string variables to factors to allow write.dta to work
combineddata$iso3 <- as.factor(combineddata$iso3)
combineddata$cluster_number <- as.factor(combineddata$cluster_number)
combineddata$uncertain_point <- as.factor(combineddata$old_uncertainty)
combineddata$buffer <- as.factor(combineddata$old_buffer)
combineddata$location_name <- as.factor(combineddata$location_name)
combineddata$location_code <- as.factor(combineddata$location_code)
combineddata$admin_level <- as.factor(combineddata$admin_level)
combineddata$shapefile <- as.factor(combineddata$shapefile)

#Remove any duplicates in the file
combineddata <- unique(combineddata)

write.dta(combineddata, paste0(outpath, "combined_codebook.dta"))


