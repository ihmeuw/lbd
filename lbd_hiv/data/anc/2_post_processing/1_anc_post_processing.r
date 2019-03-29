## ########################################################################################################################################
##
## ANC DATA CLEANING & GEOGRAPHY MATCHING
## 
## Purpose: Prepare ANC data from UNAIDS files and reports for collapse.
##            Bind all ANC data into one dataset, 
##            merge with point and polygon information in the ANC geography codebook,
##            and format and subset to prepare data for collpase.
##          
## ########################################################################################################################################

## Set-up---------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

libs <- c("data.table", "openxlsx", "raster", "foreign", "rgdal", "geosphere", "fossil", "dplyr", "rgeos", "car","plyr")
sapply(libs, require, character.only = T)
  
## Bind all UNAIDS Survaillance files together---------------------------------------------------------------------------------------------

df <- data.frame()

# Bind 2017 data to dataframe
unaid_year <- 2017
countries <- c("AGO","BDI", "CMR", "DJI", "COD","ETH", "LSO","MLI","MOZ","NGA", "RWA", "SOM", "UGA","SSD", "GAB", "TZA",
                "GHA", "ERI", "LBR", "CAF", "GMB","GNB","SWZ","MDG","GIN","NER", "ZAF","SDN","BWA","NAM","TGO", 
                "MAR","TCD","MWI","BFA","COG","BEN","CIV","ZMB")

for(i in 1:length(countries)){
  message(paste("processing", countries[i], i, " of ", length(countries)))
  surveillance_files <- list.files("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T)
  surveillance_files <- grep("SURVEILLANCE", surveillance_files, ignore.case = T, value = T)
  for(j in 1:length(surveillance_files)){
    message(paste("processing file             ", j, " of ", length(surveillance_files)))
    csv <- read.csv(surveillance_files[j])
    csv$Country <- rep(countries[i],dim(csv)[1])
    csv$data_source <- as.character(rep(unaid_year,dim(csv)[1]))
    csv$nid <- rep(317607,dim(csv)[1])
    colnames(csv)[colnames(csv)== "?..Group"] <- "Group"
    df <- rbind.fill(df, csv)
  }
}

# Bind 2016 data to dataframe  
unaid_year <- 2016
countries <- c("CPV","SLE","MRT","ZWE")
surveillance_files <- list.files("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T)
surveillance_files <- grep("SURVEILLANCE", surveillance_files, ignore.case = T, value = T)

for(i in 1:length(countries)){
  message(paste("processing", countries[i], i, " of ", length(countries)))
  temp <- grep(countries[i], surveillance_files, ignore.case = T, value = T)
  for(j in 1:length(temp)){
    message(paste("processing file             ", j, " of ", length(temp)))
    csv <- read.csv(temp[j])
    csv$Country <- rep(countries[i],dim(csv)[1])
    csv$data_source <- as.character(rep(unaid_year,dim(csv)[1]))
    csv$nid <- rep(306517,dim(csv)[1])
    colnames(csv)[colnames(csv)== "?..Group"] <- "Group"
    df <- rbind.fill(df, csv)
  }
}

# Remove duplicate data - use the data that is "In" the file if discrepency for a site
df_in <- subset(df, In == 1)
df_in <- unique(df_in)
df_not_in <- subset(df, In == 0)
df_not_in <- unique(df_not_in)
together <- rbind.fill(df_in, df_not_in)
together_dedupe <- unique(setDT(together)[order(Group,Site,Year, -In)], by = c('Group','Site','Year'))
df <- data.frame(together_dedupe)
df <- df[,!(names(df) %in% c("In"))]
rm(csv,unaid_year,temp,surveillance_files,i,j,countries, df_in, df_not_in, together, together_dedupe)
  
#Bring In Additional Data and bind to dataframe
zeros <- read.csv("<<<< FILEPATH REDACTED >>>>")
additional <- read.csv("<<<< FILEPATH REDACTED >>>>")
extraction_2018 <- read.csv("<<<< FILEPATH REDACTED >>>>")
df <- rbind.fill(df, zeros)
df <- rbind.fill(df, extraction_2018)

#Incorporate Additional Data and Discrepencies
new_data <- subset(additional, Country == "SEN" | type_of_additional %in% c(1,2))
new_data <- new_data[,c("Group","Site", "Year", "Prev", "N","Country","data_source","nid")]
df <- rbind.fill(df, new_data)
setDT(df)

replace_N <- subset(additional, type_of_additional %in% c(3,34))
setDT(replace_N)
df <- df[replace_N, on = .(Group,Site,Year), ':=' (N = i.N)]

replace_Prev <- subset(additional, type_of_additional %in% c(4,34))
setDT(replace_Prev)
df <- df[replace_Prev, on = .(Group,Site,Year), ':=' (Prev = i.Prev, data_source = "report", nid = as.numeric(i.nid))]

# Remove duplicates from TGO
df <- subset(df, !(Country == "TGO" & Group %in% c("POLULATION TOTALE","POPULATION TOTALE")))

#Remove duplicate rows
df <- unique(df)
rm(new_data, replace_N,replace_Prev)

## Merge Geography codebook and UNAIDS data -----------------------------------------------------------------------------------------------
geocodebook <- read.csv('<<<< FILEPATH REDACTED >>>>')
source('<<<< FILEPATH REDACTED >>>>/lbd_hiv/anc/2_post_processing/anc_snap.R')
geocodebook <- snap_codebook(geocodebook)
all <- merge(geocodebook, df, by.x=c("iso3", "group", "site"), by.y=c("Country", "Group", "Site"), all.x=F, all.y=T)

## Subsetting -----------------------------------------------------------------------------------------------------------------------------
# Rename blank site name in ZAF
all$site <- as.character(all$site)
all$site[all$iso3 == "ZAF" & all$site ==""] <- "KWAZULU-NATAL PROVINCE"

# Remove pseudo sites from MOZ
all <- subset(all, !(iso3 == "MOZ" & site %in% c("pseudo site", "Pseudo site", "Pseudo sites")))

# Remove COG rows where an admin 1 masks an admin 2 and remove all hospitals
all <- subset(all, !(iso3 == "COG" & Year == 2011 & point == 1)) 
all <- subset(all, !(iso3 == "COG"& Year == 2006 & site %in% c("Sibiti (%)","Nkayi (%)","Ouesso (%)"))) 
all <- subset(all, !(iso3 == "COG"& Year == 1992 & site =="Sibiti (%)"))

# Remove GNB hospitals only using admins
all <- subset(all, !(iso3 == "GNB" & point == 1))

# Remove non-ANC
all <- subset(all, !(iso3 == "CPV" & !(group %in% c("Pop f,minine restante","Pop féminine restante"))))
all <- subset(all, !(iso3 == "MAR" & !(group %in% c("Femmes"))))
all <- subset(all, !(iso3 == "MDG" & !(group %in% c("population feminine restante"))))
all <- subset(all, !(iso3 == "MRT" & !(group %in% c("Pop féminine restante"))))
all <- subset(all, !(iso3 == "NER" & !(group %in% c("Pop féminine restante"))))
all <- subset(all, !(iso3 == "SDN" & !(group %in% c("Female remaining pop"))))
all <- subset(all, !(iso3 == "STP" & !(group %in% c("Pop Fem_restante"))))
all <- subset(all, !(iso3 == "SEN" & !(group %in% c("Pop féminine restante"))))

# Remove unnecessary characters and white space in site variable
all$site <- gsub('"',"",as.character(all$site))
all$site <- gsub("\\(%\\)","", as.character(all$site))
all$site <- gsub("^\\s+|\\s+$", "", all$site)

# Remove duplicate data that represents same site but has different site name
# Also change different site names that represent same site to same name
all <- subset(all, !(iso3 == "RWA" & site == "Gikondo A"& Year != 2013))
all$site[all$iso3 == "RWA"& all$site == "Gikondo A"] <- "Gikondo"
all <- subset(all, !(iso3 == "RWA" & site == "Biryogo A"& Year != 2013))
all$site[all$iso3 == "RWA"& all$site == "Biryogo A"] <- "Biryogo"
all$site[all$iso3 == "RWA"& all$site == "CHK"] <- "Kigali CHK"
all$site[all$iso3 == "RWA"& all$site == "Ruhengeri 2"] <- "Ruhengeri"
all$site[all$iso3 == "TZA"& all$site == "Kagera (Lukole Refugee Camp)"] <- "Lukole Refugee"

# Remove outliers and RT
all <- subset(all, !(iso3 == "CIV" & site == "MATERNITE HG ABOBO SUD")) #outlier
all <- subset(all, !(iso3 == "CMR" & site == "Lolodorf" & Year == 2012)) #outlier
all <- subset(all, !(iso3 == "COD" & site == "S.Vicente" & Year == 1996)) #outlier
all <- subset(all, !(iso3 == "AGO" & site %in% c("CMI Benguela - Benguela","Mat do Lobito - Benguela") & Year == 2004)) #outlier
all <- subset(all, !(iso3 == "MOZ" & Year >= 2013)) #RT vs SS (9 site-years)
all <- subset(all, !(iso3 == "MWI" & Year == 2016)) #RT vs SS (3 site-years)

# For countries with many repeated sample sizes at a set level, 
#    change it to be the median of the others
all$N[all$iso3 == "BDI" & all$N == 300] <- 349
all$N[all$iso3 == "CMR" & all$N == 300] <- 100
all$N[all$iso3 == "COG" & all$N == 300] <- 121
all$N[all$iso3 == "LSO" & all$N == 300] <- 325
all$N[all$iso3 == "MLI" & all$N == 300] <- 299
all$N[all$iso3 == "TGO" & all$N == 300] <- 103
all$N[all$iso3 == "UGA" & all$N == 300] <- 470
all$N[all$iso3 == "ZWE" & all$N == 300] <- 333
all$N[all$iso3 == "AGO" & all$N == 500] <- 498

# Only keep variables we need
all <- all[,c("nid", "iso3", "group", "site", "point", "latitude","longitude",
              "admin_level","shapefile","location_code","Year","Prev", "N","data_source")]

## Separate geomatched and non-geomatched data --------------------------------------------------------------------------------------------
all_mapped <- subset(all,!is.na(point))

## Rename variables and change data type to prep for collapse -------------------------------------------------------------------------------
all_mapped$latitude[all_mapped$point == 0] <- NA
all_mapped$longitude[all_mapped$point == 0] <- NA
all_mapped <- transform(all_mapped, latitude = as.numeric(as.character(latitude)),
                                    longitude = as.numeric(as.character(longitude)))
all_mapped$shapefile <- as.character(all_mapped$shapefile)
all_mapped$location_code <- as.numeric(all_mapped$location_code)
all_mapped$source <- paste0("UNAIDS files - ", all_mapped$data_source)
colnames(all_mapped)[colnames(all_mapped) == "Year"] <- "year"
colnames(all_mapped)[colnames(all_mapped) == "Prev"] <- "anc_hiv_test"
colnames(all_mapped)[colnames(all_mapped) == "iso3"] <- "country"
all_mapped$N_obs <- all_mapped$N
all_mapped$anc_hiv_test <- all_mapped$anc_hiv_test *.01
all_mapped <- unique(all_mapped)
all_mapped$cluster_id <- 1:dim(all_mapped)[1]

all_mapped <- data.frame(all_mapped)

## Save files -------------------------------------------------------------------------------------------------------------------------------
module_date <- Sys.Date()
module_date <- gsub("-", "_", module_date)
saveRDS(all_mapped, file=paste0("<<<< FILEPATH REDACTED >>>>/anc_data_",module_date,".rds"))




