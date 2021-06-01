########################################################################################################################################
##
## HIV_TEST POST UBCOV EXTRACTION DATA CLEANING & GEOGRAPHY MATCHING
##
## Purpose: Prepare ANC data from UNAIDS files and reports for collapse.
##            Bind all ANC data into one dataset,
##            merge with point and polygon information in the ANC geography codebook,
##            and format and subset to prepare data for collpase.
##
## ########################################################################################################################################

## Set-up---------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())

libs <- c("data.table", "openxlsx", "raster", "foreign", "rgdal", "geosphere", "fossil", "dplyr", "rgeos", "car","plyr", "ggplot2", "gridExtra")
sapply(libs, require, character.only = T)

## Bind all UNAIDS Survaillance files together---------------------------------------------------------------------------------------------

df <- data.frame()


# Bind 2019 dataframe
unaid_year <- 2019
countries <- c("AGO","BDI", "BEN", "BFA", "BWA", "CAF", "CIV", "CMR", "COD", "COG", "COM", "CPV","DOM","HTI", "DJI", "ERI", "GAB", "GHA", "GMB", "GNB", "GNQ", "GIN","LSO", "LBR", "MAR", "MOZ",
               "MLI", "MRT", "MWI","MDG", "NAM", "NER", "NGA", "RWA", "SDN", "SEN", "SLE","SOM", "SSD", "TCD","SOM", "UGA", "TZA","TGO", "ZMB", "ZWE")


length(countries)
for(i in 1:length(countries)){
  message(paste("processing", countries[i], i, " of ", length(countries)))
  surveillance_files <- list.files(paste0("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T),
                                   pattern=".csv$", ignore.case = T, full.names = T)

  for(j in 1:length(surveillance_files)){
    csv <- read.csv(surveillance_files[j], stringsAsFactors = FALSE)
    csv$Country <- rep(countries[i],dim(csv)[1])
    csv$data_source <- as.character(rep(unaid_year,dim(csv)[1]))
    csv$nid <- rep(415978,dim(csv)[1])
    colnames(csv)[1] <- "Group"
    csv$UNAIDS_file <- NULL
    df <- rbind.fill(df, csv)
  }
}

# Bind 2018 dataframe
unaid_year <- 2018
countries <- c("AGO","BDI", "BEN", "BFA", "BWA", "CAF", "CIV", "CMR", "COD", "COG", "COM", "CPV", "DJI", "DOM","ERI","ETH", "GAB", "GHA", "GMB", "GNB", "GNQ", "GIN", "KEN", "LSO", "LBR", "MAR", "MOZ",
               "MLI","MWI","MRT", "MDG", "NAM", "NER", "NGA", "RWA","SEN", "SDN", "SLE","SOM", "SSD","SWZ", "TCD","SOM", "UGA", "TGO","TZA", "ZMB", "ZWE")

length(countries)
for(i in 1:length(countries)){
  message(paste("processing", countries[i], i, " of ", length(countries)))
  surveillance_files <- list.files(paste0("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T),
                                   pattern=".csv$", ignore.case = T, full.names = T)

  for(j in 1:length(surveillance_files)){
    csv <- read.csv(surveillance_files[j], stringsAsFactors = FALSE)
    csv$Country <- rep(countries[i],dim(csv)[1])
    csv$data_source <- as.character(rep(unaid_year,dim(csv)[1]))
    csv$nid <- rep(365412,dim(csv)[1])
    colnames(csv)[1] <- "Group"
    csv$UNAIDS_file <- NULL
    df <- rbind.fill(df, csv)
  }
}

# Bind 2017 data to dataframe
unaid_year <- 2017
countries <- c("AGO","BDI", "BEN", "BFA", "BWA", "CAF", "CIV", "CMR", "COD", "COG","DJI", "ERI", "GAB","DOM", "GHA", "GMB", "GNB", "GNQ", "GIN", "LSO", "LBR", "MAR", "MOZ",
               "MLI","MWI","MDG", "NAM", "NER", "NGA", "RWA", "SEN", "SDN","SOM", "SSD","SWZ", "TCD", "TGO", "SOM", "UGA", "TZA", "ZAF", "ZMB")

for(i in 1:length(countries)){
  message(paste("processing", countries[i], i, " of ", length(countries)))
  surveillance_files <- list.files(paste0("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T),
                                   pattern=".csv$", ignore.case = T, full.names = T)
  for(j in 1:length(surveillance_files)){
    message(paste("processing file             ", j, " of ", length(surveillance_files)))
    csv <- read.csv(surveillance_files[j], stringsAsFactors = FALSE)
    csv$Country <- rep(countries[i],dim(csv)[1])
    csv$data_source <- as.character(rep(unaid_year,dim(csv)[1]))
    csv$nid <- rep(317607,dim(csv)[1])
    csv$UNAIDS_file <- NULL
    colnames(csv)[1] <- "Group"
    #colnames(csv)[colnames(csv)== "?..Group"] <- "Group"
    df <- rbind.fill(df, csv)
  }
}

unaid_year <- 2016
countries <- c("AGO","BDI", "BFA", "BWA", "CAF", "CIV", "CMR", "COD","CPV","DJI", "ERI", "GAB", "HTI", "GHA","DOM","HTI", "GMB", "GNQ", "GIN", "LSO", "LBR", "MAR", "MOZ",
               "MLI", "MRT","MWI","MDG", "NAM", "NER", "RWA","SEN", "SDN", "SLE", "SSD","SWZ", "TCD", "TGO", "SOM", "UGA", "TZA", "ZMB", "ZAF", "ZWE")

for(i in 1:length(countries)){
  message(paste("processing", countries[i], i, " of ", length(countries)))
  surveillance_files <- list.files(paste0("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T),
                                   pattern=".csv$", ignore.case = T, full.names = T)
  for(j in 1:length(surveillance_files)){
    message(paste("processing file             ", j, " of ", length(surveillance_files)))
    csv <- read.csv(surveillance_files[j], stringsAsFactors = FALSE)
    csv$Country <- rep(countries[i],dim(csv)[1])
    csv$data_source <- as.character(rep(unaid_year,dim(csv)[1]))
    csv$nid <- rep(306517,dim(csv)[1])
    csv$UNAIDS_file <- NULL
    colnames(csv)[1] <- "Group"
    df <- rbind.fill(df, csv)
  }
}
### Cleaning site names ---------------------------------------------------------------------------------------------------------------------

df$Group[df$Country == "HTI" & df$Group == "Urban"] <- ""
df$Group[df$Country == "HTI" & df$Group == "Rural"] <- ""
df$Group[df$Country == "HTI" & df$Group == "Department"] <- ""

# LBR site differences in CSV raw data

df$Site[df$Country == "LBR" & df$Site == "Barclay.Health.Center"] <- "Barclay Health Center"
df$Site[df$Country == "LBR" & df$Site == "C.B..Dunbar.Health.Center...."] <- "C.B. Dunbar Health Center (%)"
df$Site[df$Country == "LBR" & df$Site == "C.H..Rennie.Hospital...."] <- "C.H. Rennie Hospital (%)"
df$Site[df$Country == "LBR" & df$Site == "Duport.Road.Health.Center...."] <- "Duport Road Health Center (%)"
df$Site[df$Country == "LBR" & df$Site == "Firestone.Hospital"] <- "Firestone Hospital"
df$Site[df$Country == "LBR" & df$Site == "Fish.Town.Health.Center"] <- "Fish Town Health Center"
df$Site[df$Country == "LBR" & df$Site == "Ganta.Hospital"] <- "Ganta Hospital"
df$Site[df$Country == "LBR" & df$Site == "JJ.Dossen.Hospital"] <- "JJ Dossen Hospital"
df$Site[df$Country == "LBR" & df$Site == "John.F.Kennedy.Hospital"] <- "John F Kennedy Hospital"
df$Site[df$Country == "LBR" & df$Site == "Liberia.Govt..Hospital.Tubmanburg"] <- "Liberia Govt. Hospital-Tubmanburg"
df$Site[df$Country == "LBR" & df$Site == "Liberia.Govt.Hosp.Buchanan"] <- "Liberia Govt Hosp-Buchanan"
df$Site[df$Country == "LBR" & df$Site == "Martha.Tubman.Hospital"] <- "Martha Tubman Hospital"
df$Site[df$Country == "LBR" & df$Site == "Phebe.Hospital"] <- "Phebe Hospital"
df$Site[df$Country == "LBR" & df$Site == "Redemption.Hospital"] <- "Redemption Hospital"
df$Site[df$Country == "LBR" & df$Site == "Site.22"] <- "Site 22"
df$Site[df$Country == "LBR" & df$Site == "St.Francis.Hospital"] <- "St Francis Hospital"
df$Site[df$Country == "LBR" & df$Site == "Star.of.the.Sea.Health.Center"] <- "Star of the Sea Health Center"
df$Site[df$Country == "LBR" & df$Site == "Voinjama.Health.Center"] <- "Voinjama Health Center"

# More site differences from raw data- special characters pattern in most recent collection years
df$Site[df$Country == "MDG" & df$Site == "Sambava  (%)"] <- "Sambava (%)"
df$Site[df$Country == "MDG" & df$Site == "Antsirabe  (%)"] <- "Antsirabe (%)"
df$Site[df$Country == "NGA" & df$Site == "Osun  (Iragbere)"] <- "Osun (Iragbere)"
df$Site[df$Country == "ZAF" & df$Site == "NULL"] <- "KWAZULU-NATAL PROVINCE"
df$Site[df$Country == "SOM" & df$Site == "NULL"] <- "Bosaso Central MCH"
df$Site[df$Country == "COG" & df$Site == "Cuvette.ouest"] <- "Cuvette ouest"
df$Site[df$Country == "COG" & df$Site == "Mouyondzi...."] <- "Mouyondzi (%)"
df$Site[df$Country == "COG" & df$Site == "Owando...."] <- "Owando (%)"
df$Site[df$Country == "COG" & df$Site == "Sibiti...."] <- "Sibiti (%)"
df$Site[df$Country == "CMR" & df$Site == "CMA  Tyo"] <- "CMA Tyo"
df$Site[df$Country == "BFA" & df$Site == "Dedougou...."] <- "Dedougou (%)"
df$Site[df$Country == "BFA" & df$Site == "Dori...."] <- "Dori (%)"
df$Site[df$Country == "BFA" & df$Site == "Kaya...."] <- "Kaya (%)"
df$Site[df$Country == "BFA" & df$Site == "Manga...."] <- "Manga (%)"
df$Site[df$Country == "BFA" & df$Site == "Sindou...."] <- "Sindou (%)"
df$Site[df$Country == "BDI" & df$Site == "CDS Ijenda"] <- "Ijenda"
df$Site[df$Country == "BDI" & df$Site == "CDS Kiremba"] <- "Kiremba"
df$Site[df$Country == "BDI" & df$Site == "CDS RUMONGE"] <- "H RUMONGE"
df$Site[df$Country == "BDI" & df$Site == "CDS CMC BUYENZI"] <- "CMC BUYENZI"
df$Site[df$Country == "TZA" & df$Site == "Bukoba Regional Hospital"] <- "Kagera (Bukoba Town)"

# Remove duplicate data - use the data that is "In" the file if discrepency for a site for regional sites
df_in <- subset(df, In == 1)
df_in <- unique(df_in)
df_not_in <- subset(df, In == 0)
df_not_in <- unique(df_not_in)
together <- rbind.fill(df_in, df_not_in)
together_dedupe <- unique(setDT(together)[order(Group,Site,Year, -In)], by = c('Group','Site','Year'))

# Subset regional data by "In" = 1,  assigned to a region not a population like other countries.
df_in_regions <- subset(together_dedupe, together_dedupe$Country %in% c("BEN", "CIV", "ETH", "KEN", "MOZ", "NGA", "TGO", "HTI", "ZMB", "ZWE"))
df_inout_regions <- subset(together_dedupe, !(together_dedupe$Country %in% c("BEN", "CIV", "ETH", "KEN", "MOZ", "NGA","HTI", "TGO", "ZMB", "ZWE")))
df_in_regions <- subset(df_in_regions, df_in_regions$In == 1)
df <- rbind.fill(df_in_regions, df_inout_regions)
df <- df[,!(names(df) %in% c("In"))]
rm(csv,unaid_year,surveillance_files,i,j,countries, df_in, df_not_in, together, together_dedupe, df_in_regions, df_inout_regions)


# Bring In Additional Data and bind to dataframe,
additional <- read.csv("<<<< FILEPATH REDACTED >>>>")
additional$Site <- trimws(additional$Site)


# Incorporate Additional Data and Discrepencies
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

# Remove duplicate rows
df <- unique(df)
rm(new_data, replace_N,replace_Prev)

## Merge Geography codebook and UNAIDS data -----------------------------------------------------------------------------------------------
geocodebook <- read.csv("<<<< FILEPATH REDACTED >>>>", stringsAsFactors = FALSE)
source('<<<< FILEPATH REDACTED >>>>/lbd_hiv/anc/2_post_processing/anc_snap.R')
geocodebook <- snap_codebook(geocodebook)


## Added trim/whitespace code to match more sites
geocodebook$site <- gsub(' $',"",as.character(geocodebook$site))
geocodebook$site <- trimws(geocodebook$site)

# Merge to geography information
all <- merge(geocodebook, df, by.x=c("iso3", "group", "site"), by.y=c("Country", "Group", "Site"), all.x=F, all.y=T)
all_copy <- all

## Subsetting -----------------------------------------------------------------------------------------------------------------------------
all$site <- as.character(all$site)

# Remove pseudo sites from MOZ
all <- subset(all, !(iso3 == "MOZ" & site %in% c("pseudo site", "Pseudo site", "Pseudo sites")))

# Remove COG rows where an admin 1 masks an admin 2 and remove all hospitals
all <- subset(all, !(iso3 == "COG" & Year == 2011 & point == 1)) # 23 site-years dropped
all <- subset(all, !(iso3 == "COG"& Year == 2006 & site %in% c("Sibiti (%)","Nkayi (%)","Ouesso (%)"))) # 3 site-years dropped
all <- subset(all, !(iso3 == "COG"& Year == 1992 & site =="Sibiti (%)")) # 1 site-year dropped (1 pre-2000)

# Remove GNB hospitals only using admins
all <- subset(all, !(iso3 == "GNB" & point == 1)) # 10 site-years dropped

# Remove non-ANC groups
all <- subset(all, !(iso3 == "CPV" & !(group %in% c("Pop f,minine restante","Pop fÃ©minine restante"))))
all <- subset(all, !(iso3 == "MAR" & !(group %in% c("Femmes"))))
all <- subset(all, !(iso3 == "MDG" & !(group %in% c("population feminine restante"))))
all <- subset(all, !(iso3 == "MRT" & !(group %in% c("Pop fÃ©minine restante"))))
all <- subset(all, !(iso3 == "NER" & !(group %in% c("Pop fÃ©minine restante"))))
all <- subset(all, !(iso3 == "SDN" & !(group %in% c("Female remaining pop"))))
all <- subset(all, !(iso3 == "SEN" & !(group %in% c("Pop fÃ©minine restante"))))
all <- subset(all, !(iso3 == "SOM" & !(group %in% c("Remaining females"))))
all <- subset(all, !(iso3 == "COM" & !(group %in% c("Female Population"))))
all <- subset(all, !(iso3 == "DOM" & !(group %in% c("Remanente Femenina"))))

# Removing any other weird special character things
all$site <- gsub('"',"",as.character(all$site))
all$site <- gsub("\\(%\\)","", as.character(all$site))
all$site <- gsub("\\.^\\s+|\\s+$", "", all$site)

# Remove duplicate data that represents same site but has different site name
# Also change different site names that represent same site
all <- subset(all, !(iso3 == "RWA" & site == "Gikondo A"& Year != 2013))
all$site[all$iso3 == "RWA"& all$site == "Gikondo A"] <- "Gikondo"
all <- subset(all, !(iso3 == "RWA" & site == "Biryogo A"& Year != 2013))
all$site[all$iso3 == "RWA"& all$site == "Biryogo A"] <- "Biryogo"
all$site[all$iso3 == "RWA"& all$site == "CHK"] <- "Kigali CHK"
all$site[all$iso3 == "RWA"& all$site == "Ruhengeri 2"] <- "Ruhengeri"
all$site[all$iso3 == "TZA"& all$site == "Kagera (Lukole Refugee Camp)"] <- "Lukole Refugee"


# Removing duplicates that were mapped and come from different sources
all <- subset(all, !(iso3 == "MWI" & site == "Blantyre Limbe HC" & data_source == "2017" ))
all <- subset(all, !(iso3 == "NGA" & site == "Site 2 - Lafia" & data_source == "2018" ))
all <- subset(all, !(iso3 == "NGA" & site == "Site 1 - N/ Eggon" & data_source == "2018" ))
all <- subset(all, !(iso3 == "RWA" & site == "Kirehe" & data_source == "2016" ))
all <- subset(all, !(iso3 == "RWA" & site == "Ruli" & data_source == "2016" ))
all <- subset(all, !(iso3 == "RWA" & site == "Ruli" & data_source == "2016" ))
all <- subset(all, !(iso3 == "RWA" & site == "Ruhengeri" & data_source == "2016" ))


# Remove site-years with concerning data (see "Checking weird Trends in ANC Data")
all <- subset(all, !(iso3 == "CIV" & site == "MATERNITE HG ABOBO SUD")) #outlier
all <- subset(all, !(iso3 == "CMR" & site == "Lolodorf" & Year == 2012)) #outlier
all <- subset(all, !(iso3 == "COD" & site == "S.Vicente" & Year == 1996)) #outlier
all <- subset(all, !(iso3 == "AGO" & site %in% c("CMI Benguela - Benguela","Mat do Lobito - Benguela") & Year == 2004)) #outlier
all <- subset(all, !(iso3 == "MOZ" & Year >= 2013)) #RT vs SS (9 site-years)
all <- subset(all, !(iso3 == "MWI" & Year == 2016)) #RT vs SS (3 site-years)

# Different sample sizes in 2019 vs. reports (going with sample size that is closer to trend line)
all <- subset(all, !(iso3 == "ETH" & site == "Belle HC" & Year == 2007 & N == 300))
all <- subset(all, !(iso3 == "ETH" & site == "Debate HC" & Year == 2007 & N == 308))
all <- subset(all, !(iso3 == "GHA" & site == "Nadowli" & Year == 2006 & N == 341))

# Removing Spectrum data where reports are closer to the trend line
all <- subset(all, !(iso3 == "KEN" & site == "Kisii" & Year == 1992 & data_source == "2018" ))
all <- subset(all, !(iso3 == "GHA" & site == "Nadowli" & Year == 2006 & data_source == "2018" ))

# Removing report data that was identical to  Spectrum 2019
all <- subset(all, !(iso3 == "ZAF" & site == "Namakwa DM" & Year == 2009 & data_source == "report" ))
all <- subset(all, !(iso3 == "TGO" & site == "USP KATINDI" & Year == 2011 & data_source == "report" ))
all <- subset(all, !(iso3 == "CMR" & site == "Ndop" & Year == 2002 & data_source == "report" ))
all <- subset(all, !(iso3 == "CPV" & site == "S.Vicente" & Year == 1995 & data_source == "report" ))
all <- subset(all, !(iso3 == "ERI" & site == "Akordet" & Year == 2007 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "SAINT LOUIS" & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "DAKAR" & Year == 1996 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "DAKAR" & Year == 1997 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "LOUGA" & Year == 1998 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "LOUGA" & Year == 2000 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "THIES" & Year == 1998 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "LOUGA" & Year == 2000 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "FATICK" & Year == 2004 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "KOLDA" & Year == 2004 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "KAOLACK" & Year == 1989 & data_source == "report" ))
all <- subset(all, !(iso3 == "SEN" & site == "ZIGUINCHOR" & Year == 2004 & data_source == "report" ))
all <- subset(all, !(iso3 == "LBR" & site == "Senji Health Center" & Year == 2007 & data_source == "report" ))
all <- subset(all, !(iso3 == "LBR" & site == "Senji Health Center" & Year == 2007 & data_source == "report" ))
all <- subset(all, !(iso3 == "TZA" & group == "Zanzibar" & data_source == "report" ))


## Data removal after having undue influence on trajectory of model-----------------------------------------------------------------------------------------------

# Dropped all ANC data from BDI for 2006
all <- all[!(iso3=='BDI' & year==2006)]

# Dropped all ANC data from CIV for 2008
all <- all[!(iso3=='CIV' & Year==2008)]

# Dropped the following specific site-Years
all <- all[!(iso3=='ERI' & Year %in% c(2003, 2005) & site=='Asseb')]
all <- all[!(iso3=='KEN' & Year==2000 & site %in% c('Mbale', 'Meru', 'Thika'))]
all <- all[!(iso3=='KEN' & Year==2001 & site=='Fatima')]
all <- all[!(iso3=='UGA' & Year==2011 & site=='Lira Hosp')]
all <- all[!(iso3=='BEN' & Year==2003 & site %in% c('Angaradébou', 'Kotopounga'))]
all <- all[!(iso3=='BEN' & Year==2005 & site=='Comè')]
all <- all[!(iso3=='BEN' & Year==2008 & site=='Abomey−Calavi / Godomey')]
all <- all[!(iso3=='CMR' & Year==2009 & site=='Sangmelima')]
all <- all[!(iso3=='GHA' & Year %in% c(2003, 2011) & site=='Cape Coast')]
all <- all[!(iso3=='GIN' & Year==2004 & site=='GOUECKE')]
all <- all[!(iso3=='NGA' & Year==2008 & site %in% c('Sokoto (Dogon Daji)', 'Yola', 'Abuja Bwari'))]
all <- all[!(iso3=='ETH' & Year %in% c(2000, 2001) & site=='Adama HC')]
all <- all[!(iso3=='KEN' & Year==2000 & site %in% c('Chulaimbo', 'Kisumu', 'UON4 - Baba Dogo'))]
all <- all[!(iso3=='MOZ' & Year==2000 & site %in% c('9 - Ponta-Gea', '11 - Mondlane', '12 - No. 3'))]
all <- all[!(iso3=='TZA' & Year==2000 & site %in% c('Mbeya (Itete)', 'Mbeya (Kiwanjampaka)', 'Mbeya (Kyela)', 'Mbeya(Ruanda)'))]
all <- all[!(iso3=='TZA' & Year==2007 & site %in% c('Iringa (Ngome)', 'Vwawa Hospital'))]
all <- all[!(iso3=='SLE' & Year==2006 & site %in% c('Makeni', 'Pujehun'))]
all <- all[!(iso3=='SLE' & Year==2007 & site %in% c('PCMH', 'Kissy Urban centre', 'Koidu'))]
all <- all[!(iso3=='SLE' & Year==2008 & site %in% c('Kambia', 'Mattru UBC'))]
all <- all[!(iso3=='LBR' & Year==2006 & site %in% c('Firestone Hospital', 'Phebe Hospital', 'Liberia Govt Hosp−Buchanan', 'Martha Tubman Hospital'))]
all <- all[!(iso3=='LBR' & Year==2007 & site %in% c('JJ Dossen Hospital', 'Senji Health Center', 'Karnplay Health Center', 'Voinjama Health Center', 'Fish Town Health Center'))]
all <- all[!(iso3=='GMB' & Year %in% c(2006, 2007) & site %in% c('Farafenni', 'Sibanor', 'Kuntaur',  'Basse', 'Brikama', 'Sibanor', 'Essau'))]
all <- all[!(iso3=='BFA' & Year==2008 & site == 'Sindou')]
all <- all[!(iso3=='BFA' & Year==2009 & site == 'Ouahigouya')]
all <- all[!(iso3=='CMR' & Year %in% c(2008, 2009) & site == 'Sangmelima')]
all <- all[!(iso3=='CMR' & Year==2012 & site %in% c('HR Maroua', 'HR Ebolowa', 'Doual (HD Nylon)', 'Douala (HD Congo II)', 'Pouma', 'Saa', 'Ngaoundal', 'Molokol (Mokong)'))]
all <- all[!(iso3=='GHA' & Year==2008 & site == 'Nadowli')]
all <- all[!(iso3=='NGA' & Year==2005 & site %in% c('Abuja (Garki/Wuse)', 'Mubi, Lagos (Badagry)', 'Lagos (Ikeja)'))]
all <- all[!(iso3=='NGA' & Year==2008 & site %in% c('Taraba (Zing)', 'Sokoto (Tambawal)', 'Borno (Maiduguri)', 'Kaduna (Kafanchan)', 'Plateau (Jos)', 'Abuja Bwari'))]
all <- all[!(iso3=='MDG' & Year==2018 & site == 'Sainte Marie')]
all <- all[!(iso3=='COD' & Year==2007 & site %in% c('KINGASANI', 'BUKAVU', 'MBANDAKA', 'LUBUMBASHI', 'BINZA METEO'))]
all <- all[!(iso3=='COD' & Year==2011 & site %in% c('LUKULA', 'LODJA'))]
all <- all[!(iso3=='GAB' & Year==2002 & site %in% c('SMI CS Glass 015', 'SMI CS Louis 011'))]
all <- all[!(iso3=='GIN' & Year==2004 & site %in% c('LAFOU', 'KOULE', 'LEY SARE'))]
all <- all[!(iso3=='GMB' & Year==2005 & site %in% c('Basse', 'Essau'))]
all <- all[!(iso3=='KEN' & Year==2001 & site %in% c('Mt. Elgon', 'Kitui', 'Fatima'))]
all <- all[!(iso3=='KEN' & Year==2011 & site %in% c('Lodwar', 'Tabaka'))]
all <- all[!(iso3=='SEN' & Year==2004 & site =='DAKAR')]
all <- all[!(iso3=='TCD' & Year==2011 & site %in% c('Sarh', 'Pala', 'Roi Fayçal'))]
all <- all[!(iso3=='UGA' & Year==2005 & site =='Kagadi Hosp')]
all <- all[!(iso3=='UGA' & Year==2010 & site %in% c('Kagadi Hosp', 'Mbarara Hosp'))]
all <- all[!(iso3=='UGA' & Year==2003 & site %in% c('Mbale Hosp', 'Kilembe Hosp'))]
all <- all[!(iso3=='UGA' & Year==2006 & site %in% c('Tororo Hosp', 'Nyakibaale Hosp', 'Nsambya Hospital'))]
all <- all[!(iso3=='UGA' & Year==2012 & site =='Jinja Hosp')]


# Only keep variables we need
all <- all[,c("nid", "iso3", "group", "site", "point", "latitude","longitude",
              "admin_level","shapefile","location_code","Year","Prev", "N","data_source")]

## Separate geomatched and non-geomatched data -----------------------------------------------------------------------------------------------------------
all_mapped <- subset(all,!is.na(point)) #5745 site-years dropped (4659 site-years post-2000)
locations_mapped <- all_mapped[,1:4]
locations_mapped <- unique(locations_mapped)
match_to_do <- subset(all, is.na(point) & !(iso3 %in% c("ZMB","NAM")))
match_to_do <- match_to_do[,1:4]
match_to_do <- unique(match_to_do)


## Rename variables and change data type to prep for collapse -------------------------------------------------------------------------------------------
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


## Fake Sample Size Checks -------------------------------------------------------------------------------------------------------------------------------

# For countries with many repeated sample sizes at a set level,
# change it to be the median of the others

all_mapped <- all_mapped
setDT(all_mapped)

c_all <- data.table()
y_all <- data.table()
n_all <- data.table()

# ID the unique site-years with 1 sample size only
for(c in unique(all_mapped$country)){
  print(paste0('In ', c, ' ...'))
  for(y in unique(all_mapped[country==c]$year)){
    n <- unique(all_mapped[country == c & year == y]$N)
    x <- nrow(all_mapped[country == c & year == y])
    if(x > 1 & length(n) < 2 & c != 'UGA'){
      print(paste0('Only 1 sample size (', n, ') in year: ', y, ' for ', x, ' sites'))

      c_all <- rbind(c_all, c)
      y_all <- rbind(y_all, y)
      n_all <- rbind(n_all, n)

    } else if(c == 'UGA' & y %in% 1989:2004) {
      n <- 300
      x <- nrow(all_mapped[country == c & year == y & site != 'Lacor Hosp'])
      print(paste0('Only 1 sample size (', n, ') in year: ', y, ' for ', x, ' sites'))

      c_all <- rbind(c_all, c)
      y_all <- rbind(y_all, y)
      n_all <- rbind(n_all, n)
    }
  }
}

# Bind unique site-year-sample size values
df <- cbind(c_all, y_all, n_all)
names(df) <- c('country', 'year', 'n')

## ------------------------------------------------------------------------------------------------------------------------------------------------

# Make plots of the country-sites with fake sample sizes to help in making final decisions on which to change to median & which to leave as is
plots<-list()
i<-0
for(c in unique(df$country)){
  print(paste0('In ', c, ' ...'))
  years <- unique(df[country == c]$y)
  print(years)
  for(s in unique(all_mapped[country == c & year %in% years]$site)) {
    if(c == 'UGA' & s == 'Lacor Hosp') next()
    z <- nrow(all_mapped[country == c & site == s & !(year %in% years)])
    if(z > 2) {
      tmp <- all_mapped[country == c & site == s]
      tmp <- tmp[order(year),]
      tmp[year %in% years, fake := 'yes']
      tmp[!(year %in% years), fake := 'no']
      a<-ggplot(tmp, aes(x=year, y=N, color = fake))+
        geom_point(size=3)+
        theme_classic(base_size = 12)+
        labs(title = paste0(c, ' - ', s), x = "Year", y = "Sample Size")
      i <- i + 1
      plots[[i]] <- a
    }
  }
}

pdf("<<<< FILEPATH REDACTED >>>>",
    width = 16, height = 8)
par(mfrow = c(2,3))
for(i in 1:(length(plots)-5)){
  if(i != 1 & !((i-1) %% 6 == 0)) next()
  grid.arrange(plots[[i]], plots[[i+1]], plots[[i+2]], plots[[i+3]],
               plots[[i+4]], plots[[i+5]],
               ncol=3)
}
dev.off()
## -----------------------------------------------------------------------------------------------------

# In certain cases pre-selected based on examination of the above plots, re-assign fake sample sizes as site-level median
fake_countries <-        c('BDI', 'CIV', 'CMR', 'COG',
                           'GIN', 'LSO', 'MLI', 'MOZ',
                           'NAM', 'RWA', 'TGO', 'UGA',
                           'ZWE')

all_mapped[, change_to_median := FALSE] #Will be used to ID those sites chosen for change
changed_printout <- data.table()
for(c in fake_countries){
  print(c)
  years <- unique(df[country == c]$y)
  print(years)

  # For countries where all relevant sites are changed to the median
  if(c %in% c('BDI', 'COG', 'LSO', 'MLI', 'MOZ', 'RWA', 'TGO')) {
    all_mapped[country == c & year %in% years, change_to_median := TRUE]
  }
  # For countries where only a smaller proportion of the sites should be changed to the median
  if(c == 'CIV'){
    all_mapped[country == c & year %in% years & site == 'CSR ATTINGUIE', change_to_median := TRUE]
  }
  if(c == 'ZWE'){
    all_mapped[country == c & year %in% years & site %in% c('Mnene/Musume Hosp/Mberengwa', 'Gweru City', 'Beitbridge', 'Bindura (Chipadze)', 'Banket'), change_to_median := TRUE]
  }
  # For countries where only a smaller proportion of the sites should NOT be changed to the median
  if(c == 'CMR'){
    all_mapped[country == c & year %in% years & site != 'Fondation Chantal Biya', change_to_median := TRUE]
  }
  if(c == 'GIN'){
    all_mapped[country == c & year %in% years & !(site %in% c('BOULBINET', 'KAMSAR')), change_to_median := TRUE]
  }
  if(c == 'NAM'){
    all_mapped[country == c & year %in% years & !(site %in% c('Walvisbay', 'Rundu', 'Katutura', 'Engela', 'Oshikuku', 'Oshakati', 'Katima')), change_to_median := TRUE]
  }
  if(c == 'UGA'){
    all_mapped[country == c & year %in% years & !(site %in% c('Kaabong hosp', 'Lacor Hosp', 'Matany Hosp')), change_to_median := TRUE]
  }

  for(s in unique(all_mapped[country == c & year %in% years & change_to_median == TRUE]$site)) {
    z <- nrow(all_mapped[country == c & site == s & !(year %in% years)])
    if(z > 2) {
      m <- median(all_mapped[country == c & site == s & !(year %in% years)]$N)
      f <- unique(all_mapped[country == c & site == s & year %in% years]$N)

      print(paste0('Changing fake sample size of ', f, ' in ', c, ' - ', s, ' to ', m, '. Fake N - median N == ', f-m))
      q <- data.table()
      q$year <- years
      q$fss <- f
      q$site <- s
      q$country <- c
      q$nss  <- m

      changed_printout <- rbind(changed_printout, q)

      all_mapped[country == c & site == s & (year %in% years), N := m]

    }
  }
  write.csv(changed_printout,"<<<< FILEPATH REDACTED >>>>")
}

all_mapped <- data.frame(all_mapped)

## Save files -------------------------------------------------------------------------------------------------------------------------------
module_date <- Sys.Date()
module_date <- gsub("-", "_", module_date)
saveRDS(all_mapped, file=paste0("<<<< FILEPATH REDACTED >>>>"))
