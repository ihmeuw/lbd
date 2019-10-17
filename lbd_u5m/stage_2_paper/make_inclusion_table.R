## #############################################################################
##
## STAGE 2 PAPER: Supplemental Figures - Inclusion and Exclusion tables
##
## Date: January 18, 2019
##
## #############################################################################

library(raster)
library(data.table)
library(rgdal)
library(rgeos)
library(dplyr)
library(openxlsx)
library(foreign)
library(parallel)
library(stringi)
library(gridExtra)

#update with core repo path
source("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")

#args for pulling source bias
regions <- c('soas', 'ocea+seas-mys', 'stan+mng', 'caca-mex', 'ansa+trsa-bra', 'noaf', 'mide+yem', 'cssa', 'essa-yem', 'sssa', 'wssa')
run_date <- "<<<< FILEPATH REDACTED >>>>"
plot <- T
plot_path <- "<<<< FILEPATH REDACTED >>>>"
cores <- 11


check <- readRDS("<<<< FILEPATH REDACTED >>>>")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Inclusion Table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input_nids <- c()
for(reg in regions) {
  input_data <- fread("<<<< FILEPATH REDACTED >>>>")
  input_nids <- c(input_nids, unique(input_data$nid))
}

cbh <- readRDS("<<<< FILEPATH REDACTED >>>>")
sbh <- readRDS("<<<< FILEPATH REDACTED >>>>")

meta <- fread("<<<< FILEPATH REDACTED >>>>")

stage_list <- fread("<<<< FILEPATH REDACTED >>>>")
stage_list <- stage_list[,c("location_name","iso3", "mbg_reg")]

cbh_agg <- cbh[, .(children_born = .N), by = nid]
sbh_agg <- sbh[, .(children_born = sum(ceb, na.rm =T)), by = nid]
agg <- rbind(cbh_agg, sbh_agg, fill =T)

poly_cbh <- cbh[, .(polygons = length(na.omit(as.integer(unique(location_code))))), by = nid]
point_cbh <- unique(cbh[point == 1, c("nid", "latnum", "longnum")])
point_cbh <- point_cbh[, .(points = .N), by = nid]

poly_sbh <- sbh[, .(polygons = length(na.omit(as.integer(unique(location_code))))), by = nid]
point_sbh <- unique(sbh[point == 1, c("nid", "latnum", "longnum")])
point_sbh <- point_sbh[, .(points = .N), by = nid]

poly <- rbind(poly_sbh, poly_cbh)
point <- rbind(point_sbh, point_cbh)

inclusion_table <- merge(meta, agg, by = "nid", all = T)
inclusion_table <- merge(inclusion_table, poly, by = "nid", all = T)
inclusion_table <- merge(inclusion_table, point, by = "nid", all = T)
inclusion_table <- merge(inclusion_table, stage_list, by.x = "country", by.y = "iso3", all.x = T)

inclusion_table[is.na(points), points := 0 ]
inclusion_table[is.na(polygons), polygons := 0 ]
inclusion_table <- subset(inclusion_table, select = -c(sample, rows, full_title, country))
inclusion_table <- inclusion_table[nid %in% input_nids, ]
inclusion_table <- inclusion_table[!(location_name %in% c("Brazil", "Mexico", "China", "Malaysia")), ]
#subset regions
inclusion_table <- inclusion_table[!(mbg_reg %in% c("", "soax", "ocex")),]

#some excluded surveys snuck in
inclusion_table <- inclusion_table[!(nid %in% c(1175, 336042)),]

#survey counts
nrow(inclusion_table[type == "SBH",])
nrow(inclusion_table[type == "CBH",])

#swapping file nids to survey nids to pull citations properly
inclusion_table[nid == 237958, nid := 22125]
inclusion_table[nid == 19932, nid := 19950]
inclusion_table[nid == 6825, nid := 6842]
inclusion_table[nid == 23165, nid := 23183]

# add in citation
nids <- inclusion_table$nid
ghdx_meta <- ghdx_construct_pub_table(
  nids         = nids,
  core_repo    = "<<<< FILEPATH REDACTED >>>>",
  resolve_urls = F
)

ghdx_meta <- ghdx_meta[, c("nid", "citation")]
inclusion_table <- merge(inclusion_table, ghdx_meta, by = "nid", all.x = T)

inclusion_table$citation <- gsub(pattern = "<p>|</p>", replacement = "", x=inclusion_table$citation)

#manually cite
inclusion_table[nid == 249499, citation := "Ministry of Planning and International Cooperation (Yemen), International Policy Center for Inclusive Growth, Interaction in Development (Yemen), UNICEF Yemen. Yemen National Social Protection Monitoring Survey 2012-2013."]

#format
setcolorder(inclusion_table, c("location_name", "year", "type", "children_born", "polygons", "points", "mbg_reg", "citation", "nid"))
names(inclusion_table) <- c("Country", "Year", "Data Type", "Children Born", "Polygons", "Points", "modelling region", "Citation", "GHDx ID")
setorder(inclusion_table, Country, Year)
inclusion_table <- inclusion_table[, c(1:6,8:9)]

#move pre-2000 data to excluded
moved_data <- inclusion_table[Year < 2000, ]
inclusion_table <- inclusion_table[Year >= 2000,]

stri_enc_mark(inclusion_table$Citation)
for (col in colnames(inclusion_table)){
  Encoding(inclusion_table$Citation) <- "UTF-8"
}

inclusion_table$Citation <- iconv(inclusion_table$Citation, from = 'UTF-8', to = 'ASCII//TRANSLIT')

sink("<<<< FILEPATH REDACTED >>>>", append = T)
cat("\n\n#Precision public health and child mortality, last paragraph middle\n")
cat("To produce our estimates, we collated", nrow(inclusion_table), "georeferenced household surveys and censuses")
sink()

inclusion_table <- unique(inclusion_table)

#save
write.xlsx(x = inclusion_table, file = "<<<< FILEPATH REDACTED >>>>")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exclusion Table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

excluded <- fread("<<<< FILEPATH REDACTED >>>>")

#add in citation
nids <- excluded$nid
ghdx_meta <- ghdx_construct_pub_table(
  nids         = nids,
  core_repo    = "<<<< FILEPATH REDACTED >>>>",
  resolve_urls = F
)

ghdx_meta <- ghdx_meta[, c("nid", "citation")]
exclusion_table <- merge(excluded, ghdx_meta, by = "nid", all.x =T)

exclusion_table$citation <- gsub(pattern = "<p>|</p>|<span style=\"font-size: 12px;\">|</span>", replacement = "", x=exclusion_table$citation)

exclusion_table[nid == 197909, citation := "Performance Monitoring and Accountability 2020 (PMA2020) Project. Baltimore, MD: PMA2020, Bill & Melinda Gates Institute for Population and Reproductive Health, Johns Hopkins Bloomberg School of Public Health."]
exclusion_table[nid == 244469, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Aden Nutrition Survey 2012."]
exclusion_table[nid == 244471, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Hajjah Lowland and Mountainous Ecological Zones Nutrition Survey 2012."]
exclusion_table[nid == 244473, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Taiz Mountainous and Coastal Plain Ecological Zones Nutrition Survey 2012."]
exclusion_table[nid == 246249, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Ibb Nutrition and Mortality Survey 2012."]
exclusion_table[nid == 246254, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF, Relief International. Yemen - Lahj Nutrition and Mortality Survey 2012."]
exclusion_table[nid == 246145, citation := "International Organization for Migration (IOM), International Rescue Committee (IRC), Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Abyan Nutritional Status and Mortality Survey 2012-2013."]
exclusion_table[nid == 244472, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Rayma Nutrition Survey 2012."]
exclusion_table[nid == 246209, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Dhamar Nutritional Status and Mortality Survey 2013."]
exclusion_table[nid == 246250, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Mahweet Nutrition and Mortality Survey 2013."]
exclusion_table[nid ==246246, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Hajjah Nutrition and Mortality Survey 2015."]
exclusion_table[nid ==246248, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Hodeidah Nutrition and Mortality Survey 2015."]
exclusion_table[nid ==244463, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF, Field Medical Foundation. Yemen - Aden Nutrition and Mortality Survey 2015."]
exclusion_table[nid ==244464, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Al-Baidha Nutrition and Mortality Survey 2015."]
exclusion_table[nid ==244465, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF. Yemen - Hajjah Nutrition and Mortality Survey 2015."]
exclusion_table[nid ==244467, citation := "European Commission for Humanitarian Aid and Civil Protection (ECHO), UNICEF Yemen Country Office, Ministry of Public Health and Population (Yemen). Yemen - Hodeidah Nutrition and Mortality Survey 2015."]
exclusion_table[nid ==244468, citation := "Ministry of Public Health and Population (Yemen), Yemen UNICEF, Humanitarian Aid and Development Organization (HAD). Yemen - Lahj Nutrition and Mortality Survey in Lowland and Highlands Ecological Zones 2015."]

exclusion_table <- merge(exclusion_table, stage_list, by.x = "country", by.y = "iso3")

setorder(exclusion_table, location_name, year)
exclusion_table <- subset(exclusion_table, select = -c(full_title, mbg_reg, country))

#rearrange rows
exclusion_table <- exclusion_table[, c("location_name", "year", "type", "reason", "citation", "nid")]
names(exclusion_table) <- c("Country", "Year", "Data Type", "Reason", "Citation", "NID")

exclusion_table <- exclusion_table[!(Country %in% c("Brazil", "Mexico", "China", "Malaysia")), ]

sink("<<<< FILEPATH REDACTED >>>>", append = T)
cat("\n\n#Discussion, 2nd to last\n")
cat("All input data were subject to quality checks, which resulted in the exclusion of", nrow(exclusion_table[Year >= 2000]), "surveys and censuses due to quality concerns (see Supplemental information for more details)")
sink()


#save
write.xlsx(x = exclusion_table[Year >= 2000], file = "<<<< FILEPATH REDACTED >>>>")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Source Bias Plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Making Combined Source Bias Plot - Pulling grob from reg: ")
grob_list <- list()
for(reg in regions) {
  message(reg)
  source_bias <- readRDS(sprintf("<<<< FILEPATH REDACTED >>>>"))
  source_bias <- source_bias + theme(axis.title.x=element_blank()) + 
    theme(plot.title = element_text(hjust = 0.5))
  grob_list <- c(grob_list, list(source_bias))
}

png(filename = sprintf("<<<< FILEPATH REDACTED >>>>"),
    width=8,
    height=12,
    units = "in",
    res = 300)
grid.arrange(grobs = grob_list, ncol=3, nrow=4)
dev.off()
