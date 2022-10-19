########### Merge WHO MDA data to new ESPEN IU shapefile
library(data.table)
library(raster)
library(rgdal)
library(stringr)
library(utf8)
library(sp)
library(rgeos)
library(sf)

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

IU_shp <- readOGR(<<<< FILEPATH REDACTED >>>>)
MDA <- fread(<<<< FILEPATH REDACTED >>>>)
MDA_mod <- fread(<<<< FILEPATH REDACTED >>>>)
setnames(MDA_mod, "iu_name", "iu_name_modified")
MDA_mod <- merge(MDA_mod, MDA[, c("occ_id_num", "iu_name")], by = "occ_id_num", all.x = TRUE)

shp_data <- IU_shp@data
shp_data$ADMIN1 <- str_to_title(shp_data$ADMIN1)
shp_data$IUs_NAME <- str_to_title(shp_data$IUs_NAME)

MDA_mod$iu_name_modified <- str_to_title(MDA_mod$iu_name_modified)
MDA_mod$region <- str_to_title(MDA_mod$region)

MDA_mod <- MDA_mod[iso3 %in% unique(shp_data$ADMIN0ISO3)]

# Identify IUs with the same name in a given country but different states (use the original MDA file due to some intentional duplications in MDA_mod?)
agg_st <- as.data.table(aggregate(occ_id_num ~ iso3 + iu_name, data = MDA, length))
agg_st <- agg_st[occ_id_num > 1]
agg_st[, occ_id_num := NULL]
agg_st <- merge(agg_st, MDA[, c("iso3", "iu_name", "occ_id_num")], all.x = TRUE)

multis <- MDA_mod[occ_id_num %in% unique(agg_st$occ_id_num)]
MDA_mod <- MDA_mod[!(occ_id_num %in% unique(agg_st$occ_id_num))]

MDA_mod_merged <- merge(MDA_mod, shp_data, by.x = c("iso3", "iu_name_modified"), by.y = c("ADMIN0ISO3", "IUs_NAME"), all.x = TRUE)

# Manually update IDs for some rows that failed to merge
MDA_mod_merged[occ_id_num == 1927, IU_ID := "ETH0195019271"]
MDA_mod_merged[occ_id_num == 354, IU_ID := "GHA0216521486"]
MDA_mod_merged[occ_id_num == 344, IU_ID := "GHA0216721533"]
MDA_mod_merged[occ_id_num == 355, IU_ID := "GHA0216621523"]
MDA_mod_merged[occ_id_num == 1718, IU_ID := "GHA0217021599"]
MDA_mod_merged[occ_id_num == 2748, IU_ID := "GNB0241023665"]
MDA_mod_merged[occ_id_num == 2751, IU_ID := "GNB0241023667"]
MDA_mod_merged[occ_id_num == 361, IU_ID := "GNB0241023669"]
MDA_mod_merged[occ_id_num == 2115, IU_ID := "MDG0287127887"]
MDA_mod_merged[occ_id_num == 2808, IU_ID := "MDG0286327848"]
MDA_mod_merged[occ_id_num == 721, IU_ID := "MDG0288027928"]
MDA_mod_merged[occ_id_num == 2082, IU_ID := "MDG0286527858"]
MDA_mod_merged[occ_id_num == 2078, IU_ID := "MDG0286427851"]
MDA_mod_merged[occ_id_num == 831, IU_ID := "MOZ0345333239"]
MDA_mod_merged[occ_id_num == 830, IU_ID := "MOZ0345133215"]
MDA_mod_merged[occ_id_num == 2119, IU_ID := "MOZ0344633139"]
MDA_mod_merged[occ_id_num == 2120, IU_ID := "MOZ0345033194"]
MDA_mod_merged[occ_id_num == 2500, IU_ID := "MOZ0344333108"]
MDA_mod_merged[occ_id_num == 219, IU_ID := "MOZ0345333251"]
MDA_mod_merged[occ_id_num == 2121, IU_ID := "MOZ0344933172"]
MDA_mod_merged[occ_id_num == 885, IU_ID := "NER0367235317"]
MDA_mod_merged[occ_id_num == 876, IU_ID := "NER0366835299"]
MDA_mod_merged[occ_id_num == 898, IU_ID := "NER0367435333"]
MDA_mod_merged[occ_id_num == 1257, IU_ID := "NGA0378736580"]
MDA_mod_merged[occ_id_num == 2892, IU_ID := "NGA0380436955"]
MDA_mod_merged[occ_id_num == 2156, IU_ID := "SDN53215"]
MDA_mod_merged[occ_id_num == 347, IU_ID := "SDN53217"]
MDA_mod_merged[occ_id_num == 348, IU_ID := "SDN53219"]
MDA_mod_merged[occ_id_num == 2675, IU_ID := "SEN0412540182"]
MDA_mod_merged[occ_id_num == 2738, IU_ID := "SEN0413240225"]
MDA_mod_merged[occ_id_num == 1491, IU_ID := "TZA0479646520"]
MDA_mod_merged[occ_id_num == 239, IU_ID := "GNB0241023662"]

# Now add IU_IDs to the IUs with duplicate names in a given country
multis[occ_id_num == 1143, "IU_ID" := "NGA0379836835"]
multis[occ_id_num == 2717, IU_ID := "NGA0380436957"]
multis[occ_id_num == 228, IU_ID := "NGA0380036870"]
multis[occ_id_num == 340, IU_ID := "NGA0379436736"]
multis[occ_id_num == 1294, IU_ID := "NGA0380036872"]
multis[occ_id_num == 2727, IU_ID := "NGA0378136468"]

# Append multis to MDA_mod_merged
MDA_mod_merged <- rbindlist(list(MDA_mod_merged, multis), use.names = TRUE, fill = TRUE)

# Merge MDA_mod_merged data onto IU_shp
IU_shp_temp <- merge(IU_shp, MDA_mod_merged, by = "IU_ID")
IU_shp_temp <- IU_shp_temp[!is.na(IU_shp_temp@data$occ_id_num),]

# Update shapefile paths
MDA_mod_merged[data_path == <<<< FILEPATH REDACTED >>>>, data_path := <<<< FILEPATH REDACTED >>>>]
MDA_mod_merged$data_path <- gsub(<<<< FILEPATH REDACTED >>>>, <<<< FILEPATH REDACTED >>>>, MDA_mod_merged$data_path)

# Load shapefiles for georeferenced locations in the MDA data set that are not matched to a polygon in the IU shapefile
shapes <- unique(MDA_mod_merged[is.na(IU_ID) & !is.na(data_path), data_path])
shapes <- shapes[!(shapes == "")]

shp_list <- vector(mode = "list", length = length(shapes))
for (i in 1:length(shapes)) {
  try(shp_list[i] <- readOGR(shapes[[i]]))
}

# Get list of occ_id_num to match to shapefiles
id_list <- MDA_mod_merged[is.na(IU_ID) & !(is.na(data_path) | data_path == ""),]
id_list <- id_list[!(occ_id_num %in% c(1767, 1798))] # Drop "duplicate" rows (same georeferencing and MDA history)

# Subset shapefiles
for (i in 1:length(shapes)) {
  current <- shp_list[[i]]
  current <- merge(current, id_list[data_path == shapes[[i]]], by.x = c(unique(id_list[data_path == shapes[[i]], feature_field])), by.y = "feature_code")
  current <- current[!is.na(current$occ_id_num),]
  shp_list[[i]] <- copy(current)
}

# Adjust shapefile settings to accommodate raster::erase
IU_shp_template <- copy(IU_shp_temp)
IU_shp_template <- spTransform(IU_shp_template, CRS("+proj=aea"))
IU_shp_template <- gBuffer(IU_shp_template, byid = T, width = 0)

# Add any other polygons that do not overlap with the current shapes in IU_shp_template
for (i in 1:length(shapes)) {
  print(i)
  current <- shp_list[[i]]
  current <- spTransform(current, CRS("+proj=aea"))
  current <- gBuffer(current, byid = T, width = 0)
  
  for (j in 1:length(current)) {
    print(j)
    if (length(st_intersects(st_as_sf(current[j]), st_as_sf(IU_shp_template))) == 0) {
      print("Binding...")
      IU_shp_template <- raster::bind(IU_shp_template, current[j])
    } else {
      print("Overlaps...")
    }
  }
}

# Set CRS back to original
IU_shp_template <- spTransform(IU_shp_template, CRS("+proj=longlat"))
writeOGR(IU_shp_template, <<<< FILEPATH REDACTED >>>>, "filename", driver="ESRI Shapefile")
