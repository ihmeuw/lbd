############################################################################
#process for preparing the admin 2 level populations
############################################################################


#######
# Set up Environment
#######

## Clear environment
rm(list = ls())


## Set indicator
indicator_group <- 'hiv'
indicator       <- 'epp'

## Set repos
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_hiv/")
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>", "/HIV/")
setwd(core_repo)


new_pkg_lib   <- paste0("'<<<< FILEPATH REDACTED >>>>", "/r_packages/")
if (!dir.exists(new_pkg_lib)) {dir.create(new_pkg_lib)}

.libPaths(new_pkg_lib)
test_pkg_list <- c('slackr', "RMariaDB", "DBI", "RMySQL", 'rlang')
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}
library(slackr, lib = new_pkg_lib)
library(rlang, lib = new_pkg_lib)


## Load libraries and  MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = core_repo)

## Packages
library(data.table)



### Functions
library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>/02_mortality/shared/r")


### Tables
loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

### Code
epp.list <- sort(loc.table[epp == 1, ihme_loc_id])
loc.list <- epp.list
loc.list <- c(loc.list, "ZAF", "NGA", "KEN", "ETH", "IND", "MAR", "VNM")
loc.list <- as.list(loc.list)

loc.list <- loc.list[which(nchar(loc.list) <= 3)] #only want to do this for countries not all of the GBD subnationals

shapefile_version <- '<<<< FILEPATH REDACTED >>>>'
pop_release <- "<<<< FILEPATH REDACTED >>>>"
modeling_shapefile_version <- shapefile_version
dir.create(paste0("<<<< FILEPATH REDACTED >>>>"))
year_list <- c(1969:2019)
interval_mo <- 12
rake_subnational <- F # needs to be fed in
gbd_round <- 6
decomp_step <- "step4"
age_groups <- c("a1519f", "a1519m", "a2024f", "a2024m", "a2529f", "a2529m", "a3034f", "a3034m", "a3539f", "a3539m", "a4044f", "a4044m", "a4549f", "a4549m", "a5054f", "a5054m")
n_g <- length(age_groups)


#########################################################
for (r in loc.list) {

  if (r %in% c("ETH", "ZAF", "KEN", "NGA")) {

  message(paste0(r))
  print(r)
  r
  rake_subnational <- T
  reg <- r #needs to be fed in
  #######
  # Load spatial objects
  #######
  #load regional simple raster
  simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = shapefile_version), buffer = 0.4, subset_only = FALSE, shapefile_version = shapefile_version)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon = simple_polygon_list[[2]]
  message("Building simple raster from subset_shape")
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  rm(raster_list,simple_polygon_list);gc()

  pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  #######
  # Load link table
  #######
  link_table <- get_link_table(simple_raster = simple_raster, shapefile_version = shapefile_version)



  #######
  # Load WorldPop Populations
  #######
  wps <- list()



  for (i in 1:n_g) {
    a <- age_groups[[i]]
    wps[[i]] <- load_populations_cov(reg, pop_measure = a, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
  }


  data <- wps[[1]][ ,c("pixel_id", "year", "pop")]

  for (i in 2:n_g) {
    data <- cbind(data, wps[[i]]$pop)
  }

  names(data) <- c("pixel_id", "year", age_groups)
  wps <- NULL



  #######
  # Setting up Link Table Merge
  #######
  link <- prep_link_table(link_table = link_table,
                          simple_raster = simple_raster,
                          pixel_id = pixel_id)

  cell_ids <- link_table[[2]]

  # getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
  # gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary

  connector <- get_gbd_locs(rake_subnational = rake_subnational,
                            reg = reg,
                            shapefile_version = shapefile_version)

  # merge the connector on to the link table by making sure that each cell fragment gets connected to the appropriate
  # raking geography
  link <- sub_nat_link_merge(rake_subnational,
                             link,
                             connector)

  data <- as.data.table(data)
  data[, data_id := .I] #cell_pred object ID

  cell_ids_long <- rep(cell_ids, length(year_list))

  data[,cell_id := cell_ids_long] #cell id references the country map



  #######
  # Merging on to Link Table
  #######

  data <- merge(link, data, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)
  data <- as.data.table(data)
  # space
  link <- NULL

  #######
  # Convert to fractional population
  #######

  #convert to fractional population
  data = data[,(age_groups) := lapply(age_groups, function(x) get(x)*area_fraction)]

  #######
  # Calculate Pop Scalars
  #######

  #do calculations!
  rake_geo_pop = data[, lapply(age_groups, function(x) sum(get(x), na.rm = T)), by = c('year','location_id')]
  names(rake_geo_pop) <- c("year", "location_id", age_groups)
  loc_ids = unique(connector$location_id)

  #adjust to be GBD populations
  source(paste0("<<<< FILEPATH REDACTED >>>>/get_population.R"))
  message("This script only works for the age groups needed to run EPP, if you are not using sex split 5 year bins from 15 to 54 this raking will fail")
  gbdps <- list()
  f <- grep("f", age_groups, value = TRUE)
  for (i in 1:16) {
    age_n <- age_groups[[i]]
    if (age_n %in% f) {
      sex <- 2
    } else {sex <- 1}
    age <- (((as.numeric(substr(age_n,2,3))) / 5) + 5)
    gbdps[[i]] <- get_population(age_group_id = age,
                                 location_id = loc_ids,
                                 sex = sex,
                                 year_id = unique(year_list),
                                 gbd_round_id = gbd_round,
                                 status = "best",
                                 decomp_step = decomp_step)
  }

  rake_geo_pop <- as.data.frame(rake_geo_pop)

  gbd <- gbdps[[1]][ , c("year_id", "location_id", "population")]

  for (i in 2:16) {
    gbd <- cbind(gbd, gbdps[[i]]$population)
  }

  names(gbd) <- c("year", "location_id", age_groups)

  rake_geo_pop <- as.data.table(rake_geo_pop)
  gbd <- as.data.table(gbd)
  scalar <- merge(rake_geo_pop, gbd, by = c("year", "location_id"))

  rake_geo_melted <- melt(rake_geo_pop, id.vars = c("year", "location_id"))
  gbd_pops_melted <- melt(gbd, id.vars = c("year", "location_id"))

  scalar_melted <- merge(rake_geo_melted, gbd_pops_melted, by = c("year", "location_id", "variable"))
  scalar_melted$scalar <- as.numeric(scalar_melted$value.y) / as.numeric(scalar_melted$value.x)

  data_melted <- melt(data = data, id.vars = c("ADM0_CODE", "ADM1_CODE", "ADM2_CODE", "year", "pixel_id.x", "pixel_id.y", "ID", "data_id"), measure.vars = age_groups)
  scalar_melted <- merge(scalar_melted, connector[, c("location_id", "ADM1_CODE")], by = c("location_id"))

  data_melted <- merge(data_melted, scalar_melted[,c("year", "variable", "scalar", "ADM1_CODE")], by = c("year", "variable", "ADM1_CODE"))
  data_melted$val_raked <- data_melted$value * data_melted$scalar

  data_melted <- data_melted[ ,c("ADM0_CODE", "ADM1_CODE","ADM2_CODE", "year","variable", "val_raked")]

  data_ADM2 <- dcast(data_melted, ADM2_CODE+year ~ variable, fun = sum)

  data_ADM1 <- dcast(data_melted, ADM1_CODE+year ~ variable, fun = sum)
  data_ADM0 <- dcast(data_melted, ADM0_CODE+year ~ variable, fun = sum)


  # arrange the admin_2 pops as needed for EPP
  data_ADM2 <- melt(data_ADM2, id.vars = c("year", "ADM2_CODE"))
  data_ADM2$sex_id <- 1
  data_ADM2$sex_id[which(data_ADM2$variable %in% f)] <- 2
  data_ADM2$age_group_id <- ((as.numeric(substr(data_ADM2$variable, 2,3))/5) + 5)

  data_ADM2 <- data_ADM2[ , c("age_group_id", "ADM2_CODE", "year", "sex_id", "value")]
  names(data_ADM2) <- c("age_group_id", "ADM2_CODE", "year_id", "sex_id", "population")
  data_ADM2$population <- as.numeric(data_ADM2$population)
  write.csv(data_ADM2, paste0("<<<< FILEPATH REDACTED >>>>"))


} else {

message(paste0(r))
print(r)
r
rake_subnational <- F
reg <- r #needs to be fed in
#######
# Load spatial objects
#######
#load regional simple raster
simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = shapefile_version), buffer = 0.4, subset_only = FALSE, shapefile_version = shapefile_version)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon = simple_polygon_list[[2]]
message("Building simple raster from subset_shape")
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
rm(raster_list,simple_polygon_list);gc()

pixel_id <- seegSDM:::notMissingIdx(simple_raster)
#######
# Load link table
#######
link_table <- get_link_table(simple_raster = simple_raster, shapefile_version = shapefile_version)



#######
# Load WorldPop Populations
#######
wps <- list()



for (i in 1:n_g) {
  a <- age_groups[[i]]
  wps[[i]] <- load_populations_cov(reg, pop_measure = a, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
}


data <- wps[[1]][ ,c("pixel_id", "year", "pop")]

for (i in 2:n_g) {
  data <- cbind(data, wps[[i]]$pop)
}

names(data) <- c("pixel_id", "year", age_groups)
wps <- NULL



#######
# Setting up Link Table Merge
#######
link <- prep_link_table(link_table = link_table,
                        simple_raster = simple_raster,
                        pixel_id = pixel_id)

cell_ids <- link_table[[2]]

# getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
# gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary

connector <- get_gbd_locs(rake_subnational = rake_subnational,
                          reg = reg,
                          shapefile_version = shapefile_version)

# merge the connector on to the link table by making sure that each cell fragment gets connected to the appropriate
# raking geography
link <- sub_nat_link_merge(rake_subnational,
                           link,
                           connector)

data <- as.data.table(data)
data[, data_id := .I] #cell_pred object ID

cell_ids_long <- rep(cell_ids, length(year_list))

data[,cell_id := cell_ids_long] #cell id references the africa map



#######
# Merging on to Link Table
#######

data <- merge(link, data, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)
data <- as.data.table(data)
# space
link <- NULL

#######
# Convert to fractional population
#######

#convert to fractional population
data = data[,(age_groups) := lapply(age_groups, function(x) get(x)*area_fraction)]

#######
# Calculate Pop Scalars
#######

#do calculations!
rake_geo_pop = data[, lapply(age_groups, function(x) sum(get(x), na.rm = T)), by = c('year','location_id')]
names(rake_geo_pop) <- c("year", "location_id", age_groups)
loc_ids = unique(connector$location_id)

#adjust to be GBD populations
source(paste0("<<<< FILEPATH REDACTED >>>>/get_population.R"))
message("This script only works for the age groups needed to run EPP, if you are not using sex split 5 year bins from 15 to 54 this raking will fail")
gbdps <- list()
f <- grep("f", age_groups, value = TRUE)
for (i in 1:16) {
  age_n <- age_groups[[i]]
  if (age_n %in% f) {
    sex <- 2
  } else {sex <- 1}
  age <- (((as.numeric(substr(age_n,2,3))) / 5) + 5)
  gbdps[[i]] <- get_population(age_group_id = age,
                               location_id = loc_ids,
                               sex = sex,
                               year_id = unique(year_list),
                               gbd_round_id = gbd_round,
                               status = "best",
                               decomp_step = decomp_step)
}

rake_geo_pop <- as.data.frame(rake_geo_pop)

gbd <- gbdps[[1]][ , c("year_id", "location_id", "population")]

for (i in 2:16) {
  gbd <- cbind(gbd, gbdps[[i]]$population)
}

names(gbd) <- c("year", "location_id", age_groups)

rake_geo_pop <- as.data.table(rake_geo_pop)
gbd <- as.data.table(gbd)
scalar <- merge(rake_geo_pop, gbd, by = c("year", "location_id"))

rake_geo_melted <- melt(rake_geo_pop, id.vars = c("year", "location_id"))
gbd_pops_melted <- melt(gbd, id.vars = c("year", "location_id"))

scalar_melted <- merge(rake_geo_melted, gbd_pops_melted, by = c("year", "location_id", "variable"))
scalar_melted$scalar <- as.numeric(scalar_melted$value.y) / as.numeric(scalar_melted$value.x)

data_melted <- melt(data = data, id.vars = c("ADM0_CODE", "ADM1_CODE", "ADM2_CODE", "year", "pixel_id.x", "pixel_id.y", "ID", "data_id"), measure.vars = age_groups)
names(data_melted)
data_melted$test <- data_melted$pixel_id.y == data_melted$data_id
table(data_melted$test)

data_melted <- merge(data_melted, scalar_melted[,c("year", "variable", "scalar")], by = c("year", "variable"))
data_melted$val_raked <- data_melted$value * data_melted$scalar

data_melted <- data_melted[ ,c("ADM0_CODE", "ADM1_CODE","ADM2_CODE", "year","variable", "val_raked")]

data_ADM2 <- dcast(data_melted, ADM2_CODE+year ~ variable, fun = sum)

data_ADM1 <- dcast(data_melted, ADM1_CODE+year ~ variable, fun = sum)
data_ADM0 <- dcast(data_melted, ADM0_CODE+year ~ variable, fun = sum)


# arrange the admin_2 pops as needed for EPP
data_ADM2 <- melt(data_ADM2, id.vars = c("year", "ADM2_CODE"))
data_ADM2$sex_id <- 1
data_ADM2$sex_id[which(data_ADM2$variable %in% f)] <- 2
data_ADM2$age_group_id <- ((as.numeric(substr(data_ADM2$variable, 2,3))/5) + 5)

data_ADM2 <- data_ADM2[ , c("age_group_id", "ADM2_CODE", "year", "sex_id", "value")]
names(data_ADM2) <- c("age_group_id", "ADM2_CODE", "year_id", "sex_id", "population")
data_ADM2$population <- as.numeric(data_ADM2$population)
write.csv(data_ADM2, paste0("<<<< FILEPATH REDACTED >>>>"))
}
}
