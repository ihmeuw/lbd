
rm(list = ls())

#########
## Set Up
#########
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>/mbg/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

library(dplyr)
library(data.table)
library(lme4)
library(sf)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtools)
require(nnet)
library(spdep)
library(INLA)
library(stringr)

new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>", "/r_packages/")
.libPaths(new_pkg_lib)
test_pkg_list <- c('compositions', 'dplyr')
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}
library(compositions)


'%!in%' <- function(x,y)!('%in%'(x,y))

#######
## set the HIV prevalence extraction to use for modeling ART coverage.
#######
prev_date <- "2020_04_10" #string that indicates where to save the desired outputs and diagnostic plots.
wp_version <- "wp1"

prev_info <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
mbg_run <- as.character(prev_info$Value[which(prev_info$Setting == "mbg_run")])

shp_date <- as.character(prev_info$Value[which(prev_info$Setting == "shapefile_version")])

#####
## Determine what countries we have ART data
#####

count <- list.dirs(paste0("<<<< FILEPATH REDACTED >>>>"), full.names = F, recursive = F)
exclusions <- c("AGO", "BWA", "HTI", "ERI", "IND")
countries <- count[which(count %!in% exclusions)]



#####
## Identify the appropriate file to read in for each country, in this case the newest dated file
#####
files <- list()
for (loc in countries) {
  fs <- list.files(paste0("<<<< FILEPATH REDACTED >>>>"), full.names = FALSE, recursive = FALSE)
  fs <- grep("treat_extract", fs, value = TRUE)
  if (length(fs) == 1) { files[[loc]] <- fs[[1]]}
  else {
    fs <- grep(wp_version, fs, value = TRUE)
    if (length(fs) == 0) {
      fs <- list.files(paste0("<<<< FILEPATH REDACTED >>>>"), full.names = FALSE, recursive = FALSE)
      fs <- grep("treat_extract", fs, value = TRUE)
      }
    fs <- fs[which(substr(fs, 1,2) == "20")]
    dates <- list()
    for (id in c(1:length(fs))) {
      d <- as.Date(str_replace_all(substr(fs[[id]], 1, 10), "_", "/"))
      l2 <- list()
      l2[["date"]] <- d
      l2[["file"]] <- fs[id]
      dates[[id]] <- l2
    }
    dates <- bind_rows(dates)
    files[[loc]] <- dates$file[which(dates$date == max(dates$date))]
  }
}


#######
## create a look up table to connect location code to country name
#######
loc.table <- read.csv("<<<< FILEPATH REDACTED >>>>")

loc.table <- loc.table[which(loc.table$ihme_loc_id %in% countries), ]

c_lookup <- loc.table[ ,c("ihme_loc_id", "location_name", "location_id")]
c_lookup$loc <- c_lookup$ihme_loc_id
c_lookup$ADM0_NAME <- as.character(c_lookup$location_name)

c_lookup$ADM0_NAME[which(c_lookup$ADM0_NAME == "United Republic of Tanzania")] <- "Tanzania"
c_lookup$ADM0_NAME[which(c_lookup$ADM0_NAME == "Congo")] <- "Republic of Congo"


########
## setting loop for countires to collect the data
########
DFs <- list()
raw_data <- list()
prov_data <- list()



for (c in countries) {

  print(c)
  temp <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
 temp$hiv_treat_count <- as.numeric(gsub(",","",temp$hiv_treat_count))
 temp <- temp[which(temp$age_group_category %in% c("adults", "Adults") | is.na(temp$age_group_category)), ]
 temp$site_memo <- as.character(temp$site_memo)
 if (c == "MWI") {
   temp$hiv_treat_count <- as.numeric(temp$adult_alive)
   temp <- temp[which(temp$location_code != 157), ]
 }

 if (c == "MOZ") {
   temp <- temp[which(temp$end_year > 2008), ]
 }

 if (c == "RWA") {
   temp$age_group_category <- as.character(temp$age_group_category)
   temp$age_group_category <- "adults"
 }

 if (c == "BFA") {
   temp$age_group_category <- as.character(temp$age_group_category)
   temp$age_group_category[which(temp$age_group_category == "Adults")] <- "adults"
   temp <- temp[which(temp$end_year < 2015), ]
 }

 if (c == "BEN") {
   temp$age_group_category <- as.character(temp$age_group_category)
   temp$age_group_category[which(is.na(temp$age_group_category))] <- "adults"
   temp <- temp[which(temp$end_year > 2016), ]
 }

 if (c == "BDI") {
   temp <- temp[which(temp$end_year < 2017), ]
   temp <- temp[!which(temp$ADM2_NAME == "Ntega" & temp$end_year == 2015), ]
 }

 if (c == "CIV") {
   temp <- temp[!which(temp$end_year == 2016 & temp$hiv_treat_source == "PEPFAR"), ]
 }

 if (c == "CMR") {
   temp <- temp[which(temp$end_year %in% c(2006, 2007, 2008, 2012)), ]
 }

 if (c == "ZAF") {
   temp$age_group_category[which(is.na(temp$age_group_category))] <- "adults"
   temp <- temp[!which(temp$hiv_treat_source == "PEPFAR"), ]
 }
 if (c == "NGA") {
   ondo <- grep("Ondo", temp$site_memo, value = T)
   benue <- grep("Benue", temp$site_memo, value = T)
   kaduna <- grep("Kaduna", temp$site_memo, value = T)
   plateau <- grep("Plateau", temp$site_memo, value = T)
   sites <- c(ondo, benue, kaduna, plateau)
   temp <- temp[which(temp$site_memo %in% sites), ]

 }
 if (c == "ZWE") {
   bulawayo <- grep("Bulawayo", temp$site_memo, value = T)
   harare <- grep("Harare", temp$site_memo, value = T)
   sites <- c(bulawayo, harare)
   temp <- temp[which(temp$site_memo %in% sites), ]
 }

 temp$location_code <- as.character(temp$location_code)
 temp <- temp[which(!is.na(temp$hiv_treat_count)), ]
 temp$end_year <- as.numeric(temp$end_year)


 if ("hiv_treat_source" %!in% names(temp)) { temp$hiv_treat_source <- NA }
 raw_data[[c]] <- temp[ ,c("hiv_treat_count","end_year", "age_group_category", "location_code", "hiv_treat_source")]





}

##########
## province data from unaids for use later
##########


fs <- list.files(paste0("<<<< FILEPATH REDACTED >>>>"), full.names = FALSE, recursive = FALSE)
fs <- grep(wp_version, fs, value = TRUE)
if (length(fs) == 1) {
  prov_data <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  } else {
    fs <- fs[which(substr(fs, 1,2) == "20")]
    dates <- list()
    for (id in c(1:length(fs))) {
      d <- as.Date(str_replace_all(substr(fs[[id]], 1, 10), "-", "/"))
      l2 <- list()
      l2[["date"]] <- d
      l2[["file"]] <- fs[id]
      dates[[id]] <- l2
    }
    dates <- bind_rows(dates)
    prov_data <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  }

###########
## combine extracted data into one data frame
###########

data_s <- bind_rows(raw_data, .id = "ISO")

data_s <- data_s[which(data_s$end_year < 2019), ]


##_read in geographic hierarchy
ad2 <- st_read(paste0("<<<< FILEPATH REDACTED >>>>"))
ad2$geometry <- NULL

ad2 <- ad2[ ,c("ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE", "ADM2_NAME", "ADM2_CODE")]
ad2$ADM2_CODE <- as.character(ad2$ADM2_CODE)

##_compress data to one sex
data_s <- data_s[,lapply("hiv_treat_count", function(x) sum(get(x), na.rm = TRUE)), by = c("end_year", "age_group_category", "location_code", "ISO")]


##_add the geographic hierarchy
data_s <- merge(data_s, ad2, by.x = c("location_code"), by.y = c("ADM2_CODE"), all.x = TRUE)

## read in the UNAIDS / province look up tables.

prov_lookup <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

##_calculate ADM1 and ADM0 totals
data_s$hiv_treat_count <- data_s$V1


##_pull in populations
raw_pop <- list()

for (c in countries) {
  temp <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  raw_pop[[c]] <- temp
}

##_combine pops into one data frame
pops <- as.data.table(bind_rows(raw_pop, .id = "ISO"))

##_organize things a bit
pops <- pops[,lapply(c("population"), function(x) sum(get(x))), by = c("year_id", "ADM2_CODE")]
pops$ADM2_CODE <- as.character(pops$ADM2_CODE)
pops$pop <- pops$V1

##read in prvalence data
prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

plhiv <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
plhiv_0 <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
plhiv_0$plhiv_0 <- plhiv_0$mean
plhiv_0 <- plhiv_0[,c("ADM0_NAME", "year", "plhiv_0")]
plhiv_0$ADM0_NAME[which(plhiv_0$ADM0_NAME == "Swaziland")] <- "Eswatini"


prev$ADM2_CODE <- as.character(prev$ADM2_CODE)
plhiv$ADM2_CODE <- as.character(plhiv$ADM2_CODE)

prev$prev <- prev$mean
plhiv$plhiv <- plhiv$mean

prev <- merge(prev[,c("ADM2_CODE", "ADM0_CODE", "year", "prev")], plhiv[ ,c("ADM2_CODE", "year", "plhiv")], by = c("ADM2_CODE", "year"))


########
## Make Predictions Frame
########
pred_frame <- merge(pops, prev, by.x = c("year_id", "ADM2_CODE"), by.y = c("year", "ADM2_CODE"), all.x = TRUE)



pred_frame$location_code <- pred_frame$ADM2_CODE
pred_frame$end_year <- pred_frame$year_id
pred_frame <- pred_frame[which(pred_frame$year_id > 1994), ]

########
##_find when ART started in each country
#######
unaids_art <- list()
art_start <- list()

for (c in countries) {
  f <- paste0("<<<< FILEPATH REDACTED >>>>")
  if (file.exists(f)) {
  temp <- fread(f)
  temp <- temp[,lapply(c("ART_cov_num"), function(x) sum(get(x))), by = c("year")]
  yrs <- tapply(temp$V1, temp$year, function(a)head(which(a > 0),1))
  yrs <- yrs[which(yrs > 0)]
  art_start[[c]] <- (as.numeric(names(yrs[1])))
  unaids_art[[c]] <- temp
  f <- NULL
  }
}

unaids_data <- bind_rows(unaids_art, .id = "ISO")
unaids_data$unaids_nat_art <- unaids_data$V1
unaids_data$V1 <- NULL

###########
## Add dummy data from when we know all adm2 units had 0 art patients
###########
data_s <- data_s[which(data_s$hiv_treat_count != 0), ]
data_s$old_data <- 0
for (c in countries) {

  old_data <- as.data.frame(data_s[which(data_s$ISO == c), ])
  for (id in unique(old_data$location_code)) {
    if (is.na(id)) { next }
    od <- old_data[which(old_data$location_code == id), ]
    min_yr <- min(od$end_year)
    od <- od[which(od$end_year == min_yr),]
    zero_year <- art_start[[c]] - 1
    od$end_year <- zero_year
    od$V1 <- 0
    od$hiv_treat_count <- 0
    od$old_data <- 1
    data_s <- rbind(data_s, od)
    od$end_year <- zero_year - 1
    od$old_data <- 2
    data_s <- rbind(data_s, od)
    od$end_year <- zero_year - 2
    od$old_data <- 3
    data_s <- rbind(data_s, od)
  }
}


########
## ADD pop and Prev to data
########

DF <- merge(data_s, pops, by.x = c("location_code", "end_year"), by.y = c("ADM2_CODE", "year_id"), all.x = TRUE)
DF <- merge(DF, prev, by.x = c("location_code", "end_year"), by.y = c("ADM2_CODE", "year"), all.x = TRUE)

DF <- DF[!which(is.na(DF$ADM1_NAME)), ]
DF$art_cov <- DF$hiv_treat_count / DF$plhiv

problems <- pred_frame[which(pred_frame$ADM2_CODE %!in% DF$location_code), ]
other_problems <- DF[!which(DF$location_code %in% pred_frame$ADM2_CODE), ]
problems <- unique(problems)
problems <- merge(problems, ad2, by = c("ADM2_CODE"))



c_lookup$ADM0_NAME <- as.character(c_lookup$ADM0_NAME)


pred_frame <- merge(pred_frame, ad2, by = c("ADM2_CODE"))
pred_frame$ADM0_NAME <- as.character(pred_frame$ADM0_NAME)
pred_frame$ADM0_NAME[which(pred_frame$ADM0_NAME == "Swaziland")] <- "Eswatini"
pred_frame <- merge(pred_frame, c_lookup, by = c("ADM0_NAME"))
pred_frame$end_year_scaled <- pred_frame$end_year - 1994
pred_frame <- merge(pred_frame, unaids_data, by.x = c("loc", "year_id"), by.y = c("ISO", "year"))
pred_frame <- merge(pred_frame, plhiv_0, by.x = c("ADM0_NAME", "year_id"), by.y = c("ADM0_NAME", "year"))
pred_frame$nat_cov <- pred_frame$unaids_nat_art / pred_frame$plhiv_0

pred_frame <- merge(pred_frame, prov_lookup, by = c("ADM1_CODE", "ADM1_NAME"), all.x = TRUE)
pred_frame <- merge(pred_frame, prov_data, by.x = c("full_name", "year_id"), by.y = c("full_name", "year"), all = T)

pred_frame$target_cov <- NA
pred_frame$target_cov[which(is.na(pred_frame$full_name))] <- pred_frame$nat_cov[which(is.na(pred_frame$full_name))]
pred_frame$target_cov[which(!is.na(pred_frame$full_name))] <- pred_frame$prov_cov[which(!is.na(pred_frame$full_name))]

### subbing in 2018 for locations missing 2019 coverage
pl <- pred_frame$full_name[which(is.na(pred_frame$target_cov))]

pred_frame$target_cov[which(is.na(pred_frame$target_cov) & pred_frame$full_name %in% pl & pred_frame$year_id == 2019)] <- pred_frame$target_cov[which(pred_frame$full_name %in% pl & pred_frame$year_id == 2018)]


DF <- merge(DF, unaids_data, by.x = c("ISO", "end_year"), by.y = c("ISO", "year"))
DF$ADM0_NAME <- as.character(DF$ADM0_NAME)
DF$ADM0_NAME[which(DF$ADM0_NAME == "Swaziland")] <- "Eswatini"
DF <- merge(DF, plhiv_0, by.x = c("ADM0_NAME", "end_year"), by.y = c("ADM0_NAME", "year"))


DF <- as.data.frame(DF)
pred_frame <- as.data.frame(pred_frame)



data_small_list <- list()
pred_small_list <- list()

shp <- st_read(paste0("<<<< FILEPATH REDACTED >>>>"))
shp$ADM0_NAME <- as.character(shp$ADM0_NAME)
shp$ADM0_NAME[which(shp$ADM0_NAME == "Swaziland")] <- "Eswatini"
inla.setOption(smtp = 'taucs')
inla.setOption(pardiso.license = paste0("<<<< FILEPATH REDACTED >>>>"))
#################################################################################################################################################################
c_problems <- countries[which(countries %!in% unique(DF$ISO))]

for (c in (countries)) {

message(paste0(c))


name <- unique(DF$ADM0_NAME[which(DF$ISO == c)])

c_shp <- shp[which(shp$ADM0_NAME == name), ]
c_shp$ID_nb <- seq.int(nrow(c_shp))
c_shp <- as_Spatial(c_shp)
c_nb <- poly2nb(c_shp)
nb2INLA("c_graph",c_nb)
c_nb_INLA <- inla.read.graph(filename = "c_graph")

c_nb_lookup <- c_shp@data[,c("ADM2_CODE", "ID_nb")]
c_nb_lookup$ADM2_CODE <- as.character(c_nb_lookup$ADM2_CODE)

data_small <- DF[which(DF$ISO == c), ]
data_small <- merge(data_small, c_nb_lookup, by.x = "location_code", by.y = c("ADM2_CODE"))
outcomes <- data_small[,c("location_code", "end_year", "art_cov", "hiv_treat_count")]

pred_frame_small <- pred_frame[which(pred_frame$loc == c), ]
pred_frame_small <- merge(pred_frame_small, c_nb_lookup, by.x = "location_code", by.y = c("ADM2_CODE"), all = T)
pred_frame_small <- merge(pred_frame_small, outcomes, by = c("location_code", "end_year"), all.x = TRUE)
pred_frame_small$art_cov[which(pred_frame_small$t_var == 0)] <- NA


data_yrs <- unique(data_small$t_var)


  message(paste0(c, " mod Gaus BYM on nat_cov"))

  f10 <- art_cov ~ -1 + target_cov + f(ID_nb, target_cov,
                                       model="bym2",
                                       scale.model = T,
                                       graph = c_nb_INLA,
                                       hyper = list(
                                         phi = list(prior = "pc",
                                                    param = c(0.5, 0.5),
                                                    initial = -3),
                                         prec = list(prior = "pc.prec",
                                                     param = c(1, 0.01),
                                                     initial = 5)))

  mod4_gaus_INLA <- inla(formula = f10,
                         data = pred_frame_small,
                         family = "gaussian",
                         control.compute = list(dic = F, openmp.strategy = "default",  smtp = 'taucs'),
                         num.threads = 7,
                         verbose = F,
                         control.predictor = list(link = 1))


  if (!is.numeric(mod4_gaus_INLA)) {
    pred_frame_small <- cbind(pred_frame_small, mod4_gaus_INLA$summary.fitted.values$mean)
    pred_frame_small[which(pred_frame_small$target_cov == 0),c("mod4_gaus_INLA$summary.fitted.values$mean")] <- 0 # correcting for a strange INLA thing where it puts in a small number instead of 0.
    pred_frame_small$pred4_gaus_INLA <- pred_frame_small[,c("mod4_gaus_INLA$summary.fitted.values$mean")] * pred_frame_small$plhiv
    pred_frame_small[,c("mod4_gaus_INLA$summary.fitted.values$mean")] <- NULL
    saveRDS(mod4_gaus_INLA, paste0("<<<< FILEPATH REDACTED >>>>"))
  } else {pred_frame_small$pred4_gaus_INLA <- -1 }


pred_small_list[[c]] <- pred_frame_small

pred_frame_small <- NULL
c_shp <- NULL
c_nb <- NULL
c_nb_INLA <- NULL
c_nb_lookup <- NULL
data_small <- NULL
outcomes <- NULL

}

data_final <- bind_rows(pred_small_list, .id = "ISO2")


###### removing admin 2 shapes that are just lakes, this may be temporary
lakes <- grep("Lake", unique(data_final$ADM2_NAME), value = T)
data_final <- data_final[which(data_final$ADM2_NAME %!in% lakes), ]
######


preds <- grep("pred", names(data_final), value = T)

for (var in preds) {
  data_final[which(data_final[ , c(var)] < 0), c(var)] <- 0
}
data_final <- as.data.table(data_final)
data_final$full_name[which(is.na(data_final$full_name))] <- data_final$ADM0_NAME[which(is.na(data_final$full_name))]

scalars <- data_final[,lapply(preds, function(x) sum(get(x))), by = c("ADM0_NAME", "year_id")]
names(scalars) <- c("ADM0_NAME", "year_id", paste(preds, "_tot", sep = ''))

data_final <- as.data.frame(merge(data_final, scalars, by = c("ADM0_NAME", "year_id")))



for (var in preds) {
  var_tot <- paste0(var, "_tot")
  var_scale <- paste0(var, "_scale")
  var_s <- paste0(var, "_s")
  data_final[, c(var_scale)] <- data_final$unaids_nat_art / data_final[ ,c(var_tot)]
  data_final[, c(var_scale)][which(is.na(data_final[, c(var_scale)]))] <- 0
  data_final[, c(var_scale)][which(data_final[, c(var_scale)] == Inf)] <- 1
  data_final[, c(var_s)] <- data_final[,c(var)] * data_final[, c(var_scale)]
}


codes <- as.data.table(unique(data_final[, c("ADM0_NAME","ADM1_NAME", "ADM2_NAME", "ADM2_CODE")]))
setorder(codes, "ADM0_NAME", "ADM1_NAME", "ADM2_NAME")
codes <- codes[which(!is.na(codes$ADM0_NAME)), ]
dir <- paste0("<<<< FILEPATH REDACTED >>>>")
dir.create(dir)

pdf(paste0(dir,"art_models_comp_new.pdf"), width = 11, height = 8)
for (id in 1:length(codes$ADM2_CODE)) {
  df2 <- data_final[which(data_final$location_code == codes$ADM2_CODE[[id]]), ]
  gg_count <- ggplot() +
    geom_line(data = df2, aes(x = year_id, y = (pred4_gaus_INLA ), color = "Gausian Model (bym nat_cov)", linetype = "un-scaled")) +
    geom_line(data = df2, aes(x = year_id, y = (pred4_gaus_INLA_s ), color = "Gausian Model (bym nat_cov)", linetype = "scaled")) +
    geom_point(data = df2, aes(x = year_id, y = hiv_treat_count, color = "data"), size = 3) +
    scale_color_manual(values = c("#000000","#FF0000","#FFA500", "#0000FF", "#008000")) +
    xlab("year") +
    ylab("People on ART") +
    labs(title = paste0("# ART in ad2 ", codes$ADM2_NAME[[id]], " in ad1 ", codes$ADM1_NAME[[id]], " in ", codes$ADM0_NAME[[id]])) +
    theme(legend.position = "none", plot.title = element_text(size = 12)) +
    coord_cartesian(xlim = c(1995, 2018), ylim = c(0, max(df2$pred4_gaus_INLA)))


  gg_cov <- ggplot() +
    geom_line(data = df2, aes(x = year_id, y = (pred4_gaus_INLA / plhiv), color = "Gausian Model (bym nat_cov)", linetype = "un-scaled")) +
    geom_line(data = df2, aes(x = year_id, y = (pred4_gaus_INLA_s / plhiv), color = "Gausian Model (bym nat_cov)", linetype = "scaled")) +
    geom_point(data = df2, aes(x = year_id, y = art_cov, color = "data"), size = 3) +
    scale_color_manual(values = c("#000000","#FF0000","#FFA500", "#0000FF", "#008000")) +
    xlab("year") +
    ylab("ART Coverage") +
    labs(title = "ART coverage") +
    theme(legend.position = c(-0.5,-0.20), legend.justification = c(0, 0), legend.direction = "horizontal", plot.title = element_text(size = 12)) +
    coord_cartesian(xlim = c(1995, 2018), ylim = c(0,1))

  grid.newpage()
  vp_1 <- viewport(width = 0.5, height = 0.85, x = 0.25, y = 0.575)
  vp_2 <- viewport(width = 0.5, height = 0.85, x = 0.75, y = 0.575)
  print(gg_count, vp = vp_1)
  print(gg_cov, vp = vp_2)
  print(id)
}
dev.off()



shp <- raster::shapefile(paste0("<<<< FILEPATH REDACTED >>>>"))
shp <- data.table(fortify(shp, region = "ADM2_CODE"))
shp$ADM0_NAME <- NULL
shp <- merge(ad2, shp, by.x = c("ADM2_CODE"), by.y = c("id"))
shp$ADM0_NAME <- as.character(shp$ADM0_NAME)
shp$ADM0_NAME[which(shp$ADM0_NAME == "Swaziland")] <- "Eswatini"

colors <- c('#ffffe0','#ffe4ac','#ffc879','#ffa84c','#ff8725','#ff5c03','#f12861','#cc117d','#a60383','#800080')
color_values <- scales::rescale(unique(c(seq(0, 5, length.out = 4), seq(5, 10, length.out = 4), seq(10, 25, length.out = 4))))


for (c in as.character(c_lookup$ADM0_NAME)) {
  sub_shp <- shp[which(shp$ADM0_NAME == c), ]
      pdf(file = paste0(dir,"art_maps_",c,"_nat_cov.pdf"), width = 11, height = 8)
      for (y in c(1995:2018)) {
        message(paste0(c," ",y))
        data <- data_final[which(data_final$ADM0_NAME == c & data_final$year_id == y), ]
        data <- merge(data, sub_shp, by = c("ADM2_CODE"), all.y = TRUE)
        data <- data[order(data$order),]
        gaps <- data[which(data$ADM2_CODE %in% problems$ADM2_CODE), ]
        gg_map <- ggplot() +
          geom_polygon(data = data, aes(fill = (pred4_gaus_INLA_s / plhiv), y = lat, x = long, group = group)) +
          coord_equal() +
          scale_fill_gradientn(colors = colors, values = color_values,
                               limits = c(0, 1), breaks = seq(0, 1, 0.1), labels = seq(0, 1, 0.1),
                               name = paste0('Modeled ART Coverage ', y)) +
          labs( x = NULL, y = NULL, title = NULL) +
          theme_classic(base_size = 10) +
          theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                legend.position = c(-0.1, 0.1), legend.justification = c(0, 0),
                legend.title = element_text(angle = 90, hjust = 0.5),
                plot.title = element_text(hjust = 0.5), plot.margin = unit(c(-.13, -.13, -.13, -.13), "in"),
                panel.background = element_blank()) +
          guides(fill = guide_colorbar(barwidth = 0.75, barheight = 7, nbin = 1000, title.position = 'left'))

        print(gg_map)
      }
      dev.off()
}


## Export model for use in EPP

data_final$prev_prop <- data_final$pred4_gaus_INLA_s / data_final$unaids_nat_art

export_data <- data_final[ ,c("end_year", "location_code", "prev_prop")]
names(export_data) <- c("year", "ADM2_CODE", "prev_prop")
write.csv(export_data, paste0("<<<< FILEPATH REDACTED >>>>"))







