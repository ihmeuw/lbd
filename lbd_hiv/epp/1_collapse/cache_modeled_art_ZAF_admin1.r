#### CMR admin 1 script
#### need to create an admin 1 time series for ZAF



rm(list = ls())

#########
## Set Up
#########
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf('mbg/common_inputs')
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



#####
## Determine what countries we have ART data for
#####

countries <- c("ZAF")


#######
## set the HIV prevalence extraction to use for modeling ART coverage.
#######
prev_date <- "2020_04_10" #string that indicates what prevalence is used to model coverage
wp_version <- "wp1" #string indicating what wp version the prevalence model used.

prev_info <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
mbg_run <- as.character(prev_info$Value[which(prev_info$Setting == "mbg_run")])

shp_date <- as.character(prev_info$Value[which(prev_info$Setting == "shapefile_version")])


#####
## Identify the appropriate file to read in for each country, in this case the newest dated file
#####
files <- list()
for (loc in countries) {
  fs <- list.files(paste0("<<<< FILEPATH REDACTED >>>>"), full.names = FALSE, recursive = FALSE)
  fs <- grep("thembisa", fs, value = TRUE)
  if (length(fs) == 1) { files[[loc]] <- fs[[1]]}
  else {
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

c_lookup$ADM0_NAME[which(c_lookup$ADM0_NAME == "United Republic of Tanzania")] <- "Tanzania" # the appropriate name for Tanzania is the short one.




########
## setting loop for countires to collect the data
########
DFs <- list()
raw_data <- list()
prov_data <- list()



for (c in countries) {

  temp <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  temp$hiv_treat_count <- as.numeric(temp$alive_ART)
  temp$location_code <- as.character(temp$prov_code)
  temp <- temp[which(temp$prov_code != "SA"), ] #removing data that is not actually at the admin1 level
  temp$hiv_treat_count[which(is.na(temp$hiv_treat_count))] <- 0
  temp$end_year <- as.numeric(temp$year)
  temp$age_group_category <- "adults"
  if ("hiv_treat_source" %!in% names(temp)) { temp$hiv_treat_source <- NA }
  raw_data[[c]] <- temp[ ,c("hiv_treat_count","end_year", "age_group_category", "location_code", "hiv_treat_source")]

}


###########
## combine extracted data into one data frame
###########

data_s <- bind_rows(raw_data, .id = "ISO")

data_s <- data_s[which(data_s$end_year < 2019), ] # 2019 data is not trusted at the moment because the year has not finished yet.


##_read in geographic hierarchy
ad1 <- st_read(paste0("/<<<< FILEPATH REDACTED >>>>"))
ad1$geometry <- NULL
ad1 <- ad1[ ,c("ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE")]
ad1$ADM1_CODE <- as.character(ad1$ADM1_CODE)

##_compress data to one sex
data_s <- data_s[,lapply("hiv_treat_count", function(x) sum(get(x), na.rm = TRUE)), by = c("end_year", "age_group_category", "location_code", "ISO")]


##_add the geographic hierarchy
lookup <- fread("<<<< FILEPATH REDACTED >>>>")


data_s <- merge(data_s, lookup, by.x = c("location_code"), by.y = c("UNAIDS_NAME"), all.x = TRUE)

##_calculate ADM1 and ADM0 totals
data_s$hiv_treat_count <- data_s$V1


nat_tot <- data_s[,lapply("hiv_treat_count", function(x) sum(get(x), na.rm = TRUE)), by = c("end_year", "age_group_category")]
nat_tot$nat_count <- nat_tot$V1
nat_tot$V1 <- NULL
data_s$V1 <- NULL


##_read in plhiv data

plhiv <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
plhiv_0 <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
plhiv_0$plhiv_0 <- plhiv_0$mean
plhiv_0 <- as.data.frame(plhiv_0[,c("ADM0_NAME", "year", "plhiv_0")])

plhiv$ADM1_CODE <- as.character(plhiv$ADM1_CODE)


plhiv$plhiv <- plhiv$mean

plhiv <- plhiv[ ,c("ADM1_CODE", "ADM0_CODE", "year", "plhiv")]

data_s <- merge(data_s, nat_tot, by = c("end_year"))


########
##_find when ART started in each country
#######
unaids_art <- list()
art_start <- list()

for (c in countries) {
  f <- paste0("/<<<< FILEPATH REDACTED >>>>")
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
## Add dummy data from when we know all adm2 units had 0 art patients (because UNAIDS reports 0 ART patients nationally, this is probably not totally true but it will peg the start of the roll out and improve the time trend)
###########

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
    od$hiv_treat_count <- 0
    od$nat_count <- 0
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
## Make Predictions Frame
########
pred_frame <- data_s

years_to_add <- c(1995:2018)[which(c(1995:2018) %!in% data_s$end_year)]

for (yr in years_to_add) {
  data_1 <- data_s[which(data_s$end_year == 2009), ]
  data_1$end_year <- yr
  data_1$hiv_treat_count <- NA
  data_1$nat_count <- NA
  pred_frame <- rbind(pred_frame, data_1)

}

c_lookup$ADM0_NAME <- as.character(c_lookup$ADM0_NAME)

plhiv_0 <- plhiv_0[which(plhiv_0$ADM0_NAME == "South Africa"), ]
pred_frame <- merge(pred_frame, unaids_data, by.x = c("ISO", "end_year"), by.y = c("ISO", "year"))
pred_frame <- merge(pred_frame, plhiv_0, by.x = c("end_year"), by.y = c("year"))
pred_frame$nat_cov <- pred_frame$unaids_nat_art / pred_frame$plhiv_0



pred_frame <- as.data.frame(pred_frame)


pred_small_list <- list()

shp <- st_read(paste0("<<<< FILEPATH REDACTED >>>>"))

#################################################################################################################################################################

for (c in (countries)) {

  message(paste0(c))


  name <- unique(pred_frame$ADM0_NAME[which(pred_frame$ISO == c)])

  c_shp <- shp[which(shp$ADM0_NAME == name), ]
  c_shp$ID_nb <- seq.int(nrow(c_shp))
  c_shp <- as_Spatial(c_shp)
  c_nb <- poly2nb(c_shp)
  nb2INLA("c_graph",c_nb)
  c_nb_INLA <- inla.read.graph(filename = "c_graph")

  c_nb_lookup <- c_shp@data[,c("ADM1_CODE", "ID_nb")]
  c_nb_lookup$ADM1_CODE <- as.character(c_nb_lookup$ADM1_CODE)


  pred_frame_small <- pred_frame
  pred_frame_small <- merge(pred_frame_small, c_nb_lookup, by.x = c("ADM1_CODE"), by.y = c("ADM1_CODE"), all = T)
  pred_frame_small <- merge(pred_frame_small, plhiv, by.x = c("ADM1_CODE", "end_year"), by.y = c("ADM1_CODE", "year"), all.x = TRUE)
  pred_frame_small$art_cov <- pred_frame_small$hiv_treat_count / pred_frame_small$plhiv
  n.block <- max(pred_frame_small$ID_nb)
  pred_frame_small$i.intercept <- pred_frame_small$ID_nb
  pred_frame_small$j.intercept <- pred_frame_small$ID_nb + n.block  ## see doc for iid2d

  message(paste0(c, " mod Gaus BYM on nat_cov"))

  f10 <- art_cov ~ -1 + nat_cov + f(ID_nb, nat_cov, model="bym2", scale.model = T, graph = c_nb_INLA)

  mod4_gaus_INLA <- inla(formula = f10,
                         data = pred_frame_small,
                         family = "gaussian",
                         #control.compute = list(dic = T),
                         num.threads = 7,
                         verbose = F,
                         control.predictor = list(link = 1))


  if (!is.numeric(mod4_gaus_INLA)) { pred_frame_small <- cbind(pred_frame_small, mod4_gaus_INLA$summary.fitted.values$mean)
  pred_frame_small$pred4_gaus_INLA <- pred_frame_small[,c("mod4_gaus_INLA$summary.fitted.values$mean")] * pred_frame_small$plhiv
  pred_frame_small[,c("mod4_gaus_INLA$summary.fitted.values$mean")] <- NULL
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

preds <- grep("pred", names(data_final), value = T)

for (var in preds) {
  data_final[which(data_final[ , c(var)] < 0), c(var)] <- 0
}
data_final <- as.data.table(data_final)
data_final$year_id <- data_final$end_year
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


codes <- as.data.table(unique(data_final[, c("ADM1_NAME", "location_code")]))
setorder(codes,"ADM1_NAME")
dir <- paste0("<<<< FILEPATH REDACTED >>>>")
dir.create(dir)
dir <- paste0("<<<< FILEPATH REDACTED >>>>")
dir.create(dir)

pdf(paste0(dir,Sys.Date(),"_admin1_ts.pdf"), width = 11, height = 8)
for (id in 1:length(codes$location_code)) {
  df2 <- data_final[which(data_final$location_code == codes$location_code[[id]]), ]
  gg_count <- ggplot() +
    geom_line(data = df2, aes(x = year_id, y = (pred4_gaus_INLA ), color = "Gausian Model (bym nat_cov)", linetype = "un-scaled")) +
    geom_line(data = df2, aes(x = year_id, y = (pred4_gaus_INLA_s ), color = "Gausian Model (bym nat_cov)", linetype = "scaled")) +
    geom_point(data = df2, aes(x = year_id, y = hiv_treat_count, color = "data"), size = 3) +
    scale_color_manual(values = c("#000000","#FF0000","#FFA500", "#0000FF", "#008000")) +
    xlab("year") +
    ylab("People on ART") +
    labs(title = paste0("# ART in ", codes$ADM1_NAME[[id]])) +
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





## Export


export_data <- data_final[ ,c("end_year", "ADM1_CODE", "pred4_gaus_INLA_s", "plhiv")]
names(export_data) <- c("year", "ADM1_CODE", "modeled_prov_tot", "plhiv")
export_data$modeled_prov_cov <- export_data$modeled_prov_tot / export_data$plhiv
write.csv(export_data, paste0("<<<< FILEPATH REDACTED >>>>"))

