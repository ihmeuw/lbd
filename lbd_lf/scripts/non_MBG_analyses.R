#####################################################################################################################################
### Lymphatic filariasis non-MBG analyses for small geographies
#####################################################################################################################################

### Initial setup
user <- Sys.info()[["user"]] ## Get current user name

core_repo <- paste0(<<<< FILEPATH REDACTED >>>>)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- sprintf(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))
package_list <- c(package_list, "sf")

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)

path <- paste0(<<<< FILEPATH REDACTED >>>>)
library(fasterize)
library(sf)
library(cowplot)
library(matrixStats)
library(sp)

INLA:::inla.dynload.workaround() 

modeling_shapefile_version <- "2020_02_20"
use_premade <- FALSE

output_folder <- <<<< FILEPATH REDACTED >>>>

raster_agg_factor <- 1

### Load LF data
pre_resampling_lf_data <- as.data.table(read.csv(<<<< FILEPATH REDACTED >>>>, header = TRUE))
post_resampling_lf_data <- as.data.table(read.csv(<<<< FILEPATH REDACTED >>>>, header = TRUE))

pre_resampling_lf_data <- pre_resampling_lf_data[country %in% c("ASM", "BRA", "COK", "FSM", "FJI", "PYF", "GUY", "KIR", "MDV", "MHL", "NCL", "NIU", "PLW", "WSM", "TON", "TUV", "VUT", "WLF"), ]

### Subset to non-MBG countries
admin0_shp <- rgdal::readOGR(dsn = get_admin_shapefile(admin_level = 0, version = modeling_shapefile_version))
admin1_shp <- rgdal::readOGR(dsn = get_admin_shapefile(admin_level = 0, version = modeling_shapefile_version))
admin2_shp <- rgdal::readOGR(dsn = get_admin_shapefile(admin_level = 2, version = modeling_shapefile_version))

#### Get country IDs
location_hierarchy <- as.data.table(get_location_metadata(location_set_id = 35, gbd_round_id = 7, decomp_step = "iterative"))
pre_resampling_lf_data <- merge(pre_resampling_lf_data, location_hierarchy, by.x = "country", by.y = "ihme_loc_id", all.x = TRUE)
post_resampling_lf_data <- merge(post_resampling_lf_data, location_hierarchy, by.x = "country", by.y = "ihme_loc_id", all.x = TRUE)
pre_resampling_lf_data[country %in% c("NCL", "TUV", "PYF", "NRU", "NIU", "NCL", "PYF", "WLF"), location_id := 21] ### Use population distribution for Oceania, as GBD doesn't break out these countries separately
post_resampling_lf_data[country %in% c("NCL", "TUV", "PYF", "NRU", "NIU", "NCL", "PYF", "WLF"), location_id := 21] ### Use population distribution for Oceania, as GBD doesn't break out these countries separately

post_resampling_lf_data$lf_prev <- post_resampling_lf_data$had_lf_w_resamp/post_resampling_lf_data$N

region <- "lf_non_mbg"
link_table <- paste0(<<<< FILEPATH REDACTED >>>>)

# Load simple polygon template to model over
gaul_list <- get_adm0_codes(region, shapefile_version = modeling_shapefile_version)

simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = use_premade, shapefile_version = modeling_shapefile_version)
subset_shape <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]

raster_list <- build_simple_raster_pop(subset_shape)

simple_raster <- raster_list[["simple_raster"]]
pop_raster <- raster_list[["pop_raster"]]

# Get GBD populations
locs <- get_location_metadata(location_set_id = 35, gbd_round_id = 7, decomp_step = "iterative")
loc_ids_ad0 <- locs[level == 3, location_id]
gbd_pops <- get_population(location_id = locs$location_id, year_id = 1990:2018, gbd_round_id = 7, decomp_step = "iterative")
gbd_merged <- merge(gbd_pops, locs, by = "location_id", all.x = TRUE)

####################################################################################
###### Run Niue model
rm(r)
rm(data_set_preds)

## Niue data are nationally representative
data_NIU <- post_resampling_lf_data[country == "NIU"]

### Add rows to data frame for prediction
new <- data_NIU[1]
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]

new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018

data_NIU <- rbind(data_NIU, new, fill = TRUE, use.names = TRUE)

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_NIU, control.predictor = list(compute = TRUE, link = 1), control.inla = list(int.strategy = "eb"), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_NIU$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_NIU, inv.logit(r$summary.linear.predictor[1:nrow(data_NIU), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_NIU), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_NIU), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp)], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp)], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp)], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp)], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Niue (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
print(g)
dev.off()

results_table <- setorderv(data_set_preds[is.na(had_lf_w_resamp), c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))


### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_NIU$Master_UID))
draw_preds <- matrix(nrow = nrow(data_NIU[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds_2018 <- inv.logit(draw_preds)[29,]
post_prob_2 <- c(country="NIU", post_prob_2=(length(draw_preds_2018[draw_preds_2018 < 0.02])/length(draw_preds_2018)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",")

draw_preds <- cbind("year"=predict_years, as.data.table(inv.logit(draw_preds)))
write.table(draw_preds, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_NIU[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Palau model
rm(r)
rm(data_set_preds)

data_PLW <- as.data.table(read.csv(<<<< FILEPATH REDACTED >>>>, header = TRUE))
data_PLW <- data_PLW[country == "PLW"]
data_PLW <- data_PLW[Master_UID %in% c("lit_dos_1", "lit_dos_5", "lit_dos_14")]

gaul_code <- 178
subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

new <- data_PLW[1]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_poly", "lf_prev") := list(100, NA, NA)]

data_PLW <- rbind(data_PLW, new)

formula <- had_lf_poly ~ f(year, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_PLW, control.predictor = list(compute = TRUE, link = 1), control.inla = list(int.strategy = "eb"), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, verbose = TRUE)

data_set_preds <- cbind(data_PLW, inv.logit(r$summary.linear.predictor[1:nrow(data_PLW), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_PLW), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_PLW), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_poly)], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_poly)], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_poly)], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_poly)], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Ngardmau, Palau (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

### Population data from 2005 Palau census (assume that proportion by state is constant across years)
data_set_preds$population_total <- 19907
data_set_preds$population_Ngardmau <- 166

data_set_preds$national_lf_prev <- data_set_preds$lf_prev * data_set_preds$population_Ngardmau / data_set_preds$population_total
data_set_preds$national_mean <- data_set_preds$mean * data_set_preds$population_Ngardmau / data_set_preds$population_total
data_set_preds$national_quant0.025 <- data_set_preds$quant0.025 * data_set_preds$population_Ngardmau / data_set_preds$population_total
data_set_preds$national_quant0.975 <- data_set_preds$quant0.975 * data_set_preds$population_Ngardmau / data_set_preds$population_total
g <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_poly)], aes(x = year, y = national_mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_poly)], aes(x = year, ymin = national_quant0.025, ymax = national_quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_poly)], aes(x = year, y = national_mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_poly)], aes(x = year, y = national_lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Palau (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
print(g)
dev.off()

results_table <- setorderv(data_set_preds[is.na(had_lf_poly), c("country", "year", "national_mean", "national_quant0.025", "national_quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))


### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_PLW$Master_UID))
draw_preds <- matrix(nrow = nrow(data_PLW[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds) * 166 / 19907

draw_preds_2018 <- draw_preds[29,]
post_prob_2 <- c(country="PLW", post_prob_2=(length(draw_preds_2018[draw_preds_2018 < 0.02])/length(draw_preds_2018)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

draw_preds <- cbind("year"=predict_years, as.data.table(draw_preds))
write.table(draw_preds, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_PLW[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Vanuatu model
rm(r)
rm(data_set_preds)

data_VUT <- post_resampling_lf_data[country == "VUT"]

gaul_code <- 243
subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

resampled_VUT <- data_VUT

admin2_shp <- rgdal::readOGR(dsn = get_admin_shapefile(admin_level = 2, version = modeling_shapefile_version))
admin2_shp@data$ADM1_CODE <- as.integer(levels(admin2_shp@data$ADM1_CODE))[admin2_shp@data$ADM1_CODE]
admin2_shp@data$ADM2_CODE <- as.integer(levels(admin2_shp@data$ADM2_CODE))[admin2_shp@data$ADM2_CODE]

vanuatu_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "Vanuatu", ], "SpatialPolygonsDataFrame"))
vanuatu_adm1 <- fasterize(vanuatu_df_sf, simple_raster2, field = "ADM1_CODE")
vanuatu_adm2 <- fasterize(vanuatu_df_sf, simple_raster2, field = "ADM2_CODE")

o <- over(SpatialPointsDataFrame(cbind(resampled_VUT$longitude, resampled_VUT$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_VUT), admin2_shp)
resampled_VUT$ADM1_NAME <- o$ADM1_NAME
resampled_VUT$ADM2_NAME <- o$ADM2_NAME

data_VUT_merged <- resampled_VUT

data_VUT_merged[shapefile == "lf_gadm28_adm1" & location_code == 3466, ADM1_NAME := "Shefa"]
data_VUT_merged[shapefile == "lf_gadm28_adm1" & location_code == 3465, ADM1_NAME := "Sanma"]
data_VUT_merged[shapefile == "lf_gadm28_adm1" & location_code == 3463, ADM1_NAME := "Malampa"]
data_VUT_merged[shapefile == "lf_gadm28_adm1" & location_code == 3467, ADM1_NAME := "Tafea"]
data_VUT_merged[shapefile == "lf_gadm28_adm1" & location_code == 3464, ADM1_NAME := "Penama"]
data_VUT_merged[shapefile == "lf_gadm28_adm1" & location_code == 3468, ADM1_NAME := "Torba"]
data_VUT_merged[shapefile == "lf_g2015_2009_1" & location_code == 3297, ADM1_NAME := "Penama"]

data_VUT_merged <- data_VUT_merged[!is.na(ADM1_NAME)] # Final geo-referenced data set for modeling

### INLA model for Vanuatu with adm1-level trends
data_VUT_merged$ADM1_NAME <- droplevels(data_VUT_merged$ADM1_NAME)
data_VUT_merged[ADM2_NAME %in% c("North Ambrym", "West Ambrym"), ADM1_NAME := "Ambrym Island (Malampa)"] # Treat Ambrym as a separate ADM1 due to different programmatic history

adm1 <- unique(data_VUT_merged$ADM1_NAME)
new <- data_VUT_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(adm1)), ]
new$ADM1_NAME <- rep(adm1, nrow(new)/length(adm1))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_VUT_merged <- rbind(data_VUT_merged, new)
data_VUT_merged$ADM1_NAME_int <- as.integer(data_VUT_merged$ADM1_NAME)
data_VUT_merged$year_int <- data_VUT_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01)))) + f(ADM1_NAME, model = "iid")
r <- inla(formula, family = "binomial", data = data_VUT_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_VUT_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_VUT_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_VUT_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_VUT_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_VUT_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Penama"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Penama"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Penama"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Penama"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Penama Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tafea"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tafea"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tafea"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Tafea"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tafea Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Malampa"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Malampa"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Malampa"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Malampa"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Malampa Province (sans Ambrym Island), Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Ambrym Island (Malampa)"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Ambrym Island (Malampa)"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Ambrym Island (Malampa)"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Ambrym Island (Malampa)"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Ambrym Island, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Torba"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Torba"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Torba"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Torba"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Torba Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g6 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Sanma"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Sanma"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Sanma"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Sanma"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Sanma Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g7 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Shefa"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Shefa"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Shefa"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Shefa"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Shefa Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, g6, g7, ncol = 3)


### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
preds$ADM2_NAME <- NA
predict_years <- 1990:2018

pred_indices <- which(is.na(data_VUT_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_VUT_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

for (a in predict_years) {
  print(paste0("Starting year ", a))
  pop_raster <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
  pop_raster_masked <- crop(pop_raster2, extent(subset_shape2))
  pop_raster_masked <- setExtent(pop_raster_masked, subset_shape2)
  pop_raster_masked <- raster::mask(pop_raster_masked, subset_shape2)

  for (c in unique(admin2_shp@data[admin2_shp@data$ADM0_NAME == "Vanuatu", "ADM1_NAME"])) {
    print(c)
    temp_adm1_code <- unique(admin2_shp@data[admin2_shp@data$ADM1_NAME == c, "ADM1_CODE"])
    temp_raster <- vanuatu_adm1
    temp_raster[temp_raster != temp_adm1_code] <- NA
    pop_raster_masked2 <- resample(pop_raster_masked, temp_raster)
    temp_raster2 <- setExtent(temp_raster, pop_raster_masked2)
    adm1_population <- raster::mask(pop_raster_masked2, temp_raster2)
    adm1_population <- sum(adm1_population@data@values, na.rm = TRUE)
    preds[year == a & ADM1_NAME == c, "population" := adm1_population]
  }

  print("Ambrym Island")
  temp_adm2_code <- unique(admin2_shp@data[admin2_shp@data$ADM0_NAME == "Vanuatu" & admin2_shp@data$ADM2_NAME %in% c("West Ambrym", "North Ambrym", "South East Ambrym"), "ADM2_CODE"])
  temp_raster <- vanuatu_adm2
  temp_raster[!(temp_raster %in% temp_adm2_code)] <- NA
  pop_raster_masked <- resample(pop_raster_masked, temp_raster)
  temp_raster2 <- setExtent(temp_raster, pop_raster_masked)
  adm1_population <- raster::mask(pop_raster_masked, temp_raster2)
  adm1_population <- sum(adm1_population@data@values, na.rm = TRUE)
  preds[year == a & ADM1_NAME == "Ambrym Island (Malampa)", "population" := adm1_population]
  preds[year == a & ADM1_NAME == "Malampa", "population"] <- preds[year == a & ADM1_NAME == "Malampa", "population"] - preds[year == a & ADM1_NAME == "Ambrym Island (Malampa)", "population"]
}

# National aggregates
draw_preds_pops <- draw_preds * preds$population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "Vanuatu National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- aggregate(population ~ year, data = preds, sum)$population
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

# Adm1 summaries
preds$mean <- rowMeans(draw_preds_pops) / preds$population
preds$quant0.025 <- rowQuantiles(draw_preds_pops, probs = 0.025) / preds$population
preds$quant0.975 <- rowQuantiles(draw_preds_pops, probs = 0.975) / preds$population

preds <- rbind(preds, national)

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = preds[ADM1_NAME != "Vanuatu National"], aes(x = year, y = mean, group = ADM1_NAME, color = ADM1_NAME), alpha = 0.3) + geom_line(data = preds[ADM1_NAME == "Vanuatu National"], aes(x = year, y = mean, group = ADM1_NAME), alpha = 1, size = 0.5, color = 1) + geom_ribbon(data = preds[ADM1_NAME == "Vanuatu National"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence: Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Penama"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Penama"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Penama"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Penama"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Penama Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tafea"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tafea"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tafea"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Tafea"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Tafea Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Malampa"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Malampa"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Malampa"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Malampa"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Malampa Province (sans Ambrym Island), Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Ambrym Island (Malampa)"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Ambrym Island (Malampa)"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Ambrym Island (Malampa)"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Ambrym Island (Malampa)"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Ambrym Island, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Torba"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Torba"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Torba"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Torba"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Torba Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g6 <- ggplot() + theme_classic() + geom_point(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Sanma"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Sanma"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Sanma"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Sanma"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Sanma Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g7 <- ggplot() + theme_classic() + geom_point(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Shefa"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Shefa"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Shefa"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Shefa"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Shefa Province, Vanuatu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
print(plot_grid(g1, g2, g3, g4, g5, g6, g7, ncol = 3))
dev.off()

results_table <- setorderv(preds[ADM1_NAME == "Vanuatu National", c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(<<<< FILEPATH REDACTED >>>>))


### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="VUT", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)


national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_VUT_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Marshall Islands model
rm(r)
rm(data_set_preds)

data_MHL <- pre_resampling_lf_data[country == "MHL"]

data_MHL[Master_UID %in% c("13458", "13465", "13467", "13480"), "island" := "Mejit Island"]
data_MHL[Master_UID %in% c("13461", "13466", "13468", "13481"), "island" := "Ailuk Atoll"]
data_MHL$weight <- as.numeric(data_MHL$weight)
data_MHL[shapefile == c("lf_MHL_Custom1_letoui"), c("island", "weight") := list("Mejit Island", 0.5)]
data_MHL[shapefile == c("lf_MHL_Custom2_letoui"), c("island", "weight") := list("Ailuk Atoll", 0.5)]

data_MHL <- data_MHL[!is.na(island)]

islands <- unique(data_MHL$island)
new <- data_MHL[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(islands)), ]
new$island <- rep(islands, nrow(new)/length(islands))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_poly", "lf_prev") := list(100, NA, NA)]
data_MHL <- rbind(data_MHL, new)
data_MHL$island_int <- as.integer(as.factor(data_MHL$island))
data_MHL$year_int <- data_MHL$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_poly ~ f(year_int, model = "rw1", replicate = island_int, scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_MHL, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = FALSE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_MHL$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_MHL, inv.logit(r$summary.linear.predictor[1:nrow(data_MHL), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_MHL), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_MHL), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_poly) & island == "Mejit Island"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_poly) & island == "Mejit Island"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_poly) & island == "Mejit Island"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_poly) & island == "Mejit Island"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Mejit Island, Marshall Islands (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_poly) & island == "Ailuk Atoll"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_poly) & island == "Ailuk Atoll"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_poly) & island == "Ailuk Atoll"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_poly) & island == "Ailuk Atoll"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Ailuk Atoll, Marshall Island (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, ncol = 1)


### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_MHL$Master_UID))
draw_preds <- matrix(nrow = nrow(data_MHL[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Get all-age population for entire country by year
popNumbers <- data.table(get_population(age_group_id = 22, location_id = 24, year_id = predict_years, sex_id = 3, single_year_age = F, gbd_round = 7, decomp_step = "iterative"))

# Assume that Ailuk and Mejit have represented a constant proportion of the national population across model years.
# Use the 2011 population data provided in the Marshall Islands LF elimination dossier.
national_pop_2011 <- 53158
ailuk_prop <- 339 / national_pop_2011
mejit_prop <- 348 / national_pop_2011

preds <- merge(preds, popNumbers[, c("year_id", "population")], by.x = "year", by.y = "year_id")
preds[island == "Mejit Island", island_population := mejit_prop * population]
preds[island == "Ailuk Atoll", island_population := ailuk_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$island_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$island <- "Marshall Islands National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- popNumbers$population
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

preds <- rbind(preds, national)

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = preds[island == "Marshall Islands National"], aes(x = year, y = mean, group = island), alpha = 1, size = 0.5, color = 1) + geom_ribbon(data = preds[island == "Marshall Islands National"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence: Marshall Islands (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_poly) & island == "Mejit Island"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_poly) & island == "Mejit Island"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_poly) & island == "Mejit Island"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_poly) & island == "Mejit Island"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Mejit Island, Marshall Islands (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_poly) & island == "Ailuk Atoll"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_poly) & island == "Ailuk Atoll"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_poly) & island == "Ailuk Atoll"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_poly) & island == "Ailuk Atoll"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Ailuk Atoll, Marshall Island (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, ncol = 1)
dev.off()

results_table <- setorderv(preds[island == "Marshall Islands National", c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="MHL", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_MHL[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Kiribati model
rm(r)
rm(data_set_preds)

data_KIR <- post_resampling_lf_data[country == "KIR"]

data_KIR[Master_UID %in% c(11478, 11477), "island" := "Makin"]
data_KIR[Master_UID %in% c(11472, 11479, 11474, 11473, 11481, 11469, 11480, 11490), "island" := "Kiritamiti"]
data_KIR[Master_UID %in% c(11484, 11482, 11488, 23827), "island" := "Tarawa"]
data_KIR[Master_UID %in% c(11476, 11475), "island" := "Tabuaeran"]
data_KIR[Master_UID %in% c(11494, 11491), "island" := "Teraina"]
data_KIR[Master_UID %in% c(23825), "island" := "Line Islands"]
data_KIR[Master_UID %in% c(11483, 21398, 23824, 23826), "island" := "Gilbert Islands"]
data_KIR[Master_UID %in% c(11471, 11485, 11486, 11487, 11489, 11492, 11493), island := "National"]

data_KIR <- data_KIR[island %in% c("Makin", "Kiritamiti", "Tarawa", "Tabuaeran", "Teraina")]

islands <- unique(data_KIR$island)
new <- data_KIR[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(islands)), ]
new$island <- rep(islands, nrow(new)/length(islands))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_KIR <- rbind(data_KIR, new)
data_KIR[is.na(island), island := "Other"]
data_KIR$island_int <- as.integer(as.factor(data_KIR$island))
data_KIR$year_int <- data_KIR$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "rw1", replicate = island_int, scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_KIR, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = FALSE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_KIR$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_KIR, inv.logit(r$summary.linear.predictor[1:nrow(data_KIR), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_KIR), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_KIR), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tarawa"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tarawa"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tarawa"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Tarawa"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tarawa, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Kiritamiti"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Kiritamiti"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Kiritamiti"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Kiritamiti"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Kiritamiti, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Makin"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Makin"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Makin"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Makin"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Makin, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tabuaeran"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tabuaeran"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tabuaeran"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Tabuaeran"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tabuaeran, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Teraina"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Teraina"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Teraina"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Teraina"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Teraina, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, ncol = 2)

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_KIR$Master_UID))
draw_preds <- matrix(nrow = nrow(data_KIR[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Get all-age population for entire country by year
popNumbers <- data.table(get_population(age_group_id = 22, location_id = 23, year_id = predict_years, sex_id = 3, single_year_age = F, gbd_round = 7, decomp_step = "iterative"))

# Assume that each island has represented a constant proportion of the national population across model years.
# Use the 2010 population data provided at http://www.mfed.gov.ki/sites/default/files/Revised%20Census%20Preliminary%20Report%202%20020516%20update%20%5B1306646%5D.pdf
national_pop_2010 <- 103058
tarawa_prop <- (6102 + 50182) / national_pop_2010
kiritamiti_prop <- 5586 / national_pop_2010
makin_prop <- 1798 / national_pop_2010
tabuaeran_prop <- 1960 / national_pop_2010
teraina_prop <- 1690 / national_pop_2010

preds <- merge(preds, popNumbers[, c("year_id", "population")], by.x = "year", by.y = "year_id")
preds[island == "Tarawa", island_population := tarawa_prop * population]
preds[island == "Kiritamiti", island_population := kiritamiti_prop * population]
preds[island == "Makin", island_population := makin_prop * population]
preds[island == "Tabuaeran", island_population := tabuaeran_prop * population]
preds[island == "Teraina", island_population := teraina_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$island_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "Kiribati National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- popNumbers$population
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tarawa"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tarawa"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tarawa"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Tarawa"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tarawa, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Kiritamiti"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Kiritamiti"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Kiritamiti"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Kiritamiti"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Kiritamiti, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Makin"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Makin"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Makin"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Makin"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Makin, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tabuaeran"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tabuaeran"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Tabuaeran"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Tabuaeran"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tabuaeran, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Teraina"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Teraina"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Teraina"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Teraina"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Teraina, Kiribati (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, ncol = 2)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="KIR", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_KIR[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run American Samoa model
rm(r)
rm(data_set_preds)

data_ASM <- post_resampling_lf_data[country == "ASM"]

gaul_code <- 11

subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

resampled_ASM <- data_ASM

admin2_shp <- rgdal::readOGR(dsn = get_admin_shapefile(admin_level = 2, version = modeling_shapefile_version))
admin2_shp@data$ADM1_CODE <- as.integer(levels(admin2_shp@data$ADM1_CODE))[admin2_shp@data$ADM1_CODE]
admin2_shp@data$ADM2_CODE <- as.integer(levels(admin2_shp@data$ADM2_CODE))[admin2_shp@data$ADM2_CODE]

asm_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "American Samoa", ], "SpatialPolygonsDataFrame"))
asm_adm1 <- fasterize(asm_df_sf, simple_raster2, field = "ADM1_CODE")
asm_adm2 <- fasterize(asm_df_sf, simple_raster2, field = "ADM2_CODE")

o <- over(SpatialPointsDataFrame(cbind(resampled_ASM$longitude, resampled_ASM$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_ASM), admin2_shp)
resampled_ASM$ADM1_NAME <- o$ADM1_NAME
resampled_ASM$ADM2_NAME <- o$ADM2_NAME

data_ASM_merged <- resampled_ASM

data_ASM_merged[Master_UID %in% c(17, 21, 27), ADM1_NAME := "Eastern"]
data_ASM_merged[Master_UID %in% c(13, 28, 4), ADM1_NAME := "Western"]
data_ASM_merged <- data_ASM_merged[!is.na(ADM1_NAME)] # Final geo-referenced data set for modeling

adm1s <- unique(data_ASM_merged$ADM1_NAME)
new <- data_ASM_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(adm1s)), ]
new$ADM1_NAME <- rep(adm1s, nrow(new)/length(adm1s))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_ASM_merged <- rbind(data_ASM_merged, new)
data_ASM_merged$ADM1_NAME <- droplevels(data_ASM_merged$ADM1_NAME)
data_ASM_merged$ADM1_NAME_int <- as.integer(data_ASM_merged$ADM1_NAME)
data_ASM_merged$year_int <- data_ASM_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "rw1", replicate = ADM1_NAME_int, scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_ASM_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_ASM_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_ASM_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_ASM_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_ASM_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_ASM_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Western District, American Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Eastern District, American Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Manu'a"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Manu'a"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Manu'a"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Manu'a"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Manu'a District, American Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, ncol = 1)

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_ASM_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_ASM_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Get all-age population for entire country by year
popNumbers <- data.table(get_population(age_group_id = 22, location_id = 298, year_id = predict_years, sex_id = 3, single_year_age = F, gbd_round = 7, decomp_step = "iterative"))

# Assume that each island has represented a constant proportion of the national population across model years.
# Population numbers drawn from https://www.census.gov/population/www/cen2010/cph-t/t-8tables/table1b.pdf.
national_pop_2010 <- 55519
Western_prop <- (494 + 182 + 51 + 162 + 250 + 615 + 224 + 254 + 247 + 47 + 108 + 1898 + 723 + 3195 + 1919 + 1182 + 698 + 550 + 8 + 1126 + 444 + 141 + 661 + 2450 + 193 + 965 + 299 + 7943 + 841 + 53 + 1447 + 1959) / national_pop_2010
Eastern_prop <- (0 + 524 + 495 + 54 + 646 + 96 + 920 + 18 + 855 + 359 + 2077 + 113 + 186 + 436 + 262 + 44 + 910 + 433 + 150 + 831 + 1737 + 113 + 892 + 448 + 0 + 164 + 425 + 399 + 3294 + 150 + 118 + 24 + 94 + 3656 + 75 + 297 + 2 + 405 + 684 + 48 + 74 + 640) / national_pop_2010
Manua_prop <- (162 + 117 + 183 + 153 + 176 + 172 + 5 + 175) / national_pop_2010

preds <- merge(preds, popNumbers[, c("year_id", "population")], by.x = "year", by.y = "year_id")
preds[ADM1_NAME == "Western", district_population := Western_prop * population]
preds[ADM1_NAME == "Eastern", district_population := Eastern_prop * population]
preds[ADM1_NAME == "Manu'a", district_population := Manua_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$district_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- popNumbers$population
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: American Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Western District, American Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Eastern District, American Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Manu'a"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Manu'a"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Manu'a"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Manu'a"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Manu'a District, American Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, ncol = 1)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="ASM", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_ASM_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Samoa model
rm(r)
rm(data_set_preds)

resampled_WSM <- post_resampling_lf_data[country == "WSM"]

gaul_code <- 245
subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

admin2_shp <- rgdal::readOGR(dsn = get_admin_shapefile(admin_level = 2, version = modeling_shapefile_version))
admin2_shp@data$ADM1_CODE <- as.integer(levels(admin2_shp@data$ADM1_CODE))[admin2_shp@data$ADM1_CODE]
admin2_shp@data$ADM2_CODE <- as.integer(levels(admin2_shp@data$ADM2_CODE))[admin2_shp@data$ADM2_CODE]

wsm_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "Samoa", ], "SpatialPolygonsDataFrame"))
wsm_adm1 <- fasterize(wsm_df_sf, simple_raster2, field = "ADM1_CODE")
wsm_adm2 <- fasterize(wsm_df_sf, simple_raster2, field = "ADM2_CODE")

o <- over(SpatialPointsDataFrame(cbind(resampled_WSM[!is.na(longitude), longitude], resampled_WSM[!is.na(latitude), latitude]), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_WSM[!is.na(longitude)]), admin2_shp)
resampled_WSM[!is.na(longitude), "ADM1_NAME"] <- o$ADM1_NAME
resampled_WSM[!is.na(longitude), "ADM2_NAME"] <- o$ADM2_NAME

data_WSM_merged <- resampled_WSM

# See https://mapcarta.com/ for geolocations of some villages

data_WSM_merged[Master_UID %in% c(15385, 15386, 15387, 15388, 15389, 15406, 15407), ADM1_NAME := "Atua"]
data_WSM_merged[Master_UID %in% c(15416, 15417), ADM1_NAME := "Tuamasaga"]
data_WSM_merged[Master_UID %in% c(15382, 15383), ADM1_NAME := "Gaga'emauga"]
data_WSM_merged[ADM1_NAME %in% c("Tuamasaga", "A'ana", "Aiga-i-le-Tai", "Atua", "Va'a-o-Fonoti"), island := "Upolu"]
data_WSM_merged[ADM1_NAME %in% c("Fa'asaleleaga", "Gaga'emauga", "Gagaifomauga", "Vaisigano", "Satupa'itea", "Palauli"), island := "Savai'i"]
data_WSM_merged[shapefile == "lf_ASM_Custom1_letoui", island := "Savai'i"]
data_WSM_merged[shapefile == "lf_ASM_Custom2_letoui", island := "Upolu"]
data_WSM_merged[shapefile == "lf_AMS_Custom3_letoui", island := "Upolu"]
data_WSM_merged[Master_UID %in% c(15357, 15434), island := "Upolu"]

data_WSM_merged$country <- "WSM"

data_WSM_merged <- data_WSM_merged[!is.na(island)] # Final geo-referenced data set for modeling

islands <- unique(data_WSM_merged$island)
new <- as.data.table(data_WSM_merged[1, ])
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(islands)), ]
new$island <- rep(islands, nrow(new)/2)
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_WSM_merged <- rbind(data_WSM_merged, new)
data_WSM_merged$island_int <- as.integer(as.factor(data_WSM_merged$island))
data_WSM_merged$year_int <- data_WSM_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "rw1", replicate = island_int, scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_WSM_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_WSM_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_WSM_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_WSM_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_WSM_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_WSM_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Upolu"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Upolu"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Upolu"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Upolu"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Upolu Island, Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Savai'i"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Savai'i"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Savai'i"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Savai'i"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Savai'i Island, Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, ncol = 2)

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_WSM_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_WSM_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Get all-age population for entire country by year
popNumbers <- data.table(get_population(age_group_id = 22, location_id = 27, year_id = predict_years, sex_id = 3, single_year_age = F, gbd_round = 7, decomp_step = "iterative"))

# Fit model at region level
# Use the 2011 district-level populations given at http://www.sbs.gov.ws/index.php/new-document-library?view=download&fileId=945
national_pop_2011 <- 187820
Upolu_prop <- (36735 + 62390 + 44293) / national_pop_2011
Savaii_prop <- (44402) / national_pop_2011

preds <- merge(preds, popNumbers[, c("year_id", "population")], by.x = "year", by.y = "year_id")
preds[island == "Upolu", district_population := Upolu_prop * population]
preds[island == "Savai'i", district_population := Savaii_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$district_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- popNumbers$population
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Upolu"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Upolu"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Upolu"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Upolu"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Upolu Island, Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Savai'i"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Savai'i"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & island == "Savai'i"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & island == "Savai'i"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Savai'i Island, Samoa (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, ncol = 2)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="WSM", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_WSM_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Cook Islands model
rm(r)
rm(data_set_preds)

resampled_COK <- post_resampling_lf_data[country == "COK"]

data_COK_merged <- resampled_COK

new <- data_COK_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_COK_merged <- rbind(data_COK_merged, new)

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_COK_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_COK_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_COK_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_COK_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_COK_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_COK_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds, aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds, aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds, aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Cook Islands (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, ncol = 1)

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = data_set_preds, aes(x = year, y = mean)) + geom_ribbon(data = data_set_preds, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Cook Islands (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

results_table <- setorderv(data_set_preds[is.na(had_lf_w_resamp), c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_COK_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_COK_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds_2018 <- inv.logit(draw_preds)[29,]
post_prob_2 <- c(country="COK", post_prob_2=(length(draw_preds_2018[draw_preds_2018 < 0.02])/length(draw_preds_2018)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

draw_preds <- cbind("year"=predict_years, as.data.table(inv.logit(draw_preds)))
write.table(draw_preds, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_COK_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Wallis & Futuna model
rm(r)
rm(data_set_preds)

data_WLF_merged <- pre_resampling_lf_data[country == "WLF"]

new <- data_WLF_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_poly", "lf_prev") := list(100, NA, NA)]
data_WLF_merged <- rbind(data_WLF_merged, new)

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_poly ~ f(year, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_WLF_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_WLF_merged$weight, verbose = TRUE, num.threads = 1)
r$waic$waic
summary(r)

data_set_preds <- cbind(data_WLF_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_WLF_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_WLF_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_WLF_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds, aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds, aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds, aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Wallis & Futuna (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
plot_grid(g1, ncol = 1)

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = data_set_preds, aes(x = year, y = mean)) + geom_ribbon(data = data_set_preds, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Wallis and Futuna (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) # + coord_cartesian(ylim=c(0, 1))
print(g)
dev.off()

results_table <- setorderv(data_set_preds[is.na(had_lf_poly), c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_WLF_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_WLF_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds_2018 <- inv.logit(draw_preds)[29,]
post_prob_2 <- c(country="WLF", post_prob_2=(length(draw_preds_2018[draw_preds_2018 < 0.02])/length(draw_preds_2018)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

draw_preds <- cbind("year"=predict_years, as.data.table(inv.logit(draw_preds)))
write.table(draw_preds, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_WLF_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run French Polynesia model
rm(r)
rm(data_set_preds)

resampled_PYF <- post_resampling_lf_data[country == "PYF"]

gaul_code <- 186

subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

PYF_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "French Polynesia", ], "SpatialPolygonsDataFrame"))
PYF_adm1 <- fasterize(PYF_df_sf, simple_raster2, field = "ADM1_CODE")
PYF_adm2 <- fasterize(PYF_df_sf, simple_raster2, field = "ADM2_CODE")

resampled_PYF$latitude <- as.numeric(resampled_PYF$latitude)
resampled_PYF$longitude <- as.numeric(resampled_PYF$longitude)

o <- over(SpatialPointsDataFrame(cbind(resampled_PYF$longitude, resampled_PYF$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_PYF), admin2_shp)
resampled_PYF$ADM1_NAME <- o$ADM1_NAME
resampled_PYF$ADM2_NAME <- o$ADM2_NAME

resampled_PYF[ADM1_NAME == "Îles Sous-le-Vent", ADM1_NAME := "Leeward Islands"]
resampled_PYF[ADM1_NAME == "Îles Marquises", ADM1_NAME := "Marquesas Islands"]
resampled_PYF[ADM1_NAME == "Îles du Vent", ADM1_NAME := "Windward Islands"]
resampled_PYF[Master_UID %in% c("7242", "7246", "7240", "7241", "23316"), ADM1_NAME := "Leeward Islands"]
resampled_PYF[Master_UID %in% c("7234", "23313"), ADM1_NAME := "Tuamotu-Gambier"]
resampled_PYF[Master_UID %in% c("7258", "23312"), ADM1_NAME := "Marquesas Islands"]
resampled_PYF[Master_UID %in% c("23315", "7229"), ADM1_NAME := "Austral Islands"]
resampled_PYF[shapefile == "lf_French_Polynesia_Moulia-Pelat_LF_1995_fixed" & location_code == 1, ADM1_NAME := "Leeward Islands"]
resampled_PYF[shapefile == "lf_gadm28_adm1" & location_code == 941, ADM1_NAME := "Windward Islands"]
resampled_PYF <- resampled_PYF[!is.na(ADM1_NAME)]

data_PYF_merged <- resampled_PYF

adm1s <- unique(data_PYF_merged$ADM1_NAME)
new <- data_PYF_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(adm1s)), ]
new$ADM1_NAME <- rep(adm1s, nrow(new)/length(adm1s))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_PYF_merged <- rbind(data_PYF_merged, new)
data_PYF_merged$ADM1_NAME <- droplevels(data_PYF_merged$ADM1_NAME)
data_PYF_merged$ADM1_int <- as.integer(as.factor(data_PYF_merged$ADM1_NAME))
data_PYF_merged$year_int <- data_PYF_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "rw1", replicate = ADM1_int, scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_PYF_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = FALSE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_PYF_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_PYF_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_PYF_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_PYF_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_PYF_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Leeward Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Leeward Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Leeward Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Leeward Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Leeward Islands, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Marquesas Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Marquesas Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Marquesas Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Marquesas Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Marquesas Islands, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Windward Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Windward Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Windward Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Windward Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Windward Islands, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tuamotu-Gambier"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tuamotu-Gambier"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tuamotu-Gambier"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Tuamotu-Gambier"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tuamotu-Gambier, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Austral Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Austral Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Austral Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Austral Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Austral Islands, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, ncol = 3)


### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_PYF_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_PYF_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Fit model at region level
# Use the ADM1-level populations given at http://www.ispf.pf/docs/default-source/rp2017/repart_poplegale_communes_2017_v3.pdf
national_pop_2012 <- 268270
Leeward_prop <- (34622) / national_pop_2012
Marquesas_prop <- (9264) / national_pop_2012
Windward_prop <- (200881) / national_pop_2012
Tuamotu_Gambier_prop <- (16664) / national_pop_2012
Austral_prop <- (6839) / national_pop_2012

preds$population <- national_pop_2012
preds[ADM1_NAME == "Leeward Islands", district_population := Leeward_prop * population]
preds[ADM1_NAME == "Marquesas Islands", district_population := Marquesas_prop * population]
preds[ADM1_NAME == "Windward Islands", district_population := Windward_prop * population]
preds[ADM1_NAME == "Tuamotu-Gambier", district_population := Tuamotu_Gambier_prop * population]
preds[ADM1_NAME == "Austral Islands", district_population := Austral_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$district_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- national_pop_2012
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Leeward Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Leeward Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Leeward Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Leeward Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Leeward Islands, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Marquesas Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Marquesas Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Marquesas Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Marquesas Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Marquesas Islands, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Windward Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Windward Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Windward Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Windward Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Windward Islands, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tuamotu-Gambier"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tuamotu-Gambier"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Tuamotu-Gambier"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Tuamotu-Gambier"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tuamotu-Gambier, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Austral Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Austral Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Austral Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Austral Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Austral Islands, French Polynesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, ncol = 3)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="PYF", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_PYF_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Tonga model
rm(r)
rm(data_set_preds)

resampled_TON <- pre_resampling_lf_data[country == "TON"]

data_TON_merged <- resampled_TON

new <- data_TON_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_poly", "lf_prev") := list(100, NA, NA)]
data_TON_merged <- rbind(data_TON_merged, new)

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_poly ~ f(year, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_TON_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_TON_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_TON_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_TON_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_TON_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_TON_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_point(data = data_set_preds, aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds, aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds, aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tonga (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

results_table <- setorderv(data_set_preds[is.na(had_lf_poly), c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_TON_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_TON_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds_2018 <- inv.logit(draw_preds)[29,]
post_prob_2 <- c(country="TON", post_prob_2=(length(draw_preds_2018[draw_preds_2018 < 0.02])/length(draw_preds_2018)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

draw_preds <- cbind("year"=predict_years, as.data.table(inv.logit(draw_preds)))
write.table(draw_preds, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_TON_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Tuvalu model
rm(r)
rm(data_set_preds)

resampled_TUV <- pre_resampling_lf_data[country == "TUV"]

# Population stats available from https://tuvalu.prism.spc.int/index.php/social/population.

data_TUV_merged <- resampled_TUV

new <- data_TUV_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_poly", "lf_prev") := list(100, NA, NA)]
data_TUV_merged <- rbind(data_TUV_merged, new)

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_poly ~ f(year, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_TUV_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = FALSE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_TUV_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_TUV_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_TUV_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_TUV_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_TUV_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_point(data = data_set_preds, aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds, aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds, aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Tuvalu (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

results_table <- setorderv(data_set_preds[is.na(had_lf_poly), c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_TUV_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_TUV_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds_2018 <- inv.logit(draw_preds)[29,]
post_prob_2 <- c(country="TUV", post_prob_2=(length(draw_preds_2018[draw_preds_2018 < 0.02])/length(draw_preds_2018)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

draw_preds <- cbind("year"=predict_years, as.data.table(inv.logit(draw_preds)))
write.table(draw_preds, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_TUV_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Fiji model
rm(r)
rm(data_set_preds)

resampled_FJI <- post_resampling_lf_data[country == "FJI"]

gaul_code <- 74

subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

FJI_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "Fiji", ], "SpatialPolygonsDataFrame"))
FJI_adm1 <- fasterize(FJI_df_sf, simple_raster2, field = "ADM1_CODE")
FJI_adm2 <- fasterize(FJI_df_sf, simple_raster2, field = "ADM2_CODE")

resampled_FJI$latitude <- as.numeric(resampled_FJI$latitude)
resampled_FJI$longitude <- as.numeric(resampled_FJI$longitude)

o <- over(SpatialPointsDataFrame(cbind(resampled_FJI$longitude, resampled_FJI$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_FJI), admin2_shp)
resampled_FJI$ADM1_NAME <- o$ADM1_NAME
resampled_FJI$ADM2_NAME <- o$ADM2_NAME

resampled_FJI[Master_UID %in% c("7049", "7100"), latitude := -latitude]
resampled_FJI[Master_UID %in% c("7042", "7090", "7053", "7078", "7061", "7098", "7064", "7079", "7049", "7100"), ADM1_NAME := "Eastern"]
resampled_FJI[Master_UID %in% c("7160"), ADM1_NAME := "Northern"]
resampled_FJI[Master_UID %in% c("7163"), ADM1_NAME := "Rotuma"]

resampled_FJI[Master_UID %in% c("7156"), ADM1_NAME := "Western"]
resampled_FJI[Master_UID %in% c("7153", "23311"), ADM1_NAME := "Central"]
resampled_FJI[Master_UID %in% c("7154"), ADM1_NAME := "Eastern"]
resampled_FJI[Master_UID %in% c("7155", "7175"), ADM1_NAME := "Northern"]
resampled_FJI[Master_UID %in% c("SR_C_492") & is.na(ADM1_NAME), ADM1_NAME := "Eastern"]

data_FJI_merged <- resampled_FJI

adm1s <- unique(data_FJI_merged$ADM1_NAME)
new <- data_FJI_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(adm1s)), ]
new$ADM1_NAME <- rep(adm1s, nrow(new)/length(adm1s))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_FJI_merged <- rbind(data_FJI_merged, new)
data_FJI_merged$ADM1_NAME <- droplevels(data_FJI_merged$ADM1_NAME)
data_FJI_merged$ADM1_int <- as.integer(as.factor(data_FJI_merged$ADM1_NAME))
data_FJI_merged$year_int <- data_FJI_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "rw1", replicate = ADM1_int, scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_FJI_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_FJI_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_FJI_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_FJI_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_FJI_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_FJI_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Western Division, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Central"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Central"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Central"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Central"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Central Division, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Eastern Division, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Northern"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Northern"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Northern"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Northern"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Northern Division, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Rotuma"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Rotuma"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Rotuma"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Rotuma"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Rotuma, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, ncol = 2)


### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_FJI_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_FJI_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Fit model at division level
# Use population numbers from 2017 census: https://www.statsfiji.gov.fj/index.php/census-2017

national_pop_2017 <- 886470
Western_prop <- (337041) / national_pop_2017
Central_prop <- (378284) / national_pop_2017
Eastern_prop <- (37648) / national_pop_2017
Northern_prop <- (131914) / national_pop_2017
Rotuma_prop <- (1583) / national_pop_2017

preds$population <- national_pop_2017
preds[ADM1_NAME == "Western", district_population := Western_prop * population]
preds[ADM1_NAME == "Central", district_population := Central_prop * population]
preds[ADM1_NAME == "Eastern", district_population := Eastern_prop * population]
preds[ADM1_NAME == "Northern", district_population := Northern_prop * population]
preds[ADM1_NAME == "Rotuma", district_population := Rotuma_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$district_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- national_pop_2017
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population


results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(<<<< FILEPATH REDACTED >>>>)

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Western"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Western Division, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Central"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Central"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Central"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Central"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Central Division, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Eastern"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Eastern Division, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Northern"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Northern"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Northern"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Northern"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Northern Division, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Rotuma"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Rotuma"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Rotuma"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Rotuma"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Rotuma, Fiji (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, ncol = 2)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="FJI", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_FJI_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Federated States of Micronesia model
rm(r)
rm(data_set_preds)

resampled_FSM <- post_resampling_lf_data[country == "FSM"]

gaul_code <- 78

subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

FSM_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "Micronesia", ], "SpatialPolygonsDataFrame"))
FSM_adm1 <- fasterize(FSM_df_sf, simple_raster2, field = "ADM1_CODE")
FSM_adm2 <- fasterize(FSM_df_sf, simple_raster2, field = "ADM2_CODE")

resampled_FSM$latitude <- as.numeric(resampled_FSM$latitude)
resampled_FSM$longitude <- as.numeric(resampled_FSM$longitude)

o <- over(SpatialPointsDataFrame(cbind(resampled_FSM$longitude, resampled_FSM$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_FSM), admin2_shp)
resampled_FSM$ADM1_NAME <- o$ADM1_NAME
resampled_FSM$ADM2_NAME <- o$ADM2_NAME

library(maptools)
getKMLcoordinates(<<<< FILEPATH REDACTED >>>>) # state boundaries (longitude)

resampled_FSM[Master_UID %in% c("13520", "13521"), ADM1_NAME := "Chuuk"]
resampled_FSM[Master_UID %in% c("13523"), ADM1_NAME := "Pohnpei"]
resampled_FSM[is.na(ADM1_NAME) & shapefile == "lf_gadm28_adm1" & location_code == 1854, ADM1_NAME := "Chuuk"]
resampled_FSM[is.na(ADM1_NAME) & shapefile == "lf_gadm28_adm1" & location_code == 1857, ADM1_NAME := "Yap"]
# resampled_FSM[is.na(ADM1_NAME) & latitude < 7.5 & longitude > 151, ADM1_NAME := "Chuuk"]
# resampled_FSM[is.na(ADM1_NAME) & latitude > 9.5 & longitude < 139, ADM1_NAME := "Yap"]
resampled_FSM[is.na(ADM1_NAME) & longitude <= 148, ADM1_NAME := "Yap"]
resampled_FSM[is.na(ADM1_NAME) & longitude > 148 & longitude <= 154, ADM1_NAME := "Chuuk"]
resampled_FSM[is.na(ADM1_NAME) & longitude > 154 & longitude <= 162, ADM1_NAME := "Pohnpei"]
resampled_FSM[is.na(ADM1_NAME) & longitude > 162, ADM1_NAME := "Kosrae"]

data_FSM_merged <- resampled_FSM

adm1s <- unique(data_FSM_merged$ADM1_NAME)
new <- data_FSM_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(adm1s)), ]
new$ADM1_NAME <- rep(adm1s, nrow(new)/length(adm1s))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_FSM_merged <- rbind(data_FSM_merged, new)
data_FSM_merged$ADM1_NAME <- droplevels(data_FSM_merged$ADM1_NAME)
data_FSM_merged$ADM1_int <- as.integer(as.factor(data_FSM_merged$ADM1_NAME))
data_FSM_merged$year_int <- data_FSM_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "rw1", replicate = ADM1_int, scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_FSM_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_FSM_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_FSM_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_FSM_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_FSM_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_FSM_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Chuuk"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Chuuk"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Chuuk"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Chuuk"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Chuuk, Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Kosrae"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Kosrae"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Kosrae"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Kosrae"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Kosrae, Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Pohnpei"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Pohnpei"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Pohnpei"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Pohnpei"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Pohnpei, Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Yap"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Yap"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Yap"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Yap"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Yap, Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, ncol = 2)


### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_FSM_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_FSM_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Fit model at region level
# Population stats: https://web.archive.org/web/20120629052438/http://www.fsmgov.org/info/people.html
national_pop_2000 <- 107008
Yap_prop <- (11241) / national_pop_2000
Pohnpei_prop <- (34486) / national_pop_2000
Kosrae_prop <- (7686) / national_pop_2000
Chuuk_prop <- (53595) / national_pop_2000

preds$population <- national_pop_2000
preds[ADM1_NAME == "Yap", district_population := Yap_prop * population]
preds[ADM1_NAME == "Pohnpei", district_population := Pohnpei_prop * population]
preds[ADM1_NAME == "Kosrae", district_population := Kosrae_prop * population]
preds[ADM1_NAME == "Chuuk", district_population := Chuuk_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$district_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- national_pop_2000
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Chuuk"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Chuuk"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Chuuk"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Chuuk"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Chuuk, Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Kosrae"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Kosrae"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Kosrae"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Kosrae"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Kosrae, Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Pohnpei"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Pohnpei"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Pohnpei"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Pohnpei"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Pohnpei, Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Yap"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Yap"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Yap"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Yap"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Yap, Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, ncol = 2)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="FSM", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_FSM_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run New Caledonia model
rm(r)
rm(data_set_preds)

resampled_NCL <- post_resampling_lf_data[country == "NCL"]

gaul_code <- 161

subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

NCL_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "New Caledonia", ], "SpatialPolygonsDataFrame"))
NCL_adm1 <- fasterize(NCL_df_sf, simple_raster2, field = "ADM1_CODE")
NCL_adm2 <- fasterize(NCL_df_sf, simple_raster2, field = "ADM2_CODE")

resampled_NCL$latitude <- as.numeric(resampled_NCL$latitude)
resampled_NCL$longitude <- as.numeric(resampled_NCL$longitude)

o <- over(SpatialPointsDataFrame(cbind(resampled_NCL$longitude, resampled_NCL$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_NCL), admin2_shp)
resampled_NCL$ADM1_NAME <- o$ADM1_NAME
resampled_NCL$ADM2_NAME <- o$ADM2_NAME

resampled_NCL[ADM1_NAME == "Îles Loyauté", ADM1_NAME := "Loyalty Islands"]
resampled_NCL[is.na(ADM1_NAME) | ADM1_NAME != "Loyalty Islands", ADM1_NAME := "Other"]

data_NCL_merged <- resampled_NCL

adm1s <- unique(data_NCL_merged$ADM1_NAME)
new <- data_NCL_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(adm1s)), ]
new$ADM1_NAME <- rep(adm1s, nrow(new)/length(adm1s))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_NCL_merged <- rbind(data_NCL_merged, new)
data_NCL_merged$ADM1_NAME <- droplevels(data_NCL_merged$ADM1_NAME)
data_NCL_merged$ADM1_int <- as.integer(as.factor(data_NCL_merged$ADM1_NAME))
data_NCL_merged$year_int <- data_NCL_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01)))) + f(ADM1_int, model = "iid")
r <- inla(formula, family = "binomial", data = data_NCL_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_NCL_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_NCL_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_NCL_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_NCL_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_NCL_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Loyalty Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Loyalty Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Loyalty Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Loyalty Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Loyalty Islands, New Caledonia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Other"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Other"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Other"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Other"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Other Islands, New Caledonia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, ncol = 1)

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_NCL_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_NCL_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Fit model at division level
# Use population numbers from 2010 estimates: http://www.isee.nc/component/phocadownload/category/195-donnees?download=765:la-population-annuelle-estimee

national_pop_2010 <- 247400
Loyalty_prop <- 17550 / national_pop_2010
Other_prop <- 1 - Loyalty_prop

preds$population <- national_pop_2010
preds[ADM1_NAME == "Loyalty Islands", district_population := Loyalty_prop * population]
preds[ADM1_NAME == "Other", district_population := Other_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$district_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "National"
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- national_pop_2010
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: New Caledonia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Federated States of Micronesia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Loyalty Islands"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Loyalty Islands"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Loyalty Islands"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Loyalty Islands"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Loyalty Islands, New Caledonia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Other"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Other"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Other"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Other"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Other Islands, New Caledonia (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, ncol = 1)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="NCL", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_NCL_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Maldives model
rm(r)
rm(data_set_preds)

resampled_MDV <- post_resampling_lf_data[country == "MDV"]

data_MDV_merged <- resampled_MDV

new <- data_MDV_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_MDV_merged <- rbind(data_MDV_merged, new)
data_MDV_merged$year_int <- data_MDV_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01))))
r <- inla(formula, family = "binomial", data = data_MDV_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_MDV_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_MDV_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_MDV_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_MDV_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_MDV_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp), ], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp), ], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp), ], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp), ], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Maldives (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, ncol = 1)

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = data_set_preds, aes(x = year, y = mean)) + geom_ribbon(data = data_set_preds, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Maldives (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

results_table <- setorderv(data_set_preds[is.na(had_lf_w_resamp), c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_MDV_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_MDV_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds_2018 <- inv.logit(draw_preds)[29,]
post_prob_2 <- c(country="MDV", post_prob_2=(length(draw_preds_2018[draw_preds_2018 < 0.02])/length(draw_preds_2018)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

draw_preds <- cbind("year"=predict_years, as.data.table(inv.logit(draw_preds)))
write.table(draw_preds, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_MDV_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Guyana model
rm(r)
rm(data_set_preds)

resampled_GUY <- post_resampling_lf_data[country == "GUY"]

gaul_code <- 96

subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

GUY_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "Guyana", ], "SpatialPolygonsDataFrame"))
GUY_adm1 <- fasterize(GUY_df_sf, simple_raster2, field = "ADM1_CODE")
GUY_adm2 <- fasterize(GUY_df_sf, simple_raster2, field = "ADM2_CODE")

resampled_GUY$latitude <- as.numeric(resampled_GUY$latitude)
resampled_GUY$longitude <- as.numeric(resampled_GUY$longitude)

o <- over(SpatialPointsDataFrame(cbind(resampled_GUY$longitude, resampled_GUY$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_GUY), admin2_shp)
resampled_GUY$ADM1_NAME <- o$ADM1_NAME
resampled_GUY$ADM2_NAME <- o$ADM2_NAME

data_GUY_merged <- resampled_GUY

adm1s <- unique(data_GUY_merged$ADM1_NAME)
new <- data_GUY_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(adm1s)), ]
new$ADM1_NAME <- rep(adm1s, nrow(new)/length(adm1s))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_GUY_merged <- rbind(data_GUY_merged, new)
data_GUY_merged$ADM1_NAME <- droplevels(data_GUY_merged$ADM1_NAME)
data_GUY_merged$ADM1_int <- as.integer(as.factor(data_GUY_merged$ADM1_NAME))
data_GUY_merged$year_int <- data_GUY_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.01)))) + f(ADM1_int, model = "iid")
r <- inla(formula, family = "binomial", data = data_GUY_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_GUY_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_GUY_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_GUY_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_GUY_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_GUY_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

data_set_preds[ADM1_NAME == "Upper Demerara/Berbice (Region n10)", ADM1_NAME := "Upper Demerara-Berbice"]
data_set_preds[ADM1_NAME == "East Berbice/Corentyne (Region n6)", ADM1_NAME := "East Berbice-Corentyne"]
data_set_preds[ADM1_NAME == "Mahaica Berbice (Region n5)", ADM1_NAME := "Mahaica-Berbice"]
data_set_preds[ADM1_NAME == "Cuyuni/Mazaruni (Region n7)", ADM1_NAME := "Cuyuni-Mazaruni"]
data_set_preds[ADM1_NAME == "Essequibo Islands/West Demerara (Region n3)", ADM1_NAME := "Essequibo Islands-West Demerara"]
data_set_preds[ADM1_NAME == "Demerara/Mahaica (Region n4)", ADM1_NAME := "Demerara-Mahaica"]

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Upper Demerara-Berbice"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Upper Demerara-Berbice"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Upper Demerara-Berbice"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Upper Demerara-Berbice"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Upper Demerara-Berbice, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "East Berbice-Corentyne"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "East Berbice-Corentyne"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "East Berbice-Corentyne"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "East Berbice-Corentyne"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: East Berbice-Corentyne, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Mahaica-Berbice"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Mahaica-Berbice"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Mahaica-Berbice"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Mahaica-Berbice"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Mahaica-Berbice, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Cuyuni-Mazaruni"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Cuyuni-Mazaruni"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Cuyuni-Mazaruni"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Cuyuni-Mazaruni"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Cuyuni-Mazaruni, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Essequibo Islands-West Demerara"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Essequibo Islands-West Demerara"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Essequibo Islands-West Demerara"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Essequibo Islands-West Demerara"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Essequibo Islands-West Demerara, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g6 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Demerara-Mahaica"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Demerara-Mahaica"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Demerara-Mahaica"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Demerara-Mahaica"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Demerara-Mahaica, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, g6, ncol = 2)

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_GUY_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_GUY_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# Use population numbers from WorldPop raster
pops <- data.table("ADM1_CODE" = unique(GUY_adm1), "population" = NA)
pops$population <- as.integer(pops$population)

rast <- pop_raster$WorldPop_total_global_stack.4

for (i in unique(GUY_adm1)) {
  idx <- which(values(GUY_adm1) == i)
  pops[ADM1_CODE == i, population := as.integer(sum(rast[idx]))]
}

pops <- merge(pops, unique(admin2_shp@data[, c("ADM1_NAME", "ADM1_CODE")]))

national_pop <- sum(pops$population)

national_pop_2015 <- 693363
UDB_prop <- 40430 / national_pop_2015
EBC_prop <- 109459 / national_pop_2015
MB_prop <- 58399 / national_pop_2015
CM_prop <- 7264 / national_pop_2015
EIWD_prop <- 93910 / national_pop_2015
DM_prop <- 289631 / national_pop_2015

preds$population <- national_pop_2015
preds[ADM1_NAME == "Upper Demerara-Berbice", district_population := UDB_prop * population]
preds[ADM1_NAME == "East Berbice-Corentyne", district_population := EBC_prop * population]
preds[ADM1_NAME == "Mahaica-Berbice", district_population := MB_prop * population]
preds[ADM1_NAME == "Cuyuni-Mazaruni", district_population := CM_prop * population]
preds[ADM1_NAME == "Essequibo Islands-West Demerara", district_population := EIWD_prop * population]
preds[ADM1_NAME == "Demerara-Mahaica", district_population := DM_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$district_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "National"
national$population <- as.integer(national$population)
national[, c("mean", "quant0.025", "quant0.975", "population") := NA]
national$population <- national_pop_2015
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) # + coord_cartesian(ylim=c(0, 1))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Upper Demerara-Berbice"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Upper Demerara-Berbice"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Upper Demerara-Berbice"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Upper Demerara-Berbice"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Upper Demerara-Berbice, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "East Berbice-Corentyne"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "East Berbice-Corentyne"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "East Berbice-Corentyne"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "East Berbice-Corentyne"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: East Berbice-Corentyne, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Mahaica-Berbice"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Mahaica-Berbice"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Mahaica-Berbice"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Mahaica-Berbice"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Mahaica-Berbice, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Cuyuni-Mazaruni"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Cuyuni-Mazaruni"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Cuyuni-Mazaruni"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Cuyuni-Mazaruni"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Cuyuni-Mazaruni, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Essequibo Islands-West Demerara"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Essequibo Islands-West Demerara"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Essequibo Islands-West Demerara"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Essequibo Islands-West Demerara"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Essequibo Islands-West Demerara, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g6 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Demerara-Mahaica"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Demerara-Mahaica"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM1_NAME == "Demerara-Mahaica"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM1_NAME == "Demerara-Mahaica"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Demerara-Mahaica, Guyana (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, g6, ncol = 2)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="GUY", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

write.csv(data_GUY_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))


####################################################################################
###### Run Brazil model
rm(r)
rm(data_set_preds)

resampled_BRA <- post_resampling_lf_data[country == "BRA"]

simple_polygon_list <- load_simple_polygon(gaul_list = 33, buffer = 1, tolerance = 0.4, use_premade = use_premade, shapefile_version = modeling_shapefile_version)
subset_shape <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]

raster_list <- build_simple_raster_pop(subset_shape)

simple_raster <- raster_list[["simple_raster"]]
pop_raster <- raster_list[["pop_raster"]]

gaul_code <- 33
subset_shape2 <- subset_shape[subset_shape@data$ADM0_CODE == gaul_code,]
pop_raster2 <- raster::crop(pop_raster, subset_shape2)
simple_raster2 <- raster::crop(simple_raster, subset_shape2)

BRA_df_sf <- st_as_sf(as(admin2_shp[admin2_shp@data$ADM0_NAME == "Brazil", ], "SpatialPolygonsDataFrame"))
BRA_adm1 <- fasterize(BRA_df_sf, simple_raster2, field = "ADM1_CODE")
BRA_adm2 <- fasterize(BRA_df_sf, simple_raster2, field = "ADM2_CODE")

resampled_BRA$latitude <- as.numeric(resampled_BRA$latitude)
resampled_BRA$longitude <- as.numeric(resampled_BRA$longitude)

o <- over(SpatialPointsDataFrame(cbind(resampled_BRA$longitude, resampled_BRA$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"), data = resampled_BRA), admin2_shp)
resampled_BRA$ADM1_NAME <- o$ADM1_NAME
resampled_BRA$ADM2_NAME <- o$ADM2_NAME

data_BRA_merged <- resampled_BRA

# Restrict to Pernambuco state (the only endemic state post-1999)
data_BRA_merged <- data_BRA_merged[ADM1_NAME == "Pernambuco"]
data_BRA_merged <- data_BRA_merged[(ADM2_NAME %in% c("Cabo", "Camaragibe", "Jaboatão dos Guararapes", "Moreno", "Olinda", "Recife"))]

adm2s <- unique(data_BRA_merged$ADM2_NAME)
new <- data_BRA_merged[1, ]
new <- new[rep(seq_len(nrow(new)), each = length(1990:2018)), ]
new$year <- 1990:2018
new <- new[rep(seq_len(nrow(new)), each = length(adm2s)), ]
new$ADM2_NAME <- rep(adm2s, nrow(new)/length(adm2s))
new[, c("data_collect_method", "source", "Master_UID", "diagnostic", "shapefile") := NA]
new[, c("N", "had_lf_w_resamp", "lf_prev") := list(100, NA, NA)]
data_BRA_merged <- rbind(data_BRA_merged, new)
data_BRA_merged$ADM2_NAME <- droplevels(data_BRA_merged$ADM2_NAME)
data_BRA_merged$ADM2_int <- as.integer(as.factor(data_BRA_merged$ADM2_NAME))
data_BRA_merged$year_int <- data_BRA_merged$year - 1989

inla.setOption("enable.inla.argument.weights", TRUE)

formula <- had_lf_w_resamp ~ f(year_int, model = "ar1", replicate = ADM2_int)
r <- inla(formula, family = "binomial", data = data_BRA_merged, control.predictor = list(compute = TRUE, link = 1), control.compute = list(waic = TRUE, config = TRUE, openmp.strategy = "default", smtp = "taucs"), Ntrial = N, weights = data_BRA_merged$weight, verbose = TRUE)
r$waic$waic

data_set_preds <- cbind(data_BRA_merged, inv.logit(r$summary.linear.predictor[1:nrow(data_BRA_merged), c("mean")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_BRA_merged), c("0.025quant")]))
data_set_preds <- cbind(data_set_preds, inv.logit(r$summary.linear.predictor[1:nrow(data_BRA_merged), c("0.975quant")]))
colnames(data_set_preds)[(ncol(data_set_preds) - 2):ncol(data_set_preds)] <- c("mean", "quant0.025", "quant0.975")

g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Cabo"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Cabo"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Cabo"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Cabo"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Cabo, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Jaboatão dos Guararapes"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Jaboatão dos Guararapes"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Jaboatão dos Guararapes"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Jaboatão dos Guararapes"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Jaboatão dos Guararapes, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Recife"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Recife"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Recife"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Recife"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Recife, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Olinda"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Olinda"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Olinda"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Olinda"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Olinda, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Camaragibe"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Camaragibe"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Camaragibe"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Camaragibe"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Camaragibe, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g6 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Moreno"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Moreno"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Moreno"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Moreno"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Moreno, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, g6, ncol = 2)

### Calculate country-level aggregated prevalence
draws <- inla.posterior.sample(1000, r)
preds <- data_set_preds[is.na(Master_UID)]
predict_years <- 1990:2018

pred_indices <- which(is.na(data_BRA_merged$Master_UID))
draw_preds <- matrix(nrow = nrow(data_BRA_merged[is.na(Master_UID), ]), ncol = length(draws))

for (d in 1:length(draws)) {
  draw_preds[, d] <- draws[[d]]$latent[pred_indices]
}

draw_preds <- inv.logit(draw_preds)

# 2011 city populations from "ESTIMATIVAS DA POPULAÇÃO RESIDENTE NOS MUNICÍPIOS BRASILEIROS COM DATA DE REFERÊNCIA EM 1º DE JULHO DE 2011" (in Portuguese). Instituto Brasileiro de Geografia e Estatística. 30 August 2011. Archived from the original (PDF) on 31 August 2011. Retrieved 31 August 2011.
# http://www.ibge.gov.br/home/estatistica/populacao/estimativa2011/POP2011_DOU.pdf

# 2019 city populations from "ESTIMATIVAS DA POPULAÇÃO RESIDENTE PARA OS MUNICÍPIOS E PARA AS UNIDADES DA FEDERAÇÃO BRASILEIROS COM DATA DE REFERÊNCIA EM 1º DE JULHO DE 2019"
# ftp://ftp.ibge.gov.br/Estimativas_de_Populacao/Estimativas_2019/estimativa_TCU_2019_20200116.pdf
# National estimate from UN: https://population.un.org/wpp/DVD/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/WPP2017_POP_F01_1_TOTAL_POPULATION_BOTH_SEXES.xlsx

source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)

lbd_gbd <- fread(<<<< FILEPATH REDACTED >>>>)

locs <- get_location_metadata(gbd_round_id = 7, location_set_id = 35, decomp_step = "iterative")
gbd <- get_population(decomp_step = "iterative", gbd_round_id = 7, year_id = 1990:2017, location_id = c(135, 4766))

national_pop_2019 <- 210147125
Camaragibe_prop <- 157828 / national_pop_2019
Cabo_prop <- 207048 / national_pop_2019
JabGua_prop <- 702298 / national_pop_2019
Recife_prop <- 1645727 / national_pop_2019
Olinda_prop <- 392482 / national_pop_2019
Moreno_prop <- 62784 / national_pop_2019

preds$population <- national_pop_2011
preds[ADM2_NAME == "Cabo", district_population := Cabo_prop * population]
preds[ADM2_NAME == "Camaragibe", district_population := Camaragibe_prop * population]
preds[ADM2_NAME == "Jaboatão dos Guararapes", district_population := JabGua_prop * population]
preds[ADM2_NAME == "Recife", district_population := Recife_prop * population]
preds[ADM2_NAME == "Olinda", district_population := Olinda_prop * population]
preds[ADM2_NAME == "Moreno", district_population := Moreno_prop * population]

# National aggregates
draw_preds_pops <- draw_preds * preds$district_population
annual_preds <- matrix(nrow = length(predict_years), ncol = length(draws))

for (a in 1:length(predict_years)) {
  year_indices <- which(preds$year == predict_years[a])
  for (d in 1:length(draws)) {
    annual_preds[a, d] <- sum(draw_preds_pops[year_indices, d])
  }
}

national <- preds[1, ]
national <- national[rep(seq_len(nrow(national)), each = length(1990:2018)), ]
national$year <- 1990:2018
national$ADM1_NAME <- "National"
national$population <- as.integer(national$population)
national[, c("mean", "quant0.025", "quant0.975") := NA]
national$mean <- rowMeans(annual_preds) / national$population
national$quant0.025 <- rowQuantiles(annual_preds, probs = 0.025) / national$population
national$quant0.975 <- rowQuantiles(annual_preds, probs = 0.975) / national$population

ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

# save results figure
pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 8.5, width = 11)
g <- ggplot() + theme_classic() + geom_line(data = national, aes(x = year, y = mean)) + geom_ribbon(data = national, aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
print(g)
dev.off()

pdf(paste0(<<<< FILEPATH REDACTED >>>>), height = 17, width = 22)
g1 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Cabo"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Cabo"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Cabo"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Cabo"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Cabo, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g2 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Jaboatão dos Guararapes"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Jaboatão dos Guararapes"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Jaboatão dos Guararapes"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Jaboatão dos Guararapes"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Jaboatão dos Guararapes, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g3 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Recife"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Recife"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Recife"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Recife"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Recife, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g4 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Olinda"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Olinda"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Olinda"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Olinda"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Olinda, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g5 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Camaragibe"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Camaragibe"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Camaragibe"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Camaragibe"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Camaragibe, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
g6 <- ggplot() + theme_classic() + geom_point(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Moreno"], aes(x = year, y = mean), color = 1) + geom_ribbon(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Moreno"], aes(x = year, ymin = quant0.025, ymax = quant0.975), alpha = 0.2) + geom_line(data = data_set_preds[is.na(had_lf_w_resamp) & ADM2_NAME == "Moreno"], aes(x = year, y = mean), color = 1) + geom_point(data = data_set_preds[!is.na(had_lf_w_resamp) & ADM2_NAME == "Moreno"], aes(x = year, y = lf_prev), color = 2) + xlab("Year") + ylab("Prevalence") + ggtitle("Lymphatic Filariasis Prevalence Model: Moreno, Brazil (1990-2018)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)) + coord_cartesian(ylim = c(0, 1))
plot_grid(g1, g2, g3, g4, g5, g6, ncol = 2)
dev.off()

results_table <- setorderv(national[, c("country", "year", "mean", "quant0.025", "quant0.975")], "year")
write.csv(results_table, paste0(<<<< FILEPATH REDACTED >>>>))

### Calculate country-level aggregated prevalence
national_draws <- (annual_preds/national$population)[29,]
post_prob_2 <- c(country="BRA", post_prob_2=(length(national_draws[national_draws < 0.02])/length(national_draws)))

write.table(as.data.table(t(post_prob_2)), file=paste0(<<<< FILEPATH REDACTED >>>>), append=TRUE, row.names=FALSE, sep=",", col.names=FALSE)

national_draws <- (annual_preds/national$population)
national_draws <- cbind("year"=predict_years, as.data.table(national_draws))
write.table(national_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")
saveRDS(r, file=paste0(<<<< FILEPATH REDACTED >>>>))

### Save Pernambuco estimates
Pernambuco_pop_1 <- 9557071

Pernambuco_draws <- (annual_preds/Pernambuco_pop_1)
Pernambuco_draws <- cbind("year"=predict_years, as.data.table(Pernambuco_draws))
write.table(Pernambuco_draws, file=paste0(<<<< FILEPATH REDACTED >>>>), row.names=FALSE, sep=",")

write.csv(data_BRA_merged[!is.na(lf_prev),], file=paste0(<<<< FILEPATH REDACTED >>>>))
