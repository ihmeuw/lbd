####################################################################################################
## Description: Generate descriptive plots of final model outputs
##
## Passed args: main_dir [character] -- home directory for settings and final output
##
## Requires:    Final model output ('est_all.rdata' in main_dir)
##              Shape files (shape_file)
##
## Outputs:     Plots of model predictions ('model_results.pdf' in main_dir)
##
## Run from within 'sae_central' directory!!
####################################################################################################

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)

rm(list=ls())

## Get settings ------------------------------------------------------------------------------------
main_dir <- commandArgs()[4]

source("settings.r")
get_settings(main_dir)

## Load final model output data and shapefiles -----------------------------------------------------
# Estimates
load(paste0(main_dir,"/est_all.rdata"))
est <- est[sex < 3,]
est[, sex_label := factor(sex, 2:1, c("Females", "Males"))]

# Covariates
if (!is.null(covars)) {
  load(covar_file)
  covar <- covar[, c(area_var, "year", covars), with=F]

  # applying covariate transformations
  if (!is.null(covars_trans)) {
    for (var in names(covars_trans)) {
      fun <- eval(parse(text = paste("function(x)", covars_trans[var])))
      covar[[var]] <- fun(covar[[var]])
      new_name <- paste0(gsub("\\(x\\)", "", covars_trans[var]), "_", var)
      setnames(covar, var, new_name)
      covars[covars == var] <- new_name
    }
  }
}

# Shape files
shp <- get(load(shape_file))
shp <- as.data.table(fortify(shp))
shp[, area := as.numeric(id)]

# Raw data
load(paste0(temp_dir, '/data.rdata'))
data[, year := year + min(years)]
data[, age := ages[age + 1]]
data <- rbind(data[, list(raw = sum(events)/sum(pop)), by='year,sex,age'],
              data[, list(age = 98, raw = sum(events)/sum(pop)), by='year,sex'])

## Initialize PDF ----------------------------------------------------------------------------------
pdf(paste0(main_dir, "/model_results.pdf"), width=14, height=8)

## Plot the crude and age-standardized results for the first, middle, and most recent year ---------
fdata <- est[CJ(area_var, unique(area), c(min(years), round(median(years)), max(years)), 1:2, 98),
             list(area, year, sex, sex_label, age, mx = 100000*mean), nomatch=0L]
color_values <- fdata[, quantile(mx, c(0, 0.005, 0.995, 1))]
color_values <- rescale(c(color_values[1], seq(color_values[2], color_values[3], length.out=7), color_values[4]))
fdata <- merge(fdata, shp, by="area", allow.cartesian=T)

ggplot(fdata) + facet_grid(sex_label ~ year) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=mx)) +
  scale_fill_gradientn(colours=brewer.pal(9, "OrRd"), values=color_values, name="Rate\n(per 100K)") +
  scale_x_continuous("", breaks=NULL) +
  scale_y_continuous("", breaks=NULL) +
  coord_fixed(ratio=1) +
  guides(fill = guide_colorbar(barheight=15, nbin=100)) +
  labs(title="Crude rate by sex") +
  theme_bw(base_size=12)

fdata <- est[CJ(area_var, unique(area), c(min(years), round(median(years)), max(years)), 1:2, 99),
             list(area, year, sex, sex_label, age, mx = 100000*mean), nomatch=0L]
color_values <- fdata[, quantile(mx, c(0, 0.005, 0.995, 1))]
color_values <- rescale(c(color_values[1], seq(color_values[2], color_values[3], length.out=7), color_values[4]))
fdata <- merge(fdata, shp, by="area", allow.cartesian=T)

ggplot(fdata) + facet_grid(sex_label ~ year) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=mx)) +
  scale_fill_gradientn(colours=brewer.pal(9, "OrRd"), values=color_values, name="Rate\n(per 100K)") +
  scale_x_continuous("", breaks=NULL) +
  scale_y_continuous("", breaks=NULL) +
  coord_fixed(ratio=1) +
  guides(fill = guide_colorbar(barheight=15, nbin=100)) +
  labs(title="Age-standardized rate by sex") +
  theme_bw(base_size=12)

## Plot a map of annualized rate of change by sex --------------------------------------------------
fdata <- est[CJ(area_var, unique(area), range(years), 1:2, 98:99),
             list(arc = 100*log(mean[2]/mean[1])/(year[2] - year[1])),
             by='area,sex,sex_label,age', nomatch=0L]
fdata <- merge(fdata, shp, by="area", allow.cartesian=T)

ggplot(fdata[age == 98,]) + facet_grid(~ sex_label) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=arc)) +
  scale_fill_gradientn(colours=rev(brewer.pal(5, "Spectral")), limits=max(abs(fdata$arc))*c(-1.01, 1.01), name="ARC (%)") +
  scale_x_continuous("", breaks=NULL) +
  scale_y_continuous("", breaks=NULL) +
  coord_fixed(ratio=1) +
  guides(fill = guide_colorbar(barheight=15, nbin=100)) +
  labs(title=paste0("Annualized rate of change in the crude rate by sex (", min(years), "-", max(years), ")")) +
  theme_bw(base_size=12)

ggplot(fdata[age == 99,]) + facet_grid(~ sex_label) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=arc)) +
  scale_fill_gradientn(colours=rev(brewer.pal(5, "Spectral")), limits=max(abs(fdata$arc))*c(-1.01, 1.01), name="ARC (%)") +
  scale_x_continuous("", breaks=NULL) +
  scale_y_continuous("", breaks=NULL) +
  coord_fixed(ratio=1) +
  guides(fill = guide_colorbar(barheight=15, nbin=100)) +
  labs(title=paste0("Annualized rate of change in the age-standardized rate by sex (", min(years), "-", max(years), ")")) +
  theme_bw(base_size=12)

## Plot predictions against the covariates ---------------------------------------------------------
if (!is.null(covars)) {
  fdata <- est[CJ(area_var, unique(area), c(min(year), round(median(year)), max(year)), 1:2, 98),
               list(area, year, sex, sex_label, mx = 100000*mean), nomatch=0L]
  fdata <- merge(fdata, covar, by.x=c("area", "year"), by.y=c(area_var, "year"), all.x=T)

  for (var in covars) {
    setnames(fdata, var, "var")
    p <- ggplot(fdata, aes(x=var, y=mx)) +
      facet_grid(sex_label ~ year) +
      geom_bin2d(bins=100) +
      geom_smooth(color="black", se=F) +
      scale_fill_distiller(palette="Oranges", direction=1) +
      guides(fill = guide_colorbar(barheight=15, nbin=100)) +
      labs(x=var, y="Age-standardized rate (per 100,000)", title=var) +
      theme_bw(base_size=12)
    print(p)
    setnames(fdata, "var", var)
  }
}

## Create line graphs showing national age-sex patterns --------------------------------------------
if ("natl" %in% names(geoagg_files)) {
  fdata <- merge(est[CJ("natl"), ], data, by=c("year", "sex", "age"), all.x=T)
  p <- ggplot(fdata[age < 98,]) + facet_wrap(~ year) +
    geom_ribbon(aes(x=age, ymin=100000*lb, ymax=100000*ub, fill=sex_label), alpha=0.25) +
    geom_line(aes(x=age, y=100000*mean, colour=sex_label)) +
    geom_point(aes(x=age, y=100000*raw, colour=sex_label)) +
    scale_color_discrete(name = "Sex") + scale_fill_discrete(name = "Sex") +
    labs(title="National rate by age, year, and sex",
         x="Age", y="Rate (per 100,000)") +
    theme_bw(base_size=12)
  print(p)

  p <- ggplot(fdata) + facet_wrap(~ age, scales="free") +
    geom_ribbon(aes(x=year, ymin=100000*lb, ymax=100000*ub, fill=sex_label), alpha=0.25) +
    geom_line(aes(x=year, y=100000*mean, colour=sex_label)) +
    geom_point(aes(x=year, y=100000*raw, colour=sex_label)) +
    scale_color_discrete(name = "Sex") + scale_fill_discrete(name = "Sex") +
    labs(title="National rate by year, age, and sex",
         x="Age", y="Rate (per 100,000)") +
    theme_bw(base_size=12)
  print(p)
}

## Finalize the PDF --------------------------------------------------------------------------------
dev.off()
