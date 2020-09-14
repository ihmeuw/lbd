####################################################################################################
## Description: Generate plots of model parameter estimates
##
## Passed args: main_dir [character] -- home directory for settings and final output
##
## Requires:    fitted model objects for both sexes ('model_fit_[sex].rdata' in temp_dir)
##              shape files (shape_file)
##
## Outputs:     plots of model parameter estimates ('model_fit.pdf' in main_dir)
##
## Run from within 'sae_central' directory!!
####################################################################################################

library(Matrix)
library(data.table)
library(maptools)
library(rgeos)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())

## Get settings ------------------------------------------------------------------------------------
main_dir <- commandArgs()[4]

source("settings.r")
get_settings(main_dir)

## Load fitted model objects and shapefiles --------------------------------------------------------
# model parameters
source("models/get_params.r")
mod <- rbindlist(lapply(sexes, function(sex) {
  mod <- get_params(paste0(temp_dir, "/model_fit_", sex, ".rdata"))
  mod[, sex := sex]
  mod
}))
mod[, lb := est - 2*se]
mod[, ub := est + 2*se]
mod[, sex := factor(sex, 2:1, c("Females", "Males"))]

# shape files
map <- get(load(shape_file))
num_j <- length(map)
map <- data.table(fortify(map))
map[, area := as.numeric(as.character(id))]
if (num_j != uniqueN(map$area)) stop('simplifying shape file caused a loss of areas')

## Plot hyperparameters ----------------------------------------------------------------------------
pdf(paste0(main_dir, "/model_fit.pdf"), width=14, height=8)

hyper <- mod[grepl("logit_rho|log_sigma", param),]
print(ggplot(hyper, aes(x=param, y=est, ymin=lb, ymax=ub, colour=sex)) +
  geom_pointrange(position=position_dodge(width=0.1), size=1) +
  labs(x="", y="Mean (+/- 2*SD)", title="Hyperparameters") +
  theme_bw(base_size=12) + theme(axis.text.x = element_text(angle=45, hjust=1)))
rm(hyper)

## Plot fixed effects ------------------------------------------------------------------------------
fe <- mod[param == "B",]
fe[, param := c("int", covars, covars_as), sex]
fe[, param := factor(param, levels=c("int", covars, covars_as))]
fe[, covar := as.numeric(param != "int")]
print(ggplot(fe, aes(x=param, y=est, ymin=lb, ymax=ub, colour=sex)) +
        facet_wrap(~ covar, scales="free") +
        geom_pointrange(position=position_dodge(width=0.1), size=1) +
        labs(x="", y="Mean (+/- 2*SD)", title="Fixed effects") +
        theme_bw(base_size=12) +
        theme(axis.text.x = element_text(angle=45, hjust=1),
              strip.background=element_blank(), strip.text=element_blank()))
rm(fe)

## Random age-time intercept -----------------------------------------------------------------------
if (model %in% c(1, 2, 3, 5, 6, 7)) { # in most models, this is a single (interacted) random effect
  re1 <- mod[param == "re1",]
  re1[, age := rep(ages, length(years)), by='sex']
  re1[, year := rep(years, each = length(ages)), by='sex']

  print(ggplot(re1, aes(x=year, y=est, colour=factor(age))) +
          facet_grid(sex ~ .) + geom_line() +
          labs(x="Year", y="Mean", title="Random effects: age-year intercept") +
          theme_bw(base_size=12))
  print(ggplot(re1, aes(x=age, y=est, colour=factor(year))) +
          facet_grid(sex ~ .) + geom_line() +
          labs(x="Age", y="Mean", title="Random effects: age-year intercept") +
          theme_bw(base_size=12))
  rm(re1)

} else if (model == 4) { # model 4 has two (non-interacted) random effects for age and time
  re1 <- mod[param == "re1",]
  re1[, age := ages]
  print(ggplot(re1, aes(x=age, y=est)) + facet_grid(sex ~ .) + geom_line() +
          labs(x="Age", y="Mean", title="Random effects: age intercept") +
          theme_bw(base_size=12))
  rm(re1)

  re2 <- mod[param == "re2",]
  re2[, year := years]
  print(ggplot(re2, aes(x=year, y=est)) + facet_grid(sex ~ .) + geom_line() +
          labs(x="Year", y="Mean", title="Random effects: year intercept") +
          theme_bw(base_size=12))
  rm(re2)

} else if (model %in% c("1b", "2b")) { # model 1b and 2b have only a time-effect
  re1 <- mod[param == "re1",]
  re1[, year := years]
  print(ggplot(re1, aes(x=year, y=est)) + facet_grid(sex ~ .) + geom_line() +
          labs(x="Year", y="Mean", title="Random effects: year intercept") +
          theme_bw(base_size=12))
  rm(re1)

}

## Random area intercept ---------------------------------------------------------------------------
re2 <- mod[param == "re2"]
if (model == 4) re2 <- mod[param == "re3"] # the random effect numbering is different for model 4
re2[, area := 1:num_j - 1L, sex]
ggplot(re2, aes(x=area, y=est)) + facet_grid(sex ~ .) + geom_point() +
  labs(x="Area", y="Mean", title="Random effects: area intercept") +
  theme_bw(base_size=12)

re2 <- merge(map, re2, by="area", allow.cartesian=T)
ggplot(re2) + facet_grid(~ sex) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=est)) +
  scale_fill_gradientn(colours=rev(brewer.pal(5, "Spectral")), name="RE", limits=re2[, max(abs(est))]*c(-1.01,1.01)) +
  scale_x_continuous("", breaks=NULL) +
  scale_y_continuous("", breaks=NULL) +
  coord_fixed(ratio=1) +
  labs(title="Random effects: area intercept") +
  theme_bw(base_size=12)
rm(re2)

## Random area-time slope --------------------------------------------------------------------------
if (model %in% c("1", "1b", "2", "2b", "6", "7")) {
  re3 <- mod[param == "re3"]
  re3[, area := 1:num_j - 1L, sex]
  print(ggplot(re3, aes(x=area, y=est)) + facet_grid(sex ~ .) + geom_point() +
          labs(x="Area", y="Mean", title="Random effects: area-time slope") +
          theme_bw(base_size=12))

  re3 <- merge(map, re3, by="area", allow.cartesian=T)
  print(ggplot(re3) + facet_grid(~ sex) +
          geom_polygon(aes(x=long, y=lat, group=group, fill=est)) +
          scale_fill_gradientn(colours=rev(brewer.pal(5, "Spectral")), name="RE", limits=re3[, max(abs(est))]*c(-1.01,1.01)) +
          scale_x_continuous("", breaks=NULL) +
          scale_y_continuous("", breaks=NULL) +
          coord_fixed(ratio=1) +
          labs(title="Random effects: area-time slope") +
          theme_bw(base_size=12))
  rm(re3)
}

## Random area-age slope ---------------------------------------------------------------------------
if (model %in% c("1", "2", "6")) {
  re4 <- mod[param == "re4"]
  re4[, area := 1:num_j - 1L, sex]
  print(ggplot(re4, aes(x=area, y=est)) + facet_grid(sex ~ .) + geom_point() +
          labs(x="Area", y="Mean", title="Random effects: area-age slope") +
          theme_bw(base_size=12))

  re4 <- merge(map, re4, by="area", allow.cartesian=T)
  print(ggplot(re4) + facet_grid(~ sex) +
          geom_polygon(aes(x=long, y=lat, group=group, fill=est)) +
          scale_fill_gradientn(colours=rev(brewer.pal(5, "Spectral")), name="RE", limits=re4[, max(abs(est))]*c(-1.01,1.01)) +
          scale_x_continuous("", breaks=NULL) +
          scale_y_continuous("", breaks=NULL) +
          coord_fixed(ratio=1) +
          labs(title="Random effects: area-age slope") +
          theme_bw(base_size=12))
  rm(re4)
}

## Random area-time intercept ----------------------------------------------------------------------
if (model %in% c("2b", 6, 7)) {
  re5 <- mod[param == "re5",]
  re5[, area := rep(1:num_j - 1L, length(years)), sex]
  re5[, year := rep(years, each=num_j), sex]
  print(ggplot(re5, aes(x=area, y=est, colour=sex)) +
          facet_wrap(~ year) + geom_point(size=0.3) +
          labs(x="Area", y="Mean", title="Random effects: area-year intercept") +
          guides(colour=guide_legend(override.aes=list(size=3))) +
          theme_bw(base_size=12))
  rm(re5)
}

## Random area-age intercept -----------------------------------------------------------------------
if (model == 6) {
  re6 <- mod[param == "re6",]
  re6[, area := rep(1:num_j - 1L, length(ages)), sex]
  re6[, age := rep(ages, each=num_j), sex]
  print(ggplot(re6, aes(x=area, y=est, colour=sex)) +
          facet_wrap(~ age) + geom_point(size=0.3) +
          labs(x="Area", y="Mean", title="Random effects: area-age intercept") +
          guides(colour=guide_legend(override.aes=list(size=3))) +
          theme_bw(base_size=12))
  rm(re6)
}

## Random area-time-age intercept ------------------------------------------------------------------
if (model == 2) {
  re5 <- mod[param == "re5",]
  re5[, area := rep(1:num_j - 1L, length(ages)*length(years)), sex]
  re5[, age := rep(rep(ages, each = num_j), length(years)), sex]
  re5[, year := rep(years, each = num_j * length(ages)), sex]
  for (this_year in years) {
    print(ggplot(re5[year == this_year,], aes(x=area, y=est)) +
            facet_grid(age ~ sex) + geom_point(size=0.5) +
            labs(x="Area", y="Mean", title=paste0("Random effects: area-age-year intercept (", this_year, ")")) +
            theme_bw(base_size=12))
  }
  rm(re5)
}

dev.off()
