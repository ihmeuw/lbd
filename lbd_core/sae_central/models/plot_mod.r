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
library(gridExtra)
library(dplyr)
library(stringr)
library(gtools)

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
if (model %in% c(1, 2, 3, 5, 6, "6_c", 7)) { # in most models, this is a single (interacted) random effect
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

} else if (model %in% c(4, "4_c")) { # model 4 has two (non-interacted) random effects for age and time
  re1 <- mod[param == "re1",]
  re1[, age := ages, sex]
  print(ggplot(re1, aes(x=age, y=est)) + facet_grid(sex ~ .) + geom_line() +
          labs(x="Age", y="Mean", title="Random effects: age intercept") +
          theme_bw(base_size=12))
  rm(re1)

  re2 <- mod[param == "re2",]
  re2[, year := years, sex]
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
if (model %in% c(4, "4_c")) re2 <- mod[param == "re3"] # the random effect numbering is different for model 4
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
if (model %in% c("1", "1b", "2", "2b", "6", "6_c", "7")) {
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
if (model %in% c("1", "2", "6", "6_c")) {
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
if (model %in% c("2b", 6, "6_c", 7)) {
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
if (model %in% c(6, "6_c")) {
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


## Completeness posterior plots
if (grepl("_c", model)) {
  
  # Load some specific libraries
  library(tidyr)
  library(stringr)
  library(gridExtra)
  library(tidybayes, lib.loc = "<<<< FILEPATH REDACTED >>>>")
  
  ## First map in logit space -----------------------------------------------------------
  # Add mapping to age groups, area, year to posteriors
  
  load(paste0(temp_dir, "/data.rdata"))
  comp_map <- unique(data[, .(area, year, age, C)][, age_group := ifelse(age < 3, 0, 1)])[order(age_group, area, year)]
  
  if (admin1_completeness) {
    # Make sure areas are admin1
    load(geoagg_files[["admin1"]])
    area_admin1 <- unique(weights[, .(area, admin1)])
    comp_map <- unique(merge(comp_map, area_admin1)[, area := admin1][, .(area, year, age_group, C)])[order(age_group, area, year)]
  }
  
  pi_length <- length(unique(comp_map$C))
  # comp_map <- data.table(readRDS(completeness_prior_file))[, year := as.integer(year - years[1])][order(age_group, area, year)]
  # comp_map <- comp_map[, C := 1:nrow(comp_map) - 1L][, .(area, year, C, age_group)]
  # 
  # Make sure this aligns if completeness if fixed for adults or children
  if (!calc_child_comp) {
    comp_post <- mod[param == "logit_pi"][, C := seq(pi_length / 2 + 1, pi_length) - 1L, by = sex]
  } else if (!calc_adult_comp) {
    comp_post <- mod[param == "logit_pi"][, C := seq(1, pi_length / 2) - 1L, by = sex]
  } else {
    comp_post <- mod[param == "logit_pi"][, C := 1:pi_length - 1L, by = sex]
  }
  
  comp_post <- merge(comp_map, comp_post)
  comp_post[, year := year + years[1]]
  comp_post <- comp_post[, .(C, area, year, age_group, sex, est, lb, ub)]
  comp_post[, type := "post"]
  
  completeness_prior <- data.table(readRDS(completeness_prior_file))[order(age_group, area, year)] 
  completeness_prior[, C := 1:nrow(completeness_prior) - 1L]
  completeness_prior[, ub:= mean + 1.96*sd]
  completeness_prior[, lb := mean - 1.96*sd]
  setnames(completeness_prior, "mean", "est")
  completeness_prior <- completeness_prior[C %in% unique(comp_post$C), .(C, area, year, age_group, est, lb, ub)]
  completeness_prior <- rbind(copy(completeness_prior[, sex := "Males"]), copy(completeness_prior[, sex := "Females"]))
  completeness_prior[, type := "prior"]
  
  comp_post <- rbind(comp_post, completeness_prior)
  comp_post[, age_group := factor(age_group, levels = c(0, 1), labels = c("Child", "Adult"))]
  comp_post[, type := factor(type, levels = c("prior", "post"))]
  comp_post[, `:=`(est = inv.logit(est), lb = inv.logit(lb), ub = inv.logit(ub))]
  
  # Make plots
  int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
  upper <- ceiling(quantile(comp_post$ub, 0.99))
  lower <- floor(quantile(comp_post$lb, 0.01))
  areas <- unique(comp_post$area)
  gg_comp_post <- 
    lapply(areas, function(a) {
      comp_post[area == a] %>% 
        ggplot(aes(x = year, y = est, ymin = lb, ymax = ub, colour = type)) +
        geom_pointrange(position = position_dodge(width = 0.75)) +
        facet_grid(age_group ~ sex) +
        scale_x_continuous(breaks = int_breaks) +
        coord_cartesian(ylim = c(lower, upper)) + 
        labs(x="Year", y="Mean (+/- 2*SD)", title= paste("Model estimate of completeness in area", a)) +
        guides(colour=guide_legend(override.aes=list(size=1))) +
        theme_bw()
    })
  
  # Make graph of priors, batched by 4 areas to one sheet
  dir.create(paste0(main_dir, "/completeness_plots/"))
  pdf(paste0(main_dir, "/completeness_plots/prior_post.pdf"), width=14, height=8)
  indices <- data.table(value1 = seq(1, length(areas), 4), value2 = c(seq(4, length(areas), 4), length(areas)))
  apply(indices, 1, function(ind){
    rows = ceiling(length(ind[[1]]:ind[[2]]) / 2) 
    do.call(grid.arrange, gg_comp_post[ind[[1]]:ind[[2]]])
  })
  dev.off()
  
  ## Now map draws ----------------------------------------------------
  # Combine female and male draws
  pi_males   <- readRDS(paste0(main_dir, "/temp_dir/completeness_draws_1"))
  pi_males[, sex := 0]
  pi_females <- readRDS(paste0(main_dir, "/temp_dir/completeness_draws_2"))
  pi_females[, sex := 1]
  
  pi_post <-
    rbind(pi_males, pi_females) %>%
    pivot_longer(cols = paste0("V", 1:1000), names_to = "draw") %>%
    dplyr::mutate(draw = as.numeric(str_remove(draw, "V"))) %>%
    data.table()
  
  pi_post <- merge(pi_post, comp_map, by = "C")
  pi_post[, age_group := factor(age_group,
                                levels = c(0, 1),
                                labels = c("Child", "Adult"))]
  pi_post[, sex := factor(sex,
                          levels = c(0, 1),
                          labels = c("Males", "Females"))]
  pi_post[, year := year + years[1]]
    
  pi_post <- pi_post[, .(area, year, sex, age_group, draw, post = value)]
    
    # Grab the priors 
    # Eventually change to save priors to input 
  country = str_remove(main_dir, "<<<< FILEPATH REDACTED >>>>") %>% substr(1, 3)
  if (!calc_child_comp) {
    
    # For models of only adult completeness
    adult_prior <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
    adult_prior[, age_group := 1]
    pi_prior <- rbind(copy(adult_prior[, sex := 0]), copy(adult_prior[, sex := 1]))
    
  } else if (!calc_adult_comp) {
   
    # For models of only child completeness
    child_prior <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
    child_prior[, age_group := 0]
    pi_prior <- rbind(copy(child_prior[, sex := 0]), copy(child_prior[, sex := 1]))
    
  } else {
    
    # Running both adult and child completeness
    child_prior <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>/final_child_completeness.RDS"))
    child_prior[, age_group := 0]
    adult_prior <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>/final_adult_completeness.RDS"))
    adult_prior[, age_group := 1]
    # Combine child and adult priors, since we fit seperately by sex both have same prior
    pi_prior <- rbind(rbind(child_prior, adult_prior)[, sex := 0], rbind(child_prior, adult_prior)[, sex := 1])
  }
   
  pi_prior[, age_group := factor(age_group, levels = c(0, 1), labels = c("Child", "Adult"))]
  pi_prior[, sex := factor(sex, levels = c(0, 1), labels = c("Males", "Females"))]
  pi_prior <- pi_prior[, .(area, year, sex, age_group, draw, prior = final_complete)]
  
  # Merge prior and posterior
  final_pi <- merge(pi_prior, pi_post, by = c("area", "year", "sex", "age_group", "draw"))
  
  # Make plot of draws as well as normal approximation
  int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
  library(tidybayes, lib.loc = "<<<< FILEPATH REDACTED >>>>")
  
  gg_comp_post <- 
    lapply(unique(final_pi$area), function(a) {
      final_pi[area == a] %>% 
        ggplot() +
        stat_slab(aes(x = year, y = prior), side = "left", color = NA, fill = "grey80", size = 0.1, alpha = 0.7) +
        stat_slab(aes(x = year, y = post), side = "left", alpha = 1, fill = "grey10", color = "grey10", size = 0.5) +
        facet_grid(age_group ~ sex) + 
        lims(y = c(0, 1)) + 
        scale_x_continuous(breaks = int_breaks) +
        labs(x="Year", y="Posterior completeness", title= paste("Model posterior of completeness in area", a)) +
        guides(colour=guide_legend(override.aes=list(size=1))) +
        theme_bw()
    })
  
  # Make graph of priors and posteriors, batched by 4 areas to one sheet
  pdf(paste0(main_dir, "/completeness_plots/prior_posterior_draws.pdf"), width=14, height=8)
  indices <- data.table(value1 = seq(1, num_j, 4), value2 = c(seq(4, num_j, 4), num_j))
  apply(indices, 1, function(ind){
    rows = ceiling(length(ind[[1]]:ind[[2]]) / 2) 
    do.call(grid.arrange, gg_comp_post[ind[[1]]:ind[[2]]])
  })
  
  dev.off()
  
  # Make a scatter
  dist_pi <- 
    final_pi[, .(mean_prior = mean(prior),
                 mean_post  = mean(post),
                 ci_prior   = quantile(prior, 0.975) - quantile(prior, 0.025),
                 ci_post    = quantile(post, 0.975) - quantile(post, 0.025)), 
             by = c("area", "year", "sex", "age_group")]
  

  source('<<<< FILEPATH REDACTED >>>>/lbd_hiv/mbg/functions/crosswalk_functions.r')
    
  gg_mean  <- 
    dist_pi %>% 
    ggplot(aes(x = mean_prior, y = mean_post)) + 
    geom_point(alpha = 0.3, shape = 1) +
    geom_abline(slope = 1, intercept = 0) + 
    stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse = TRUE,
                     xpos = .01,
                     ypos = 0.99) +  
    facet_grid(age_group ~ sex) + 
    lims(y = c(0, 1), x = c(0, 1)) + 
    theme_bw() + 
    labs(x = "Prior mean", y = "Posterior mean")
  
  lim <- round(max(c(dist_pi$ci_post, dist_pi$ci_prior)), 1)
  gg_ci <- 
    dist_pi %>% 
    ggplot(aes(x = ci_prior, y = ci_post)) + 
    geom_point(alpha = 0.3, shape = 1) +
    geom_abline(slope = 1, intercept = 0) + 
    stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse = TRUE,
                     xpos = quantile(dist_pi$ci_prior, .01),
                     ypos = lim) +  
    facet_grid(age_group ~ sex) + 
    lims(y = c(0, lim), x = c(0, lim)) + 
    theme_bw() + 
    labs(x = "Prior CI length", y = "Posterior CI length")
  
  pdf(paste0(main_dir, "/completeness_plots/prior_posterior_mean_compare.pdf"), width=14, height=8)
  print(gg_mean)
  print(gg_ci)
  dev.off()
}

message("Done!")