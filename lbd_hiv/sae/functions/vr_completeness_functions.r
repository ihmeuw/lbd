## Function needed for running completeness prep code


#####################################################
## Grab most recent date from folder
###################################################

most_recent_date <- function(dir, date_format = "%Y_%m_%d", out_format = "%Y_%m_%d", file_pattern = NULL) {
  date_pattern <- gsub("%y|%m|%d", "[[:digit:]]{2}", gsub("%Y", "[[:digit:]]{4}", date_format))
  dates <- dir(dir, pattern = date_pattern)
  if (!is.null(file_pattern)) dates <- grep(file_pattern, dates, value = T)
  dates <- gsub(paste0("(.*)(", date_pattern, ")(.*)"), "\\2", dates)
  dates <- as.Date(dates, date_format)
  format(max(dates), out_format)
}

#############################################################
## This function calculates the smoothed number of recorded under 5 deaths by area
############################################################

completeness_shrink_u5m_deaths <- function(country, 
                                           years, 
                                           pop_file,
                                           shape,
                                           adjmat_file) {
  
  # Get VR under-5 deaths across all years by area
  source(paste0("<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/functions/vr_model_prep_functions.r"))
  vr_u5m <- vr_prep_all_cause_mortality(country, "adm2", years)
  vr_u5m <- vr_u5m[age < 5, .(vr_u5m = sum(deaths)), by = .(area)]
  
  # Make sure no areas have VR deaths less than 5
  if ((nrow(filter(vr_u5m, vr_u5m > 0)) != length(shape[["uid"]]))){
    message(paste0("There are ", nrow(filter(vr_u5m, vr_u5m == 0)),
                   " areas with 0 registered under 5 deaths out of ", length(shape[["link_id"]]), " (",
                   round(100*nrow(filter(vr_u5m, vr_u5m == 0))/length(shape[["link_id"]]), 1), "%)"))
    
    # Get population 
    load(pop_file)
    pop <- pop[age < 5, .(pop = sum(pop)), by = .(area)]
    pop_missing <- round(pop[area %in% vr_u5m[vr_u5m == 0, area], pop])
    
    message("Admin names are ", vr_pull_loc_meta(paste0(country, "_adm2")) %>% 
              filter(uid %in% vr_u5m[vr_u5m == 0, area]) %>% 
              pull(adm_name) %>% paste(collapse = ", "), " with populations over all years of ", 
            paste(pop_missing, collapse = ", "))
  }
  
  # Grab population over all years
  load(pop_file)
  pop <- pop[age < 5, .(pop = sum(pop)), by = .(area)]
  
  # Make a data table of unadjusted VR deaths for each area, asign ID to be used for shrinkage model
  unadjusted_vr <- left_join(vr_u5m, pop, by = "area") %>% mutate(ID = area + 1)
  
  # Grab adjacency matrix for model 
  adj_mat <- get(load(adjmat_file))
  
  # Fit shrinkage model using BYM2 spatial random effect, with PC priors assigned from Riebler et al.
  # The prior states the probability that the residuale relative risk is smaller than 2 is 0.99 prob
  # Note, this needs to run using old version of INLA that doesn't require rounding counts
  formula <- vr_u5m ~ 1 + f(ID, 
                            model = "bym2", 
                            graph = adj_mat,
                            scale.model = TRUE,
                            constr = TRUE,
                            hyper = list(
                              phi = list(prior = "pc",
                                         param = c(0.5, 0.5),
                                         initial = 1),
                              prec = list(prior = "pc.prec",
                                          param = c(1, 0.01),
                                          initial = 5))) 
  
  # Fit spatial model 
  mod <- inla(formula,
              family = "poisson", 
              data = unadjusted_vr,
              offset = log(pop))
  
  # Add fixed effects to random effect and get estimate of mortality rate for each area, only select area rate
  shrink_mod <- 
    data.table(mod$summary.random$ID[unadjusted_vr$ID,]) %>% 
    mutate(intercept = mod$summary.fixed[[1]]) %>% 
    mutate(rate = exp(intercept + mean)) %>% 
    mutate(area = ID - 1) %>% # Had to change index for INLA
    select(area, rate)
  
  # Join this to unadjusted VR to get new shrinkage estimates
  adjusted_vr <- 
    unadjusted_vr %>% 
    left_join(shrink_mod, by = "area") %>% 
    mutate(adjusted = rate*pop) %>%
    mutate(dif = adjusted / vr_u5m) %>% 
    select(area, vr_u5m, adjusted, dif, pop)
  
  # Return model and adjusted vr
  return(list(model = mod, adjusted_vr = adjusted_vr))
}


##########################################
##Function used to find K when adjusting completeness to match national completeness, modified from logit raking
#########################################

LogitFindK_mort <- function(gbdval, pixelval, weightval, MaxIter = 40, MaxJump = 20, FunTol = 1e-5, approx_0_1){
  
  # Logit raking won't work with any values of 0 or 1 in cell_pred
  # Adjust values slightly to avoid -Inf or Inf in NewPixelVal
  if (approx_0_1) {
    pixelval[pixelval == 0] <- 1e-10
    pixelval[pixelval == 1] <- 1-(1e-10)
  }
  
  NumIter <- ceiling(-log2(FunTol / MaxJump))
  NumIter <- 50
  
  if(NumIter > MaxIter){
    stop(paste("Maximum number of iterations is less than projected iterations required:", NumIter / MaxIter))
  }
  
  CurrentError <- EvalDiff_mort(gbdval, pixelval, weightval)
  if (CurrentError < 0){
    Range <- c(0, MaxJump)
  } else {
    Range <- c(-MaxJump, 0)
  }
  
  a <- Range[1]
  b <- Range[2]
  F_a <- EvalDiff_mort(gbdval, NewPixelVal(a, pixelval), weightval)
  F_b <- EvalDiff_mort(gbdval, NewPixelVal(b, pixelval), weightval)
  
  if (F_a * F_b > 0) {
    stop("Your estimates are WAY too far away from GBD")
  } else {
    i <- 1
    c <- (a + b) / 2
    F_c <- EvalDiff_mort(gbdval, NewPixelVal(c, pixelval), weightval)
    Success <- (abs(F_c) <= FunTol) 
    while (!Success & i < NumIter){
      if (sign(F_c) == sign(F_a)){
        a <- c
        F_a <- F_c
      } else {
        b <- c
        F_b <- F_c
      }
      c <- (a + b) / 2
      F_c <- EvalDiff_mort(gbdval, NewPixelVal(c, pixelval), weightval)
      Success <- (abs(F_c) <= FunTol)
      i <- i + 1
    }
    if (Success){
      return(c)
    } else {
      return(sprintf("Need more iterations, current output: K = %g, F(K) = %g",c, F_c))
    }
  }
}


NewPixelVal <- function(K, vals){
  return(ilogit(logit(vals)+K))
}


#vals: matrix
#weightval: vector
NewEst_mort <- function(vals, weightval){
  
  vals = (1 / vals) * weightval
  vals = apply(vals, 2, sum)
  
  return(mean(vals))
}


EvalDiff_mort <- function(gbdval, vals, weightval){
  return(gbdval - NewEst_mort(vals, weightval))
}


logit <- function(x) {
  log(x/(1-x))
}


ilogit <- function(x) {
  exp(x)/(1+exp(x))
}

#####################################################################
## Function for visualizing the completeness
#####################################################################

# Visualize prior completeness for adults and children
# Plotting tools for completeness

library(tidybayes, lib.loc = "<<<< FILEPATH REDACTED >>>>")
library(fitdistrplus)

plot_completeness_prior <- function(final_completeness, 
                                    prior_dist,
                                    save_file) {
  
  years = sort(unique(final_completeness$year))
  areas = sort(unique(final_completeness$area))
  
  gg_graphs_logit <- 
    lapply(areas, function(a) {
      message("working on area ", a)
      
      # Simulate from prior districution
      area_prior <-
        prior_dist[area == a, .(normal_prior = inv.logit(rnorm(n_draws, mean = mean, sd = sd)),
                                draw = 1:n_draws), by = c("area", "year")]
      
      # Merge on area specific completeness
      area_completeness <- merge(final_completeness[area == a], area_prior, by = c("area", "year", "draw"))
      
      # Make plot of draws as well as normal approximation
      int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
      
      gg_complete <- 
        area_completeness %>% 
        ggplot() +
        stat_halfeye(aes(x = year, y = final_complete), side = "left") +
        stat_slab(aes(x = year, y = normal_prior), side = "left", color = "red", fill = NA, linetype = "dotted") +
        labs(y = "completeness", title = paste("Completeness for area", a)) +
        scale_x_continuous(breaks = int_breaks) +
        theme_classic() + 
        coord_cartesian(ylim = c(0, 1.1)) 
    })
  
  # Make indeces for graphing these plots 
  total_areas <- length(areas)
  indices <- data.table(value1 = seq(1, total_areas, 10), value2 = c(seq(10, total_areas, 10), total_areas))
  
  # Make graph of priors, batched by 10 areas to one sheet
  pdf(file = save_file, height = 10, width = 18)
  apply(indices, 1, function(ind){
    rows = ceiling(length(ind[[1]]:ind[[2]]) / 2) 
    do.call(grid.arrange, gg_graphs_logit[ind[[1]]:ind[[2]]])
  })
  
  dev.off()
}

plot_dist_comparison <- function(final_completeness,
                                 prior_dist,
                                 save_file) { 
  
  years = sort(unique(final_completeness$year))
  areas = sort(unique(final_completeness$area))
  
  dist_diff <- 
    rbindlist(mclapply(areas, mc.cores = 10, function(a) {
      message("working on area ", a)
      
      # Simulate from prior districution
      area_prior <-
        prior_dist[area == a, .(normal_prior = inv.logit(rnorm(n_draws, mean = mean, sd = sd)),
                                draw = 1:n_draws), by = c("area", "year")]
      
      # Merge on area specific completeness
      area_completeness <- merge(final_completeness[area == a], area_prior, by = c("area", "year", "draw"))
      
      
      dist_diff <-  
        area_completeness[, .(mean = mean(final_complete),
                              mean_logit = mean(normal_prior),
                              lower = sort(final_complete)[25],
                              lower_logit = sort(normal_prior)[25],
                              upper = sort(final_complete)[975],
                              upper_logit = sort(normal_prior)[975]),
                          by = c("area", "year")]
    }))
  
  gg_mean <- 
    dist_diff %>% 
    mutate(logit_diff = round(100*(mean_logit - mean) / mean, 2)) %>% 
    ggplot() + 
    geom_histogram(aes(x = logit_diff, y = ..density..), position = "dodge", alpha = 0.3) + 
    geom_density(aes(x = logit_diff, y = ..density..)) + 
    theme_bw() + 
    labs(x = "Percent difference", title = "Percent difference between mean and distribution mean")
  
  gg_lower <- 
    dist_diff %>% 
    mutate(logit_diff = round(100*(lower_logit - lower) / lower, 2)) %>% 
    ggplot() + 
    geom_histogram(aes(x = logit_diff, y = ..density..), position = "dodge", alpha = 0.3) + 
    geom_density(aes(x = logit_diff, y = ..density..)) + 
    theme_bw() + 
    labs(x = "Percent difference", title = "Percent difference between lower UI bound and distribution lower UI bound")
  
  gg_upper <- 
    dist_diff %>% 
    mutate(logit_diff = round(100*(upper_logit - upper) / upper, 2)) %>%
    ggplot() + 
    geom_histogram(aes(x = logit_diff, y = ..density..), position = "dodge", alpha = 0.3) + 
    geom_density(aes(x = logit_diff, y = ..density..)) + 
    theme_bw() + 
    labs(x = "Percent difference", title = "Percent difference between upper UI bound and distribution upper UI bound")
  
  gg_ci <- 
    dist_diff %>% 
    mutate(ci_dif = upper - lower) %>% 
    mutate(ci_dif_logit = upper_logit - lower_logit) %>% 
    mutate(logit_diff = round(100*(ci_dif_logit - ci_dif) / ci_dif, 2)) %>% 
    ggplot() + 
    geom_histogram(aes(x = logit_diff, y = ..density..), position = "dodge", alpha = 0.3) + 
    geom_density(aes(x = logit_diff, y = ..density..)) + 
    theme_bw() + 
    labs(x = "percentage difference", title = "Percent difference in length of UI by distribution length of UI")
  
  pdf(file = save_file, height = 10, width = 18)
  print(gg_mean)
  print(gg_lower)
  print(gg_upper)
  print(gg_ci)
  dev.off()
}
