clean_model_results_table <- function(rd   = run_date,
                                      regs = Regions,
                                      ages = 0,
                                      nm   = '',
                                      indic = indicator,
                                      ig = indicator_group,
                                      stackers = stacked_fixed_effects,
                                      coefs.sum1 = as.logical(coefs_sum1),
                                      tmb = fit_with_tmb){
  
  str_match <- stringr::str_match
  
  require(magrittr)
  
  sharedir <- <<<< FILEPATH REDACTED >>>>
  
  # make loopvars
  lv <- expand.grid(regs,ages)
  
  # grab formatted model fit objects
  stacker_names <- strsplit(stackers, " + ", fixed=T)[[1]]
  mods <- model_fit_table(lv=lv,rd=rd,nullmodel=nm, indicator = indic, indicator_group = ig,
                          coefs.sum1 = coefs.sum1, stacker_name_vec = stacker_names,
                          use_stacking_covs = use_stacking_covs, use_gp = use_gp,
                          spde_prior = spde_prior, fit_with_tmb = tmb)
  
  # add region column
  all_mods <- lapply(names(mods), function(n) {
    
    mod <- mods[[n]] %>% as.data.table
    
    r <- str_match(n,"(.*)_")[1,2]
    a <- str_match(n, "_(.*)")[1,2]
    mod[, region := r]
    mod[, age := a]
    return(mod)
  })
  
  all_mods <- rbindlist(all_mods)
  colorder <- c("region", "age",
                names(all_mods)[!(names(all_mods) %in% c("region", "age"))])
  setcolorder(all_mods, colorder)
  
  write.csv(all_mods, <<<< FILEPATH REDACTED >>>>,
            row.names = F)
  
  return(mods)
  
}
