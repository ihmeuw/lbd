cw_water <- function(mydat) {
	cw_dat <- fread('<<<< FILEPATH REDACTED >>>>') #wrtten in write_cw.R
  regions <- fread('<<<< FILEPATH REDACTED >>>>')
	
	cw_reg_dat <- cw_dat %>%
	  group_by(reg) %>%
	  summarize(well_imp = sum(well_imp),
	            well_unimp = sum(well_unimp),
	            spring_imp = sum(spring_imp),
	            spring_unimp = sum(spring_unimp),
	            piped = sum(piped),
	            piped_imp = sum(piped_imp))

	results <- list()
	for (i in unique(mydat$iso3)) {
		message(i)
		cw_sub <- subset(cw_dat, iso3 == i)
		cw_reg <- subset(cw_reg_dat, reg == subset(regions, iso3 == i)$region)
		if (nrow(cw_sub) == 0) {
			cw_sub[1,] <- 0
		}
		if (nrow(cw_reg) == 0) {
		  cw_reg[1,] <- 0
		}
		#piped ratio
	  if(cw_sub$sources < 5 | 
	      (cw_sub$piped == 0) |
	     (cw_sub$piped_imp == 0)){
  	  ipiped_pct <- cw_reg$piped_imp/(cw_reg$piped + cw_reg$piped_imp)
	  }else{
	    ipiped_pct <- cw_sub$piped_imp/(cw_sub$piped + cw_sub$piped_imp)
	  }
		#well ratio
		if(cw_sub$sources < 5 | 
		   (cw_sub$well_imp == 0) |
		   (cw_sub$well_unimp == 0)){
		  iwell_pct <- cw_reg$well_imp/(cw_reg$well_imp + cw_reg$well_unimp)
		}else{
		  iwell_pct <- cw_sub$well_imp/(cw_sub$well_imp + cw_sub$well_unimp)
		}
		#spring ratio
		if(cw_sub$sources < 5 | 
		   (cw_sub$spring_imp == 0) |
		   (cw_sub$spring_unimp == 0)){
		  ispring_pct <- cw_reg$spring_imp/(cw_reg$spring_imp + cw_reg$spring_unimp)
		}else{
		  ispring_pct <- cw_sub$spring_imp/(cw_sub$spring_imp + cw_sub$spring_unimp)
		}
		
		subdat <- mydat %>%
				 filter(iso3 == i)
		if(length(ipiped_pct) == 0){
		  ipiped_pct <- NA
		}
		if(length(iwell_pct) == 0){
		  iwell_pct <- NA
		}
		if(length(ispring_pct) == 0){
		  ispring_pct <- NA
		}

		if (is.na(ipiped_pct)) {
			ipiped_pct <- 1
			subdat$piped_cw <- 0
		}

		if (is.na(iwell_pct)) {
			iwell_pct <- 1
			subdat$well_cw <- 0
		}

		if (is.na(ispring_pct)) {
			ispring_pct <- 1
			subdat$spring_cw <- 0
		}

		subdat <- subdat %>%
				 mutate(network = piped + (ipiped_pct*piped_cw),
				   piped = piped + piped_cw + piped_imp,
				 		unimp = unimp + well_unimp + 
				 				(well_cw * (1 - iwell_pct)) +
				 				spring_unimp + (spring_cw * (1 - ispring_pct)),
				 		surface = surface) %>%
				 mutate(imp = well_imp + (well_cw * iwell_pct) +
				 			  spring_imp + (spring_cw * ispring_pct) +
				 			  imp) %>%
				 rename(N = total_hh) %>%
				 select(nid, iso3, survey_series, 
				 		lat, long, shapefile, location_code,
				 		year_start, int_year, year_median, sum_old_N,
				 		N, sum_of_sample_weights,
				 		network, piped, imp, unimp, surface, row_id)
		results[[length(results) + 1]] <- subdat

	}

	results <- do.call(rbind, results)

	return(results)
}

cw_sani <- function(mydat) {
  cw_dat <- fread('<<<< FILEPATH REDACTED >>>>')
  regions <- fread('<<<< FILEPATH REDACTED >>>>')
  
  cw_reg_dat <- cw_dat %>%
    group_by(reg) %>%
    summarize(latrine_imp = sum(latrine_imp),
              latrine_unimp = sum(latrine_unimp),
              flush_imp = sum(flush_imp),
              flush_unimp = sum(flush_unimp))
  
  results <- list()
  for (i in unique(mydat$iso3)) {
    message(i)
    cw_sub <- subset(cw_dat, iso3 == i)
    cw_reg <- subset(cw_reg_dat, reg == subset(regions, iso3 == i)$region)
    if (nrow(cw_sub) == 0) {
      cw_sub[1,] <- 0
    }
    if (nrow(cw_reg) == 0) {
      cw_reg[1,] <- 0
    }
    #flush ratio
    if(cw_sub$sources < 5 | 
       (cw_sub$flush_imp == 0) |
       (cw_sub$flush_unimp == 0)){
      iflush_pct <- cw_reg$flush_imp/(cw_reg$flush_imp + cw_reg$flush_unimp)
    }else{
      iflush_pct <- cw_sub$flush_imp/(cw_sub$flush_imp + cw_sub$flush_unimp)
    }
    #latrine ratio
    if(cw_sub$sources < 5 | 
       (cw_sub$latrine_imp == 0) |
       (cw_sub$latrine_unimp == 0)){
      ilatrine_pct <- cw_reg$latrine_imp/(cw_reg$latrine_imp + cw_reg$latrine_unimp)
    }else{
      ilatrine_pct <- cw_sub$latrine_imp/(cw_sub$latrine_imp + cw_sub$latrine_unimp)
    }
    
		subdat <- mydat %>%
				 filter(iso3 == i)
		if(length(ilatrine_pct) == 0){
		  ilatrine_pct <- NA
		}
		if(length(iflush_pct) == 0){
		  iflush_pct <- NA
		}

		if (is.na(ilatrine_pct)) {
			ilatrine_pct <- 1
			subdat$latrine_cw <- 0
		}

		if (is.na(iflush_pct)) {
			iflush_pct <- 1
			subdat$flush_cw <- 0
		}

		subdat <- subdat %>%
		  mutate(piped = (iflush_pct * flush_cw) + flush_imp + s_piped,
		         imp = imp + 
		           (ilatrine_pct * latrine_cw) + latrine_imp,
		         unimp = unimp + 
		           (latrine_cw * (1 - ilatrine_pct)) + latrine_unimp +
		           (flush_cw * (1 - iflush_pct)) + flush_unimp,
		         od = od) %>%
		  rename(N = total_hh) %>%
		  select(nid, iso3, survey_series, 
		         lat, long, shapefile, location_code,
		         year_start, int_year, year_median, sum_old_N,
		         N, sum_of_sample_weights, network, piped,
		         imp, unimp, od, row_id)
		results[[length(results) + 1]] <- subdat

	}

	results <- do.call(rbind, results)

	return(results)
}