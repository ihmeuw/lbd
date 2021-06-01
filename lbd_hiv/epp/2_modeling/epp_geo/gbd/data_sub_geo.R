################################################################################
## Purpose: These functions substitute data into the already formed GBD geography
##          epp model.  It is used to replace certain kinds of data that GBD
##          wants to be able to togle on and off.
##        Also note the 0 multiplier on migration. We are putting the migration
##        responsibilty on the popadust setting rather than a direct migration
##        model.  That is managed here.
################################################################################


### Code
#' @title set.anc.prior
#' @description helper function that sets the anc bias prior based on location.  Not directly used in our model.
#'
#' @param loc the location whose EPP data object is being
#' @param subpop the sub population whose ANC bias prior is getting set, ususally there is only one per loc, but not always.
#'
#' @return The same EPP data object as before but with updated anc bias priors.
#'
#' @export

set.anc.prior <- function(loc, subpop) {

	  gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
		table <- fread(anc.bias.file)
		if ((loc %in% unique(table$ihme_loc_id) & loc == subpop) | (grepl("IND", loc) & subpop %in% gen.pop.dict)) {
			ancbias.pr.mean <<- table[ihme_loc_id == loc, mean]
			ancbias.pr.sd <<- table[ihme_loc_id == loc, sd]
		# } else if (grepl("IND", loc) & subpop %in% gen.pop.dict) {
		# 	ancbias.pr.mean <<- 0.1821865
		# 	ancbias.pr.sd <<- 0.05083865
		} else {
			ancbias.pr.mean <<- 0.0
			ancbias.pr.sd <<- 0.001
		}
  # }
	print(paste0("ANC bias prior set to N(", ancbias.pr.mean, ", ", ancbias.pr.sd, ")"))
}


#' @title extend.trans.params
#' @description helper function that restructures the transition parameters in the EPP object to allow for the larger specification from GBD.
#'
#' @param dt An EPP data object whose transition parametrs will be expanded
#' @param start.year the first year of the simulation, should be set eslewhere
#' @param stop.year the last year of the simulation, should be set elsewhere
#'
#'
#' @return The same EPP data object as before but with expanded transition parameters.
#'
#' @export

extend.trans.params <- function(dt, start.year, stop.year) {
	for (n in names(dt)) {
		## create time series of on/off art mortality and cd4 progression
		n.years <- length(start.year:stop.year)
		# mortality
		single.year.cd4artmort <- attr(dt[[n]], "eppfp")$cd4artmort
		time.series.cd4artmort <- do.call("rbind", rep(list(single.year.cd4artmort), n.years))
		attr(dt[[n]], "eppfp")$cd4artmort <- time.series.cd4artmort
		#cd4 progression
		single.year.cd4prog <- attr(dt[[n]], "eppfp")$cd4prog
		time.series.cd4prog <- do.call("rbind", rep(list(single.year.cd4prog), n.years))
		attr(dt[[n]], "eppfp")$cd4prog <- time.series.cd4prog

		attr(dt[[n]], "eppfp")$mortyears <- n.years
		attr(dt[[n]], "eppfp")$cd4years <- n.years

	}
	return(dt)
}



#' @title sub.prev
#' @description helper function that substitutes out the hhs prevalence data in the EPP data object for hhs data managed centrally by GBD.
#' This is done to improve replicability and understanding of how out data inputs are impacting our model outputs.
#'
#' @param dt An EPP data object whose hhs prevalence data will be replaced
#' @param loc the location of the EPP model, used to find the hhs prevalence surveys
#'
#'
#' @return The same EPP data object as before but with new hhs data.
#'
#' @export


sub.prev <- function(loc, dt, likelihood_type, prev_date, prev_type){
    gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
    if (length(dt) == 1) {
        gen.pop.i <- 1
    } else {
        gen.pop.i <- which(names(dt) %in% gen.pop.dict)
    }
    code <- get_adm0_codes(loc) # yes you want gaul codes here because this is based on the direct output from a GAUL run, this needs to be updated
    surv.path <- ad0_prev_path
    data4 <- fread(surv.path)[ADM0_CODE == code]
    data4[,c("ADM0_CODE", "V1") := NULL]

    if (likelihood_type == "multi_norm") {
      load(paste0("<<<< FILEPATH REDACTED >>>>"))
      attr(dt[[gen.pop.i]], 'likdat')$means <- means
      attr(dt[[gen.pop.i]], 'likdat')$cov_mat <- mat

    }

    if (nrow(data4) > 0) {
        data4[,used := TRUE]
        data4[prev == 0,used := FALSE]
        data4[,W.hhs := qnorm(prev)]
        data4[,v.hhs := 2*pi*exp(W.hhs^2)*se^2]
        data4[,sd.W.hhs := sqrt(v.hhs)]
        data4[,idx := year - (start.year - 1)]

        ## logit version
        data4[,W.logit.hhs := logit(prev)]
        data4[,sd.logit.hhs := (se)/(prev*(1 - prev))]



        data4 <- data4[order(data4$year),]
        attr(dt[[gen.pop.i]], 'eppd')$hhs <- as.data.frame(data4[, 1:5])
        if (all(!attr(dt[[gen.pop.i]], "eppd")$anc.used)) {
            attr(dt[[gen.pop.i]], 'likdat')$hhslik.dat <- data4
            attr(dt[[gen.pop.i]], 'likdat')$lastdata.idx <- max(attr(dt[[gen.pop.i]], 'likdat')$hhslik.dat$idx)
            attr(dt[[gen.pop.i]], 'likdat')$firstdata.idx <- min(attr(dt[[gen.pop.i]], 'likdat')$hhslik.dat$idx)
        } else {
            attr(dt[[gen.pop.i]], 'likdat') <- fnCreateLikDat(attr(dt[[gen.pop.i]], 'eppd'), floor(attr(dt[[gen.pop.i]], "eppfp")$proj.steps[1]), no.anc = no.anc)
        }
    } else {
        print(paste0("No surveys for ",loc))
    }
    return(dt)
}

#' @title sub.anc
#' @description helper function that adds backcast anc data into the EPP data object for use in the GBD models.  The LBD models do not use
#' backcast anc data so this function should not get called in the LBD model at all.  I am inducing a warning whenever this gets called to
#' hopefully help know what is happening.
#'
#' @param dt An EPP data object whose anc prevalence data will be replaced
#' @param loc the location of the EPP model, used to find the anc data
#'
#'
#' @return The same EPP data object as before but with new anc data.
#'
#' @export


sub.anc <- function(loc, dt) {
# Make adjustments to ANC coming from PJNZ files ** add more **
    ## Prep EPP data
  warning("The sub.anc function got called and backcast anc data is being added to your model")
    # Choose subpopulation for substitution
    gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
    if (length(names(dt)) > 1) {
    	gen.pop <- names(dt)[names(dt) %in% gen.pop.dict]
    } else {
    	gen.pop <- 1
    }
    eppd <- attr(dt[[gen.pop]], "eppd")
	if (grepl("ZAF", loc) | grepl("SWZ", loc)) {
	# Collapse up to single provincial ANC site
		# Extract first year of data and use that site as provincial site
		first.year <- min(as.integer(colnames(eppd$anc.prev)[sapply(colnames(eppd$anc.prev), function(col) {
			!all(is.na(eppd$anc.prev[, col]))
		})]))
		prov.sites <- which(!is.na(eppd$anc.prev[, as.character(first.year)]))
		for (i in 1:length(prov.sites)) {
			prov.site <- prov.sites[i]
			row.lower <- prov.sites[i] + 1
			row.upper <- ifelse(i == length(prov.sites), nrow(eppd$anc.prev), prov.sites[i + 1] - 1)
			eppd$anc.used[row.lower:row.upper] <- F
			# Sum administrative units to provincial level
			site.prev <- eppd$anc.prev[row.lower:row.upper,]
			site.n <- eppd$anc.n[row.lower:row.upper,]
			site.pos <- site.prev * site.n
			sub.pos <- colSums(site.pos, na.rm = T)
			sub.n <- colSums(site.n, na.rm = T)
			sub.prev <- sub.pos / sub.n
			# Append to provincial site
			for (c in colnames(eppd$anc.prev)) {
				if (is.na(eppd$anc.prev[prov.site, c])) {
					eppd$anc.prev[prov.site, c] <- sub.prev[c]
					eppd$anc.n[prov.site, c] <- sub.n[c]
				}
			}
		}
	} else {
	# Add imputed data
		# Read ANC data from back cast
		anc.dir <- "<<<< FILEPATH REDACTED >>>>"
		recent <- "180529"#max(as.integer(list.files(anc.dir)))
		anc.path <- paste0(anc.dir, recent, "/data/", loc, ".csv")
		anc.dt <- fread(anc.path)
		anc.dt[, clinic := gsub("[^[:alnum:] ]", "",clinic)] # For differences in naming like added special characters
		anc.dt <- anc.dt[order(clinic)]
		# Add draw level data from ANC backcast
	    for (cl in unique(anc.dt$clinic)) {
			clinic.idx <- which(grepl(gsub(" ", "", cl),gsub(" ", "",  rownames(eppd$anc.prev))))
			sub.dt <- anc.dt[clinic == cl]
			sub.dt[pred == "Data", (paste0("draw_", i)) := mean]
			for (y in unique(anc.dt[pred == "Data"]$year_id)) {
				eppd$anc.prev[clinic.idx, as.character(y)] <- sub.dt[year_id == y, get(paste0("draw_", i))]
				eppd$anc.n[clinic.idx, as.character(y)] <- sub.dt[year_id == y, n]
			}
		}
	}

	# Reformat EPP object with updated data
    attr(dt[[gen.pop]], "eppd") <- eppd

    set.list.attr <- function(obj, attrib, value.lst)
    mapply(function(set, value){ attributes(set)[[attrib]] <- value; set}, obj, value.lst)
 	attr(dt[[gen.pop]], "likdat") <- fnCreateLikDat(eppd, anchor.year = floor(attr(dt[[gen.pop]], "eppfp")$proj.steps[1]), no.anc = no.anc)

    return(dt)
}


#' @title sub.pop.params
#' @description function that swaps the GBD population parameters into the EPP data object.  Quick note, this is where migration is set to 0 so
#' that popadjust handels the migration instead of our modeled migration.
#'
#' @param epp.subp a part of an EPP data object.  This function is called inside prep_epp_data and that function ids this particular data structure
#' @param loc the location of the EPP model, used to find the demography data needed
#' @param mig a boolian indicating if GBD migration "TRUE" or SPECTRUM migration "FALSE" is to be used.  This is kind of vestigal becasue all
#' migration gets set to 0 for logical consistency with the admin 2 models where that greatly improved performance in small areas.
#'
#'
#' @return The same EPP data object as before but with new population and migration data.
#'
#' @export

sub.pop.params <- function(epp.subp, loc, mig = FALSE) {
  	## Load central functions
    loc.id <- loc.table[ihme_loc_id == loc, location_id]
    years <- epp.subp[[1]]$year
	  path <- paste0("<<<< FILEPATH REDACTED >>>>")
	  if (file.exists(path)) {
		  in.pop <- fread(path)[year_id %in% years]
		  in.pop <- in.pop[which(in.pop$age_group_id != 15), ]
	  } else {
		if (!"get_population" %in% ls()) {
   			source(paste0("<<<< FILEPATH REDACTED >>>>/get_population.R"))
    	}
	}

	# add in missing years in the future
	max.pop <- copy(in.pop[year_id == max(year_id)])
	missing.dt <- rbindlist(lapply(setdiff(years, unique(in.pop$year_id)), function(year) {
		dt <- copy(max.pop)
		dt[, year_id := year]
	}))
	bound.pop <- rbind(in.pop, missing.dt)
	both.pop <- bound.pop[, .(population = sum(population)), by = .(age_group_id, year_id)]
    pop15 <- both.pop[age_group_id == 8, population] / 5
    pop50 <- both.pop[age_group_id == 14, population] / 5
    pop15to49 <- both.pop[, .(population = sum(population)), by = .(year_id)]$population

    if (mig == TRUE) {
      nat_mig <- fread(nat_mig_path)
      nat_mig <- nat_mig[which(nat_mig$ihme_loc_id == loc), ]
      # add in missing years in the future
      max.mig <- copy(nat_mig[year_id == max(year_id)])
      missing.dt <- rbindlist(lapply(setdiff(years, unique(nat_mig$year_id)), function(year) {
        dt <- copy(max.mig)
        dt[, year_id := year]
      }))
      bound.mig <- rbind(nat_mig, missing.dt)
      names(bound.mig) <- c("V1","year_id", "ihme_loc_id", "mig")
      netmigr <- bound.mig$mig


      for (pop in names(epp.subp)) {
        temp.dt <- epp.subp[[pop]]
        if (pop == "subpops") {
          for (subpop in names(temp.dt)) {
            epp.subp[[pop]][[subpop]]$pop15to49 <- pop15to49
            epp.subp[[pop]][[subpop]]$pop15 <- pop15
            epp.subp[[pop]][[subpop]]$pop50 <- pop50
            epp.subp[[pop]][[subpop]]$netmigr <- (0 * netmigr) ####################################### this allows for a scalar to be applied to the migration, that scalar should be 0 for a final run becasue we are putting the migration responsibility on the popadjust setting, my migration model was junk
          }
        } else {
          epp.subp[[pop]]$pop15to49 <- pop15to49
          epp.subp[[pop]]$pop15 <- pop15
          epp.subp[[pop]]$pop50 <- pop50
          epp.subp[[pop]]$netmigr <- (0 * netmigr) ####################################### this allows for a scalar to be applied to the migration, that scalar should be 0 for a final run becasue we are putting the migration responsibility on the popadjust setting, my migration model was junk
        }
      }
    } else {

    for (pop in names(epp.subp)) {
    	temp.dt <- epp.subp[[pop]]
        if (pop == "subpops") {
            for (subpop in names(temp.dt)) {
            	epp.subp[[pop]][[subpop]]$pop15to49 <- pop15to49
            	epp.subp[[pop]][[subpop]]$pop15 <- pop15
            	epp.subp[[pop]][[subpop]]$pop50 <- pop50

            }
        } else {
            	epp.subp[[pop]]$pop15to49 <- pop15to49
            	epp.subp[[pop]]$pop15 <- pop15
            	epp.subp[[pop]]$pop50 <- pop50

        }
    }
    }
    return(epp.subp)
}


#' @title calc.expand.pop
#' @description function that expands the GBD population numbers using Incidence Rate Ratios from literature to determine what the age and sex
#' breakdown of the HIV+ population is so we can get the best estimates of HIV+ mortality parameters and CD4 Progression parameters.
#'
#' @param loc the location of the EPP model, used to find the demography data needed
#' @param sex.agg a boolian indicating if exes are to be aggregated.
#'
#'
#' @return an expanded version of the population.
#'
#' @export

calc.expand.pop <- function(loc, sex.agg = T, test_IRR = F) {

	loc.id <- loc.table[ihme_loc_id == loc, location_id]

	## Load central function
	age.table <- fread("<<<< FILEPATH REDACTED >>>>")
	IRR2 <- fread(age_irr_path)
	IRR2 <- IRR2[age < 55,]

	sex_IRR <- fread(sex_irr_path)
	sex_IRR <- sex_IRR[epidemic_class == "GEN",]
	sex_IRR[,year := year + start.year - 1]

	missing_years <- c()
	if (sex_IRR[,max(year)] < stop.year)
	  missing_years <- (sex_IRR[,max(year)] + 1):stop.year
	replace_IRR <- sex_IRR[order(year)][rep(nrow(sex_IRR), times = length(missing_years))]
	if (length(missing_years) > 0)
	  replace_IRR[,year := missing_years]
	sex_IRR <- rbind(sex_IRR, replace_IRR)

	sex_IRR[,sex := 2]

	male_IRR <- copy(sex_IRR)
	male_IRR[,FtoM_inc_ratio := 1.0]
	male_IRR[,sex := 1]

	sex_IRR <- rbind(sex_IRR, male_IRR)

	## Read in population for age-sex structure for aggregation
	path <- paste0("<<<< FILEPATH REDACTED >>>>")
	if (file.exists(path)) {
		in.pop <- fread(path)[year_id %in% start.year:stop.year]
	} else {
		if (!"get_population" %in% ls()) {
			source(paste0("<<<< FILEPATH REDACTED >>>>/get_population.R"))
		}
	}
	pop <- merge(in.pop, age.table[, .(age_group_id, age_group_name_short)], by = "age_group_id", all.x = T)
	setnames(pop,
	c("year_id", "sex_id", "age_group_name_short", "population"),
	c("year", "sex", "age", "value")
	)
	pop[, (setdiff(names(pop), c("year", "sex", "age", "value"))) := NULL]


	pop$age <- strtoi(pop$age)
	pop[(age - 5) %%  10 != 0, age := as.integer(age - 5)]
	pop[,value := as.numeric(value)]

	pop1 <- data.table(aggregate(value ~ sex + age + year,pop,FUN = sum))[order(sex,age)]
	missing_years <- c()
	if (pop1[,max(year)] < stop.year)
	  missing_years <- (pop1[,max(year)] + 1):stop.year
	replace_pop <- pop1[rep(which(pop1[,year] == pop1[,max(year)]), times = length(missing_years))]
	replace_years <- rep(1:(stop.year - pop1[,max(year)]), each = length(which(pop1[,year] == pop1[,max(year)])))
	replace_pop[,year := year + replace_years]
	pop1 <- rbind(pop1, replace_pop)
	IRR <- runif(16, IRR2$lower, IRR2$upper)
	IRR2[,IRR := IRR]

	IRR2[,IRR := IRR2[,IRR]/IRR2[age == 25 & sex == 1,IRR]]

	combined_IRR <- merge(sex_IRR, IRR2, by = 'sex', allow.cartesian = TRUE)
	combined_IRR[,comb_IRR := FtoM_inc_ratio * IRR]

	if (test_IRR == F) {pop2 <- merge(pop1, combined_IRR, by = c('sex', 'age', 'year'))}
	if (test_IRR == T) {
	  test_IRRs <- read.csv("<<<< FILEPATH REDACTED >>>>")
	  test_IRRs$loc <- as.character(test_IRRs$loc)
	  test_IRRs <- test_IRRs[which(test_IRRs$loc == loc), ]
	  test_IRRs$comb_IRR <- test_IRRs$IRR
	  test_IRRs$IRR <- NULL
	  test_IRRs$age <- (test_IRRs$age_group_id - 5)*5
	  test_IRRs$comb_IRR[which(test_IRRs$mean == 0)] <- 1 #fix the 0/0 issue when there is no incience.
	  test_IRRs <- as.data.table(test_IRRs)
	  pop2 <- merge(pop1, test_IRRs, by.x = c('sex', 'age', 'year'), by.y = c("sex_id", "age", "year_id"))
	}

	pop2[,wt := comb_IRR*value]

	sex_agg <- pop2[,.(wt = sum(wt)),by = .(year, age)]

	total <- pop2[,.(total = sum(wt)),by = .(year)]
	pop2 <- merge(pop2, total, by = c('year'))
	pop2[,ratio := wt/total]


	sex_agg <- merge(sex_agg, total, by = c('year'))
	sex_agg[,ratio := wt/total]

	if (sex.agg) {
		out.pop <- sex_agg
	} else {
		out.pop <- pop2
	}

	return(out.pop)
}

#' @title sub.on.art
#' @description function that substitutes the on art HIV + mortality parameters after weighting by our best estimate of the age and sex structure
#' of the HIV+ population in an area.
#'
#' @param dt an EPP data object that will have its HIV+ on art mortality parameters replaced
#' @param loc the location of the EPP model, used to find the apropriate mortality parameters
#' @param k indicates what draw this is.
#'
#'
#' @return The same EPP data object as before but with new on art HIV+ mortality parameters.
#'
#' @export

sub.on.art <- function(dt, loc, k, mean_test = FALSE, tirr = F) {

  mortart <- fread(paste0(mortart_path, loc,"_HIVonART.csv"))
	mortart <- melt(mortart,
	                id = c("durationart", "cd4_category", "age", "sex","cd4_lower",
	                       "cd4_upper"))
	setnames(mortart, c("variable","value"), c("drawnum","draw"))
	mortart <- mortart[,drawnum := substr(drawnum, 5,8)]
	mortart <- mortart[order(durationart,cd4_category,age,sex,drawnum)]
	mortart <- mortart[,c("durationart", "cd4_category", "age", "sex","cd4_lower",
	                      "cd4_upper", "drawnum", "draw"), with = F]
	mortart_read <- data.table(dcast(mortart,durationart+cd4_category+age+sex~drawnum, value.var = 'draw'))
	for (i in 1:1000) {
	  j <- i + 4
	  setnames(mortart_read, j, paste0("draw",i))
	}
	overs <- names(mortart_read)[5:1004]
	mortart_read[,drawmean := rowMeans(mortart_read[,overs, with = F])]
	mortart_read <- melt(mortart_read, id = c("durationart", "cd4_category", "age", "sex"))
	setnames(mortart_read, c("variable","value","cd4_category"),c("draw","mort","cd4"))
	mortart_read <- mortart_read[age != "55-100",]
	mortart_read <- mortart_read[,draw := substr(draw,5,8)]
	mortart_read$draw <- as.numeric(mortart_read$draw)


	mortart <- mortart_read[draw == k,]
	if (mean_test == TRUE) {mortart <- mortart_read[is.na(draw),] }
	mortart[,age := as.integer(sapply(strsplit(mortart[,age],'-'), function(x) {x[1]}))]
	mortart[,sex := as.integer(sex)]
	cd4_cats <- unique(mortart[,cd4])
	durat_cats <- unique(mortart[,durationart])
	cd4_vars <- expand.grid(durationart = durat_cats, cd4 = cd4_cats)
	expanded_pop <- calc.expand.pop(loc, sex.agg = F, test_IRR = tirr)
	n <- nrow(expanded_pop)
	expanded_pop <- expanded_pop[rep(1:n, times = length(cd4_cats)*length(durat_cats))]
	expanded_pop <- expanded_pop[order(year, sex, age)]
	expanded_pop <- cbind(expanded_pop, cd4_vars[rep(1:(nrow(cd4_vars)), times = n),])
	combined_mort <- merge(expanded_pop, mortart, by = c('durationart', 'cd4', 'sex', 'age'))
	mortart <- combined_mort[,.(mort = sum(ratio*mort)), by = .(durationart, cd4, year)]

	mortart <- mortart[cd4 == "ARTGT500CD4", cat := 1]
	mortart <- mortart[cd4 == "ART350to500CD4", cat := 2]
	mortart <- mortart[cd4 == "ART250to349CD4", cat := 3]
	mortart <- mortart[cd4 == "ART200to249CD4", cat := 4]
	mortart <- mortart[cd4 == "ART100to199CD4", cat := 5]
	mortart <- mortart[cd4 == "ART50to99CD4", cat := 6]
	mortart <- mortart[cd4 == "ARTLT50CD4", cat := 7]
	mortart[,risk := -1*log(1 - mort)]
	mortart <- mortart[,c("risk","cat","durationart","year"), with = F]
	mortart <- mortart[, setattr(as.list(risk), 'names', cat), by = c("year","durationart")]
	mortart <- mortart[order(year, durationart)]
	# mortart <- mortart[age!="55-100",]

	mortart1 <- mortart[durationart == "LT6Mo",]
	mortart2 <- mortart[durationart == "6to12Mo",]
	mortart3 <- mortart[durationart == "GT12Mo",]

	alpha1gbd <- as.matrix(data.frame(mortart1[,c("1","2","3","4","5","6", "7"), with = F]))
	alpha2gbd <- as.matrix(data.frame(mortart2[,c("1","2","3","4","5","6", "7"), with = F]))
	alpha3gbd <- as.matrix(data.frame(mortart3[,c("1","2","3","4","5","6", "7"), with = F]))

	alpha1 <- as.vector(t(alpha1gbd))
	alpha2 <- as.vector(t(alpha2gbd))
	alpha3 <- as.vector(t(alpha3gbd))
	for (n in names(dt)) {
	    attr(dt[[n]], 'eppfp')$cd4artmort[,1] <- alpha1
	    attr(dt[[n]], 'eppfp')$cd4artmort[,2] <- alpha2
	    attr(dt[[n]], 'eppfp')$cd4artmort[,3] <- alpha3

	}
	return(dt)
}

#' @title sub.off.art
#' @description function that substitutes the off art HIV + mortality parameters after weighting by our best estimate of the age and sex structure
#' of the HIV+ population in an area.
#'
#' @param dt an EPP data object that will have its HIV+ off art mortality parameters replaced
#' @param loc the location of the EPP model, used to find the apropriate mortality parameters
#' @param k indicates what draw this is.
#'
#'
#' @return The same EPP data object as before but with new off art HIV+ mortality parameters.
#'
#' @export

sub.off.art <- function(dt, loc, k, mean_test = FALSE, tirr = F) {

  # Off-ART Mortality
	mortnoart <- fread(paste0(mortnoart_path,loc,"_mortality_par_draws.csv"))
	mortnoart[,draw := rank(-mort,ties.method = "first"),by = c("age","cd4")]
	mortnoart <- mortnoart[order(age,cd4,draw)]
	mortnoart_read <- mortnoart[,c("age","cd4","draw","mort"), with = F]
	mortnoart <- mortnoart_read[draw == k,]
	if (mean_test == TRUE) {
	  temp <- mortnoart_read[, lapply(c("mort"), function(x) mean(get(x), na.rm = T)), by = c("age", "cd4")]
	  temp$draw <- "m"
	  temp$mort <- temp$V1
	  temp$V1 <- NULL
	  mortnoart <- temp
	}

	mortnoart[,age := as.integer(sapply(strsplit(mortnoart[,age],'-'), function(x) {x[1]}))]
	mortnoart[,risk := -1*log(1 - mort)/0.1]
	mortnoart[,prob := 1 - exp(-1*risk)]

	cd4_cats <- unique(mortnoart[,cd4])
	cd4_vars <- data.table(cd4 = cd4_cats)
	sex_agg <- calc.expand.pop(loc, sex.agg = T, test_IRR = tirr)
	expanded_pop <- sex_agg[rep(1:nrow(sex_agg), times = length(cd4_cats))]
	expanded_pop <- expanded_pop[order(year, age)]
	expanded_pop <- cbind(expanded_pop, cd4_vars[rep(1:(nrow(cd4_vars)), times = nrow(sex_agg)),])
	combined_mu <- merge(expanded_pop, mortnoart, by = c('cd4', 'age'))
	mortnoart <- combined_mu[,.(prob = sum(ratio*prob)), by = .(cd4, year)]

	mortnoart <- mortnoart[cd4 == "GT500CD4", cat := 1]
	mortnoart <- mortnoart[cd4 == "350to500CD4", cat := 2]
	mortnoart <- mortnoart[cd4 == "250to349CD4", cat := 3]
	mortnoart <- mortnoart[cd4 == "200to249CD4", cat := 4]
	mortnoart <- mortnoart[cd4 == "100to199CD4", cat := 5]
	mortnoart <- mortnoart[cd4 == "50to99CD4", cat := 6]
	mortnoart <- mortnoart[cd4 == "LT50CD4", cat := 7]
	mortnoart[,risk := -1*log(1 - prob)]
	mortnoart <- mortnoart[,.(year,risk,cat)]
	mortnoart <- mortnoart[, setattr(as.list(risk), 'names', cat), by = .(year)]
	mortnoart <- mortnoart[order(year)]
	mortnoart <- mortnoart[,c("1","2","3","4","5","6", "7"), with = F]
	mortnoart <- data.frame(mortnoart)
	mugbd <- as.matrix(mortnoart)
	mu <- as.vector(t(mugbd))
	for (n in names(dt)) {
        attr(dt[[n]], 'eppfp')$cd4artmort[,0] <- mu
    }
	return(dt)
}

#' @title sub.cd4.prog
#' @description function that substitutes the CD4 progression parameters after weighting by our best estimate of the age and sex structure
#' of the HIV+ population in an area.
#'
#' @param dt an EPP data object that will have its CD4 progression parameters replaced
#' @param loc the location of the EPP model, used to find the apropriate parameters
#' @param k indicates what draw this is.
#'
#'
#' @return The same EPP data object as before but with new CD4 progression parameters.
#'
#' @export



sub.cd4.prog <- function(dt, loc, k, mean_test = FALSE, tirr = F) {
	progdata <- fread(paste0(progdata_path, loc, "_progression_par_draws.csv"))
	progdata <- progdata[order(age,cd4,draw)]
	progdata_read <- progdata[,c("age","cd4","draw","prog"), with = F]
	progdata_read <- progdata_read[,lambda := 1/prog]

	progdata <- progdata_read[draw == k,]

	if (mean_test == TRUE) {
	  temp <- progdata_read[, lapply(c("prog", "lambda"), function(x) mean(get(x), na.rm = T)), by = c("age", "cd4")]
	  temp$draw <- "m"
	  temp$prog <- temp$V1
	  temp$V1 <- NULL
	  temp$lambda <- temp$V2
	  temp$V2 <- NULL
	  progdata <- temp
	}


	progdata[,risk := -1*log(1 - prog)/0.1]
	progdata[,prob := 1 - exp(-1*risk)]
	progdata[,age := as.integer(sapply(strsplit(progdata[,age],'-'), function(x) {x[1]}))]

	cd4_cats <- unique(progdata[,cd4])
	cd4_vars <- data.table(cd4 = cd4_cats)

	sex_agg <- calc.expand.pop(loc, sex.agg = T, test_IRR = tirr)
	expanded_pop <- sex_agg[rep(1:nrow(sex_agg), times = length(cd4_cats))]
	expanded_pop <- expanded_pop[order(year, age)]
	expanded_pop <- cbind(expanded_pop, cd4_vars[rep(1:(nrow(cd4_vars)), times = nrow(sex_agg)),])
	combined_prog <- merge(expanded_pop, progdata, by = c('cd4', 'age'))
	progdata <- combined_prog[,.(prob = sum(ratio*prob)), by = .(cd4, year)]

	progdata <- progdata[cd4 == "GT500CD4", cat := 1]
	progdata <- progdata[cd4 == "350to500CD4", cat := 2]
	progdata <- progdata[cd4 == "250to349CD4", cat := 3]
	progdata <- progdata[cd4 == "200to249CD4", cat := 4]
	progdata <- progdata[cd4 == "100to199CD4", cat := 5]
	progdata <- progdata[cd4 == "50to99CD4", cat := 6]
	progdata[,risk := -1*log(1 - prob)]
	progdata <- progdata[,.(year,risk,cat)]
	progdata <- progdata[, setattr(as.list(risk), 'names', cat), by = .(year)]
	progdata <- progdata[order(year)]
	progdata <- progdata[,c("1","2","3","4","5","6"), with = F]
	progdata <- as.matrix(progdata)
	for (n in names(dt)) {
        attr(dt[[n]], 'eppfp')$cd4prog <- progdata
    }
    return(dt)
}



#' @title sub.art
#' @description function that substitutes GBD extrapolated ART data past the end of UNAIDS reporting.  This happens in all locations if the art.sub setting
#' is put to TRUE
#'
#' @param epp.input a portion of the EPP data object that will have the art data in it. This is called from prep_epp_data and that function sets this up.
#' @param loc the location of the EPP model, used to find the apropriate art data
#'
#'
#' @return The same EPP data object as before but with GBD extrapolated art data.
#'
#' @export

sub.art <- function(epp.input, loc) {
    if (grepl("KEN", loc) & loc.table[ihme_loc_id == loc, level] == 5) {
        temp.loc <- loc.table[location_id == loc.table[ihme_loc_id == loc, parent_id], ihme_loc_id]
    } else {
        temp.loc <- loc
    }

	for (c.year in c('UNAIDS_2019','UNAIDS_2018','UNAIDS_2017', 'UNAIDS_2016', 'UNAIDS_2015', '140520')) {
	  art.path <- paste0(art.dir, c.year, "/", temp.loc, "_Adult_ART_cov.csv")
	  if (file.exists(art.path)) {
	    art.dt <- fread(art.path)
	    message(paste0("using ART data from ", c.year))
	    break;
	  }
	}

	art.dt <- fread(art.path)
	art.dt[is.na(art.dt)] <- 0
	art.dt[, type := ifelse(ART_cov_pct > 0, "P", "N")]
	years <- epp.input$epp.art$year
	if (max(years) > max(art.dt$year)) {
		max.dt <- copy(art.dt[year == max(year)])
		missing.years <- setdiff(years, art.dt$year)
		add.dt <- rbindlist(lapply(missing.years, function(cyear) {
			copy.dt <- copy(max.dt)
			copy.dt[, year := cyear]
		}))
		art.dt <- rbind(art.dt, add.dt)
	}
	epp.input$epp.art$m.isperc <- art.dt[year %in% years & sex == 1, type]
	epp.input$epp.art$f.isperc <- art.dt[year %in% years & sex == 1, type]
	art.dt[, ART_cov_val := ifelse(ART_cov_pct > 0, ART_cov_pct, ART_cov_num)]
	epp.input$epp.art$m.val <- art.dt[year %in% years & sex == 1, ART_cov_val]
	epp.input$epp.art$f.val <- art.dt[year %in% years & sex == 2, ART_cov_val]
	return(epp.input)
}



#' @title ken_ne_input
#' @description function that manages data processing for the NE KEN 2015 PJNZ.  hopefully we can abandon this as we move to the 2018 or more recent files.
#' I am keeping this here in case we need to repilcate this process in the future.
#'
#' @param pjnz the file path to a pjnz file that will be used as the basis for constructing the EPP data object.  in this case the KEN NE PJNZ file for 2015
#' lacked an XML so it was nearly unusable.  This function cobbeles together data from other sources and constructs a replacement EPP data ojbect.
#'
#' @return an EPP data object that replaces the corrupted PJNZ.
#'
#' @export

ken_ne_input <- function(pjnz){

  ## ep1
  ep1file <- grep(".ep1", unzip(pjnz, list = TRUE)$Name, value = TRUE)
  con <- unz(pjnz, ep1file)
  ep1 <- scan(con, "character", sep = "\n")
  close(con)
  #removing end of line characters
  ep1 <- gsub("^,*|(?<=,),|,*$", "", ep1, perl = T)

  #grep here only needed on windows
  country.idx <- grep("COUNTRY", ep1)
  firstprojyr.idx <-  which(sapply(ep1, substr, 1, 11) == "FIRSTPROJYR")
  lastprojyr.idx <-  which(sapply(ep1, substr, 1, 10) == "LASTPROJYR")
  popstart.idx <- which(ep1 == "POPSTART") + 1
  popend.idx <- which(ep1 == "POPEND") - 1

  country <- as.character(read.csv(text = ep1[country.idx], header = FALSE, as.is = TRUE)[2])
  country.code <- as.integer(read.csv(text = ep1[country.idx], header = FALSE)[3])

  start.year <- as.integer(read.csv(text = ep1[firstprojyr.idx], header = FALSE)[2])
  stop.year <- as.integer(read.csv(text = ep1[lastprojyr.idx], header = FALSE)[2])
  epp.pop <- setNames(read.csv(text = ep1[popstart.idx:popend.idx], header = FALSE, as.is = TRUE),
                      c("year", "pop15to49", "pop15", "pop50", "netmigr"))

  ## ep4
  ep4file <- grep(".ep4", unzip(pjnz, list = TRUE)$Name, value = TRUE)
  con <- unz(pjnz, ep4file)
  ep4 <- scan(con, "character", sep="\n")
  close(con)
  ep4<- gsub("^,*|(?<=,),|,*$", "", ep4, perl=T)

  cd4lim.idx <- which(sapply(ep4, substr, 1, 12) == "CD4LOWLIMITS")
  lambda.idx <- which(sapply(ep4, substr, 1, 6) == "LAMBDA")
  cd4init.idx <- which(sapply(ep4, substr, 1, 13) == "NEWINFECTSCD4")
  mu.idx <- which(sapply(ep4, substr, 1, 3) == "MU_")
  alpha1.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA1")
  alpha2.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA2")
  alpha3.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA3")
  infectreduc.idx <- which(sapply(ep4, substr, 1, 11) == "INFECTREDUC")
  artstart.idx <- which(ep4 == "ARTSTART")+1
  artend.idx <- which(ep4 == "ARTEND")-1

  cd4lim <- as.integer(read.csv(text=ep4[cd4lim.idx], header=FALSE)[-1])
  cd4init <- as.matrix(read.csv(text=ep4[cd4init.idx], header=FALSE, row.names=1))
  lambda <- as.matrix(read.csv(text=ep4[lambda.idx], header=FALSE, row.names=1))
  mu <- as.matrix(read.csv(text=ep4[mu.idx], header=FALSE, row.names=1))
  alpha1 <- as.matrix(read.csv(text=ep4[alpha1.idx], header=FALSE, row.names=1))
  alpha2 <- as.matrix(read.csv(text=ep4[alpha2.idx], header=FALSE, row.names=1))
  alpha3 <- as.matrix(read.csv(text=ep4[alpha3.idx], header=FALSE, row.names=1))
  infectreduc <- as.numeric(read.csv(text=ep4[infectreduc.idx], header=FALSE)[2])

  temp.epp.art <- read.csv(text=ep4[artstart.idx:artend.idx], header=FALSE, as.is=TRUE)
  epp.art.vars <- c("year", "m.isperc", "m.val", "f.isperc", "f.val", "cd4thresh", "m.perc50plus", "f.perc50plus", "perc50plus", "1stto2ndline")

  if (!ncol(temp.epp.art) == length(epp.art.vars)) {
    temp.epp.art$V9 <- (temp.epp.art$V8 + temp.epp.art$V7) / 2
    temp.epp.art$V10 <- 0
  }
  epp.art <- setNames(temp.epp.art, epp.art.vars)

  specpop.idx <- grep("SPECPOP", ep4)
  art.specpop <- read.csv(text=ep4[specpop.idx], header=FALSE,
                          colClasses=c("NULL", "character", "numeric", "integer"),
                          col.names=c(NA, "specpop", "percelig", "yearelig"))
  art.specpop$percelig <- art.specpop$percelig/100


  cd4median.start.idx <- which(ep4 == "CD4MEDIAN_START")+1
  cd4median.end.idx <- which(ep4 == "CD4MEDIAN_END")-1
  if(length(cd4median.start.idx) > 0)
    epp.pop$cd4median <- read.csv(text=ep4[cd4median.start.idx:cd4median.end.idx], header=FALSE, colClasses=c("NULL", "numeric"))[[1]]
  else
    epp.pop$cd4median <- 0

  hivp15yr.start.idx <- which(ep4 == "HIVPOS_15YEAROLDS")+1
  hivp15yr.end.idx <- which(ep4 == "HIVPOS_15YEAROLDS_END")-1
  if(length(hivp15yr.start.idx) > 0)
    epp.pop$hivp15yr <- read.csv(text=ep4[hivp15yr.start.idx:hivp15yr.end.idx], header=FALSE, colClasses=c("NULL", "numeric"))[[1]]
  else
    epp.pop$hivp15yr <- 0

  art15yr.start.idx <- which(ep4 == "HIVPOS_15YEAROLDSART")+1
  art15yr.end.idx <- which(ep4 == "HIVPOS_15YEAROLDSART_END")-1
  if(length(art15yr.start.idx) > 0)
    epp.art$art15yr <- read.csv(text=ep4[art15yr.start.idx:art15yr.end.idx], header=FALSE, colClasses=c("NULL", "numeric"))[[1]]
  else
    epp.art$art15yr <- 0

  artdropout.start.idx <- which(ep4 == "ARTDROPOUTRATE")+1
  artdropout.end.idx <- which(ep4 == "ARTDROPOUTRATE_END")-1
  if(length(artdropout.start.idx) > 0)
    epp.art$artdropout <- read.csv(text=ep4[artdropout.start.idx:artdropout.end.idx], header=FALSE, colClasses=c("NULL", "numeric"))[[1]]
  else
    epp.art$artdropout <- 0

  hivp15yr.cd4dist.idx <- which(ep4 == "HIVPOS15_CD4")+1
  if(length(hivp15yr.cd4dist.idx) > 0)
    hivp15yr.cd4dist <- as.numeric(read.csv(text=ep4[hivp15yr.cd4dist.idx], header=FALSE))
  else
    hivp15yr.cd4dist <- rep(0, length(cd4lim))

  art15yr.cd4dist.idx <- which(ep4 == "HIVPOS15ART_CD4")+1
  if(length(art15yr.cd4dist.idx) > 0)
    art15yr.cd4dist <- as.numeric(read.csv(text=ep4[art15yr.cd4dist.idx], header=FALSE))
  else
    art15yr.cd4dist <- rep(0, length(cd4lim))


  ## Between 2014 and 2015, Spectrum switched from passing CD4 stage duration to passing
  ## CD4 stage progression rate. This change is unmarked in the .ep4 file.

  if(mean(lambda[,1]) > 1)
    lambda <- lambda
  else
    lambda <- 1/lambda

  epidemic.start <- 1975

  eppin <- list(start.year       = start.year,
                stop.year        = stop.year,
                epidemic.start   = epidemic.start,
                epp.pop          = epp.pop,
                cd4lowlim        = cd4lim,
                cd4initperc      = cd4init,
                cd4stage.dur     = lambda,
                cd4mort          = mu,
                artmort.less6mos = alpha1,
                artmort.6to12mos = alpha2,
                artmort.after1yr = alpha3,
                infectreduc      = infectreduc,
                epp.art          = epp.art,
                art.specpop      = art.specpop,
                hivp15yr.cd4dist = hivp15yr.cd4dist,
                art15yr.cd4dist  = art15yr.cd4dist)
  class(eppin) <- "eppin"
  attr(eppin, "country") <- country
  attr(eppin, "country.code") <- country.code

  return(eppin)
}

