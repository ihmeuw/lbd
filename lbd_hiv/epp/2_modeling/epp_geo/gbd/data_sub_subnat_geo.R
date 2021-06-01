################################################################################
## Purpose: These functions substitute data for the admin 2 models into the GBD
## geography level EPP data object that was created by the higher level model.
#
################################################################################



#' @title sub.prev.params.adm2
#' @description function that substitutes the the LBD prevalence estimates for the specific ADM 2 being modeled into the EPP data object
#'
#' @param eppd is the EPP prevalence data used in EPP this function is called in the prep_epp_data_adm2 function which handels queueing up the eppd
#' @param loc the location of the higher level EPP model that this ADM 2 model fits into, used to find the apropriate parameters
#' @param subnat the ADM2_CODE of the admin 2 which we are modeling at this moment
#' @param prev_date the date of the prevalence run which is to be substituted into the EPP object.  This is the data of when the prevalence data was prepared,
#' not the date of the MBG model run being used.  This is becasue there is a possibility that we might make more than one processed prevalence set from a given
#' MBG run.
#'
#'
#' @return An EPP data object but with the higher level prevalence data replaced with the ADM2 prevalence estimates from MBG.
#'
#' @export

sub.prev.params.adm2 <- function(eppd, loc, subnat, prev_date){

  hhs <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"))
  hhs <- hhs[which(hhs$ADM2_CODE == subnat), ]


  eppd[[1]]$anc.used <- FALSE
  eppd[[1]]$anc.prev <- NULL
  eppd[[1]]$anc.n <- NULL
  eppd[[1]]$ADM2_CODE <- subnat
  eppd[[1]]$hhs <- hhs


  return(eppd)
}


#' @title sub.pop.params.adm2
#' @description function that substitutes the the raked World Pop populations for the admin 2 getting modeled in place of the national populations that
#' were natively in the EPP data object
#'
#' @param epp.subp is the EPP pop data used in EPP this function is called in the prep_epp_data_adm2 function which handels queueing up the epp.subp
#' @param loc the location of the higher level EPP model that this ADM 2 model fits into, used to find the apropriate parameters
#' @param subnat the ADM2_CODE of the admin 2 which we are modeling at this moment
#' @param migration_run the name of the admin 2 to admin 2 migration model run that is being used for this EPP run, VESTIGAL SINCE WE FORCE ALL MIGRATION
#' TO BE 0 AND THEN ALLOW THE POPADJUST SETTING TO ACCOUNT FOR IT.  NEED TO EXPLORE REMOVING THIS FROM THE MODEL
#' @param migration_type the tyep of the admin 2 to admin 2 migration model run that is being used for this EPP run, VESTIGAL SINCE WE FORCE ALL MIGRATION
#' TO BE 0 AND THEN ALLOW THE POPADJUST SETTING TO ACCOUNT FOR IT.  NEED TO EXPLORE REMOVING THIS FROM THE MODEL
#' @param shapefile_version The shapefile version we are using for this model run.  this is used to make sure we are getting the correct aggregated pops
#'
#' @return An EPP data object but with the higher level pop data replaced with the ADM2 pop estimates from World Pop that have been raked to GBD totals.
#'
#' @export

sub.pop.params.adm2 <- function(epp.subp, loc, subnat, migration_run, migration_type, shapefile_version) {

  pops <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  pops <- pops[which(pops$ADM2_CODE == subnat & pops$year_id > 1969), ]
  pops <- pops[which(pops$age_group_id != 15), ]
  pops <- as.data.table(pops)
  use.pop <- pops[, .(population = sum(population)), by = .(age_group_id, year_id)]
  sub_n_pop15 <- use.pop[age_group_id == 8, population] / 5
  sub_n_pop50 <- use.pop[age_group_id == 14, population] / 5
  sub_n_pop15to49 <- use.pop[, .(population = sum(population)), by = .(year_id)]$population

  #########################################################################################################

  replace <- as.data.frame(cbind(sort(unique(pops$year_id)), sub_n_pop15to49))
  replace <- cbind(replace, sub_n_pop15)
  replace <- cbind(replace, sub_n_pop50)
  replace$netmigr <- 0 #<- cbind(replace, mig$netmigr) This is set to 0
epp.subp$subpops[[1]] <- replace
names(epp.subp$subpops[[1]]) <- names(epp.subp$total)
names(epp.subp$subpops) <- subnat
return(epp.subp)
}


#' @title sub.art.adm2
#' @description function that adjusts the number of art doses available for the model to distribute to a number that is more reasonable
#'
#' @param epp.input is the EPP parameters data used in EPP which included the number of ART patients. this function is called in the prep_epp_data_adm2 function
#' which handels queueing up the epp.input
#' @param loc the location of the higher level EPP model that this ADM 2 model fits into, used to find the apropriate parameters
#' @param subnat the ADM2_CODE of the admin 2 which we are modeling at this moment
#' @param art_dist_type how we are distributing art, there are several model types to choose from.
#' @param prev_data the prevalence extraction to be used in this EPP run, equal coverage splitting is based on this data.
#'
#' @return An EPP data object but with the higher level ART data scaled down to our estimate for the admin 2
#'
#' @export

sub.art.adm2 <- function(epp.input, loc, subnat, art_dist_type = "Coverage", prev_date) {
  art_data <- epp.input$epp.art

  if (art_dist_type == "Coverage") {
  art.dist <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  }
  if (art_dist_type == "Modeled") {
    warning("you are using modeled admin 2 art this is only possible in MOZ, MWI, CIV, TZA, UGA, RWA, ZAF, and NAM (and the NAM data sucks)")
    art.dist <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  }
  if (art_dist_type == "Hybrid") {
    warning("you are using hybrid admin 2 art this is out of date currecntly and your model is pointless, dispare and go fix it")
    art.dist <- fread("<<<< FILEPATH REDACTED >>>>")
  }

  art.dist <- art.dist[which(art.dist$ADM2_CODE == subnat), ]

  art_data <- as.data.frame(merge(art_data, art.dist[ , c("year", "prev_prop")], by = c("year")))
  art_data$prev_prop[which(is.na(art_data$prev_prop))] <- 0 #deals with a divide by 0 problem coming out of the modeled art proportions.  (if there is 0 art in a coutry the proportion is NA)
  art_data$m.val <- art_data$m.val * art_data$prev_prop
  art_data$f.val <- art_data$f.val * art_data$prev_prop

  epp.input$epp.art <- art_data
  return(epp.input)

}

#' @title fnCreateEPPSubpops_adm2
#' @description function that organizes the EPP data ojbects as needed to run EPP
#'
#' @param epp.input is the EPP parameters
#' @param epp.subpops is the population info
#' @param epp.data is the prevalence info that the model is fit to
#' @param no.anc indicates if anc data should be kept or stripped from the EPP object
#' @param stop.year is the last year the model will be run
#'
#' @return An EPP data object for the admin 2 of interest that is structured as needed.
#'
#' @export

fnCreateEPPSubpops_adm2 <- function(epp.input, epp.subpops, epp.data, no.anc = no.anc, stop.year){

  epp.subpop.input <- list()

  for (subpop in names(epp.subpops$subpops)) {
    epp.subpop.input[[subpop]] <- epp.input
    epp.subpop.input[[subpop]]$epp.pop <- epp.subpops$subpops[[subpop]]
    epp.subpop.input[[subpop]]$epp.pop$cd4median <- epp.input$epp.pop$cd4median[which((epp.input$epp.pop$year) <= max(epp.subpop.input[[subpop]]$epp.pop$year))]
    epp.subpop.input[[subpop]]$epp.pop$hivp15yr <- epp.input$epp.pop$hivp15yr[which((epp.input$epp.pop$year) <= max(epp.subpop.input[[subpop]]$epp.pop$year))]

    epp.art <- epp.input$epp.art

    epp.subpop.input[[subpop]]$epp.art <- epp.art

    if (!is.null(attr(epp.subpops$subpops[[subpop]], "epidemic.start")))
      epp.subpop.input[[subpop]]$epidemic.start <- attr(epp.subpops$subpops[[subpop]], "epidemic.start")
  }

  return(epp.subpop.input)
}



#' @title fnCreateEPPSubpops_adm2
#' @description function that extends the EPP data object tohave the decimal year time steps and interpotates the appropriate inputs.
#'
#' @param epp.input The EPP input data that will end up betting spread over the decimal year time steps
#' @param dt numeric, the fraction of a year in each time step
#' @param proj.start the first time step, usually half a year after the start year because EPP reports out mid year estimates
#' @param proj.end the last time step of the model, usually the stop year plus 0.5 NEED TO FIX THIS FOR MORT TRACKIGNG HERE AND IN THE NATIONAL MODEL
#' @param tsEpidemicStart the time step when the iota prior initiates the epidemic
#' @param cd4stage.weights weights given to various CD4 statges
#' @param art1yr.weight
#' @param num.knots number of knots in the spline parameterization of the force of infection
#' @param no.anc to remove ANC data or not to remove ANC data
#' @param gbd.mort to use GBD mortality parameters or not to use GBD mortality parameters.
#'
#' @return An EPP data object for the admin 2 of interest that is structured as needed with the apropriate interpolation fo the inputs to the decimal year time steps.
#' NOTE: THERE IS A SCALAR OF 1 APPLIED TO THE HIV FREE MORT IN THIS FUCNTION, THIS CAN BE CHANGED FOR SENSITIVITY ANALYSIS
#'
#' @export

fnCreateEPPFixPar_adm2 <- function(epp.input,
                                   dt = 0.1,
                                   proj.start = start.year + 0.5,
                                   proj.end = stop.year + 0.5,
                                   tsEpidemicStart = epp.input$epidemic.start + dt*ceiling(1/(2*dt)),
                                   cd4stage.weights = c(1.3, 0.6, 0.1, 0.1, 0.0, 0.0, 0.0),
                                   art1yr.weight = 0.1,
                                   num.knots = 7,
                                   no.anc,
                                   gbd.mort = FALSE, gbd_mort_run = gbd_mort_run){

  shp_v <- epp.input$shp
  locat <- epp.input$loc
  subnat <- epp.input$subnat

  #########################
  ##  Mortality  inputs  ##
  #########################
 if (gbd.mort == TRUE) {
  path <- paste0("<<<< FILEPATH REDACTED >>>>")
  pops <- fread(path)
  age_groups <- unique(pops$age_group_id)

  gbd_q <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
  gbd_q <- gbd_q[which(gbd_q$ihme_loc_id == loc & gbd_q$estimate_stage_id == 7), ]
  gbd_q <- gbd_q[ ,c("age_group_id", "sex_id", "year_id", "mean")]


  p <- pops[which(pops$ADM2_CODE == subnat), ]

  p <- merge(p, gbd_q, by = c("age_group_id", "sex_id", "year_id"))

  p$weighted <- p$population * p$mean

  collapsed <- p[,lapply(c("weighted", "population"), function(x) sum(get(x), na.rm = T)), by = c("year_id", "ADM2_CODE")]
  collapsed$mx <- (collapsed$V1/collapsed$V2)*1 ################################################ This is a mutiplier that can be used to test the model's sensitivity to HIV free mortality levels

 }



  #########################
  ##  Population inputs  ##
  #########################
  if (gbd.mort == TRUE) {
  epp.pop <- merge(epp.input$epp.pop, epp.input$epp.art[c("year", "art15yr")], all.x = TRUE)
  epp.pop <- merge(epp.pop, collapsed[ , c("year_id", "mx")], by.x = c("year"), by.y = c("year_id"))
  epp.pop$art15yr[is.na(epp.pop$art15yr)] <- 0
  proj.steps <- seq(proj.start, proj.end, dt)
  epp.pop.ts <- data.frame(pop15to49  = approx(epp.pop$year + 0.4, epp.pop$pop15to49,     proj.steps, rule = 2)$y,
                           age15enter = approx(epp.pop$year + 0.4, dt * epp.pop$pop15,    proj.steps, rule = 2)$y,
                           age50exit  = approx(epp.pop$year + 0.4, dt * epp.pop$pop50,    proj.steps, rule = 2)$y,
                           netmigr    = approx(epp.pop$year + 0.4, dt * epp.pop$netmigr,  proj.steps, rule = 2)$y,
                           hivp15yr   = approx(epp.pop$year + 0.4, dt * epp.pop$hivp15yr, proj.steps, rule = 2)$y,
                           art15yr    = approx(epp.pop$year + 0.4, dt * epp.pop$art15yr,  proj.steps, rule = 2)$y,
                           mx         = approx(epp.pop$year + 0.4, dt * epp.pop$mx,       proj.steps, rule = 2)$y)

  proj.years <- floor(proj.steps)

  epp.pop.ts$netmigr <- epp.pop.ts$netmigr * with(subset(epp.pop, epp.pop$year %in% proj.years), rep(ifelse(netmigr != 0, netmigr * table(proj.years) * dt / tapply(epp.pop.ts$netmigr, floor(proj.steps), sum), 0), times = table(proj.years)))
  epp.pop.ts$hivp15yr <- epp.pop.ts$hivp15yr * with(subset(epp.pop, epp.pop$year %in% proj.years), rep(ifelse(hivp15yr != 0, hivp15yr * table(proj.years) * dt / tapply(epp.pop.ts$hivp15yr, floor(proj.steps), sum), 0), times = table(proj.years)))
  epp.pop.ts$art15yr <- epp.pop.ts$art15yr * with(subset(epp.pop, epp.pop$year %in% proj.years), rep(ifelse(art15yr != 0, art15yr * table(proj.years) * dt / tapply(epp.pop.ts$art15yr, floor(proj.steps), sum), 0), times = table(proj.years)))
  epp.pop.ts$age50rate <- (epp.pop.ts$age50exit/epp.pop.ts$pop15to49)/dt

  epp.pop.ts$mx <- epp.pop.ts$mx*dt
  } else {
      epp.pop <- merge(epp.input$epp.pop, epp.input$epp.art[c("year", "art15yr")], all.x = TRUE)
      epp.pop$art15yr[is.na(epp.pop$art15yr)] <- 0
      proj.steps <- seq(proj.start, proj.end, dt)
      epp.pop.ts <- data.frame(pop15to49  = approx(epp.pop$year + 0.4, epp.pop$pop15to49,     proj.steps, rule = 2)$y,
                               age15enter = approx(epp.pop$year + 0.4, dt * epp.pop$pop15,    proj.steps, rule = 2)$y,
                               age50exit  = approx(epp.pop$year + 0.4, dt * epp.pop$pop50,    proj.steps, rule = 2)$y,
                               netmigr    = approx(epp.pop$year + 0.4, dt * epp.pop$netmigr,  proj.steps, rule = 2)$y,
                               hivp15yr   = approx(epp.pop$year + 0.4, dt * epp.pop$hivp15yr, proj.steps, rule = 2)$y,
                               art15yr    = approx(epp.pop$year + 0.4, dt * epp.pop$art15yr,  proj.steps, rule = 2)$y)

      proj.years <- floor(proj.steps)

      epp.pop.ts$netmigr <- epp.pop.ts$netmigr * with(subset(epp.pop, epp.pop$year %in% proj.years), rep(ifelse(netmigr != 0, netmigr * table(proj.years) * dt / tapply(epp.pop.ts$netmigr, floor(proj.steps), sum), 0), times = table(proj.years)))
      epp.pop.ts$hivp15yr <- epp.pop.ts$hivp15yr * with(subset(epp.pop, epp.pop$year %in% proj.years), rep(ifelse(hivp15yr != 0, hivp15yr * table(proj.years) * dt / tapply(epp.pop.ts$hivp15yr, floor(proj.steps), sum), 0), times = table(proj.years)))
      epp.pop.ts$art15yr <- epp.pop.ts$art15yr * with(subset(epp.pop, epp.pop$year %in% proj.years), rep(ifelse(art15yr != 0, art15yr * table(proj.years) * dt / tapply(epp.pop.ts$art15yr, floor(proj.steps), sum), 0), times = table(proj.years)))
      epp.pop.ts$age50rate <- (epp.pop.ts$age50exit/epp.pop.ts$pop15to49)/dt
      epp.pop.ts$mx <- c(1.0 - (epp.pop.ts$pop15to49[-1] - (epp.pop.ts$age15enter - epp.pop.ts$age50exit + epp.pop.ts$netmigr)[-length(proj.steps)]) / epp.pop.ts$pop15to49[-length(proj.steps)], NA) / dt
      epp.pop.ts[length(proj.steps), "mx"] <- epp.pop.ts[length(proj.steps) - 1, "mx"]
  }

  ##################
  ##  ART inputs  ##
  ##################

  epp.art <- epp.input$epp.art

  ## number of persons who should be on ART at end of timestep
  artnum.ts <- with(subset(epp.art, m.isperc == "N"), approx(year + 1 - dt, (m.val + f.val)*(1.0 - perc50plus/100), proj.steps, rule = 2))$y  # offset by 1 year because number on ART are Dec 31
  ## !! Currently only dealing with numbers on ART, assuming no increase in coverage after last number

  epp.art$artelig.idx <- match(epp.art$cd4thresh, c(1, 2, 500, 350, 250, 200, 100, 50))
  artelig.idx.ts <- approx(epp.art$year, epp.art$artelig.idx, proj.steps, "constant", rule = 2)$y

  epp.art$specpop.percelig <- rowSums(with(epp.input$art.specpop, mapply(function(percent, year) rep(c(0, percent), c(year - min(epp.art$year), max(epp.art$year) - year + 1)), percelig, yearelig)))
  specpop.percelig.ts <- approx(epp.art$year + 0.4, epp.art$specpop.percelig, proj.steps, "constant", rule = 2)$y

  cd4prog <- 1/colMeans(epp.input$cd4stage.dur[c(2:3,6:7),])
  cd4init <- 0.01*colMeans(epp.input$cd4initperc[c(2:3,6:7),])
  cd4artmort <- cbind(colMeans(epp.input$cd4mort[c(2:3,6:7),]),
                      colMeans(epp.input$artmort.less6mos[c(2:3,6:7),]),
                      colMeans(epp.input$artmort.6to12mos[c(2:3,6:7),]),
                      colMeans(epp.input$artmort.after1yr[c(2:3,6:7),]))

  relinfectART <- 1.0 - epp.input$infectreduc


  ###########################
  ##  r-spline parameters  ##
  ###########################
  inner.knots <- num.knots
  if (inner.knots == 1) {
    ord <- 1
    rvec.spldes <- matrix(1, nrow = length(proj.steps))
  } else {
    if (inner.knots == 2) {
      ord <- 2
    } else {
      ord <- 4
    }

    numKnots <- ord + inner.knots
    proj.dur <- diff(range(proj.steps))
    rvec.knots <- seq(min(proj.steps) - (numKnots - inner.knots) / 2 * (proj.dur / (inner.knots - 1)),
                      max(proj.steps) + (numKnots - inner.knots) / 2 * (proj.dur / (inner.knots - 1)),
                      proj.dur / (inner.knots - 1))
    rvec.spldes <- splines::splineDesign(rvec.knots, proj.steps, ord = ord, outer.ok = T)
  }

    val <- list(proj.steps      = proj.steps,
                tsEpidemicStart = tsEpidemicStart,
                dt              = dt,
                epp.pop.ts      = epp.pop.ts,
                artnum.ts       = artnum.ts,
                artelig.idx.ts  = artelig.idx.ts,
                specpop.percelig.ts = specpop.percelig.ts,
                cd4prog         = cd4prog,
                cd4init         = cd4init,
                cd4artmort      = cd4artmort,
                relinfectART    = relinfectART,
                hivp15yr.cd4dist = epp.input$hivp15yr.cd4dist,
                art15yr.cd4dist = epp.input$art15yr.cd4dist,
                numKnots        = inner.knots,
                rvec.spldes     = rvec.spldes,
                iota            = 0.0025)

  class(val) <- "eppfp"
  return(val)
}

#' @title calc.expand.pop.adm2
#' @description performs the same function as calc.expand.pop but for admin 2 units
#'
#' @param loc the high level location that the admin 2 model fits into.
#' @param sex.agg do we aggregate sexes or not
#' @param subnat the ADM2_CODE for the admin 2 unit that we are modeling.
#' @param shapefile_version the shapefile version we are modeling on, this determines the relationships between different data inputs since the ADM2_CODEs can change
#'
#' @return expands the population structure using IRR to best estimate the HIV+ motality and CD4 progression parameters
#'
#' @export


calc.expand.pop.adm2 <- function(loc, sex.agg = T, subnat, shapefile_version, test_IRR = F) {

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
  in.pop <- read.csv(path)
  in.pop <- in.pop[which(in.pop$ADM2_CODE == subnat), ]
  in.pop <- in.pop[which(in.pop$year_id %in% start.year:stop.year), ]
  in.pop <- as.data.table(in.pop)
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
    test_IRRs$comb_IRR[which(test_IRRs$mean == 0)] <- 1 #fix the 0/0 issue
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

#' @title sub.on.art.adm2
#' @description performs the same function as sub.on.art but for admin 2 units
#'
#' @param dt the EPP data object that will have its on art HIV+ mortality parameters replaced
#' @param loc the high level location that the admin 2 model fits into.
#' @param k the draw number to use for selecting the on art HIV+ mortality
#' @param subnat the ADM2_CODE for the admin 2 unit that we are modeling.
#' @param shapefile_version the shapefile version we are modeling on, this determines the relationships between different data inputs since the ADM2_CODEs can change
#'
#' @return EPP object with new on art HIV+ mortality parameters
#'
#' @export

sub.on.art.adm2 <- function(dt, loc, k, subnat, shapefile_version, mean_test = FALSE, tirr = FALSE) {
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
  expanded_pop <- calc.expand.pop.adm2(loc, sex.agg = F, subnat, shapefile_version, test_IRR = tirr)
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

#' @title sub.off.art.adm2
#' @description performs the same function as sub.off.art but for admin 2 units
#'
#' @param dt the EPP data object that will have its off art HIV+ mortality parameters replaced
#' @param loc the high level location that the admin 2 model fits into.
#' @param k the draw number to use for selecting the off art HIV+ mortality
#' @param subnat the ADM2_CODE for the admin 2 unit that we are modeling.
#' @param shapefile_version the shapefile version we are modeling on, this determines the relationships between different data inputs since the ADM2_CODEs can change
#'
#' @return EPP object with new off art HIV+ mortality parameters
#'
#' @export

sub.off.art.adm2 <- function(dt, loc, k, subnat, shapefile_version, mean_test = FALSE, tirr = FALSE) {

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
  sex_agg <- calc.expand.pop.adm2(loc,sex.agg = T, subnat, shapefile_version, test_IRR = tirr)
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


#' @title sub.cd4.prog.adm2
#' @description performs the same function as sub.cd4.prog but for admin 2 units
#'
#' @param dt the EPP data object that will have its on CD$ progression parameters replaced
#' @param loc the high level location that the admin 2 model fits into.
#' @param k the draw number to use for selecting the CD4 progression
#' @param subnat the ADM2_CODE for the admin 2 unit that we are modeling.
#' @param shapefile_version the shapefile version we are modeling on, this determines the relationships between different data inputs since the ADM2_CODEs can change
#'
#' @return EPP object with new CD4 progression parameters
#'
#' @export

sub.cd4.prog.adm2 <- function(dt, loc, k, subnat, shapefile_version, mean_test = FALSE, tirr = FALSE) {
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

  sex_agg <- calc.expand.pop.adm2(loc, sex.agg = T, subnat, shapefile_version, test_IRR = tirr)
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
