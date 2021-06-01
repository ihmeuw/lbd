
###########################################################################
####  Function to read inputs from Spectrum to EPP (.ep1, .ep3, .ep4)  ####
###########################################################################
ind.prepare.epp.fit <- function(filepath, proj.end=2016.5, sub.anc = T, sub.art.cov=T){

  ##File paths
  sub.art.cov.path <- paste0("<<<< FILEPATH REDACTED >>>>")
  sub.anc.path <- paste0("<<<<FILEPATH REDACTED>>>>>")
  sub.anc.ibbs.path <- paste0("<<<<FILEPATH REDACTED>>>>>")

  loc_id <- loc.table[ihme_loc_id == loc, location_id]

  # old loc id map
  input.loc.map <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))[, .(ihme_loc_id, gbd15_spectrum_loc)]
  input.loc.map[ ihme_loc_id=="IND_4871", gbd15_spectrum_loc:= "IND_TL"]
  iso3_init <- input.loc.map[ihme_loc_id == loc, gbd15_spectrum_loc]

  # Load get_population
  if(!"get_population" %in% ls()) {
    source(paste0("<<<< FILEPATH REDACTED >>>>/get_population.R"))
  }

  ep1 <- scan(paste(filepath, ".ep1", sep=""), "character", sep="\n")
  ep1 <<- ep1[3:length(ep1)]
  ep4 <- scan(paste(filepath, ".ep3", sep=""), "character", sep="\n")
  ep4 <<- ep4[3:length(ep4)]

  ## epp
  eppd <- ind.read.epp.data(paste(filepath, ".xml", sep=""))
  if(sub.anc) {
    ## Prevalence
    # Prep ANC prevalence data for insertion
    gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
    sub.anc.dt <- fread(sub.anc.path, header = T)[gbd15_spectrum_loc == iso3_init]
    gen.pop <- names(eppd)[names(eppd) %in% gen.pop.dict]
    year.range <- colnames(eppd[[gen.pop]]$anc.prev)
    # Add years where needed
    for(add.year in setdiff(year.range, names(sub.anc.dt))) {
      sub.anc.dt[, (add.year) := NA]
    }
    # Convert to matrix
    sub.anc.matrix <- as.matrix(sub.anc.dt[, year.range, with=F]) / 100  # Convert from per 100 to per 1
    rownames(sub.anc.matrix) <- sub.anc.dt[["Site Name"]]
    # Sub in prevalence
    eppd[[gen.pop]]$anc.prev <- sub.anc.matrix

    ## Sample size
    # Prep sample size for insertion
    sub.anc.matrix[!is.na(sub.anc.matrix)] <- 400
    eppd[[gen.pop]]$anc.n <- sub.anc.matrix

    ## Used True/False
    eppd[[gen.pop]]$anc.used <- rep(TRUE, nrow(sub.anc.matrix))

    ##### Substitue the subpop anc with ibbs prev rate and sample size
    sub.anc.ibbs.dt <- fread(sub.anc.ibbs.path, header = T)[ihme_loc_id == loc_id]
    for (sub.pop in unique(sub.anc.ibbs.dt$subpop)){
    subpop.anc.ibbs.dt <- sub.anc.ibbs.dt[subpop==sub.pop]
      for (year.id in unique(subpop.anc.ibbs.dt$year)){
      year.range <- colnames(eppd[[sub.pop]]$anc.prev)
      ## add prev rate
      add.matrix <- t(rep(NA, length(year.range)))
      colnames(add.matrix) <- year.range
      rownames(add.matrix) <- "IBBS"
      add.matrix[, as.character(year.id)] <- subpop.anc.ibbs.dt[subpop==sub.pop & year==year.id]$prev/100
      eppd[[sub.pop]]$anc.prev <- rbind(eppd[[sub.pop]]$anc.prev, add.matrix)
      ## add sample size n
      n.add.matrix <- t(rep(NA, length(year.range)))
      colnames(n.add.matrix) <- year.range
      rownames(n.add.matrix) <- "IBBS"
      n.add.matrix[, as.character(year.id)] <- subpop.anc.ibbs.dt[subpop==sub.pop & year==year.id]$n
      eppd[[sub.pop]]$anc.n <- rbind(eppd[[sub.pop]]$anc.n, n.add.matrix)
      ## add used index
      eppd[[sub.pop]]$anc.used <- rep(TRUE, nrow(eppd[[sub.pop]]$anc.n))
      }
    }
  }

  epp.subp <- ind.read.epp.subpops(paste(filepath, ".xml", sep=""))
  ### add in GBD adult pop to scale the subpops
  # if(iso3_init %in% c("IND_AP", "IND_TL")) {

    loc_id <- input.loc.map[gbd15_spectrum_loc==iso3_init, ihme_loc_id]

    location_id <- loc.table[ihme_loc_id==loc_id, location_id]
    gbd16.pop <- get_population(year_id = -1, location_id = location_id, sex_id = 3, age_group_id = 24, gbd_round_id = 5)
    epp.adult.pop <- data.table(year_id=epp.subp$total$year, epp_pop=epp.subp$total$pop15to49)
    merged.pop <- merge(epp.adult.pop, gbd16.pop[,.(year_id, population)], by="year_id", all.x=T)
    merged.pop[,scale:=population/epp_pop]
    scale <- merged.pop[,scale]
    scale <- c(scale, rep(scale[length(scale)], (nrow(epp.subp$total)- length(scale))))

    for (pop in names(epp.subp)) {
    temp.dt <- epp.subp[[pop]]
          if(pop=="subpops") {
            for (subpop in names(temp.dt)){
               temp.dt1 <- temp.dt[[subpop]]
             for (col_id in c("pop15to49", "pop15", "pop50", "netmigr")){
              epp.subp[[pop]][[subpop]][,col_id] <- temp.dt1[,col_id]*scale
             }
            }
          } else {
              for (col_id in c("pop15to49", "pop15", "pop50", "netmigr")){
                epp.subp[[pop]][,col_id] <- temp.dt[,col_id]*scale
              }
            }
    }
 # }
  ###########################
  epp.input <- ind.read.epp.input(filepath)

  if(sub.art.cov) {
    art.cov.dt <- fread(sub.art.cov.path)
    year.range <- unique(epp.input[["epp.art"]][["year"]])
    subset.art.cov <- art.cov.dt[year %in% year.range]
    for(sex in c("m", "f")) {
      sex_id <- ifelse(sex == "m",  1, 2)
      epp.input[["epp.art"]][[paste0(sex, ".val")]] <- round(subset.art.cov[sex == sex_id, ART_cov_num])
    }
  }

    ind_totals <- list(eppd.ind = eppd, epp.subp.ind = epp.subp, epp.input.ind = epp.input)
    return(ind_totals)
}


ind.read.epp.input <- function(ep.path){

  ## ep1
  # ep1 <- scan(paste(ep.path, ".ep1", sep=""), "character", sep="\n")
  # ep1 <- ep1[3:length(ep1)]

  firstprojyr.idx <-  which(sapply(ep1, substr, 1, 11) == "FIRSTPROJYR")
  lastprojyr.idx <-  which(sapply(ep1, substr, 1, 10) == "LASTPROJYR")
  popstart.idx <- which(ep1 == "POPSTART")+1
  popend.idx <- which(ep1 == "POPEND")-1

  start.year <- as.integer(read.csv(text=ep1[firstprojyr.idx], header=FALSE)[2])
  stop.year <- as.integer(read.csv(text=ep1[lastprojyr.idx], header=FALSE)[2])

  # if (iso3 %in% c('AGO', 'SSD'))
  #   start.year <- 1985

  epp.pop <- setNames(read.csv(text=ep1[popstart.idx:popend.idx], header=FALSE, as.is=TRUE),
                      c("year", "pop15to49", "pop15", "pop50", "netmigr"))

  ## ep4
  # ep4 <- scan(paste(ep.path, ".ep4", sep=""), "character", sep="\n")
  # ep4 <- ep4[3:length(ep4)]

  cd4lim.idx <- which(sapply(ep4, substr, 1, 12) == "CD4LOWLIMITS")
  lambda.idx <- which(sapply(ep4, substr, 1, 6) == "LAMBDA")
  cd4init.idx <- which(sapply(ep4, substr, 1, 13) == "NEWINFECTSCD4")
  mu.idx <- which(sapply(ep4, substr, 1, 2) == "MU")
  alpha1.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA1")
  alpha2.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA2")
  alpha3.idx <- which(sapply(ep4, substr, 1, 6) == "ALPHA3")
  infectreduc.idx <- which(sapply(ep4, substr, 1, 11) == "INFECTREDUC")
  artstart.idx <- which(ep4 == "ARTSTART")+1
  artend.idx <- which(ep4 == "ARTEND")-1

  cd4lim <- as.integer(read.csv(text=ep4[cd4lim.idx], header=FALSE)[-1])
  cd4init <- as.matrix(read.csv(text=ep4[cd4init.idx], header=FALSE, row.names=1))
  lambda <- as.matrix(read.csv(text=ep4[lambda.idx], header=FALSE, row.names=1))
  # row.names(progdata) <- row.names(lambda)

  cd4init <- t(as.matrix(c(100,0)))
  cd4init[,1] <- 100
  cd4init[,2] <- 0

  #lambda[,1:6] <- 0.1*progdata[,1:6]
  #lambda <- 0.1 * progdata[,1:6]
  # lambda <- progdata

  #lambda <- lambda*ratio

  mu <- as.matrix(read.csv(text=ep4[mu.idx], header=FALSE, row.names=1))
  alpha1 <- as.matrix(read.csv(text=ep4[alpha1.idx], header=FALSE, row.names=1))
  alpha2 <- as.matrix(read.csv(text=ep4[alpha2.idx], header=FALSE, row.names=1))
  alpha3 <- as.matrix(read.csv(text=ep4[alpha3.idx], header=FALSE, row.names=1))

  # row.names(mugbd) <- row.names(mu)
  # row.names(alpha1gbd) <- row.names(alpha1)
  # row.names(alpha2gbd) <- row.names(alpha2)
  # row.names(alpha3gbd) <- row.names(alpha3)

  #mu <- 10*mugbd
  # mu <- mugbd
  #mu <- mu*ratio
  # alpha1 <- alpha1gbd
  # alpha2 <- alpha2gbd
  # alpha3 <- alpha3gbd

  #alpha1 <- alpha1*ratio
  #alpha2 <- alpha2*ratio
  #alpha3 <- alpha3*ratio

  infectreduc <- as.numeric(read.csv(text=ep4[infectreduc.idx], header=FALSE)[2])

  epp.art <- read.csv(text=ep4[artstart.idx:artend.idx], header=FALSE, as.is=TRUE)
  if (length(names(epp.art)) == 6) {
    names(epp.art) <- c("year", "m.isperc", "m.val", "f.isperc", "f.val", "cd4thresh")
    # epp.art <- setNames(read.csv(text=ep4[artstart.idx:artend.idx], header=FALSE, as.is=TRUE),
    #                     c("year", "m.isperc", "m.val", "f.isperc", "f.val", "cd4thresh"))
    epp.art$m.perc50plus <- 0
    epp.art$f.perc50plus <- 0
    epp.art$perc50plus <- 0

  } else {
    names(epp.art) <- c("year", "m.isperc", "m.val", "f.isperc", "f.val", "cd4thresh", "m.perc50plus", "f.perc50plus", "perc50plus", "")
  }
  first.art.i <- min(which(epp.art$m.val > 0))
  epp.art$m.isperc[1:(first.art.i-1)] <- epp.art$m.isperc[first.art.i]


  eppin <- list(start.year       = start.year,
                stop.year        = stop.year,
                epp.pop          = epp.pop,
                cd4lowlim        = cd4lim,
                cd4initperc      = cd4init,
                cd4stage.dur     = lambda,
                cd4mort          = mu,
                artmort.less6mos = alpha1,
                artmort.6to12mos = alpha2,
                artmort.after1yr = alpha3,
                infectreduc      = infectreduc,
                epp.art          = epp.art)
  class(eppin) <- "eppin"

  return(eppin)
}


############################################################################
####  Function to read prevalence data used in EPP fitting (from .xml)  ####
############################################################################

library(XML)

ind.read.epp.data <- function(epp.xml){
  obj <- xmlTreeParse(epp.xml)

  r <- xmlRoot(obj)[[1]]
  eppSetChildren.idx <- which(xmlSApply(r, xmlAttrs) == "eppSetChildren")
  country <- xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetCountry")]][[1]])

  epp.data <- list() # declare list to store output
  attr(epp.data, "country") <- country
  used.indices <- NULL
  for(eppSet.idx in 1:xmlSize(r[[eppSetChildren.idx]])){

    tmp.eppSet <- r[[eppSetChildren.idx]][[eppSet.idx]][[1]]
    n.iter <- 1
    if (xmlSize(xmlSApply(tmp.eppSet, xmlAttrs)) > 0 & xmlSize(which(xmlSApply(tmp.eppSet, xmlAttrs) == "siteNames")) == 0) {
      newSetChildren.idx <- which(xmlSApply(tmp.eppSet, xmlAttrs) == "eppSetChildren")
      n.iter <- xmlSize(tmp.eppSet[[newSetChildren.idx]])
    }
    for (eppSet.idx.2 in 1:n.iter) {
      eppSet <- tmp.eppSet
      if (xmlSize(xmlSApply(tmp.eppSet, xmlAttrs)) > 0 & xmlSize(which(xmlSApply(tmp.eppSet, xmlAttrs) == "siteNames")) == 0)
        eppSet <- tmp.eppSet[[newSetChildren.idx]][[eppSet.idx.2]][[1]]

      if (xmlSize(eppSet) == 0) {

        tmp_len <- 0
        tmp_index <- 0
        obj_size <- 0
        while(obj_size == 0) {
          tmp_index <- tmp_index + 1
          if (tmp_index != eppSet.idx & !(tmp_index %in% used.indices)) {
            print(paste('tmp_index',tmp_index))
            eppSet <- r[[eppSetChildren.idx]][[tmp_index]][[1]]
            tmp_len <- length(which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo"))
            if (tmp_len > 0) {
              tmp_eppSet <- eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo")]]
              set_location <- ifelse(xmlSize(tmp_eppSet) == 1, 1, which(xmlSApply(tmp_eppSet, xmlAttrs) == "epp2011.core.sets.ProjectionSet"))
              tmp_eppSet <- tmp_eppSet[[set_location]]
              obj_size <- xmlSize(tmp_eppSet)
            }
          }
        }
        used.indices <- c(used.indices, tmp_index)
        eppSet <- eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo")]]
        set_location <- ifelse(xmlSize(eppSet) == 1, 1, which(xmlSApply(eppSet, xmlAttrs) == "epp2011.core.sets.ProjectionSet"))
        eppSet <- eppSet[[set_location]]
      }
      eppName <- xmlToList(eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "name")]][["string"]])
      print(eppName)
      ##  ANC data  ##

      siteNames.idx <- which(xmlSApply(eppSet, xmlAttrs) == "siteNames")
      siteSelected.idx <- which(xmlSApply(eppSet, xmlAttrs) == "siteSelected")
      survData.idx <- which(xmlSApply(eppSet, xmlAttrs) == "survData")
      survSampleSizes.idx <- which(xmlSApply(eppSet, xmlAttrs) == "survSampleSizes")

      siteNames <- xmlSApply(eppSet[[siteNames.idx]][[1]], xmlSApply, xmlToList, FALSE)
      siteIdx <- as.numeric(xmlSApply(eppSet[[siteNames.idx]][[1]], xmlAttrs)) ## 0 based
      anc.used <- as.logical(xmlSApply(eppSet[[siteSelected.idx]][[1]], xmlSApply, xmlToList, FALSE))

      nsites <- length(siteNames)
      nANCyears <- max(as.integer(xmlSApply(eppSet[[survData.idx]][["array"]][[1]][[1]], xmlAttrs))) + 1

      ## ANC prevalence
      anc.prev <- matrix(NA, nsites, nANCyears)
      rownames(anc.prev) <- siteNames
      colnames(anc.prev) <- 1985+0:(nANCyears-1)
      for(clinic.idx in 1:nsites){
        clinic <- eppSet[[survData.idx]][["array"]][[clinic.idx]][[1]]
        prev <- as.numeric(xmlSApply(clinic, xmlSApply, xmlToList, FALSE))
        idx <- as.integer(xmlSApply(clinic, xmlAttrs)) + 1
        anc.prev[clinic.idx, idx] <- prev
      }
      anc.prev[is.na(anc.prev)] <- 0.0 ## NOTE: appears that if value is 0.0, the array index is omitted from XML file, might apply elsewhere.
      anc.prev[anc.prev == -1] <- NA
      anc.prev <- anc.prev/100

      ## ANC sample sizes
      anc.n <- matrix(NA, nsites, nANCyears)
      rownames(anc.n) <- siteNames
      colnames(anc.n) <- 1985+0:(nANCyears-1)
      for(clinic.idx in 1:nsites){
        clinic <- eppSet[[survSampleSizes.idx]][["array"]][[clinic.idx]][[1]]
        n <- as.numeric(xmlSApply(clinic, xmlSApply, xmlToList, FALSE))
        idx <- as.integer(xmlSApply(clinic, xmlAttrs)) + 1
        anc.n[clinic.idx, idx] <- n
      }
      anc.n[anc.n == -1] <- NA

      # anc.prev[,colnames(anc.prev) <= 2000] <- NA
      # anc.n[,colnames(anc.n) <= 2000] <- NA
      ##  HH surveys  ##

      hhsUsed.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyIsUsed")
      hhsHIV.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyHIV")
      hhsSampleSize.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveySampleSize")
      hhsSE.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyStandardError")
      hhsYear.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyYears")

      nhhs <- max(xmlSize(eppSet[[hhsYear.idx]][[1]]),
                  xmlSize(eppSet[[hhsHIV.idx]][[1]]),
                  xmlSize(eppSet[[hhsSE.idx]][[1]]),
                  xmlSize(eppSet[[hhsSampleSize.idx]][[1]]),
                  xmlSize(eppSet[[hhsUsed.idx]][[1]]))

      hhs <- data.frame(year = rep(NA, nhhs), prev = rep(NA, nhhs), se = rep(NA, nhhs), n = rep(NA, nhhs), used = rep(NA, nhhs))

      hhs$year[as.integer(xmlSApply(eppSet[[hhsYear.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[hhsYear.idx]][[1]], xmlSApply, xmlToList, FALSE))
      hhs$prev[as.integer(xmlSApply(eppSet[[hhsHIV.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[hhsHIV.idx]][[1]], xmlSApply, xmlToList, FALSE))/100
      hhs$se[as.integer(xmlSApply(eppSet[[hhsSE.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[hhsSE.idx]][[1]], xmlSApply, xmlToList, FALSE))/100
      hhs$n[as.integer(xmlSApply(eppSet[[hhsSampleSize.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[hhsSampleSize.idx]][[1]], xmlSApply, xmlToList, FALSE))
      hhs$used[as.integer(xmlSApply(eppSet[[hhsUsed.idx]][[1]], xmlAttrs))+1] <- as.logical(xmlSApply(eppSet[[hhsUsed.idx]][[1]], xmlSApply, xmlToList, FALSE))

      hhs <- subset(hhs, !is.na(prev))

      epp.data[[eppName]] <- list(country=country,
                                  region=eppName,
                                  anc.used=anc.used,
                                  anc.prev=anc.prev,
                                  anc.n=anc.n,
                                  hhs=hhs)
    }
  }


  class(epp.data) <- "eppd"

  return(epp.data)
}


################################################################################
####  Function to read subpopulation sizes used in EPP fitting (from .xml)  ####
################################################################################

ind.read.epp.subpops <- function(epp.xml){

  obj <- xmlTreeParse(epp.xml)
  r <- xmlRoot(obj)[[1]]
  eppSetChildren.idx <- which(xmlSApply(r, xmlAttrs) == "eppSetChildren")
  country <- xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetCountry")]][[1]])

  epp.pops <- list() # declare list to store output
  attr(epp.pops, "country") <- country

  workset.startyear <- as.integer(xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetStartYear")]][[1]]))
  workset.endyear <- as.integer(xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetEndYear")]][[1]]))
  ## workset.popbaseyear <- as.integer(xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetPopBaseYear")]][[1]])) # not sure what this is...

  pop15to49.idx <- which(xmlSApply(r, xmlAttrs) == "pop15to49")
  pop15.idx <- which(xmlSApply(r, xmlAttrs) == "pop15")
  pop50.idx <- which(xmlSApply(r, xmlAttrs) == "pop50")
  netMigration.idx <- which(xmlSApply(r, xmlAttrs) == "netMigration")

  epp.pops$total <-  data.frame(year = workset.startyear:workset.endyear,
                                pop15to49 = 0,
                                pop15 = 0,
                                pop50 = 0,
                                netmigr = 0)
  epp.pops$total$pop15to49[as.integer(xmlSApply(r[[pop15to49.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(r[[pop15to49.idx]][[1]], xmlSApply, xmlToList))
  epp.pops$total$pop15[as.integer(xmlSApply(r[[pop15.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(r[[pop15.idx]][[1]], xmlSApply, xmlToList))
  epp.pops$total$pop50[as.integer(xmlSApply(r[[pop50.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(r[[pop50.idx]][[1]], xmlSApply, xmlToList))
  epp.pops$total$netmigr[as.integer(xmlSApply(r[[netMigration.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(r[[netMigration.idx]][[1]], xmlSApply, xmlToList))

  epp.pops$subpops <- list()
  used.indices <- NULL
  for(eppSet.idx in 1:xmlSize(r[[eppSetChildren.idx]])){

    tmp.eppSet <- r[[eppSetChildren.idx]][[eppSet.idx]][[1]]
    n.iter <- 1
    if (xmlSize(xmlSApply(tmp.eppSet, xmlAttrs)) > 0 & xmlSize(which(xmlSApply(tmp.eppSet, xmlAttrs) == "siteNames")) == 0) {
      newSetChildren.idx <- which(xmlSApply(tmp.eppSet, xmlAttrs) == "eppSetChildren")
      n.iter <- xmlSize(tmp.eppSet[[newSetChildren.idx]])
    }
    for (eppSet.idx.2 in 1:n.iter) {
      eppSet <- tmp.eppSet
      if (xmlSize(xmlSApply(tmp.eppSet, xmlAttrs)) > 0 & xmlSize(which(xmlSApply(tmp.eppSet, xmlAttrs) == "siteNames")) == 0)
        eppSet <- tmp.eppSet[[newSetChildren.idx]][[eppSet.idx.2]][[1]]

      if (xmlSize(eppSet) == 0) {

        tmp_len <- 0
        tmp_index <- 0
        obj_size <- 0
        while(obj_size == 0) {
          tmp_index <- tmp_index + 1
          if (tmp_index != eppSet.idx & !(tmp_index %in% used.indices)) {
            # print(paste('tmp_index',tmp_index))
            eppSet <- r[[eppSetChildren.idx]][[tmp_index]][[1]]
            tmp_len <- length(which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo"))
            if (tmp_len > 0) {
              tmp_eppSet <- eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo")]]
              set_location <- ifelse(xmlSize(tmp_eppSet) == 1, 1, which(xmlSApply(tmp_eppSet, xmlAttrs) == "epp2011.core.sets.ProjectionSet"))
              tmp_eppSet <- tmp_eppSet[[set_location]]
              obj_size <- xmlSize(tmp_eppSet)
            }
          }
        }
        used.indices <- c(used.indices, tmp_index)

        eppSet <- eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo")]]
        set_location <- ifelse(xmlSize(eppSet) == 1, 1, which(xmlSApply(eppSet, xmlAttrs) == "epp2011.core.sets.ProjectionSet"))
        eppSet <- eppSet[[set_location]]
      }

      eppName <- xmlToList(eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "name")]][["string"]])
      # print(eppName)

      pop15to49.idx <- which(xmlSApply(eppSet, xmlAttrs) == "pop15to49")
      pop15.idx <- which(xmlSApply(eppSet, xmlAttrs) == "pop15")
      pop50.idx <- which(xmlSApply(eppSet, xmlAttrs) == "pop50")
      netMigration.idx <- which(xmlSApply(eppSet, xmlAttrs) == "netMigration")

      subp <- data.frame(year = workset.startyear:workset.endyear,
                         pop15to49 = 0,
                         pop15 = 0,
                         pop50 = 0,
                         netmigr = 0)
      subp$pop15to49[as.integer(xmlSApply(eppSet[[pop15to49.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[pop15to49.idx]][[1]], xmlSApply, xmlToList))
      subp$pop15[as.integer(xmlSApply(eppSet[[pop15.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[pop15.idx]][[1]], xmlSApply, xmlToList))
      subp$pop50[as.integer(xmlSApply(eppSet[[pop50.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[pop50.idx]][[1]], xmlSApply, xmlToList))
      subp$netmigr[as.integer(xmlSApply(eppSet[[netMigration.idx]][[1]], xmlAttrs))+1] <- as.numeric(xmlSApply(eppSet[[netMigration.idx]][[1]], xmlSApply, xmlToList))

      epp.pops$subpops[[eppName]] <- subp
    }
  }
  class(epp.pops) <- "eppsubp"

  return(epp.pops)
}

#######################################################
#### Function to extend India by x number of years ####
######################################################

extend_india <- function(ind_totals){

  sub.art.cov.path <- paste0("<<<< FILEPATH REDACTED >>>>")
  sub.art.cov.path.update <- paste0("<<<< FILEPATH REDACTED >>>>")

  eppd <- ind_totals$eppd.ind
  epp.subp <- ind_totals$epp.subp.ind
  epp.input <- ind_totals$epp.input.ind
  subpops <- names(eppd)
  pop.cols <- colnames(epp.input$epp.pop)[!colnames(epp.input$epp.pop) %in% "year"]

  ##Function to decide between linear, exponential or no extrapolation for populations
  get_best_adjustment <- function(dt.ext,outcome,ds){
      nostop <- TRUE

      m1 <- lm(get(outcome) ~ year, data = ds)
      m2 <- lm(get(outcome) ~ year + I(year^2), data = ds)

      if(summary(m2)$coefficients[3,4] <= 0.05){
        dt.ext[,outcome] <- predict(m2,newdata=dt.ext)
        predict(m2)
        nostop <- FALSE
      }

      if(summary(m1)$coefficients[2,4] <= 0.05 & nostop==TRUE){
        dt.ext[,outcome] <- predict(m1,newdata=dt.ext)

      }

      return(dt.ext)

    }

  ##Iterate through the subpopulations and add empty columns for additional years  for ANC data
  for (subpop in subpops){
      subpop.dt.prev = eppd[[subpop]]$anc.prev
      subpop.dt.n = eppd[[subpop]]$anc.n
      years.original <- as.integer(colnames(subpop.dt.prev))

        if (!floor(proj.end) %in% years.original){
            add.years <- seq(max(years.original)+1,floor(proj.end))
            col.names <- add.years
            add.years <- matrix(,nrow(subpop.dt.prev), length(add.years),)
            colnames(add.years) <- col.names
            anc.prev <- cbind(subpop.dt.prev ,add.years)
            anc.n <- cbind(subpop.dt.n,add.years)

        }

      eppd[[subpop]]$anc.prev <-  anc.prev
      eppd[[subpop]]$anc.n <- anc.n

    }


  if(!max(years.original) %in% floor(proj.end)){
    for (subpop in subpops){
        subpop.pop = epp.subp$subpops[[subpop]]
        ds <- subpop.pop[subpop.pop$year > floor(proj.end)-8,]
               dt.ext <- data.frame()

               for(yeardiff in seq(max(subpop.pop$year)+1,floor(proj.end))){
                 ext.base <- subpop.pop[subpop.pop$year == max(subpop.pop$year),]
                 extension <- ext.base
                 extension$year <- yeardiff
                 dt.ext <- rbind(dt.ext,extension)
               }

               for(outcome in colnames(subpop.pop)[!colnames(subpop.pop) %in% "year"]){
                 dt.ext[,outcome] <- get_best_adjustment(dt.ext=dt.ext,outcome=outcome,ds=ds)[,outcome]
               }

               epp.subp$subpops[[subpop]] <- rbind(subpop.pop, dt.ext)

    }


    epp.pop.total = epp.subp$total
    ds <- epp.pop.total[epp.pop.total$year > floor(proj.end)-8,]
    dt.ext <- data.frame()
    for(yeardiff in seq(max(epp.pop.total$year)+1,floor(proj.end))){
        ext.base <- epp.pop.total[epp.pop.total$year == max(epp.pop.total$year),]
        extension <- ext.base
        extension$year <- yeardiff
        dt.ext <- rbind(dt.ext,extension)
    }

    for(outcome in pop.cols){
      dt.ext[,outcome] <- get_best_adjustment(dt.ext=dt.ext,outcome=outcome,ds=ds)[,outcome]
    }

    epp.subp$total <- rbind(epp.pop.total, dt.ext)


    epp.input.pop = epp.input$epp.pop
    ds <- epp.input.pop[epp.input.pop$year > floor(proj.end)-8,]
    dt.ext <- data.frame()
    for(yeardiff in seq(max(epp.input.pop$year)+1,floor(proj.end))){
      ext.base <- epp.input.pop[epp.input.pop$year == max(epp.input.pop$year),]
      extension <- ext.base
      extension$year <- yeardiff
      dt.ext <- rbind(dt.ext,extension)
    }

    for(outcome in pop.cols){
      dt.ext[,outcome] <- get_best_adjustment(dt.ext=dt.ext,outcome=outcome,ds=ds)[,outcome]
    }

    epp.input$epp.pop <- rbind(epp.input$epp.pop,dt.ext)


}

##Input object and ART coverage by sex/year
  epp.input$stop.year <- stop.year

  ##Extend m.val and f.val if necessary using an exponential trend
  if(!stop.year %in% epp.input$epp.art$year) {
    dt.ext <- data.frame()
    epp.art = epp.input$epp.art
    for(yeardiff in seq(max(epp.art$year)+1,stop.year)){
      ext.base <- epp.art[epp.art$year == max(epp.art$year),]
      extension <- ext.base
      extension$year <- yeardiff
      dt.ext <- rbind(dt.ext,extension)
    }

    m1.m <- lm(m.val ~ year , data = epp.input$epp.art[epp.input$epp.art$year > stop.year-5,])
    dt.ext$m.val <- predict(m1.m, newdata = dt.ext)
    m1.f <- lm(f.val ~ year  , data = epp.input$epp.art[epp.input$epp.art$year > stop.year-5,])
    dt.ext$f.val <- predict(m1.f, newdata = dt.ext)
    epp.input$epp.art <- rbind(epp.input$epp.art,dt.ext)

}


  ind_totals$eppd.ind <- eppd
  ind_totals$epp.subp.ind <- epp.subp
  ind_totals$epp.input.ind <- epp.input


  return(ind_totals)

}





###################
####  Example  ####
###################

## sa.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/South Africa 2014/SouthAfrica_p_2014June10-H-c"
## sa.epp.input <- read.epp.input(sa.path)
## sa.eppd <- read.epp.data(paste(sa.path, ".xml", sep=""))
## sa.eppsubp <- read.epp.subpops(paste(sa.path, ".xml", sep=""))
