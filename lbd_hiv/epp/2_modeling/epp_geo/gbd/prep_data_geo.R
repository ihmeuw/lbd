################################################################################
## Purpose: This is the set of main wrapper functions that collect the main epp
##          data inputs and put them into the format need to run epp.  There are
##          a few other substitutions that can be specified by other functions
##          in the sub_data* scripts but those assume you are opperating on the
##          output of the prep_epp_data() function.  That function should repurn
##          a list which has the nessecary data stored as attributes.
################################################################################

#' @title prep_epp_data
#' @description This is the main data packaging function for EPP.  It controls the reading in and formatting of constituent parts of the EPP
#' model object that is ues to run the EPP model
#'
#' @param loc the location id for the GBD geography being modeled.  This is the ISO 3 code for national models and the iso3_somenumber
#' for the GBD subnational models
#' @param proj.end the decimal year of the last modeled time step
#' @param stop_collapse Boolian indicating if the data for geographic subpopulations in the SPECTRUm files should be compbined or maintainted
#' separately.
#' @param gbd_pop a boolian indicating if EPP should use the SPECTRUM populations 'FALSE' or the GBD populations 'TRUE'
#' @param art_sub a boolian indicating if EPP should use the native SPECTRUM ART data 'FALSE' or the GBD future projected data 'TRUE'.  To be
#' clear, this only impacts ART in the later years when both the SPECTRUM file and GBD are using projected ART numbers, not actual programatic
#' data
#' @param num_knots the number of knots in the spline parameterization of the force of infection vector
#' @param no.anc boolian indicating if ANC data is used 'FALSE' or not used 'TRUE' in the EPP fitting process
#' @param run.name the run name that us being prepared.  this is important because some of the early intermediary files are saved for later use
#' in preparing the admin 2 EPP objects
#' @param gbd.mort boolian indicating if HIV free mortality should be calculated from GBD estimates 'TRUE' or as a remainder of population
#' differences 'FALSE' this second option is what was originally programed into EPP
#' @param mig bollian indicating if the migration should come from GBD 'TRUE' or SPECTRUM 'FALSE'  This is possibly vestigal in that we have
#' migration set to 0 later in the data processing.  This is becasue estimateing it for admin 2 was very difficult so we went the route of
#' the popadjust setting.
#' @param shapefile_version is the version of the admin 2 shapes that we are going to be using to build out the geopgarphic aspects of the
#' models. This relates to the prevalence inputs as well.
#'
#' @return a list with one geographic subpopulation in it.  the actuall EPP data objects are attributes to that list.
#'
#' @export



prep_epp_data <- function(loc, proj.end = 2019.5,
                          stop_collapse = FALSE, gbd_pop = TRUE, art_sub = TRUE, num_knots = 7,
                          no.anc, run.name, i = 1, gbd.mort = gbd.mort, mig = mig, shapefile_version = shapefile_version, gbd_mort_run = gbd_mort_run, mort_draw = mort_draw) {

    if (stop_collapse | grepl("NGA_", loc) | grepl("ZAF_", loc) | grepl("IND_", loc)) {
        collapse <- F
    } else {
        collapse <- T#loc.table[ihme_loc_id == loc, collapse_subpop]
    }
    if (collapse) {
        print("Collapsing subpopulations")
        epp_totals <- collapse_epp(loc)
        saveRDS(epp_totals, file = paste0("<<<< FILEPATH REDACTED >>>>"))
        eppd <- epp_totals$eppd.tot
        epp.subp <- epp_totals$epp.subp.tot
        if (gbd_pop) {
            epp.subp <- sub.pop.params(epp.subp, loc, mig = mig)
        }
        epp.input <- epp_totals$epp.input.tot
        if (art_sub) {
            epp.input <- sub.art(epp.input, loc)
        }
        if (loc == "MRT") {
          epp.input$epp.art$art15yr[which(epp.input$epp.art$year == 2001)] <- epp.input$epp.art$art15yr[which(epp.input$epp.art$year == 2000)] #correcting an error in the MRT 2019 spectrum file.  The MRT 2017 file has a 0 here which is more consistent with the overall time trend.
        }
        if (loc == "NGA") {
          epp.input$epp.art$art15yr[which(epp.input$epp.art$year == 2000)] <- epp.input$epp.art$art15yr[which(epp.input$epp.art$year == 1999)] #correcting an error in the NGA 2019 spectrum file.  The NGA 2017 file has smoother time trend here that can be aproximated by copying the value from 1999.
          epp.input$epp.art$art15yr[which(epp.input$epp.art$year == 2001)] <- epp.input$epp.art$art15yr[which(epp.input$epp.art$year == 1999)]
          epp.input$epp.art$art15yr[which(epp.input$epp.art$year == 2002)] <- epp.input$epp.art$art15yr[which(epp.input$epp.art$year == 1999)]
        }


        epp_totals$eppd.tot <- eppd
        epp_totals$epp.subp.tot <- epp.subp
        epp_totals$epp.input.tot <- epp.input
        epp_totals$mort_draw <- mort_draw
        saveRDS(epp_totals, file = paste0("<<<< FILEPATH REDACTED >>>>"))
        epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd)
        epp.subp.input[[loc]]$reg <- loc
        val <- setNames(vector("list", length(eppd)), names(eppd))
        set.list.attr <- function(obj, attrib, value.lst) mapply(function(set, value) {
            attributes(set)[[attrib]] <- value
            set
        }, obj, value.lst)

        #Adding the use case for extracting the EPP objects with no anc data
        if (no.anc == TRUE & exists("anc.used", where = eppd[[loc]])) {
        for (counter in 1:length(eppd[[loc]]$anc.used)) {eppd[[loc]]$anc.used[counter] <- FALSE}
        eppd[[loc]]$anc.n <- NULL
        eppd[[loc]]$anc.prev <- NULL}

        val <- set.list.attr(val, "eppd", eppd)
        val <- set.list.attr(val, "likdat", lapply(eppd, fnCreateLikDat, anchor.year = epp.input$start.year, no.anc = no.anc))
        val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar, proj.end = proj.end, num.knots = num_knots, no.anc = no.anc, gbd.mort = gbd.mort, loc = loc, subnat = "national", shapefile_version = shapefile_version, gbd_mort_version = gbd_mort_run))
        val <- set.list.attr(val, "country", attr(eppd, "country"))
        val <- set.list.attr(val, "region", names(eppd))
        val <- set.list.attr(val, "mort_draw", mort_draw)

        prepped.dt <- val

       } else {
        ## Get location of file
        # Kenya counties
        if (grepl("KEN", loc) & loc.table[ihme_loc_id == loc, level] == 5) {
            if (loc.table[ihme_loc_id == loc, parent_id] == 44797) {
                temp.loc <- "KEN_44795" # This is because the North Eastern province doesn't have an XML file
            } else {
               temp.loc <- loc.table[location_id == loc.table[ihme_loc_id == loc, parent_id], ihme_loc_id]
            }
        } else {
            temp.loc <- loc
        }


        for (c.year in c('2018','2017','2016','2015','2013')) {
          if (c.year == 2016 | c.year == 2017 | c.year == 2018) {
            dir <- paste0("<<<<FILEPATH REDACTED>>>>>")
          } else {
            dir <- paste0("<<<<FILEPATH REDACTED>>>>>")
          }

        if (file.exists(paste0(dir, temp.loc, ".PJNZ")))
        break;
        }

        ## South Africa location map
        if (grepl("ZAF", loc)) {
          c.year <- 2017#loc.table[ihme_loc_id == temp.loc, unaids_recent]
          if (c.year == 2016 | c.year == 2017 | c.year == 2018) {
            dir <- paste0("<<<<FILEPATH REDACTED>>>>>")
          } else {
            dir <- paste0("<<<<FILEPATH REDACTED>>>>>")
          }
            zaf.dict <- list("MP" = "ZAF_487", "GP" = "ZAF_484", "KZN" = "ZAF_485",
                             "WC" = "ZAF_490", "EC" = "ZAF_482", "LP" = "ZAF_486",
                             "FS" = "ZAF_483", "NW" = "ZAF_488", "NC" = "ZAF_489")
            pjnz <- paste0(dir, "ZAF.PJNZ")
            nat.dt <- prepare_epp_fit(pjnz, proj.end = proj.end, num_knots = num_knots, no.anc = no.anc)
            prepped.dt <- list()
            prepped.dt[[loc]] <- nat.dt[[names(which(zaf.dict == loc))]]
        ## Special case for India
        } else if (grepl("IND", loc)) {
            filepath <- paste0(dir, loc)
            ind_totals_no_extension <- ind.prepare.epp.fit(filepath, proj.end = proj.end)
            ind_totals <- extend_india(ind_totals = ind_totals_no_extension )
            eppd <- ind_totals$eppd.ind
            epp.subp <- ind_totals$epp.subp.ind
            epp.input <- ind_totals$epp.input.ind

            ## Fill in missing data (epp_totals vs. ind_totals)
            epp.input$epidemic.start <- 1985
            pop.dim <- length(epp.input$epp.pop$year)

            # epp.pop
            epp.input$epp.pop$cd4median <- rep(0, times = pop.dim)
            epp.input$epp.pop$hivp15yr <- rep(0, times = pop.dim) # can pull from spectrum estimates (prevelance among 15 year olds)

            # Expand progression parameters to matrix
            prog.names <- c("cd4stage.dur", "cd4mort", "artmort.less6mos", "artmort.6to12mos", "artmort.after1yr")
            for (param in prog.names) {
            epp.input[[param]] <- matrix(rep(epp.input[[param]], times = 8), nrow = 8, byrow = T)
            }

            # cd4initperc
            std.cd4initperc <- matrix(data = c(64.3, 35.7, 0, 0, 0, 0, 0, 60.7, 39.3, 0, 0, 0, 0, 0, 58.5, 41.5, 0, 0, 0, 0, 0, 55.2, 44.8, 0, 0, 0, 0, 0, 64.3, 35.7, 0, 0, 0, 0, 0, 60.7, 39.3, 0, 0, 0, 0, 0, 58.5, 41.5, 0, 0, 0, 0, 0, 55.2, 44.8, 0, 0, 0, 0, 0),
                                        nrow = 8, ncol = 7, byrow = T)
            row.names(std.cd4initperc) <- c("NEWINFECTSCD4_M_15_24", "NEWINFECTSCD4_M_25_34", "NEWINFECTSCD4_M_35_44", "NEWINFECTSCD4_M_45_54", "NEWINFECTSCD4_F_15_24", "NEWINFECTSCD4_F_25_34", "NEWINFECTSCD4_F_35_44", "NEWINFECTSCD4_F_45_54")
            colnames(std.cd4initperc) <- paste0("V", 2:8)
            epp.input$cd4initperc <- std.cd4initperc

            # epp.art
            art.dim <- nrow(epp.input$epp.art)
            epp.input$epp.art$'1stto2ndline'  <- rep(0, times = art.dim)
            epp.input$epp.art$art15yr  <- rep(0, times = art.dim)
            epp.input$epp.art$artdropout  <- rep(0, times = art.dim)

            # art.specpop
            specpop.sub <- data.frame(specpop = c("PW", "TBHIV", "DC", "FSW", "MSM"),
                                      percelig = c(0, 0, 0, 0, 0),
                                      yearelig = c(2015, 2015, 2015, 2015, 2015))
            epp.input$art.specpop <- specpop.sub

            # hivp15yr.cd4dist
            epp.input$hivp15yr.cd4dist <- c(0.056, 0.112, 0.112, 0.070, 0.140, 0.230, 0.280)

            # art15yr.cd4dist
            epp.input$art15yr.cd4dist <- c(0.00, 0.00, 0.11, 0.23, 0.23, 0.14, 0.29)

            epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd)
            val <- setNames(vector("list", length(eppd)), names(eppd))
            set.list.attr <- function(obj, attrib, value.lst) mapply(function(set, value) {
                attributes(set)[[attrib]] <- value
                set
            }, obj, value.lst)


            val <- set.list.attr(val, "eppd", eppd)
            val <- set.list.attr(val, "likdat", lapply(eppd, fnCreateLikDat, anchor.year = epp.input$start.year, no.anc = no.anc))
            val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar, proj.end = proj.end, num.knots = num_knots, no.anc = no.anc))
            val <- set.list.attr(val, "country", attr(eppd, "country"))
            val <- set.list.attr(val, "region", names(eppd))

            prepped.dt <- val
        } else {

            pjnz <- paste0(dir, temp.loc, ".PJNZ")
            prepped.dt <- prepare_epp_fit(pjnz, proj.end = proj.end, num_knots = num_knots, no.anc = no.anc)
        }
        ## Screen ANC data for Kenya
        if (grepl("KEN", loc) & loc.table[ihme_loc_id == loc, level] == 5) {
            # ANC
            ken.anc.path <- paste0("<<<<FILEPATH REDACTED>>>>>")
            ken.anc <- fread(ken.anc.path)
            county.sites <- ken.anc[ihme_loc_id == loc, site]
            prov.sites <- row.names(attr(prepped.dt[[1]], "eppd")$anc.prev)
            keep.index <- which(prov.sites %in% county.sites)
            attr(prepped.dt[[1]], "eppd")$anc.used[] <- FALSE
            if (length(keep.index) > 0) {
                attr(prepped.dt[[1]], "eppd")$anc.used[keep.index] <- TRUE
            }

            # ART
            prop.path <- paste0("<<<<FILEPATH REDACTED>>>>>")
            prop.dt <- fread(prop.path)
            prop <- prop.dt[ihme_loc_id == loc, prop_pepfar]
            attr(prepped.dt[[1]], "eppfp")$artnum.ts <- attr(prepped.dt[[1]], "eppfp")$artnum.ts * prop
        }
    }
    ## Sample from multiple start years
    range <- 6
    for (subpop in names(prepped.dt)) {
        epi.start <- attr(prepped.dt[[subpop]], "eppfp")$tsEpidemicStart
        attr(prepped.dt[[subpop]], "eppfp")$tsEpidemicStart <- max(1970.5, sample(seq((epi.start - (range / 2)),(epi.start + (range / 2)), 0.1), 1))
    }
    return(prepped.dt)
}




#' @title mode
#' @description helper function that finds the mode of a list
#'
#' @param x some sort of vector
#'
#' @return the most frequent value in x
#'
#' @export


Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

#' @title collapse_epp
#' @description helper function does the work of reading in the data from the PJNZ files or file for agiven location.
#'
#' @param loc the location whose PJNZ file or files are to be read.
#'
#' @return a list with 3 objects, eppd.tot, epp.subp.tot, and epp.input.tot
#'
#' @export


collapse_epp <- function(loc) {
    pjnz.list <- c()

    for (c.year in c('2019', '2018','2017','2016','2015','2013')) {
      if ( loc == "ZAF") { c.year <- 2017 } # there are issues with the 2018 ZAF file and there is no 2019 ZAF file
      if ( loc == "TGO") { c.year <- 2018 } # there are issues with the 2019 TGO file.
      print(paste0("Checking UNAIDS data for ", c.year))
      if (c.year == 2016 | c.year == 2017 | c.year == 2018 | c.year == 2019) {
         dir <- paste0("<<<<FILEPATH REDACTED>>>>>")
          } else {
        dir <- paste0("<<<<FILEPATH REDACTED>>>>>")
        }
    if (dir.exists(dir)) {
        pjnz.list <- list.files(dir, pattern = "PJNZ", full.names = T)
        file.list <- grep(loc, pjnz.list, value = T)
        #if (loc == "NGA") file.list <- c()
    } else {
        one.up <- paste(head(unlist(tstrsplit(dir, "/")), -1), collapse = "/")
        dir.list <- list.files(one.up, pattern = loc, full.names = T)
        file.list <- unlist(lapply(dir.list, function(dir) {
            list.files(dir, pattern = "PJNZ", full.names = T)
        }))
    }
    if (length(file.list) == 0) {

        loc.name <- loc.table[ihme_loc_id == loc, location_name]
        loc.name <- gsub(" ", "", gsub("[^[:alnum:] ]", "", loc.name))
        file.list <- grep(loc.name, pjnz.list, value = T)
    }

    if (length(file.list) > 0)

    break;

}
  if (loc == "KEN_44797") { file.list <- paste0("<<<< FILEPATH REDACTED >>>>")
  c.year <- 2018}
  if (length(file.list) == 0) {
      print("No data files found")
    } else {
      print(paste0("Using ", c.year, " PJNZ file"))
    }


    ## eppd
    eppd.list <- lapply(file.list, function(file) {
        pjnz <- file
        eppd <- read_epp_data(pjnz, c.year)
    })

    eppd.list <- unlist(eppd.list,recursive = FALSE )

    eppd.tot <- eppd.list[1]
    subpop.tot <- loc
    names(eppd.tot) <- subpop.tot

    # region
    eppd.tot[[subpop.tot]]$region <- subpop.tot

    #country
    attr(eppd.tot,"country") <- eppd.tot[[1]]$country

    # anc.used (append)
    eppd.tot[[subpop.tot]]$anc.used <- unlist(lapply(eppd.list, function(eppd) {
        subpop <- names(eppd)
        anc.used <- eppd$anc.used
    }))

    # anc.prev (append)
    eppd.tot[[subpop.tot]]$anc.prev <- do.call(rbind, lapply(eppd.list, function(eppd) {
        subpop <- names(eppd)
        anc.prev <- eppd$anc.prev
    }))

    # anc.n (append)
    eppd.tot[[subpop.tot]]$anc.n <- do.call(rbind, lapply(eppd.list, function(eppd) {
        subpop <- names(eppd)
        anc.n <- eppd$anc.n
    }))

    # hhs (append) ** be careful "not used TRUE"
    hhs.temp <- data.table(do.call(rbind, lapply(eppd.list, function(eppd) {
        subpop <- names(eppd)
        hhs <- eppd$hhs
    })))
    hhs.temp <- hhs.temp[used == TRUE]
    hhs.temp[, pos := n * prev]
    hhs.sum <- hhs.temp[, lapply(.SD, sum), by = .(year)]
    hhs.sum[, prev := pos / n]
    hhs.sum[, se := ((prev * (1 - prev)) / n)**0.5]
    hhs.sum[, used := NULL]
    hhs.sum[, pos := NULL]
    hhs.sum[, used := TRUE]
    eppd.tot[[subpop.tot]]$hhs <- as.data.frame(hhs.sum)
    # eppd.tot[[subpop.tot]]$hhs <- as.data.frame(hhs.temp)

    ## epp.subp
    epp.subp.list <- lapply(file.list, function(file) {
        pjnz <- file
        epp.subp <- read_epp_subpops(pjnz)
    })

    epp.subp.tot <- list()


    # total
    total.temp <- data.table(do.call(rbind, lapply(epp.subp.list, function(epp.subp) {
        anc.n <- epp.subp$total
    })))
    total.sum <- total.temp[, lapply(.SD, sum), by = .(year)]
    epp.subp.tot$total <- as.data.frame(total.sum)
    epp.subp.tot$subpops[[subpop.tot]] <- as.data.frame(total.sum)

    ## epp.input
    if (length(file.list) > 1) {
        epp.input.list <- lapply(file.list, function(file) {
            pjnz <- file
            epp.input <- read_epp_input(pjnz)
        })
        epp.input.tot <- epp.input.list[[1]]
        attr(epp.input.tot,"country") <- subpop.tot

        # start.year (check for difference)
        start.years <- unlist(lapply(epp.input.list, function(epp.input) {
            start.year <- epp.input$start.year
        }))
        length(unique(start.years)) == 1

        # stop.year (check for difference)
        stop.years <- unlist(lapply(epp.input.list, function(epp.input) {
            stop.year <- epp.input$stop.year
        }))
        length(unique(stop.years)) == 1

        # epidemic.start (check for difference)
        epidemic.starts <- unlist(lapply(epp.input.list, function(epp.input) {
            epidemic.start <- epp.input$epidemic.start
        }))
        epp.input.tot$epidemic.start <- min(epidemic.starts)

        # epp.pop (sum and mean)
        epp.pop.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
            epp.pop <- epp.input$epp.pop
        })))
        epp.pop.sum <- epp.pop.temp[, lapply(.SD, sum), by = .(year)]
        epp.pop.mean <- epp.pop.temp[, lapply(.SD, mean), by = .(year)]
        epp.pop.comb <- cbind(epp.pop.sum[, .(year, pop15to49, pop15, pop50, netmigr)], epp.pop.mean[, .(cd4median, hivp15yr)])
        epp.input.tot$epp.pop <- epp.pop.comb

        # cd4lowlim (check for difference)
        cd4lowlim.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
            cd4lowlim <- epp.input$cd4lowlim
        })))

        # cd4initperc (check for difference)
        cd4initperc.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
            cd4initperc <- epp.input$cd4initperc
        })))

        # cd4stage.dur (check for difference)
        cd4stage.dur.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
            cd4stage.dur <- epp.input$cd4stage.dur
        })))

        # cd4mort, artmort.less6mos, artmort.6to12mos, artmort.after1yr (leave the same)

        # infectreduc (check for difference)
        infectreducs <- unlist(lapply(epp.input.list, function(epp.input) {
            infectreduc <- epp.input$infectreduc
        }))
        length(unique(infectreducs)) == 1

         # epp.art (sum and mean) ** beware of percentages!!! also not sure whether 1stto2ndline is count or percent
        epp.art.temp <- rbindlist(lapply(epp.input.list, function(epp.input) {
            epp.art <- epp.input$epp.art
        }), fill = T)
        epp.art.temp[is.na(m.isperc), m.isperc := "N"]
        epp.art.temp[is.na(f.isperc), f.isperc := "N"]# this is where the issues in CIV comes from.  ART data that is coverage gets labeled as numbers.
        if ("P" %in% unique(c(epp.art.temp$m.isperc, epp.art.temp$f.isperc))) {
            # Add prevalence
            epp.prev <- unlist(lapply(file.list, function(pjnz) {
              spu <- read_spu(pjnz)$prev
                mean.spu <- rowMeans(spu)
            }))
            epp.prev.subset <- epp.prev[names(epp.prev) %in% paste0(unique(epp.art.temp$year))]
            if ( loc == "ETH_44856" & loc.table$unaids_recent[which(loc.table$ihme_loc_id == loc)] == 2018) {epp.art.temp <- epp.art.temp[c(1:55), ]} #there is an error in the prevalence measures for the concentrated epidemic that makes processing the art coverage impossible, but that error only impacts the estimates for 2022

            if (nrow(epp.art.temp) != length(epp.prev.subset)) {
                stop("ART collapse problem")
            }
            epp.art.temp[, prev := epp.prev.subset]

            # Add population
            pop <- epp.pop.temp[year %in% unique(epp.art.temp$year)]
            epp.art.temp <- cbind(epp.art.temp, pop[, .(pop15to49)])
            epp.art.temp[, c("m.val", "f.val") := .(as.numeric(m.val), as.numeric(f.val))]
            epp.art.temp[m.isperc == "P", m.val := (weighted.mean(m.val, w = pop15to49 * prev) / .N), by = .(year)]
            epp.art.temp[f.isperc == "P", f.val := (weighted.mean(f.val, w = pop15to49 * prev) / .N), by = .(year)]
            epp.art.temp[, c("pop15to49", "prev") := NULL]
        }
        epp.art.hold <- epp.art.temp[1:length(min(epp.art.temp$year):max(epp.art.temp$year)), .(m.isperc, f.isperc)]
        epp.art.temp[, m.isperc := NULL]
        epp.art.temp[, f.isperc := NULL]
        epp.art.sum <- epp.art.temp[, lapply(.SD, sum), by = .(year)]
        epp.art.mean <- epp.art.temp[, lapply(.SD, mean), by = .(year)]
        epp.art.mode <- epp.art.temp[, lapply(.SD, Mode), by = .(year)]
        epp.art.comb <- cbind(epp.art.sum[, .(year, m.val, f.val, artdropout)],
                                epp.art.mode[, .(cd4thresh)],
                                epp.art.mean[, c("m.perc50plus", "f.perc50plus", "perc50plus", "1stto2ndline", "art15yr"), with = F],
                                epp.art.hold)
        epp.art.order <- epp.art.comb[, c("year", "m.isperc", "m.val", "f.isperc", "f.val", "cd4thresh", "m.perc50plus", "f.perc50plus", "perc50plus", "1stto2ndline", "art15yr"), with = F]
        epp.input.tot$epp.art <- as.data.frame(epp.art.order)

        # art.specpop (check for difference)
        art.specpop.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
            art.specpop <- epp.input$art.specpop
        })))

        # hivp15yr.cd4dist (check for difference)
        hivp15yr.cd4dist.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
            hivp15yr.cd4dist <- epp.input$hivp15yr.cd4dist
        })))

        # art15yr.cd4dist (check for difference)
        art15yr.cd4dist.temp <- data.table(do.call(rbind, lapply(epp.input.list, function(epp.input) {
            art15yr.cd4dist <- epp.input$art15yr.cd4dist
        })))
    } else {
        epp.input.tot <- read_epp_input(file.list[1])
    }


    # epidemic.type (check for difference)
    # epidemic.types <- unlist(lapply(epp.input.list, function(epp.input) {
    #   epidemic.type <- epp.input$epidemic.type
    # }))
    # length(unique(epidemic.types)) == 1
    epp_totals <- list(eppd.tot = eppd.tot, epp.subp.tot = epp.subp.tot, epp.input.tot = epp.input.tot )
    return(epp_totals)

}





####### code from EPPASM

prepare_rhybrid <- function(fp,
                            tsEpidemicStart = fp$ss$time_epi_start + 0.5,
                            rw_start = fp$rw_start,
                            rw_trans = fp$rw_trans,
                            rw_dk = fp$rw_dk){

  if (is.null(rw_start)) {rw_start <- 2003}

  if (is.null(rw_trans)) {rw_trans <- 5}

  if (is.null(rw_dk)) {rw_dk <- 5}

  fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - tsEpidemicStart))]

  rt <- list()

  rt$proj.steps <- fp$proj.steps

  rt$rw_start <- rw_start
  rt$rw_trans <- rw_trans

  switch_idx <- max(which(fp$proj.steps <= rw_start))
  rt$rlogistic_steps <- fp$proj.steps[1:switch_idx]
  rt$rw_steps <- fp$proj.steps[switch_idx:length(fp$proj.steps)]

  rt$n_rw <- ceiling((max(rt$proj.steps) - rw_start) / rw_dk)
  rt$rw_dk <- rw_dk
  rt$rw_knots <- seq(rw_start, rw_start + rt$rw_dk * rt$n_rw, by = rt$rw_dk)
  rt$rw_idx <- findInterval(rt$rw_steps[-1], rt$rw_knots)

  rt$n_param <- 4 + rt$n_rw  # 4 parameters for rlogistic

  ## Linearly interpolate between 0 and 1 over the period (rw_start, rw_start + rw_trans)
  ## Add a small value to avoid R error in approx() if rw_trans = 0
  rt$rw_transition <- approx(c(rw_start, rw_start + rw_trans + 0.001), c(0, 1), rt$rw_steps[-1], rule = 2)$y

  rt$dt <- 0.1

  rt$eppmod <- "rhybrid"
  fp$rt <- rt

  if (!exists("eppmod", fp))
    fp$eppmod <- "rhybrid"
  fp$iota <- NULL

  return(fp)
}
