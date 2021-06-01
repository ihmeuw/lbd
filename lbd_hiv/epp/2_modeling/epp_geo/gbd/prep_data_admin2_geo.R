################################################################################
## Purpose: This is the set of main wrapper functions that prepare the EPP data
##          objects for the admin 2 models.  It reads in the data from a higher
##          level PJNZ file and then replaces much of that data with data from
##          the admin 2 being modeled.
################################################################################

#' @title prep_epp_data_adm2
#' @description This is the main data packaging function for EPP at teh admin 2 level.  It controls the reading in and formatting of
#' constituent parts of the EPP model object that is ues to run the EPP model
#'
#' @param loc the location id for the GBD geography that the admin 2 being modeled is part of.  This is the ISO 3 code for national models
#' and the iso3_somenumber for the GBD subnational models
#' @param subnat This is the ADM2_CODE for the admin 2 being modeled
#' @param run.name the run name that is happening, used for retreaving the national data
#' @param i the transition parameters draw that is being modeled.
#' @param proj.end the decimal year of the last modeled time step
#' @param stop_collapse Boolian indicating if the data for geographic subpopulations in the SPECTRUm files should be compbined or maintainted
#' separately.
#' @param gbd_pop a boolian indicating if EPP should use the SPECTRUM populations 'FALSE' or the GBD populations 'TRUE'
#' @param art_sub a boolian indicating if EPP should use the native SPECTRUM ART data 'FALSE' or the GBD future projected data 'TRUE'.  To be
#' clear, this only impacts ART in the later years when both the SPECTRUM file and GBD are using projected ART numbers, not actual programatic
#' data
#' @param num_knots the number of knots in the spline parameterization of the force of infection vector
#' @param no.anc boolian indicating if ANC data is used 'FALSE' or not used 'TRUE' in the EPP fitting process
#' @param prev_date the date of the prevalence model that will be used in this admin 2 run.
#' @param migration_run A VESTIGAL SETTING SINCE WE SET ALL MIGRATION = 0
#' @param migration_type ANNOTHER VESTIGAL SETTING SINCE ALL MIGRATION IS SET TO 0
#' @param stop.year the last year of the model
#' @param art.type which art splitting method to use.
#' @param gbd.mort boolian indicating if HIV free mortality should be calculated from GBD estimates 'TRUE' or as a remainder of population
#' differences 'FALSE' this second option is what was originally programed into EPP
#' @param shapefile_version is the version of the admin 2 shapes that we are going to be using to build out the geopgarphic aspects of the
#' models. This relates to the prevalence inputs as well.
#'
#' @return a list with one geographic subpopulation in it.  the actuall EPP data objects are attributes to that list.
#'
#' @export


prep_epp_data_adm2 <- function(loc, subnat, run.name, i, proj.end = 2019.5,
                               stop_collapse = FALSE, gbd_pop = TRUE, art_sub = FALSE,
                               num_knots = 7, no.anc, prev_date = prev_date, migration_run = migration_run,
                               migration_type = migration_type, stop.year = stop.year, gbd.mort = gbd.mort,
                               art_type = art_type, shapefile_version = shapefile_version, gbd_mort_run = gbd_mort_run, likelihood_type = likelihood_type) {

    dt_raw <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
    mort_draw <- dt_raw$mort_draw
    eppd <- dt_raw$eppd.tot
    eppd <- sub.prev.params.adm2(eppd, loc, subnat, prev_date)
    epp.subp <- dt_raw$epp.subp.tot
    epp.subp <- sub.pop.params.adm2(epp.subp, loc, subnat, migration_run, migration_type, shapefile_version)
    epp.input <- dt_raw$epp.input.tot
    epp.input <- sub.art.adm2(epp.input, loc, subnat, art_dist_type = art_type, prev_date)
    epp.subp.input <- fnCreateEPPSubpops_adm2(epp.input, epp.subp, eppd, stop.year)
    epp.subp.input[[1]]$loc <- loc
    epp.subp.input[[1]]$subnat <- subnat
    epp.subp.input[[1]]$shp <- shapefile_version
    val <- setNames(vector("list", length(eppd)), subnat)
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
        val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar_adm2, proj.start = (start.year + 0.5), proj.end = (stop.year + 0.5), num.knots = num_knots, no.anc = no.anc, gbd.mort = gbd.mort, gbd_mort_run = gbd_mort_run))
        val <- set.list.attr(val, "country", attr(eppd, "country"))
        val <- set.list.attr(val, "region", names(eppd))
        val <- set.list.attr(val, "mort_draw", mort_draw)

        prepped.dt <- val
return(prepped.dt)
}

val  <- NULL

