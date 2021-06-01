################################################################################
## Purpose: This function checks the input files that are not part of the PJNZ
##          for internal consistency and that they are up to date.  It also prepares
##          some plots of these inputs to help with troubleshooting later.
################################################################################

check_epp_inputs <- function(config){

  ## shapes
  shp_version <- config$Value[which(config$Setting == "shapefile_version")]
  shp <- st_read(paste0("<<<< FILEPATH REDACTED>>>>"))
  shp$geometry <- NULL

  ## loc_list
  ad0_codes <- get_adm0_codes(loc.list)
  shp <- shp[which(shp$ADM0_CODE %in% ad0_codes), ]
  ad2_list <- shp$ADM2_CODE

  ## year_list
  start_year <- as.numeric(as.character(config$Value[which(config$Setting == "start.year")]))
  stop_year <- as.numeric(as.character(config$Value[which(config$Setting == "stop.year")]))
  year_list <- c(start_year:stop_year)

  ## admin 2 pops
  gbd_info_pops <- paste0("<<<< FILEPATH REDACTED>>>>")
  gbd_info_pops <- fread(gbd_info_pops)
  pops <- list.files(paste0("<<<< FILEPATH REDACTED>>>>"))
  for (loc in loc.list) {
    p <- fread(paste0("<<<< FILEPATH REDACTED>>>>"))
    ad2s <- shp$geo_id[which(shp$ADM0_CODE == get_adm0_codes(loc))]
    problem_years <- as.list(year_list[which(year_list %!in% p$year_id)])
    if (length(problem_years) == 0) {print(paste0(loc, " admin 2 populations cover all nessecary years"))} else {print(paste0(loc, " admin 2 populations are missing ", problem_years))
                                                                                                                  break}
    problem_admins <- as.list(ad2s[which(ad2s %!in% p$ADM2_CODE)])
    if (length(problem_admins) == 0) {print(paste0(loc, " admin 2 populations cover all nessecary admin 2 units"))} else {print(paste0(loc, " admin 2 populations are missing in ", paste0(problem_admins), " admin 2 units"))
                                                                                                                          break}
  }

  ## admin 2 prev
  prev_date <- config$Value[which(config$Setting == "prev_date")]
  prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  prev0 <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  problem_adm0 <- as.list(as.character(ad0_codes[which(ad0_codes %!in% prev0$ADM0_CODE)]))
  if (length(problem_adm0) == 0) {print(paste0("All admin 0 prevalence is present"))} else {print(paste0(" missing admin 0 prevalence in ", problem_adm0))
                                                                                              break}

  if (length(problem_adm0) != 0) {
  for (i in 1:length(problem_adm0)) {
    if (nchar(problem_adm0[[i]]) == 2) {problem_adm0[[i]] <- paste0("0", problem_adm0[[i]])}
  }}


  problem_adm2 <- as.list(as.character(ad2_list[which(ad2_list %!in% prev$ADM2_CODE)])) # this would mean that there is a miss match with the shapefile versions
  na_prev_admin2 <- as.list(as.character(prev$ADM2_CODE[which(is.na(prev$prev))])) # this would mean there is an issue in the mbg modeling process and it is outputting an NA prevalence for a given admin 2 unit.
  problem_adm2 <- c(problem_adm2, na_prev_admin2)
  problem_adm2 <- problem_adm2[which(substr(problem_adm2, (nchar(problem_adm2) - 2), nchar(problem_adm2)) %!in% problem_adm0)]
  if (length(problem_adm2) == 0) {print(paste0("All admin 2 prevalence is present in places where we have admin 0 prevalence"))} else {print(paste0(" missing admin 2 prevalence in ", unique(problem_adm2)))
                                                                                                                                        break}

  #admin2 art, years and shapes and consistent with the correct prev run
  prev_date <- config$Value[which(config$Setting == "prev_date")]
  art_type <- as.character(config$Value[which(config$Setting == "art_type")])
  if (art_type == "Coverage") {
  art_cov <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))



  problem_adm2_cov <- as.character(ad2_list[which(ad2_list %!in% art_cov$ADM2_CODE)]) # this would mean that there is a miss match with the shapefile versions
  na_prev_admin2_cov <- as.character(art_cov$ADM2_CODE[which(art_cov$prev_prop == 0)]) # this would mean there is an issue in the mbg modeling process and it is outputting an NA prevalence for a given admin 2 unit.
  problem_adm2_cov <- c(problem_adm2_cov, na_prev_admin2_cov)
  problem_adm2_cov <- problem_adm2_cov[which(substr(problem_adm2_cov, (nchar(problem_adm2_cov) - 2), nchar(problem_adm2_cov)) %!in% problem_adm0)]
  if (length(problem_adm2_cov) == 0) {print(paste0("All admin 2 art proportions for equal coverage are present in places where we have admin 0 prevalence"))} else {print(paste0("missing admin 2 art equal coverage info in ", unique(problem_adm2_cov)))
                                                                                                                                                                    break}
  }

  if (art_type == "Modeled") {
    art_mod <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

    problem_adm2_mod <- as.character(ad2_list[which(ad2_list %!in% art_mod$ADM2_CODE)]) # this would mean that there is a miss match with the shapefile versions
    problem_adm2_mod <- problem_adm2_mod[which(substr(problem_adm2_mod, (nchar(problem_adm2_mod) - 2), nchar(problem_adm2_mod)) %!in% problem_adm0)]
    if (length(problem_adm2_mod) == 0) {print(paste0("All admin 2 art proportions for modeled coverage are present in places where we have admin 0 prevalence"))} else {print(paste0("missing admin 2 art modeled coverage info in ", unique(problem_adm2_mod)))
                                                                                                                                                                        break}
  }

  #gbd_migration
  gbd_mig_version <- config$Value[which(config$Setting == "gbd_mig_run")]
  nat_migs <- fread(paste0("<<<< FILEPATH REDACTED>>>>"))
  problem_nats <- loc.list[which(loc.list %!in% nat_migs$ihme_loc_id)]
  problem_years <- year_list[which(year_list %!in% nat_migs$year_id)]
  if (length(problem_nats) == 0) {print("all top level models have GBD migration estimates")} else {print(paste0(problem_nats, " is missing a GBD migration estimate"))
                                                                                                    break}
  if (length(problem_years) == 0) {print("all years have GBD migration estimates")} else {print(paste0(problem_years, " is missing a GBD migration estimate"))
                                                                                          break}

  #gbd_migration
  gbd_mort_version <- config$Value[which(config$Setting == "gbd_mort_run")]
  nat_morts <- readRDS(paste0("<<<< FILEPATH REDACTED>>>>"))
  problem_nats <- loc.list[which(loc.list %!in% nat_morts$ihme_loc_id)]
  problem_years <- year_list[which(year_list %!in% nat_morts$year_id)]
  if (length(problem_nats) == 0) {print("all locs have GBD HIV free mortality estimates")} else {print(paste0(problem_nats, " is missing a GBD HIV free mortality estimate"))
                                                                                                  break}
  if (length(problem_years) == 0) {print("all years have GBD HIV free mortality estimates")} else {print(paste0(problem_years, " is missing a GBD HIV free mortality estimate"))
                                                                                                    break}

}
























