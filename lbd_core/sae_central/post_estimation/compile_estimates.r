####################################################################################################
## Description: Compile point estimates, CI, and standard errors for all sexes, years, and
##              geographic levels. This code also runs basic consistency checks on output -- if any
##              of these checks fail, a report about the identified issues is saved instead of the
##              final output file.
##
## Passed args: main_dir [character] -- home directory for settings and final output
##
## Requires:    point estimates and CI ('[temp_dir]/est_[level]_[year]_[sex].rdata')
##
## Outputs:     If all checks pass:
##                compiled estimates ('[main_dir]/est_all.rdata')
##
##              Otherwise:
##                issues in the estimates ('[main_dir]/est_all_issues.txt')
##
## Run from within 'sae_central' directory!!
####################################################################################################

library(data.table)

rm(list=ls())

## Get settings ------------------------------------------------------------------------------------
main_dir <- commandArgs()[4]

source("settings.r")
get_settings(main_dir)

# Function to check estimates
check <- function(est) {
  issues <- NULL

# check: all estimates are >= 0
  if (est[, sum(mean < 0 | lb < 0 | ub < 0)] > 0) issues <- c(issues, paste("Not all estimates are positive"))

# check: point estimates within confidence intervals
  if (est[, sum(!between(mean, lb, ub))] > 0) issues <- c(issues, paste("Point estimates are outside confidence bounds"))

# check: there are no duplicates
  if (nrow(est) != nrow(unique(est[, list(level, area, year, sex, age)]))) issues <- c(issues, "There are duplicated rows")

# check: all expected level/area/year/sex/age combinations are present
  all <- CJ(level = c(area_var, names(geoagg_files)), sex = c(sexes, 3), year = years, age = c(fread(age_std_file)$age, 98, 99))
  all <- merge(all, unique(est[, list(level, area)]), by = "level", allow.cartesian = T)
  if (!is.null(geoagg_files)) {
    for (this_level in names(geoagg_files)) {
      load(geoagg_files[this_level])
      all <- all[level != this_level | year %in% unique(weights$year),]
      rm(weights)
    }
  }
  setkeyv(all, key(est))
  if (class(all.equal(all, est[, names(all), with=F])) == "character") issues <- c(issues, "There are missing level/area/year/sex/age combinations")
  rm(all)

# check: value for both sexes combined is between value for males and females
  if (length(sexes) == 2) {
    temp <- est[age != 99, list(level, area, year, sex, age, value=mean)]
    temp <- dcast.data.table(temp, level + area + year + age ~ sex, value.var="value")
    setnames(temp, as.character(1:3), paste0("value_", 1:3))
    temp <- na.omit(temp)
    if (temp[, mean(between(value_3, 0.99999*pmin(value_1, value_2), 1.00001*pmax(value_1, value_2))) < 1]) issues <- c(issues, "estimates for both sexes combined is not always between the male and female estimates")
    rm(temp)
  }

  issues
}

## Compile estimates -------------------------------------------------------------------------------
# load and combine estimates for all years, sexes, and levels
files <- dir(temp_dir, pattern = "est_")
est <- rbindlist(lapply(files, function(f) {
  load(paste0(temp_dir, "/", f))
  est
}))

# format
setkeyv(est, c("level", "area", "year", "sex", "age"))
setcolorder(est, c(key(est), "mean", "lb", "ub", "se"))

# check for issues
issues <- check(est)

# if there are no problems, save the estimates, otherwise save the list of issues
file.remove(paste0(main_dir, "/est_all_issues.csv"), showWarnings = F) # in case this a resub, remove any old issues file
if (length(issues) == 0) save(est, file = paste0(main_dir, "/est_all.rdata"))
if (length(issues) > 0) write.table(issues, file = paste0(main_dir, "/est_all_issues.csv"), row.names = F, col.names = F)
