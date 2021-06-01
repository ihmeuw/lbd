
waitformodelstofinish <- function(sleeptime=100,
                                  path =  paste0("<<<<FILEPATH REDACTED>>>>>", indicator_group, '/', indicator, '/output/', run_date),
                                  rd   = run_date,
                                  lv   = loopvars,
                                  showfiles = TRUE,
                                  showcluster = FALSE){

  n_finished <- length(grep("fin_", list.files(path)))

  lv <- data.table(lv)
  names(lv) <- c("reg", "holdout")
  lv$reg <- as.character(lv$reg)

  lv[, file := paste0(path, "/", "fin__bin0_", reg, "_", holdout)]

  while(n_finished != nrow(lv)){

    n_finished <- length(grep("fin_", list.files(path)))

    message('\n====================================================================================')
    message(sprintf('=====================      Run Date: %s      ======================',rd))
    message(paste0('\nAt ',Sys.time(),' .... ',n_finished,' Models have written output.'))
    if(showfiles){
      message('\nCurrently missing models:')
      for(i in 1:nrow(lv))
        if(file.exists(lv[i, file]) == F)
          message(paste('Region =',lv[i,1],'| Holdout =',lv[i,2]))
    }
    n_cluster <- system("qstat | grep job_ | wc | awk '{print $1}'", intern = T)
    message(paste0('\nFuthermore, there are still ', n_cluster, ' jobs running on the cluster.'))
    if(showcluster){
      system("qstat -r | grep jobname | grep -oP \"(?<=job_)(.*)\"")
    }
    message('\n====================================================================================')
    message('====================================================================================')
    message("\n")
    Sys.sleep(sleeptime)
  }

  unlink(lv$file) # Clean up by deleting extra files once done with the loop
}

waitforpostesttofinish <- function(sleeptime = 100,
                                   indicator = indicator,
                                   indicator_group = indicator_group,
                                   run_date = run_date,
                                   strata,
                                   showfiles = TRUE) {

  path <- paste0("<<<<FILEPATH REDACTED>>>>>", indicator_group, '/', indicator, '/output/', run_date, "/temp_post_est/")

  lv <- data.table(strata)
  names(lv) <- "stratum"

  lv[, filename := paste0(path, stratum, "_post_est_list.RData")]

  n_finished <- 0
  n_needed <- length(strata)

  while (n_finished != n_needed) {
    lv[, exists := file.exists(filename)]
    n_finished <- nrow(lv[exists == T])

    message('\n====================================================================================')
    message('\nRunning post-estimation in parallel')
    message(paste0('Run Date: ', run_date))
    message(paste0('\nAt ',Sys.time(),' .... ',n_finished,' out of ', n_needed, ' strata have written output.'))
    if (showfiles) {
      message(paste0('\nCurrently missing strata:'))
      for(i in 1:nrow(lv[exists == F])) {
        message(subset(lv, exists == F)[i, stratum])
      }
    }
    message('\n====================================================================================')
    message('====================================================================================')
    Sys.sleep(sleeptime)
  }
}

waitforresultstable <- function(st = strata,
                                rd = run_date,
                                indic = indicator,
                                ig = indicator_group,
                                baseline_year = 2000,
                                measure = 'prevalence') {

  message("Waiting for results tables to finish...")
  r_left <- length(st)
  u_left <- length(st)
  while(length(r_left) > 0 & length(u_left) > 0) {
    message(paste0("\nCurrent time: ", Sys.time()))
    str_match <- stringr::str_match
    results_dir_r <- paste0("<<<<FILEPATH REDACTED>>>>>", ig, '/',indic,'/output/', rd, '/table_', baseline_year, '/')
    results_dir_u <- paste0("<<<<FILEPATH REDACTED>>>>>", ig, '/',indic,'/output/', rd, '/table_', baseline_year, '_unraked/')

    r_regs_done <- list.files(results_dir_r) %>% str_match(., paste0(measure, "_(.*).csv")) %>% .[,2] %>% .[!is.na(.)]
    u_regs_done <- list.files(results_dir_u) %>% str_match(., paste0(measure, "_(.*).csv")) %>% .[,2] %>% .[!is.na(.)]

    r_left <- st[!(st %in% r_regs_done)]
    u_left <- st[!(st %in% u_regs_done)]

    message(paste0(c("  Raked tables remaining: ", paste(r_left, collapse = ", "))))
    message(paste0(c("  Unraked tables remaining: ", paste(u_left, collapse = ", "))))

    Sys.sleep(60)

  }
}

waitforaggregation <- function(rd = run_date,
                               indic = indicator,
                               ig = indicator_group,
                               ages,
                               regions,
                               holdouts,
                               raked,
                               dir_to_search = NULL) {

  if (is.null(dir_to_search)) {
    dir_to_search <- paste0("<<<<FILEPATH REDACTED>>>>>",indicator_group,"/",indicator,"/output/",run_date,"/")
  }

  lv <- expand.grid(regions, holdouts, ages, raked) %>% as.data.table
  names(lv) <- c("reg", "holdout", "age", "raked")
  lv[, file := paste0(dir_to_search, "fin_agg_", reg, "_", holdout, "_", age, "_", raked)]
  n_left <- nrow(lv)

  message("Waiting for aggregation to finish...")

  while(n_left > 0) {
    message(paste0("\n=====================================",
                   "\nCurrent time: ", Sys.time()))

    lv[, file_exists := file.exists(file)]
    n_left <- nrow(lv[file_exists == F])
    lv_left <- lv[file_exists==F]

    if (n_left > 0) {
      message(paste0("Aggregation jobs remaining: ", n_left))
      for (i in 1:nrow(lv_left)) {
        message(paste0("  ",
                       "Region: ", lv_left[i, reg], " | ",
                       "Holdout: ", lv_left[i, holdout], " | ",
                       "Age: ", lv_left[i, age], " | ",
                       "Raked: ", lv_left[i, raked]))
      }
    } else {
      break
    }

    Sys.sleep(60)

  }
}

#' @title Combine aggregation
#' @description Combine aggregation objects across region
#' @param run_date, indicator, indicator_group for this run
#' @param ages: single value or vector of ages
#' @param regions: vector of regions used
#' @param holdouts: vector of holdouts used, e.g. just 0 or c(1,2,3,4,5,0)
#' @param raked: vector of raked values, e.g. just T, just F, or c(T,F)
#' @param dir_to_search: which directory to search in (defaults to share directory)
#' @param delete_region_files: logical. Should we delete the region-specific intermediate files?
#' @param merge_hierarchy_list: logical. Do you want to merge the sp_hierarchy_list onto your admin tables?
#' @param check_for_dupes PARAM_DESCRIPTION, Default: F
#' @return rdata files for each combo of age/holdout/raked
#'   each with admin_0, admin_1, admin_2 data table objects & the sp_hierarchy_list object
#'   that maps them to names of admin units
combine_aggregation <- function(rd = run_date,
                                indic = indicator,
                                ig = indicator_group,
                                ages,
                                regions,
                                holdouts,
                                raked,
                                dir_to_search = NULL,
                                delete_region_files = T,
                                merge_hierarchy_list = F,
                                check_for_dupes = F) {

  # Combine aggregation objects across region

  # Args:
  #   run_date, indicator, indicator_group for this run
  #   ages: single value or vector of ages
  #   regions: vector of regions used
  #   holdouts: vector of holdouts used, e.g. just 0 or c(1,2,3,4,5,0)
  #   raked: vector of raked values, e.g. just T, just F, or c(T,F)
  #   dir_to_search: which directory to search in (defaults to share directory)
  #   delete_region_files: logical. Should we delete the region-specific intermediate files?
  #   merge_hierarchy_list: logical. Do you want to merge the sp_hierarchy_list onto your admin tables?

  # Outputs:
  #   rdata files for each combo of age/holdout/raked
  #   each with admin_0, admin_1, admin_2 data table objects & the sp_hierarchy_list object
  #   that maps them to names of admin units

  if (is.null(dir_to_search)) {
    dir_to_search <- paste0("<<<<FILEPATH REDACTED>>>>>",ig,"/",indic,"/output/",run_date,"/")
  }

  message("Combining aggregation results...")

  for(rake in raked) {
    for (holdout in holdouts) {
      for (age in ages) {
        message(paste0("\nWorking on age: ", age, " | holdout: ", holdout, " | raked: ", rake))

        # Set up lists
        ad0 <- list()
        ad1 <- list()
        ad2 <- list()
        sp_h <- list()


        for (reg in regions) {
          message(paste0("  Region: ", reg))

          load(paste0(dir_to_search, indic, "_", ifelse(rake == 'raked', "raked", ifelse(rake == 'raked_c', "raked_c", "unraked")),
                      "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))

          if(merge_hierarchy_list == T) {
            # Prepare hierarchy list for adm0
            ad0_list <- subset(sp_hierarchy_list, select = c("ADM0_CODE", "ADM0_NAME", "region")) %>% unique

            # Prepare hierarchy list for adm1
            ad1_list <- subset(sp_hierarchy_list,
                               select = c("ADM0_CODE", "ADM1_CODE", "ADM0_NAME", "ADM1_NAME", "region")) %>%
              unique

            # Merge
            admin_0 <- merge(ad0_list, admin_0, by = "ADM0_CODE", all.y = T)
            admin_1 <- merge(ad1_list, admin_1, by = "ADM1_CODE", all.y = T)
            admin_2 <- merge(sp_hierarchy_list, admin_2, by = "ADM2_CODE", all.y = T)
            rm(ad0_list, ad1_list)
          }
          if(check_for_dupes){
            adms <- get_adm0_codes(reg)
            sp_hier <- get_sp_hierarchy()
            include_ad0 <- sp_hier$ADM0[ADM0_CODE %in% adms, ADM0_CODE]
            include_ad1 <- sp_hier$ADM1[ADM0_CODE %in% adms, ADM1_CODE]
            include_ad2 <- sp_hier$ADM2[ADM0_CODE %in% adms, ADM2_CODE]

            ad0[[reg]] <- admin_0[ADM0_CODE %in% include_ad0]
            ad1[[reg]] <- admin_1[ADM1_CODE %in% include_ad1]
            ad2[[reg]] <- admin_2[ADM2_CODE %in% include_ad2]
            sp_h[[reg]] <- sp_hierarchy_list
          } else{
            ad0[[reg]] <- admin_0
            ad1[[reg]] <- admin_1
            ad2[[reg]] <- admin_2
            sp_h[[reg]] <- sp_hierarchy_list
          }

          rm(admin_0, admin_1, admin_2, sp_hierarchy_list)
        }

        # Get to long format & save
        message("  Combining...")
        admin_0 <- rbindlist(ad0)
        admin_1 <- rbindlist(ad1)
        admin_2 <- rbindlist(ad2)
        sp_hierarchy_list <- rbindlist(sp_h)

        message("  Saving combined file...")
        save(admin_0, admin_1, admin_2, sp_hierarchy_list,
             file = paste0(dir_to_search, indic, "_",
                           ifelse(rake == 'raked', "raked", ifelse(rake == 'raked_c', "raked_c", "unraked")),
                           "_admin_draws_eb_bin", age, "_",
                           holdout, ".RData"))
      }
    }
  }

  if (delete_region_files == T) {
    # Make sure all full files are written
    combos <- expand.grid(ifelse(rake == 'raked', "raked", ifelse(rake == 'raked_c', "raked_c", "unraked")), ages, holdouts)
    files_to_check <- sapply(1:nrow(combos), function(i) {
      paste0(dir_to_search, indic, "_", combos[i, 1],
             "_admin_draws_eb_bin", combos[i, 2], "_", combos[i,3], ".RData")
    })

    if (all(file.exists(files_to_check))) {
      message("All anticipated combined files were created successfully.  Deleting intermediate files...")
      combos <- expand.grid(ifelse(rake == 'raked', "raked", ifelse(rake == 'raked_c', "raked_c", "unraked")), ages, regions, holdouts)
      files_to_delete <- sapply(1:nrow(combos), function(i) {
        paste0(dir_to_search, indic, "_", combos[i, 1],
               "_admin_draws_eb_bin", combos[i, 2], "_", combos[i,3], "_", combos[i, 4], ".RData")
      })
      unlink(files_to_delete)
    } else {
      warning("Did not delete intermediate files - not all output files created successfully!")
    }
  }

  # Finally, delete the "fin" files
  fin_files_to_delete <- list.files(dir_to_search, pattern = "fin_agg_", full.names=T)
  unlink(fin_files_to_delete)
}


insertRaster <- function (raster, new_vals, idx = NULL) {


  # calculate cell index if not provided
  if (is.null(idx)) idx <- cellIdx(raster)

  # check the index makes superficial sense
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))

  # create results raster
  n <- ncol(new_vals)
  raster_new <- raster::brick(replicate(n,
                                        raster[[1]],
                                        simplify = FALSE))
  names(raster_new) <- colnames(new_vals)

  # update the values
  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }

  return (raster_new)

}

condSim <- function (vals, weights = NULL, group = NULL, fun = NULL, ...) {
  # given a matrix of pixel-level prevalence samples `prev`
  # where each rows are pixels and columns are draws, a vector
  # of corresponding pixel populations `pop`, and an optional pixel
  # grouping factor `group`, return draws for the total deaths in each
  # group, or overall if groups are not specified

  # get dimensions of vals
  ncell <- nrow(vals)
  ndraw <- ncol(vals)

  # capture function as a string
  fun_string <- deparse(substitute(fun))

  # check fun accepts a

  # check dimensions of weights and group, set to 1 if not specified
  if (is.null(weights)) {
    weights <- rep(1, ncell)
  } else {
    if (length(weights) != ncell) {
      stop (sprintf('number of elements in weights (%i) not equal to number of cells in vals (%i)',
                    length(weights),
                    ncell))
    }
  }

  if (is.null(group)) {
    group <- rep(1, length(weights))
  } else {
    if (length(group) != ncell) {
      stop (sprintf('number of elements in group (%i) not equal to number of cells in vals (%i)',
                    length(group),
                    ncell))
    }
  }

  # otherwise, get the levels in group and create a matrix of results
  levels <- unique(na.omit(group))
  nlevel <- length(levels)

  ans <- matrix(NA,
                ncol = ndraw,
                nrow = nlevel)
  rownames(ans) <- levels

  # loop through levels in group, getting the results
  for (lvl in 1:nlevel) {

    # get an index o pixels in the level
    idx <- which(group == levels[lvl])

    # by default, calculate a weighted sum
    if (is.null(fun)) {

      # get draws and add to results
      # exception for if area has 1 cell, transpose matrix so it conforms (RB)
      if(all(dim(t(vals[idx, ]))==c(1,ndraw))){
        ans[lvl, ] <- weights[idx] %*% t(vals[idx, ])
      } else {
        ans[lvl, ] <- weights[idx] %*% vals[idx, ]
      }

    } else {

      # otherwise, apply function to each column
      ans[lvl, ] <- apply(vals[idx, ], 2, fun, weights = weights[idx], ...)

    }

  }

  # if only one level, make this a vector
  if (nlevel == 1) ans <- as.vector(ans)

  # return result
  return (ans)

}

defaultOptions <- function (resolution = 5,
                            location = loc(),
                            cores = switch(loc(), oxford = 60, seattle = 60),
                            spacetime = TRUE,
                            start = Sys.time()) {

  # built in check for if they selected more than one country

  # if these options are not already set, set them at these default values

  # list all visible objects in this environment
  object_names <- ls()

  # get them in a named list
  objects <- lapply(object_names,
                    get,
                    envir = environment())
  names(objects) <- object_names

  # empty vector of the options to keep
  keep <- c()

  # loop through identifying those that don't already exist
  for (i in 1:length(objects)) {

    # if doesn't exist yet...
    if (is.null(options(names(objects)[i])[[1]])) {

      # add it to the list
      keep <- c(keep, i)

    }
  }

  # keep the new options
  objects <- objects[keep]

  if (length(objects) > 0) {
    # notify the user
    message('\nThe following options were not defined and have been set to their defaults:\n')
    for(i in 1:length(objects)) {
      message(sprintf('    %s = %s',
                      names(objects)[i],
                      paste(objects[[i]], collapse = ', ')))
    }
    message('\n')
  }

  # add to options
  options(objects)
}


# for cleaning clutter
# https://stackoverflow.com/a/11625075
rmlike <- function(...) {
  names <- sapply(
    match.call(expand.dots = FALSE)$..., as.character)
  names = paste(names,collapse="|")
  Vars <- ls(1)
  r <- Vars[grep(paste("^(",names,").*",sep=""),Vars)]
  rm(list=r,pos=1)
}



reportTime <- function () {

  # get start time stored in options
  start <- options()$start

  # if it exists
  if (!is.null(start)) {

    # get elapsed time
    now <- Sys.time()
    diff <- difftime(now, start)

    # format nicely and report
    diff_string <- capture.output(print(diff, digits = 3))
    diff_string <- gsub('Time difference of', 'Time elapsed:', diff_string)
    message (diff_string)
  }
}


subtractYears <- function (Date, years) {
  # given a vector of dates in Date format,
  # and a number of years (positive numeric, can be non-integer)
  # subtract the required number of years and return

  stopifnot(years >= 0)
  stopifnot(inherits(Date, 'Date'))

  n_Date <- length(Date)
  n_years <- length(years)

  if (n_years != n_Date) {
    if (n_years == 1) {
      years <- rep(years, n_Date)
    } else {
      stop('Date and years have different lengths')
    }
  }

  Date <- as.POSIXlt(Date)

  month <- round(12 * years %% 1)
  year <- floor(years)

  Date$year <- Date$year - year
  Date$mon <- Date$mon - month

  Date <- as.Date(Date)

  return (Date)

}

toYear <- function(Date) format(Date, '%Y')

targetYear <- function (period, period_end, width = 60) {
  # get the target year, based on the period, period size (in months) and end
  # date of the final period

  # convert date to cmc
  period_end <- Date2cmc(period_end)

  # get number of months preceeding
  months <- width / 2 + width * (period - 1)

  # subtract
  cmc <- period_end - months + 1

  # format as a year
  ans <- toYear(cmc2Date(cmc))

  return (ans)

}

# convert strings or numerics to a sequential numeric index
# pass any number of vectors of the same length, get a numeric ID for the
# unique ones
idx <- function (...) {
  list <- list(...)

  # get their lengths
  lengths <- lapply(list, length)
  if (!isTRUE(do.call(all.equal, lengths))) {
    stop('vectors do not appear to have the same length')
  }
  # combine them
  x <- do.call(paste, list)
  #get numeric index
  match(x, unique(x))
}

getPaths <- function(pathfile = '~/Z/ABRAID/prevalence modelling/under five mortality/paths_for_nick.csv') {
  # get file paths for the key datasets
  # path points to a csv file containing named paths, the function returns a dataframe
  # of named filepaths
  paths <- read.csv(pathfile,
                    row.names = 1)
  data.frame(t(paths),
             stringsAsFactors = FALSE)
}


# given the admin level, a GAUL code and the admin1 and 2 shapefiles,
# return an SPDF with the relevant area
getPoly <- function (level, code, admin1, admin2) {

  # get admin level
  if (level == 1) {
    admin <- admin1
  } else {
    admin <- admin2
  }

  # find the reight area
  idx <- match(code, admin$GAUL_CODE)

  # if it's valid
  if (length(idx) == 1) {
    ans <- admin[idx, ]
  } else {
    # handle errors
    warning (paste0("something's wrong on row ", i))
    ans <- NULL
  }

  # return result
  return (ans)
}

# functions to determine proportion of women in the given age group
matAgeRate <- function (age_group,
                        groups = c('15-19', '20-24', '25-29', '30-34')) {
  # given a vector `age_group` reporting the age group to which mothers belong,
  # return the proportion falling in each of the age groups in `groups`.

  # check it's a character vector
  stopifnot(class(age_group) == 'character')

  # count the number in all age groups
  counts <- table(age_group)

  # add on any groups not represented a 0s
  missing_groups <- which(!(groups %in% names(counts)))

  if (length(missing_groups) > 0) {
    dummy_counts <- rep(0, length(missing_groups))
    names(dummy_counts) <- groups[missing_groups]
    counts <- c(counts, dummy_counts)
  }

  # get proportions
  props <- counts / sum(counts)

  # find the ones we want
  idx_keep <- match(groups, names(props))

  # and return these
  return (props[idx_keep])

}


# functions to determine proportion of women in the given age group
matAgeParRate <- function (ceb,
                           age_group,
                           groups = c('15-19', '20-24', '25-29', '30-34')) {
  # given vectors `ceb` and `age_group` reporting the number of children ever
  # born to mothers and the age groups to which they belong,
  # return the proportion of births falling in each of the age groups in
  # `groups`.

  # check it's a character vector
  stopifnot(class(ceb) %in% c('numeric', 'integer'))
  stopifnot(class(age_group) == 'character')

  # count the number of births in all age groups
  counts <- tapply(ceb, age_group, sum)

  # add on any groups not represented a 0s
  missing_groups <- which(!(groups %in% names(counts)))

  if (length(missing_groups) > 0) {
    dummy_counts <- rep(0, length(missing_groups))
    names(dummy_counts) <- groups[missing_groups]
    counts <- c(counts, dummy_counts)
  }

  # get proportions
  props <- counts / sum(counts)

  # find the ones we want
  idx_keep <- match(groups, names(props))

  # and return these
  return (props[idx_keep])

}

# function to tabulate maternal age rate
tabMatAgeRate <- function (age_group,
                           cluster_id,
                           groups = c('15-19', '20-24', '25-29', '30-34')){
  # given a vector `age_group` reporting the age group to which mothers belong,
  # and a vector `cluster_id` giving the cluster to which each mother belongs
  # return a matrix - with number of rows equal to the number of unique elements
  # in `cluster_id` and number of columns equal to the length of `groups` -
  # giving the proportion falling in each of the age groups in `groups`.

  # get the grouped data as a list
  ans <- tapply(age_group,
                cluster_id,
                matAgeRate,
                groups)

  # combine into a matrix
  ans <- do.call(rbind, ans)

  # and return
  return (ans)

}


# function to tabulate maternal age rate
tabMatAgeParRate <- function (ceb,
                              age_group,
                              cluster_id,
                              groups = c('15-19', '20-24', '25-29', '30-34')){
  # given vectors `ceb` and `age_group` reporting the number of children ever
  # born to mothers and the age groups to which they belong,
  # return a matrix - with number of rows equal to the number of unique elements
  # in `cluster_id` and number of columns equal to the length of `groups` -
  # giving the proportion of births falling in each of the age groups in
  # `groups`.

  # get unique cluster ids
  cluster_ids <- unique(cluster_id)

  # create dummy results matrix
  ans <- matrix(NA,
                nrow = length(cluster_ids),
                ncol = length(groups))

  rownames(ans) <- cluster_ids
  colnames(ans) <- groups

  # loop through the cluster ids calculating the rates
  for (i in 1:length(cluster_ids)) {

    # get index for the cluster
    idx_cluster <- which(cluster_id == cluster_ids[i])

    ans[i, ] <- matAgeParRate(ceb = ceb[idx_cluster],
                              age_group = age_group[idx_cluster],
                              groups = groups)

  }

  # and return
  return (ans)

}

# get predictions from predictive INLA glm models
predGLM <- function (result,
                     intercept = '(Intercept)',
                     fixed_continuous = NULL,
                     fixed_group = NULL,
                     random_group = 'cluster_id') {

  # starting means and variances
  pred_mean <- pred_var <- 0

  # add intercept terms
  if (!is.null(intercept)) {
    pred_mean <- pred_mean + result$summary.fixed[intercept, 'mean']
    pred_var <- pred_var + result$summary.fixed[intercept, 'sd'] ^ 2
  }

  # add continuous fixed effects terms
  if (!is.null(fixed_continuous)) {
    for (name in fixed_continuous) {

      # get coefficients
      coef_mean <- result$summary.fixed[name, 'mean']
      coef_var <- result$summary.fixed[name, 'sd'] ^ 2

      # get covariate
      cov <- result$model.matrix[, name]

      pred_mean <- pred_mean + cov * coef_mean
      pred_var <- pred_var + (cov ^ 2) * coef_var

    }
  }

  # add discrete fixed effects terms
  if (!is.null(fixed_group)) {
    for (name in fixed_group) {

      # find group members
      idx <- grep(sprintf('^%s*', name), rownames(result$summary.fixed))

      # get clean names
      names <- rownames(result$summary.fixed)[idx]
      names_clean <- gsub(name, '', names)

      # get coefficients
      coef_mean <- result$summary.fixed[idx, 'mean']
      coef_var <- result$summary.fixed[idx, 'sd'] ^ 2

      # get covariate
      cov <- result$model.matrix[, names]

      pred_mean <- pred_mean + as.vector(cov %*% coef_mean)
      pred_var <- pred_var + as.vector((cov ^ 2) %*% coef_var)

    }
  }

  # add discrete random effects terms
  if (!is.null(random_group)) {
    for (name in random_group) {

      # get coefficients
      coef_level <- result$summary.random[[name]][, 'ID']
      coef_mean <- result$summary.random[[name]][, 'mean']
      coef_var <- result$summary.random[[name]][, 'sd'] ^ 2

      # get covariate
      cov <- result$.args$data[, name]

      # match up the levels
      length(cov)
      length(coef_level)
      idx <- match(cov, coef_level)

      pred_mean <- pred_mean + coef_mean[idx]
      pred_var <- pred_var + coef_var[idx]

    }
  }

  # return result
  ans <- data.frame(mean = pred_mean,
                    sd = sqrt(pred_var))

  return (ans)

}

getCols <- function (df, period = 1) {
  # subset results matrix
  df[, grep(sprintf('^p%s_', period),
            colnames(df))]
}

rnormMatrix <- function (n, mean, sd) {
  # sample random normals with matrix parameters
  # returna an array as a result

  # coerce to matrix
  mean <- as.matrix(mean)
  sd <- as.matrix(sd)

  # get & check dimensions
  ncol <- ncol(mean)
  nrow <- nrow(mean)
  stopifnot(ncol(sd) == ncol)
  stopifnot(nrow(sd) == nrow)

  # convert to vector
  mean <- as.vector(mean)
  sd <- as.vector(sd)

  # sample
  draws <- rnorm(n * length(mean), mean, sd)

  # reshape
  ans <- array(draws, dim = c(nrow, ncol, n))

  return (ans)
}

accumulate <- function (mean, sd, months = c(1, 11, 24, 24), nsample = 1000) {
  # given matrices of means and standard deviations of
  # logit component death probabilities, each column giving
  # consecutive and adjacent age bins and rows giving the records,
  # calculate the logit of the cumulative mortality probability across
  # the age bins by Monte Carlo sampling

  # generate random draws from these logits
  draws <- rnormMatrix(nsample, mean, sd)

  # convert to draws of survival probabilities
  draws_p <- 1 - plogis(draws)

  # raise to power of number of months per bin
  for (i in 1:dim(draws_p)[2]) {
    draws_p[, i, ] <- draws_p[, i, ] ^ months[i]
  }

  # loop through age bins accumulating them
  for (i in 2:dim(draws_p)[2]) {
    draws_p[, i, ] <- draws_p[, i, ] * draws_p[, i - 1, ]
  }

  # convert back to logit mortality probabilities
  draws_l <- qlogis(1 - draws_p)

  # calculate the means and standard deviations of these logits
  l_mean <- apply(draws_l, c(1, 2), mean)
  l_sd <- apply(draws_l, c(1, 2), sd)

  # return as a list
  return (list(mean = l_mean,
               sd = l_sd))

}


# functions to centre and scale matrices, columnwise

getCentreScale <- function (x, exclude = NULL, na.rm = TRUE) {
  # get dataframe of centreing and scaling values to convert x
  # to the standard normal. exclude is an optional character vector
  # giving column names to exclude from scaling

  # get means and SDs for all columns
  df <- data.frame(name = colnames(x),
                   mean = colMeans(x, na.rm = na.rm),
                   sd = apply(x, 2, sd, na.rm = na.rm))
  rownames(df) <- NULL

  # replace any zero standard deviations with 1
  # to avoid divide-by-zero errors
  df$sd[df$sd == 0] <- 1

  # if any named covariates are to be excluded, set mean to 0 and sd to 1
  if (!is.null(exclude)) {
    idx <- match(exclude, df$name)
    df$mean[idx] <- 0
    df$sd[idx] <- 1
  }

  return (df)
}

centreScale <- function (x, df, inverse = FALSE) {
  # apply pre-calculated centreing/scaling to matrix x,
  # with fixed dataframe of means/sds df
  # or uncentre/unscale if inverse = TRUE

  # get the centreing/scaling dataframe if not available
  if (is.null(df))
    df <- getCentreScale(x)

  # get index to match up values with column names
  names <- colnames(x)
  idx <- match(names, df$name)

  if (any(is.na(idx))) {
    stop ('could not match up column names with the values in df')
  }

  df <- df[idx, ]

  # apply transformations
  if (!inverse) {
    # move to standard normal

    # centre
    x <- sweep(x, MARGIN = 2, STATS = df$mean, FUN = '-')
    # scale
    x <- sweep(x, MARGIN = 2, STATS = df$sd, FUN = '/')

  } else {
    # inverse case, move from standard normal to original

    # unscale
    x <- sweep(x, MARGIN = 2, STATS = df$sd, FUN = '*')
    # uncentre
    x <- sweep(x, MARGIN = 2, STATS = df$mean, FUN = '+')

  }

  return (x)

}

startLog <- function(file = 'full_run.log') {
  # create a logfile and start logging
  con <- file(file)
  sink(con,
       split = TRUE)
  sink(con,
       type = 'message')

  # report session info
  message(sprintf('# starting log at %s\n', Sys.time()))
  message('# session info:\n')
  print(sessionInfo())
  message('\n# run log:\n')
}

stopLog <- function() {
  # report session info
  message('\n# session info:\n')
  print(sessionInfo())
  message(sprintf('\n# stopping log at %s', Sys.time()))
  # stop logging to the logfile
  sink()
  sink(type = 'message')
}

list2paths <- function (x) {

  # get node and branch names
  nodes <- unlist(x, use.names = FALSE)
  branches <- unlist(getBranches(x),
                     use.names = FALSE)

  # append the node name
  paths <- file.path(branches, nodes)

  # remove extra slashes
  paths <- gsub('/+', '/', paths)

  paths
}

getBranches <- function(x, parent = "") {
  # loop through levels of a tree-like list,
  # flattening names of levels into
  # slash-separated character vector

  # get element names
  names <- names(x)

  # if unnamed, repeat the parent names
  if (is.null(names)) {
    names <- rep.int(parent, length(x))
  } else{
    names <- paste0(parent, names)
  }

  # if there are more branches ahead
  if (is.list(x)) {

    # add a slash
    names <- paste0(names, "/")

    # define a function to loop through them
    getTwigs <- function(i) {
      getBranches(x[[i]], names[i])
    }

    # get the names from these
    names <- lapply(1:length(x), getTwigs)

  }

  return (names)

}

unpack <- function (tmp) {
  # unpack a nested list where the outer list is unnamed, but each inner
  # list contains a single named element.
  # Assign this element to the parent environment and delete the list
  # it is called on fromt he parent ennvironment

  # get name of list in calling environment
  tmp_name <- deparse(substitute(tmp))

  # get calling environment
  pf <- parent.frame()

  # unpack into a single list
  tmp <- unlist(tmp, recursive = FALSE)

  # loop through assigning to calling environment
  for (i in 1:length(tmp)) {
    assign(names(tmp)[i], tmp[[i]], envir = pf)
  }

  # remove object from calling environment
  rm(list = tmp_name, envir = pf)

  # return nothing
  return (invisible(NULL))
}

prepLines <- function (rate,
                       year,
                       ylab = '',
                       title = '',
                       xlim = c(1999, 2016),
                       ylim = c(0, max(rate)),
                       line_years = c(2000, 2005, 2010, 2015)) {
  # set up the base plot for the custom line chart

  plot(rate ~ year,
       type = 'n',
       ylab = '',
       xlab = '',
       axes = FALSE,
       xlim = xlim,
       ylim = ylim)

  for (ly in line_years) {
    lines(x = rep(ly, 2),
          y = ylim,
          lwd = 3,
          lty = 3,
          col = grey(0.8))
  }

  axis(1,
       lty = 0,
       col.axis= grey(0.4),
       line = -1)

  axis(2,
       las = 2,
       cex.axis = 0.8,
       col = grey(0.4),
       col.axis= grey(0.4))

  title(main = title,
        col.main = grey(0.35),
        cex.main = 1.2,
        line = 0.5)

  title(ylab = ylab,
        col.lab = grey(0.4),
        cex.lab = 1.2)
}

addLines <- function (rate,
                      year,
                      country,
                      countries = NULL,
                      col = grey(0.7),
                      size = 1,
                      border = grey(0.4)) {

  # given vectors: rate (y axis), year (x axis)
  # and country (grouping factor), make a nice line
  # plot with lollipop-like likes with border colour
  # 'border' and fill colour 'col' (repeated if length one)
  # If 'countries' is NULL, all countries are plotted,
  # otherwise only those named in this character vector

  # check inputs
  stopifnot(all.equal(length(rate),
                      length(year),
                      length(country)))


  # sort countries
  all_countries <- sort(unique(country))
  if (is.null(countries)) {
    countries <- all_countries
  } else {
    stopifnot(all(countries %in% all_countries))
  }
  n_ctry <- length(countries)

  # expand col and bg if needed
  if (length(col) == 1) {
    col <- rep(col, n_ctry)
  } else {
    stopifnot(length(col) == n_ctry)
  }

  # loop through countries
  for (i in 1:n_ctry) {

    ctry <- countries[i]

    idx_ctry <- which(country == ctry)

    # dark grey outline
    lines(rate[idx_ctry] ~ year[idx_ctry],
          col = border,
          lwd = 7.5 * size)
    points(rate[idx_ctry] ~ year[idx_ctry],
           col = border,
           cex = 1 * size,
           pch = 16)

    # coloured foreground
    lines(rate[idx_ctry] ~ year[idx_ctry],
          col = col[i],
          lwd = 6 * size)

    points(rate[idx_ctry] ~ year[idx_ctry],
           col = col[i],
           pch = 16,
           cex = 0.85 * size)

  }
}

addLabels <- function (rate,
                       year,
                       country,
                       countries = NULL,
                       col = grey(0.7),
                       gap = diff(range(rate)) / 60,
                       cex = 0.7,
                       adj = 0,
                       xpd = NA,
                       ...) {

  # add country names on RHS
  # arguments as before, with gap to define spacing between labels.
  # dots are passed to text
  require (plotrix)

  # check inputs
  stopifnot(all.equal(length(rate),
                      length(year),
                      length(country)))

  # sort countries
  all_countries <- sort(unique(country))
  if (is.null(countries)) {
    countries <- all_countries
  } else {
    stopifnot(all(countries %in% all_countries))
  }
  n_ctry <- length(countries)

  # expand col and bg if needed
  if (length(col) == 1) {
    col <- rep(col, n_ctry)
  } else {
    stopifnot(length(col) == n_ctry)
  }

  # keep only those for latest year
  max_year <- max(as.numeric(year))
  year_idx <- which(year == max_year)
  rate <- rate[year_idx]
  year <- year[year_idx]
  country <- country[year_idx]

  # keep only those for countries requires
  ctry_idx <- which(country %in% countries)
  rate <- rate[ctry_idx]
  year <- year[ctry_idx]
  country <- country[ctry_idx]

  # order by countries
  ctry_o <- match(countries, country)
  rate <- rate[ctry_o]
  year <- year[ctry_o]
  country <- country[ctry_o]

  # get a good gap
  y_pos <- plotrix::spreadout(rate, gap)
  text(x = max_year + 1,
       y = y_pos,
       labels = country,
       col = col,
       cex = cex,
       adj = adj,
       xpd = xpd,
       ...)

}


splitGeoNames <- function (geo) {
  # split country/years (in format 'country_year') out of rownames
  # of a geostatistical conditional simulation object and add as columns
  splits <- strsplit(rownames(geo), split = '_')
  ctry <- sapply(splits, '[', 1)
  year <- sapply(splits, '[', 2)
  geo <- data.frame(iso3 = ctry,
                    year = year,
                    geo,
                    stringsAsFactors = FALSE)
  rownames(geo) <- NULL
  return (geo)
}

getEst <- function (df, iso3, year) {
  # for each dataframe of estimates, line up the country
  # and year, add object name as prefix to the column nmes and
  # return as a dataframe

  # get object name
  prefix <- deparse(substitute(df))

  # get index
  iso3_year_target <- paste(iso3, year, sep = '_')
  iso3_year_df <- paste(df$iso3, df$year, sep = '_')
  idx <- match (iso3_year_target, iso3_year_df)

  # subset, add prefix to column names and return
  df <- df[idx, -(1:2)]
  colnames(df) <- paste(prefix, colnames(df), sep = '_')
  return (df)
}

elogit <- function (y, n) log ( (y + 0.5) / (n - y + 0.5) )

## Combine input data and covariates layers for models run by region
## Necessary for saving to Shiny tool
## Save "covs", "tv_*", and "df" to new combined snapshot in model_image_history
combine_region_image_history <- function(indicator, indicator_group, run_date, fixed_effects) {


  load(paste0("<<<<FILEPATH REDACTED>>>>>"))
  new_cov_list <- list()
  for(cov in names(cov_list)[!grepl("gaul_code", names(cov_list))]) {
    pull_raster_covs <- function(region) {
      load(paste0("<<<<FILEPATH REDACTED>>>>>"))
      cov_raster <- cov_list[[cov]]
      return(cov_raster)
    }
    region_layers <- lapply(Regions, pull_raster_covs)
    combined_layers <- do.call(raster::merge, region_layers)
    names(combined_layers) <- gsub("layer", cov, names(combined_layers))
    if(length(names(combined_layers))==1) new_cov_list[[cov]] <- combined_layers
    if(length(names(combined_layers))!=1) assign(paste0('tv_', cov), combined_layers)
  }
  #cov_list <- do.call(raster::brick, new_cov_list)
  cov_list <- brick(new_cov_list)

  # Combine input data
  pull_df <- function(region) {
    load(paste0("<<<<FILEPATH REDACTED>>>>>"))
    return(df)
  }
  df <- lapply(Regions, pull_df)
  df <- do.call(rbind.fill, df)

  covs <- cov_list

  save(list = c('df','covs',grep('^tv_*', ls(), value = TRUE)), file = paste0("<<<<FILEPATH REDACTED>>>>>"))

}


# These 3 functions are used as a summstats argument via a configuration option. DO NOT DELETE
cirange = function(x){
  z=quantile(x,probs=c(.025,.975),na.rm=T)
  return(z[2]-z[1])
}
lower = function(x) quantile(x,probs=.025,na.rm=T)
upper = function(x) quantile(x,probs=.975,na.rm=T)


cleanup_inla_scratch <- function(run_date) {

  if(keep_inla_files==FALSE) {

    # Clean up INLA intermediate directories unless user has specified to keep them.
    inla_working_dir <- '<<< FILEPATH REDACTED >>>>'
    inla_working_dirs <- list.dirs(inla_working_dir, recursive = FALSE)
    inla_working_dirs <- inla_working_dirs[grepl(run_date, inla_working_dirs)]
    for(inla_dir in inla_working_dirs) {
      unlink(inla_dir, recursive=TRUE)
    }

  }

  if(keep_inla_files==TRUE) {

    message('Keeping INLA intermediate files because keep_inla_files==TRUE in config.')
    message(paste0('Files stored here: <<< FILEPATH REDACTED >>>>'))

  }

}

region_violin <- function(indicator,
                          indicator_group,
                          run_date,
                          output_file) {

  results_dir <- paste0("<<<<FILEPATH REDACTED>>>>>", indicator_group, '/', indicator, '/output/', run_date)
  preds <- raster(paste0(results_dir,'/',indicator,'_mean_raster.tif'))
  extract_region <- function(region) {
    input_data <- fread(paste0(results_dir,'/input_data_bin0_',region,'_0.csv'))
    input_data <- input_data[weight==1 & N >= 20,]
    preds_at_points <- extract(preds, input_data[, c('longitude','latitude'), with = F])
    input_data <- input_data[, pred := preds_at_points]
    input_data <- input_data[!is.na(pred),]
    input_data <- input_data[, data := get(indicator)/N]
    input_data <- input_data[, c('latitude','longitude','pred','data'), with=F]
    input_data <- melt(input_data, id.vars = c('longitude','latitude'), measure.vars = c('data','pred'))
    input_data <- input_data[, region := region]
  }
  input_data <- rbindlist(lapply(c('essa','wssa','cssa','sssa','name'), extract_region))
  pdf(output_file)
  library(ggplot2)
  ggplot(data=input_data) + geom_violin(aes(x = variable, y = value, fill = region)) + facet_wrap(~region)
  dev.off()

}

# Timer Functions ---------------------------------------------------------

## Functions to time within the child stacking regions
##
## General usage:
##    require(tictoc)
##
##    tic("Step 1")
##    **your code here**
##    toc(log = T)
##
##    ticlog <- tic.log(format = F)
##    generate_time_log(ticlog)
##
##  Returns: data table with two columns
##     "step": names of events (e.g. "Step 1")
##     "time": time elapsed (as text: Xh Xm Xs)
##
##  Note: can nest tic/toc pairs

generate_time_log <- function(ticlog) {

  require(magrittr)
  require(data.table)
  # Functions in functions
  strip_time <- function(x) {
    sec <- as.numeric(x$toc - x$tic)
    time <- format_time(sec)
    name <- x$msg

    df <- c(name, time) %>%
      t %>%
      as.data.table

    names(df) <- c("step", "time")

    return(df)
  }

  format_time <- function(run_time) {
    run_time <- round(as.numeric(run_time),0)

    hours <- run_time %/% 3600
    remainder <- run_time %% 3600
    minutes <- remainder %/% 60
    seconds <- remainder %% 60

    run_time <- paste0(hours, "h ", minutes, "m ", seconds, "s")
    return(run_time)
  }

  df_out <- lapply(ticlog, strip_time) %>% rbindlist

  return(df_out)

}


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

  sharedir <- paste0("<<<<FILEPATH REDACTED>>>>>", ig, '/', indic, '/output/', rd, '/')

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

  write.csv(all_mods, paste0(sharedir, indic, "_model_results_table.csv"),
            row.names = F)

  return(mods)

}

## ~~~~~~~~~~~ make table of model results ~~~~~~~~~~~~~~~~~

model_fit_table <- function(lv=loopvars[loopvars[,3]==0,],
                            rd=run_date,
                            nullmodel='',
                            indicator = indicator,
                            indicator_group = indicator_group,
                            holdout = 0,
                            coefs.sum1,
                            use_gp = use_gp,
                            spde_prior = spde_prior,
                            use_stacking_covs = use_stacking_covs,
                            stacker_name_vec,
                            fit_with_tmb){
  # load models
  require(INLA)
  message(sprintf('Pulling together table for %s models',rd))
  tlist=list()
  sharedir <- paste0("<<<<FILEPATH REDACTED>>>>>", indicator_group, '/', indicator, '/output/', rd, '/')
  for(i in 1:nrow(lv)){
    reg <- lv[i,1]
    age <- lv[i,2]
    message(sprintf('%i of %i loopvars. %s %i',i,nrow(lv),lv[i,1],lv[i,2]))

    # Recreate the INLA data stack
    pathaddin <- paste0('_bin', age, '_', reg, '_', holdout)
    f=paste0("<<<<FILEPATH REDACTED>>>>>")

    if(!file.exists(f)){
      message('FAILED TO OPEN')
    } else {
      load(f)
    }

    stacker_name_vec <- intersect(stacker_name_vec, colnames(df))
    if(length(stacker_name_vec) > 0){
      df = df[,paste0(stacker_name_vec) := lapply(stacker_name_vec,
                                                  function(x) get(paste0(x,'_cv_pred')))]
    }

    # Get the spde
    input_data <- build_mbg_data_stack(df = df,
                                       fixed_effects = all_fixed_effects,
                                       mesh_s = mesh_s,
                                       mesh_t = mesh_t,
                                       use_ctry_res = FALSE,
                                       use_nugget = FALSE,
                                       stacker_names = stacker_name_vec,
                                       exclude_cs    = stacker_name_vec,
                                       spde_prior = spde_prior)

    spde <- input_data[[2]]

    # Now get & transform model fit
    message('::::loading in INLA fit\n')
    f=paste0(sharedir, indicator, '_model_eb_bin', age, "_", reg, "_0.RData")

    if(!file.exists(f)){
      message('FAILED TO OPEN')
    } else {
      load(f)
    }

    if(fit_with_tmb == TRUE){
      tlist[[paste0(reg, "_", age)]] <- fitted_param_table_tmb(res_fit)
    }  else{
      ## columns we'll show return
      keep.cols <- c('0.025quant', '0.5quant', '0.975quant')

      ## other hyperparmas
      hyps <- summary(res_fit)$hyperpar[-(1:2), keep.cols] ## first two rows are
      ## theta1/range, theta2/sd

      if(as.logical(use_gp)){
        if(eval(parse(text=spde_prior))$type=="pc") {
          ## extract values from the fit directly
          range <- res_fit$summary.hyperpar[1,keep.cols]
          nom.var <- res_fit$summary.hyperpar[2,keep.cols]^2
        } else {
          ## now we extract what we need from the fit to get transformed spatial params
          res.field <- INLA::inla.spde2.result(res_fit, 'space', spde, do.transf=TRUE)

          ## nominal range at 0.025, 0.5, 0.975 quantiles
          range   <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.range.nominal[[1]])
          nom.var <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.variance.nominal[[1]])
        }
        spat.hyps <- rbind(range, nom.var)
        rownames(spat.hyps) <- c('Nominal Range', 'Nominal Variance')
        colnames(spat.hyps) <- keep.cols
      }

      ## fixed effects from coefs.sum1
      if(as.logical(coefs.sum1) & as.logical(use_stacking_covs)){
        fixed.sum1 <- res_fit$summary.random$covar
        fixed.sum1$ID <- NULL
        rownames(fixed.sum1) <- stacker_name_vec
        fixed.sum1 <- fixed.sum1[, keep.cols]
      }else{
        fixed.sum1 <- NULL
      }

      ## all other coefs (e.g. intercept and raw covs)
      fixed <- summary(res_fit)$fixed
      if(is.null(nrow(fixed))){
        fixed <- matrix(fixed, ncol = length(fixed)) ## in the event of one row, convert numeric back to data.frame
        rownames(fixed) <- rownames( summary(res_fit)$fixed )
        colnames(fixed) <- colnames( summary(res_fit)$fixed )
      }
      fixed <- fixed[, keep.cols]

      ## combine the two types of 'fixed' results
      fixed <- rbind(fixed, fixed.sum1)

      ## combine them all and just keep three quantiles
      all.res <- rbind(fixed,
                       if(use_gp){spat.hyps}else{NULL},
                       hyps)

      ## rename GPRandom rho for time
      all.res <- as.data.table(all.res, keep.rownames = T)
      setnames(all.res, "rn", "parameter")
      if(use_gp) all.res[parameter == "GroupRho for space", parameter := "GPRandom rho for time"]
      all.res

      tlist[[paste0(reg, "_", age)]] <- all.res
    }

  }

  return(tlist)
}

check_config <- function(cr = core_repo) {

  # TODO: update package to use data()
  # data(must_haves)
  must_haves <- read.csv(paste0(cr, '/mbg_central/share_scripts/common_inputs/config_must_haves.csv'), header = F, stringsAsFactors = F)$V1

  message("\nRequired covariates: ")
  for(confs in must_haves){
    if (exists(confs)) {
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'use_gp') {
      message("You are missing a 'use_gp' argument in your config. Defaulting it to TRUE")
      use_gp <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'use_stacking_covs') {
      message("You are missing a 'use_stacking_covs' argument in your config. Defaulting it to TRUE")
      use_stacking_covs <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'use_raw_covs') {
      message("You are missing a 'use_raw_covs' argument in your config. Defaulting it to FALSE")
      use_raw_covs <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'fit_with_tmb') {
      message("You are missing a 'fit_with_tmb' argument in your config. Defaulting it to FALSE")
      fit_with_tmb <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'gbd_fixed_effects_measures') {
      message("You are missing a 'gbd_fixed_effects_measures' argument in your config. Defaulting it to 'covariate' for all elements of gbd_fixed_effects")
      gbd_fixed_effects_measures <<- paste(rep("covariate", length(strsplit(gbd_fixed_effects, split = " \\+ ")[[1]])), collapse = " + ")
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'gbd_fixed_effects_age') {
      message("You are missing a 'gbd_fixed_effects_age' argument in your config. Defaulting to '2 3 4 5'")
      gbd_fixed_effects_age <<- '2 3 4 5'
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'z_list') {
      message("You are missing a 'z_list' argument in your config. Defaulting it to 0")
      z_list <<- 0
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'zcol') {
      message("You are missing a 'zcol' argument in your config. Defaulting it to z_column_default_blank")
      zcol <<- 'z_column_default_blank'
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'summstats') {
      message("You are missing a 'summstats' argument in your config. Defaulting to c('mean','lower','upper','cirange')")
      summstats <<- c('mean','lower','upper','cirange')
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == 'scale_gaussian_variance_N') {
      message("You are missing a 'scale_gaussian_variance_N' argument in your config. Defaulting to TRUE")
      scale_gaussian_variance_N <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "nugget_prior") {
      message("You are missing a 'nugget_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      nugget_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "ctry_re_prior") {
      message("You are missing a 'ctry_re_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      ctry_re_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "nid_re_prior") {
      message("You are missing a 'nid_re_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      nid_re_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "use_nid_res") {
      message("You are missing a 'use_nid_res' argument in your config. Defaulting to FALSE")
      use_nid_res <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "rho_prior") {
      message("You are missing a 'rho_prior' argument in your config. Defaulting to 'list(prior = 'normal', param = c(0, 0.1502314))'")
      rho_prior <<- "list(prior = 'normal', param = c(0, 1/(2.58^2)))"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "gp_sum_to_zero") {
      message("You are missing a 'gp_sum_to_zero' argument in your config. Defaulting to FALSE")
      gp_sum_to_zero <<- FLASE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "use_s2_mesh") {
      message("You are missing a 'use_s2_mesh' argument in your config. Defaulting to FALSE")
      use_s2_mesh <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "s2_mesh_params") {
      message("You are missing a 's2_mesh_params' argument in your config. Defaulting to c(50, 500, 1000)")
      s2_mesh_params <<- "c(25, 500, 1000)"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "sparse_ordering") {
      message("You are missing a 'sparse_ordering' argument in your config. Defaulting to TRUE")
      sparse_ordering <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "modeling_shapefile_version") {
      message("You are missing a 'modeling_shapefile_version' argument in your config. Defaulting to 'current'")
      modeling_shapefile_version <<- "current"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "raking_shapefile_version") {
      message("You are missing a 'raking_shapefile_version' argument in your config. Defaulting to 'current'")
      raking_shapefile_version <<- "current"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "subnational_raking") {
      message("You are missing a 'subnational_raking' argument in your config. Defaulting to TRUE")
      subnational_raking <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "check_cov_pixelcount") {
      message("You are missing a 'check_cov_pixelcount' argument in your config. Defaulting to FALSE")
      check_cov_pixelcount <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "gbd_fixed_effects_constraints") {
      message("You are missing a 'gbd_fixed_effects_constraints' argument in your config. Defaulting to FALSE")
      gbd_fixed_effects_constraints <<- "c(0)"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "fixed_effects_constraints") {
      message("You are missing a 'fixed_effects_constraints' argument in your config. Defaulting to FALSE")
      fixed_effects_constraints <<- "c(0)"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "memory") {
      message("You are missing a 'memory' argument in your config. Defaulting to 10G")
      memory <<- 10

    } else if (confs == "singularity_version") {
      message("You are missing a 'singularity_version' argument in your config. Defaulting to 'default'")
      singularity_version <<- "default"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "queue") {
      message("You are missing a 'queue' argument in your config. Defaulting to 'long.q', unless you have use_geos_nodes to TRUE, which will override this to geospatial.q")
      queue <<- "long.q"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "run_time") {
      message("You are missing a 'run_time' argument in your config. Defaulting to 16 days ('16:00:00:00')")
      run_time <<- "16:00:00:00"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "countries_not_to_rake") {
      message("You are missing a 'countries_not_to_rake' argument in your config. Defaulting to ESH+GUF")
      countries_not_to_rake <<- "ESH+GUF"
      message(paste0("  ", confs, ": ", get(confs)))

    } else if (confs == "countries_not_to_subnat_rake") {
      message("You are missing a 'countries_not_to_subnat_rake' argument in your config. Defaulting to PHL+NGA+PAK+ETH+KEN")
      countries_not_to_subnat_rake <<- "PHL+NGA+PAK+ETH+KEN"
      message(paste0("  ", confs, ": ", get(confs)))

    } else if (confs == "rake_countries") {
      message("You are missing a 'rake_countries' argument in your config. Defaulting to TRUE")
      rake_countries <<- TRUE
      message(paste0("  ", confs, ": ", get(confs)))

    } else if (confs == "use_space_only_gp") {
      message("You are missing a 'use_space_only_gp' argument in your config. Defaulting to FALSE")
      use_space_only_gp <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "st_gp_int_zero") {
      message("You are missing a 'st_gp_int_zero' argument in your config. Defaulting to FALSE")
      st_gp_int_zero <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "s_gp_int_zero") {
      message("You are missing a 's_gp_int_zero' argument in your config. Defaulting to FALSE")
      s_gp_int_zero <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "use_time_only_gmrf") {
      message("You are missing a 'use_time_only_gmrf' argument in your config. Defaulting to FALSE")
      use_time_only_gmrf <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "time_only_gmrf_type") {
      message("You are missing a 'time_only_gmrf_type' argument in your config. Defaulting to FALSE")
      time_only_gmrf_type <<- "rw2"
      message(paste0('  ', confs, ': ', get(confs)))

    } else if (confs == "spde_prior") {
      message("You are missing a 'spde_prior' argument in your config. Defaulting to 'list(type='pc')'")
      spde_prior <<- "list(type='pc')"
      message(paste0('  ', confs, ': ', get(confs)))

    } else {
      stop(paste0(confs, " is missing, add it to your config"))
    }
  }


  ## Test for subnational random effect
  if(exists("use_subnat_res", envir = .GlobalEnv)) {
    stopifnot(exists("subnat_country_to_get", envir = .GlobalEnv))
    # stopifnot(length(eval(parse(text = subnat_country_to_get))) == 1)
  } else {
    use_subnat_res <<- FALSE
    subnat_country_to_get <<- FALSE
  }


  message("\nAdditional config arguments: ")
  extras <- config$V1[!(config$V1 %in% must_haves)]
  for (extra in extras) message(paste0('  ', extra, ': ', get(extra)))

  ## print out shapefile info
  m.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = modeling_shapefile_version))
  r.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = raking_shapefile_version))
  message("\n\n\nSHAPEFILE VERSION INFORMATION: ")
  message(sprintf("\n--MODELING SHAPEFILE VERSION: %s -- which contains %s codes", m.sf.info$shpfile_date, toupper(m.sf.info$shpfile_type)))
  message(sprintf("\n--RAKING SHAPEFILE VERSION:   %s -- which contains %s codes\n", r.sf.info$shpfile_date, toupper(r.sf.info$shpfile_type)))
}


#' @title Easy eval-parse
#' @description Allows for easily eval-parsing through a config dataset
#' @param data The data.table
#' @param column The column with the string call
#' @return Evaluated call
#' @export
#' @rdname ez_evparse
ez_evparse <- function(data, column) {
  return(eval(parse(text = data[, column, with = FALSE])))
}






#' @title Set up config
#' @description Setting up configuration variables for an MBG run
#' @param repo Location where you've cloned the MBG repository for your indicator.
#' @param core_repo Location where you've cloned the lbd_core repository. Not necessary in the package version.
#' @param indicator_group Category of indicator, i.e. "education"
#' @param indicator Specific outcome to be modeled within indicator category, i.e. "edu_0"
#' @param config_name Name of configuration file in the indicator folder, Default: NULL
#' @param config_file Full path to configuration file that overrides \code{config_name}, Default: NULL
#' @param covs_name Name of covariates configuration file, Default: NULL
#' @param covs_file Full path to covariates configuration file that overrides \code{covs_name}, Default: NULL
#' @param post_est_only Set up only for post estimation? Default: FALSE
#' @param run_date Run date, Default: ''
#' @param push_to_global_env Should the config parameters be pushed to the global environment? Default: TRUE
#' @param run_tests Run the assertion tests? This will run the tests and error out if there's an
#' inconsistent config parameter. Default: TRUE
#' @param return_list Return a list result or just the config? Default FALSE
#' @return Depends on return_list. If FALSE (default) returns just the MBG config (a list). If True, returns a
#' named list of configs, where "config" is the usual MBG config, and "fixed_effects_config" is the config info
#' of the fixed effects
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   config <- load_config(repo = core_repo,
#'     indicator_group = indicator_group,
#'     indicator = indicator,
#'     config_name = 'config_training',
#'     covs_name = 'covs_training')
#' }
#' }
#' @rdname set_up_config
#' @importFrom assertthat is.flag is.string is.number
#' @export
set_up_config <- function(repo,
                          core_repo = repo,
                          indicator_group,
                          indicator,
                          config_name = NULL,
                          config_file = NULL,
                          covs_name = NULL,
                          covs_file = NULL,
                          post_est_only = FALSE,
                          run_date = "",
                          push_to_global_env = TRUE,
                          run_tests = TRUE,
                          return_list = FALSE) {

  ###### Block 1: Equivalent to load_config ######

  print("[1/6] Load the configs")

  ####### Logic checking for model config #######
  ## Make sure only one of config_name or config_file are not null
  if (!is.null(config_name) & !is.null(config_file)) {
    stop("You must specify just one of config_name or config_file, not both", call. = FALSE)
  }

  ## Pull config from indicator repo
  if (is.null(config_name) & is.null(config_file)) {
    message("Pulling config from default folder, since config_name and config_file are NULL")
    if (post_est_only == FALSE)
      config <- data.table::fread(paste0(repo, "/", indicator_group, "/config_", indicator, ".csv"), header = FALSE)
    ## If running analysis on existing model, use config from that model's outputs folder
    if (post_est_only == TRUE)
      config <- data.table::fread(paste0("<<<<FILEPATH REDACTED>>>>>"))
  }

  ## Pull by specific config name
  if (!is.null(config_name) & is.null(config_file)) {
    message("Pulling config from specified name")
    config <- data.table::fread(paste0(repo, "/", indicator_group, "/", config_name, ".csv"), header = FALSE)
  }
  ## Pull specified config file
  if (is.null(config_name) & !is.null(config_file)) {
    message("Pulling config from specified filepath")
    config <- data.table::fread(config_file, header = FALSE)
  }

  ####### Logic checking for covariates config #######
  ## Make sure only one of covs_name or covs_file are not null
  if (!is.null(covs_name) & !is.null(covs_file)) {
    stop("You must specify just one of covs_name or covs_file, not both", call. = FALSE)
  }

  ## Covs not pulled
  if (is.null(covs_name) & is.null(covs_file)) {
    message("Not pulling covs since covs_name and covs_file are NULL")
    covs <- NULL
  }

  ## Pull by specific covs name
  if (!is.null(covs_name) & is.null(covs_file)) {
    message("Pulling covs from specified name")
    covs <- read_covariate_config(paste0(repo, "/", indicator_group, "/", covs_name, ".csv"))
  }

  ## Pull specified covs file
  if (is.null(covs_name) & !is.null(covs_file)) {
    message("Pulling covs from specified filepath")
    covs <- read_covariate_config(covs_file)
  }

  ## For parsimony, let's make sure that the config column names are V1 and V2
  config <- data.table(config)
  if(colnames(config)[1] != "V1" & colnames(config)[2] != "V2") {
    warning("Renaming config column names to V1 and V2. Please verify that 'config' is properly built")
    colnames(config) <- c("V1", "V2")
  }


  # If a covariate .csv file exists, use that instead
  if (!is.null(covs)) {

    # Grab fixed effects & measures (and gbd fixed effects & measures) from CSV if present

    # After update to data.table 1.11.4, 'T' and 'F' are not read in as logical,
    ## but as characters, which we need to remedy here.
    ## We are assuming that the 'covs.csv' has an 'include' and 'gbd' column here
    covs[, `:=`(gbd, as.logical(gbd))]
    covs[, `:=`(include, as.logical(include))]
    covs <- subset(covs, include == T)  # Use only those where include flag set to true
    fe <- subset(covs, gbd == F)
    update_fixed_effect_config_with_missing_release(fe)
    gbd <- subset(covs, gbd == T)
    gbd[measure != "output", `:=`(measure, "covariate")]
    fixed_effects <- paste0(fe$covariate, collapse = " + ")
    fixed_effects_measures <- paste0(fe$measure, collapse = " + ")
    gbd_fixed_effects <- paste0(gbd$covariate, collapse = " + ")
    gbd_fixed_effects_measures <- paste0(gbd$measure, collapse = " + ")

    if(!("constraint" %in% names(covs))){
      fixed_effects_constraints <- paste0("c(", paste(rep(0, nrow(fe)), collapse=", "), ")")
      gbd_fixed_effects_constraints <- paste0("c(", paste(rep(0, nrow(gbd)), collapse=", "), ")")
    }
    else{
      fixed_effects_constraints <- paste0("c(", paste(unname(fe$constraint), collapse=", "), ")")
      gbd_fixed_effects_constraints <- paste0("c(", paste(unname(gbd$constraint), collapse=", "), ")")
    }

    # Remove any other versions from original config and
    # override with covariates config outputs
    all_varz <- c(
      "fixed_effects", "fixed_effects_measures", "fixed_effects_constraints",
      "gbd_fixed_effects", "gbd_fixed_effects_measures", "gbd_fixed_effects_constraints"
    )
    for(varz in all_varz) {
      if(!(varz %in% colnames(config))) {
        config <- rbindlist(list(config, data.table(V1 = varz, V2 = get(varz))))
      } else {
        config[V1 == varz, V2:= get(varz)]
      }
    }

  }


  ###### Block 2: Add fields in config that are not in the default set ######

  print("[2/6] Add fields that are in the default config set but not in user's config")

  ## Load in the default config dataset
  if (.in.package()) {
    data("default_config_values", package = packageName())
  } else {
    default_config_values <- data.table::fread(file.path(core_repo, '/mbg_central/share_scripts/common_inputs/config_values.csv'), header = TRUE, stringsAsFactors = FALSE)
  }

  ## Now, go through each of the values in `config_values` and
  ## add on all the fields that are not in the user-specified config
  config <- set_default_config_values(config, default_config_values)


  ###### Block 3: Extra parameters in config ######

  print("[3/6] Add fields that are in user's config but not in the default config set")
  message("\nAdditional covariates: ")
  extras <- config$V1[!(config$V1 %in% names(default_config_values))]
  for (extra in extras) {
    message(paste0("  ", extra, ": ", config[V1 == extra, V2] ))
  }

  ###### Block 4: Print out shapefile info from config. Resolve 'current' to fixed version date ######

  print("[4/6] Print out shapefile info from config")

  ## get the shapefile info
  m.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = config[V1 == 'modeling_shapefile_version', V2]))
  r.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = config[V1 == 'raking_shapefile_version', V2]))

  ## replace shapefile version in config (and env variable) with the actual version date
  ## if a specific date was already set, nothing changes.
  ## if 'current' had been selected, then it will be replaced by the version date that is currently symlinked to 'current'
  config[V1 == 'modeling_shapefile_version', V2 := m.sf.info$shpfile_date]
  config[V1 == 'raking_shapefile_version', V2 := r.sf.info$shpfile_date]

  ## print out the shapefile info
  message("\n\n\nSHAPEFILE VERSION INFORMATION: ")
  message(sprintf("\n--MODELING SHAPEFILE VERSION: %s -- which contains %s codes", m.sf.info$shpfile_date, toupper(m.sf.info$shpfile_type)))
  message(sprintf("\n--RAKING SHAPEFILE VERSION:   %s -- which contains %s codes\n", r.sf.info$shpfile_date, toupper(r.sf.info$shpfile_type)))


  ###### Block 5: Run tests on all the configuration variables loaded ######
  if(run_tests) {
    print("[5/6] Running simple type-assertion tests on config parameters")
    if (.in.package()) {
      data("config_tests", package = packageName())
    } else {
      config_tests <- data.table::fread(paste0(core_repo, '/mbg_central/share_scripts/common_inputs/config_tests.csv'), header = TRUE, stringsAsFactors = FALSE)
    }

    ## Test for params only in the config_tests list of params
    for (param in sort(config[, V1])) {
      cat(paste0("Testing config parameter: ", param, " "))
      if(param %in% config_tests$variable) {
        test_call_1 <- config_tests[variable == param, test_call]
        test_call_2 <- config_tests[variable == param, extra_test1]
        test_call_3 <- config_tests[variable == param, extra_test2]

        if(test_call_1 != "") {
          ## For a string in the config file, the eval-parse combo will
          ## fail to evaluate it, and so we build in this exception for that
          tryCatch(
            get(test_call_1)(ez_evparse(config[V1 == param, ], "V2")),
            error = function(e) {
              if(attributes(e)$class[[1]] == 'simpleError') {
                if (test_call_1 == "is.string") {
                  message(paste0("Assertion on ", param, " errored out because it's tested as a string. Please check for the real type manually"))
                } else {
                  stop(sprintf("%s errored with message: %s", test_call_1, geterrmessage()))
                }
              }
            }
          )
        }
        if(test_call_2 != ""  ) {
          tryCatch(
            assertthat::assert_that(eval(parse(text = test_call_2))),
            error = function(e) {
              stop(paste0("The following test failed: ", test_call_2) )
            }
          )
        }
        if(test_call_3 != ""  ) {
          tryCatch(
            assertthat::assert_that(eval(parse(text = test_call_3))),
            error = function(e) {
              stop(paste0("The following test failed: ", test_call_3) )
            }
          )
        }
        cat(" OK. \n")
      }
    }

    ## Stop if using z or poly aggregation strategies without TMB
    if(as.logical(config[V1 == "poly_ag", "V2"]) | config[V1 == "zcol_ag", "V2"] != "NULL") {
      if(!as.logical(config[V1 == "fit_with_tmb", "V2"])) {
        stop("Must use TMB when incorporating polygon or aggregated z data")
      }
      if(as.logical(config[V1 == "makeholdouts", "V2"])) {
        stop("There is aggregated data and functionality for making holdouts is
             not yet implemented. Set makeholdouts to FALSE.")
      }
      if(as.logical(config[V1 == "test", "V2"])) {
        stop("Testing with aggregated data not yet implemented. Set test to FALSE.")
      }
    }

  } else {
    warning("[5/6] Skipping over type-assertion")
  }




  ###### Final Block :  Assign all the covariates to the environment if desired ######

  if(push_to_global_env) {
    print("[6/6] Pushing config parameters into global environment")
    for (param in config[, V1]) {
      assign(param, config[V1 == param, V2], envir = globalenv())
    }
    if (!is.null(covs)) {
      assign("fixed_effects_config", fe, envir = globalenv())
      assign("gbd_fixed_effects_config", gbd, envir = globalenv())
    }
    # Processing of z config arguments
    if (zcol_ag != "NULL") {
      if (exists("z_ag_mat_file")) {
        assign("z_ag_mat", read.csv(z_ag_mat_file, header=F), envir = globalenv())
      } else {
        assign("z_ag_mat", NULL, envir = globalenv())
        assign("zcol_ag_id", NULL, envir = globalenv())
      }
      if (exists("z_map_file")) {
        assign("z_map", read_z_mapping(z_map_file), envir = globalenv())
      } else {
        assign("z_map", NULL, envir = globalenv())
      }
    } else {
      assign("zcol_ag", NULL, envir = globalenv())
      assign("z_map", NULL, envir = globalenv())
      assign("z_ag_mat", NULL, envir = globalenv())
      assign("zcol_ag_id", NULL, envir = globalenv())
    }
  } else {
    print("[6/6] Config parameters not passed into global environment")
  }

  ## Return the config data.table
  message("Saving out config...")
  if (return_list) {
    return(list("config" = config, "fixed_effects_config" = fe))
  } else {
    return(config)
  }
}


#' @title lonlat3D
#' @description  This function takes in a vector of longitude and a vector of
#' latitude and returns coordinates on the S2 sphere (globe living in
#' 3D) in (x, y, z) coords on a sphere with radius 1
#'
#' @param lon numeric vector of longitude coords
#' @param lat numeric vector of latitude coords
#'
#' @return 3 column numeric matrix where each row is a (x,y,z) of the
#'   transformed (long, lat) coords

lonlat3D <- function(lon,lat){
  cbind(cos((lon/180)*pi)*cos((lat/180)*pi),
        sin((lon/180)*pi)*cos((lat/180)*pi),
        sin((lat/180)*pi))
}


#' @title Get output regions
#' @description This function takes in a directory where mbg modeling has stored
# outputs (*cell_pred* objects) and infers the regions specified in
# the model
#'
#' @description Determines modeling regions from written output dir objects
#'
#' @param in_dir directory path containing completed mbg cell_pred objects
#'
#' @return A vector string of region names

get_output_regions <- function(in_dir) {
  return(unique(stringr::str_match(list.files(in_dir, pattern = paste0('_cell_draws_eb_')),
                                   '_cell_draws_[^_]+_[^_]+_(.*)_[0-9].RData')[,2]))
}


#' @title rank_draws
#' @description function to transform values to ranks by draw
#' @param df data frame with columns for each draw
#' @param 'high' means a high value is rank 1, 'low' means a
#' low value should be rank 1
#' @param columns: vector of column names that contain the draws
#' @export

rank_draws <- function(df, ordr, columns) {
  library(data.table)

  # Ensure data.table class
  df <- as.data.table(df)
  # Separate into df with draws and the rest
  df1 <- df[, setdiff(names(df), columns), with = FALSE]
  df2 <- df[, columns, with = FALSE]

  # Transform draws to ranks
  if (ordr == 'high') {
    descending <- TRUE
  } else {
    descending <- FALSE
  }
  df2 <- as.data.table(apply(df2, 2, order, decreasing = descending))

  # Summarize ranks
  df2 <- df2[, `:=`(median=apply(df2, 1, median),
                    lci=apply(df2, 1, quantile, 0.025),
                    uci=apply(df2, 1, quantile, 0.975))]

  # Bind to make a single df
  df <- cbind(df1, df2)
  return(df)

}

#' Return value if provided in function, else return global of same name.
#'
#' Returns the value named \code{name} from the calling function's environment. If that
#' results in an error (e.g., beause the value was not provided) then return an identically
#' named value from the global environment. If no such value exists in the global
#' environment then error.
#'
#' @param name character name of the value to return.
#'
#' @return the value.
#'
#' @export
use_global_if_missing <- function(name) {
  tryCatch(
    {
      return(get(name, pos = parent.frame(1)))
    },
    error = function(e) {
      if (name %in% names(.GlobalEnv)) {
        return(.GlobalEnv[[name]])
      } else {
        stop(sprintf("Variable %s not provided to function and not available in global environment", name))
      }
    }
  )
}
