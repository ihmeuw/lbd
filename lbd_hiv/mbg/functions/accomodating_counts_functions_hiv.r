# These are to be functions that will bue used after the fractional raking has taken place. This script will get called from launch.




############# edit of post_load_combine_save to accomodate the counts data created in the fractional raking
#####################################################
#####################################################

post_load_combine_save_counts <- function(regions    = strata,
                                          summstats  = c('mean','cirange','upper','lower'),
                                          raked      = c('raked','unraked', 'raked_c'),
                                          rf_table   = TRUE,
                                          run_summ   = TRUE,
                                          indic      = indicator,
                                          ig         = indicator_group,
                                          sdir       = sharedir){

  message(paste0("indic: ", indic))
  message(paste0("ig: ", ig))

  rake_addin <- character()
  if ("unraked" %in% raked) {
    lookup_dir <- paste0(sprintf('<<<< FILEPATH REDACTED >>>>'))
    ur <- length(grep(paste0(indic, ".*unraked.*raster.tif"), list.files(lookup_dir)))
    if (ur > 0) rake_addin <- c(rake_addin, unraked = "_unraked")
    if (ur == 0) rake_addin <- c(rake_addin, unraked = "")
  }

  if ("raked" %in% raked) {
    rake_addin <- c(rake_addin, raked = "_raked")
  }

  if ("raked_c" %in% raked) {
    rake_addin <- c(rake_addin, raked = "_raked_c")
  }
  # loop through and combine all rasters
  message("\nCombining rasters...")
  for (rake in rake_addin) {
    message(names(rake_addin)[which(rake_addin == rake)])
    rr <- rake
    for (ss in summstats) {
      message(paste0('  ',ss))
      rlist <- list()
      for (reg in regions) {
        message(paste0('    ',reg))
        rlist[[reg]] <-
          brick(sprintf('<<<< FILEPATH REDACTED >>>>'))
      }
      if (length(rlist) > 1) rlist <- do.call(raster::merge,unname(rlist)) else rlist <- rlist[[1]]
      if (ss == 'cirange') ssname = 'range' else ssname = ss # naming convention
      save_post_est(rlist, 'raster',
                    paste0(ssname,rr,'_raster'),
                    indic)
    }
  }

  # do rf also
  if (rf_table) {
    message('RF table')
    rflist <- list()
    for (reg in regions) {
      rflist[[reg]] <-
        read.csv(sprintf('<<<< FILEPATH REDACTED >>>>'))
    }
    save_post_est(do.call(rbind.fill,rflist),'csv','rf', indic)
  }

  # make a run summary graph
  if (run_summ) {
    graph_run_summary(run_date = run_date,
                      indicator_group = ig,
                      indicator = indic)
  }
}


######################## edit of combine_aggregation to accomodate the counts data created in the fractional raking
#####################################################################
#####################################################################

combine_aggregation_counts <- function(rd = run_date,
                                       indic = indicator,
                                       ig = indicator_group,
                                       ages,
                                       regions,
                                       holdouts,
                                       raked,
                                       dir_to_search = NULL,
                                       delete_region_files = T,
                                       merge_hierarchy_list = F) {

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
    dir_to_search <- paste0('<<<< FILEPATH REDACTED >>>>')
  }

  message("Combining aggregation results...")

  for (rake in raked) {
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

          if (merge_hierarchy_list == T) {
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

          ad0[[reg]] <- admin_0
          ad1[[reg]] <- admin_1
          ad2[[reg]] <- admin_2
          sp_h[[reg]] <- sp_hierarchy_list

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
  fin_files_to_delete <- list.files(dir_to_search, pattern = "fin_agg_", full.names = T)
  unlink(fin_files_to_delete)
}




################################  edit of summarize_admins to accomodate the counts data created in the fractional raking
########################################################################
########################################################################

summarize_admins_counts <- function(ind = indicator,
                                    ig = indicator_group,
                                    summstats = c("mean", "lower", "upper", "cirange"),
                                    raked    = c("raked", "raked_c", "unraked"),
                                    ad_levels = c(0,1,2),
                                    file_addin = NULL,
                                    ...) {

  sharedir       <- sprintf('<<<< FILEPATH REDACTED >>>>')
  input_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  output_dir <- paste0('<<<< FILEPATH REDACTED >>>>')
  dir.create(output_dir, recursive = T, showWarnings = F)

  # Convert raked to character
  rr <- raked

  # If file_addin present, use it
  if (!is.null(file_addin)) file_addin <- paste0("_", file_addin)
  if (is.null(file_addin)) file_addin <- ""

  # Summarize and save admin preds
  for (rake in rr) {
    load(paste0(input_dir, ind, "_", rake, "_admin_draws_eb_bin0_0.RData"))
    sp_hierarchy_list <- mutate_if(sp_hierarchy_list, is.factor, as.character)
    sp_hierarchy_list <- mutate_at(sp_hierarchy_list, grep('_CODE', names(sp_hierarchy_list), value = T), as.numeric)

    for (ad in ad_levels) {
      message(paste0("Summarizing ", ind, ": admin ", ad, " (", rake, ")"))
      ad_summary_table <- make_admin_pred_summary(admin_pred = get(paste0("admin_", ad)),
                                                  sp_hierarchy_list,
                                                  summary_stats = summstats,
                                                  ...)
      fwrite(ad_summary_table,
             file = paste0(output_dir, ind, "_admin_", ad, "_", rake, file_addin, "_summary.csv"))
    }
  }
}



