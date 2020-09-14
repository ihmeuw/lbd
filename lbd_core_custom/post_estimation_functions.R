
## summarize_admins ################################################

#' Function to summarize admin_pred objects
#'
#' This is a wrapper for `make_admin_pred_summary()`
#'
#' @param ind indicator
#' @param ig indicator_group
#' @param summstats Summary statistics (functions) to compute.
#'                  Order will be the order they are written in csv
#'                  This is passed to `make_admin_pred_summary()`
#' @param raked Raked (T), unraked (F), or both (`c(T,F)`)?
#' @param ad_levels Admin levels to summarize (0, 1, and 2 by default)
#' @return Writes csv files to `'<<<< FILEPATH REDACTED >>>>'`
#' @examples
#' summarize_admins(summstats = c("mean", "lower", "upper", "cirange"),
#'                  ad_levels = c(0,1,2),
#'                  raked     = c(T,F))

summarize_admins <- function(ind = indicator,
                             ig = indicator_group,
                             summstats = c("mean", "lower", "upper", "cirange"),
                             raked = c(T,F),
                             ad_levels = c(0,1,2),
                             file_addin = NULL,
                             measure = 'prevalence',
                             metrics = c('rates', 'counts'),
                             eti = 'none',
                             ...) {
  
  sharedir       <- '<<<< FILEPATH REDACTED >>>>'
  input_dir <- paste0(sharedir, "/output/", run_date, "/")
  output_dir <- paste0(input_dir, "/pred_derivatives/admin_summaries/")
  dir.create(output_dir, recursive = T, showWarnings = F)
  
  # Convert raked to character
  rr <- character()
  if (T %in% raked) rr <- c(rr, paste0("raked_", measure))
  if (F %in% raked) rr <- c(rr, "unraked")
  
  # If file_addin present, use it
  if (!is.null(file_addin)) file_addin <- paste0("_", file_addin)
  if (is.null(file_addin)) file_addin <- ""
  
  # Summarize and save admin preds
  for (et in eti){
    for (metric in metrics) {
      for (rake in rr) {
        load(paste0(input_dir, ind, "_", ifelse(et == 'none', '', paste0(et, '_')), rake, ifelse(metric == "counts", "_c", ""), "_admin_draws_eb_bin0_0.RData"))
        sp_hierarchy_list <- mutate_if(sp_hierarchy_list, is.factor, as.character)
        sp_hierarchy_list <- mutate_at(sp_hierarchy_list, grep('_CODE', names(sp_hierarchy_list), value = T), as.numeric)
        
        for (ad in ad_levels) {
          message(paste0("Summarizing ", ind, ": admin ", ad, " (", rake, ifelse(et == 'none', '', paste0(' ', et)), ' ', metric, ")"))
          ad_summary_table <- make_admin_pred_summary(admin_pred = get(paste0("admin_", ad)),
                                                      sp_hierarchy_list,
                                                      summary_stats = summstats,
                                                      ...)
          fwrite(ad_summary_table,
                 file = paste0(output_dir, ind, ifelse(et == 'none', '', paste0('_', et)), ifelse(metric == "counts", "_c", ""), "_admin_", ad, "_", rake, file_addin, "_summary.csv"))
        }
      }
    }
  }
}
  
summarize_admins_tridalys <- function(ind = indicator,
                             ig = indicator_group,
                             summstats = c("mean", "lower", "upper", "cirange"),
                             raked = c(T,F),
                             ad_levels = c(0,1,2),
                             file_addin = NULL,
                             measure = 'prevalence',
                             metrics = c('rates', 'counts'),
                             input_dir,
                             output_dir,
                             ...) {
  
  sharedir       <- '<<< FILEPATH REDACTED >>>'
  dir.create(output_dir, recursive = T, showWarnings = F)
  
  # Convert raked to character
  rr <- character()
  if (T %in% raked) rr <- c(rr, paste0("raked_", measure))
  if (F %in% raked) rr <- c(rr, "unraked")
  
  # If file_addin present, use it
  if (!is.null(file_addin)) file_addin <- paste0("_", file_addin)
  if (is.null(file_addin)) file_addin <- ""
  
  # Summarize and save admin preds
  for (metric in metrics) {
    for (rake in rr) {
      load(paste0(input_dir, ind, "_", rake, ifelse(metric == "counts", "_c", ""), "_admin_draws_eb_bin0_0.RData"))
      sp_hierarchy_list <- mutate_if(sp_hierarchy_list, is.factor, as.character)
      sp_hierarchy_list <- mutate_at(sp_hierarchy_list, grep('_CODE', names(sp_hierarchy_list), value = T), as.numeric)
      
      for (ad in ad_levels) {
        message(paste0("Summarizing ", ind, ": admin ", ad, " (", rake, ")"))
        ad_summary_table <- make_admin_pred_summary(admin_pred = get(paste0("admin_", ad)),
                                                    sp_hierarchy_list,
                                                    summary_stats = summstats,
                                                    ...)
        fwrite(ad_summary_table,
               file = paste0(output_dir, ind, ifelse(metric == "counts", "_c", ""), "_admin_", ad, "_", rake, file_addin, "_summary.csv"))
      }
    }
  }
}

## make_admin_pred_summary ###################################################

#' Take an admin pred object and an sp_hierarchy list and generate a sorted,
#' cleaned table for a given set of summary statistics
#'
#' @param admin_pred the admin draws object (e.g. admin_1, etc)
#' @param sp_hierarchy_list spatial hierarchy object from the admin draws RData files
#' @param summary_stat summary statistic to apply
#' @param ... any other arguments to pass to the `summary_stats` functions
#'            (note that currently will be passed to all functions)
#' @return data table of summary stats and admin info
#' @examples
#' # load("indicator_raked_admin_draws_eb_bin0.RData") # Replace with your admin draws
#' make_admin_pred_summary(admin_2,
#'                         sp_hierarchy_list,
#'                         c("mean", "cirange", "upper", "lower"))

make_admin_pred_summary    <- function(admin_pred,
                                       sp_hierarchy_list,
                                       summary_stats = 'mean',
                                       ...){
  
  ### Get set up
  str_match <- stringr::str_match
  
  # Split up your input data
  ad_code <- subset(admin_pred, select = grep("ADM[0-9]_CODE", names(admin_pred))) %>%
    as.data.table()
  year <- subset(admin_pred, select = "year") %>%
    as.data.table()
  pop <- subset(admin_pred, select = c("pop", grep("ADM[0-9]_CODE", names(admin_pred), value = TRUE), 'year')) %>%
    as.data.table()
  draws <- subset(admin_pred, select = grep("V[0-9]*", names(admin_pred))) %>%
    as.data.table()
  
  # Get the admin level
  ad_code_name <- names(ad_code)
  ad_level <- as.numeric(str_match(ad_code_name,"ADM([0-9])_CODE")[,2])
  
  # Get all admin levels
  all_admin_levels <- as.numeric(str_match(names(sp_hierarchy_list), "ADM([0-9])_CODE")[,2])
  all_admin_levels <- unique(all_admin_levels)[!is.na(unique(all_admin_levels))]
  
  ### Make summary
  summ_list <- lapply(summary_stats, function(ss) {
    apply(draws, 1, ss, ...)})
  
  if (length(summ_list) > 1) {
    output_df <- as.data.table(cbind(year, ad_code, do.call(cbind, summ_list)))
    names(output_df)[grep("V[0-9]*", names(output_df))] <- summary_stats
  } else if (length(summ_list) == 1) {
    output_df <- as.data.table(cbind(year, ad_code, summ_list[[1]]))
    names(output_df)[grep("V[0-9]*", names(output_df))] <- summary_stats
  }
  
  ### Add on identifiers
  
  # Drop smaller admin levels
  drop_levels <- all_admin_levels[all_admin_levels > ad_level]
  
  if (length(drop_levels) > 0) {
    for (d in drop_levels) {
      drop_cols <- names(sp_hierarchy_list)[grepl(paste0("ADM",d), names(sp_hierarchy_list))]
      sp_hierarchy_list <- subset(sp_hierarchy_list, select = !(names(sp_hierarchy_list) %in% drop_cols))
    }
  }
  sp_hierarchy_list <- unique(sp_hierarchy_list)
  
  output_df <- merge(sp_hierarchy_list, output_df, all.y = T, all.x = F) %>%
    merge(pop)
  
  # Clean up & sort
  order_vars <- c(paste0("ADM", 0:ad_level, "_NAME"), "year")
  setorderv(output_df, order_vars)
  
  # Get col order
  ad_col_order <- sort(names(output_df)[grep("AD*_*", names(output_df))])
  cols_to_order <- c(ad_col_order, "region", "year", summary_stats)
  cols_to_order <- cols_to_order[cols_to_order %in% names(output_df)]
  if (length(cols_to_order) < length(names(output_df))) {
    other_cols <- names(output_df)[!(names(output_df) %in% cols_to_order)]
    cols_to_order <- c(cols_to_order, other_cols)
  }
  
  setcolorder(output_df, cols_to_order)
  
  #add on GBD region, superregion, and ISO3
  region_list <- fread('<<<< FILEPATH REDACTED >>>>') %>%
    select('iso3', 'spr_reg_nm', 'reg_name', 'gadm_geoid') %>%
    rename('gbd_super_region' = 'spr_reg_nm', 'gbd_region' = 'reg_name', 'ADM0_CODE' = 'gadm_geoid')
  
  output_df <- merge(output_df, region_list)
  
  return(output_df)
  
}
