#' @title Make admin prediction summary
#' @description Take an admin pred object and an sp_hierarchy list and generate a sorted,
#' cleaned table for a given set of summary statistics
#'
#' @param admin_pred the admin draws object (e.g. admin_1, etc)
#' @param sp_hierarchy_list spatial hierarchy object from the admin draws RData files
#' @param summary_stat summary statistic to apply
#' @param ... any other arguments to pass to the `summary_stats` functions
#'            (note that currently will be passed to all functions)
#' @return data table of summary stats and admin info
#' @examples
#' \dontrun{
#' # load("indicator_raked_admin_draws_eb_bin0.RData") # Replace with your admin draws
#' make_admin_pred_summary(
#'   admin_2,
#'   sp_hierarchy_list,
#'   c("mean", "cirange", "upper", "lower")
#' )
#' }
#' @export
make_admin_pred_summary <- function(admin_pred,
                                    sp_hierarchy_list,
                                    summary_stats = "mean",
                                    ...) {

  ### Get set up
  str_match <- stringr::str_match

  # Split up your input data
  ad_code <- subset(admin_pred, select = grep("ADM[0-9]_CODE", names(admin_pred)))
  year <- subset(admin_pred, select = "year")
  pop <- subset(admin_pred, select = "pop")
  draws <- subset(admin_pred, select = grep("V[0-9]*", names(admin_pred)))

  # Get the admin level
  ad_code_name <- names(ad_code)
  ad_level <- as.numeric(str_match(ad_code_name, "ADM([0-9])_CODE")[, 2])

  # Get all admin levels
  all_admin_levels <- as.numeric(str_match(names(sp_hierarchy_list), "ADM([0-9])_CODE")[, 2])
  all_admin_levels <- unique(all_admin_levels)[!is.na(unique(all_admin_levels))]

  ### Make summary
  summ_list <- lapply(summary_stats, function(ss) {
    apply(draws, 1, ss, ...)
  })

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
      drop_cols <- names(sp_hierarchy_list)[grepl(paste0("ADM", d), names(sp_hierarchy_list))]
      sp_hierarchy_list <- subset(sp_hierarchy_list, select = !(names(sp_hierarchy_list) %in% drop_cols))
    }
  }
  sp_hierarchy_list <- unique(sp_hierarchy_list)

  output_df <- merge(sp_hierarchy_list, output_df, all.y = T, all.x = F)
  
  # Pass other, custom cols (such as aggregated stacker predictions) through to output
  non_draw_cols <- grep("V", names(admin_pred), invert = T, value = T)
  if(length(non_draw_cols[!(non_draw_cols %in% c(names(output_df), "year", "pop"))]) > 0){
    output_df <- merge(output_df, admin_pred[, non_draw_cols[!(non_draw_cols %in% c("pop"))], with = F], all = TRUE)
  }

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

  return(output_df)
}
