#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param indicator PARAM_DESCRIPTION
#' @param indicator_group PARAM_DESCRIPTION
#' @param run_date PARAM_DESCRIPTION
#' @param save_objs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname prep_postest
#' @export
prep_postest <- function(indicator,
                         indicator_group,
                         run_date,
                         reg,
                         age,
                         sex,
                         save_objs) {
  
  # Save a list of objects in a standard location for parallel scripts to pull from
  main_dir <- paste0(<<<< FILEPATH REDACTED >>>> "/output/", run_date, "/")
  temp_dir <- paste0(main_dir, "temp_post_est/")
  temp_file <- paste0(temp_dir, "post_est_temp_objs_bin", age, "_sex", sex, "_", reg, ".RData")
  dir.create(temp_dir, showWarnings = F)
  save(list = save_objs, file = temp_file)
}