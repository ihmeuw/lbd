#' @title Returns the stage master list
#'
#' @description Loads from a pre-saved data file.
#'
#' @export
load_stage_list <- function() {
  if (.in.package()) {
    data("stage_master_list")
    stage_master_list
  } else {
    data.table::fread(<<<< FILEPATH REDACTED >>>>)
  }
}
