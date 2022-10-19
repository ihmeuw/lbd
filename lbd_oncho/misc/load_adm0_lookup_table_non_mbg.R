#' @title Load the GAUL lookup table
#'
#' @description Loads the most recent version of the lookup table that links
#'   ADM0 codes with other identifiers such as GBD location IDs, MBG modeling
#'   regions, and ISO codes, among others.
#'
#' @return Returns a data.frame of the ADM0 lookup table. If package data.table
#'   is loaded, also converts the lookup table to a data.table
#' @export
load_adm0_lookup_table <- function() {
  # Define consistent path to the ADM0 lookup table
  lookup_table_filepath <- <<<< FILEPATH REDACTED >>>>
  # If data.table is loaded, read in as a data.table
  if ("package:data.table" %in% search()) {
    lookup_table <- fread(lookup_table_filepath)
  } else {
    # Otherwise, read in as a data.frame
    lookup_table <- read.csv(lookup_table_filepath)
  }
  # Set ISO codes as lowercase for easy lookup
  lookup_table$iso3 <- tolower(lookup_table$iso3)
  return(lookup_table)
}
