## ##########################################################
##  Purpose: Check for results, by location, of a given process
##          This will loop through the list of locations and check for the existence of files within the directory
##          e.g. if locations = "AFG","CAN" and prefix = "draws" and postfix = ".csv", we would look for "drawsAFG.csv", "drawsCAN.csv", etc.
##          Sleeping for X seconds each time until they are written (if they ever get written)
##  Options:
##    locations: list of locations to look through
##    check_dir: what directory to look into
##    prefix: what is the prefix of the filename (e.g. "draws_")
##    postfix: what is the filename after the location (e.g. "_HIVMortonART.csv")
##    sleep: How many seconds to sleep between checks?

check_loc_results <- function(locations, check_dir, prefix = "", postfix, sleep = 300) {
  counter <- 0
  time_counter <- 0
  while (counter == 0) {
    inner_counter <- 0
    missing_list <- ""
    for (loc in locations) {
      # Check duration draws (from no-ART process)
      if (file.exists(paste0(check_dir, "/", prefix, loc, postfix))) {
        inner_counter <- inner_counter + 1
      } else {
        missing_list <- paste0(missing_list, " ", loc)
      }
    }

    if (length(locations) == inner_counter) {
      print("All results are present")
      counter <- 1
    } else {
      print(paste0("Have ", inner_counter, " results: expecting ", length(locations), " ", Sys.time()))
      if (inner_counter > (length(locations) * .75)) {
        print(paste0("Still Missing: ", missing_list))
      }
      time_counter <- time_counter + 1
      Sys.sleep(sleep)
    }
  }
}


check_loc_results_kick <- function(locations, check_dir, prefix = "", postfix, sleep = 300) {
  counter <- 0
  time_counter <- 0
  while (counter == 0) {
    inner_counter <- 0
    missing_list <- ""
    for (loc in locations) {
      # Check duration draws (from no-ART process)
      if (file.exists(paste0(check_dir, "/", prefix, loc, postfix))) {
        inner_counter <- inner_counter + 1
      } else {
        missing_list <- paste0(missing_list, " ", loc)
      }
    }

    if (length(locations) == inner_counter) {
      print("All results are present")
      counter <- 1
      next
    } else {
      print(paste0("Have ", inner_counter, " results: expecting ", length(locations), " ", Sys.time()))
      if (inner_counter > (length(locations) * .75)) {
        print(paste0("Still Missing: ", missing_list))
      }
      time_counter <- time_counter + 1
      Sys.sleep(sleep)
    }
  }
}
