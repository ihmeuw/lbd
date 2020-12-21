# HEADER ------------------------------------------------------------------
# Purpose: Functions to speed up shapefile processing
#
# Details: These functions are meant to allow synchronization of our main
#          shapefile directory, where all shapefiles are stored in .shp
#          format, and a new directory where shapefiles are stored in RDS
#          format.  Loading an RDS file is approximately 40,000 times faster
#          than loading a shapefile with readOGR() in R.
#
#          An effort has been made to creating a robust one-way
#          synchronization environment with file locking to help avoid
#          corruption of the rds-formatted shapefiles, but this method
#          may result in some errors that would be absent from the original
#          readOGR()-based method (e.g. if synchronization fails).
#
#          None of this code should modify the main shapefile directory in
#          any way, so the objects of reference should be unchanged.
#
# Sample usage:
#  synchronize_shapefile_directories(cores = 12) #To ensure library up to date
#  shp <- fast_load_shapefile("ALB_2005_MICS.rds")
#**************************************************************************


#####################################################################
## Functions to read and write shapefiles in .rds format
#####################################################################

# Synchronize the two shapefile directories
# Only new or more recently modified files will be changed (for speed)
synchronize_shapefile_directories <- function(cores = NULL,
                                              verbose = F,
                                              in_dir = '<<<< FILEPATH REDACTED >>>>',
                                              out_dir = '<<<< FILEPATH REDACTED >>>>') {

  # This function performs a one-way sync from in_dir --> out_dir
  # Run this before using any of the fast shapefile functions

  message("Checking to see if shapefile synchronization needed...")

  # Load one function from stringr
  str_match <- stringr::str_match

  # Make list of input files from in_dir with last modified time
  in_files <- list.files(in_dir, pattern = '.shp$', full.names = T)
  in_files <- as.data.table(in_files) %>% setnames(., "in_files", "filename_in")
  in_files[, last_modified_in := file.mtime(filename_in)]
  in_files[, shape := str_match(filename_in, "//(.*).shp$")[,2]]

  # Make list of output files from out_dir with last modified time
  out_files <- list.files(out_dir, pattern = '.rds$', full.names = T)
  out_files <- as.data.table(out_files) %>% setnames(., "out_files", "filename_out")
  out_files[, last_modified_out := file.mtime(filename_out)]
  out_files[, shape := str_match(filename_out, "//(.*).rds$")[,2]]

  out_files <- subset(out_files, !grepl("_lockfile", shape))

  # create a single master list
  in_files <- merge(in_files, out_files, all.x = T, all.y = F, by = "shape")

  # Get list of files in the in_dir but not the out_dir
  files_to_synchronize <- subset(in_files, is.na(filename_out) | last_modified_out < last_modified_in)
  files_to_synchronize <- unique(files_to_synchronize$filename_in)

  if(length(files_to_synchronize) > 0) {

    # Simplify shapes
    message(paste0("\nSynchronizing ", length(files_to_synchronize), " shapefiles..."))

    if (is.null(cores)) {
      results <- lapply(files_to_synchronize, factory(save_shapefile_as_rds),
                        in_dir = in_dir, out_dir = out_dir, verbose = verbose)
    } else if (!is.null(cores)) {
      # Set multithreading to serial for `mclapply()`:
      set_serial_threads()
      results <- mclapply(files_to_synchronize, factory(save_shapefile_as_rds),
                          in_dir = in_dir, out_dir = out_dir, verbose = verbose,
                          mc.cores = cores)
      # Return to multithreading (if any):
      set_original_threads()
    }

    # Make & format a results table
    pull_result <- function(x) x[[1]] %>% as.character

    out_table <- results %>%
                  lapply(., function(x) as.character(x[[1]])) %>%
                  unlist %>%
                  as.data.table %>%
                  setnames(., c("."), c("result")) %>%
                  .[, shape := str_match(files_to_synchronize, "//(.*).shp$")[,2]] %>%
                  setcolorder(., c("shape", "result"))

    out_table <- subset(out_table, result != "success")

    if (nrow(out_table) > 0) {
      warning(paste0("The following shapefiles failed to synchronize:\n   ",
               paste(out_table$shape, collapse = "\n   "),
               "\nSee the output of this function for a summary of errors."))
      return(out_table)
    } else {
      message(paste0("Successfully synchronized ", length(files_to_synchronize), " shapes."))
      return_table <- as.data.table(c(paste0("Successfully synchronized ",
                                            length(files_to_synchronize),
                                            " shapes."),
                                    files_to_synchronize))
      return(return_table)
    }

  } else {
    message("No files to synchronize")
    return(as.data.table("No files to synchronize"))
  }
}

# Function to convert from .shp to RDS
save_shapefile_as_rds <- function(shapefile_path,
                                  verbose = F,
                                  in_dir = '<<<< FILEPATH REDACTED >>>>',
                                  out_dir = '<<<< FILEPATH REDACTED >>>>') {

  # Load one function from stringr
  str_match <- stringr::str_match

  shape <- str_match(shapefile_path, "//(.*).shp$")[2]

  lock_file(shape)

  if(verbose == T) {
    message(paste0("     ", shape))
  }

  in_dir_no_slash <- gsub('/$','',in_dir)

  the_shape <- try(readOGR(dsn = in_dir_no_slash, layer = shape), silent = T)

  if (is.error(the_shape)) {
    unlock_file(shape)
    return(the_shape)
  } else {
    saveRDS(the_shape, file = paste0(out_dir, shape, ".rds"))
    unlock_file(shape)
    return("success")
  }
}

# Simple function to quickly load a shapefile
fast_load_shapefile <- function(shape,
                                fast_shapefile_dir = '<<<< FILEPATH REDACTED >>>>') {

  wait_for_lock(shape)
  return(readRDS(paste0(fast_shapefile_dir, shape, ".rds")))
}

#####################################################################
## Functions to implement a locking system
#####################################################################

# Lock a shapefile (e.g. while editing / updating)
lock_file <- function(shp,
                      fast_shapefile_dir = '<<<< FILEPATH REDACTED >>>>') {

  lockfile <- paste0(fast_shapefile_dir, ".",shp,"_lockfile.rds")

  new_lock <- data.table(user = Sys.info()['user'],
                         time = Sys.time(),
                         shape = shp)

  if (check_for_lock(shp) == T) wait_for_lock(shp = shp)

  if (file.exists(lockfile)) {
    locks <- readRDS(lockfile)
    locks <- rbind(locks, new_lock)
  } else if (!file.exists(lockfile)) {
    locks <- new_lock
  }

  invisible(saveRDS(locks, file = lockfile))
}

# Wrapper for check_for_lock that checks every `interval` seconds up to `tries` attempts
#   If unable to obtain the lock in that time period, will throw a `stop()`
#   Useful in case someone is trying to modify a shapefile that you want to read/write
wait_for_lock <- function(shp,
                        fast_shapefile_dir = '<<<< FILEPATH REDACTED >>>>',
                        tries = 20,
                        interval = 30) {

  for (i in 1:tries) {
     if (check_for_lock(shp) == F) invisible(return())
     Sys.sleep(interval)
  }

  lockfile <- paste0(fast_shapefile_dir, ".",shp,"_lockfile.rds")
  locks <- readRDS(lockfile)
  locks[shape == shp, user]

  stop(paste0("shape ", shp,
              " locked by user ", locks[shape == shp, user],
              " at ", locks[shape == shp, time]))
}

# Check to see if the file is locked
check_for_lock <- function(shp,
                           fast_shapefile_dir = '/<<<< FILEPATH REDACTED >>>>') {

  lockfile <- paste0(fast_shapefile_dir, ".",shp,"_lockfile.rds")

  if (!file.exists(lockfile)) return(FALSE)

  locks <- readRDS(lockfile)
  if (shp %in% locks$shape) {
      return(TRUE)
    } else if (!(shp %in% locks$shape)) {
      return (FALSE)
    }
}

# Unlock a file (e.g. if done with updating it)
unlock_file <- function(shp,
                        fast_shapefile_dir = '<<<< FILEPATH REDACTED >>>>') {

  lockfile <- paste0(fast_shapefile_dir, ".",shp,"_lockfile.rds")

  if (!(file.exists(lockfile))) return(invisible())

  locks <- readRDS(lockfile)

  if (locks[shape == shp, user] != Sys.info()["user"]) {
    stop(paste0("this file is locked by another user (",
                locks[shape == shp, user], ")"))
  } else {
    locks <- subset(locks, shape != shp)
    if(nrow(locks) == 0) {
      unlink(lockfile)
    } else {
      invisible(saveRDS(locks, file = lockfile))
    }
  }
  return(invisible())
}

#####################################################################
## Error handling functions
#####################################################################

is.error <- function(x) inherits(x, "try-error")

# Additional error handling functions from
#   https://stackoverflow.com/questions/4948361/

factory <- function(fun)
    function(...) {
        warn <- err <- NULL
        res <- withCallingHandlers(
            tryCatch(fun(...), error=function(e) {
                err <<- conditionMessage(e)
                NULL
            }), warning=function(w) {
                warn <<- append(warn, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        list(res, warn=warn, err=err)
    }

.has <- function(x, what) !sapply(lapply(x, "[[", what), is.null)
hasWarning <- function(x) .has(x, "warn")
hasError <- function(x) .has(x, "err")
isClean <- function(x) !(hasError(x) | hasWarning(x))
value <- function(x) sapply(x, "[[", 1)
cleanv <- function(x) sapply(x[isClean(x)], "[[", 1)


#####################################################################
## Functions relating to offical world shapefiles
#####################################################################
## get_admin_shape_dir
#' @title Return path to admin shapes directory
#'
#' @description
#' Returns path to the official LBD administrative shape file directory. This
#' actually includes non-shape data that is important for mapping, notably
#' the standard link table and standard id raster.
#'
get_admin_shape_dir <- function(version = "current") {
    paste0("<<<< FILEPATH REDACTED >>>>", version, "/")
}

## get_admin_shapefile
#' @title Return path to world admin shapefile at specified admin level
#'
#' @description
#' Returns path to world admin shapefile given \code{admin_level}. stop()s if
#' no file exists at that admin_level. Defaults to returning the ".shp" file
#' path, but will substitute any other file \code{suffix} provided.
#'
#' @param admin_level Valid admin level we have a shapefile for. Current 0/1/2.
#'
#' @param suffix '.shp' by default, provide any other suffix to e.g., get the .dbf file
#' associated with the admin shapefile.
#'
#' @param raking boolean, default F. If TRUE pulls subnational raking shapefile.
#'
#' @examples
#' get_admin_shapefile(2)
#' get_admin_shapefile(2, suffix = '.dbf')
#'

get_admin_shapefile <- function(admin_level = 0, suffix = ".shp", type = "admin", version = 'current', raking = F) {
    if (raking) type = "raking"  # backwards compatibility

    base_dir <- get_admin_shape_dir(version)

    if (type == "admin" ) {
        path <- paste0(base_dir, "lbd_standard_admin_", admin_level, suffix)
    } else if (type == "raking") {
        path <- paste0(base_dir, "lbd_standard_raking", suffix)
    } else if (type == "disputed_mask") {
        path <- paste0(base_dir, "lbd_disputed_mask", suffix)
    } else {
        stop(paste("Unknown admin shapefile type '", type, "'"))
    }

    if (!file.exists(path)) {
        stop(paste("Could not locate admin shapefile (", path, ")"))
    }
    return(path)
}

## is_admin_shapefile_string
#' @title Return logical indicating if string is an admin shapefile version string.
#'
#' @description
#' Checks strings to see if they are equal to "current" or match a YYYY_MM_DD format.
#' If so, returns TRUE. Else FALSE. No checking is done to ensure that the date string
#' has a corresponding subdirectory in the admin_shapefiles directory.
#'
#' @param s String the string to check.
#'
is_admin_shapefile_string <- function(s) {
  if (s == "current") {
      TRUE
  } else if (length(grep("^\\d{4}_\\d{2}_\\d{2}$", s))) {
      TRUE # "2019_02_27" or another admin shapefile release date
  } else {
    FALSE
  }
}



## detect_adm_shapefile_date_type() ################################################

#' Detects whether the admin shapefile you are using is gaul or gadm
#' and returns the version/date of the shapefile even if the 'current'
#' version was specified
#'
#' @param shpfile_path path to admin shapefile we want to learn about
#'
#' @return two element named list:
#'   1) list element 'shpfile_type' contains string: either 'gadm' or 'gaul'
#'   2) list element 'shpfile_date' contains the actual date of the
#'   shapefile, even if version='current' was used in shpfile_path
#'
#' @examples
#' detect_shapefile_type(get_admin_shapefile(version='current'))

detect_adm_shapefile_date_type <- function(shpfile_path = get_admin_shapefile(version = modeling_shapefile_version)){

  ## resolve the symlink to the path if version='current'
  if(grepl(pattern = 'current', x = shpfile_path)){
    resolve.command <- paste("readlink -f", shpfile_path)
    full.path <- system(resolve.command, intern = TRUE)
  }else{
    full.path <- shpfile_path ## this is faster than resolving if it's an option...
  }

  ## grab the date from the full path
  sf.date <- strsplit(full.path, '/')[[1]][6]

  ## determine if shpfile.date pertains to gaul or gadm

  ## assume versions dated on and after Sept. 1, 2018
  ## (transition.date) are GADM, while dates before transition.date
  ## are GAUL unless otherwise coded in as an exception
  transition.date <- "2018_09_01"
  gaul.exceptions <- c() ## dates after transition.date that are actually GAUL
  gadm.exceptions <- c('2018_08_01') ## dates before transition.date are actually GADMw

  ## determine gaul or gadm
  if(sf.date >= transition.date){
    sf.type <- 'gadm'
  } else{
    sf.type <- 'gaul'
  }

  if(sf.date %in% gaul.exceptions) sf.type <- 'gaul'
  if(sf.date %in% gadm.exceptions) sf.type <- 'gadm'

  return(list(shpfile_type = sf.type,
              shpfile_date = sf.date))
}
