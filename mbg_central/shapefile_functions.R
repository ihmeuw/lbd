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
                                              in_dir = '<<<< FILE PATH REDACTED >>>>',
                                              out_dir = '<<<< FILE PATH REDACTED >>>>') {
  
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
      results <- mclapply(files_to_synchronize, factory(save_shapefile_as_rds), 
                          in_dir = in_dir, out_dir = out_dir, verbose = verbose,
                          mc.cores = cores)
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
                                  in_dir = '<<<< FILE PATH REDACTED >>>>',
                                  out_dir = '<<<< FILE PATH REDACTED >>>>') {

  # Load one function from stringr
  str_match <- stringr::str_match

  shape <- str_match(shapefile_path, "//(.*).shp$")[2]

  lock_file(shape)
  
  if(verbose == T) {
    message(paste0("     ", shape))
  }

  the_shape <- try(readOGR(dsn = in_dir, layer = shape), silent = T)

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
                                fast_shapefile_dir = '<<<< FILE PATH REDACTED >>>>') {

  wait_for_lock(shape)
  return(readRDS(paste0(fast_shapefile_dir, shape, ".rds")))
}

#####################################################################
## Functions to implement a locking system
#####################################################################

# Lock a shapefile (e.g. while editing / updating)
lock_file <- function(shp,
                      fast_shapefile_dir = '<<<< FILE PATH REDACTED >>>>') {
  
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
                          fast_shapefile_dir = '<<<< FILE PATH REDACTED >>>>',
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
                           fast_shapefile_dir = '<<<< FILE PATH REDACTED >>>>') {

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
                        fast_shapefile_dir = '<<<< FILE PATH REDACTED >>>>') {

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