# Some functions to load R packages from specific directories and R functions
# from one or more repo directories

## load_R_packages() ---------------------------------------------------------->
#' @title A function to load a list of R packages
#'
#' @description
#' \code{load_R_packages} loads a list of R packages depending on where the
#' script is being run.
#'
#' @details
#' WARNING: The logic for this function is dependent on the SINGULARITY_NAME
#'          environmental being written into the Renviron in the Singularity
#'          image so that it can be detected from within RStudio. This should
#'          match the name of the image and will therefore be an exact match for
#'          the SINGULARITY_NAME environmental variable created by Singularity
#'          for consistency. Unfortunately, this variable is not accessible from
#'          RStudio, which is why we go to the trouble of hard coding it into
#'          the Renviron.
#'
#' WARNING: In addition to the SINGULARITY_NAME this function also assumes that
#'          all packages built into LBD Singularity images live at:
#'          /usr/local/R-<version>/library
#'
#' This function is an attempt to address the tangle that is loading packages
#' from different locations. This was much worse in the past when we were
#' loading packages for different versions of R installed on different cluster
#' nodes running different OS's as well as RStudio and R both under Singularity
#' images. After the team moved to R and RStudio only running within Singularity
#' containers, the process of determining where to load packages from became
#' much simpler.
#'
#' The following Scenario's are currently supported by this function:
#'
#' Scenario | How is R Run?            | Load Packages From Directory
#' ---------|------------------------------------------------------------------
#'     A    | R in LBD image           | '<<<< FILEPATH REDACTED >>>>'
#'     B    | RStudio in LBD image     | '<<<< FILEPATH REDACTED >>>>'
#'     C    | RStudio in non-LBD image | '<<<< FILEPATH REDACTED >>>>'
#'
#' This function tries to identify each of these scenarios and then load
#' packages from the correct place. Scenario A and B are when running R from
#' within a Singularity container spun ' up from an image that *has been* built
#' by the LBD team. We want to check for this specifically so we know 1) that
#' the required packages are also built into the image and we know where they
#' are located and 2) that special ' packages (like TMB) are there.
#'
#' Scenario C would be for those users that continue to use RStudio from
#' images that LBD does not construct or support
#' ('<<<< FILEPATH REDACTED >>>>') and only has a limited number of packages
#' built in to the image, which is why we rely on the directory of R packages
#' outside of the image.
#'
#' Once this function determines where it is being run, it will attempt to load
#' the list of packages from the appropriate directory also setting
#' '.libPaths()' as necessary.
#'
#' Packages are loaded with calls to 'base::library()' so we get errors if
#' packages cannot be loaded.
#'
#' @param package_list Character vector defining one or more package names.
#'
#' @return None
#'
#' @family MBG setup functions
#'
#' @seealso This function is used by:
#' \code{\link{mbg_setup}}
#'
load_R_packages <- function(package_list) {
  # This central scripts are intended to be run on the cluster, and there are
  # have been issues with packages being overwritten by Windows in the past,
  # so bail if we aren't on a Linux machine as a very basic check. Could also
  # have this check as part of 'mbg_setup()', but specifically loading packages
  # under Windows has really been the issue in the past we want to test for here
  # and we may also want to use this function stand alone outside of
  # 'mbg_setup()'
  if(Sys.info()[1] != 'Linux') {
    stop("MBG central scripts should be run on the IHME cluster only...\nExiting!")
  }
  # Not an LBD RStudio running within Singularity container
  if(is_rstudio(check_singularity = TRUE) & !is_lbd_singularity()) {
    dir_name <- paste(R.version$major, R.version$minor, sep = ".")
    package_lib <- paste('<<<< FILEPATH REDACTED >>>>', dir_name, sep = "/")
    # Create directory if it doesn't exist (e.g., because R was upgraded)
    if (!dir.exists(package_lib)) {
        stop(paste0("You're using a version of R for which no LBD packages have been installed. Please ask the lbd_core team to install packages for ", dir_name))
    }
    message("In RStudio within non-LBD Singularity container...")
  # A Singularity container specifically from an LBD image
  # If the SINGULARITY_NAME environmental has been setup properly (as above) it
  # doesn't matter if this is RStudio or not. LBD Singularity images always have
  # their packages installed at the same location.
  } else if(is_lbd_singularity()) {
    # All packages are in the '/usr/local/R-<version>/library' directory in the
    # LBD Singularity image. In the past, other package directories have snuck
    # in so we enforce this here.
    package_lib <- paste(R.home(), 'library', sep = '/')
    if(!dir.exists(package_lib)) {
      stop(paste0("Could not find expected LBD Singularity image package directory: '",
                   package_lib, "'...\nExiting!"))
    }
    message("In LBD Singularity container...")
  # Outside of Singularity container
  } else stop("Only LBD Singularity images or non-LBD Singularity RStudio images currently supported...\nExiting!")

  # Now set `.libPaths()` and load the packages from the 'package_lib' directory
  message(paste0("Loading packages from location: '", package_lib, "'"))
  .libPaths(package_lib)
  # Load all of the packages with invisible `lapply()` call so it isn't echoed
  invisible(lapply(package_list, library,
                                 lib.loc        = package_lib,
                                 character.only = TRUE))
}

## source_functions() --------------------------------------------------------------->
#' @title Sources a list of functions
#'
#' @description
#' \code{source_functions} sources a list of functions with
#' with complete path names
#'
#' @details
#' A standalone function that will source in functions from a provided character
#' vector of functions names with complete paths.
#'
#' @param functions Character vector of functions names including COMPLETE path
#'
#' @return None
#'
#' @family MBG setup functions
#'
#' @seealso This function is used by:
#'   \code{\link{load_mbg_functions}}
#'
#' @examples
#' mbg_functions <- c('setup.R','mbg_functions.R', 'misc_functions.R')
#' source_functions(paste0(
#'                   '<<<< FILEPATH REDACTED >>>>'))
#'
source_functions <- function(functions) {
  if(length(functions) == 0) {
    stop(paste0("'functions' argument empty...\nExiting!"))
  } else {
    for(funct in functions) {
      message(paste0("--> Loading: '", funct, "'"))
      source(funct)
    }
  }
}

## load_mbg_functions() ------------------------------------------------------->
#' @title Loads all of the other MBG functions from a specified repo (directory)
#'
#' @description
#' \code{load_mbg_functions} finds all of the *.R files with 'functions' in the
#' naming convention and sources them with a message.
#' @details
#' Makes sure that the provided repo directory exists and if so, attempts to
#' load scripts of R functions from within that directory. Note that this
#' function will *only* load R scripts with 'functions' in the naming convention
#' of the R script file name it finds after searching recursively from the
#' directory provided for 'repo'.
#'
#' @param repo Path to repository directory to attempt to load functions from
#'
#' @return None
#'
#' @family MBG setup functions
#'
#' @seealso This function is used by:
#' \code{\link{mbg_setup}}
#'
load_mbg_functions <- function(repo) {
  message(paste0("\nChecking directory '", repo, "' for R scripts with 'functions' in the filename"))
  if(!dir.exists(repo)) stop(paste0("Directory '", repo, "' does not exist...\nExiting!"))
  functs <- list.files(path       = repo,
                       recursive  = TRUE,
                       pattern    = 'functions',
                       full.names = TRUE)
  if(length(functs) == 0) {
    stop(paste0("Found no R scripts in '", repo, "' with 'functions' in the filename...\nExiting!"))
  } else source_functions(functs)
}


## fix_raster_tmpdir() -------------------------------------------------------->
#' @title
#' Set temporary directory used by the raster package and clear old files.
#'
#' @description
#' \code{fix_raster_tmpdir()} loads and configures the raster package.
#'
#' @details
#' By default the raster package uses /tmp to store temporary files by default.
#' This is problematic as IHME machines are not configured to have a large
#' amount of /tmp space, and multiple users will quickly fill the directory
#' leading to a non-functioning computer. This function does two things: set the
#' temporary directory to '<<<< FILEPATH REDACTED >>>>'and also delete all files a week or older
#' in that directory owned by whomever is running the function.
#'
#' @return NULL
#'
#' @seealso This is called by:
#' \code{\link{mbg_setup}}
#'
fix_raster_tmpdir <- function() {
    library(raster)
    # Give the user a message about what is happening
    message("Loading and configuring the raster package...")
    # set temporary file dir.
    raster::rasterOptions(tmpdir = '<<<< FILEPATH REDACTED >>>>')
    # delete files older than 1 week.
    if (!interactive()) raster::removeTmpFiles(h = 24 * 7)
}

## mbg_setup() ---------------------------------------------------------------->
#' @title Loads R packages and sources functions for MBG modeling
#'
#' @description
#' \code{mbg_setup()} Loads R packages and sources functions for MBG modeling
#' with helper functions in 'MBG setup functions' family
#'
#'
#' @details
#' This function makes use of some other helper functions in the 'MBG setup
#' functions' family to load all the necessary packages in a provided vector
#' as well as functions for the core_repo and optionally additional repos such
#' as indicator repos and/or development repos. Functions will be read in from
#' repos in the order they appear in the 'repos' vector argument. The examples
#' below show how this could be run with one or multiple repos.
#'
#' @param package_list A vector of R package names to attempt to load
#' @param repos A vector of paths to repository directories. At a
#'   minimum, this should contain the 'core_repo', but may also contain an
#'   indicator repo and one or more additional repos. Functions will be loaded
#'   from each repo in the order they appear in this vector. If any path has
#'   'lbd_core' as the last directory in the path, only it's subdirectory
#'   'mbg_central' will be recursively searched for functions to load as is the
#'   custom. A single, complete path (or object storing a single path) can also
#'   be provided as a string.
#'
#' @return None
#'
#' @family MBG setup functions
#'
#' @seealso This function uses:
#' \code{\link{load_R_packages}}
#' \code{\link{load_mbg_functions}}
#' \code{\link{is_lbd_core_repo}}
#'
#' @examples
#' mbg_setup(package_list = 'ggplot2', repos = '<<<< FILEPATH REDACTED >>>>/lbd_core')
#' mbg_setup(package_list = package_list, repos = c(core_repo, indic_repo))
#'
mbg_setup <- function(package_list, repos) {
  # check 'package_list' and 'repos' list to make sure they are of correct type
  # and actually contain something
  if(!is.character(package_list) | length(package_list) == 0) {
    stop("Package list not characer vector and/or empty...\nExiting!")
  }
  # Similarly, check repos list and source functions
  if(!is.character(repos) | length(repos) == 0) {
    stop("Vector of repo paths not character vector and/or empty...\nExiting!")
  }

  # load packages
  load_R_packages(package_list)

  # source in functions
  for(repo in repos) {
    # if you are pointing to the lbd_core repo (or a fork) only search the
    # '/mbg_central' subdirectory for functions to source
    if(is_lbd_core_repo(repo)) repo <- paste0(repo, '/mbg_central')
    # load the functions from the repo
    load_mbg_functions(repo)
  }

  # Loading and configuring the raster package
  fix_raster_tmpdir()

  # If we are in an LBD Singularity image, load "setthreads.so" and set to
  # serial if we are running RStudio
  load_setthreads()
}
