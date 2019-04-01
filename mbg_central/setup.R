# Some functions to load R packages from specific directories and R functions
# from one or more repo directories
#
## is_singularity() ----------------------------------------------------------->
#' @title Tests if within a Singularity container
#'
#' @description
#' \code{is_singularity} returns TRUE if within a Singularity container, FALSE
#' otherwise.
#'
#' @details
#' This function can determine if it is being run in a Singularity container.
#' If the 'image_dir' argument is left blank, it will return TRUE if it is
#' within any Singularity container. If a path is provided for 'image_dir', this
#' function will test if the container was spun up from an image that lives in
#' the specified directory. This will help us determine if the container was
#' spun up from a specific team's image, which can be important to distinguish
#' where load packages from or if one should attempt to fit a model with a
#' special package like TMB.
#'
#' @param image_dir A path to a directory containing Singularity images
#    [default = ""]. Must be a valid path.
#'
#' @return TRUE/FALSE
#'
#' @family MBG setup functions
#'
#' @seealso This function is used by:
#'   \code{\link{is_lbd_singularity}}
#'
#' @examples
#' # return TRUE if in any Singularity container / FALSE otherwise
#' is_singularity()
#' # TRUE only if container is spun up from image in "<<<< FILE PATH REDACTED >>>>"
#' is_singularity(image_dir = "<<<< FILE PATH REDACTED >>>>")
#' # Will `stop()` because an invalid path was supplied
#' is_singularity(image_dir = "/path/to/nowhere")
#'
is_singularity <- function(image_dir = "") {
  if('SINGULARITY_NAME' %in% names(Sys.getenv())) {
    # if we only want to test if this is an image or not
    if(image_dir == "") return(TRUE)
    # if this is a container and the 'group_image' is not blank, find out if the
    # container was spun up from an image in the specified directory
    else {
      if(!dir.exists(image_dir)) stop(paste0("Invalid file path: '", image_dir, "'\n...Exiting!"))
      image_name <- system('echo $SINGULARITY_NAME', intern = TRUE)
      if(image_name %in% list.files(image_dir)) return(TRUE)
      else return(FALSE)
    }
  }
  else return(FALSE)
}

## is_lbd_singularity() ------------------------------------------------------->
#' @title Tests if within LBD specific Singularity image
#'
#' @description
#' \code{is_lbd_singularity} convenience function to check if a Singularity
#' container was spun up specifically from a LBD image.
#'
#' @details
#' This could be used in multiple places and only wanted to change the path
#' to where the Singularity images live (if it ever changes) in one place and
#' not many.
#'
#' @return TRUE/FALSE
#'
#' @family MBG setup functions
#'
#' @seealso This function is used by:
#' \code{\link{load_R_packages}}
#'
is_lbd_singularity <- function() {
  is_singularity(image_dir = "<<<< FILE PATH REDACTED >>>>")
}

## is_rstudio() --------------------------------------------------------------->
#' @title Tests if within RStudio
#'
#' @description \code{is_rstudio} returns TRUE if within RStudio, FALSE otherwise.
#'
#' @details
#' This function can tell if it is being run from within RStudio. If the
#' 'check_singularity' argument is FALSE (default) it will return TRUE if it is
#' in any RStudio session. If the argument is set to TRUE, it will return TRUE
#' only if 'singularity' is found within 'LD_LIBRARY_PATH'. Unfortunately, we
#' have to do it this way because rstudio-server seems to mask all of the normal
#' system environmental variables, so we can't check for SINGULARITY_NAME in the
#' environment as we do above in \code{is_singularity}.
#'
#' @param check_singularity Logical to check if RStudio is running in a
#'   Singularity container [default = FALSE]
#'
#' @return TRUE/FALSE
#'
#' @family MBG setup functions
#'
#' @seealso This function is used by:
#'   \code{\link{load_R_packages}}
#'
#' @examples
#' # return TRUE if in RStudio / FALSE otherwise
#' is_rstudio()
#' # TRUE only if RStudio lives within a Singularity container
#' is_rstudio(check_singularity = FALSE)
#'
is_rstudio <- function(check_singularity = FALSE) {
  if('RSTUDIO=1' %in% system('env', intern = TRUE)) {
    if(!check_singularity) return(TRUE)
    else {
      if(grepl('singularity', Sys.getenv()['LD_LIBRARY_PATH'])) return(TRUE)
      else return(FALSE)
    }
  } else return(FALSE)
}

## load_R_packages() ---------------------------------------------------------->
#' @title A function to load a list of R packages
#'
#' @description
#' \code{load_R_packages} loads a list of R packages depending on where the
#' script is being run.
#'
#' @details
#' This function is an attempt to address the tangle that is loading packages
#' from different locations. We have different package directories because the
#' packages were build under different environments depending on which cluster
#' node the package is built for, or under a Singularity environment. For
#' example, attempting to load a package build on a prod node will fail to load
#' on a geos node, unless the package is pure R (no shared lib (*.so) to load).
#'
#' There are currently the following scenarios for loading packages
#'
#' Scenario  | How is R Run?           | Load Packages From Directory
#' ----------|----------------------------------------------------------------
#'     A     | RStudio in Singularity  | <<<< FILE PATH REDACTED >>>>/singularity_packages
#'     B     | Cmd-line in Singularity | <<<< FILE PATH REDACTED >>>>/library
#'     C     | Cmd-line on geos node   | <<<< FILE PATH REDACTED >>>>/geos_packages
#'     D     | Cmd-line on prod node   | <<<< FILE PATH REDACTED >>>>/packages
#'
#' This function tries to identify each of these scenarios and then load
#' packages from the correct place. Scenarios C & D are outside of a Singularity
#' container and might be changed to a single option if the OS on the geos nodes
#' is changed to match the prod nodes.
#'
#' Scenario A is currently when running RStudio within a Singularity container
#' spun up from an image *not* built by the LBD team
#' (<<<< FILE PATH REDACTED >>>>/rstudio).
#'
#' Running within RStudio presents some additional difficulty. As discussed
#' above, it is harder to detect if RStudio is being run from within a
#' Singularity container built by a specific team. See details for
#' \code{is_singularity} and \code{is_rstudio}. Once RStudio is moved into an
#' LBD Singularity image, we will likely have to address this.
#'
#' Scenario B is when running R from within a Singularity container spun up from
#' an image that *has been* built by the LBD team. We want to check for this
#' specifically so we know 1) that the required packages are also built into the
#' image and we know where they are located and 2) that special packages
#' (like TMB) are there.
#'
#' Once this function determines where it is being run, it will attempt to load
#' the list of packages from the appropriate directory also setting
#' '.libPaths()' as necessary.
#'
#' Packages are loaded with calls to 'base::library()' so we get errors if
#' packages can't be loaded.
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
  # RStudio running within Singularity container
  if(is_rstudio(check_singularity = TRUE)) {
    dir_name <- paste(R.version$major, R.version$minor, sep = ".")
    # e.g., '<<<< FILE PATH REDACTED >>>>.../singularity_packages/3.5.1/'
    package_lib <- paste('<<<< FILE PATH REDACTED >>>>/singularity_packages', dir_name, sep = "/")
    # Create directory if it doesn't exist (e.g., because R was upgraded)
    if (!file.exists(package_lib)) {
        stop(paste0("You're using a version of R for which no LBD packages have been installed. Please ask the lbd_core team to install packages for ", dir_name))
    }
    message("In RStudio within Singularity container...")
  }
  # Not in RStudio, but in a Singularity container specifically from an LBD image
  else if(is_lbd_singularity()) {
    # All packages are in the '<<<< FILE PATH REDACTED >>>>/R-<version>/library' directory in the
    # LBD Singularity image. In the past, other package directories have snuck in
    # so we enforce this here.
    package_lib <- paste(system('echo $R_HOME', intern = TRUE), 'library', sep = '/')
    if(!dir.exists(package_lib)) {
      stop(paste0("Could not find expected LBD Singularity image package directory: '",
                   package_lib, "'...\nExiting!"))
    }
    message("In LBD Singularity container...")
  # Outside of Singularity container
  } else {
    nodename <- Sys.info()[[4]]
    if(grepl('geos', nodename)) {
      package_lib <- '<<<< FILE PATH REDACTED >>>>/geos_packages'
      message(paste0("On geos node: ", nodename))
    } else {
      package_lib <- '<<<< FILE PATH REDACTED >>>>/packages'
      message(paste0("On prod node: ", nodename))
    }
  }
  # Now set `.libPaths()` and load the packages from the 'package_lib' directory
  message(paste0("Loading packages from location: '", package_lib, "'"))
  .libPaths(package_lib)
  # save the package list from 'lapply' to a variable so it isn't echoed
  packages <- lapply(package_list, library,
                                   lib.loc        = package_lib,
                                   character.only = TRUE)
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
#'                   <<<< FILE PATH REDACTED >>>>/lbd_core/mbg_central/',
#'                   mbg_functions))
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
#' @title Loads all of the other mbg functions from a specified repo (directory)
#'
#' @description
#' \code{load_mbg_functions} finds all of the *.R files with 'functions' in the
#' naming convention and sources them with a message.
#'
#' @details
#' Makes sure that the provided repo directory exists and if so, attempts to
#' load scrips of R functions from within that directory. Note that this
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

## is_lbd_core_repo() --------------------------------------------------------->
#' @title
#' A function to detect if a directory path is the central lbd_core repo or
#' (likely) a fork of that repo.
#'
#' @description
#' \code{is_lbd_core_repo} returns TRUE/FALSE if final subdirectory is
#' 'lbd_core'
#'
#' @details
#' A function to detect if a directory path is the central lbd_core repo or
#' likely a fork of that repo. This is necessary for \code{mbg_setup} to know
#' whether or not to search only the 'mbg_central' subdirectory or not.
#'
#' @param path A path
#'
#' @return TRUE/FALSE
#'
#' @family MBG setup functions
#'
#' @seealso This function is used by:
#' \code{\link{mbg_setup}}
#'
is_lbd_core_repo <- function(path) {
  dir_names <- strsplit(path, '/')[[1]]
  dir_names <- dir_names[which(dir_names != "")]
  if(dir_names[length(dir_names)] == 'lbd_core') return(TRUE)
  else return(FALSE)
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
#' temporary directory to <<<< FILE PATH REDACTED >>>>/geospatial-tempfiles (a location
#' agreed to by IHME infrastructure) and also delete all files a week or older
#' in that directory owned by whomever is running the function.
#'
#' @return NULL
#'
#' @seealso This is called by:
#' \code{\link{mbg_setup}}
#'
fix_raster_tmpdir <- function() {
    require(raster)
    # set temporary file dir. INFR does not regularly delete files from here
    raster::rasterOptions(tmpdir = "<<<< FILE PATH REDACTED >>>>/geospatial-tempfiles")
    # delete files older than 1 week. INFR may take over this duty at some point
    if (!interactive()) raster::removeTmpFiles(h = 24 * 7)
}

## mbg_setup() ---------------------------------------------------------------->
#' @title Loads R packages and sources functions for MBG modeling
#'
#' @description
#' \code{mbg_setup()} Loads R packages and sources functions for MBG modeling
#' with helper functions in 'MBG setup functions' family
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
#' mbg_setup(package_list = 'ggplot2', repos = '<<<< FILE PATH REDACTED >>>>/lbd_core')
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
    stop("Vector of repo paths not characer vector and/or empty...\nExiting!")
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

  fix_raster_tmpdir()
}
