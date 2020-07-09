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
#'   [default = ""]. Must be a valid path.
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
#' # TRUE only if container is spun up from image in "<<<< FILEPATH REDACTED >>>>"
#' is_singularity(image_dir = "<<<< FILEPATH REDACTED >>>>")
#' # Will \code{stop()} because an invalid path was supplied
#' is_singularity(image_dir = "/path/to/nowhere")
#'
is_singularity <- function(image_dir = "") {
  if('SINGULARITY_NAME' %in% names(Sys.getenv())) {
    # if we only want to test if this is an image or not
    if(image_dir == "") return(TRUE)
    # if this is a container and the 'group_image' is not blank, find out if the
    # container was spun up from an image in the specified directory
    else {
      if(!dir.exists(image_dir)) stop(paste0("Invalid file path: '", image_dir, "'...\nExiting!"))
      image_name <- Sys.getenv('SINGULARITY_NAME')
      if(file.exists(paste0(image_dir, "/", image_name))) return(TRUE)
      else return(FALSE)
    }
  } else return(FALSE)
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
  is_singularity(image_dir = "<<<< FILEPATH REDACTED >>>>") ||
  is_singularity(image_dir = "<<<< FILEPATH REDACTED >>>>") ||
  is_singularity(image_dir = "<<<< FILEPATH REDACTED >>>>")
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
  if(Sys.getenv("RSTUDIO") == 1) {
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
#'          <<<< FILEPATH REDACTED >>>>
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
#'     A    | R in LBD image           | <<<< FILEPATH REDACTED >>>>
#'     B    | RStudio in LBD image     | <<<< FILEPATH REDACTED >>>>
#'     C    | RStudio in non-LBD image | <<<< FILEPATH REDACTED >>>>
#'

#' This function tries to identify each of these scenarios and then load
#' packages from the correct place. Scenario A and B are when running R from
#' within a Singularity container spun up from an image that *has been* built
#' by the LBD team. We want to check for this specifically so we know 1) that
#' the required packages are also built into the image and we know where they
#' are located and 2) that special packages (like TMB) are there.
#'
#' Scenario C would be for those users that continue to use RStudio from
#' images that LBD does not construct or support
#' (<<<< FILEPATH REDACTED >>>>) and only has a limited number of packages
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
    dir_name <- '<<<< FILEPATH REDACTED >>>>'
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
    # All packages are in the '<<<< FILEPATH REDACTED >>>>' directory in the
    # LBD Singularity image. In the past, other package directories have snuck
    # in so we enforce this here.
    package_lib <- '<<<< FILEPATH REDACTED >>>>'
    if(!dir.exists(package_lib)) {
      stop(paste0("Could not find expected LBD Singularity image package directory: '",
                   package_lib, "'...\nExiting!"))
    }
    message(paste0("In LBD Singularity container from image: '", Sys.getenv("SINGULARITY_NAME"), "'"))
  # Outside of Singularity container
  } else stop("Only LBD Singularity images or non-LBD Singularity RStudio images currently supported...\nExiting!")

  # Now set \code{.libPaths()} and load the packages from the 'package_lib'
  # directory
  message(paste0("Loading packages from location: '", package_lib, "'"))
  .libPaths(package_lib)
  # Load all of the packages with invisible \code{lapply()} call so it isn't
  # echoed
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
#'                   '<<<< FILEPATH REDACTED >>>>',
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
#' @title Loads all of the other MBG functions from a specified repo (directory)
#'
#' @description
#' \code{load_mbg_functions} finds all of the *.R files with 'functions' in the
#' naming convention and sources them with a message.
#'
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
  
  ## Coerce relative to absolute path
  repo <- normalizePath(repo)
  
  message(paste0("\nChecking directory '", repo, "' for R scripts with 'functions' in the filename"))
  if(!dir.exists(repo)) stop(paste0("Directory '", repo, "' does not exist...\nExiting!"))
  functs <- list.files(path       = repo,
                       recursive  = TRUE,
                       pattern    = c('functions.+(R|r)'),
                       full.names = TRUE)
  ## Drop testing functions from being sourced
  functs <- functs[!grepl('testing', functs)]
  functs <- functs[!grepl('LBDCore', functs)]
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
#' temporary directory to <<<< FILEPATH REDACTED >>>> (a location
#' agreed to by IHME infrastructure) and also delete all files a week or older
#' in that directory owned by whomever is running the function.
#'
#' @return NULL
#'
#' @seealso This is called by:
#' \code{\link{mbg_setup}}
#'
fix_raster_tmpdir <- function() {
    # Give the user a message about what is happening
    message("Loading and configuring the raster package...")
    # set temporary file dir. INFR does not regularly delete files from here
    raster_tmp_dir <- paste("<<<< FILEPATH REDACTED >>>>", Sys.getenv("USER"), sep="/")
    if(!dir.exists(raster_tmp_dir)) dir.create(raster_tmp_dir)
    raster::rasterOptions(tmpdir = raster_tmp_dir)
    # delete files older than 25 days (maximum days in geospatial.q)
    if (!interactive()) raster::removeTmpFiles(h = 24 * 25)
}

## load_setthreads() ---------------------------------------------------------->
#' @title Loads a pre-compiled executable for setting MKL/OMP threads
#'
#' @description
#' \code{load_setthreads()} \code{dyn.loads()} the pre-compiled setthreads.so if
#' in an LBD Singularity image.
#'
#' @details
#' This function will first check to see if it is being run in an LBD '
#' Singularity image, and if so, \code{dyn.load()} the pre-compiled shared
#' library "setthreads.so" from the specified location. If it also detects that
#' we are running in RStudio, then it will go ahead and set MKL_NUM_THREADS and
#' OMP_NUM_THREADS to 1 since:
#'   * Since RStudio cannot see environmental variables outside of those set in
#'     the Renviron, we have no way to fire up RStudio in a container with a
#'     specified number threads
#'   * We do not want to bake in serial MKL/OMP execution into the Renviron
#'   * The user can always use \code{setthreads()} to set the number of MKL/OMP
#'     threads to whatever they choose afterward.
#'   * This is the safest thing to do considering issues with \code{mclapply()}
#'     and potentially other functions
#'   * Attempt to set MKL DYNAMIC to FALSE and OMP_NESTED to TRUE as is the 
#'     recommendation for programs with OpenMP parallel regions which may also
#'     call multi-threaded MKL processes. This is just a precaution in case 
#'     the user later decides to enable multi-threading for either/both OpenMP
#'     or MKL. Older LBD images don't have the "setMKLdynamic" or "setOMPnested"
#'     functions built into "setthreads.so", in which case these won't be set.
#' This function won't do anything if it is run outside of an LBD Singularity
#' container.
#'
#' WARNING: The "OMP_NESTED=true" and "MKL_DYNAMIC=false" that are normally set
#'          when an LBD Singularity container is launched through various shell
#'          scripts (e.g. 'singR.sh' and 'shell_sing.sh') may not be set if 
#'          running RStudio from an LBD image that is R3.5.0 or older, or if 
#'          the image is not an LBD image. If these settings are unable to be
#'          made, you may see some unexpected behavior when running programs
#'          with both MKL and OMP parallel regions (like TMB) where both have
#'          been set to multicore. For more recent LBD images, this script will
#'          set those below.
#'
#' @return None
#'
#' @family MBG setup functions
#'
#' @seealso Threads for MKL and OpenMP as well as MKL DYNAMIC and OMP NESTED
#' are set here automatically and can be set to different values using these
#' related functions:
#' \code{\link{setmklthreads()}}
#' \code{\link{setompthreads()}}
#' \code{\link{setmkldynamic()}}
#' \code{\link{setompnested()}}
#'
load_setthreads <- function() {
  # This shared library will only exist in an LBD Singularity image
  if(is_lbd_singularity()) {
    # but let's check anyway
    if(!file.exists('/opt/compiled_code_for_R/setthreads.so')) {
      stop(paste0("'setthreads.so' not found in /opt/compiled_code_for_R...\nExiting!"))
    }
    # Load it and make sure that the two functions for setting threads are
    # available
    dyn.load("/opt/compiled_code_for_R/setthreads.so")
    if(!is.loaded("setMKLthreads")) {
      stop("C function 'setMKLthreads' not loaded...\nExiting!")
    }
    if(!is.loaded("setOMPthreads")) {
      stop("C function 'setOMPthreads' not loaded...\nExiting!")
    }
    # If we are running RStudio in this LBD container, set the threads to 1
    if(is_rstudio()) {
      # Tell the user about it:
      message("Setting MKL/OMP to serial for RStudio")
      invisible(.C("setMKLthreads", as.integer(1)))
      invisible(.C("setOMPthreads", as.integer(1)))
      # The 'setthreads.so' built into early images did not have the
      # 'setMKLdynamic' or 'setOMPnested' functions included. In RStudio, we can
      # get by without them, but to be fully consistent with qsubed jobs, we
      # want to set these if we can.
      # If these functions are available, setting to the recommended default of
      # MKL dynamic = FALSE (0) and OpenMP nested parallelism = TRUE (1)
      if(is.loaded("setMKLdynamic")) invisible(.C("setMKLdynamic", as.integer(0)))
      if(is.loaded("setOMPnested"))  invisible(.C("setOMPnested", as.integer(1)))
    }
  } else message("WARNING: Not and LBD Singularity image. No thread adjustment made.")
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
#' mbg_setup(package_list = 'ggplot2', repos = '<<<< FILEPATH REDACTED >>>>')
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

## set_serial_threads() ------------------------------------------------------->
#' @title Sets MKL/OMP threads to serial
#'
#' @description
#' \code{set_serial_threads()} Uses functions in the "setthreads.so" shared
#' library built into LBD Singularity images to set MKL/OMP threads to serial.
#'
#' @details
#' Uses functions in the "setthreads.so" shared ' library built into LBD
#' Singularity images to set MKL/OMP threads to serial. The shared library
#' should exist in ' the LBD Singularity image and should have already been
#' loaded by \code{mbg_setup()}. This function checks to make sure that it is
#' loaded and if not, attempts to use \code{load_setthreads()} to do so. A
#' warning is generated if this is run outside of an LBD Singularity image and
#' no thread adjustment is done.
#'
#' This can be used to set any multi-threaded operations to serial before a call
#' to a process that forks (like \code{mclapply()}) which is known to have
#' issues with threaded processes. Then \code{set_original_threads()} can be
#' used to set the threads back to their original setting.
#'
#' @return None
#'
#' @family MBG setup functions
#'
#' @seealso This function depends on:
#' \code{\link{load_setthreads}}
#' and is usually coupled with:
#' \code{\link{set_original_threads}}
#'
#' @examples
#' set_serial_threads()
#' mclapply(...)
#' set_original_threads()
#'
set_serial_threads <- function() {
  # This shared library will only exist in a Singularity image
  if(is_lbd_singularity()) {
    # "setthreads.so" should already be loaded, but if it isn't, let's try
    # to load it again
    if(!"setthreads" %in% names(getLoadedDLLs())) load_setthreads()
    # Give the user a message and set the threads to 1
    message("Setting MKL/OMP to serial")
    invisible(.C("setMKLthreads", as.integer(1)))
    invisible(.C("setOMPthreads", as.integer(1)))
  } else message("WARNING: Not and LBD Singularity image. No thread adjustment made.")
}

## set_original_threads() ------------------------------------------------------->
#' @title Sets MKL/OMP threads back to original setting
#'
#' @description
#' \code{set_original_threads()} Uses functions in the "setthreads.so" shared
#' library built into LBD Singularity images to set MKL/OMP back to an original
#' setting.
#'
#' @details
#' Uses functions in the "setthreads.so" shared library built into LBD
#' Singularity images to set MKL/OMP threads back to an original setting. The
#' shared library should exist in the LBD Singularity image and should have
#' already been loaded by \code{mbg_setup()}. This function checks to make sure
#' that it is loaded and if not, attempts to use \code{load_setthreads()} to do
#' so. A warning is generated if this is run outside of an LBD Singularity image
#' and no thread adjustment is done.
#'
#' This function is meant to be used after \code{set_threads_serial()} was used
#' to set MKL/OMP operations to serial before a call to a process that forks
#' (like \code{mclapply()}) which is known to have issues with threaded
#' processes. This function will try and retrieve the original thread settings
#' from the MKL_NUM_THREADS or OMP_NUM_THREADS environmental variable and again
#' set the threads back to those values. If those environmentals are not set,
#' it sets them to 1.
#'
#' @return None
#'
#' @family MBG setup functions
#'
#' @seealso This function depends on:
#' \code{\link{load_setthreads}}
#' and is usually coupled with:
#' \code{\link{set_serial_threads}}
#'
#' @examples
#' set_serial_threads()
#' mclapply(...)
#' set_original_threads()
#'
set_original_threads <- function() {
  # This shared library will only exist in a Singularity image
  if(is_lbd_singularity()) {
    # "setthreads.so" should already be loaded, but if it isn't, let's try
    # to load it again
    if(!"setthreads" %in% names(getLoadedDLLs())) load_setthreads()

    # See if there was an MKL setting based on environmental originally
    if('MKL_NUM_THREADS' %in% names(Sys.getenv())) {
      message("Setting MKL threads back to original setting")
      mkl_threads <- as.integer(Sys.getenv('MKL_NUM_THREADS'))
    } else {
       message("WARNING: No MKL_NUM_THREADS environmental set in LBD image. Setting MKL to serial.")
      mkl_threads <- 1
    }
    invisible(.C("setMKLthreads", as.integer(mkl_threads)))

    # See if there was an OMP setting based on environmental originally
    if('OMP_NUM_THREADS' %in% names(Sys.getenv())) {
      message("Setting OMP threads back to original setting")
      omp_threads <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
    } else {
      message("WARNING: No OMP_NUM_THREADS environmental set in LBD image. Setting OMP to serial.")
      omp_threads <- 1
    }
    invisible(.C("setOMPthreads", as.integer(omp_threads)))
  } else message("WARNING: Not and LBD Singularity image. No thread adjustment made.")
}

## is_integer() ------------------------------------------------------------>
#' @title Tests if object is an integer
#'
#' @description
#' \code{is_integer()} Tests if object is an integer since \code{is.integer()}
#' does not.
#'
#' @details
#' We need to be able to test if an object is an integer or not. Unfortunately
#' R's \code{is.integer()} does not actually tell you if an object is an integer
#' (`is.integer(5)` returns FALSE for example).  This function will return FALSE
#' if something that is not an integer is passed in, including something like a
#' string that R likes to coerce into a NA which when compared with something
#' else resulting in NA instead of TRUE or FALSE, ex.:
#' `x <- 'a'; x == as.integer(x)`
#'
#' @param i An object
#'
#' @return TRUE/FALSE
#'
#' @family MBG setup functions
#'
#' @seealso
#' \code{\link{setmklthreads}}
#' \code{\link{setompthreads}}
#'
#' @examples
#' is_integer(1.2)            # returns FALSE
#' is_integer(1.)             # returns TRUE
#' is_integer(-1)             # returns TRUE
#' is_integer('a')            # returns FALSE
#' is_integer(5)              # returns TRUE
#' is_integer(c(1, 1.2, 3))   # returns TRUE FALSE TRUE
#' is_integer(c('a', 1.2, 3)) # returns FALSE
#'
is_integer <- function(i) {
  # Give a FALSE if R gives a warning or an error. If R tries to coerce it isn't
  # an integer.
  tryCatch(i == as.integer(i),
    warning = function(w) FALSE,
    error   = function(e) FALSE
  )
}

## setmklthreads() ------------------------------------------------------------>
#' @title Sets number of MKL threads
#'
#' @description
#' \code{setmklthreads()} Uses a function in the "setthreads.so" shared library
#' built into LBD Singularity images to set the number of MKL threads to a user
#' defined value
#'
#' @details
#' Uses a function in the "setthreads.so" shared library built into LBD
#' Singularity images to set the number of MKL threads to a user defined value.
#' The shared library should exist in the LBD Singularity image and should have
#' already been loaded by \code{mbg_setup()}. This function checks to make sure
#' that it is loaded and if not, attempts to use \code{load_setthreads()} to do
#' so. A warning is generated if this is run outside of an LBD Singularity image
#' and no thread adjustment is done.
#'
#' @param threads A single value indicating the number of threads to use
#'   [default = 1]
#'
#' @return None
#'
#' @family Mutlti-threading Functions
#'
#' @seealso This function depends on:
#' \code{\link{load_setthreads()}}
#' \code{\link{is_integer()}}
#' Is used by:
#' \code{\link{set_serial_threads()}}
#' \code{\link{set_original_threads()}}
#' And is related to:
#' \code{\link{get_mkl_threads()}}
#' OMP thread setting is done with \code{\link{setompthreads}}
#'
#' @examples
#' setmklthreads(2) # sets the MKL threads to 2
#'
setmklthreads <- function(threads = 1) {
  # This shared library will only exist in a Singularity image
  if(is_lbd_singularity()) {
    # "setthreads.so" should already be loaded, but if it isn't, let's try
    # to load it again
    if(!"setthreads" %in% names(getLoadedDLLs())) load_setthreads()

    # Verify that only an integer has been passed in for the number of threads
    if(!is_integer(threads) | threads <= 0) {
      stop("Integer thread values >= 0 only...\nExiting!")
    } else {
      if(threads > 1) {
        message("WARNING: Enabling multi-threaded MKL...")
        message("         Use caution not to oversubscribe the node or when using")
        message("         with forked processes like 'mclapply()'")
      }
    }
    invisible(.C("setMKLthreads", as.integer(threads)))
  } else message("WARNING: Not and LBD Singularity image. No MKL thread adjustment made.")
}

## setmkldynamic() ------------------------------------------------------------>
#' @title Enables MKL to dynamically change the number of OpenMP threads
#'
#' @description
#' \code{setmkldynamic()} Uses a function in the "setthreads.so" shared library
#' built into LBD Singularity images to enable/disable MKL's ability to change
#' the number of OpenMP threads dynamically.
#'
#' @details
#' Uses a function in the "setthreads.so" shared library built into LBD
#' Singularity images to enable/disable MKL's ability to change the number of
#' OpenMP threads dynamically. Since many of the packages we use in R rely on
#' OpenMP we normally want to disable MKL dynamic, which is done by default here
#' (0 is a FALSE, i.e. requests disabling dynamic adjustment, see: 
#' https://software.intel.com/en-us/mkl-developer-reference-c-mkl-set-dynamic).
#' The shared library should exist in the LBD Singularity image and should have
#' already been loaded by \code{mbg_setup()}. This function checks to make sure
#' that it is loaded and if not, attempts to use \code{load_setthreads()} to do
#' so. A warning is generated if this is run outside of an LBD Singularity image
#' and no MKL dynamic adjustment is done. This function is normally used along
#' with \code{setompnested()} as follows:
#' \code{setmkldynamic(enable = FALSE)} and \code{setompnested(enable = TRUE)}
#' as described here: https://software.intel.com/en-us/articles/recommended-settings-for-calling-intel-mkl-routines-from-multi-threaded-applications
#'
#' @param enable A single logical indicating whether or not to enable MKL
#'   dynamic. [default = FALSE]
#'
#' @return None
#'
#' @family Mutlti-threading Functions
#'
#' @seealso This function depends on:
#' \code{\link{load_setthreads()}}
#' And is related to:
#' \code{\link{setompnested()}}
#' \code{\link{setmklthreads()}}
#' \code{\link{setompthreads()}}
#'
#' @examples
#' setmkldynamic(enable = FALSE) # disables MKL dynamic
#'
setmkldynamic <- function(enable = FALSE) {
  # Only logicals allowed for our only argument
  if(!is.logical(enable)) stop("Logical values only to enable/disable MKL dynamic...\nExiting!")

  # This shared library will only exist in a Singularity image
  if(is_lbd_singularity()) {
    # "setthreads.so" should already be loaded, but if it isn't, let's try
    # to load it again
    if(!"setthreads" %in% names(getLoadedDLLs())) load_setthreads()

    # The 'setthreads.so' built into early images did not have a 'setMKLdynamic'
    # function included. Still can get by without it in some instances, but
    # let's at least let the user know.
    if(!is.loaded("setMKLdynamic")) {
      message(paste0("WARNING: 'setmkldynamic' unavailable in LBD image: '", 
                     Sys.getenv("SINGULARITY_NAME"), "'"))
    } else {
      if(enable) {
        message("WARNING: Enabling MKL dynamic...")
        message("         May cause performance issues when used with OpenMP")
        message("         enabled programs.")
      }
      invisible(.C("setMKLdynamic", as.integer(enable)))
    }
  } else message("WARNING: Not and LBD Singularity image. No MKL dynamic adjustment made.")
}

## setompthreads() ------------------------------------------------------------>
#' @title Sets number of OMP threads
#'
#' @description
#' \code{setompthreads()} Uses a function in the "setthreads.so" shared library
#' built into LBD Singularity images to set the number of OMP threads to a user
#' defined value
#'
#' @details
#' Uses a function in the "setthreads.so" shared library built into LBD
#' Singularity images to set the number of OMP threads to a user defined value.
#' The shared library should exist in the LBD Singularity image and should have
#' already been loaded by \code{mbg_setup()}. This function checks to make sure
#' that it is loaded and if not, attempts to use \code{load_setthreads()} to do
#' so. A warning is generated if this is run outside of an LBD Singularity image
#' and no thread adjustment is done.
#'
#' @param threads A single value indicating the number of threads to use
#'   [default = 1]
#'
#' @return None
#'
#' @family Mutlti-threading Functions
#'
#' @seealso This function depends on:
#' \code{\link{load_setthreads()}}
#' \code{\link{is_integer()}}
#' Is used by:
#' \code{\link{set_serial_threads()}}
#' \code{\link{set_original_threads()}}
#' And is related to:
#' \code{\link{get_omp_threads()}}
#' MKL thread setting is done with \code{\link{setmklthreads}}
#'
#' @examples
#' setompthreads(2) # sets the OMP threads to 2
#'
setompthreads <- function(threads = 1) {
  # This shared library will only exist in a Singularity image
  if(is_lbd_singularity()) {
    # "setthreads.so" should already be loaded, but if it isn't, let's try to
    # load it again
    if(!"setthreads" %in% names(getLoadedDLLs())) load_setthreads()

    # Verify that only an integer has been passed in for the number of threads
    if(!is_integer(threads) | threads <= 0) {
      stop("Integer thread values >= 0 only...\nExiting!")
    } else {
      if(threads > 1) {
        message("WARNING: Enabling multi-threaded OMP...")
        message("         Use caution not to oversubscribe the node or when using")
        message("         with forked processes like 'mclapply()'")
      }
    }
    invisible(.C("setOMPthreads", as.integer(threads)))
  } else message("WARNING: Not and LBD Singularity image. No OMP thread adjustment made.")
}


## setompnested() ------------------------------------------------------------>
#' @title Enables OpenMP nested parallelism
#'
#' @description
#' \code{setompnested()} Uses a function in the "setthreads.so" shared library
#' built into LBD Singularity images to enable/disable OpenMP nested
#' parallelism.
#'
#' @details
#' Uses a function in the "setthreads.so" shared library built into LBD
#' Singularity images to enable/disable OpenMP nested parallelism. Since many of
#' the packages we use in R rely on both OpenMP and MKL, where an OpenMP
#' parallel region can call the MKL and threading is enabled on both, we want to
#' enable OMP nested parallelism, which is done by default here (!= 0 is a TRUE,
#' i.e. requests enabling OMP nested parallelism, see:
#' https://msdn.microsoft.com/en-us/library/sk3zt8e1.aspx)
#' The shared library should exist in the LBD Singularity image and should have
#' already been loaded by \code{mbg_setup()}. This function checks to make sure
#' that it is loaded and if not, attempts to use \code{load_setthreads()} to do
#' so. A warning is generated if this is run outside of an LBD Singularity image
#' and no OMP nesting is done. This function is normally used along with
#' \code{setmkldynamic()} as follows:
#' \code{setmkldynamic(enable = FALSE)} and \code{setompnested(enable = TRUE)}
#' as described here: https://software.intel.com/en-us/articles/recommended-settings-for-calling-intel-mkl-routines-from-multi-threaded-applications
#'
#' @param enable A single logical indicating whether or not to enable OpenMP
#'   nested parallelism [default = TRUE]
#'
#' @return None
#'
#' @family Mutlti-threading Functions
#'
#' @seealso This function depends on:
#' \code{\link{load_setthreads()}}
#' And is related to:
#' \code{\link{setmkldynamic()}}
#' \code{\link{setmklthreads()}}
#' \code{\link{setompthreads()}}
#'
#' @examples
#' setompnested(enable = TRUE) # disables OMP nesting
#'
setompnested <- function(enable = TRUE) {
  # Only logicals allowed for our only argument
  if(!is.logical(enable)) stop("Logical values only to enable/disable OpenMP nested parallelism...\nExiting!")

  # This shared library will only exist in a Singularity image
  if(is_lbd_singularity()) {
    # "setthreads.so" should already be loaded, but if it isn't, let's try
    # to load it again
    if(!"setthreads" %in% names(getLoadedDLLs())) load_setthreads()

    # The 'setthreads.so' built into early images did not have a 'setOMPnested'
    # function included. Still can get by without it in some instances, but
    # let's at least let the user know.
    if(!is.loaded("setOMPnested")) {
      message(paste0("WARNING: 'setompnested' unavailable in LBD image: '", 
                     Sys.getenv("SINGULARITY_NAME"), "'"))
    } else {
      if(!enable) {
        message("WARNING: Disabling OMP nesting...")
        message("         May cause performance issues when used with MKL")
        message("         enabled programs.")
      }
      invisible(.C("setOMPnested", as.integer(enable)))
    }
  } else message("WARNING: Not and LBD Singularity image. No OpenMP nested parallelism adjustment made.")
}

## get_mkl_threads() ---------------------------------------------------------->
#' @title Finds number of threads set for MKL operations
#'
#' @description
#' \code{get_mkl_threads()} Uses environmental variable "MKL_NUM_THREADS" used
#' in LBD Singualirty images to determine how many threads have been assigned
#' for MKL operations.
#'
#' @details
#' The OpenMP function \code{omp_get_num_threads()} will report how many threads
#' have been set for OpenMP operations. Unfortunately, there is no
#' "mkl_get_num_threads()" function in the MKL library, so we have to rely on
#' our MKL_NUM_THREADS environmental variable to find out how many threads have
#' been assigned for MKL operations. Fortunately, we can guarantee that
#' MKL_NUM_THREADS has been set in LBD Singularity containers spun up by either
#' "shell_sing.sh" or "singR.sh".
#'
#' @return number of threads assigned to MKL
#'
#' @family Mutlti-threading Functions
#'
#' @seealso This function is used by:
#' \code{\link{get_total_threads()}}
#'
get_mkl_threads <- function() {
  if(is_lbd_singularity() & !is_rstudio()) {
    if('MKL_NUM_THREADS' %in% names(Sys.getenv())) {
      mkl_threads <- as.integer(Sys.getenv('MKL_NUM_THREADS'))
    } else {
      message("Warning: Unexpectedly, no MKL_NUM_THREADS environmental set in LBD image.")
      mkl_threads <- 1
    }
    return(mkl_threads)
  } else return(1)
}

## get_omp_threads() ---------------------------------------------------------->
#' @title Finds number of threads set for OpenMP operations
#'
#' @description
#' \code{get_omp_threads()} Uses environmental variable "OMP_NUM_THREADS" used
#' in LBD Singualirty images to determine how many threads have been assigned
#' for OpenMP operations.
#'
#' @details
#' The OpenMP function \code{omp_get_num_threads()} will report how many threads
#' have been set for OpenMP operations. Unfortunately, there is no
#' "mkl_get_num_threads()" function in the MKL library, so we have to rely on
#' our OMP_NUM_THREADS environmental variable to find out how many threads have
#' been assigned for OpenMP operations as we do for MKL. Fortunately, we can
#' guarantee that OMP_NUM_THREADS has been set in LBD Singularity containers
#' spun up by either "shell_sing.sh" or "singR.sh".
#'
#' @return number of threads assigned to OpenMP
#'
#' @family Mutlti-threading Functions
#'
#' @seealso This function is used by:
#' \code{\link{get_total_threads()}}
#'
get_omp_threads <- function() {
  if(is_lbd_singularity() & !is_rstudio()) {
    if('OMP_NUM_THREADS' %in% names(Sys.getenv())) {
      omp_threads <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
    } else {
      message("Warning: Unexpectedly, no OMP_NUM_THREADS environmental set in LBD image.")
      omp_threads <- 1
    }
    return(omp_threads)
  } else return(1)
}

## get_total_threads() -------------------------------------------------------->
#' @title Finds total number of threads available
#'
#' @description
#' \code{get_total_threads()} Uses functions \code{get_mkl_threads()} and
#' \code{get_omp_threads()} to determine total cores allocated to the job.
#'
#' @details
#' In an LBD image, the total number of threads are broken up over OpenMP
#' operations (OMP_NUM_THREADS) and MKL operations (MKL_NUM_THREADS). Since we
#' set both OpenMP and MKL to serial for things like \code{mclapply()}, we can
#' at least try and parallelize as much as possible over forked process or other
#' parallel processes. We need to know how many total cores have been allocated
#' to the job to do this.
#'
#' @return total number of threads available for any operation
#'
#' @family Mutlti-threading Functions
#'
#' @seealso This function depends on:
#' \code{\link{get_mkl_threads()}}
#' \code{\link{get_total_threads()}}
#' This function is used by:
#' \code{\link{get_max_forked_threads()}}
#'
#' @examples Set to serial operation, use \code{get_total_threads()} to
#' determine how many total threads are available, run the parallel operation,
#' then set the threads back to the original setting:
#' set_serial_threads()
#' cores <- get_total_threads()
#' model <- fit_glmnet(...)
#' set_original_threads()
#'
get_total_threads <- function() {
  return(get_mkl_threads() * get_omp_threads())
}

## get_max_forked_threads() ----------------------------------------------------->
#' @title Determines how many threads to use in forked applications
#'
#' @description
#' \code{get_max_forked_threads()} Determines how many threads to give forked
#' applications (like \code{mclapply()}).
#'
#' @details
#' When using \code{mclapply()} we should only parallelize over the smaller of
#' the two values: the number of total cores we have available or the number of
#' objects being parallelized over. By first determining this, we can make sure
#' not to use more cores than have been allocated for the job or spawn a bunch
#' of forked \code{mclapply()} processes unnecessarily that then just sleep. In
#' other words, it makes no sense to launch more processes than are being
#' than you have operations to run in parallel, or launch more processes than
#' there are cores available for your job.
#'
#' In order for this to work properly, we assume that OpenMP and MKL operations
#' have been set to serial by \code{set_serial_threads()} before a forked
#' operation is launched, which it should be to avoid hanging the operation.
#'
#' @param threads If NULL, \code{get_total_threads()} is used to determine
#'   how many cores are available for the process. This assumes that all other
#'   multi-core applications have been set to serial with
#'   \code{set_serial_threads()}. An integer value of threads may also be passed
#'   in [default = NULL].
#' @param nobjs Number of objects the forked process will parallelize over
#'
#' @return number of threads to use in forked operation
#'
#'
#' @family Mutlti-threading Functions
#'
#' @seealso This function depends on:
#' \code{\link{get_total_threads()}}
#' Related functions:
#' \code{\link{set_serial_threads()}}
#' \code{\link{set_original_threads()}}
#'
#' @examples Set to serial operation, use \code{get_max_forked_threads()} to
#' determine how many threads to use, run the forked operation, then set the
#' threads back to the original setting:
#' set_serial_threads()
#' cores <- get_max_forked_threads(nobjs = length(folds))
#' model <- mclapply(folds, ...)
#' set_original_threads()
#'
get_max_forked_threads <- function(threads = NULL, nobjs) {
  if(is.null(threads)) threads <- get_total_threads()
  if(!is_integer(threads) | threads <= 0) {
    stop("Integer thread values >= 0 only...\nExiting!")
  }
  return(min(threads, nobjs))
}
