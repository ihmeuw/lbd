####################################################################################################
## Description: Define a function for submitting jobs to the cluster
##
## Inputs:      code [character] -- filepath to the code to run.
##              name [character] -- the job name. If NULL, an argument will be constructed from
##                the name of the code and the arguments.
##              arguments [vector] -- a vector of arguments to pass.
##              hold [numeric vector] -- a vector of job IDs to hold on.
##              slots [numeric] -- the number of slots to assign to this job.
##              intel [logical] -- should the job be submitted specifically to an intel cluster?
##              proj [character] -- the cluster project to run.
##              sgeoutput [character] -- file path where sgeoutput files should be written (if NULL
##                then no sgeoutput files are written).
##              skip_if_exists -- file path(s) to check for output file. If the output file(s) are
##                found no job will be submitted.
##              singularity_opts -- pass in a named list of environmental variables that the shell_sing.sh script knows about: 'SET_OMP_THREADS'
##                and 'SET_MKL_THREADS'. Passing in other environmental names in the list will result in an error. If left as 'NULL' and the
##                shell_sing.sh script will use the default setting of SET_OMP_THREADS = 1 and SET_MKL_THREADS = 1. For example SET_OMP_THREADS=1
##                and SET_MKL_THREADS=4 can be achieved by passing in singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=4)
##
## Outputs:     the job ID assigned to the submitted job
####################################################################################################

qsub <- function(code, name = NULL, arguments = NULL, hold = NULL,
                 slots = 1, intel = F, proj = NULL, sgeoutput = NULL,
                 skip_if_exists = NULL, singularity_opts = NULL) {

  # Define shell, use default singularity image as a default
  corerepo  <- "<<<< FILEPATH REDACTED >>>>"
  shell <- "<<<< FILEPATH REDACTED >>>>"

  # Set singularity opts, make sure they are defined correctly
  if(is.null(singularity_opts)) singularity_opts <- list(SET_OMP_THREADS = 1, SET_MKL_THREADS = 1)
  if(!is.list(singularity_opts)) {
    message("Expected a list for 'singularity_opts' argument")
    stop("Exiting!")
  } else if (!all(names(singularity_opts) %in% c("SET_OMP_THREADS", "SET_MKL_THREADS"))) {
    message("SET_OMP_THREADS and SET_MKL_THREADS were unexpected in the singularity_opts argument, see documentation")
    stop("Exiting!")
  }
  sing_environment <- paste("-v sing_image=default", paste0(" -v ", names(singularity_opts), "=", singularity_opts, collapse = ""))

  # check for skip_if_exists file (if not NULL)
  if (!is.null(skip_if_exists)) {
    if (sum(!file.exists(skip_if_exists)) == 0) return(1) # '1' here is a placeholder and signifies no job submitted
  }

  # format job name, if you don't have a custom name
  if (is.null(name)){
    name <- gsub(".r$", "", tail(strsplit(code, "/")[[1]], 1))
    name <- paste0(c(name, arguments[-1]), collapse="_") # skip first argument as this is always main_dir
  }

  # format arguments and hold IDs
  if(!is.null(arguments)) arguments <- paste(arguments, collapse=" ")
  if(!is.null(hold) & length(hold)>1) hold <- paste(hold, collapse=",")

  # make sure the sgeoutput folder exists
  if (!is.null(sgeoutput)) {
    if (!file.exists(sgeoutput)) {
      sgeoutput <- NULL
    }
  }

  # make error and output directories
  if(!is.null(sgeoutput)) {
    if (!dir.exists(paste0(sgeoutput, "/output"))) dir.create(paste0(sgeoutput, "/output"))
    if (!dir.exists(paste0(sgeoutput, "/errors"))) dir.create(paste0(sgeoutput, "/errors"))
  }

  # check if we're on the dev cluster, and if so set intel to F (note: prod = cn200-cn453; dev = cn454-cn519)
  if (substr(Sys.info()["nodename"], 3, 5) >= 454) intel <- F

  # construct and submit qsub command and return the job ID
  x <- paste("qsub -cwd",
             if (!is.null(sgeoutput)) paste0("-o ", sgeoutput, "/output/ -e ", sgeoutput, "/errors/"),
             if (slots > 1) paste("-pe multi_slot", slots),
             if (!is.null(proj)) paste("-P", proj),
             if (!is.null(proj)) if (proj == "proj_geo_nodes") "-l geos_node=TRUE",
             "-N", name,
             sing_environment,
             if (!is.null(hold)) paste("-hold_jid", hold),
             if (intel) "-l hosttype=intel",
             shell, code,
             if (!is.null(arguments)) arguments)
  id <- system(x, intern = T)
  return(as.numeric(strsplit(id, " ")[[1]][3]))
}
