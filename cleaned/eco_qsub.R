############
## Define qsub function
############

## jobname:      Name for job, must start with a letter
## code:         Filepath to code file
## hold:         Comma-separated list of jobnames to hold the job on
## pass:         List of arguments to pass on to receiving script
## submit:       Should we actually submit this job?
## log:          Should this job create a log in FILEPATH - output and errors?
## proj:         What is the ADDRESS to be used?
## user:         Which FILEPATH are the errors and outputs going to?
## arr_len:      If not NULL, how many array jobs?
## q:            What ADDRESS to get in?
## archive:      Does the job need FILEPATH access?
## m_mem_free:   Maximum memory allotment - required; use G, M, K rather than g, m, k; ie - 10G v 10g
## fthreads:     Maximum thread allotment - required
## max_run_time: Maximum job runtime - err larger than smaller; format is HH:MM:SS

## In code_dir, there must be an r_shell, python_shell, and stata_shell available for use.

## Example submission:
# setwd(code_dir)
# pn <- "nothing"
# jobids <- NULL
# for (cc in c("AFG","USA","VIR")) {
#   if (start<=6 & end >=6) qsub(paste("ADDRESS", cc, sep="_"), paste0(code_dir,"/FILEPATH.r"), hold = paste(pn, collapse=","), pass=list(cc), slots=4, submit=!test)
#   jobids[count] <- paste("ADDRESS", cc, sep="_")
# }
#
# if (start<=6 & end >=6) qsub("ADDRESS", paste0(code_dir,"/FILEPATH.R"),  hold=paste(jobids, collapse=","), submit =!test) 

qsub <- function(jobname, 
                 code, 
                 hold=NULL, 
                 pass=NULL, 
                 submit=F,
                 log=T,
                 intel=F,
                 proj = "ADDRESS",
                 user=NULL, 
                 arr_len=NULL,
                 q = "ADDRESS",
                 archive = F,
                 m_mem_free,
                 fthreads,
                 max_run_time = NULL,
                 error_file = "errors",
                 output_file = "output",
                 priority = "0"
                 ) {
  # user <- Sys.getenv("USER") # Default for linux user grab. "USERNAME" for Windows
  # choose appropriate shell script
  # if(grepl(".r", code, fixed=T) | grepl(".R", code, fixed=T)) shell <- "r_shell.sh" else if(grepl(".py", code, fixed=T)) shell <- "python_shell.sh" else shell <- "stata_shell.sh"
  shell <- "FILEPATH.sh"
  sing_true <-  "-v sing_image=default "
  queue <- paste0("-q ", q)
  # set up jobs to hold for
  if (!is.null(hold)) {
    hold.string <- paste(" -hold_jid \"", hold, "\"", sep="")
  }
  # set up arguments to pass in
  if (!is.null(pass)) {
    pass.string <- ""
    for (ii in pass) pass.string <- paste(pass.string, " \"", ii, "\"", sep="")
  }
  # construct the command
  sub <- paste("qsub",
               paste0(" -l m_mem_free=", m_mem_free, " "),
               paste0(" -l fthread=", fthreads, " "),
               if (!is.null(max_run_time)) paste0(" -l h_rt=", max_run_time, " "),
               if (archive==T) paste0(" -l archive=TRUE "),
               if(log==F) " -e /FILEPATH -o /FILEPATH ",  # don't log (if there will be many log files)
               if(log==T && !is.null(user)) paste0(" -e /FILEPATH/", error_file, " -o /FILEPATH/", output_file, " "),
               if(intel==T) paste0(" -l hosttype=intel "),
               if(proj != "") paste0(" -P ",proj," "),
               if(priority!="0") paste0(" -p ",priority," "),
               sing_true,
               queue,
               if (!is.null(hold)) hold.string,
               " -N ", jobname, " ",
               shell, " ",
               code, " ",
               if (!is.null(pass)) pass.string,
               if (!is.null(arr_len)) paste0(" -t 1:", arr_len, " "),
               sep="")
  # submit the command to the system
  if (submit) {
    system(sub)
  } else {
    cat(paste("\n", sub, "\n\n "))
    flush.console()
  }
}
