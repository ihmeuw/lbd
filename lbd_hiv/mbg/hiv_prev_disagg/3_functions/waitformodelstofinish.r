#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sleeptime PARAM_DESCRIPTION, Default: 100
#' @param path PARAM_DESCRIPTION, Default: paste0(<<<< FILEPATH REDACTED >>>> "/output/", run_date)
#' @param rd PARAM_DESCRIPTION, Default: run_date
#' @param lv PARAM_DESCRIPTION, Default: loopvars
#' @param showfiles PARAM_DESCRIPTION, Default: TRUE
#' @param showcluster PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname waitformodelstofinish
#' @export
waitformodelstofinish <- function(sleeptime = 100,
                                  path = paste0(<<<< FILEPATH REDACTED >>>> "/output/", run_date),
                                  rd = run_date,
                                  lv = loopvars,
                                  showfiles = TRUE,
                                  showcluster = FALSE,
                                  age=0,
                                  sex=0) {
  n_finished <- length(grep("fin_", list.files(path)))
  
  lv <- data.table(lv)
  names(lv) <- c("reg", "holdout", "age", "sex")
  lv$reg <- as.character(lv$reg)
  
  lv[, file := paste0(path, "/", "fin__bin", age, "_sex", sex, "_", reg, "_", holdout)]
  
  while (n_finished != nrow(lv)) {
    n_finished <- length(grep("fin_", list.files(path)))
    
    message("\n====================================================================================")
    message(sprintf("=====================      Run Date: %s      ======================", rd))
    message(paste0("\nAt ", Sys.time(), " .... ", n_finished, " Models have written output."))
    if (showfiles) {
      message("\nCurrently missing models:")
      for (i in 1:nrow(lv))
        if (file.exists(lv[i, file]) == F) {
          message(paste("Region =", lv[i, 1], "| Holdout =", lv[i, 2], "| Age =", lv[i, 3], "| Sex =", lv[i, 4]))
        }
    }
    n_cluster <- system("qstat | grep job_ | wc | awk '{print $1}'", intern = T)
    message(paste0("\nFuthermore, there are still ", n_cluster, " jobs running on the cluster."))
    if (showcluster) {
      system("qstat -r | grep jobname | grep -oP \"(?<=job_)(.*)\"")
    }
    message("\n====================================================================================")
    message("====================================================================================")
    message("\n")
    Sys.sleep(sleeptime)
  }
  
  unlink(lv$file) # Clean up by deleting extra files once done with the loop
}






#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sleeptime PARAM_DESCRIPTION, Default: 100
#' @param path PARAM_DESCRIPTION, Default: paste0(<<<< FILEPATH REDACTED >>>> "/output/", run_date)
#' @param rd PARAM_DESCRIPTION, Default: run_date
#' @param lv PARAM_DESCRIPTION, Default: loopvars
#' @param showfiles PARAM_DESCRIPTION, Default: TRUE
#' @param showcluster PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname waitformodelstofinish
#' @export
waitformodelstofinish_prediction <- function(sleeptime = 100,
                                             path = paste0(<<<< FILEPATH REDACTED >>>> "/output/", run_date, '/prediction_chunks/'),
                                             rd = run_date,
                                             lv = loopvars,
                                             showfiles = TRUE,
                                             showcluster = FALSE) {
  n_finished <- length(grep("fin_", list.files(path)))
  
  lv <- data.table(lv)
  names(lv) <- c("chunk")
  
  lv[, file := paste0(path, "/", "fin__bin",age,'_',reg,'_',holdout, '_', chunk)]
  
  while (n_finished != nrow(lv)) {
    n_finished <- length(grep("fin_", list.files(path)))
    
    message("\n====================================================================================")
    message(sprintf("=====================      Run Date: %s      ======================", rd))
    message(paste0("\nAt ", Sys.time(), " .... ", n_finished, " prediction chunks have written output."))
    if (showfiles) {
      message("\nCurrently missing chunks:")
      for (i in 1:nrow(lv))
        if (file.exists(lv[i, file]) == F) {
          message(paste("Region =", reg, "| Holdout =", holdout, "| Chunk =", lv[i,1]))
        }
    }
    n_cluster <- system("qstat | grep job_ | wc | awk '{print $1}'", intern = T)
    message(paste0("\nFuthermore, there are still ", n_cluster, " jobs running on the cluster."))
    if (showcluster) {
      system("qstat -r | grep jobname | grep -oP \"(?<=job_)(.*)\"")
    }
    message("\n====================================================================================")
    message("====================================================================================")
    message("\n")
    Sys.sleep(sleeptime)
  }
  
  unlink(lv$file) # Clean up by deleting extra files once done with the loop
}
