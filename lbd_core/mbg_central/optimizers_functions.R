#' @title Run BRT Optimizer in Python
#'
#' @description Tune hyperparameters for BRT and save tuning results with shrinkage, intersection depth, and the number of trees.
#'
#' @param funcs.file_path directory containing the optimizer script
#'
#' @param funcs.file file name of optimizer script file
#'
#' @param bounds.file_path directory containing the search space boundaries csv
#'
#' @param bounds.file file name of search space boundaries csv file
#'
#' @param data.file_path model data file path.
#'
#' @param data.loc model data file sub dir -- i.e. -- data.loc <- paste0('/data_train_', jobnum, '.csv'),
#'
#' @param optimizer choice of non-parametric model to optimize with, either: 'rf' (Random Forest), 'gp' (Gaussian Process), 'brt' (Boosted Regression Tree)
#'
#' @param learner choice of BRT model, either: brtR, brtC, xgbR, xgbC -- R=regression, C=classification, xgb=use XGBoost BRT, BRT=use sklearn brt
#'
#' @param cv_folds number of cross-validation folds / set of parameters
#'
#' @param n_calls number of draws to make to inform posterior distribution over parameter space
#'
#' @param jobnum The resulting file is saved to: data.file_path + "best_pars_" + jobnum + ".csv"
#'
#' @param col_start numeric index. Assumes passed array has covariate data from here until the last column
#'
#' @param return_cmd_instead_of_running if TRUE returns the string to run the program instead of calling system(cmd). You probably don't want this unless you're testing or want to embed this command in a qsub call.
#' @return The string of the qsub command, or nothing if the job is chosen to be submitted
#'
#' @export
run_optimizerPy <- function(funcs.file_path,
                            funcs.file,
                            bounds.file_path,
                            bounds.file,
                            data.file_path,
                            data.loc,
                            optimizer,
                            learner,
                            cv_folds,
                            n_calls,
                            jobnum,
                            col_start,
                            return_cmd_instead_of_running = FALSE) {
  func_path <- paste0(funcs.file_path, funcs.file)
  bounds_path <- paste0(bounds.file_path, bounds.file)
  data_path <- paste0(data.file_path, data.loc)

  if (!file.exists(func_path)) stop(paste("optimizer file", func_path, "does not exist. Check the path"))
  if (!file.exists(bounds_path)) stop(paste("bounds file", bounds_path, "does not exist. Check the path"))
  if (!file.exists(data_path)) stop(paste("model data file", data_path, "does not exist. Check the path"))

  # optimizer should be in "gp", "rf" or "brt"

  cmd <- paste0("FILEPATH")
  if (return_cmd_instead_of_running) {
    return(cmd)
  } else {
    system(cmd)
  }
}
