#' @title Read INLA prior for TMB
#'
#' @description Read in a prior specification that is suited for INLA and make 
#'   it TMB readable.
#'
#' @param prior_string character, character vec of length 1 specifying priors
#'
#' @return List specifying a TMB prior, containing three elements:
#'   - logNormal: Is the prior lognormal (1) or not (0)?
#'   - par1: The first shape parameter. In the lognormal case, the mean
#'   - par2: The second shape parameter. In the lognormal case, the variance
#'
read_inla_prior <- function(prior_string){
  prior_list <- eval(parse(text=prior_string[1]))
  if(!(prior_list$prior %in% c("normal", "loggamma"))){
    stop("TMB implementation only supports normal and loggamma priors.")
  }
  return(list(
    logNormal = ifelse(prior_list$prior == "normal", 1, 0),
    par1 = prior_list$param[1],
    par2 = prior_list$param[2]
  ))
}