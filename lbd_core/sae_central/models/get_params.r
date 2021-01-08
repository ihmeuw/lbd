####################################################################################################
## Description: Define a function for extracting model parameters
##
## Arguments:   model [character] -- file path to the saved model object (output from sdreport(),
##                expected to be named 'out')
##
## Output:      a data.table with four columns:
##              - param: parameter name
##              - part:  part of the model ('fixed' or 'random')
##              - est:   point estimate
##              - se:    standard error
####################################################################################################

get_params <- function(model) {

  # load the model
  load(model)

  # pull out the point estimates and SEs for the fixed effects
  par.fixed <- out$par.fixed
  se.fixed <- sqrt(diag(out$cov.fixed))

  # pull out the point estimates and SEs for the random effects
  par.random <- out$par.random
  se.random <- sqrt(out$diag.cov.random)

  # combine all elements into a data.table and return
  data.table(id = 1:length(c(par.fixed, par.random)),
             param = c(names(par.fixed), names(par.random)),
             part = c(rep("fixed", length(par.fixed)), rep("random", length(par.random))),
             est = c(par.fixed, par.random),
             se = c(se.fixed, se.random))
}
