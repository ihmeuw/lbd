####################################################################################################
## Description: Define functions for collapsing draws:
##
##    'collapse_draws':
##      Generic function to calculate point estimates, confidence intervals, and standard errors
##
## Inputs:
##          draws - a data.table of draws of some variable.
##          id_vars - id variables to collapse over, these also become the data.table keys.
##                    By default is set to "area", "year", "sex", and "age.
##
## Outputs: data.table with point estimates, confidence intervals, and standard errors of the given
##          variable for all specified id variables.
####################################################################################################


collapse_draws <- function(draws, id_vars=c("level", "area", "year", "sex", "age")) {
  est <- draws[, as.list(c(mean(value),
                           quantile(value, c(0.025, 0.975), type=5),
                           sd(value))),
               by=id_vars]

  setnames(est, c(id_vars, c("mean", "lb", "ub", "se")))
  setkeyv(est, id_vars)
  est
}
