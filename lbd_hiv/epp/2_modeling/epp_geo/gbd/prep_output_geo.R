################################################################################
## Purpose: This is the set of functions that convert the model that is output 
##          by the IMIS process to create the full set of epp outputs.  Since 
##          IMIS runs epp hundreds of thousands of times to create the fit, it 
##          only evaluates the prevalence output of a given epp run to assing 
##          a likelihood to the force of infection that created that prevalence
##          output.  These fucntions here are used to take the posterior 
##          distribution of those force of infection vectors and use them to 
##          simulate the posterior distributions of the other outputs such as
##          incidence, prevalence, and art coverage.
################################################################################

#' @title new.inf
#' @description helper function that extract the new infections from an EPP model post fitting.
#'
#' @param mod the model fit being used for the simulation of new infections.
#' @param fp the model object that holds the parameters used for the new infections, in this case explicitly the spline parameters for the force of infection
#' the proj_steps, and the relative infectivity of someone on ART
#'
#' @return none explicitly, this operates in the simmod function
#'
#' @export

new.inf <- function(mod, fp) {
  attr(mod, "rvec")[fp$proj.steps %% 1 == 0.5] * (rowSums(mod[,-1,1]) + fp$relinfectART * rowSums(mod[,-1,-1])) / rowSums(mod) * mod[,1,1]
}


#' @title suscept.pop
#' @description helper function that extracts the susceptible pop from an EPP model post fitting.
#'
#' @param mod the model fit being used for the simulation of susceptible population.
#'
#' @return none explicitly, this operates in the simmod function
#'
#' @export
suscept.pop <- function(mod) {
  mod[,1,1]
}

#' @title plwh
#' @description helper function that extracts the people living with hiv from an EPP model post fitting.
#'
#' @param mod the model fit being used for the simulation of people living with hiv.
#'
#' @return none explicitly, this operates in the simmod function
#'
#' @export
plwh <- function(mod) {
  rowSums(mod[,-1,])
}


#' @title art
#' @description helper function that extracts the number of art patients from an EPP model post fitting.
#'
#' @param mod the model fit being used for the simulation of number of art patients.
#'
#' @return none explicitly, this operates in the simmod function
#'
#' @export
art <- function(mod) {
  rowSums(mod[,-1,-1])
}


#' @title total.pop
#' @description helper function that extracts the total population from an EPP model post fitting.
#'
#' @param mod the model fit being used for the simulation of a population.
#'
#' @return none explicitly, this operates in the simmod function
#'
#' @export
total.pop <- function(mod) {
  rowSums(mod)
}



#' @title nat.draws
#' @description helper function extracts the national draws from an EPP fit.
#'
#' @param result NEED TO FIGURE THIS OUT.
#'
#' @return FIGURE THIS OUT!!!!!!
#'
#' @export
nat.draws <- function(result) {
  nat.inf <- Reduce('+', lapply(result, function(x){x$new.inf}))
  nat.suscept <- Reduce('+', lapply(result, function(x){x$suscept.pop}))

  nat.incid <- nat.inf/nat.suscept

  nat.plwh <- Reduce('+', lapply(result, function(x){x$plwh}))
  nat.total.pop <- Reduce('+', lapply(result, function(x){x$total.pop}))

  nat.art <- Reduce('+', lapply(result, function(x){x$art}))

  nat.art.cov <- nat.art / nat.plwh
  nat.prev <- nat.plwh/nat.total.pop
  output <- list(prev = nat.prev, incid = nat.incid, art = nat.art.cov, art_num = nat.art, pop = nat.total.pop)
  return(output)
}

#' @title simfit.gbd
#' @description this is the main function that takes in an EPP fit and turns it into meaningful outputs like prevalence and incidence.
#'
#' @param fit the output of the IMIS functions that do the actual fitting process
#' @param random_walk a boolian indicating if the random walk specification for r was used in fitting or not.
#' @param no.anc a boolian indicating if there was anc data in the fitting process or not.
#'
#' @return outputs a new fit object with meaningful disease quantities like prevalence, incidence, and mortality
#'
#' @export

# Simulate fit and output incidence and prevalence
simfit.gbd <- function(fit, random_walk = F, no.anc = no.anc){
  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp, no.anc = no.anc))
  # Random-walk projection method
  if (random_walk) {
    if (exists("eppmod", where = fit$fp) && fit$fp$eppmod == "rtrend")
      stop("Random-walk projection is only used with r-spline model")

    fit$rvec.spline <- sapply(fit$param, "[[", "rvec")
    firstidx <- (fit$likdat$firstdata.idx - 1) / (fit$fp$dt + 1) # which(fit$fp$proj.steps == fit$fp$tsEpidemicStart) ## assume SD only depends on rvec of years with ANC/survey data
    lastidx <- (fit$likdat$lastdata.idx - 1) / (fit$fp$dt + 1)

    firstidx.mean <- min(which(fit$fp$artnum.ts > 0)) # which(fit$fp$proj.steps == fit$fp$tsEpidemicStart) 

    ### Get the year of prev data that starts to decrease in the recent years
    hhs.dt <- fit$likdat$hhslik.dat
    if (nrow(fit$likdat$hhslik.dat) != 0) {
      hhs.dt <- hhs.dt[order(hhs.dt$year),]
      hhs.dt$dif <- c(0, diff( hhs.dt$prev))
      proj.year <- hhs.dt[hhs.dt$dif <= 0,]$year[nrow(hhs.dt[hhs.dt$dif <= 0,]) - 1]
      proj.year.idx <- ifelse(length(proj.year) > 0, which(fit$fp$proj.steps == proj.year), 0)
    } else {proj.year.idx <- 0}
    ## replace rvec with random-walk simulated rvec
    fit$param <- lapply(fit$param, function(par){
      ### compare the prev year with the lowest rvec year  
      rvec.min.idx <- which(par$rvec == min(par$rvec))
      first_projidx <- ifelse(proj.year.idx > rvec.min.idx, proj.year.idx, rvec.min.idx)

      par$rvec <- sim_rvec_rwproj(par$rvec, firstidx, lastidx, first_projidx, firstidx.mean, fit$fp$dt); par}) 
  }
  
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list = par))
  mod.list <- lapply(fp.list, simmod)
  
  if ("mvec" %in% names(attributes(mod.list[[1]]))) {fit$mvec <- sapply(mod.list, attr, "mvec")}
  if ("svec" %in% names(attributes(mod.list[[1]]))) {fit$svec <- sapply(mod.list, attr, "svec")}
  fit$rvec <- sapply(mod.list, attr, "rvec")
  fit$prev <- sapply(mod.list, prev)
  fit$incid <- mapply(incid, mod = mod.list, fp = fp.list)
  fit$popsize <- sapply(mod.list, rowSums)
  # fit$pregprev <- mapply(epp::fnPregPrev, mod.list, fp.list)

  ## GBD additions
  fit$new.inf <- mapply(new.inf, mod = mod.list, fp = fp.list)
  fit$suscept.pop <- sapply(mod.list, suscept.pop)
  fit$plwh <- sapply(mod.list, plwh)
  fit$total.pop <- sapply(mod.list, total.pop)
  fit$art <- sapply(mod.list, art)
  return(fit)
}
