#################
####  Prior  ####
#################

ldinvgamma <- function(x, alpha, beta){
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  return(log.density)
}

## r-hybrid prior parameters
rlog_pr_mean <- c(log(0.35), log(0.09), log(0.2), 1993)
rlog_pr_sd <- c(0.5, 0.3, 0.5, 5)

rw_prior_shape <- 300
rw_prior_rate <- 1.0
rw_prior_sd <- 0.06

rlogistic <- function(t, p){
  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}

lprior_iota <- function(par, fp){
  
  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }
  
  
  if(exists("logitiota", fp) && fp$logitiota)
    ldinvlogit(par)  # Note: parameter is defined on range logiota.unif.prior, so no need to check bound
  else
    dunif(par, logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE)
}

sample_iota <- function(n, fp){
  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }
  if(exists("logitiota", fp) && fp$logitiota)
    return(logit(runif(n)))
  else
    runif(n, logiota.unif.prior[1], logiota.unif.prior[2])
}






## r-spline prior parameters
expand.iota <- TRUE
if(expand.iota) {
  logiota.unif.prior <- c(log(1e-14), log(0.000025))  
} else {
  logiota.unif.prior <- c(log(1e-14), log(0.0025)) 
}
tau2.prior.rate <- 0.5
invGammaParameter <- 0.001   #Inverse gamma parameter for tau^2 prior for spline
muSS <- 1/11.5               #1/duration for r steady state prior

## r-trend prior parameters
t0.unif.prior <- c(1970, 1990)
t1.unif.prior <- c(10, 30)
logr0.unif.prior <- c(1/11.5, 10)
rtrend.beta.pr.sd <- 0.2

vinfl.prior.rate <- 1/0.015

lprior <- function(theta, fp, no.anc = TRUE){
 if (no.anc == TRUE){
   if(!exists("eppmod", where = fp)) {fp$eppmod <- "rspline"}  # backward compatibility
   
   if(fp$eppmod == "rhybrid"){
     epp_nparam <- fp$rt$n_param+1
     lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE)) +
       sum(dnorm(theta[4+1:fp$rt$n_rw], 0, rw_prior_sd, log=TRUE))
     lpr <- lpr + lprior_iota(theta[fp$rt$n_param+1], fp)
     return(lpr)
   } else if(fp$eppmod == "rspline"){
     nk <- fp$numKnots
     tau2 <- exp(theta[nk+2])
     
     return(sum(dnorm(theta[3:nk], 0, sqrt(tau2), log=TRUE)) +
              dunif(theta[nk+1], logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE) + 
              #dnorm(theta[nk+2], ancbias.pr.mean, ancbias.pr.sd, log=TRUE) +
              ldinvgamma(tau2, invGammaParameter, invGammaParameter) + log(tau2) +  # + log(tau2): multiply likelihood by jacobian of exponential transformation
              dexp(exp(theta[nk+3]), vinfl.prior.rate, TRUE) + theta[nk+3])         # additional ANC variance
   } else if(fp$eppmod == "rtrend"){ # rtrend
     
     return(dunif(theta[1], t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
              dunif(theta[2], t1.unif.prior[1], t1.unif.prior[2], log=TRUE) +
              dunif(theta[3], logr0.unif.prior[1], logr0.unif.prior[2], log=TRUE) +
              sum(dnorm(theta[4:7], 0, rtrend.beta.pr.sd, log=TRUE)) +
              #dnorm(theta[8], ancbias.pr.mean, ancbias.pr.sd, log=TRUE) +
              dexp(exp(theta[8]), vinfl.prior.rate, TRUE) + theta[8])   # additional ANC variance
   }
 } else {
   if(!exists("eppmod", where = fp)) {fp$eppmod <- "rspline"}  # backward compatibility
   
   
   if(fp$eppmod == "rspline"){
     nk <- fp$numKnots
     tau2 <- exp(theta[nk+3])
     
     return(sum(dnorm(theta[3:nk], 0, sqrt(tau2), log=TRUE)) +
              dunif(theta[nk+1], logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE) + 
              dnorm(theta[nk+2], ancbias.pr.mean, ancbias.pr.sd, log=TRUE) +
              ldinvgamma(tau2, invGammaParameter, invGammaParameter) + log(tau2) +  # + log(tau2): multiply likelihood by jacobian of exponential transformation
              dexp(exp(theta[nk+4]), vinfl.prior.rate, TRUE) + theta[nk+4])         # additional ANC variance
   } else { # rtrend
     
     return(dunif(theta[1], t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
              dunif(theta[2], t1.unif.prior[1], t1.unif.prior[2], log=TRUE) +
              dunif(theta[3], logr0.unif.prior[1], logr0.unif.prior[2], log=TRUE) +
              sum(dnorm(theta[4:7], 0, rtrend.beta.pr.sd, log=TRUE)) +
              dnorm(theta[8], ancbias.pr.mean, ancbias.pr.sd, log=TRUE) +
              dexp(exp(theta[9]), vinfl.prior.rate, TRUE) + theta[9])   # additional ANC variance
   }
 }
 
 
}

################################
####                        ####
####  HH survey likelihood  ####
####                        ####
################################

fnPrepareHHSLikData <- function(hhs, anchor.year = 1970L){
  hhs$W.hhs <- qnorm(hhs$prev)
  hhs$v.hhs <- 2*pi*exp(hhs$W.hhs^2)*hhs$se^2
  hhs$sd.W.hhs <- sqrt(hhs$v.hhs)
  hhs$idx <- hhs$year - (anchor.year - 1)
  hhs$W.logit.hhs <- logit(hhs$prev)
  hhs$sd.logit.hhs <- (hhs$se)/(hhs$prev*(1 - hhs$prev))
  

  hhslik.dat <- subset(hhs, used)
  return(hhslik.dat)
}

fnHHSll <- function(qM, hhslik.dat, likelihood_type, prev_type, means = NULL, cov_mat = NULL){
  if (likelihood_type == "norm") {
    if (prev_type == "prev") {
      means_norm <- hhslik.dat$prev
      var_param <- hhslik.dat$se
    } else if (prev_type == "probit") {
      means_norm <- hhslik.dat$W.hhs
      var_param <- hhslik.dat$sd.W.hhs
    } else if (prev_type == "logit") {
      means_norm <- hhslik.dat$W.logit.hhs
      var_param <- hhslik.dat$sd.logit.hhs
    }
  return(sum(dnorm(means_norm, qM[hhslik.dat$idx], var_param, log = TRUE)))
    } else if (likelihood_type == "multi_norm") {
    return(dmvnorm(means, qM[hhslik.dat$idx], cov_mat, log = TRUE))
  }
}


###########################
####                   ####
####  Full likelihood  ####
####                   ####
###########################

fnCreateLikDat <- function(epp.data, anchor.year=1970L, no.anc){
if(no.anc==TRUE){# adding this use case for running EPP with no ANC data
  likdat <- list(hhslik.dat = fnPrepareHHSLikData(epp.data$hhs, anchor.year=anchor.year))
  likdat$lastdata.idx <- max(likdat$hhslik.dat$idx)
  likdat$firstdata.idx <- min(likdat$hhslik.dat$idx)
  } else {
  likdat <- list(anclik.dat = anclik::fnPrepareANCLikelihoodData(epp.data$anc.prev, epp.data$anc.n, epp.data$anc.used, anchor.year=anchor.year),
                 hhslik.dat = fnPrepareHHSLikData(epp.data$hhs, anchor.year=anchor.year))
  likdat$lastdata.idx <- max(unlist(likdat$anclik.dat$anc.idx.lst), likdat$hhslik.dat$idx)
  likdat$firstdata.idx <- min(unlist(likdat$anclik.dat$anc.idx.lst), likdat$hhslik.dat$idx)
  }
  
  return(likdat)
}

create_rvec <- function(theta, rt){
  if(rt$eppmod == "rhybrid"){
    
    par <- theta[1:4]
    par[3] <- exp(par[3])
    rvec_rlog <- rlogistic(rt$rlogistic_steps, par)
    
    th_rw <- theta[4+1:rt$n_rw]
    
    diff_rlog <- diff(rlogistic(rt$rw_steps, par))
    diff_rw <- rt$dt * th_rw[rt$rw_idx] / sqrt(rt$rw_dk)
    diff_rvec <- (1 - rt$rw_transition) * diff_rlog + rt$rw_transition * diff_rw
    rvec_rw <- cumsum(c(rvec_rlog[length(rvec_rlog)], diff_rvec))
    
    rvec <- c(rvec_rlog, rvec_rw)
    
    return(exp(rvec))
  }
  else
    stop(paste(rt$eppmod, "is not impmented in create_rvec()"))
}


transf_iota <- function(par, fp){
  
  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }
  
  if(exists("logitiota", fp) && fp$logitiota)
    exp(invlogit(par)*diff(logiota.unif.prior) + logiota.unif.prior[1])
  else
    exp(par)  
}




fnCreateParam <- function(theta, fp, no.anc){

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"
  if(fp$eppmod == "rhybrid"){
    
    
    epp_nparam <- fp$rt$n_param+1
    param <- list()
    param$rvec <- create_rvec(theta[1:fp$rt$n_param], fp$rt)
    param$iota <- transf_iota(theta[fp$rt$n_param+1], fp)
    param$v.infl <- exp(theta[fp$rt$n_param+2])
    return(param)
    
    
  } else if(fp$eppmod == "rspline"){
    u <- theta[1:fp$numKnots]
    beta <- numeric(fp$numKnots)
    beta[1] <- u[1]
    if(fp$numKnots >= 2) {
      beta[2] <- u[1]+u[2]
        if(fp$numKnots >= 3) {
          for(i in 3:fp$numKnots)
            beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]  
             
        }
    }
    
    # return(list(rvec = ifelse(all(beta > 0), as.vector(fp$rvec.spldes %*% beta), -abs(as.vector(fp$rvec.spldes %*% beta))),
    if (no.anc==TRUE) {
      return(list(rvec = as.vector(fp$rvec.spldes %*% beta),
                  iota = exp(theta[fp$numKnots+1]),
                  v.infl = exp(theta[fp$numKnots+3])))
    } else {
      return(list(rvec = as.vector(fp$rvec.spldes %*% beta),
                  iota = exp(theta[fp$numKnots+1]),
                  ancbias = theta[fp$numKnots+2],
                  v.infl = exp(theta[fp$numKnots+4])))
    }
    
    
  } else { # rtrend
    if (no.anc == TRUE) {
      return(list(tsEpidemicStart = fp$proj.steps[which.min(abs(fp$proj.steps - theta[1]))], # t0
                  rtrend = list(tStabilize = theta[1]+theta[2],  # t0 + t1
                                r0 = exp(theta[3]),              # r0
                                beta = theta[4:7]),
                  v.infl = exp(theta[8])))
    } else {
    return(list(tsEpidemicStart = fp$proj.steps[which.min(abs(fp$proj.steps - theta[1]))], # t0
                rtrend = list(tStabilize = theta[1]+theta[2],  # t0 + t1
                              r0 = exp(theta[3]),              # r0
                              beta = theta[4:7]),
                ancbias = theta[8],
                v.infl = exp(theta[9])))
    }
  }
}

ll <- function(theta, fp, likdat, no.anc, likelihood_type = likelihood_type, prev_type = prev_type){
  
  
  
  param <- fnCreateParam(theta, fp, no.anc)

  fp <- update(fp, list=param)
  
  if(!exists("eppmod", where=fp) || fp$eppmod == "rspline" || fp$eppmod == "rhybrid"){
    if(min(fp$rvec)<0 || max(fp$rvec)>20){ # Test positivity of rvec
      return(-Inf)
      }
  }
  
  mod <- simmod(fp)
  
  if (prev_type == "prev") {
    qM.all <- prev(mod)
  } else if (prev_type == "probit") {
    qM.all <- qnorm(prev(mod))
  } else if (prev_type == "logit") {
    qM.all <- logit(prev(mod))
  }
  
  if(no.anc==FALSE) {
    ll.anc <- log(anclik::fnANClik(qM.all+fp$ancbias, likdat$anclik.dat, fp$v.infl))
    }
  if (likelihood_type == "norm") {
  ll.hhs <- fnHHSll(qM.all, likdat$hhslik.dat, likelihood_type, prev_type)
  } else if (likelihood_type == "multi_norm") {
    
    ll.hhs <- fnHHSll(qM.all, likdat$hhslik.dat, likelihood_type, prev_type, means = likdat$means, cov_mat = likdat$cov_mat)
  }
  

  if(exists("equil.rprior", where=fp) && fp$equil.rprior){
    rvec.ann <- fp$rvec[fp$proj.steps %% 1 == 0.5]
    equil.rprior.mean <- muSS/(1-pnorm(qM.all[likdat$lastdata.idx]))
    equil.rprior.sd <- sqrt(mean((muSS/(1-pnorm(qM.all[likdat$lastdata.idx - 10:1])) - rvec.ann[likdat$lastdata.idx - 10:1])^2))  # empirical sd based on 10 previous years
    ll.rprior <- sum(dnorm(rvec.ann[(likdat$lastdata.idx+1L):length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts year after last data
  } else {
    ll.rprior <- 0
  }
  if(no.anc==TRUE) {
    return(ll.hhs+ll.rprior)
  } else {
    return(ll.anc+ll.hhs+ll.rprior)
  }
}


##########################
####  IMIS functions  ####
##########################

sample.prior <- function(n, fp, no.anc = no.anc){
 
  if(no.anc == TRUE){
  if(!exists("eppmod", where = fp)) {fp$eppmod <- "rspline"}  # backward compatibility
  
   if(fp$eppmod == "rhybrid") {
     nparam <- fp$rt$n_param+2
     mat <- matrix(NA, n, nparam)
     
     logiota.unif.prior <- log(c(1e-13, 0.0025))
     r0logiotaratio.unif.prior <- c(-25, -5)
     
     logit <- function(p) log(p/(1-p))
     invlogit <- function(x) 1/(1+exp(-x))
     ldinvlogit <- function(x){v <- invlogit(x); log(v) + log(1-v)}
    
     
     ldsamp_iota <- lprior_iota
     
     
     mat[,1:4] <- t(matrix(rnorm(4*n, rlog_pr_mean, rlog_pr_sd), 4))
     mat[,4+1:fp$rt$n_rw] <- rnorm(n*fp$rt$n_rw, 0, rw_prior_sd)  
     mat[,fp$rt$n_param+1] <- sample_iota(n, fp)
     mat[,fp$rt$n_param+2] <- log(rexp(n, vinfl.prior.rate))      # v.infl
     
   }

  if(fp$eppmod == "rspline"){
       
    mat <- matrix(NA, n, fp$numKnots+3)

    ## sample penalty variance
    tau2 <- rexp(n, tau2.prior.rate)                                                # variance of second-order spline differences
    if(fp$numKnots == 1) {
      mat[,1] <- rnorm(n, 0.5, 1) 
    } else {
      mat[,1] <- rnorm(n, 1.5, 1)                                                   # u[1]
    }
                                                         
    if(fp$numKnots > 1) {
      mat[,2:fp$numKnots] <- rnorm(n*(fp$numKnots-1), 0, sqrt(tau2))                # u[2:numKnots]
    }                                                     
    mat[,fp$numKnots+1] <-  runif(n, logiota.unif.prior[1], logiota.unif.prior[2])  # iota
    mat[,fp$numKnots+2] <- log(tau2)                                                # tau2
    mat[,fp$numKnots+3] <- log(rexp(n, vinfl.prior.rate))                           # v.infl

  } 
  if(fp$eppmod == "rtrend") { # r-trend

    mat <- matrix(NA, n, 8)

    mat[,1] <- runif(n, t0.unif.prior[1], t0.unif.prior[2])        # t0
    mat[,2] <- runif(n, t1.unif.prior[1], t1.unif.prior[2])        # t1
    mat[,3] <- runif(n, logr0.unif.prior[1], logr0.unif.prior[2])  # r0
    mat[,4:7] <- rnorm(4*n, 0, rtrend.beta.pr.sd)                  # beta
    mat[,8] <- log(rexp(n, vinfl.prior.rate))                      # v.infl
  }
 } 
  if (no.anc == FALSE) {
   if(!exists("eppmod", where = fp))  # backward compatibility
     fp$eppmod <- "rspline"
   
   if(fp$eppmod == "rspline"){
     
     mat <- matrix(NA, n, fp$numKnots+4)
     
     ## sample penalty variance
     tau2 <- rexp(n, tau2.prior.rate)                  # variance of second-order spline differences
     if(fp$numKnots == 1) {
       mat[,1] <- rnorm(n, 0.5, 1) 
     } else {
       mat[,1] <- rnorm(n, 1.5, 1)                                                   # u[1]
     }
     
     if(fp$numKnots > 1) {
       mat[,2:fp$numKnots] <- rnorm(n*(fp$numKnots-1), 0, sqrt(tau2))                # u[2:numKnots]
     }                                                     
     mat[,fp$numKnots+1] <-  runif(n, logiota.unif.prior[1], logiota.unif.prior[2])  # iota
     mat[,fp$numKnots+2] <-  rnorm(n, ancbias.pr.mean, ancbias.pr.sd)                # ancbias parameter
     mat[,fp$numKnots+3] <- log(tau2)                                                # tau2
     mat[,fp$numKnots+4] <- log(rexp(n, vinfl.prior.rate))                           # v.infl
     
   } 
   if(fp$eppmod == "rtrend") { # r-trend
     
     mat <- matrix(NA, n, 9)
     
     mat[,1] <- runif(n, t0.unif.prior[1], t0.unif.prior[2])        # t0
     mat[,2] <- runif(n, t1.unif.prior[1], t1.unif.prior[2])        # t1
     mat[,3] <- runif(n, logr0.unif.prior[1], logr0.unif.prior[2])  # r0
     mat[,4:7] <- rnorm(4*n, 0, rtrend.beta.pr.sd)                  # beta
     mat[,8] <- rnorm(n, ancbias.pr.mean, ancbias.pr.sd)            # ancbias parameter
     mat[,9] <- log(rexp(n, vinfl.prior.rate))                      # v.infl
   }
 }
  return(mat)
}

prior <- function(theta, fp, no.anc){
  if(is.vector(theta))
    return(exp(lprior(theta, fp, no.anc)))
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(lprior(theta[i,], fp, no.anc))))))
}

likelihood <- function(theta, fp, likdat, no.anc, likelihood_type, prev_type){
  if(is.vector(theta))
    return(exp(ll(theta, fp, likdat, no.anc, likelihood_type, prev_type)))
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(ll(theta[i,], fp, likdat, no.anc, likelihood_type, prev_type))))))
}
