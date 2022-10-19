##########################################################################################
############################  Oncho age and diagnostic crosswalk  ########################
##########################################################################################
## This version fits both ss and nod data concurrently, with separate splines for each
## type of diagnostic test, and an indicator for ss observations

##################################################################################
###### I. Setup
##################################################################################
source(<<<< FILEPATH REDACTED >>>>)
load_from_parallelize()
message(paste0("Using ", core_repo))

user <- Sys.info()[["user"]] ## Get current user name
indic_repo <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>), header = FALSE)))
package_list <- c(package_list, "sf")
path <- <<<< FILEPATH REDACTED >>>>

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

library(jsonlite)
mbg_setup(package_list = package_list, repos = core_repo)

library(ggthemes)
library(fasterize)
library(fda)
library(subplex)
library(spdep)
library(boot)

if (perform_bootstrap) {
  message(paste0("Boostrap replicate ", replicate))
  data_final <- bootstrap_samples[[replicate]]
} else {
  message("Running a crosswalk analysis using full data set (bootstrap replicate `0`)")
}

##################################################################################
###### V. Define crosswalk functions.
##################################################################################
#### Calculate the prevalence-age curve in logit space; can put something other than splines in here
applySpline <- function(vec, Basis, avec) {
  AgeSpline <- fd(vec, Basis)
  return(predict(AgeSpline, avec))
}

### Apply constraints to shape of spline
applyBasisConstraints <- function(betas) {
  alphas <- copy(betas)
  
  return(alphas)
}

#### Calculate prevalence
calculatePrevalence <- function(AgeSpline, scale, diag_fe) {
  return((ilogit(logit(0.00001) + AgeSpline + scale + diag_fe))) # Specify near-zero intercept
}

fit_age_crosswalk <- function(Data, test, iter, use_freq = TRUE, age_threshold = 49.5, restart_run = NULL, burnin = NULL) {
  message(paste0("Fitting age crosswalk model for diagnostic type ", test))
  write.csv(Data, file=<<<< FILEPATH REDACTED >>>>)
  
  #### Set up age vector for spline fits
  avec <- (0:94)
  
  #### Grab unique study-site levels
  USiteTemp <- Data$cohort
  USite <- unique(USiteTemp)
  
  ### Sum up counts by age for identification of spline knot locations
  NumByAge_ss <- as.data.table(cbind(age = avec, num = (0 * avec)))
  Data_ss <- Data[diagnostic == "ss"]
  for (i in 1:nrow(Data_ss)) {
    # message(i)
    Start <- Data_ss$age_start[i]
    End <- Data_ss$age_end[i]
    Test <- Data_ss$sample_size[i]
    
    # Efit this to pull age distribution by location and year
    PGroup <- Data_ss[i, which(colnames(Data_ss) == paste0("age_prop.", Start)):which(colnames(Data_ss) == paste0("age_prop.", End))]
    PIndvRel <- PGroup / sum(PGroup)
    NumAge <- data.table(age = Start:End, num = t(Test * PIndvRel))
    colnames(NumAge)[2] <- "num"
    for (j in Start:End) {
      NumByAge_ss[age == j, "num"] <- NumByAge_ss[age == j, num] + NumAge[age == j, num]
    }
  }
  
  NumByAge_nod <- as.data.table(cbind(age = avec, num = (0 * avec)))
  Data_nod <- Data[diagnostic == "nod"]
  for (i in 1:nrow(Data_nod)) {
    # message(i)
    Start <- Data_nod$age_start[i]
    End <- Data_nod$age_end[i]
    Test <- Data_nod$sample_size[i]
    
    # Efit this to pull age distribution by location and year
    PGroup <- Data_nod[i, which(colnames(Data_nod) == paste0("age_prop.", Start)):which(colnames(Data_nod) == paste0("age_prop.", End))]
    PIndvRel <- PGroup / sum(PGroup)
    NumAge <- data.table(age = Start:End, num = t(Test * PIndvRel))
    colnames(NumAge)[2] <- "num"
    for (j in Start:End) {
      NumByAge_nod[age == j, "num"] <- NumByAge_nod[age == j, num] + NumAge[age == j, num]
    }
  }
  
  avec <- c(0, (0:94) + 0.5)
  
  PerByAge_ss <- data.table("age" = NumByAge_ss$age, "perc" = (NumByAge_ss$num / sum(NumByAge_ss$num)))
  CUMSUM_ss <- data.table("age" = PerByAge_ss$age, "cum_sum" = cumsum(PerByAge_ss$perc))
  
  PerByAge_nod <- data.table("age" = NumByAge_nod$age, "perc" = (NumByAge_nod$num / sum(NumByAge_nod$num)))
  CUMSUM_nod <- data.table("age" = PerByAge_nod$age, "cum_sum" = cumsum(PerByAge_nod$perc))
  
  SplineNum <- 4
  BREAKS_ss <- c(3, 6, 65)
  BREAKS_nod <- c(3, 6, 65)
  for (i in 1:(SplineNum - 2)) {
    BREAKS_ss <- c(BREAKS_ss, CUMSUM_ss[max(which(CUMSUM_ss$cum_sum < (i / (SplineNum - 1)))), age])
    BREAKS_nod <- c(BREAKS_nod, CUMSUM_nod[max(which(CUMSUM_nod$cum_sum < (i / (SplineNum - 1)))), age])
  }
  
  BREAKS_ss <- sort(unique(BREAKS_ss))
  BREAKS_nod <- sort(unique(BREAKS_nod))
  
  SplineNum <- length(BREAKS_ss) + 2
  
  print(BREAKS_ss)
  print(BREAKS_nod)
  
  #### Create bspline basis and calculate numbers of parameters
  Basis_ss <- create.bspline.basis(rangeval = range(0, 94.5), breaks = BREAKS_ss)
  Basis_nod <- create.bspline.basis(rangeval = range(0, 94.5), breaks = BREAKS_nod)
  
  #### Calculate negative log likelihood
  NegLL <- function(vec) return(-LL(vec))
  
  #### Calculate likelihood
  LL <- function(vec) {
    ## Skin snip
    alphas_ss <- applyBasisConstraints(vec[1:SplineNum])
    AgeSpline_ss <- applySpline(alphas_ss, Basis_ss, avec)
    
    ## Nodules
    alphas_nod <- applyBasisConstraints(vec[(SplineNum + 1):(2 * SplineNum)])
    AgeSpline_nod <- applySpline(alphas_nod, Basis_nod, avec)
    
    #### Set up study-level scaling factors
    Scales <- vec[((2 * SplineNum) + 1):(length(vec) - 1)]
    
    #### Calculate likelihood for each study-site
    ll <- 0
    for (i in 1:length(USite)) {
      AgeLoc_ss <- data.table("age" = avec, "prev" = calculatePrevalence(AgeSpline_ss, Scales[i], vec[length(vec)])) # add skin snip fixed effect
      colnames(AgeLoc_ss)[2] <- "prev"
      
      AgeLoc_nod <- data.table("age" = avec, "prev" = calculatePrevalence(AgeSpline_nod, Scales[i], 0))
      colnames(AgeLoc_nod)[2] <- "prev"
      
      #### Constrain prevalence to be flat beyond a threshold age
      AgeLoc_ss[age > age_threshold, prev := AgeLoc_ss[age == age_threshold, prev]]
      AgeLoc_nod[age > age_threshold, prev := AgeLoc_nod[age == age_threshold, prev]]
      locs <- which(USiteTemp == USite[i])
      
      #### For each location, calculate the likelihood of the data given the prevalance-by-age curve formed from the parameters
      for (j in 1:length(locs)) {
        Diagnostic <- Data$diagnostic[locs[j]]
        Start <- Data$age_start[locs[j]]
        End <- Data$age_end[locs[j]]
        Test <- Data$sample_size[locs[j]]
        Pos <- Data$cases[locs[j]]
        PGroup <- as.vector(as.matrix(Data[locs[j], which(colnames(Data) == paste0("age_prop.", Start)):which(colnames(Data) == paste0("age_prop.", End))]))
        
        if (Diagnostic == "ss") {
          P <- sum(AgeLoc_ss[(age - 0.5) %in% Start:End, prev] * PGroup / sum(PGroup))
        } else if (Diagnostic == "nod") {
          P <- sum(AgeLoc_nod[(age - 0.5) %in% Start:End, prev] * PGroup / sum(PGroup))
        }
        
        ll <- ll + dbinom(Pos, Test, P, log = TRUE)
      }
    }
    
    return(ll)
  }
  
  #### Fit model
  InitialGuess <- c(logit(c(rep(0.5, SplineNum * 2))), logit(rep(0.5, length(USite))), logit(0.5))
  
  if (use_freq) { #### Fit using maximum likelihood
    ## Set upper and lower starting limits
    LLim_start <- c(rep(-5, SplineNum * 2), rep(-5, length(USite)), -5)
    ULim_start <- c(rep(5, SplineNum * 2), rep(5, length(USite)), 5)
    
    Optim <- NULL
    while(is.null(Optim)) {
      starting_values <- matrix(nrow = 1, ncol = length(LLim_start))
      for (i in 1:length(LLim_start)) {
        starting_values[i] <- runif(1, min = LLim_start[[i]], max = ULim_start[[i]])
      }
      try(Optim <- optim(starting_values, NegLL, control = list(maxit = final_iter, trace = 6), method = "L-BFGS-B"))
    } 
    return(list("BREAKS_ss" = BREAKS_ss, "BREAKS_nod" = BREAKS_nod, "Optim" = Optim, "avec" = avec, "Data" = Data))
  } else { #### Optionally fit model using Bayesian MCMC
    LLim <- c(rep(-10, SplineNum), rep(-15, length(USite)))
    ULim <- c(rep(10, SplineNum), rep(15, length(USite)))
    
    LLim_start <- c(rep(-5, SplineNum), rep(-5, length(USite)))
    ULim_start <- c(rep(5, SplineNum), rep(5, length(USite)))
    
    starting_values <- matrix(nrow = 3, ncol = length(LLim))
    for (j in 1:3) {
      for (i in 1:length(LLim)) {
        starting_values[j, i] <- runif(1, min = LLim_start[[i]], max = ULim_start[[i]])
      }
    }
    
    setup <- createBayesianSetup(LL, lower = LLim, upper = ULim, best = InitialGuess)
    
    if (is.null(restart_run)) {
      settings <- list(iterations = iter,  message = TRUE, startValue = starting_values)
      out <- runMCMC(setup, sampler = "DEzs", settings = settings)
    } else {
      if (!is.null(burnin)) {
        x = getSample(restart_run, start = burnin)
        meansPost = apply(x, 2, mean)
        rangePost = apply(x, 2, range)
        newZ = matrix(runif(500 * ncol(x), rangePost[1,], rangePost[2,]), ncol = ncol(x), byrow = T)
        settings = list(Z = newZ, startValue = x[(nrow(x) - 2):nrow(x),], iterations = iter,  message = TRUE)
        out <- runMCMC(setup,  sampler = "DEzs", settings = settings)
      } else {
        settings <- list(iterations = iter, message = TRUE)
        out <- runMCMC(restart_run, sampler = "DEzs", settings = settings)
      }
    }
    
    return(list("BREAKS" = BREAKS, "avec" = avec, "Data" = Data, "mcmc" = out, "age_threshold" = age_threshold))
  }
}


##################################################################################
###### VI. Infer age and diagnostic crosswalk model
##################################################################################
# Set options for all age crosswalk models
# age_threshold <- 64.5
age_threshold <- 99

if (as.logical(load_and_continue)) {
  load(paste0(<<<< FILEPATH REDACTED >>>>))
  oncho_crosswalk_model <- fit_age_crosswalk(Data = data_final, test = diagnostic_test, iter = final_iter, burnin = burnin_iter, use_freq = use_freq, age_threshold = age_threshold, restart_run = oncho_crosswalk_model)
  
  save(oncho_crosswalk_model, file=<<<< FILEPATH REDACTED >>>>)
} else {
  oncho_crosswalk_model <- fit_age_crosswalk(Data = data_final, test = diagnostic_test, iter = initial_iter, use_freq = use_freq, age_threshold = age_threshold)
  save(oncho_crosswalk_model, file=<<<< FILEPATH REDACTED >>>>)
}

message("Run finished!")
message(paste0("Outputs have been saved to ", <<<< FILEPATH REDACTED >>>>))
