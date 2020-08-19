##########################################################################################
############################  LF age and diagnostic crosswalk  ###########################
##########################################################################################

##################################################################################
###### I. Setup
##################################################################################

source(<<<< FILEPATH REDACTED >>>>)
load_from_parallelize()
message(paste0("Using ", core_repo))

user <- Sys.info()[["user"]] ## Get current user name
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- sprintf(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))
package_list <- c(package_list, "sf")
path <- paste0(<<<< FILEPATH REDACTED >>>>)

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

library(jsonlite)
mbg_setup(package_list = package_list, repos = core_repo)

library(ggthemes)
library(fasterize)
library(fda)
library(subplex)
library(spdep)
library(boot)

message(paste0("Boostrap replicate ", replicate))
data_final <- bootstrap_samples[[replicate]]
  
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
  
  # Enforce monotonicity (always increasing) between specific knots (mf only)
  if (diagnostic_test == "mf") {
    for (a in 2:2) {
      alphas[a] <- betas[1] + sum(exp(betas[2:a]))
    }
  }
  
  return(alphas)
}

#### Calculate prevalence
calculatePrevalence <- function(AgeSpline, scale) {
  return((ilogit(logit(0.00001) + AgeSpline + scale))) # Specify near-zero intercept
}

fit_age_crosswalk <- function(Data, test, iter, use_freq = TRUE, age_threshold = 49.5, restart_run = NULL, burnin = NULL) {
  message(paste0("Fitting age crosswalk model for diagnostic type ", test))
  write.csv(Data, file=paste0(<<<< FILEPATH REDACTED >>>>))
  
  #### Set up age vector for spline fits
  avec <- (0:94)
  
  #### Grab unique study-site levels
  USiteTemp <- Data$cohort
  USite <- unique(USiteTemp)
  
  ### Sum up counts by age for identification of spline knot locations
  NumByAge <- as.data.table(cbind(age = avec, num = (0 * avec)))
  for (i in 1:nrow(Data)) {
    Start <- Data$age_start[i]
    End <- Data$age_end[i]
    Test <- Data$sample_size[i]
    
    # Efit this to pull age distribution by location and year
    PGroup <- Data[i, which(colnames(Data) == paste0("age_prop.", Start)):which(colnames(Data) == paste0("age_prop.", End))]
    PIndvRel <- PGroup / sum(PGroup)
    NumAge <- data.table(age = Start:End, num = t(Test * PIndvRel))
    colnames(NumAge)[2] <- "num"
    for (j in Start:End) {
      NumByAge[age == j, "num"] <- NumByAge[age == j, num] + NumAge[age == j, num]
    }
  }
  
  avec <- c(0, (0:94) + 0.5)
  
  NumByAge_mod <- copy(NumByAge)
  
  PerByAge <- data.table("age" = NumByAge_mod$age, "perc" = (NumByAge_mod$num / sum(NumByAge_mod$num)))
  CUMSUM <- data.table("age" = PerByAge$age, "cum_sum" = cumsum(PerByAge$perc))
  
  SplineNum <- 4
  BREAKS <- c(3, 6, 65)
  for (i in 1:(SplineNum - 2)) {
    BREAKS <- c(BREAKS, CUMSUM[max(which(CUMSUM$cum_sum < (i / (SplineNum - 1)))), age])
  }
  
  BREAKS <- sort(unique(BREAKS))
  
  SplineNum <- length(BREAKS) + 2
  
  print(BREAKS)
  
  #### Create bspline basis and calculate numbers of parameters
  Basis <- create.bspline.basis(rangeval = range(0, 94.5), breaks = BREAKS)
  ParNum <- SplineNum + length(USite)
  
  #### Calculate negative log likelihood
  NegLL <- function(vec) return(-LL(vec))
  
  #### Calculate likelihood
  LL <- function(vec) {
    alphas <- applyBasisConstraints(vec[1:(SplineNum)])
    AgeSpline <- applySpline(alphas, Basis, avec)
    
    #### Set up study-level scaling factors
    Scales <- vec[(SplineNum + 1):(length(vec))]
    
    #### Calculate likelihood for each study-site
    ll <- 0
    for (i in 1:length(USite)) {
      AgeLoc <- data.table("age" = avec, "prev" = calculatePrevalence(AgeSpline, Scales[i]))
      colnames(AgeLoc)[2] <- "prev"
      
      #### Constrain prevalence to be flat beyond a threshold age
      AgeLoc[age > age_threshold, prev := AgeLoc[age == age_threshold, prev]]
      
      locs <- which(USiteTemp == USite[i])
      
      #### For each location, calculate the likelihood of the data given the prevalance-by-age curve formed from the parameters
      for (j in 1:length(locs)) {
        Start <- Data$age_start[locs[j]]
        End <- Data$age_end[locs[j]]
        Test <- Data$sample_size[locs[j]]
        Pos <- Data$cases[locs[j]]
        PGroup <- as.vector(as.matrix(Data[locs[j], which(colnames(Data) == paste0("age_prop.", Start)):which(colnames(Data) == paste0("age_prop.", End))]))
        P <- sum(AgeLoc[(age - 0.5) %in% Start:End, prev] * PGroup / sum(PGroup))
        
        ll <- ll + dbinom(Pos, Test, P, log = TRUE)
      }
    }
    
    return(ll)
  }
  
  #### Fit model
  InitialGuess <- c(logit(c(rep(0.5, SplineNum))), logit(rep(0.5, length(USite))))
  
  if (use_freq) { #### Fit using maximum likelihood
    LLim_start <- c(rep(-5, SplineNum), rep(-5, length(USite)))
    ULim_start <- c(rep(5, SplineNum), rep(5, length(USite)))
    
    Optim <- NULL
    while(is.null(Optim)) {
      starting_values <- matrix(nrow = 1, ncol = length(LLim_start))
      for (i in 1:length(LLim_start)) {
        starting_values[i] <- runif(1, min = LLim_start[[i]], max = ULim_start[[i]])
      }
      try(Optim <- optim(starting_values, NegLL, control = list(maxit = final_iter, trace = 6), method = "L-BFGS-B"))
    } 
    return(list("BREAKS" = BREAKS, "Optim" = Optim, "avec" = avec, "Data" = Data))
  }
}


##################################################################################
###### VI. Infer ICT age crosswalk model
##################################################################################

# Set options for all age crosswalk models
if (diagnostic_test == "mf") {
  age_threshold <- 99.5
} else { #ICT/FTS
  age_threshold <- 64.5
}

if (as.logical(load_and_continue)) {
  load(paste0(<<<< FILEPATH REDACTED >>>>))
  ICT_crosswalk_model <- fit_age_crosswalk(Data = data_final, test = diagnostic_test, iter = final_iter, burnin = burnin_iter, use_freq = use_freq, age_threshold = age_threshold, restart_run = ICT_crosswalk_model)
  
  save(ICT_crosswalk_model, file=paste0(<<<< FILEPATH REDACTED >>>>))
} else {
  ICT_crosswalk_model <- fit_age_crosswalk(Data = data_final, test = diagnostic_test, iter = initial_iter, use_freq = use_freq, age_threshold = age_threshold)
  save(ICT_crosswalk_model, file=paste0(<<<< FILEPATH REDACTED >>>>))
}

message("Run finished!")
message(paste0("Outputs have been saved to <<<< FILEPATH REDACTED >>>>"))
