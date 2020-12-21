.libPaths('<<<< FILEPATH REDACTED >>>>')
library(PresenceAbsence)

auc2 <- function (DATA,
                  st.dev = TRUE,
                  which.model = 1,
                  na.rm = FALSE) {
  if (is.logical(st.dev) == FALSE) {
    stop ("'st.dev' must be of logical type")
  }
  if (is.logical(na.rm) == FALSE) {
    stop ("'na.rm' must be of logical type")
  }
  if (sum(is.na(DATA)) > 0) {
    if (na.rm == TRUE) {
      NA.rows <- apply(is.na(DATA), 1, sum)
      warning (length(NA.rows[NA.rows > 0]), " rows ignored due to NA values")
      DATA <- DATA[NA.rows == 0, ]
    } else {
      return (NA)
    }
  }
  if (length(which.model) != 1) {
    stop ("this function will only work for a single model, 'which.model' must be of length one")
  }
  if (which.model < 1 || round(which.model) != which.model) {
    stop ("'which.model' must be a positive integer")
  }
  if (which.model + 2 > ncol(DATA)) {
    stop ("'which.model' must not be greater than number of models in DATA")
  }
  DATA <- DATA[, c(1, 2, which.model + 2)]
  DATA[DATA[, 2] > 0, 2] <- 1
  OBS <- DATA[, 2]
  PRED <- DATA[, 3]
  if (length(OBS[OBS == 1]) == 0 || length(OBS[OBS == 1]) == 
      nrow(DATA)) {
    if (st.dev == FALSE) {
      return (NaN)
    } else {
      return (data.frame(AUC = NaN, AUC.sd = NaN))
    }
  }
  rm(DATA)
  PRED.0 <- PRED[OBS == 0]
  PRED.1 <- PRED[OBS == 1]
  N <- length(PRED)
  n0 <- as.double(length(PRED.0))
  n1 <- as.double(length(PRED.1))
  R <- rank(PRED, ties.method = "average")
  R0 <- R[OBS == 0]
  R1 <- R[OBS == 1]
  U <- n0 * n1 + (n0 * (n0 + 1)) / 2 - sum(R0)
  AUC <- U / (n0 * n1)
  
  rm(PRED)
  rm(OBS)
  if (st.dev == FALSE) {
    return (AUC = AUC)
  } else {
    RR0 <- rank(PRED.0, ties.method = "average")
    RR1 <- rank(PRED.1, ties.method = "average")
    pless.0 <- (R0 - RR0) / n1
    pless.1 <- (R1 - RR1) / n0
    var.0 <- var(pless.0)
    var.1 <- var(pless.1)
    var.AUC <- (var.0 / n0) + (var.1 / n1)
    st.dev.AUC <- var.AUC ^ 0.5
    return (data.frame(AUC = AUC, AUC.sd = st.dev.AUC))
  }
}

rmse <- function(truth, prediction)
  # root mean squared error of prediction from true probability
{
  sqrt(mean(abs(prediction - truth) ^ 2))
}

devBern <- function (truth, prediction)
  # predictive deviance from a bernoulli distribution
{
  -2 * sum(dbinom(truth, 1, prediction, log = TRUE))
}

calcStats <- function(df) {
  
  # if any elements of df are NAs, return NAs
  if (any(is.na(df))) {
    
    results <- c(deviance = NA,
                 rmse = NA,
                 kappa = NA,
                 auc = NA,
                 sens = NA,
                 spec = NA,
                 pcc = NA,
                 kappa_sd = NA,
                 auc_sd = NA,
                 sens_sd = NA,
                 spec_sd = NA,
                 pcc_sd = NA,
                 thresh = NA)
    
  } else {
    
    # add an id column (needed for PresenceAbsence functions)
    df <- data.frame(id = 1:nrow(df), df)
    
    # ~~~~~~~~~~~
    # copntinuous probability metrics
    
    # bernoulli deviance
    dev <- devBern(df[, 2], df[, 3])
    
    # root mean squared error
    rmse <- rmse(df[, 2], df[, 3])
    
    # auc (using my safe version - see above)
    auc <- auc2(df, st.dev = TRUE)
    
    # ~~~~~~~~~~~  
    # discrete classification metrics
    
    # calculate the 'optimum' threshold - one at which sensitivity == specificity
    opt <- optimal.thresholds(df, threshold = 101, which.model = 1, 
                              opt.methods = 3)
    
    # create confusiuon matrix at this threshold
    confusion <- cmx(df, threshold = opt[1, 2])
    
    # kappa (using threshold)
    kappa <- Kappa(confusion, st.dev = TRUE)
    
    # sensitivity and specificity using threshold
    sens <- sensitivity(confusion, st.dev = TRUE)
    spec <- specificity(confusion, st.dev = TRUE)
    
    # proportion correctly classified using threshold
    pcc <- pcc(confusion, st.dev = TRUE)
    
    # create results vector
    results <- c(deviance = dev,
                 rmse = rmse,
                 kappa = kappa[, 1],
                 auc = auc[, 1],
                 sens = sens[, 1],
                 spec = spec[, 1],
                 pcc = pcc[, 1],
                 kappa_sd = kappa[, 2],
                 auc_sd = auc[, 2],
                 sens_sd = sens[, 2],
                 spec_sd = spec[, 2],
                 pcc_sd = pcc[, 2],
                 thresh = opt[1, 2],
                 bf = bag_fraction,
                 cv = cv_folds,
                 exp_num = experiment_num)
    
    
  }
  
  # and return it
  return (results)
}
