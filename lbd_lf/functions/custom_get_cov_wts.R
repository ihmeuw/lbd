###########################################################################################
# Getting covariate weights for model runs (w/ custom bar graph output)
###########################################################################################

get.cov.wts <- function(rd, ## run_date
                        ind, ## indicator
                        ind_gp, ## indicator_group
                        reg,
                        age = 0,
                        holdout = 0,
                        vallevel = "", ## validation level if used in qsub
                        stackers_used = c("gam", "gbm", "lasso", "ridge", "enet")) {
  ## ##########################################
  ## load the workspaces and objects we need ##
  ## ##########################################

  pathaddin <- paste0(<<<< FILEPATH REDACTED >>>>)

  this_config <- fread(<<<< FILEPATH REDACTED >>>>)
  stacker_list <- this_config[V1 == "stacked_fixed_effects", V2]
  stackers_used <- strsplit(stacker_list, " ")
  stackers_used <- stackers_used[[1]][stackers_used[[1]] != "+"]

  ## to get the covariate names and data used to fit the model
  ## we also get the data that went into INLA
  ## we get these things from the loaded df object
  load(<<<< FILEPATH REDACTED >>>>)
  fit.data <- fit.data <- as.data.table(df)

  fixed_effects <- this_config[V1 == "fixed_effects", V2][2]
  coefs_sum1 <- as.logical(this_config[V1 == "coefs_sum1", V2])
  selected_fixed_effects <- strsplit(fixed_effects, " ")[[1]][strsplit(fixed_effects, " ")[[1]] != "+"] %>% unique()
  fe_cols <- sapply(selected_fixed_effects, FUN = function(x) grep(x, names(fit.data), value = T)) %>%
    unlist() %>%
    unique()
  fe_cols <- fe_cols[!(fe_cols %in% c("latitude", "longitude"))]
  X <- fit.data[, fe_cols, with = F]

  ## we have now reconstructed the design matrix used in the stackers
  rc <- colnames(X)

  ## we also need to load in the stacker fits
  load(<<<< FILEPATH REDACTED >>>>)

  smo <- child_models

  ## and we need to load the INLA fit
  load(<<<< FILEPATH REDACTED >>>>)
  fit <- res_fit
  waic <- fit$waic$waic

  ## first we strip the design matrix down to just extracted raw covariate columns
  col.idx <- numeric(length(rc))
  for (i in 1:length(rc)) {
    col.idx[i] <- which(colnames(X) == rc[i])
  }
  X <- as.data.frame(X)[, col.idx] ## this is now in order of rc

  ## second we make a matrix to hold the p-values/importances
  imp.mat <- as.data.frame(matrix(ncol = length(rc), nrow = length(smo) + 1))
  colnames(imp.mat) <- rc
  rownames(imp.mat) <- c(names(smo), "INLA COMBINED")

  ## #####################################################################################
  ## now, for each of the models, add the pvalues to the corresponding rows and columns ##
  ## #####################################################################################

  ## ~~~~~
  ## gam ~
  ## ~~~~~
  gam <- smo[["gam"]]
  if (!is.null(gam) & class(gam)[1] != "try-error") {
    index <- grep("gam", row.names(imp.mat))
    smoothed <- summary(gam)$s.table[, 4] ## smoothed table
    for (i in 1:length(smoothed)) {
      names(smoothed)[i] <- substr(names(smoothed)[i], 3, (nchar(names(smoothed)[i]) - 1))
    }
    unsmoothed <- summary(gam)$p.table[, 4] ## parametric table
    all.p <- c(smoothed, unsmoothed)

    ## now match names and stick into our pval.matrix
    for (i in 1:length(rc)) {
      idx <- which(rc[i] == names(all.p))
      imp.mat[index, i] <- all.p[idx]
    }

    ##  convert from p-val to `importance`
    imp.mat[index, ] <- -log(imp.mat[index, ])
    imp.mat[index, ][imp.mat[index, ] == Inf] <- max(imp.mat[index, ][imp.mat[index, ] != Inf])
    imp.mat[index, ] <- imp.mat[index, ] / sum(imp.mat[index, ], na.rm = T)
  } else {
    ## leave the row full of NAs
  }

  ## ~~~~~
  ## gbm ~
  ## ~~~~~
  gbm <- smo[["gbm"]]
  if (!is.null(gbm) & class(gbm)[1] != "try-error") {
    index <- grep("gbm", row.names(imp.mat))
    rel.inf <- gbm$contributions
    for (i in 1:length(rc)) {
      idx <- which(rc[i] == rownames(rel.inf))
      imp.mat[index, i] <- rel.inf$rel.inf[idx]
    }

    ## convert to scaled importance
    imp.mat[index, ] <- imp.mat[index, ] / sum(imp.mat[index, ])
  } else {
    ## leave the row full of NAs
  }

  ## all the penalized regression DO NOT give SD or p-vals...
  ## I 'standardize' the coefs as per a suggestion in this thread:
  ## https://stats.stackexchange.com/questions/14853/variable-importance-from-glmnet
  ## to do this, we first make sure the design matrix is correctly ordered
  if (sum(names(X) != rc)) {
    X <- X[, match(rc, names(X))] ## reorder so everything
    ## matches... but it should be in
    ## order already
  }

  ## ~~~~~~~
  ## lasso ~
  ## ~~~~~~~
  lasso <- smo[["lasso"]]
  if (!is.null(lasso) & class(lasso)[1] != "try-error") {
    index <- grep("lasso", row.names(imp.mat))
    l <- lasso$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(lasso$lambda == l)

    sds <- apply(na.omit(X), 2, sd)
    unscaled.coefs <- lasso$beta[, l.idx]
    ## drop gaul code (country fixed effects) if necessary
    if (sum(grepl("gaul_code", names(unscaled.coefs))) > 0) unscaled.coefs <- unscaled.coefs[-which(grepl("gaul_code", names(unscaled.coefs)))]
    ## ensure correct ordering
    unscaled.coefs <- unscaled.coefs[match(names(unscaled.coefs), rc)]
    ## scaled
    unscaled.coefs <- unscaled.coefs[is.na(names(unscaled.coefs)) == F]
    scaled.coefs <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    for (i in 1:length(rc)) {
      idx <- which(rc[i] == names(scaled.coefs))
      imp.mat[index, i] <- scaled.coefs[idx]
    }

    ## convert to scaled importance
    imp.mat[index, ] <- imp.mat[index, ] / sum(imp.mat[index, ])
  } else {
    ## leave the row full of NAs
  }

  ## ~~~~~~~
  ## ridge ~
  ## ~~~~~~~
  ridge <- smo[["ridge"]]
  if (!is.null(ridge) & class(ridge)[1] != "try-error") {
    index <- grep("ridge", row.names(imp.mat))
    l <- ridge$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(ridge$lambda == l)

    sds <- apply(X, 2, sd)
    unscaled.coefs <- ridge$beta[, l.idx]
    ## drop gaul code (country fixed effects) if necerssary
    if (sum(grepl("gaul_code", names(unscaled.coefs))) > 0) unscaled.coefs <- unscaled.coefs[-which(grepl("gaul_code", names(unscaled.coefs)))]
    ## ensure correct ordering
    unscaled.coefs <- unscaled.coefs[match(names(unscaled.coefs), rc)]
    unscaled.coefs <- unscaled.coefs[(which(is.na(names(unscaled.coefs)) == F))]
    ## scaled
    scaled.coefs <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    for (i in 1:length(rc)) {
      idx <- which(rc[i] == names(scaled.coefs))
      imp.mat[index, i] <- scaled.coefs[idx]
    }

    ## convert to scaled importance
    imp.mat[index, ] <- imp.mat[index, ] / sum(imp.mat[index, ])
  } else {
    ## leave the row full of NAs
  }

  ## ~~~~~~
  ## enet ~
  ## ~~~~~~
  enet <- smo[["enet"]]
  if (!is.null(enet) & class(enet)[1] != "try-error") {
    index <- grep("enet", row.names(imp.mat))
    l <- enet$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(enet$lambda == l)

    sds <- apply(X, 2, sd)
    unscaled.coefs <- enet$beta[, l.idx]
    ## drop gaul code (country fixed effects) if necessary
    if (sum(grepl("gaul_code", names(unscaled.coefs))) > 0) unscaled.coefs <- unscaled.coefs[-which(grepl("gaul_code", names(unscaled.coefs)))]
    ## ensure correct ordering
    unscaled.coefs <- unscaled.coefs[match(names(unscaled.coefs), rc)]
    ## scaled
    scaled.coefs <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    for (i in 1:length(rc)) {
      idx <- which(rc[i] == names(scaled.coefs))
      imp.mat[index, i] <- scaled.coefs[idx]
    }

    ## convert to scaled importance
    imp.mat[index, ] <- imp.mat[index, ] / sum(imp.mat[index, ])
  } else {
    ## leave the row full of NAs
  }

  ## ##########################################
  ## Now we propagate through the INLA coefs ##
  ## ##########################################

  ## to account for different scaling in INLA we also create SDs of the covariates that go into INLA
  inla.X <- as.data.frame(matrix(
    ncol = length(smo),
    nrow = nrow(X)
  ))
  for (i in 1:ncol(inla.X)) {
    if (!is.null(smo[[i]]) & class(smo[[i]])[1] != "try-error") {
      inla.X[, i] <- df[, grep(sprintf("%s_cv_pred", names(smo)[i]), colnames(df)), with = F]
    }
  }

  inla.sds <- apply(inla.X, 2, sd)
  if (as.logical(coefs_sum1) == F) {
    inla.coefs <- summary(fit)$fixed[, 1]

    ## make sure they are in the right order
    inla.coef.idx <- numeric(length(smo))
    for (i in 1:(length(smo))) {
      inla.coef.idx[i] <- match(names(smo)[i], names(inla.coefs))
    }
    scaled.inla.coefs <- inla.coefs[inla.coef.idx] * inla.sds
  } else {
    inla.coefs <- fit$summary.random$covar$mean
    scaled.inla.coefs <- inla.coefs * inla.sds
  }

  ## add the scaled coefs as a column
  imp.mat <- cbind(imp.mat, c(scaled.inla.coefs, NA))
  colnames(imp.mat)[ncol(imp.mat)] <- "SCALED.INLA.COEFS"

  ## make the INLA weighted rel. impotance
  for (i in 1:(ncol(imp.mat) - 1)) {
    imp.mat[nrow(imp.mat), i] <- sum(scaled.inla.coefs * imp.mat[1:length(smo), i], na.rm = TRUE)
  }

  ## finally scale INLA weights
  imp.mat[nrow(imp.mat), ] <- abs(imp.mat[nrow(imp.mat), ])
  imp.mat[nrow(imp.mat), ] <- imp.mat[nrow(imp.mat), ] / sum(imp.mat[nrow(imp.mat), ], na.rm = TRUE)

  ## save to general output directory
  save(imp.mat, waic,
    file = <<<< FILEPATH REDACTED >>>>
  )

  return(imp.mat)
}
