## #################################################################
## ~~~~~~~~~~~ make table of INLA model results ~~~~~~~~~~~~~~~~~ ##
## #################################################################
## takes in standard model run inputs outputs a table of fixed effect,
## spatio-temporal hyperparameter, and random effects parameter
## summaries.
## note: takes a little while since it has to recreate the
## SPDE INLA object since neither we nor INLA saved that object

model_fit_table_list <- function(regions, rd=run_date, holdout = 0,
                                 age = 0,
                                 ind= indicator,
                                 ind_gp = indicator_group,
                                 sharedir = <<<< FILEPATH REDACTED >>>>){
  ## load models
  require(INLA)
  message(sprintf('Pulling together results for %s models',rd))

  tlist=list()

  for(rr in regions){
    message(sprintf('::on region %s',rr))
    reg  <-  rr

    message("::::loading in pre-INLA objects to get spde")
    pathaddin  <-  paste0(<<<< FILEPATH REDACTED >>>>)
    load(<<<< FILEPATH REDACTED >>>>)

    modnames = c('gam','gbm','ridge','enet','lasso')

    full_raster_list <- cov_list

    for(mm in modnames){
      if(min(na.omit(values(full_raster_list[[mm]][[1]]))) < 0){
        message(sprintf("un-logiting: %s", mm))
        full_raster_list[[mm]] <- ilogit(full_raster_list[[mm]])
      }
    }

    ## for stacking, overwrite the columns matching the model_names so
    ## that we can trick inla into being our stacker
    df = df[,paste0(child_model_names) := lapply(child_model_names,
                                                 function(x) get(paste0(x,'_cv_pred')))]

    ## Create SPDE INLA stack
    input_data <- build_mbg_data_stack(df = df,
                                       fixed_effects = all_fixed_effects,
                                       mesh_s = mesh_s,
                                       mesh_t = mesh_t,
                                       use_ctry_res = use_country_res,
                                       use_nugget = use_inla_nugget)

    spde <- input_data[[2]]
    ## this is what we need!

    message('::::loading in INLA fit\n')
    f <-  sprintf(<<<< FILEPATH REDACTED >>>>)
    res_fit <- readRDS(f)

    ## now we extract what we need from the fit to get transformed spatial params
    res.field <- inla.spde2.result(res_fit, 'space', spde, do.transf=TRUE)

    ## nominal range at 0.025, 0.5, 0.975 quantiles
    range   <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.range.nominal[[1]])
    nom.var <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.variance.nominal[[1]])
    spat.hyps <- rbind(range, nom.var)
    rownames(spat.hyps) <- c('Nominal Range', 'Nominal Variance')

    ## other hyperparmas
    hyps <- summary(res_fit)$hyperpar[-(1:2), ] ## first two rows are
    ## theta1, theta2 which
    ## we have in range and
    ## nom.var

    colnames(spat.hyps) <- colnames(hyps)[3:5]
    ## fixed effects
    fixed <- summary(res_fit)$fixed[,1:6]

    ## combine them all and just keep three quantiles

    all.res <- rbind(fixed[, 3:5],
                     spat.hyps,
                     hyps[, 3:5])
    tlist[[rr]] <- all.res
  }
  return(tlist)
}
