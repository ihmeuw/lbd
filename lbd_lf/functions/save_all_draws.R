### This script generates and saves all INLA draws to file, for use by parallel runs of predict_mbg(). This is intended to reduce
### peak memory demands when predicting out to a large area and timescale

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save_all_draws <- function (model_fit, samples, rd, outputdir) {
    message(paste0("Generating ", samples, " posterior draws from INLA model fit..."))
    
    ### Generate draws
    draws <- inla.posterior.sample(samples, model_fit)
    
    # Save full inla draws object above
    run_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
    message(paste0("Now saving posterior draws to ", run_dir))
    dir.create(run_dir, showWarnings = FALSE)
    save(list = "draws", file = paste0(<<<< FILEPATH REDACTED >>>>))
}
