#' @title Checks for covariate issues
#'
#' @description This function looks for issues with covariates that will cause issues down the line and tries to fix them.
#' The two issues looked for here are uniform covariates and pixel coverage. Uniform Covariates: If one of the extracted variables
#' does not vary over the data area. Pixel Coverage: If one of the covariates is missing in large parts of your model region,
#' these areas will be NA in your results later on. This often happens if you are using a modeled covariate in a new or
#' partially new region which was not previously modelled. This function takes in all parallel model objects that would need to
#' change if a covariate was removed for the code below in the parallel mode to work.
#'
#' @param cc object to change: cs_covs
#' @param afe object to change: all_fixed_effects
#' @param afeb object to change: all_fixed_effects_brt
#' @param acl object to change: all_cov_layers
#' @param tc object to change: the_covs
#' @param check_uniform boolean, do a check for uniform covariates, defaults to TRUE
#' @param check_pixelcount boolean, do a check for covariates with too few pixel (ie too low of geographic coverage in the region), defaults to TRUE
#' @param check_pixelcount_thresh for a pixelcount check, what proportion of the maximum observed pixel coverage is needed to keep the covariate in? must be between 0 and 1. defaults to 0.95
#'
#' @return the necessary objects as a named list that will later get assigned
#' to the environment within the parallel_model.R script
#' @export
check_for_cov_issues   <- function(cc                      = cs_covs,
                                   afe                     = all_fixed_effects,
                                   afeb                    = all_fixed_effects_brt,
                                   fe                      = fixed_effects,
                                   tc                      = the_covs,
                                   acl                     = all_cov_layers,
                                   check_uniform           = TRUE,
                                   check_pixelcount        = TRUE,
                                   check_pixelcount_thresh = 0.95,
                                   drop_nonvarying_covs    = TRUE){
  
  # check that threshold is between 0 and 1
  if(check_pixelcount_thresh < 0 | check_pixelcount_thresh > 1){
    stop('check_pixelcount_thresh must be between 0 and 1')
  }
  
  # make a few useful objects to track names and dropped covs
  covs <- as.character(cc$cs_df$name)
  covs <- covs[!grepl('gaul_code', covs)]
  dropcovs <- c()
  
  if (drop_nonvarying_covs) { # run through extracted data and find any non-varying covariates
    if(check_uniform == TRUE){
      message('Checking for uniform covariates ... ')
      for(covname in covs){
        if (all( abs(na.omit(cc$covs[[covname]]) - mean(na.omit(cc$covs[[covname]]))) == 0 )){
          message(sprintf('WARNING: %s did not vary in your data and is being removed as a covariate', covname))
          dropcovs <- c(dropcovs,covname)
        }
      }
      if(length(dropcovs) == 0){
        message('  All clear')
      }
    }
  }
  
  # check to see if any of the covariates have significant area missing, if so drop them 
  # this typically happens if you use a modelled surface as a covariate in a region that was missing one or more countries
  # count the number of pixels in each covariate, make sure they are within 95% non-NA of the max covered covariate
  if(check_pixelcount == TRUE){
    message('Checking for pixel counts ... ') 
    px_cnt <- c()
    for(covname in covs){
      dimtoget <- max(dim(acl[[covname]])[3]) # grab just the most recent (or synoptic) layer, assuming they are all the same
      px_cnt   <- c(px_cnt, length(cellIdx(acl[[covname]][[dimtoget]]))) 
    }
    threshdrops <- covs[px_cnt/max(px_cnt) < check_pixelcount_thresh]
    if(length(threshdrops) == 0){
      message('  All clear')
    } else {
      message(sprintf('WARNING: the following covariates had < %i per cent pixel coverage in this region and are being dropped:\n %s', 
                      round(check_pixelcount_thresh,2)*100, paste(threshdrops,collapse=', ')))
      dropcovs <- c(dropcovs, threshdrops)
    }
  }
  
  # make sure dropcovs is unique
  dropcovs <- unique(dropcovs)
  
  if (drop_nonvarying_covs) { # if any non-varying covariates were detected, then remove them from any objects
    if(length(dropcovs) > 0){
      for(dc in dropcovs){
        
        # drop from the cs_covs object, which itself is a list
        cc$covs[, (dc) := NULL]
        cc$cs_df <- cc$cs_df[!as.character(cc$cs_df$name) %in% dc,]
        
        # drop from the covariate layers list of rasters
        acl <- replace(acl, dc, NULL)
      }
    } else {
      message('No non-varying covariates were detected in the data. Yay!')
    }
  }
  
  # redo the all fixed effects pseudo-formula strings
  fe   <- paste0(format_covariates(fe)  [!format_covariates(fe)   %in% dropcovs], collapse = ' + ')
  afe  <- paste0(format_covariates(afe) [!format_covariates(afe)  %in% dropcovs], collapse = ' + ')
  afeb <- paste0(format_covariates(afeb)[!format_covariates(afeb) %in% dropcovs], collapse = ' + ')
  tc   <- tc[!tc %in% dropcovs]
  
  # return the objects as a named list that will later get assigned to the environment within the parallel_model.R script
  return(list( cs_covs               = cc,
               all_fixed_effects     = afe,
               all_fixed_effects_brt = afeb,
               fixed_effects         = fe,
               the_covs              = tc,
               all_cov_layers        = acl))
}