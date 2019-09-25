#################################################################################
### Takes in raked or raw draw-level estimates and makes stat summary rasters
## Inputs:
# draw_level_cell_pred: Cells by Draws matrix which is output from predict_mbg() or from rake_predictions()
# mask: Should be the simple_raster
# return_as_raster: If TRUE returns as raster, else as table
# summary_stat: ie mean, cirange, quantile, sd
## Outputs: Summary table or raster of the cell_pred table put in
#################################################################################
make_cell_pred_summary    <- function(draw_level_cell_pred,
                                      mask                 = simple_raster,
                                      return_as_raster     = TRUE,
                                      summary_stat         = 'mean',
                                      ...){

  # make summary
  summ <- apply(draw_level_cell_pred, 1, summary_stat, ...)

  # put it in a raster
  if(return_as_raster){
    yrs = dim(draw_level_cell_pred)[1]/length(cellIdx(mask))
    message(sprintf('Making a RasterBrick with %i layers',yrs))
    summ <- insertRaster(mask,  matrix(summ,  ncol = yrs))
  }


  return(summ)

}



#################################################################################
### Saves post-estimation output files in the proper directory
## Inputs:
## Outputs:
#################################################################################
save_post_est   <- function(x,
                            filetype,
                            filename,
                            indic = indicator){

  output_dir <- paste0("<<<< FILEPATH REDACTED >>>>")


  dir.create(output_dir, showWarnings = FALSE)

  filetype = tolower(filetype)
  if (!filetype %in% c('rdata','raster','csv'))
    stop('filetype argument has to be either rdata or raster or csv')


  if (filetype == 'raster')
    writeRaster(
      x,
      file = paste0(output_dir, '/', indic,'_',filename),
      format = 'GTiff',
      overwrite = TRUE
    )

  if (filetype == 'rdata')
    save(
      x,
      file = paste0(output_dir, '/', indic,'_',filename,'.RData'),
      compress = TRUE
    )

  if (filetype == 'csv')
    write.csv(
      x,
      file = paste0(output_dir, '/', indic,'_',filename,'.csv')
    )


}



#################################################################################
### Pull cell preds
## Inputs:
## Outputs:
#################################################################################

load_cell_preds <- function(indicator_group,
                            indicator,
                            rd = run_date,
                            region,
                            agebin,
                            u5m=FALSE,
                            other='',
                            ageasindic=TRUE){


  if (u5m) {
    if (ageasindic == FALSE) {
      load(paste0("<<<< FILEPATH REDACTED >>>>",
                  indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'NA.RData')) # the 0 are no holdout
    } else {
      load(paste0("<<<< FILEPATH REDACTED >>>>",
                  indicator,'_age',agebin,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'.RData'))
    }
    cell_pred <- cptmp
  } else {
    load(paste0("<<<< FILEPATH REDACTED >>>>",
                indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0.RData')) # the 0 are no holdouts
  }
  return(cell_pred)
}





prep_postest <- function(indicator,
                         indicator_group,
                         run_date,
                         save_objs) {

  # Save a list of objects in a standard location for parallel scripts to pull from
  main_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  temp_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  temp_file <- paste0(temp_dir, "post_est_temp_objs.RData")
  dir.create(temp_dir, showWarnings = F)
  save(list = save_objs, file = temp_file)
}

post_load_combine_save <- function(regions    = strata,
                                   summstats  = c('mean','cirange','upper','lower'),
                                   raked      = c('raked','unraked'),
                                   rf_table   = TRUE,
                                   run_summ   = TRUE,
                                   indic      = indicator,
                                   ig         = indicator_group,
                                   sdir       = sharedir){

  message(paste0("indic: ", indic))
  message(paste0("ig: ", ig))

  rake_addin <- character()
  if ("unraked" %in% raked) {
    lookup_dir <- paste0(sprintf("<<<< FILEPATH REDACTED >>>>"))
    ur <- length(grep(paste0(indic, ".*unraked.*raster.tif"), list.files(lookup_dir)))
    if (ur > 0) rake_addin <- c(rake_addin, unraked = "_unraked")
    if (ur == 0) rake_addin <- c(rake_addin, unraked = "")
  }

  if ("raked" %in% raked) {
    rake_addin <- c(rake_addin, raked = "_raked")
  }

  # loop through and combine all rasters
  message("\nCombining rasters...")
  for (rake in rake_addin) {
    message(names(rake_addin)[which(rake_addin == rake)])
    rr <- rake
    for (ss in summstats) {
      message(paste0('  ',ss))
      rlist <- list()
      for (reg in regions) {
        message(paste0('    ',reg))
        rlist[[reg]] <-
          brick(sprintf("<<<< FILEPATH REDACTED >>>>"))
      }
      if (length(rlist) > 1) rlist <- do.call(raster::merge,unname(rlist)) else rlist <- rlist[[1]]
      if (ss == 'cirange') ssname = 'range' else ssname = ss # naming convention
      save_post_est(rlist, 'raster',
                    paste0(ssname,rr,'_raster'),
                    indic)
    }
  }

  # do rf also
  if (rf_table) {
    message('RF table')
    rflist <- list()
    for (reg in regions) {
      rflist[[reg]] <-
        read.csv(sprintf("<<<< FILEPATH REDACTED >>>>"))
    }
    save_post_est(do.call(rbind.fill,rflist),'csv','rf', indic)
  }

  # make a run summary graph
  if (run_summ) {
    graph_run_summary(run_date = run_date,
                      indicator_group = ig,
                      indicator = indic)
  }
}


clean_after_postest <- function(indicator,
                                indicator_group,
                                run_date,
                                strata,
                                delete_region_rasters = F) {

  # Delete intermediate files that we no longer need
  main_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  temp_dir <- paste0("<<<< FILEPATH REDACTED >>>>")

  # Deleting - be careful!
  unlink(temp_dir, recursive = T)

  # If desired, insert code here to delete other temporary objects (eg. region-specific rasters)
  grep_string <- paste0(indicator, "_(", paste(strata, collapse = "|"), ").*_raster.tif")
  region_rasters <- grep(grep_string,list.files(main_dir), value = T) %>%
    paste0(main_dir, .)
  if (delete_region_rasters == T) unlink(region_rasters)
}
