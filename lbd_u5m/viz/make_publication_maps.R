## #############################################################################
## 
## MAKE PUBLICATION-READY MAPS
## 
## Purpose: Pull output data and create a large series of publication-ready maps
## NOTE: Please run this on the geospatial R Singularity image!
## 
## #############################################################################

library(data.table)

## SET PROGRAM INPUTS HERE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_date <- '<<<< REDACTED >>>>'
out_dir <- paste0('<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()),'/')
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
figs_only <- TRUE
shp_version <- '<<<< REDACTED >>>>'
include_non_gbd <- FALSE


## SOURCE ALL MAPPING FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

repo_path <- '<<<< FILEPATH REDACTED >>>>'
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))


## PROGRAM EXECUTION STARTS HERE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mapping_data_full <- prep_u5m_data_for_mapping(
  run_date      = run_date,
  core_repo     = '<<<< FILEPATH REDACTED >>>>',
  u5m_repo      = repo_path,
  shp_version   = shp_version, 
  start_year    = 2000,
  end_year      = 2017,
  resume        = TRUE,
  cache_results = TRUE,
  include_non_gbd = include_non_gbd
)

if(figs_only){
  ## Make ONLY plots needed for publication
  ad2_topics_17 <- c(
    'under5_mean_raked','neonatal_mean_raked','infant_mean_raked',
    'under5_lower_raked','under5_upper_raked','nn_u5_ratio'
  )
  for(fig_topic in ad2_topics_17){
    message(sprintf("Plotting %s at adm2 for 2017",fig_topic))
    fig <- make_publication_map(
      full_data       = mapping_data_full,
      topic           = fig_topic,
      year            = 2017,
      geo_level       = 'adm2',
      title           = '',
      subtitle        = NULL,
      col_scale_title = 'detect',
      col_breaks      = 'detect',
      col_labels      = 'detect',
      col_scale       = 'detect',
      col_transformation = 'identity',
      add_projection  = TRUE,
      add_masks       = FALSE
    )
    out_file <- paste0(out_dir,'/',fig_topic,'_adm2_2017.png')
    png(out_file, height=5.5, width=12, units='in', res=1000)
    print(fig)
    dev.off()
  }

  # Plot the raster file in 2017
  message("Working on raster plot...")
  fig <- make_publication_map(
    full_data       = mapping_data_full,
    topic           = 'under5_mean_raked',
    year            = 2017,
    geo_level       = 'rast',
    title           = '',
    subtitle        = NULL,
    col_scale_title = 'detect',
    col_breaks      = 'detect',
    col_labels      = 'detect',
    col_scale       = 'detect',
    col_transformation = 'identity',
    add_projection  = TRUE,
    add_masks       = TRUE
  )
  out_file <- paste0(out_dir,'/under5_mean_raked_rast_2017.png')
  png(out_file, height=5.5, width=12, units='in', res=1000)
  print(fig)
  dev.off()

  # Plot the one admin2 plot in 2000
  message("Working on 2000 plot...")
  fig <- make_publication_map(
    full_data       = mapping_data_full,
    topic           = 'under5_mean_raked',
    year            = 2000,
    geo_level       = 'adm2',
    title           = '',
    subtitle        = NULL,
    col_scale_title = 'detect',
    col_breaks      = 'detect',
    col_labels      = 'detect',
    col_scale       = 'detect',
    col_transformation = 'identity',
    add_projection  = TRUE,
    add_masks       = FALSE
  )
  out_file <- paste0(out_dir,'/under5_mean_raked_adm2_2000.png')
  png(out_file, height=5.5, width=12, units='in', res=1000)
  print(fig)
  dev.off()

  # Plot death counts for 2017 at all admin levels
  message("Working on death counts plots...")
  for(ad_lev in c('rast','adm2','adm1','adm0')){
    message(sprintf('  - %s',ad_lev))
    fig <- make_publication_map(
      full_data       = mapping_data_full,
      topic           = "under5_deathcounts_mean_raked",
      year            = 2017,
      geo_level       = ad_lev,
      title           = '',
      subtitle        = NULL,
      col_scale_title = 'detect',
      col_breaks      = 'detect',
      col_labels      = 'detect',
      col_scale       = 'detect',
      col_transformation = 'log10',
      add_projection  = TRUE,
      add_masks       = FALSE # False for all death counts plots
    )
    out_file <- sprintf('%s/under5_deathcounts_%s_2017.png',out_dir,ad_lev)
    png(out_file, height=5.5, width=12, units='in', res=1000)
    print(fig)
    dev.off()
  }

} else {
  ## Make all possible plots
  for(adm in c('adm1','adm2','adm0')){
    # Get all unique topics
    topics <- names(mapping_data_full$dt[[adm]])
    topics <- unique(sapply(topics, function(x) substr(x, 1, nchar(x)-5)))
    years <- c(2000:2017, 'aroc')
    for(topic in topics[3:length(topics)]){
      if(str_detect(topic,'deathcounts')) c_trans <- 'log10' else c_trans <- 'identity'
      for(yr in years){
        topic_yr <- paste0(topic, '_', yr)
        if( topic_yr %in% names(mapping_data_full$dt[[adm]]) ){
          tryCatch(
            {
              out_file <- paste0(out_dir, topic_yr, '_', adm, '.png')
              fig <- make_publication_map(
                full_data       = mapping_data_full,
                topic           = topic,
                year            = yr,
                geo_level       = adm,
                title           = '',
                subtitle        = NULL,
                col_scale_title = 'detect',
                col_breaks      = 'detect',
                col_labels      = 'detect',
                col_scale       = 'detect',
                col_transformation = c_trans,
                add_projection  = TRUE,
                add_masks       = FALSE
              )  
              png(out_file, height=5.5, width=12, units='in', res=1000)
              print(fig)
              dev.off()
              message(sprintf("  Finished with %s", topic_yr))
            },
            error = function(e){message(sprintf("  %s failed.",topic_yr))},
            warning = function(w){}
          )
        }
      }
    }
    message(sprintf("** FINISHED WITH %s **", adm))
  }

  for(rast_topic in names(mapping_data_full$rast$annual)){
    if(str_detect(rast_topic,'deathcounts')) c_trans <- 'log10' else c_trans <- 'identity'
    for(yr in 2000:2017){
      out_file <- paste0(out_dir, rast_topic,'_',yr,'_rast.png')
      fig <- make_publication_map(
        full_data       = mapping_data_full,
        topic           = rast_topic,
        year            = yr,
        geo_level       = 'rast',
        title           = '',
        subtitle        = NULL,
        col_scale_title = 'detect',
        col_breaks      = 'detect',
        col_labels      = 'detect',
        col_scale       = 'detect',
        col_transformation = c_trans,
        add_projection  = TRUE,
        add_masks       = TRUE
      )  
      png(out_file, height=5.5, width=12, units='in', res=1000)
      print(fig)
      dev.off()
      message(sprintf("  Finished with raster %s", rast_topic))
    }  
  }
}

message("FINISHED WITH ALL MAPS!! *TADA*")
