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

admin0s_to_map <- get_adm0_codes('IND')
mapping_bounds <- c('S'=3.5, 'N'=39, 'W'=66, 'E'=105)

mapping_topics <- c('under5_mean_raked','infant_mean_raked','neonatal_mean_raked')
levels <- c('adm1','adm2','rast')
years <- c(2000, 2005, 2010, 2017)

for(year in years){
  for(level in levels){
    for(topic in mapping_topics){
      message(sprintf("Plotting %s at %s for %s...",topic,level,year))
      fig <- make_publication_map(
        full_data       = mapping_data_full,
        topic           = topic,
        year            = year,
        geo_level       = level,
        title           = '',
        subtitle        = NULL,
        adm0_restriction = admin0s_to_map,
        mapping_bounds   = mapping_bounds,
        col_scale_title = 'detect',
        col_breaks      = 'detect',
        col_labels      = 'detect',
        col_scale       = 'detect',
        col_transformation = 'identity',
        add_projection  = FALSE,
        add_masks       = FALSE,
        add_disputed    = FALSE,
        add_bra_mex     = FALSE,
        point_size      = 0.35
      )
      out_file <- sprintf('%s/%s_%s_%s.png', out_dir, topic, level, year)
      png(out_file, height=6, width=8, units='in', res=600)
      print(fig)
      dev.off()      
    }
  }
}
