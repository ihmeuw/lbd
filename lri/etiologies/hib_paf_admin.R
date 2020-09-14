################################################################################
# Calculate Hib PAFs given LBD Hib estimates and Hib vaccine coverage
#
# Katie Welgan
# 16 April 2020
################################################################################


# (1) Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set user arguments
repo                        <- '<<<< FILEPATH REDACTED >>>>'
lri_run_date                <- '2020_01_10_15_18_27'
vax_run_date                <- '2019_12_17_16_03_24_raked'
modeling_shapefile_version  <- '2019_09_10'

# set indicator arguments
indicator_groups         <- list('lri', 'vaccine')
indicators               <- list('has_lri', 'hib3_cov')
run_dates                <- list(lri_run_date, vax_run_date)
measures                 <- list('incidence', '')
suffixes                 <- list('_eb_bin0_0', '_eb_bin0_0')
rks                      <- list(T, T)

# load mbg packages and functions
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>', header=FALSE)))
source(paste0(repo, 'lbd_core/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = repo)

# load custom functions
source(paste0(repo, 'lri/etiologies/functions/format_draw_objects.R'))
source(paste0(repo, 'lri/etiologies/functions/calculate_paf.R'))

# define share_dir
share_dir <- '<<<< FILEPATH REDACTED >>>>'

# create output dir
outdir <- paste0(share_dir, '/etiologies/')
dir.create(outdir, showWarnings = FALSE)

# get modelling regions
indicator <- 'has_lri'
region <- get_output_regions(share_dir)


# (2) Pull draw-level LRI and vaccine aggregates ------------------------------------------------------------

# pull estimate aggregates
dt <- format_admin_results(ind_gp = indicator_groups,
                           ind = indicators,
                           rd = run_dates,
                           measure = measures,
                           suffix = suffixes,
                           rk = rks)

# clean up
dt[, name := NULL]


# Get GBD estimates ------------------------------------------------------------------

# Get scaling factors
message('Getting GBD base PAF draws')
base_paf_draws <- read.csv('<<<< FILEPATH REDACTED >>>>')
ve[, V1 := NULL]

# remove niger state in Nigeria for simplicity (we're just pulling national-level VE for now)
ve <- ve[location_id != 25344]

# add admin 1 and 2 codes to GBD results
ve <- merge(sp_hierarchy_list[, c('ADM0_NAME', 'ADM0_CODE', 'ADM1_CODE', 'ADM2_CODE')], ve,
            by.x = 'ADM0_NAME', by.y = 'location_name', allow.cartesian = T)

# clean
setnames(ve, 'year_id', 'year')
ve <- ve[, c('ADM0_NAME', 'ADM0_CODE', 'ADM1_CODE', 'ADM2_CODE', 'year', 'pred_ve')]
ve <- ve[year >= 2000 & year <= 2017]


## Split draw-level aggregates ------------------------------------------------------------

# loop over admin levels
for (a in 0:2) {
  message('Performing rotavirus adjustment for admin level: ', a)

  # create a copy and subset
  dt_admin <- copy(dt)
  dt_admin <- dt_admin[agg_level == paste0('ADM', a)]
  dt_admin[, agg_level := NULL]

  # create admin column name objects
  ad_code <- paste0('ADM', a, '_CODE')
  all_ad_codes <- paste0('ADM', 0:a, '_CODE')
  all_ad_names <- paste0('ADM', 0:a, '_NAME')

  # merge MBG results with with GBD estimates
  setnames(dt_admin, 'code', c(ad_code))
  dt_admin <- merge(dt_admin,
                    unique(ve[, c(ad_code, 'year', 'pred_ve'), with = FALSE]),
                    by = c(ad_code, 'year'), allow.cartesian = T)

  # multiply to get split estimate and clean up
  message('~>adjusting rota incidence with vaccine coverage')
  dt_admin[, rota := lapply(.SD, get_rota_inc,
                            mbg_diaInc = dia_inc,
                            mbg_vax = vax,
                            gbd_ve = pred_ve),
           .SDcols = 'rota_inc', by = c('year', 'draw', ad_code)]
  dt_admin[, c('dia_inc', 'vax', 'pred_ve', 'rota_inc') := NULL]

  # get mean, upper, and lower
  message('~~>calculating and saving mean, upper, and lower estimates')
  dt_admin[, mean := lapply(.SD, mean, na.rm = T),
           .SDcols = 'rota', by = c(ad_code, 'year')]
  dt_admin[, upper := lapply(.SD, upper),
           .SDcols = 'rota', by = c(ad_code, 'year')]
  dt_admin[, lower := lapply(.SD, lower),
           .SDcols = 'rota', by = c(ad_code, 'year')]

  # save summary and free up space
  ad_summary <- merge(unique(sp_hierarchy_list[, c(all_ad_codes, all_ad_names), with = FALSE]),
                      unique(dt_admin[, c(ad_code, 'year', 'mean', 'upper', 'lower'), with = FALSE]),
                      by = ad_code, allow.cartesian = T)
  write.csv(ad_summary, paste0(outdir, 'vax_adjusted_rotavirus_incidence_adm',
                               tolower(a), '_summary.csv'))
  rm(ad_summary)

  # reshape draw object wide
  message('~~~>reshaping wide')
  dt_admin[, c('mean', 'upper', 'lower') := NULL]
  dt_admin <- dcast(dt_admin, ... ~ draw, value.var = 'rota')

  # check for NAs
  message('TESTING: Percent of NA rows per column is: ', mean(is.na(dt_admin[, V1])), '%')

  # save draw objects
  outpath <- paste0(outdir, 'vax_adjusted_rotavirus_incidence_adm',
                    tolower(a), '_draws.rds')
  message('-- finished making rotavirus PAFs across admin draws. now saving at \n', outpath)
  saveRDS(dt_admin, outpath)
  rm(dt_admin)

} # End loop over admins

# done!
message('Finished creating updated rotavirus PAFS for admin-level draws!')