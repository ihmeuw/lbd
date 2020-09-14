##########################################################
# Collate covariates
###########################################################

# (1) Setup -----------------------------------------------
rm(list = ls())
in_dir <- '<<<< FILEPATH REDACTED >>>>'
regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name-ESH', 'dia_sssa',
             'dia_mcaca', 'dia_s_america_n', 'dia_s_america_s', 'dia_central_asia',
             'dia_se_asia', 'dia_malay', 'dia_south_asia-IND', 'dia_mid_east', 'dia_essa', 'IND','MNG')

# (2) Collate ----------------------------------------------
covs <- data.frame(covariate = NULL, measure = NULL, gbd = NULL, region = NULL, included = NULL)
for (reg in regions){
  reg_file <- fread(paste0(in_dir, '/covs_27_', reg, '.csv')) %>%
    filter(include == TRUE) %>%
    select(-V1, -include)
  reg_file$region <- reg
  reg_file$included <- 'TRUE'
  covs <- rbind(covs, reg_file)
}
covs <- as.data.table(covs) %>%
  dcast(covariate + measure + gbd ~ region, value.var = 'included')
covs[is.na(covs)] <- 'FALSE'

# (3) Save -------------------------------------------------
write.csv(covs, paste0(in_dir, 'covariate_summary.csv'))
