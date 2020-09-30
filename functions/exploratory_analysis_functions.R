# ---------------------------------------------------------------------------------------------
# Functions useful for exploratory analysis of ORS, RHF, ORT, and diarrhea
# ---------------------------------------------------------------------------------------------


# -------------------------------------------------------------------
# Read in and clean data

get_data <- function(indi, run_date, indicator_group = 'ort', ad = 'admin_0') {
  
  # load data
  mydat <- fread('<<< FILEPATH REDACTED >>>')
  
  # clean data
  mydat <- mydat[!is.na(mean)]
  mydat[, cirange := NULL]
  mydat[, indicator := indi]
  
  # end function
  return(mydat)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Calculate fold difference between max and min admin unit

min_max_diff <- function(x) {
  
  # calculate fold difference
  diff <- max(x)/min(x)
  
  # end function
  return(diff)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Calculate relative inequity

rel_dev <- function(x, y) {
  
  # calculate fold difference
  ineq <- (x - y) / y
  
  # end function
  return(ineq)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Calculate absolute inequity

abs_dev <- function(x, y) {
  
  # calculate fold difference
  ineq <- x - y
  
  # end function
  return(ineq)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Correlation coefficient function

# get mean spearman statistic
get_cor_coef <- function(x, y) {
  if (!is.na(mean(x)) | !is.na(mean(y))) { 
    return(cor.test(x, y, method = 'spearman')[['estimate']][['rho']])
  } else {
    return(NA)
  }
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# General deaths averted for comparative risk assessment

# get deaths averted
get_deaths_averted <- function(observed_rate, counterfactual_rate, observed_deaths, rr) {
  
  paf_obs <- (observed_rate*(rr - 1))/(observed_rate*(rr - 1) + 1)
  
  paf_ctf <- (counterfactual_rate*(rr - 1))/(counterfactual_rate*(rr - 1) + 1)
  
  averted <- observed_deaths*paf_ctf - observed_deaths*paf_obs
  
  return(averted)
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Get PAF

get_paf <- function(rate, rr) {
  
  paf <- (rate*(rr - 1))/(rate*(rr - 1) + 1)
  
  return(paf)
  
}

# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Deaths averted across four indicators

# get total deaths averted
get_total_deaths_averted <- function(observed_deaths, 
                                     paf_O1, paf_O2, paf_O3, paf_O4,
                                     paf_C1, paf_C2, paf_C3, paf_C4,
                                     only_keep_pos = FALSE) {
  
  paf_obs <- 1 - (1 - paf_O1) * (1 - paf_O2) * (1 - paf_O3) * (1 - paf_O4)
  
  paf_ctf <- 1 - (1 - paf_C1) * (1 - paf_C2) * (1 - paf_C3) * (1 - paf_C4)
  
  averted <- observed_deaths * paf_ctf - observed_deaths * paf_obs
  
  if (only_keep_pos) averted <- max(0, averted)
  
  return(averted)
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Deaths averted across two indicators

# get total deaths averted
get_half_deaths_averted <- function(observed_deaths, 
                                    paf_O1, paf_O2,
                                    paf_C1, paf_C2,
                                    only_keep_pos = FALSE) {
  
  paf_obs <- 1 - (1 - paf_O1) * (1 - paf_O2)
  
  paf_ctf <- 1 - (1 - paf_C1) * (1 - paf_C2)
  
  averted <- observed_deaths * paf_ctf - observed_deaths * paf_obs
  
  if (only_keep_pos) averted <- max(0, averted)
  
  return(averted)
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# General number of deaths averted by to a given risk factor per 1,000 (rate)

# get ratio averted
get_rate_averted <- function(deaths_averted, population) {
  
  rate <- deaths_averted/population*1000
  
  return(rate)
  
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# General ratio of change in deaths attributable to a given risk factor

# get ratio averted
get_ratio_averted <- function(deaths_averted, change_in_deaths) {
  
  ratio <- deaths_averted/change_in_deaths
  
  return(ratio)
  
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------
# Get correlation at the admin and cell level

make_corr <- function(ind_gp, ind, rd,
                      matrix_pred_name = NULL,
                      type,
                      measure="mortality",
                      skip_cols = NULL,
                      year_list = c(2000:2017),
                      raked,
                      extra_file_tag = '',
                      shapefile_version = 'current') {
  
  # define directories
  share_dir <- '<<< FILEPATH REDACTED >>>'
  output_dir <- paste0(share_dir, "/pred_derivatives/corr/")
  dir.create(output_dir, showWarnings=F, recursive = T)
  
  ## load in regions used in this model and run_date
  regions <- get_output_regions(in_dir = share_dir)
  
  ## get number of years
  num_yrs <- length(year_list)
  
  if('cell' %in% type){
    message('Working on CELL level')
    for(ii in 1:length(regions)){
      message(sprintf('- On region: %s', regions[ii]))
      
      cell_pred_1 <- get_cell_pred_for_aroc(ind_gp,
                                            ind[[1]],
                                            rd[[1]],
                                            regions[ii],
                                            measure,
                                            matrix_pred_name,
                                            skip_cols, rk = raked,
                                            shapefile_version = shapefile_version)
      cell_pred_1 <- cell_pred_1[, 1:12]
      
      cell_pred_2 <- get_cell_pred_for_aroc(ind_gp,
                                            ind[[2]],
                                            rd[[2]],
                                            regions[ii],
                                            measure,
                                            matrix_pred_name,
                                            skip_cols, rk = raked,
                                            shapefile_version = shapefile_version)
      cell_pred_2 <- cell_pred_2[, 1:12]
      
      # check to make sure the nubmer of rows match and that the index ids match
      if (nrow(cell_pred_1) != nrow(cell_pred_2)) stop('The number of rows in the two cell preds do not match. Make sure the same shapefiles were used for modeling.')
      if (sum(unique(cell_pred_1[, 2]) == unique(cell_pred_2[, 2])) != length(unique((cell_pred_1[, 2])))) stop('The row indexes of your two cell pred objects do not match. Stopping.')
      
      ## now we have some year column and some id column and the rest are draw columns.
      ## we can now calculate correlation by idx across time
      num_draws <- ncol(cell_pred_1) - 2
      num_idx <- length(unique(cell_pred_1[, 2])) ## second col is 'idx'
      idx_names <- unique(cell_pred_1[, 2])
      
      ## for each draw, calculate the rates of change between years and
      ## then the weighted total AROC across all years
      message('-- making correlation across years for each draw')
      corr_draws <- matrix(ncol = num_draws + 1, nrow = num_idx)
      corr_draws[, 1] <- idx_names
      
      # loop over draws 
      for(dd in 3:ncol(cell_pred_1)){
        
        # create correlation matrix
        corr_mat <- matrix(ncol = 1, nrow = num_idx)
        
        # lapply over indxes
        get_idx_correlation <- function(xx) {
          corr_val <- get_cor_coef(subset(cell_pred_1[cell_pred_1[, 'idx'] == idx_names[[xx]], ], select = dd), 
                                   subset(cell_pred_2[cell_pred_2[, 'idx'] == idx_names[[xx]], ], select = dd))
          corr_mat[xx, ] <- corr_val
        }
        mclapply(1:length(idx_names), get_idx_correlation, mc.cores = detectCores())
        
        # add to draw matrix
        corr_vec <- as.vector(corr_mat)
        corr_draws[, dd-1] <- corr_vec
        
      }
      
      # finish up and save
      message(sprintf('TESTING: Percent of NA rows per column is: %f%%', mean(is.na(corr_draws[, 1]))))
      message('-- finished making correlation across draws. now saving')
      saveRDS(object = corr_draws,
              file = sprintf('%s/%s_%s_corr_cell_draw_matrix_%s%s.RDs',
                             output_dir, ind[[1]], measure, regions[ii],
                             extra_file_tag))
    } # Close region loop
  }
  
  if('admin' %in% type){
    
    message('Working on ADMIN level')
    ## load the admin objects
    ## try two different locations until we standardize
    file_1 <- '<<< FILEPATH REDACTED >>>'
    file_2 <- '<<< FILEPATH REDACTED >>>'
    load(file_1)
    admins_1 <- list(admin_0, admin_1, admin_2)
    load(file_2)
    admins_2 <- list(admin_0, admin_1, admin_2)
    rm(admin_0, admin_1, admin_2)
    
    ## load spatial admin hierarchy
    admins <- get_sp_hierarchy(shapefile_version = shapefile_version)
    
    ## this contains admin_0, admin_1, and admin_2 matrix_draw objects
    ## col1 is year, col2 is ADM*_CODE, final column is pop. all other columns are draws
    
    for(aa in 0:2){ ## for each admin type
      message(sprintf('- On admin: %i', aa))
      
      cell_pred_1 <- admins_1[[aa+1]]
      cell_pred_1 <- as.data.table(cell_pred_1)
      cell_pred_1[, c('gam', 'gbm', 'enet') := NULL]
      cell_pred_2 <- admins_2[[aa+1]]
      cell_pred_2 <- as.data.table(cell_pred_2)
      cell_pred_2[, c('gam', 'gbm', 'enet') := NULL]
      
      # check to make sure the nubmer of rows match and that the index ids match
      if (nrow(cell_pred_1) != nrow(cell_pred_2)) stop('The number of rows in the two cell preds do not match. Make sure the same shapefiles were used for modeling.')
      if (sum(unique(cell_pred_1[, 2]) == unique(cell_pred_2[, 2])) != nrow(unique((cell_pred_1[, 2])))) stop('The row indexes of your two cell pred objects do not match. Stopping.')
      
      # get ids
      num_idx <- length(unique(cell_pred_1$ADM0_CODE))
      idx_names <- unique(cell_pred_1$ADM0_CODE)
      num_draws <- ncol(cell_pred_1) - 2
      
      ## for each draw, calculate the correlation between model outputs
      message('-- making correlation between years for each draw')
      corr_draws <- matrix(ncol = num_draws + 1, nrow = num_idx)
      corr_draws[, 1] <- idx_names
      
      # loop over draws 
      for (dd in 1:(num_draws - 1)) {
        
        # create correlation matrix
        corr_mat <- matrix(ncol = 1, nrow = num_idx)
        
        # loop over indxes
        for (xx in 1:length(idx_names)) {
          
          # get correlation coefficient
          corr_val <- get_cor_coef(cell_pred_1[get(paste0('ADM', aa, '_CODE')) == idx_names[[xx]], get(paste0('V', dd))], 
                                   cell_pred_2[get(paste0('ADM', aa, '_CODE')) == idx_names[[xx]], get(paste0('V', dd))])
          corr_mat[xx, ] <- corr_val
        }
        
        # add to draw matrix
        corr_vec <- as.vector(corr_mat)
        corr_draws[, dd+1] <- corr_vec
        
      }
      
      message(sprintf('TESTING: Percent of NA rows per column is: %f%%', mean(is.na(corr_draws[, 1]))))
      message('-- finished making AROC across draws. now saving')
      saveRDS(object = corr_draws,
              file = sprintf('%s/%s_%s_corr_adm%i_draw_matrix%s.RDs',
                             output_dir, ind[1], measure, aa,
                             extra_file_tag))
      
    }
  } 
}


# -------------------------------------------------------------------