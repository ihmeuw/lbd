jobnum <- commandArgs()[6]
opt_type <- commandArgs()[7]
lrnr_type <- commandArgs()[8]
bounds_version <- commandArgs()[9]
print(bounds_version)
experiment_version <- commandArgs()[10]

code_dir <- '<<<< FILEPATH REDACTED >>>>'
in_dir <- '<<<< FILEPATH REDACTED >>>>'
out_dir <- paste0(in_dir, lrnr_type, '/', experiment_version)
dir.create(out_dir)
dir.create(paste0(out_dir, '/stats_train_output/'))
dir.create(paste0(out_dir, '/sstats_test_output/'))
dir.create(paste0(out_dir, '/stats_output/'))
dir.create(paste0(out_dir, '/model_output/'))

library(seegSDM, lib.loc = '<<<< FILEPATH REDACTED >>>>')
library(gbm)
source(paste0(code_dir, '/calcStats.R'))
library(data.table)

#################################################################
## CONSTRUCT PYTHON SYSTEM COMMAND TO RUN 'optimizers.py'

run_optimizerPy <- function(funcs.file_path,
                            funcs.file,
                            bounds.file_path,
                            bounds.file,
                            data.file_path,
                            data.loc,
                            optimizer,
                            learner,
                            cv_folds,
                            n_calls,
                            jobnum,
                            col_start,
                            col_response)
{
  system(
    paste0('<<<< FILEPATH REDACTED >>>>', 
           funcs.file_path, funcs.file,
           bounds.file_path, bounds.file, 
           data.file_path, data.loc,
           optimizer, 
           learner,
           cv_folds,
           n_calls,
           jobnum,
           col_start,
           col_response)
  )
}
#################################################################

# save bounds file for a record
bnds_file <- read.csv(paste0(code_dir, '/space_bounds', bounds_version, '.csv'))
write.csv(bnds_file, paste0(out_dir, '/space_bounds_settings.csv'))

# read in covariates and data
dat_all <- read.csv(paste0(in_dir, 'data/', jobnum, '.csv'))

## MAKE .CSV OF OPTIMIZED HPARS

# split to make OoS hold-out data
indx <- sample(nrow(dat_all), 0.80*nrow(dat_all))
data_train <- dat_all[indx,]
data_test <- dat_all[-indx,]
write.csv(data_train, file = paste0(out_dir, '/data_train_', jobnum, '.csv'), row.names=F)

# get col names for covs
first_cov <- toString(which(names(dat_all) == 'access2'))
response_col <- toString(which(names(dat_all) == 'rate'))

# run optimization in python via system
run_optimizerPy(funcs.file_path <- code_dir,
                funcs.file <- '/optimizers.py ',
                bounds.file_path <- code_dir,
                bounds.file <- paste0('/space_bounds', bounds_version, '.csv '),
                data.file_path <- paste0(out_dir, ' '),
                data.loc <- paste0(out_dir, '/data_train_', jobnum, '.csv '),
                optimizer <- paste0(opt_type, ' '),
                learner <- paste0(lrnr_type, 'R '),
                cv_folds <- '10 ',
                n_calls <- '150 ',
                jobnum <- paste0(jobnum, ' '),
                col_start <- paste0(first_cov, ' '),
                col_response <- response_col)
