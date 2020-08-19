#### Normal Run Launcher ####
# set repo
rm(list = ls())

setwd("<<<< FILEPATH REDACTED >>>>")

# User Inputs
## set node preference
nodes <- 'geos'
## ihme username
user <- "<<<< USERNAME REDACTED >>>>"
# script
script <- '<<<< FILEPATH REDACTED >>>>/mbg/wash/09_fit_statistics.R'
# Define nodes
r_shell <- '<<<< FILEPATH REDACTED >>>>/mbg/mbg_central/share_scripts/shell_sing.sh'
# define comp params
cores = 2

comp_params <- read.csv('00_comp_params.csv')

regions <- c('w_imp_cr')

groups <- c('sani', 'water')


    for (rr in regions) {

        hours <- 12
        ram <- 128
        jname <- paste(rr, "fit_stats", sep = "_")
        sys.sub <- paste0("qsub -now no -l fthread=", cores, " -l m_mem_free=", ram,
            "G -P proj_geo_nodes -q geospatial.q -l archive=TRUE -l h_rt=", hours,
            ":00:00:00 -N ", jname,
          " -e <<<< FILEPATH REDACTED >>>> -o <<<< FILEPATH REDACTED >>>>")
        args <- paste(rr)
        print(paste(rr, ram, hours))
        print(system(paste(sys.sub, r_shell, 1, script, args)))
        Sys.sleep(0.5)
    }
