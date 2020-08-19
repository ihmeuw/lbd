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
script <- '<<<< FILEPATH REDACTED >>>>/mbg/wash/99_make_ui_ras.R'
# Define nodes
r_shell <- '<<<< FILEPATH REDACTED >>>>/mbg/wash/00_r_shell.sh'
# define comp params
cores = 2

comp_params <- read.csv('00_comp_params.csv')

regions <- c('dia_s_america-per', 'dia_chn_mng')

groups <- c('w_piped', 'w_imp_cr', 'w_unimp_cr', 's_piped', 's_imp_cr',
  's_unimp_cr')

for (gg in groups) {
    for (rr in regions) {
        hours <- 12
        ram <- 128
        jname <- paste(gg, rr, "ui_ras", sep = "_")
        sys.sub <- paste0("qsub -now no -l fthread=", cores, " -l m_mem_free=", ram,
            "G -P proj_geo_nodes -q geospatial.q -l archive=TRUE -l h_rt=", hours,
            ":00:00:00 -N ", jname,
          " -e <<<< FILEPATH REDACTED >>>> -o <<<< FILEPATH REDACTED >>>>")
        args <- paste(gg, rr)
        print(paste(gg, rr, ram, hours))
        print(system(paste(sys.sub, r_shell, 1, script, args)))
        Sys.sleep(0.5)
    }
}
