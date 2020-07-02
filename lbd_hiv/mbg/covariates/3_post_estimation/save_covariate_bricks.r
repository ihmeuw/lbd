##############################################################
# Script to save covariate bricks once models have been run
##############################################################

## Set repo
core_repo  <- paste0("/lbd_core/")
if (!dir.exists(core_repo)) core_repo <- "/lbd_core/"
setwd(core_repo)

## Load libraries and  MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = core_repo)

table <-
  fread(paste0("<<<< FILEPATH REDACTED >>>>")) %>%
  mutate(cov_name = case_when(indicator == "multiple_partners_year" & data_tag == "_WN" ~ "partners_year_wn",
                              indicator == "multiple_partners_year" & data_tag == "_MN" ~ "partners_year_mn",
                              TRUE ~ indicator)) %>% 
  dplyr::select(indicator, cov_name, run_date) %>%
  as.matrix()

## Save covariate bricks
apply(table, 1, function(x) {
  indicator <- x[[1]]
  cov_name <- x[[2]]
  run_date <- x[[3]]

  save_standard_covariate(cov_name = cov_name,
                          ig = "hiv",
                          run_date = run_date,
                          ind = indicator,
                          raked = F,
                          year_list = 2000:2018,
                          measure = "mean")

})
