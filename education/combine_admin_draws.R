rm(list=ls())
## Create run date in correct format
run_date <- as.character(commandArgs()[4])
indicator <- as.character(commandArgs()[5])
raked <- as.character(commandArgs()[6])
## Set repo locations and indicator group
user <- Sys.info()['user']
message(user)
message(sessionInfo())

core_repo <- paste0("<<<< FILEPATH REDACTED >>>>")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")
indicator_group <- 'education'
rerun_run_date <- NA

sharedir       <- sprintf("<<<< FILEPATH REDACTED >>>>")
commondir      <- sprintf("<<<< FILEPATH REDACTED >>>>")

package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

year_list <- 2000:2017
## Set strata as character vector of each strata (in my case, just stratifying by region whereas U5M stratifies by region/age)
strata <- c('ansa+pry+sur+guy+guf', 'caca', 'bra+ury', 'soas', 'seas', 'eaas', 'ocea', 'stan', 'mide', 'noaf', 'wssa', 'sssa', 'essa', 'cssa')
regions <- strata
r <- T
raked <- 'raked'

combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                    ages = 0, 
                    regions = regions,
                    holdouts = 0,
                    raked = r,
                    delete_region_files = F)

load("<<<< FILEPATH REDACTED >>>>")

sum_admin_probs <- function(admin_level, target_type, goal_threshold, indicator, indicator_group, run_date, raked) {
  
  admin <- get(paste0('admin_', admin_level))
  
  classify_pixel <- function(x,...) {
    if(target_type == 'less') {
      x[!is.na(x) & x > goal_threshold] <- 0
      x[!is.na(x) & x != 0 & x < goal_threshold] <- 1
    }
    if(target_type == 'greater') {
      x[!is.na(x) & x < goal_threshold] <- 0
      x[!is.na(x) & x != 0 & x > goal_threshold] <- 1
    }
    return(x)
  }
  
  adm_code <- paste0('ADM', admin_level, '_CODE')
  
  ## Admin probabilities
  message(paste0('Calculating ', adm_code, ' probability that ', indicator, ' is ', target_type, ' than ', goal_threshold, '...'))
  admin_probs <- admin[, lapply(.SD, classify_pixel), by=c(adm_code,'year'), .SDcols=grep("^V", names(admin))]
  admin_probs <- admin_probs[, value := apply(.SD, 1, mean), by=c(adm_code,'year'), .SDcols=grep("^V", names(admin_probs))]
  admin_probs <- admin_probs[, c(grep("^ADM", names(admin_probs), value = TRUE),'year','value'), with=FALSE]
  write.csv(admin_probs, paste0("<<<< FILEPATH REDACTED >>>>"))

  
  ## Admin counts (if prop)
  if(grepl('prop', indicator)) {
    message(paste0('Calculating ', adm_code, ' counts...'))
    admin_summary <- admin[, lapply(.SD,  function(x) x*pop), by=c(adm_code,'year'), .SDcols=grep("^V", names(admin))]
    admin_summary <- admin_summary[, lower := apply(.SD, 1, quantile, c(.025), na.rm=T), .SDcols=grep("^V", names(admin_summary))]
    admin_summary <- admin_summary[, value := apply(.SD, 1, mean), .SDcols=grep("^V", names(admin_summary))]
    admin_summary <- admin_summary[, upper := apply(.SD, 1, quantile, c(.975), na.rm=T), .SDcols=grep("^V", names(admin_summary))]
    admin_summary <- admin_summary[, c(grep("^ADM", names(admin_summary), value = TRUE),'year','value','upper','lower'), with=FALSE]
    write.csv(admin_summary, "<<<< FILEPATH REDACTED >>>>")
  }

  # ## Admin summaries
  message(paste0('Calculating ', adm_code, ' mean/upper/lower...'))
  admin_summary <- admin[, lower := apply(.SD, 1, quantile, c(.025), na.rm=T), .SDcols=grep("^V", names(admin))]
  admin_summary <- admin_summary[, value := apply(.SD, 1, mean), .SDcols=grep("^V", names(admin_summary))]
  admin_summary <- admin_summary[, upper := apply(.SD, 1, quantile, c(.975), na.rm=T), .SDcols=grep("^V", names(admin_summary))]
  admin_summary <- admin_summary[, c(grep("^ADM", names(admin_summary), value = TRUE),'year','value','upper','lower'), with=FALSE]
  write.csv(admin_summary, "<<<< FILEPATH REDACTED >>>>")
  
}

lapply(c(0,1, 2), sum_admin_probs,
       goal_threshold = 1,
       target_type = 'greater',
       indicator_group = indicator_group,
       indicator = indicator,
       run_date = run_date,
       raked = raked)
