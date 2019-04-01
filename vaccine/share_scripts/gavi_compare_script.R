###############################################################################
###############################################################################
## Gavi FCE comparison script
##
## Purpose: Generate comparison plots for MBG vs Gavi FCE project estimates
###############################################################################
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<<<< FILEPATH REDACTED >>>>'
indic_repo         <- '<<<< FILEPATH REDACTED >>>>'

## sort some directory stuff
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific code begins here #######################################################

# Set run_date option
run_date <- NULL # Model run date

# Output location
output_dir <- paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/",
                   run_date, "/gavi_comparisons/")

dir.create(output_dir)

### Load all GAVI data
uga <- "<<<< FILEPATH REDACTED >>>>/dpt3_cov_binomial_car_spline_source_preds.rdata"
zmb <- "<<<< FILEPATH REDACTED >>>>/dpt3_cov_binomial_car_spline_source_preds.rdata"
tcd <- "<<<< FILEPATH REDACTED >>>>/dpt3_binomial_car_spline_source_preds.rdata"
cmr <- "<<<< FILEPATH REDACTED >>>>/dpt3_cov_binomial_car_spline_preds.rdata"
moz <- "<<<< FILEPATH REDACTED >>>>/dpt3_cov_binomial_car_spline_source_preds.rdata"

##load and rename GAVI workspace
load(uga)
uga_gavi <- preds
rm(preds)

load(zmb)
zmb_gavi <- preds
rm(preds)

load(tcd)
tcd_gavi <- preds
rm(preds)

load(cmr)
cmr_gavi <- preds
rm(preds)

load(moz)
moz_gavi <- preds
rm(preds)


# load mbg data 
mbg_data_admin1 <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/",run_date,"/pred_derivatives/admin_summaries/dpt3_cov_admin_1_raked_summary.csv"))
mbg_data_admin2 <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/",run_date,"/pred_derivatives/admin_summaries/dpt3_cov_admin_2_raked_summary.csv"))

uga_mbg <- mbg_data_admin1[ADM0_NAME == "Uganda",]
zmb_mbg <- mbg_data_admin2[ADM0_NAME == "Zambia",]
tcd_mbg <- mbg_data_admin1[ADM0_NAME == "Chad",]
cmr_mbg <- mbg_data_admin2[ADM0_NAME == "Cameroon",]
moz_mbg <- mbg_data_admin2[ADM0_NAME == "Mozambique",]

# Rename GAVI tables
rename_gavi <- function(dt){
  names <- names(dt)
  names[[4]] <- "year"
  names[[7]] <- "gavi_mean"
  names(dt) <- names
  return(dt)
}

uga_gavi <- rename_gavi(uga_gavi)
zmb_gavi <- rename_gavi(zmb_gavi)
tcd_gavi <- rename_gavi(tcd_gavi)
cmr_gavi <- rename_gavi(cmr_gavi)
moz_gavi <- rename_gavi(moz_gavi)

compare_uga <- merge(uga_mbg, uga_gavi[level=="dist112",], by.x = c("ADM1_NAME", "year"), by.y=c("area_name","year"))
compare_zmb <- merge(zmb_mbg, zmb_gavi[level=="dist72",], by.x = c("ADM2_NAME", "year"), by.y=c("area_name","year"))
compare_tcd <- merge(tcd_mbg, tcd_gavi[level=="region15",], by.x = c("ADM1_NAME", "year"), by.y=c("area_name","year"))
compare_cmr <- merge(cmr_mbg, cmr_gavi[level=="dept58",], by.x = c("ADM2_NAME", "year"), by.y=c("area_name","year"))
compare_moz <- merge(moz_mbg, moz_gavi[level=="dist148",], by.x = c("ADM2_NAME", "year"), by.y=c("area_name","year"))

## Add variable to color points that have uncertainty bounds that fall outside the 1:1
compare_uga[,r:=1]
compare_uga[((gavi_mean - lower) < 0) | ((gavi_mean - upper) > 0),r:=0]

compare_zmb[,r:=1]
compare_zmb[((gavi_mean - lower) < 0) | ((gavi_mean - upper) > 0),r:=0]

compare_tcd[,r:=1]
compare_tcd[((gavi_mean - lower) < 0) | ((gavi_mean - upper) > 0),r:=0]

compare_cmr[,r:=1]
compare_cmr[((gavi_mean - lower) < 0) | ((gavi_mean - upper) > 0),r:=0]

compare_moz[,r:=1]
compare_moz[((gavi_mean - lower) < 0) | ((gavi_mean - upper) > 0),r:=0]

#for time trends
compare_uga[,z:=1]
compare_uga[((ub-lower) < 0) | ((lb - upper) > 0),z:=0]

compare_zmb[,z:=1]
compare_zmb[((ub - lower) < 0) | ((lb - upper) > 0),z:=0]

compare_tcd[,z:=1]
compare_tcd[((ub-lower) < 0) | ((lb - upper) > 0),z:=0]

compare_cmr[,z:=1]
compare_cmr[((ub-lower) < 0) | ((lb - upper) > 0),z:=0]

compare_moz[,z:=1]
compare_moz[((ub-lower) < 0) | ((lb - upper) > 0),z:=0]

## Uganda
# Gavi vs MBG
out_loc <- paste0(output_dir,"UGA.png")
png(out_loc, res = 600, height = 5, width = 7, units = "in")
ggplot(compare_uga,aes(x=gavi_mean,y=mean)) + 
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0, 1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.01, size=.05) +
  geom_point(alpha=0.2, color = "black")+
  geom_abline(colour="red") + 
  ggtitle("UGA MBG vs GAVI") + 
  labs(x = "GAVI Admin 1 mean DPT3 Coverage",y = "MBG Admin 1 mean DPT3 Coverage") +
  theme_bw()
dev.off()

## Zambia
# Gavi vs MBG
out_loc <- paste0(output_dir,"ZMB.png")
png(out_loc, res = 600, height = 5, width = 7, units = "in")
ggplot(compare_zmb,aes(x=gavi_mean,y=mean)) + 
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0, 1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.01, size=.05) +
  geom_point(alpha=0.2, color = "black") + 
  geom_abline(colour="red") + 
  ggtitle("ZMB MBG vs GAVI") + 
  labs(x = "GAVI Admin 2 mean DPT3 Coverage",y = "MBG Admin 2 mean DPT3 Coverage") +
  theme_bw()
dev.off()

## Chad  
# Gavi vs MBG
out_loc <- paste0(output_dir,"TCD.png")
png(out_loc, res = 600, height = 5, width = 7, units = "in")
ggplot(compare_tcd,aes(x=gavi_mean,y=mean)) + 
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0, 1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.01, size=.1) +
  geom_point(alpha=0.2, color = "black") + 
  geom_abline(colour="red") + 
  ggtitle("TCD MBG vs GAVI") + 
  labs(x = "GAVI Admin 1 mean DPT3 Coverage",y = "MBG Admin 1 mean DPT3 Coverage") +
  theme_bw()
dev.off()

## Cameroon
# Gavi vs MBG
out_loc <- paste0(output_dir,"CMR.png")
png(out_loc, res = 600, height = 5, width = 7, units = "in")
ggplot(compare_cmr,aes(x=gavi_mean,y=mean)) + 
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0, 1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.01, size=.1) +
  geom_point(alpha=0.2, color = "black") + 
  geom_abline(colour="red") + 
  ggtitle("CMR MBG vs GAVI") + 
  labs(x = "GAVI Admin 2 mean DPT3 Coverage",y = "MBG Admin 2 mean DPT3 Coverage") +
  theme_bw()
dev.off()

## Mozambique
out_loc <- paste0(output_dir,"MOZ.png")
png(out_loc, res = 600, height = 5, width = 7, units = "in")
ggplot(compare_moz,aes(x=gavi_mean,y=mean)) + 
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0, 1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.01, size=.05) +
  geom_point(alpha=0.2, color = "black") + 
  geom_abline(colour="red") + 
  ggtitle("MOZ MBG vs GAVI") + 
  labs(x = "GAVI Admin 2 mean DPT3 Coverage",y = "MBG Admin 2 mean DPT3 Coverage") +
  theme_bw()
dev.off()
