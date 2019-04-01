###############################################################################
###############################################################################
## Plot GBD vs MBG national-level estimates
##
## Purpose: Create scatterplots of GBD vs uncalibrated MBG estimates at the
##          national level
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

## Script-specific code begins here ##########################################

rd <- NULL # Change to model run date
ig <- "vaccine"

for (ind in c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout")) {

  # Load admin0 (national level) estimates
  ad0_raked <- fread(paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", rd, "/pred_derivatives/admin_summaries/",
                          ind, "_admin_0_raked_summary.csv"))
  ad0_unraked <- fread(paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", rd, "/pred_derivatives/admin_summaries/",
                          ind, "_admin_0_unraked_summary.csv"))

  ad0_raked <- subset(ad0_raked, select = c("ADM0_CODE", "ADM0_NAME", "year", "mean", "lower", "upper"))
  ad0_unraked <- subset(ad0_unraked, select = c("ADM0_CODE", "ADM0_NAME", "year", "mean", "lower", "upper"))

  setnames(ad0_raked, 
           c("mean", "upper", "lower"), 
           c("mean_raked", "upper_raked", "lower_raked"))

  setnames(ad0_unraked, 
         c("mean", "upper", "lower"), 
         c("mean_unraked", "upper_unraked", "lower_unraked"))

  ad0_all <- merge(ad0_raked, ad0_unraked)
  ad0_all[, diff := abs(mean_raked - mean_unraked)]
  ad0_all[diff > 0.25, lab := paste0(ADM0_NAME, " ", year)]
  ad0_all[is.na(lab), lab := ""]

  ad0_all <- subset(ad0_all, ADM0_NAME != "Ma'tan al-Sarra")

  out_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", rd, "/number_plugging/rf_plots/")
  dir.create(out_dir, showWarnings = F)
  
  if (ind == "dpt3_cov") plot_title <- "DPT3 Coverage"
  if (ind == "dpt1_cov") plot_title <- "DPT1 Coverage"
  if (ind == "dpt1_3_abs_dropout") plot_title <- "DPT1-3 Absolute Dropout"
  
  gg <- ggplot(data = ad0_all, aes(x=mean_raked, y = mean_unraked)) +
    geom_point(alpha = 0.2) +
    geom_errorbar(aes(ymin = lower_unraked, ymax = upper_unraked), alpha = 0.1) +
    theme_bw() +
    geom_abline(slope = 1, color = "red") +
    geom_text_repel(aes(label = lab), size = 1.8, segment.size = 0.2) +
    coord_equal() +
    xlim(0,1) +
    ylim(0,1) +
    labs(x = "GBD estimate", 
         y = "MBG estimate",
         title = plot_title)
  
  png(file = paste0(out_dir, "rf_plot_", ind, ".png"),
      height = 6,
      width = 6,
      units = "in",
      res = 300)
    
  print(gg)
  dev.off()
}

