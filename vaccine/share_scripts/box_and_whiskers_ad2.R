###############################################################################
###############################################################################
## Box and whisker plotting code
##
## Purpose: Produce box-and-whisker plot of admin 2 estimates by country
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

library("beeswarm")
library("ggbeeswarm")

rd <- NULL # Change this to reflect the run date
ind <- "dpt3_cov"
ig <- "vaccine"
ad_level <- 2
raked <- T
yr <- 2016

in_dir <- paste0("<<<< FILEPATH REDACTED >>>>/",ig,"/", ind, "/output/", rd, "/pred_derivatives/admin_summaries/")
out_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", rd, "/number_plugging/")

ad2_file <- paste0(in_dir, ind, "_admin_", ad_level, ifelse(raked, "_raked_", "_unraked_"), "summary.csv")
ad2_df <- fread(ad2_file)

ad0_file <- paste0(in_dir, ind, "_admin_0", ifelse(raked, "_raked_", "_unraked_"), "summary.csv")
ad0_df <- fread(ad0_file)

# Subset to year of interest
ad2_df <- subset(ad2_df, year == yr)
ad0_df <- subset(ad0_df, year == yr)

# Drop Ma'tan al-Sarra
ad2_df <- subset(ad2_df, ADM0_CODE != 40762)
ad0_df <- subset(ad0_df, ADM0_CODE != 40762)

# Truncate long names
ad0_df[ADM0_NAME == "Central African Republic", ADM0_NAME := "CAR"]
ad0_df[ADM0_NAME == "Democratic Republic of the Congo", ADM0_NAME := "DR Congo"]
ad0_df[ADM0_NAME == "Sao Tome and Principe", ADM0_NAME := "STP"]
ad0_df[ADM0_NAME == "United Republic of Tanzania", ADM0_NAME := "Tanzania"]

ad2_df[ADM0_NAME == "Central African Republic", ADM0_NAME := "CAR"]
ad2_df[ADM0_NAME == "Democratic Republic of the Congo", ADM0_NAME := "DR Congo"]
ad2_df[ADM0_NAME == "Sao Tome and Principe", ADM0_NAME := "STP"]
ad2_df[ADM0_NAME == "United Republic of Tanzania", ADM0_NAME := "Tanzania"]

# Sort ad0_df
setorderv(ad0_df, "mean")
ad0_df[, ADM0_RANK := .I]

ad2_df <- merge(ad2_df, subset(ad0_df, select = c("ADM0_CODE", "ADM0_RANK")))

gg_violin <- ggplot(data = ad2_df, aes(x = forcats::fct_reorder(ADM0_NAME, ADM0_RANK), y = mean)) +
  #geom_point(alpha = 0.1) +
  geom_violin(trim = TRUE, scale = "width", fill = "gray90", color = "black") +
  geom_point(data = ad0_df, aes(x = forcats::fct_reorder(ADM0_NAME, ADM0_RANK), y = mean),
             color = "red", shape = 8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,0.1), labels = scales::percent) +
  labs(x = "Country", y = "DPT3 Coverage at the Second Administrative Level (2016)")

gg_box <- ggplot(data = ad2_df, aes(x = forcats::fct_reorder(ADM0_NAME, ADM0_RANK), y = mean)) +
  #geom_point(alpha = 0.1) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,0.1), labels = scales::percent) +
  geom_point(data = ad0_df, aes(x = forcats::fct_reorder(ADM0_NAME, ADM0_RANK), y = mean),
             color = "red", shape = 8) +
  labs(x = "Country", y = "DPT3 Coverage at the Second Administrative Level (2016)")

png(file = paste0(out_dir, ind, "_ad2_boxplot_by_country_2016.png"),
    width = 11,
    height = 5.5,
    units = "in",
    res = 300)
print(gg_box)
dev.off()

pdf(file = paste0(out_dir, ind, "_ad2_boxplot_by_country_2016.pdf"),
    width = 11,
    height = 5.5)
print(gg_box)
dev.off()

png(file = paste0(out_dir, ind, "_ad2_violin_plot_by_country_2016.png"),
    width = 11,
    height = 5.5,
    units = "in",
    res = 300)
print(gg_violin)
dev.off()

pdf(file = paste0(out_dir, ind, "_ad2_violin_plot_by_country_2016.pdf"),
    width = 11,
    height = 5.5)
print(gg_violin)
dev.off()


