###############################################################################
###############################################################################
## MBG "Number Plugging" script
##
## Purpose: This script performs calculations required for the DPT manuscript
##          using modeled results, and writes those to .csv and .txt files
###############################################################################
###############################################################################

# Create directories
for(indicator in postest_indicators){
  input_dir =  paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', indicator, '/output/', run_date)
  output_dir =  paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', indicator, '/output/', run_date, "/number_plugging/")
  dir.create(output_dir, showWarnings = FALSE)  
  out_file <- paste0(output_dir, indicator, "_number_plugging.txt")
  if(file.exists(out_file)) unlink(out_file)
}
# Define a convenience function to write to a file
write_np <- function(text, 
                     ind = indicator,
                     ig = indicator_group,
                     rd = run_date) {
  out_file <- paste0('<<<< FILEPATH REDACTED >>>>/', ig, '/', ind, '/output/', rd, "/number_plugging/", ind, "_number_plugging.txt")
  if (!file.exists(out_file)) file.create(out_file)
  fileConn<-file(out_file, open = "a")
  writeLines(text, fileConn)
  close(fileConn)
}

for (indicator in postest_indicators) {
  write_np("##########################################################################", ind = indicator)
  write_np(paste0("## Number plugging for ", indicator, " | run_date: ", run_date), ind = indicator)
  write_np("##########################################################################", ind = indicator)
}

### Calculating from pre-pointpoly processed input data

## Key Points -------------------------------

# Number of countries with >95% probability that all admin2s are >80% in 2016
all_admin_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/dpt3_cov/output/", 
                       run_date, "/pred_derivatives/target_probs/",
                       "dpt3_cov_all_admins_2016_absolute_greater_0.8_admin_target_probs.csv"))

all_admin_df <- subset(all_admin_df, ADM0_CODE != 40762) # exclude Ma'tan al-Sarra
total_countries <- length(unique(all_admin_df$ADM0_NAME))
total_met_goal <- sum(all_admin_df$probability >= 0.95)

write_np(paste0("\nIn 2016, ", total_met_goal, " of ", total_countries, " countries ",
               "had >80% coverage in all admin2 units with high certainty (>95%):",
               paste(unique(all_admin_df[probability>=0.95]$ADM0_NAME), collapse = ",")), ind = "dpt3_cov")

# Number of countries with >95% probability that all admin2s are >80% in 2020
all_admin_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/dpt3_cov/output/", 
                       run_date, "/pred_derivatives/target_probs/",
                       "dpt3_cov_all_admins_2020_absolute_greater_0.8_admin_target_probs.csv"))

all_admin_df <- subset(all_admin_df, ADM0_CODE != 40762) # exclude Ma'tan al-Sarra
total_countries <- length(unique(all_admin_df$ADM0_NAME))
total_met_goal <- sum(all_admin_df$probability >= 0.95)

write_np(paste0("\nIn 2020, ", total_met_goal, " of ", total_countries, " countries ",
               "projected to have >80% coverage in all admin2 units with high certainty (>95%):\n"), ind = "dpt3_cov")

sink(paste0('<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/', run_date, "/number_plugging/dpt3_cov_number_plugging.txt"), append = T)
print(subset(all_admin_df, probability >= 0.95, 
             select = c("ADM0_NAME", "ADM0_CODE", "probability")))
sink()

## Methods ----------------------------------

# Number of Covariates
config <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/", run_date, "/config.csv"))

covs <- config[56,2]
covs <- strsplit(as.character(covs[[1]])," + ", fixed=TRUE)[[1]]
GBD_covs <- config[58,2]
GBD_covs <- strsplit(as.character(GBD_covs[[1]])," + ", fixed=TRUE)[[1]]
total_covs <- c(covs, GBD_covs)
number_covs <- uniqueN(total_covs)

write_np(paste0("Number of unique covariates used: ", number_covs), ind="dpt3_cov")
write_np(paste0("Number of unique covariates used: ", number_covs), ind="dpt1_cov")

## Results -----------------------------------

#########################################################################################
### DPT3 ################################################################################
#########################################################################################

# Load data sets and investigate differences
admin2_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt3_cov_admin_2_raked_summary.csv"))
admin2_diff_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt3_cov_admin_2_raked_diff_2000-2016.csv"))
first_year <- min(unique(admin2_df$year)) 
last_year <- max(unique(admin2_df$year))

# Look at biggest absolute differences subnationally
admin2_last_year <- subset(admin2_df, year == last_year)
admin2_last_year[, max := max(mean, na.rm = T), by = ADM0_CODE]
admin2_last_year[, min := min(mean, na.rm = T), by = ADM0_CODE]
ad2_maxmin <- subset(admin2_last_year, max == mean | min == mean)
ad2_maxmin[, diff := max - min]

write_np("Biggest difference subnationally:", ind = "dpt3_cov")
biggest_difference <- subset(ad2_maxmin, diff == max(ad2_maxmin$diff, na.rm = T))
sink(paste0('<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/', run_date, "/number_plugging/dpt3_cov_number_plugging.txt"), append = T)
print(biggest_difference)
sink()

fwrite(ad2_maxmin, file = paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/", run_date, "/number_plugging/dpt3_ad2_maxmin_2016.csv"))

# Additional custom analysis for NGA, ETH, DRC -- populous countries not meeting GVAP goals 

for (country in c("Nigeria", "Ethiopia", "Democratic Republic of the Congo")) {
  
  # admin 2 means
  first_year_adm2_df <- subset(admin2_df, ADM0_NAME == country & year == first_year)
  last_year_adm2_df <- subset(admin2_df, ADM0_NAME == country & year == last_year)

  # admin 0 mean
  admin0_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt3_cov_admin_0_raked_summary.csv"))
  first_year_df <- subset(admin0_df, ADM0_NAME == country & year == first_year)
  last_year_df <- subset(admin0_df, ADM0_NAME == country & year == last_year)

  write_np(paste0("Average admin0 coverage in ", country, " changed from ", 
                 round(first_year_df$mean * 100, 1), "% in ", first_year, " to ",
                 round(last_year_df$mean * 100, 1), "% in ", last_year, "."), ind = "dpt3_cov")

  ad2_diff_country <- subset(admin2_diff_df, ADM0_NAME == country)
  ad2_diff_country <- subset(ad2_diff_country, (complete.cases(ad2_diff_country)))
  best_change <- as.numeric(subset(ad2_diff_country, mean == max(mean), select = ADM2_CODE))
  worst_change <- as.numeric(subset(ad2_diff_country, mean == min(mean),select = ADM2_CODE))

  best_change_df <- rbind(first_year_adm2_df[ADM2_CODE == best_change],
                          last_year_adm2_df[ADM2_CODE == best_change])

  worst_change_df <- rbind(first_year_adm2_df[ADM2_CODE == worst_change],
                           last_year_adm2_df[ADM2_CODE == worst_change])

  write_np(paste0("Best change in ", country, " seen in ", unique(best_change_df$ADM2_NAME),
                 " where coverage changed from ", round(best_change_df[year == first_year]$mean * 100, 1), 
                 "% (95% UI: ", round(best_change_df[year == first_year]$lower * 100, 1), "% - ", 
                 round(best_change_df[year == first_year]$upper * 100, 1), "%)", 
                 " in ", first_year, " to ", round(best_change_df[year == last_year]$mean * 100, 1),
                 "% (95% UI: ", round(best_change_df[year == last_year]$lower * 100, 1), "& - ", 
                 round(best_change_df[year == last_year]$upper * 100, 1), "%)", 
                 " in ", last_year, "."), ind = "dpt3_cov")

  write_np(paste0("Worst change in ", country, " seen in ", unique(worst_change_df$ADM2_NAME),
               " where coverage changed from ", round(worst_change_df[year == first_year]$mean * 100, 1), 
               "% (95% UI: ", round(worst_change_df[year == first_year]$lower * 100, 1), "% - ", 
               round(worst_change_df[year == first_year]$upper * 100, 1), "%)", 
               " in ", first_year, " to ", round(worst_change_df[year == last_year]$mean * 100, 1),
               "% (95% UI: ", round(worst_change_df[year == last_year]$lower * 100, 1), "& - ", 
               round(worst_change_df[year == last_year]$upper * 100, 1), "%)", 
               " in ", last_year, "."), ind = "dpt3_cov")

}

# Target probabilities for NGA, ETH, DRC
prob_above_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/",run_date,"/pred_derivatives/admin_summaries/dpt3_cov_admin_2_raked_p_0.8_or_better_summary.csv"))

for (country in c("Nigeria", "Ethiopia", "Democratic Republic of the Congo")) {
  country_above_df <- subset(prob_above_df, ADM0_NAME == country & year == last_year)
  print(country_above_df[p_above >= 0.95])
  print(country_above_df[p_above <= 0.05])
  message(country)
  high_prob <- nrow(country_above_df[p_above >= 0.95])
  low_prob <- nrow(country_above_df[p_above <= 0.05])
  write_np(paste0("\n", country), ind = "dpt3_cov")
  write_np(paste0("High > 95% prob coverage > 80%: ", high_prob, " out of ", nrow(country_above_df), " (", round((high_prob / nrow(country_above_df))*100, 1), "%)"), ind = "dpt3_cov")
  write_np(paste0("Low < 5% prob coverage > 80%: ", low_prob, " out of ", nrow(country_above_df), " (", round((low_prob / nrow(country_above_df))*100, 1), "%)"), ind = "dpt3_cov")
}


# Look at variation compared to national averages
admin2_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt3_cov_admin_2_raked_summary.csv"))
admin2_df <- subset(admin2_df, year == max(admin2_df$year))

admin0_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt3_cov_admin_0_raked_summary.csv"))
admin0_df <- subset(admin0_df, year == max(admin0_df$year))
setnames(admin0_df, "mean", "nat_mean")
admin0_df <- subset(admin0_df, select = c("ADM0_CODE", "year", "nat_mean"))

admin2_df <- merge(admin2_df, admin0_df)
admin2_df[, ratio := mean / nat_mean]

ratio_df <- lapply(unique(admin2_df$ADM0_CODE), function (ad0_code) {
  c_df <- admin2_df[ADM0_CODE == ad0_code]
  return(data.table(ADM0_NAME = unique(c_df$ADM0_NAME),
                    ADM0_CODE = unique(c_df$ADM0_CODE),
                    min_ratio = min(c_df$ratio, na.rm = T),
                    max_ratio = max(c_df$ratio, na.rm = T)))
})

ratio_df <- rbindlist(ratio_df)
ratio_df[order(min_ratio)]
ratio_df <- subset(ratio_df, ADM0_NAME != "Ma'tan al-Sarra")

numer <- nrow(ratio_df[min_ratio < 0.75])
denom <- nrow(ratio_df)

write_np(paste0("For DPT3, ", numer, " out of ", denom, " (", round((numer/denom)*100, 2), "%)",
                "countries have at least one ad2 < 75% of the national average."), ind = "dpt3_cov")


numer <- nrow(ratio_df[min_ratio < 0.5])
denom <- nrow(ratio_df)

write_np(paste0("For DPT3, ", numer, " out of ", denom, " (", round((numer/denom)*100, 2), "%)",
                "countries have at least one ad2 < 50% of the national average."), ind = "dpt3_cov")

ratio_df[, diff_ratio := max_ratio-min_ratio]
ratio_df <- ratio_df[rev(order(diff_ratio))]
fwrite(ratio_df, file = paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt3_cov/output/", 
                               run_date, "/number_plugging/dpt3_cov_ad2_ratios_by_country.csv"))

#########################################################################################
### DPT1 ################################################################################
#########################################################################################

# Load data sets and investigate differences
admin2_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt1_cov_admin_2_raked_summary.csv"))
admin2_diff_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt1_cov_admin_2_raked_diff_2000-2016.csv"))
first_year <- min(unique(admin2_df$year)) 
last_year <- max(unique(admin2_df$year))

# Look at biggest absolute differences subnationally
admin2_last_year <- subset(admin2_df, year == last_year)
admin2_last_year[, max := max(mean, na.rm = T), by = ADM0_CODE]
admin2_last_year[, min := min(mean, na.rm = T), by = ADM0_CODE]
ad2_maxmin <- subset(admin2_last_year, max == mean | min == mean)
ad2_maxmin[, diff := max - min]

biggest_difference <- subset(ad2_maxmin, diff == max(ad2_maxmin$diff, na.rm = T))
sink(paste0('<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/', run_date, "/number_plugging/dpt1_cov_number_plugging.txt"), append = T)
print(biggest_difference)
sink()

fwrite(ad2_maxmin, file = paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/", run_date, "/number_plugging/dpt1_ad2_maxmin_2016.csv"))

# Additional custom analysis for NGA, ETH, DRC -- populous countries not meeting GVAP goals 
for (country in c("Nigeria", "Ethiopia", "Democratic Republic of the Congo")) {
  ##adm2
  first_year_adm2_df <- subset(admin2_df, ADM0_NAME == country & year == first_year)
  last_year_adm2_df <- subset(admin2_df, ADM0_NAME == country & year == last_year)

  #admin 0 mean
  admin0_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt1_cov_admin_0_raked_summary.csv"))
  first_year_df <- subset(admin0_df, ADM0_NAME == country & year == first_year)
  last_year_df <- subset(admin0_df, ADM0_NAME == country & year == last_year)

  write_np(paste0("\n", country), ind = "dpt1_cov")

  write_np(paste0("Average admin0 coverage in ", country, " changed from ", 
                 round(first_year_df$mean * 100, 1), "% in ", first_year, " to ",
                 round(last_year_df$mean * 100, 1), "% in ", last_year, "."), ind = "dpt1_cov")

  ad2_diff_country <- subset(admin2_diff_df, ADM0_NAME == country)
  ad2_diff_country <- subset(ad2_diff_country, (complete.cases(ad2_diff_country)))
  best_change <- as.numeric(subset(ad2_diff_country, mean == max(mean), select = ADM2_CODE))
  worst_change <- as.numeric(subset(ad2_diff_country, mean == min(mean),select = ADM2_CODE))

  best_change_df <- rbind(first_year_adm2_df[ADM2_CODE == best_change],
                          last_year_adm2_df[ADM2_CODE == best_change])

  worst_change_df <- rbind(first_year_adm2_df[ADM2_CODE == worst_change],
                           last_year_adm2_df[ADM2_CODE == worst_change])

  write_np(paste0("Best change in ", country, " seen in ", unique(best_change_df$ADM2_NAME),
                 " where coverage changed from ", round(best_change_df[year == first_year]$mean * 100, 1), 
                 "% (95% UI: ", round(best_change_df[year == first_year]$lower * 100, 1), "% - ", 
                 round(best_change_df[year == first_year]$upper * 100, 1), "%)", 
                 " in ", first_year, " to ", round(best_change_df[year == last_year]$mean * 100, 1),
                 "% (95% UI: ", round(best_change_df[year == last_year]$lower * 100, 1), "& - ", 
                 round(best_change_df[year == last_year]$upper * 100, 1), "%)", 
                 " in ", last_year, "."), ind = "dpt1_cov")

  write_np(paste0("Worst change in ", country, " seen in ", unique(worst_change_df$ADM2_NAME),
               " where coverage changed from ", round(worst_change_df[year == first_year]$mean * 100, 1), 
               "% (95% UI: ", round(worst_change_df[year == first_year]$lower * 100, 1), "% - ", 
               round(worst_change_df[year == first_year]$upper * 100, 1), "%)", 
               " in ", first_year, " to ", round(worst_change_df[year == last_year]$mean * 100, 1),
               "% (95% UI: ", round(worst_change_df[year == last_year]$lower * 100, 1), "& - ", 
               round(worst_change_df[year == last_year]$upper * 100, 1), "%)", 
               " in ", last_year, "."), ind = "dpt1_cov")

}

# Target probabilities for NGA, ETH, DRC
prob_above_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/",run_date,"/pred_derivatives/admin_summaries/dpt1_cov_admin_2_raked_p_0.8_or_better_summary.csv"))

for (country in c("Nigeria", "Ethiopia", "Democratic Republic of the Congo")) {
  country_above_df <- subset(prob_above_df, ADM0_NAME == country & year == last_year)
  print(country_above_df[p_above >= 0.95])
  print(country_above_df[p_above <= 0.05])
  message(country)
  high_prob <- nrow(country_above_df[p_above >= 0.95])
  low_prob <- nrow(country_above_df[p_above <= 0.05])
  write_np(paste0("\n", country), ind = "dpt1_cov")
  write_np(paste0("High > 95% prob coverage > 80%: ", high_prob, " out of ", nrow(country_above_df), " (", round((high_prob / nrow(country_above_df))*100, 1), "%)"), ind = "dpt1_cov")
  write_np(paste0("Low < 5% prob coverage > 80%: ", low_prob, " out of ", nrow(country_above_df), " (", round((low_prob / nrow(country_above_df))*100, 1), "%)"), ind = "dpt1_cov")
}

# Look at variation compared to national averages
admin2_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt1_cov_admin_2_raked_summary.csv"))
admin2_df <- subset(admin2_df, year == max(admin2_df$year))

admin0_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/", run_date, "/", 
                          "pred_derivatives/admin_summaries/dpt1_cov_admin_0_raked_summary.csv"))
admin0_df <- subset(admin0_df, year == max(admin0_df$year))
setnames(admin0_df, "mean", "nat_mean")
admin0_df <- subset(admin0_df, select = c("ADM0_CODE", "year", "nat_mean"))

admin2_df <- merge(admin2_df, admin0_df)
admin2_df[, ratio := mean / nat_mean]

ratio_df <- lapply(unique(admin2_df$ADM0_CODE), function (ad0_code) {
  c_df <- admin2_df[ADM0_CODE == ad0_code]
  return(data.table(ADM0_NAME = unique(c_df$ADM0_NAME),
                    ADM0_CODE = unique(c_df$ADM0_CODE),
                    min_ratio = min(c_df$ratio, na.rm = T),
                    max_ratio = max(c_df$ratio, na.rm = T)))
})

ratio_df <- rbindlist(ratio_df)
ratio_df[order(min_ratio)]
ratio_df <- subset(ratio_df, ADM0_NAME != "Ma'tan al-Sarra")

numer <- nrow(ratio_df[min_ratio < 0.75])
denom <- nrow(ratio_df)

write_np(paste0("For DPT1, ", numer, " out of ", denom, " (", round((numer/denom)*100, 2), "%)",
                "countries have at least one ad2 < 75% of the national average."), ind = "dpt1_cov")

ratio_df[, diff_ratio := max_ratio-min_ratio]
ratio_df <- ratio_df[rev(order(diff_ratio))]
fwrite(ratio_df, file = paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_cov/output/", 
                               run_date, "/number_plugging/dpt1_cov_ad2_ratios_by_country.csv"))

#########################################################################################
### RELATIVE DROPOUT ####################################################################
#########################################################################################

# Load data sets and investigate differences
admin2_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_3_rel_dropout/output/", run_date, "/", 
                        "pred_derivatives/admin_summaries/dpt1_3_rel_dropout_admin_2_raked_summary.csv"))
admin2_diff_df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_3_rel_dropout/output/", run_date, "/", 
                        "pred_derivatives/admin_summaries/dpt1_3_rel_dropout_admin_2_raked_diff_2000-2016.csv"))
first_year <- min(unique(admin2_df$year)) 
last_year <- max(unique(admin2_df$year))

# Look at biggest absolute differences subnationally
admin2_last_year <- subset(admin2_df, year == last_year)
admin2_last_year[, max := max(mean, na.rm = T), by = ADM0_CODE]
admin2_last_year[, min := min(mean, na.rm = T), by = ADM0_CODE]
ad2_maxmin <- subset(admin2_last_year, max == mean | min == mean)
ad2_maxmin[, diff := max - min]

biggest_difference <- subset(ad2_maxmin, diff == max(ad2_maxmin$diff, na.rm = T))
sink(paste0('<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_3_rel_dropout/output/', run_date, "/number_plugging/dpt1_3_rel_dropout_number_plugging.txt"), append = T)
print(biggest_difference)
sink()

fwrite(ad2_maxmin, file = paste0("<<<< FILEPATH REDACTED >>>>/vaccine/dpt1_3_rel_dropout/output/", run_date, "/number_plugging/dpt1_3_rel_dropout_ad2_maxmin_2016.csv"))

# Look at differences over the study period
first_year_df <- subset(admin2_df, year == first_year)
last_year_df <- subset(admin2_df, year == last_year)

first_year_mean <- mean(first_year_df$mean, na.rm = T)
last_year_mean <- mean(last_year_df$mean, na.rm = T)

ad2_diff_country <- subset(admin2_diff_df)
ad2_diff_country <- subset(ad2_diff_country, (complete.cases(ad2_diff_country)))
best_change <- as.numeric(subset(ad2_diff_country, mean == max(mean), select = ADM2_CODE))
worst_change <- as.numeric(subset(ad2_diff_country, mean == min(mean),select = ADM2_CODE))

best_change_df <- rbind(first_year_df[ADM2_CODE == best_change],
                        last_year_df[ADM2_CODE == best_change])

worst_change_df <- rbind(first_year_df[ADM2_CODE == worst_change],
                         last_year_df[ADM2_CODE == worst_change])

write_np(paste0("Best change in relative dropout seen in ", unique(best_change_df$ADM2_NAME),
               " where coverage changed from ", round(best_change_df[year == first_year]$mean * 100, 1), 
               "% (95% UI: ", round(best_change_df[year == first_year]$lower * 100, 1), "% - ", 
               round(best_change_df[year == first_year]$upper * 100, 1), "%)", 
               " in ", first_year, " to ", round(best_change_df[year == last_year]$mean * 100, 1),
               "% (95% UI: ", round(best_change_df[year == last_year]$lower * 100, 1), "& - ", 
               round(best_change_df[year == last_year]$upper * 100, 1), "%)", 
               " in ", last_year, "."), ind = "dpt1_3_rel_dropout")

write_np(paste0("Worst change in relative dropout seen in ", unique(worst_change_df$ADM2_NAME),
             " where coverage changed from ", round(worst_change_df[year == first_year]$mean * 100, 1), 
             "% (95% UI: ", round(worst_change_df[year == first_year]$lower * 100, 1), "% - ", 
             round(worst_change_df[year == first_year]$upper * 100, 1), "%)", 
             " in ", first_year, " to ", round(worst_change_df[year == last_year]$mean * 100, 1),
             "% (95% UI: ", round(worst_change_df[year == last_year]$lower * 100, 1), "& - ", 
             round(worst_change_df[year == last_year]$upper * 100, 1), "%)", 
             " in ", last_year, "."), ind = "dpt1_3_rel_dropout")

#########################################################################################
### ADDITIONAL ADMIN-LEVEL ANALYSES #####################################################
#########################################################################################

## set admin levels to loop over
sum_metric_regions <-c( "ad1","ad2","country")
admin_levels <- c("admin_0","admin_1","admin_2")

for(indicator in postest_indicators){
  message("indicator: ", indicator)
  input_dir =  paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', indicator, '/output/', run_date)
  output_dir =  paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', indicator, '/output/', run_date, "/number_plugging/")
  
    ### Create new simplified tables for summary metrics for admin0, 1 and 2
    for (reg in sum_metric_regions){
      if (file.exists(paste0(input_dir, "/summary_metrics/",reg,"_metrics.csv"))) {

      sum_metric <- fread(paste0(input_dir, "/summary_metrics/",reg,"_metrics.csv"))
      ## Subset to columns of interest
      cleaned_sum <- sum_metric[,c("Year","Mean Err.","RMSE","Corr.","95% Cov.")]
      ## Reassign names for final table if desired
      names(cleaned_sum) <- c("Year","Mean Error","RMSE","Corr.","95% Cov.")

      write.csv(cleaned_sum, paste0(output_dir,reg,"_summary.csv"))
    
    }
  }

  
  ### Identify 10 highest and lowest regions for admin 0, 1 and 2
  for (ad_level in admin_levels){
    i <- 10
    #reg <- "admin_1"
    summary <- fread(paste0(input_dir, "/pred_derivatives/admin_summaries/",indicator,"_",ad_level,"_raked_summary.csv"))
    ## set to most recent year
    top_year <- max(unique(summary$year))
    #set year to numeric
    summary[, year := as.numeric(year)]
    
    year_summary<- summary[year==top_year & !is.na(mean),]
    year_summary <- year_summary[order(-mean)]
    top <- head(year_summary,i)
    bottom <- tail(year_summary,i)
    
    #export to csv
    write.csv(top, paste0(output_dir,ad_level,"_",top_year,"_highest.csv"))
    write.csv(bottom, paste0(output_dir,ad_level,"_",top_year,"_lowest.csv"))
  }
  
  ### Identify greatest change in mean vaccine coverage from 2000-2016
  write_np("\nIncreases and decreases in coverage:", ind = indicator)
  summary <- fread(paste0(input_dir, "/pred_derivatives/admin_summaries/",indicator,"_",ad_level,"_raked_diff_2000-2016.csv"))

  # Check if CIs overlap
  summary[,ci_overlap := ifelse((mean > 0 & lower > 0) | (mean < 0 & upper < 0), F, T)]
    
  total_ad2s <- length(unique(summary$ADM2_CODE))

  # report out the numbers
  sig_decline_table <- summary[mean < 0 & ci_overlap == F,]
  sig_decline_table <- sig_decline_table[order(mean)]
  decline_ad2s <- length(unique(sig_decline_table$ADM2_CODE))
  write_np(paste0(decline_ad2s, " out of ", total_ad2s, " (",
                 round ((decline_ad2s/total_ad2s)*100,1), 
                 "%) admin 2s had decreasing coverage with P>95%"),
           ind = indicator)
  write.csv(sig_decline_table, paste0(output_dir,ad_level,"_", unique(summary$year),"_decline.csv"))

  sig_increase_table <- summary[mean > 0 & ci_overlap == F,]
  sig_increase_table <- sig_increase_table[order(-mean)]
  increase_ad2s <- length(unique(sig_increase_table$ADM2_CODE))
  write_np(paste0(increase_ad2s, " out of ", total_ad2s, " (",
                 round ((increase_ad2s/total_ad2s)*100,1), 
                 "%) admin 2s had increasing coverage with P>95%"),
           ind = indicator)
  write.csv(sig_increase_table, paste0(output_dir,ad_level,"_",unique(summary$year),"_increase.csv"))
  
  #### Code to generate file for countries with at least 1 admin 2 unit below 80% coverage in 2016
  summary <- fread(paste0(input_dir, "/pred_derivatives/admin_summaries/",indicator,"_admin_2_raked_summary.csv"))
  top_year <- max(unique(summary$year))
  #convert year to numeric
  summary[,year:=as.numeric(year)]
  
  recent_summary <- summary[year==top_year, ]
  recent_summary[,num_adm2 := .N, by=.(ADM0_NAME)]
  
  ## Admin2's below 80% w/o confidence interval
  below_80_no_CI <- recent_summary[mean < .80,]
  ## Admin2's below 80% w/ confidence interval
  below_80_w_CI <- recent_summary[upper < .80,]
  
  ## Identify total number of admin2's below 80% and calculate % of total admin2s per country
  below_80_no_CI <- below_80_no_CI[,.N, by=.(ADM0_NAME,  ADM0_CODE, year,num_adm2)]
  names(below_80_no_CI) <- c("ADM0_NAME",   "ADM0_CODE","year","total_num_ad2", "num_ad2_below_80")
  below_80_no_CI[,percent_below_80 := num_ad2_below_80 / total_num_ad2]
  
  below_80_w_CI <- below_80_w_CI[,.N, by=.(ADM0_NAME, ADM0_CODE,  year,num_adm2)]
  names(below_80_w_CI) <- c("ADM0_NAME",   "ADM0_CODE", "year","total_num_ad2","num_ad2_below_80")
  below_80_w_CI[,percent_below_80 := num_ad2_below_80 / total_num_ad2]
  
  write.csv(below_80_w_CI, paste0(output_dir,"ad0_w_ad2_below_80_including_CI_",top_year,".csv"))
  write.csv(below_80_no_CI, paste0(output_dir,"ad0_w_ad2_below_80_without_CI_",top_year,".csv"))
  
  labels <- c("admin0s w/ admin2s below 80% without confidence interval","admin0s w/ admin2s below 80% including confidence interval")
  vals <- c(nrow(below_80_no_CI),nrow(below_80_w_CI))
  output_summary <- data.frame(labels, vals)
  write.csv(output_summary, paste0(output_dir,"ad0_w_ad2_below_80_count_",top_year,".csv"))

  # Look at changes in ad2s draw-wise
  load(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/", indicator, "/output/", run_date, "/", indicator, "_raked_admin_draws_eb_bin0_0.RData"))
  ad2_first_year_draws <- admin_2[year == min(unique(admin_2$year))]
  ad2_last_year_draws <- admin_2[year == max(unique(admin_2$year))]

  # Check to make sure all aligned
  if (all.equal(ad2_first_year_draws$ADM2_CODE, ad2_last_year_draws$ADM2_CODE)) {
    draw_cols <- names(ad2_first_year_draws)[grepl("V[0-9]+", names(ad2_first_year_draws))]
    ad2_diff <- as.matrix(subset(ad2_last_year_draws, select = draw_cols)) - as.matrix(subset(ad2_first_year_draws, select = draw_cols))
    ad2_increases <- colSums(ad2_diff > 0, na.rm = T) / colSums(!is.na(ad2_diff))
    ad2_decreases <- colSums(ad2_diff < 0, na.rm = T) / colSums(!is.na(ad2_diff))

    write_np(paste0("Estimated ", indicator, " increased in ", round(mean(ad2_increases)*100,1), 
                    "% (", round(quantile(ad2_increases, 0.025)*100, 1), " - ", 
                    round(quantile(ad2_increases, 0.975)*100, 1), 
                    "%) of second-level administrative units in Africa from 2000-2016"),
             ind = indicator)

    write_np(paste0("Estimated", indicator, " decreased in ", round(mean(ad2_decreases)*100,1), 
                    "% (", round(quantile(ad2_decreases, 0.025)*100, 1), " - ", 
                    round(quantile(ad2_decreases, 0.975)*100, 1), 
                    "%) of second-level administrative units in Africa from 2000-2016"), 
             ind = indicator)

  # Draw-wise calculation of GVAP targets
  pct_met_gvap_threshold <- function(x) {
    x <- x[!is.na(x)]
    return(sum(x>=0.8)/length(x))
  }

  ad2_gvap_draws <- subset(merge(ad2_last_year_draws, sp_hierarchy_list), select = c("ADM0_NAME", draw_cols))
  ad2_gvap_draws[ , (draw_cols) := lapply(.SD, pct_met_gvap_threshold), by = ADM0_NAME, .SDcols = draw_cols]
  ad2_gvap_draws <- unique(ad2_gvap_draws)
  ad2_gvap_draws <- ad2_gvap_draws[ADM0_NAME != "Ma'tan al-Sarra"]
  ad2_gvap_matrix <- as.matrix(subset(ad2_gvap_draws, select = draw_cols))
  n_countries_met_gvap <- colSums(ad2_gvap_matrix == 1)

  write_np(paste0("For indicator ", indicator, ", estimated number of countries meeting GVAP target (by draw):", round(mean(n_countries_met_gvap),1), 
                "% (", round(quantile(n_countries_met_gvap, 0.025), 1), " - ", 
                round(quantile(n_countries_met_gvap, 0.975), 1), 
                "%) of second-level administrative units in Africa from 2000-2016"), 
         ind = indicator)

  ad2_gvap_country_probs <- Reduce(merge, list(ad2_gvap_draws[, .(mean = rowMeans(.SD)), .SDcols=draw_cols, by=ADM0_NAME],
                                               ad2_gvap_draws[, .(lower = as.list(quantile(.SD, 0.025))), .SDcols=draw_cols, by=ADM0_NAME],
                                               ad2_gvap_draws[, .(upper = as.list(quantile(.SD, 0.975))), .SDcols=draw_cols, by=ADM0_NAME]))
  
  # Write a table of the % of ad2s meeting GVAP target, by country, in 2016
  fwrite(ad2_gvap_country_probs, 
         file = paste0("<<<< FILEPATH REDACTED >>>>/vaccine/", indicator, "/output/", run_date, "/number_plugging/ad2_gvap_country_probs.csv"))
                    

  } else {
    write_np("Order of admin draws not correct for ad2 draw-level change comparisons", ind = indicator)
  }
}

## Last, generate a .csv of the proportion of districts meeting the GVAP threshold of >= 80% in each country for each year

# Reload draw object
load(paste0("<<<< FILEPATH REDACTED >>>>/vaccine/", indicator, "/output/", run_date, "/", indicator, "_raked_admin_draws_eb_bin0_0.RData"))
admin_2 <- merge(admin_2, sp_hierarchy_list, by = "ADM2_CODE")

draw_cols <- names(admin_2)[grepl("V[0-9]+", names(admin_2))]
year_list <- unique(admin_2$year)

admin_2[, (draw_cols) := lapply(.SD, function(x) as.numeric(x >= 0.8)), .SDcols = draw_cols]

missing_admins <- unique(admin_2[!complete.cases(admin_2)]$ADM2_CODE)
if (length(missing_admins) > 0) {
  warning(paste0("Dropping ", length(missing_admins), " missing ADM2_CODES:\n  ",
                 paste(missing_admins, collapse = "\n  ")))
  
  admin_2 <- subset(admin_2, !(ADM2_CODE %in% missing_admins)) 
}

# Add placeholder for total sum
admin_2[, n_admins := 1]

# Collapse to year/country
admin_2 <- admin_2[, lapply(.SD, sum), by = c("year", "ADM0_CODE"), .SDcols = c(draw_cols, "n_admins")]

ad2_estimates <- admin_2[, .(ADM0_CODE = ADM0_CODE,
                             year = year,
                             mean = rowMeans(.SD),
                             median = apply(.SD, 1, quantile, 0.5),
                             lower = apply(.SD, 1, quantile, 0.025),
                             upper = apply(.SD, 1, quantile, 0.975),
                             n_admins = n_admins),
                         .SDcols = draw_cols]

ad2_estimates <- merge(ad2_estimates, unique(subset(sp_hierarchy_list, select = c("ADM0_CODE", "ADM0_NAME"))))
ad2_estimates[, pct_mean := mean / n_admins]
ad2_estimates[, pct_median := median / n_admins]
ad2_estimates[, pct_lower := lower / n_admins]
ad2_estimates[, pct_upper := upper / n_admins]

setcolorder(ad2_estimates, c("ADM0_CODE", "ADM0_NAME", "year", "n_admins",
                             "mean", "median", "lower", "upper", 
                             "pct_mean", "pct_median", "pct_lower", "pct_upper"))

 fwrite(ad2_estimates, 
         file = paste0("<<<< FILEPATH REDACTED >>>>/vaccine/", indicator, "/output/", run_date, "/number_plugging/ad2_pct_districts_GVAP_target.csv"))
                    