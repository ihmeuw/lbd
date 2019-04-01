###############################################################################
###############################################################################
## Plot administrative-level changes
##
## Purpose: Produce scatter plots of admin2 changes over time by country
###############################################################################
###############################################################################

# Source directly from within master script

indicator <- "dpt3_cov"

library(ggplot2)
library(data.table)
library(tidyr)

in_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", 
                 indicator, "/output/", run_date, "/")

out_dir <- paste0(in_dir, "number_plugging/")

ad2_df <- fread(paste0(in_dir, "pred_derivatives/admin_summaries/",
                       indicator, "_admin_2_raked_summary.csv"))

#Subset to years of interest
ad2_df <- subset(ad2_df, year == max(year_list) | year == min(year_list))

# Convenience
ad2_df[year==max(year_list), year := 1]
ad2_df[year==min(year_list), year := 0]

# Convert to wide
ad2_df <- ad2_df %>% 
  gather(metric, value, mean:cfb) %>% 
  unite(metric_year, metric, year) %>%
  spread(metric_year, value) %>%
  as.data.table

# Grab populations
load(paste0(in_dir, indicator, "_raked_admin_draws_eb_bin0_0.RData"))
pops <- subset(admin_2, year == max(year_list), select = c(ADM2_CODE, pop))
ad2_df <- merge(ad2_df, pops, by = "ADM2_CODE", all.x= T, all.y = F)

# Figure out which had significant increases or decreases
decreases <- fread(paste0(in_dir, "/number_plugging/admin_2_",
                          max(year_list), " - ", min(year_list), "_decline.csv"))

increases <- fread(paste0(in_dir, "/number_plugging/admin_2_",
                          max(year_list), " - ", min(year_list), "_increase.csv"))

sig_ad2s <- as.numeric(unique(c(increases$ADM2_CODE, decreases$ADM2_CODE)))

ad2_df[, sig := ifelse(ADM2_CODE %in% sig_ad2s, T, F)]

# Drop Ma'tan al Sarra
ad2_df<-subset(ad2_df, ADM0_NAME != "Ma'tan al-Sarra")

# Add in some probabilities of > 80% in 2015
ad2_80p_target <- fread(paste0(in_dir, "/pred_derivatives/admin_summaries/",
                               indicator, "_admin_2_raked_p_0.8_or_better_summary.csv"))
ad2_80p_target <- subset(ad2_80p_target, year == max(year_list),
                         select = c("ADM2_CODE", "p_above"))

ad2_df <- merge(ad2_df, ad2_80p_target, by= "ADM2_CODE")
ad2_df[, p_above_cut := cut(p_above, breaks = seq(0,1,by=0.1))]

# For faceted plots
ad2_df[ADM0_NAME == "Democratic Republic of the Congo", ADM0_NAME := "DRC"]
ad2_df[ADM0_NAME == "Sao Tome and Principe", ADM0_NAME := "STP"]
ad2_df[ADM0_NAME == "United Republic of Tanzania", ADM0_NAME := "Tanzania"]
ad2_df[ADM0_NAME == "Central African Republic", ADM0_NAME := "CAR"]

# Use carto color scheme
carto_colors <- c("#5F4690","#38A6A5","#73AF48","#E17C05","#94346E","#6F4070","#994E95")

carto_diverging <- c("#009392","#39b185","#9ccb86","#e9e29c","#eeb479","#e88471","#cf597e")

big_colors <- c("#d34a53","#36dee6","#862b1b","#45bc8d","#d74d82","#64b764","#ad74d6","#92b440",
                "#5f6ed3","#caa331","#573586","#c2a957","#628bd5","#c8772c","#cd76c7","#627123",
                "#a54381","#c6673c","#ac455f","#d17160")


# Create a plot of first vs last year estimates (all countries together)

png(filename=paste0(out_dir, "ad2_", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=300)

gg <- ggplot() +
  geom_point(data = ad2_df,
             aes(x = mean_0, y = mean_1, size = pop, color = p_above, shape = sig),
             alpha = 0.6) +
  geom_abline(alpha = 0.6) +
  geom_abline(slope = 0, intercept = 0.8, 
              color = "black", alpha = 0.6,
              linetype = "dotted") +
  scale_size_area() +
  scale_color_gradientn(colors = rev(carto_diverging)) +
  theme_classic() +
  scale_shape_manual(values = c(1,19)) +
  labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
       y = paste0("Mean DPT3 Coverage: ", max(year_list)),
       size = "Number of children \n< 5 years of age",
       color = paste0("Probability of \nDPT3 coverage > 80% \n(", max(year_list), ")"),
       shape = "High (>95%) probability \nof true change 2000-2015") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))

print(gg)

dev.off()

# Create a plot of first vs last year estimates (faceted by country)

png(filename=paste0(out_dir, "ad2_by_adm0_", min(year_list), "-", max(year_list), ".png"), 
    type="cairo",
    units="in", 
    width=12, 
    height=14, 
    pointsize=12, 
    res=300)

gg <- ggplot() +
  geom_point(data = ad2_df,
             aes(x = mean_0, y = mean_1, size = pop, color = p_above, shape = sig),
             alpha = 0.6) +
  geom_abline(alpha = 0.6) +
  geom_abline(slope = 0, intercept = 0.8, 
              color = "black", alpha = 0.6,
              linetype = "dotted") +
  scale_size_area() +
  scale_color_gradientn(colors = rev(carto_diverging), limits = c(0,1)) +
  facet_wrap(~ADM0_NAME, ncol = 6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(1,19)) +
  labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
       y = paste0("Mean DPT3 Coverage: ", max(year_list)),
       size = "Number of children \n< 5 years of age",
       color = paste0("Probability of \nDPT3 coverage > 80% \n(", max(year_list), ")"),
       shape = "High (>95%) probability \nof true change 2000-2015") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))

print(gg)

dev.off()

# Create plots of first vs last year estimates (one for each country)

for (adname in unique(ad2_df$ADM0_NAME)) {
  country_df <- subset(ad2_df, ADM0_NAME == adname)
  
  png(filename=paste0(out_dir, "ad2_by_adm0_", adname, "_", min(year_list), "-", max(year_list), ".png"), 
      type="cairo",
      units="in", 
      width=8, 
      height=6, 
      pointsize=12, 
      res=300)
  
  gg <- ggplot() +
    geom_point(data = country_df,
               aes(x = mean_0, y = mean_1, size = pop, color = p_above, shape = sig),
               alpha = 0.6) +
    geom_abline(alpha = 0.6) +
    geom_abline(slope = 0, intercept = 0.8, 
                color = "black", alpha = 0.6,
                linetype = "dotted") +
    scale_size_area() +
    scale_color_gradientn(colors = rev(carto_diverging), limits = c(0,1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_shape_manual(values = c(1,19)) +
    labs(x = paste0("Mean DPT3 Coverage: ", min(year_list)),
         y = paste0("Mean DPT3 Coverage: ", max(year_list)),
         size = "Number of children \n< 5 years of age",
         color = paste0("Probability of \nDPT3 coverage > 80% \n(", max(year_list), ")"),
         shape = "High (>95%) probability \nof true change 2000-2015") +
    scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,1))
  
  print(gg)
  
  dev.off()
  
}
