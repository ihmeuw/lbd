####################################################################################################
## Description:   Make plots of national-level temporal trends and data.
##
## Inputs:        Pre-collapse microdata (<<<< FILEPATH REDACTED >>>>.RDS)
##            
##                Collapsed MBG input data (<<<< FILEPATH REDACTED >>>>)
##                Unraked admin0 predictions (<<<< FILEPATH REDACTED >>>>)
##                Raked admin0 predictions (<<<< FILEPATH REDACTED >>>>)
##
## Output:        PDF of plots (<<<< FILEPATH REDACTED >>>>)
####################################################################################################

require(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
require(RColorBrewer)

admin0_data_and_estimates_plots <- function(data_date, run_date, indicator, indicator_group = "hiv", weight = "pweight") {

## Load data and estimates -------------------------------------------------------------------------
# national estimates
pred_dir <- paste0("<<<< FILEPATH REDACTED >>>>", "/pred_derivatives/admin_summaries/")
pred <- fread(paste0(pred_dir, indicator, "_admin_0_unraked_summary.csv"))

inla <- fread(paste0(pred_dir, indicator, "_admin_0_unraked_summary.csv"))
inla <- inla[, list(type = "unraked", ADM0_CODE, ADM0_NAME, year, mean, lower, upper)]


if (file.exists(paste0(pred_dir, indicator, "_admin_0_stackers.csv"))) {
  stackers <- fread(paste0(pred_dir, indicator, "_admin_0_stackers.csv"))
  stackers <- melt(stackers, id.vars = c("ADM0_CODE", "year"), value.name = "mean", variable.name = "type")
  levels(stackers$type) <- paste("stackers:", levels(stackers$type))
} else {
  stackers <- inla[type == "stacker",] # empty data frame, but with the correct columns
  stackers[, type := factor(type)]
}

pred <- rbind(inla, stackers, fill=T)
pred[, ADM0_NAME := ADM0_NAME[1], ADM0_CODE]
pred[, type := factor(type, levels = c("unraked", rev(levels(stackers$type))))]
rm(pred_dir, inla, stackers)

#input data WANT TO EVENTUALLY MOVE AWAY FROM DATA TAG
input_data <- fread(paste0("<<<< FILEPATH REDACTED >>>>", "/input_data.csv"))
input_data <- input_data[input_data$year %in% pred$year, list(nid, country, source, year, male_circumcision = male_circumcision / N)]

# microdata (collapse to national) DATA date is different
data_date <- as.Date(data_date, "%Y_%m_%d")
micro <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
micro <- micro[nid %in% input_data$nid,]
micro <- micro[between(age_year, 15, 49) & !is.na(male_circumcision) & !is.na(pweight),]
micro <- micro[, list(int_year = floor(median(int_year, na.rm = T)),
                      male_circumcision = weighted.mean(male_circumcision, pweight)),
               by = 'nid,country,survey_name']
input_data[nid %in% micro$nid, type := "Microdata"]

# report data
report <- fread(paste0("/lbd_hiv/1_report_extraction/male_circumcision/circumcision_extraction.csv"))
report <- report[nid %in% input_data$nid & !nid %in% micro$nid,]
input_data[nid %in% report$nid & !nid %in% micro$nid, type := "Report"]
report <- report[start_age ==  15 & end_age ==  49 & location_type ==  "National",]
report <- report[, list(nid, country, survey_name, int_year, male_circumcision = male_circumcision / 100)]

# Annual number of voluntary medical male circumcision in eastern and southern Africa by country, 2008-2016
location_name <- c(replicate(9,"Botswana"), replicate(9,"Ethiopia"), replicate(9,"Kenya"), replicate(9, "Lesotho"),
             replicate(9, "Malawi"), replicate(9, "Mozambique"), replicate(9, "Namibia"), replicate(9, "Rwanda"), 
             replicate(9,"South Africa"), replicate(9,"Swaziland"), replicate(9,"Uganda"), replicate(9,"Tanzania"),
             replicate(9,"Zambia"), replicate(9,"Zimbabwe"))
country <- c(replicate(9,"BWA"), replicate(9,"ETH"), replicate(9,"KEN"), replicate(9, "LSO"),
                   replicate(9, "MWI"), replicate(9, "MOZ"), replicate(9, "NAM"), replicate(9, "RWA"), 
                   replicate(9,"ZAF"), replicate(9,"SWZ"), replicate(9,"UGA"), replicate(9,"TZA"),
                   replicate(9,"ZMB"), replicate(9,"ZWE"))

year <- rep(2008:2016, 14)

#Country data 
circumcisions <- c(0, 5424, 5773, 14661, 38005, 46793, 30033, 15722, 24042,
0, 769, 2689, 7542, 11961, 16393, 11831, 9744, 10306,
11663, 80719, 139905, 159196, 151517, 190580, 193576, 207014, 243447,
0,0,0,0,10835, 37655, 36245, 25966, 34157,
589, 1234, 1296, 11811, 21250, 40835, 80419, 108672, 129975,
0,100, 7633, 29592, 135000, 146046, 240507, 198340, 253079, 
0,224, 1763, 6123, 4863, 1182, 4165, 18459, 27340,
0, 0, 1694, 25000, 138711, 116029, 173191, 138216, 137218,
190, 9168, 131117, 296726, 422009, 514991, 482474, 485552, 497186,
1110, 4336, 18869, 13791, 9977, 10105, 12289, 12952, 17374, 
0, 0, 21072, 77756, 368490, 801678, 878109, 556546, 411459,
0, 1033, 18026, 120261, 183480, 329729, 573845, 435302, 548390,
2758, 17180, 61911, 85151, 173992, 294466, 315168, 222481, 311792,
0, 2801, 11176, 36603, 40755, 112084, 209125, 188732, 205784)
df <- data.frame(country, year, circumcisions)
df$country <- as.character(df$country)


## Make plots --------------------------------------------------------------------------------------

# get GAUL to iso mapping
gaul <- fread("<<<< FILEPATH REDACTED >>>>")
gaul <- gaul[level ==  3, list(location_id = loc_id, GAUL_CODE)]

source("<<<< FILEPATH REDACTED >>>>")
loc <- get_location_metadata(location_set_id = 1)
loc <- merge(loc[location_type ==  "admin0",], gaul)

dir.create(paste0("<<<< FILEPATH REDACTED >>>>"))
# loop over countries and make plots
pdf(paste0("<<<< FILEPATH REDACTED >>>>", "/male_circumcision_admin0_data_and_results.pdf"), width = 14, height = 8)
for (cc in pred[order(ADM0_NAME), unique(ADM0_CODE)]) {
  iso <- loc[GAUL_CODE ==  cc, ihme_loc_id]
  name <- loc[GAUL_CODE ==  cc, location_name]
  
  # main data and estimates plot
  
  plot_colors <- brewer.pal(nlevels(pred$type), "Set2")
  
  p1 <- ggplot() +
    geom_line(data = pred[ADM0_CODE ==  cc,], aes(x = year, y = mean)) +
    scale_color_manual(values = plot_colors) +
    scale_size_manual(values = c(1.5, 1.5, rep(0.5, length(plot_colors) - 2))) +
    scale_x_continuous(limits = range(pred$year) + c(-0.5, 0.5), expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 3)) +
    labs(y = "Male Circumcision Prevalence", title = name) +
    theme_bw() + theme(legend.position = "bottom", legend.direction = "horizontal",
                       legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0, "cm"),
                       axis.title.x = element_blank(), plot.margin = unit(rep(0.5, 4), "cm"))
  
  p2 <- ggplot() +
    geom_ribbon(data = pred[ADM0_CODE ==  cc & type == "unraked",], aes(x = year, ymin = lower, ymax = upper, fill = type), color = NA, alpha = 0.2, show.legend = F) +
    geom_line(data = pred[ADM0_CODE ==  cc & type == "unraked"], aes(x = year, y = mean, color = type), size = 1.5, show.legend = F) +
    geom_point(data = input_data[country ==  iso,], aes(x = year, y = male_circumcision), color = "gray40", shape = 1, size = 1, position = position_jitter(width = 0.2)) +
    geom_point(data = micro[country ==  iso,], aes(x = int_year, y = male_circumcision), color = "black", shape = 19, size = 5) +
    geom_point(data = report[country ==  iso,], aes(x = int_year, y = male_circumcision), color = "black", shape = 19, size = 5) +
    scale_color_manual(values = plot_colors[1]) +
    scale_fill_manual(values = plot_colors[1]) +
    scale_x_continuous(limits = range(pred$year) + c(-0.5, 0.5), expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 3)) +
    labs(x = "Year", y = "Male Cirumcision Prevalence") +
    theme_bw() + theme(plot.margin = unit(rep(0.5, 4), "cm"))

  # data table
  tab <- input_data[country ==  iso, list(source = unique(source), years = paste(sort(unique(year)), collapse = "-")), nid][order(years),]
  if (nrow(tab) ==  0) tab <- data.table(nid = "", source = "", years = "")
  tab <- tableGrob(tab, rows = NULL,
                   theme = ttheme_default(core = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 8)),
                                          colhead = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 8, fontface = "bold"))))
  
  scaleFUN <- function(x) round(as.numeric(x), digits=2)
  p3 <- ggplot(data = df[df$country == iso,], aes(x = year, y = cumsum(circumcisions))) + 
    geom_smooth(color= "grey", se = FALSE) +
    geom_ribbon(aes(ymin = 0, ymax = predict(loess(cumsum(circumcisions) ~ year))),
                alpha = 0.3,fill = 'grey') +
    scale_x_continuous(limits = c(2000, 2015) + c(-0.5, 0.5), expand = c(0, 0)) +
    xlab("Year") +
    ylab("Voluntary circumcisions") +
    theme_bw() + theme(legend.position = c(0.1, 0.9), plot.margin = unit(rep(0.5, 4), "cm"))
  
  # plot chart and table into one object
  if (iso %in% c("BTA", "ETH", "KEN", "LSO", "MWI", "MOZ", "NAM", "RWA", "ZAF", "SWZ", "UGA", "TZA", "ZMB", "ZWE")){
  grid.newpage()
  grid.draw(arrangeGrob(p1, p2, p3, tab, layout_matrix=matrix(c(rep(c(1, 1, 2, 2, 3), 3), 4, 4, 4, 4, 4), nrow = 5)))
  } else {
    grid.newpage()
    grid.draw(arrangeGrob(p1, p2, tab, layout_matrix=matrix(c(rep(c(1, 1, 2, 2, 2), 3), 3, 3, 3, 3, 3), nrow = 5)))
}
}
dev.off()

return("Plots saved!")
}
