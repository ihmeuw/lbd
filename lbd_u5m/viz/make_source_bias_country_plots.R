library(raster)
library(data.table)
library(rgdal)
library(rgeos)
library(dplyr)
library(openxlsx)
library(foreign)
library(parallel)
library(stringi)
library(gridExtra)

#update with core repo path
source("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")

#args for pulling source bias
regions <- c('soas', 'ocea+seas-mys', 'stan+mng', 'caca-mex', 'ansa+trsa-bra', 'noaf', 'mide+yem', 'cssa', 'essa-yem', 'sssa', 'wssa')
run_date <- "<<<< REDACTED >>>>"
plot <- F
plot_path <- "<<<< FILEPATH REDACTED >>>>"
cores <- 11

#get source bias from survey-level effects
source_bias_table <- pull_u5m_source_bias_adjustment(regions = regions, 
                                                     run_date = run_date,
                                                     plot = plot,
                                                     plot_path = plot_path,
                                                     cores = cores)

x <- source_bias_table
source_bias <- source_bias_table[, .(nid, source_bias_adj)]

#data that feeds directly into the model
train <- mclapply(X=regions, FUN = function(x) {fread(sprintf("<<<< FILEPATH REDACTED >>>>", run_date, x))}, mc.cores = cores)
train2 <- do.call("rbind", c(train, list(fill = T)))

train_agg <- train2[, .(q_ab = sum(died) / sum(N)), by=c("nid", "year", "country", "ab")]
train_agg2 <- train_agg[, .(q = 1 - prod(1 - q_ab)), by=c("nid", "year", "country")]

full_data <- merge(train_agg2, source_bias, by="nid")

pdf("<<<< FILEPATH REDACTED >>>>")
for(coun in sort(unique(full_data$country))){
  plot_data <- full_data[country == coun,]
  plot_data$nid <- as.factor(plot_data$nid)
  gg <- ggplot(data = plot_data) +
    geom_point(aes(x = year, y = q, color = source_bias_adj, shape = nid))+
    ggtitle(coun) +
    theme_minimal() +
    scale_color_gradient2(low="red", mid = "green", high = "blue") +
    labs(color = "Bias Adjustment") +
    ylab("mean weighted haz") +
    scale_shape_manual(values=1:nlevels(plot_data$nid))
  
  print(gg)
}
dev.off()

##############################################################################
#random effect distribution
types <- unique(train2[, c("nid", "data_type")])
d_plot <- merge(source_bias_table, types, by = "nid")
d_plot <- unique(d_plot[, c("nid", "data_type", "random_effect")])

cbh <- d_plot[data_type == "cbh",]
sbh <- d_plot[data_type == "sbh",]
summary(cbh$random_effect)
summary(sbh$random_effect)

gg <- ggplot() +
  geom_histogram(data = d_plot, aes(x = random_effect)) +
  facet_wrap(facets = ~data_type)
  
pdf("~/test.pdf")
print(gg)
dev.off()



#############################################################################
# survey pairs, 5 years back

#look at rows dropped
y <- train2
train2 <- train2[year <= svyyr,]

train_agg <- train2[, .(q_ab = sum(died) / sum(N), N = sum(N)), by=c("nid", "year", "country", "ab")]
train_agg2 <- train_agg[, .(q = 1 - prod(1 - q_ab), N = sum(N)), by=c("nid", "year", "country")]

train_agg2 <- train_agg2[!(country %in% c("CHN", "BRA", "MEX")),]

final_data = data.table(long_avg = numeric(),
                        short_avg = numeric(),
                        year_diff = integer(),
                        long_survey = integer(),
                        short_survey = integer())

# loop through each country
for(coun in unique(train_agg2$country)) {
  message(coun)
  #order nid by latest year
  coun_data <- train_agg2[country == coun,]
  niddict <- unique(coun_data[, c("nid", "year")])
  niddict <- setorder(niddict, -year)
  niddict[, max_year := max(year), by=nid]
  niddict2 <- unique(niddict[, c("nid", "max_year")])
  
  #find survey pairs to compare
  nid_list <- list()
  for(n in 1:length(niddict2$nid)){
    start_row <- n+1
    m_year <- niddict2[n,]$max_year
    nid <- niddict2[n,]$nid
    x <- niddict2[start_row:nrow(niddict2),]
    x <- x[max_year < m_year]
    if(nrow(x) > 0){
      nid_temp <-  x$nid
      nid_list <- c(nid_list, list(nid_temp))
    }
  }
  
  if(length(nid_list) == 0) {
    next
  }
  
  #loop through each survey pair within a country
  for(i in 1:length(nid_list)){
    long_survey <- niddict2[i, ]$nid
    long_data <- setorder(coun_data[nid == long_survey], -year)
    long_max_year <- long_data[1]$year
    
    for(j in nid_list[[i]]){
      short_survey <- j
      short_data <- setorder(coun_data[nid == short_survey], -year)
      short_max_year <- short_data[1]$year - 1
      
      years <- seq(short_max_year, short_max_year - 4)
      
      if(short_max_year <= 2002 | nrow(short_data) < 5) {
        next
      }
      
      if(nrow(long_data[year == short_max_year]) == 0) {
        next
      }
      
      q_adj <- long_data[year == short_max_year,]$q - short_data[year == short_max_year,]$q
      short_data[, q := q + q_adj]
      
      long_avg <- weighted.mean(long_data[year %in% years]$q, long_data[year %in% years]$N, na.rm = T)
      short_avg <- weighted.mean(short_data[year %in% years]$q, short_data[year %in% years]$N, na.rm = T)
      year_diff <- long_max_year - short_max_year
      
      temp_data <- data.table(long_avg = long_avg,
                              short_avg = short_avg,
                              year_diff = year_diff,
                              long_survey = long_survey,
                              short_survey = j)
      
      final_data <- rbind(final_data, temp_data)
    }
  }
}
#find adjustment to bring shorter survey to longer survey based on most recent year of shorter survey
#adjust all shorter values
#calculate 5-yr mean for shorter and longer
#find end_year longer - end_year_shorter
#plot points, year on x, difference on y

final_data[, avg_diff := long_avg / short_avg]
final_data <- unique(final_data)

gg <- ggplot(data = final_data) +
  geom_jitter(aes(y = avg_diff, x = year_diff, alpha = .5)) +
  theme_bw() +
  xlab("Years between start of surveys") +
  ylab("Difference in average 5q0 over 5 years between survey pairs, adjusted q in latest shared year") +
  geom_hline(yintercept = 0, col="red") +
  geom_smooth(aes(x = year_diff, y = avg_diff), method = "loess") +
  theme(legend.position = "none")
  

png("<<<< FILEPATH REDACTED >>>>",
    width=8,
    height=8,
    units = "in",
    res = 300)
print(gg)
dev.off()
