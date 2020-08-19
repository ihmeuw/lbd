######## Combine region- and rundate-specific model results into global rasters
library(raster)
library(data.table)

#### Setup
run_dates <- <<<< FILEPATH REDACTED >>>>
regions <- c("lf_endem_afr", "lf_s_se_asia", "lf_hispaniola")
viz_location <- <<<< FILEPATH REDACTED >>>>
year_list <- c(2000:2017)
measures <- c("cases", "prev")
measures_output <- c("count", "prev")
summ_stats <- c("lower", "mean", "upper")

#### Merge regional rasters by year, measure and summary stat, and save to file
for (a in 1:length(year_list)) {
  print(year_list[[a]])
  for (b in 1:length(measures)) {
    print(measures[[b]])
    for (c in 1:length(summ_stats)) {
      print(summ_stats[[c]])
      template <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
      for (i in 2:length(regions)) {
        template <- raster::merge(template, raster(paste0(<<<< FILEPATH REDACTED >>>>)))
      }
      writeRaster(template, paste0(<<<< FILEPATH REDACTED >>>>), format = "GTiff")
      print("Done")
    }
  }
}

#### Merge and save posterior prob. < 2% rasters
rm(template)
template <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
for (a in 2:length(run_dates)) {
  template <- raster::merge(template, raster(paste0(<<<< FILEPATH REDACTED >>>>)))
}
writeRaster(template, paste0(<<<< FILEPATH REDACTED >>>>), format = "GTiff")

#### Merge and save posterior prob. < 2% CSVs by admin

all <- data.table("id" = integer(), "absolute_goal_prob" = numeric())
for (a in 1:length(run_dates)) {
  adm0 <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
  adm1 <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
  adm2 <- fread(paste0(<<<< FILEPATH REDACTED >>>>))

  setnames(adm0, "ADM0_CODE", "id")
  setnames(adm1, "ADM1_CODE", "id")
  setnames(adm2, "ADM2_CODE", "id")

  all <- rbindlist(list(all, as.data.table(adm0[, c("id", "absolute_goal_prob")]), as.data.table(adm1[, c("id", "absolute_goal_prob")]), as.data.table(adm2[, c("id", "absolute_goal_prob")])))
}

setnames(all, "absolute_goal_prob", "value")
setorder(all, "id")
fwrite(all, file = <<<< FILEPATH REDACTED >>>>)
