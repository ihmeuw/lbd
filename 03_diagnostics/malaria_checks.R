# (1) check sums in summary csvs
admin <- 0

yll <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(mean_yll = mean, upper_yll = upper, lower_yll = lower)

yld <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(mean_yld = mean, upper_yld = upper, lower_yld = lower)

daly <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(mean_daly = mean, upper_daly = upper, lower_daly = lower)

loc_info <- grep(names(daly), pattern = 'ADM', value = T)

sum <- merge(yll, yld, by = c(loc_info,'region','year')) %>%
  merge(daly, by = c(loc_info,'region','year'))

#rates
sum[,yll_yld := mean_yll + mean_yld]

plot <- ggplot(data = sum, aes(x = mean_daly, y = yll_yld)) + geom_point() + geom_abline(slope = 1, intercept = 0, color = 'red') + ggtitle(paste0('ADM', admin))
plot(plot)


# (2) check sums in admin draws files
#load Yll rates
run_date <- '2019_07_24'
load('<<< FILEPATH REDACTED >>>')
yll_rate_ad0 <- admin_0
yll_rate_ad1 <- admin_1
yll_rate_ad2 <- admin_2

#load yld rates
load('<<< FILEPATH REDACTED >>>')
yld_rate_ad0 <- admin_0
yld_rate_ad1 <- admin_1
yld_rate_ad2 <- admin_2

#load calculated dalys
load('<<< FILEPATH REDACTED >>>')
daly_rate_ad0 <- admin_0
daly_rate_ad1 <- admin_1
daly_rate_ad2 <- admin_2


#sum to dalys rates
admin_0 <- merge(yll_rate_ad0, yld_rate_ad0, by = c('year','ADM0_CODE','region','pop'))
for (n in 1:100){
  admin_0[,paste0('V', n) :=  rowSums(.SD, na.rm=T), .SDcols= c(paste0('V',n, '.x'), paste0('V',n, '.y'))]
}

admin_1 <- merge(yll_rate_ad1, yld_rate_ad1, by = c('year','ADM1_CODE','region','pop'))
for (n in 1:100){
  admin_1[,paste0('V', n) :=  rowSums(.SD, na.rm=T), .SDcols= c(paste0('V',n, '.x'), paste0('V',n, '.y'))]
}

admin_2 <- merge(yll_rate_ad2, yld_rate_ad2, by = c('year','ADM2_CODE','region','pop'))
for (n in 1:100){
  admin_2[,paste0('V', n) :=  rowSums(.SD, na.rm=T), .SDcols= c(paste0('V',n, '.x'), paste0('V',n, '.y'))]
}

# first, check that the above sum works:
for (n in 1:100){
  yld_ind <- n + 4
  yll_ind <- n + 104
  daly_ind <- n + 204
  
  yll <- admin_2[,..yll_ind] %>%
    as.matrix()
  yld <- admin_2[,..yld_ind] %>%
    as.matrix()
  daly <- admin_2[,..daly_ind] %>%
    as.matrix()
  sum <- yll + yld
  dt <- data.frame(daly, sum)
  mod <- lm(daly ~ sum, data=dt)
  print(mod$coefficients[2])
}


#take means for each admin
daly <- daly_rate_ad2[,daly_mean :=  rowMeans(.SD, na.rm=T), .SDcols= c(paste0('V',1:100))]
yld <- yld_rate_ad2[,yld_mean :=  rowMeans(.SD, na.rm=T), .SDcols= c(paste0('V',1:100))]
yll <- yll_rate_ad2[,yll_mean :=  rowMeans(.SD, na.rm=T), .SDcols= c(paste0('V',1:100))]

merge <- merge(yll, yld, by = c('ADM2_CODE','year')) %>%
  merge(daly, by = c('ADM2_CODE', 'year') )

merge[,sum := yld_mean + yll_mean]
plot <- ggplot(data = merge, aes(x = sum, y = daly_mean)) + geom_point() + geom_abline(slope = 1, intercept = 0, color = 'red')
plot(plot)


raster <- raster('<<< FILEPATH REDACTED >>>')
missing <- raster

missing[missing != 0] <- 1
missing[missing == 0] <- 100
spplot(missing)
