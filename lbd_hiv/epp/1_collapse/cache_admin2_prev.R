#######
# Set up Environment
#######


## Clear environment
rm(list = ls())

## Set indicator
indicator_group <- 'hiv'
indicator       <- 'epp'

## Set repos
core_repo  <- paste0("/<<<< FILEPATH REDACTED >>>>/", "/lbd_core/")
indic_repo <- paste0("/<<<< FILEPATH REDACTED >>>>/", "/lbd_hiv/")
code.dir <- paste0("/<<<< FILEPATH REDACTED >>>>/", "/HIV/")
setwd(core_repo)

## Load libraries and  MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "share_scripts/common_inputs/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = core_repo)
library(MASS)


######################################################
output_date <- "<<<<FILEPATH REDACTED>>>>>"
input_run   <- "<<<<FILEPATH REDACTED>>>>>"
dir.create(paste0("<<<< FILEPATH REDACTED >>>>",output_date,"/"))
dir.create(paste0("<<<< FILEPATH REDACTED >>>>",output_date,"/"))
######################################################




load(paste0("<<<< FILEPATH REDACTED >>>>"))
means_0 <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
means_2 <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))


codes_0 <- unique(means_0$ADM0_CODE)




## Probit case
for (c in codes_0) {
  means <- means_0[which(means_0$ADM0_CODE == c),]
  means <- qnorm(means$mean)


  data <- admin_0[which(admin_0$ADM0_CODE == c), ]

  data$ADM0_CODE <- NULL
  data$pop <- NULL
  transposed <- melt(data, id.vars = c("year"))
  transposed <- dcast(transposed, variable ~ year)
  transposed$variable <- NULL

  cov_mat <- cov(qnorm(as.matrix(transposed)))

  mat <- cov_mat
  save(means, mat, file = paste0("<<<< FILEPATH REDACTED >>>>"))
}


## Un-Transformed case
for (c in codes_0) {
  means <- means_0[which(means_0$ADM0_CODE == c),]
  means <- means$mean


  data <- admin_0[which(admin_0$ADM0_CODE == c), ]

  data$ADM0_CODE <- NULL
  data$pop <- NULL
  transposed <- melt(data, id.vars = c("year"))
  transposed <- dcast(transposed, variable ~ year)
  transposed$variable <- NULL

  cov_mat <- cov(as.matrix(transposed))

  mat <- cov_mat
  save(means, mat, file = paste0("<<<< FILEPATH REDACTED >>>>"))
}


## logit case
for (c in codes_0) {
  means <- means_0[which(means_0$ADM0_CODE == c),]
  means <- logit(means$mean)


  data <- admin_0[which(admin_0$ADM0_CODE == c), ]

  data$ADM0_CODE <- NULL
  data$pop <- NULL
  transposed <- melt(data, id.vars = c("year"))
  transposed <- dcast(transposed, variable ~ year)
  transposed$variable <- NULL

  cov_mat <- cov(logit(as.matrix(transposed)))

  mat <- cov_mat
  save(means, mat, file = paste0("<<<< FILEPATH REDACTED >>>>"))
}




admin_0 <- as.data.table(admin_0)

admin_0$prev <-  rowMeans(admin_0[ ,V1:V1000])
admin_0$var  <- -1
admin_0$se   <- -1
admin_0$n    <- -1


for (r in 1:nrow(admin_0)) {
  message(r)
  cy <- unlist(as.list(admin_0[r, V1:V1000]))
  var_cy <- var(cy)
  var_cy
  se_cy <- sqrt(var_cy)
  se_cy
  p_cy <- admin_0[r, prev]
  p_cy
  n_cy <- (p_cy * (1 - p_cy)) / (var_cy)
  n_cy
  low_cy <- as.numeric(quantile(cy, probs = c(0.025), na.rm = TRUE))
  high_cy <- as.numeric(quantile(cy, probs = c(0.975), na.rm = TRUE))
  ci_cy <- high_cy - low_cy
  admin_0[r, 1005] <- var_cy
  admin_0[r, 1006] <- se_cy
  admin_0[r, 1007] <- n_cy
}




vars <- c(2,1, 1004, 1006, 1007)
admin_0 <- as.data.frame(admin_0)
hhs <- admin_0[ , vars]
hhs$used <- TRUE
write.csv(hhs, paste0("<<<< FILEPATH REDACTED >>>>"))







codes_2 <- unique(means_2$ADM2_CODE)

## Probit case
for (c in codes_2) {
  means <- means_2[which(means_2$ADM2_CODE == c),]
  means <- qnorm(means$mean)


  data <- admin_2[which(admin_2$ADM2_CODE == c), ]

  data$ADM2_CODE <- NULL
  data$pop <- NULL
  transposed <- melt(data, id.vars = c("year"))
  transposed <- dcast(transposed, variable ~ year)
  transposed$variable <- NULL

  cov_mat <- cov(qnorm(as.matrix(transposed)))

  mat <- cov_mat
  save(means, mat, file = paste0("<<<< FILEPATH REDACTED >>>>"))
}

## Un-Transformed case
for (c in codes_2) {
  means <- means_2[which(means_2$ADM2_CODE == c),]
  means <- means$mean


  data <- admin_2[which(admin_2$ADM2_CODE == c), ]

  data$ADM2_CODE <- NULL
  data$pop <- NULL
  transposed <- melt(data, id.vars = c("year"))
  transposed <- dcast(transposed, variable ~ year)
  transposed$variable <- NULL

  cov_mat <- cov(as.matrix(transposed))

  mat <- cov_mat
  save(means, mat, file = paste0("<<<< FILEPATH REDACTED >>>>"))
}

## Logit case
for (c in codes_2) {
  means <- means_2[which(means_2$ADM2_CODE == c),]
  means <- logit(means$mean)


  data <- admin_2[which(admin_2$ADM2_CODE == c), ]

  data$ADM2_CODE <- NULL
  data$pop <- NULL
  transposed <- melt(data, id.vars = c("year"))
  transposed <- dcast(transposed, variable ~ year)
  transposed$variable <- NULL

  cov_mat <- cov(logit(as.matrix(transposed)))

  mat <- cov_mat
  save(means, mat, file = paste0("<<<< FILEPATH REDACTED >>>>"))
}



admin_2 <- merge(admin_2, sp_hierarchy_list, by = c("ADM2_CODE"))
admin_2 <- as.data.table(admin_2)
admin_2$prev <-  rowMeans(admin_2[ ,V1:V1000])
admin_2$var  <- -1
admin_2$se   <- -1
admin_2$n    <- -1
for (r in 1:nrow(admin_2)) {
  message(r)
  cy <- unlist(as.list(admin_2[r, V1:V1000]))
    var_cy <- var(cy)
    var_cy
    se_cy <- sqrt(var_cy)
    se_cy
    p_cy <- admin_2[r, prev]
    p_cy
    n_cy <- (p_cy * (1 - p_cy)) / (var_cy)
    n_cy
    low_cy <- as.numeric(quantile(cy, probs = c(0.025), na.rm = TRUE))
    high_cy <- as.numeric(quantile(cy, probs = c(0.975), na.rm = TRUE))
    ci_cy <- high_cy - low_cy
    admin_2[r, 1011] <- var_cy
    admin_2[r, 1012] <- se_cy
    admin_2[r, 1013] <- n_cy
}

vars <- c("ADM2_CODE","year", "prev", "se", "n")
admin_2 <- as.data.frame(admin_2)
hhs <- admin_2[ , vars]
hhs$used <- TRUE

write.csv(hhs, paste0("<<<< FILEPATH REDACTED >>>>"))


########### ADM2 art proportions work arround
load(paste0("<<<< FILEPATH REDACTED >>>>"))
admin_2 <- as.data.table(admin_2)
vars2 <- c("year", "ADM2_CODE", "mean")
vars0 <- c("year", "ADM0_CODE", "mean")

admin_2$mean <- rowMeans(admin_2[ ,V1:V1000])
admin_0$mean <- rowMeans(admin_0[ ,V1:V1000])

admin_2 <- as.data.frame(admin_2)
admin_0 <- as.data.frame(admin_0)

admin_2 <- admin_2[ , vars2]
admin_0 <- admin_0[ , vars0]

admin_2 <- merge(admin_2, sp_hierarchy_list, by = c("ADM2_CODE"))

art_props <- merge(admin_2, admin_0, by = c("ADM0_CODE", "year"))
art_props$prev_prop <- art_props$mean.x / art_props$mean.y




write.csv(art_props, paste0("<<<< FILEPATH REDACTED >>>>"))


### set Admin 2 exclusions, these are the admin 2 units where the upper bound of the PLHIV CI is always below 1 in MBG modeled output.

plhiv <- as.data.frame(fread(paste0("<<<< FILEPATH REDACTED >>>>")))

plhiv_prob <- plhiv[which(plhiv$upper < 1), ]


n_years <- length(unique(plhiv$year))

table <- as.data.frame(table(plhiv_prob$ADM2_CODE))
ex <- as.list(table$Var1[which(table$Freq == n_years)])

saveRDS(ex, paste0("<<<< FILEPATH REDACTED >>>>"))
