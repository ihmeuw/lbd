
###########
## SETUP ##
###########
sys <- ifelse(Sys.info()[1] == 'Linux', 'cluster', 'local')
if(sys == 'local'){
  require(data.table)
}else{
  core_repo <- 'FILEPATH'
  commondir      <- paste('FILEPATH')
  package_list <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))
  source(paste0(core_repo, '/mbg_central/setup.R'))
  mbg_setup(package_list = package_list, repos = core_repo)
}

###############
## LOAD FILE ##
###############
date <- '20181018'
date <- '20181022'
if(sys == 'local') { ## ie on mac
  dat <- fread(sprintf('FILEPATH', date))
} else { ## ie on cluster
  dat <- fread(sprintf('FILEPATH'))
}
colnames(dat)

#######################
## FORMAT AND FILTER ##
#######################

## 1 ## split out parallel model jobs
non.par.runs <- subset(dat, !grepl("parallel_model", dat[,qsubbed_command]))
dat <- subset(dat, grepl("parallel_model", dat[,qsubbed_command]))

## 2 ## add region as a column 
dat[, region := unlist(lapply(strsplit(qsubbed_command,split=" "), function(x){x[3]}))]

## 3 ## split out TMB jobs from INLA jobs
tmb.runs  <- subset(dat, fit_with_tmb == TRUE)
dat <- subset(dat, fit_with_tmb != TRUE)

## 4 ## determine region size (pixels)
## this will take a while...
## but we need to get pixel counts from every region in dat
reg.pix <- data.table(region = sort(unique(dat$region)), pixels = 0.0)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~ setup  custom regions and older regions from other teams: ~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
gaul_ref <- list()
gaul_ref[["s2_essa"]]   <- c(170, 152, 205, 270, 133, 150, 220, 257, 253, 271, 58, 43, 235, 142)
gaul_ref[["s2_afr_horn"]]   <- c(226, 74, 61013, 79, 77, 70, 40760, 6, 102, 269)
gaul_ref[["s2_cssa"]]   <- c(89, 76, 214, 59, 68, 49, 8, 40762)
gaul_ref[["s2_wssa"]]   <- c(50, 45, 243, 221, 217, 182, 181, 159, 155, 144, 105, 90, 106, 94, 47, 66, 42, 29)
gaul_ref[["s2_name"]]   <- c(169, 145, 248, 268, 40765, 4)
gaul_ref[["s2_sssa"]]  <- c(172, 227, 35)
gaul_ref[["s2_mcacaf"]]  <- c(71, 72, 99, 123, 209, 211, 108, 191, 75, 180, 162, 111, 103, 61, 28)
gaul_ref[["s2_s_america"]]  <- c(57, 107, 73, 86, 194, 233, 260, 263, 195, 37, 33)
gaul_ref[["s2_central_asia"]]  <- c(250, 261, 239, 138)
gaul_ref[["s2_chn_mng"]]  <- c(167, 193, 67, 216, 147295, 2, 230)
gaul_ref[["s2_se_asia"]]  <- c(171, 44, 264, 139, 153, 240)
gaul_ref[["s2_malay"]]  <- c(196, 116, 242, 192)
gaul_ref[["s2_south_asia"]]  <- c(154, 175, 188, 115, 40781, 231, 31, 23, 15)
gaul_ref[["s2_mid_east"]]  <- c(117, 1, 118, 130, 238, 21)
gaul_ref[["s2_oceania"]]  <- c(262, 225, 83, 157, 163, 135, 212, 245, 5)
gaul_ref[["asia"]] <- c(116, 44, 139, 231, 154, 171, 153, 196, 220, 240, 242, 264, 147295, 67)
gaul_ref[["asia_png"]] <- c(116, 44, 139, 231, 154, 171, 153, 196, 220, 240, 242, 264, 147295, 67, 192)
gaul_ref[["china"]] <- c(147295)
gaul_ref[["asia_png_no_china"]] <- c(116, 44, 139, 231, 154, 171, 153, 196, 220, 240, 242, 264, 67, 192)
gaul_ref[["oceania"]] <- c(5, 83, 163, 135, 157, 192, 225, 245, 262, 212)
gaul_ref[["asia_w_oc"]] <- c(116, 44, 139, 231, 154, 171, 153, 196, 220, 240, 242, 264, 147295, 67, 5, 83, 163, 135, 157, 192, 225, 245, 262, 212)
gaul_ref[["latin_america"]] <- c(28, 71, 72, 99, 107, 108, 123, 209, 233, 211, 33, 73, 195, 37, 194, 260, 57, 61,75, 103, 111, 162, 180, 191, 263)
gaul_ref[["central_america"]] <- c(57, 61, 103, 111, 162, 180, 191,  75, 263)
gaul_ref[["south_america_carr"]] <- c(28, 71,72,  99, 107, 108, 123, 209, 233, 211,  33,  73, 195,  37, 194, 260)
gaul_ref[["south_asia"]] <- c(23, 31, 115, 175, 188)
gaul_ref[["middle_east"]] <- c(1, 21, 117, 118, 130, 138, 167, 187, 267, 238, 239, 250, 261, 269)
gaul_ref[["south_america"]] <- c(33,37,57,73,107,195,194,233,263)
gaul_ref[["lf_endem_afr"]] <- c(29,42,45,47,66,90,94,106,105,144,155,159,181,182,214,217,221,243,
                                8,49,59,68,76,89,50,6,74,
                                43,58,70,77,79,133,150,152,170,205,226,257,253,270,271)

#' @title Pull custom modeling regions
#' @description Define modeling regions that are not simple combinations of the
#'   default MBG regions (in other words, regions that are not combinations of
#'   four-letter MBG regions such as "wssa" or "seas+ocea" or ISO-3 codes such
#'   as 'ZAF' or 'CHN').
#' @param custom_region character vector of custom modeling regions
#' @return Returns a named list of custom modeling regions with associated
#'    "standard" (non-custom) modeling regions that can be directly interpreted
#'    by get_adm0_codes().
pull_custom_modeling_regions <- function(custom_regions){
  custom_regions <- tolower(custom_regions)

  # FULL LIST OF ALL REFERENCE REGIONS
  # If you need to add a new custom region, add it to this list
  ref_reg_list <- list(
    'africa'          = 'noaf+essa+wssa+cssa+sssa-yem',
    'middle_east'     = 'mide+stan-pak',
    'eastern_europe'  = "blr+est+lva+ltu+mda+ukr",
    'latin_america'   = 'caca+trsa+ansa',
    'south_asia'      = 'soas+chn_d2+pak',
    'central_america' = 'caca',
    'south_america'   = 'ansa+trsa',
    # se_asia was historically inclusive of East Asia + SE Asia
    'se_asia'         = 'eaas+seas+ocea+png', 
    
    #data coverage regions
    'africa_dcp' = 'noaf+essa+wssa+cssa+sssa+yem',
    'middle_east_dcp' = 'mide+stan-yem-pak+tur+isr+leb',
    'latin_america_dcp' = 'latin_america+cub',
    'south_asia_dcp' = 'south_asia-mdv-syc',
    'se_asia_dcp' = 'eaas+seas+png+idn+phl+tls+mys+twn',

      'stage1' = 'noaf+essa+wssa+cssa+sssa-yem',
    # ONLY stage 2 countries (not inclusive of Stage 1)
    'stage2' = 'ansa+caca+stan+eaas+mide+ocea+soas+seas+trsa+yem',
    'stage3' = 'all-stage1-stage2',

    'vax_soas' = 'soas+pak-syc',
    'vax_seas' = 'seas+ocea-asm-fji-kir-wsm-ton',
    'vax_eaas' = 'eaas',
    'vax_caeu' = 'arm+aze+geo+kgz+mda+tjk+tkm+ukr+uzb',
    'vax_crbn' = 'cub+dma+dom+grd+hti+jam+lca+vct+vir',
    'vax_ctam' = 'blz+col+cri+slv+gtm+hnd+mex+nic+pan+ven',
    'vax_ansa' = 'ansa-col-ven',
    'vax_trsa' = 'trsa',
    'vax_name' = 'noaf+mide+afg+omn',
    'vax_cssa' = 'cssa',
    'vax_essa' = 'essa+syc',
    'vax_sssa' = 'sssa',
    'vax_wssa' = 'wssa'
    
  )
  # Warn if there are any custom regions not in the reference list
  missing_regions <- custom_regions[ !(custom_regions %in% names(ref_reg_list)) ]
  if( length(missing_regions) > 0 ){
    message(paste0('WARNING: The following custom regions are not defined: ',
                   paste(missing_regions, collapse=','))
            )
  }
  # Return a named list of all custom regions
  custom_regions_list <- ref_reg_list[ names(ref_reg_list) %in% custom_regions ]
    return(custom_regions_list)
}


## NOTE: This part takes a while and should be done on the cluster
modeling_shapefile_version = 'current'
simple_raster_list <- simple_polygon_list <- NULL
missing.regs <- NULL
for(ii in 1:nrow(reg.pix)){
  reg <- reg.pix[ii, region]
  message(sprintf('on region %i out of %i: %s', ii, nrow(reg.pix), reg))

  ## get the adm0 codes and load the region shapefile
  ## 1) try get_adm_codes
  gaul_list    <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  ## 2) try custom coded region from gaul_ref
  if(length(gaul_list) == 0){ ## try getting from gaul_ref
    if(reg %in% names(gaul_ref)){
      gaul_list <- gaul_ref[[reg]]
    }
  }
  ## 3) if we still haven't found it, out of luck
  if(length(gaul_list) == 0){ 
    missing.regs <- c(missing.regs, reg)
    message(sprintf('--could not find adm0 list for region:', reg))
  }else{
    
    polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
                                        shapefile_version = modeling_shapefile_version)
    subset_shape       <- polygon_list[[1]]
    raster_list        <- build_simple_raster_pop(subset_shape)
    simple_raster      <- raster_list[['simple_raster']] ## get pixels from this

    ## get pixel count
    reg.pix[ii, pixels := ncell(simple_raster)]

    ## save things b/c it's slow
    simple_raster_list[[reg]] <- raster_list
    simple_polygon_list[[reg]] <- simple_polygon_list
  } 
}

## check what percent of model runs we've matched regions for
found.reg <- setdiff(reg.pix[, region], missing.regs)
message(sprintf('The number of runs we have pixel counts for is: %i out of %i', sum(dat$region %in% found.reg), nrow(dat)))
message(sprintf('The percent of runs we have pixel counts for is: %f', mean(dat$region %in% found.reg)))

## look at indicators submitting jobs for the unmatched regions
unmatched.dat <- subset(dat, region %in% missing.regs)
remaining.regs <- unique(cbind(unmatched.dat$indicator_group, unmatched.dat$region))
remaining.regs[order(remaining.regs[, 1]), ]


## turns out that we really wanted the number of non-NA pixels in simple raster!
## here we get that isntead of ncell(simple_raster)
for(ii in 1:nrow(reg.pix)){
  reg <- reg.pix[ii, region]
  if(reg %in% found.reg){
    reg.pix[ii, pixels := length(seegSDM:::notMissingIdx(simple_raster_list[[reg]]$simple_raster))]
  }
}


## merge on pixel counts to the matched portion 
dat <- subset(dat, region %in% found.reg)
dat <- merge(dat, reg.pix, by = 'region')

saveRDS(dat, file = sprintf('FILEPATH'))

###########################################################
## LOOK INTO RELATIONSHIP BETWEEN PIXELS, DRAWS, MAX MEM ##
###########################################################

## LOAD IN THE PROCESSED RUNS TO FIT THE MODEL

## the oct 18 extract has many more runs, so we use it to fit
train.date <- '20181018'
train.dat <- readRDS(sprintf('FILEPATH'))
train.dat[, log_ram_gb := log(ram_gb)]

## convert pixels into pixel-years
num.yrs  <- unlist(lapply(sapply(parse(text = train.dat[, year_list]), eval), length))
train.dat[, pixels := pixels * num.yrs]

## first, plot max mem by pixel and max mem by draw
par(mfrow = c(1, 2))
plot(train.dat$pixels, train.dat$ram_gb / 1000)
plot(train.dat$samples, train.dat$ram_gb / 1000)

## second, facet pixels vs mem by samples
require(ggplot2)
gg <- ggplot(train.dat, aes(x = pixels / 10000, y = ram_gb / 1000)) + xlab('pixel*yrs/10000') + ylab('RAM (GB)')
gg.ram <- gg + geom_point(aes(color = samples)) 
gg + geom_point() + facet_wrap( ~ samples, nrow = 1)
png('~/Desktop/SlotAllocation/ram_use.png', width = 10, height = 8, res = 300, units = 'in')
print(gg.ram)
dev.off()

## try a linear regression model
mod1 <- lm(data = train.dat, ram_gb ~ pixels + samples)
summary(mod1)

## try a linear regression model with interaction
mod2 <- lm(data = train.dat, ram_gb~ pixels + samples + pixels * samples)
summary(mod2)

## try a linear regression model with interaction and log(ram)
mod3 <- lm(data = train.dat, log(ram_gb)~ pixels + samples + pixels * samples)
summary(mod3)


## try to determine some better transformations of our primary variables

#############
## try ACE ##
#############
require(acepack)
ace.res <- ace(x = train.dat[, .(pixels, samples)], y = log(train.dat$ram_gb))

## plot some stuff from transforms
par(mfrow = c(2, 3))
plot(mod2$fitted.values, train.dat$ram_gb, main = 'linear regression (ram=samples:pixels)') ## untransformed inputs
plot(rowSums(ace.res$tx), ace.res$y, main = 'linear regression with transforms') ## transformed inputs
plot(1, 1, pch = '', xlim = c(0, 1), ylim = c(0, 1));abline(a = 0, b = 1);abline(a = 1, b = -1)
plot(ace.res$x[1, ], ace.res$tx[, 1], main = 'pixels transform', xlab = 'orig', ylab = 'tranform(orig)' )
plot(ace.res$x[2, ], ace.res$tx[, 2], main = 'samples transform', xlab = 'orig', ylab = 'tranform(orig)' )
plot(ace.res$y, ace.res$ty, main = 'log(ram) transform', xlab = 'orig', ylab = 'tranform(orig)' )


##############
## try AVAS ##
##############
require(acepack)
avas.res <- avas(x = train.dat[, .(pixels, samples)], y = log(train.dat$ram_gb))

library(MASS)
library(Hmisc)

ace.r <- areg.boot(log_ram_gb~ pixels + samples, B = 100, data = train.dat)
f     <- Function(ace.r, ytype='inverse')
pred  <- f$log_ram_gb(predict(ace.r,train.dat))

## plot some stuff from transforms
## png('~/Desktop/SlotAllocation/ram_ace.png', width = 12, height = 8, res = 300, units = 'in')
par(mfrow = c(2, 3))
plot(mod2$fitted.values, train.dat$ram_gb, main = 'linear regression (ram=samples:pixels)') ## untransformed inputs
plot(rowSums(avas.res$tx), avas.res$y, main = 'linear regression with transforms')    ## transformed inputs
points(rowSums(avas.res$tx), pred, col = "red")
points(predict(ace.r, train.dat), pred, col = "blue")
plot(1, 1, pch = '', xlim = c(0, 1), ylim = c(0, 1));abline(a = 0, b = 1);abline(a = 1, b = -1)
plot(avas.res$x[1, ], avas.res$tx[, 1], main = 'pixels transform', xlab = 'orig', ylab = 'tranform(orig)' )
plot(avas.res$x[2, ], avas.res$tx[, 2], main = 'samples transform', xlab = 'orig', ylab = 'tranform(orig)' )
plot(avas.res$y, avas.res$ty, main = 'log(ram) transform', xlab = 'orig', ylab = 'tranform(orig)' )
## dev.off()


## ##################################################################
## now we can compare predicted RAM use against slots! see how much
## we'd be saving
## ##################################################################

## LOAD IN THE PROCESSED RUNS TO ASSESS HOW WELL OUR PREDICTIONS COMPARE TO SLOT REQUESTS

## the oct 22 date also added some slots and things, so we can see how use compared to requests
## and we can see how predicted use comapares to requests
test.date <- '20181022'
test.dat <- readRDS(sprintf('/Users/azimmer/Desktop/SlotAllocation/profiled_mem_data_%s.RDS', test.date))
test.dat[, log_ram_gb := log(ram_gb)]

## convert pixels into pixel-years
num.yrs  <- unlist(lapply(sapply(parse(text = test.dat[, year_list]), eval), length))
test.dat[, pixels := pixels * num.yrs]

## first, figure out how much RAM was requested
slot.ram <- numeric(nrow(train.dat))

lbd.machines <- which(grepl('lbd-cluster', test.dat$hostname))
## first 20 are at 17.9 GB/slot. 21-48 are at 12 GB/slot
lbd.num <- as.numeric(unlist(lapply(strsplit(test.dat$hostname[lbd.machines], "-"), function(x){ substr(x[3], 2, 3)})))
lbd.first.20 <- lbd.machines[which(lbd.num <= 20)]
lbd.after.20 <- lbd.machines[which(lbd.num > 20)]
slot.ram[lbd.first.20] <- 17.9
slot.ram[lbd.after.20]  <- 12
## geos machines are same as lbd-cluster{1..20}
slot.ram[which(grepl('geos-app', test.dat$hostname))] <- 17.9

## now, onto the c2 nodes
slot.ram[which(grepl('c2', test.dat$hostname))] <- 9.14

## all the other cn nodes are at 2.5, 3, or 4.
slot.ram[which(grepl('cn', test.dat$hostname))] <- 2.5

test.dat[, slot.ram := slot.ram]
test.dat[, ram.req := slot.ram * slots]



## we want to add in some buffer on based on residuals of the avas transforms
## since we want to conservatively select slots
res.sd <- sd(ace.r$residuals)

## now we can predict RAM use in the test data and compare against requested RAM
test.dat[, ram.pred0sd := exp(f$log_ram_gb(predict(ace.r,train.dat) + 0 * res.sd))] ## no sd buffer
test.dat[, ram.pred1sd := exp(f$log_ram_gb(predict(ace.r,train.dat) + 1 * res.sd))] ## 1 sd buffer
test.dat[, ram.pred2sd := exp(f$log_ram_gb(predict(ace.r,train.dat) + 2 * res.sd))] ## 2 sd buffer
test.dat[, ram.pred3sd := exp(f$log_ram_gb(predict(ace.r,train.dat) + 3 * res.sd))] ## 3 sd buffer



## MAKE SOME PLOTS TO COMPARE PREDICTIONS, SLOT REQUESTS, and ACTUAL MEM USE


png('~/Desktop/SlotAllocation/prediction1_comparison.png', units = 'in', res = 300, height = 9, width = 12)
par(mfrow = c(4, 3))

## no buffer - i.e. best fit
plot(test.dat[, ram.req], test.dat[, ram_gb] / 1000, main = sprintf("Req vs Used: %.1f%% over", 100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.req])),
     xlab = 'requested', ylab = 'max used');abline(a = 0, b = 1, col='red');legend('topright', legend = c('0 SD buffer'))
plot(test.dat[, ram.pred0sd] / 1000, test.dat[, ram_gb] / 1000, main = sprintf("Pred vs Used: %.1f%% over", 100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.pred0sd] / 1000)),
     xlab = 'predicted', ylab = 'max used');abline(a = 0, b = 1, col='red');legend('topright', legend = c('0 SD buffer'))
plot(test.dat[, ram.req], test.dat[, ram.pred0sd] / 1000, main = sprintf("Req vs Pred: %.1f%% over", 100 * mean(test.dat[, ram.pred0sd] / 1000 > test.dat[, ram.req])),
     xlab = 'requested', ylab = 'predicted');abline(a = 0, b = 1, col='red');legend('topright', legend = c('0 SD buffer'))

## 1 sd buffer
plot(test.dat[, ram.req], test.dat[, ram_gb] / 1000, main = sprintf("Req vs Used: %.1f%% over", 100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.req])),
     xlab = 'requested', ylab = 'max used');abline(a = 0, b = 1, col='red');legend('topright', legend = c('1 SD buffer'))
plot(test.dat[, ram.pred1sd] / 1000, test.dat[, ram_gb] / 1000, main = sprintf("Pred vs Used: %.1f%% over", 100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.pred1sd] / 1000)),
     xlab = 'predicted', ylab = 'max used');abline(a = 0, b = 1, col='red');legend('topright', legend = c('1 SD buffer'))
plot(test.dat[, ram.req], test.dat[, ram.pred1sd] / 1000, main = sprintf("Req vs Pred: %.1f%% over", 100 * mean(test.dat[, ram.pred1sd] / 1000 > test.dat[, ram.req])),
     xlab = 'requested', ylab = 'predicted');abline(a = 0, b = 1, col='red');legend('topright', legend = c('1 SD buffer'))

## 2 sd buffer
plot(test.dat[, ram.req], test.dat[, ram_gb] / 1000, main = sprintf("Req vs Used: %.1f%% over", 100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.req])),
     xlab = 'requested', ylab = 'max used');abline(a = 0, b = 1, col='red');legend('topright', legend = c('2 SD buffer'))
plot(test.dat[, ram.pred2sd] / 1000, test.dat[, ram_gb] / 1000, main = sprintf("Pred vs Used: %.1f%% over", 100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.pred2sd] / 1000)),
     xlab = 'predicted', ylab = 'max used');abline(a = 0, b = 1, col='red');legend('topright', legend = c('2 SD buffer'))
plot(test.dat[, ram.req], test.dat[, ram.pred2sd] / 1000, main = sprintf("Req vs Pred: %.1f%% over", 100 * mean(test.dat[, ram.pred2sd] / 1000 > test.dat[, ram.req])),
     xlab = 'requested', ylab = 'predicted');abline(a = 0, b = 1, col='red');legend('topright', legend = c('2 SD buffer'))

## 3 sd buffer
plot(test.dat[, ram.req], test.dat[, ram_gb] / 1000, main = sprintf("Req vs Used: %.1f%% over", 100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.req])),
     xlab = 'requested', ylab = 'max used');abline(a = 0, b = 1, col='red');legend('topright', legend = c('3 SD buffer'))
plot(test.dat[, ram.pred3sd] / 1000, test.dat[, ram_gb] / 1000, main = sprintf("Pred vs Used: %.1f%% over", 100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.pred3sd] / 1000)),
     xlab = 'predicted', ylab = 'max used');abline(a = 0, b = 1, col='red');legend('topright', legend = c('3 SD buffer'))
plot(test.dat[, ram.req], test.dat[, ram.pred3sd] / 1000, main = sprintf("Req vs Pred: %.1f%% over", 100 * mean(test.dat[, ram.pred3sd] / 1000 > test.dat[, ram.req])),
     xlab = 'requested', ylab = 'predicted');abline(a = 0, b = 1, col='red');legend('topright', legend = c('3 SD buffer'))

dev.off()

## save relevant objects for use in mbg_central
percent.over <- c(100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.pred0sd] / 1000),
                  100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.pred1sd] / 1000),
                  100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.pred2sd] / 1000),
                  100 * mean(test.dat[, ram_gb] / 1000 > test.dat[, ram.pred3sd] / 1000))
save(res.sd, percent.over, f, ace.r, file = '~/Documents/GitRepos/lbd_core/mbg_central/util/slot_estimation_fits.RData')
