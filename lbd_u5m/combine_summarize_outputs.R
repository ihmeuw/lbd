## ####################################################################################################
## ####################################################################################################
## AUTHOR:    Roy Burstein
## DATE:      2018
## PURPOSE:   Run after postestimation to merge together rasters, and to summarize and 
##            append admin aggregations. Save output into format lucas expects. 
## ADDL NOTE: follwoing Lucas's standard conventions for naming outputs: 
##            https://hub.ihme.washington.edu/pages/viewpage.action?spaceKey=GT&title=MBG+Output+Maps
## ####################################################################################################
## ####################################################################################################

doplotting <- TRUE

if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

mask <- raster('<<<< FILEPATH REDACTED >>>>')
mask[is.na(mask)] <- 20
mask[mask==1] <- NA
mask[mask==20] <-1


# Drop GUF and ESH for now
drop_codes <- c(69,94)

# Create country mask
all_adm0 <- unique(unlist(lapply(Regions, function(x) suppressMessages(get_adm0_codes(x)))))
all_adm0 <- all_adm0[not(all_adm0 %in% drop_codes)]
# Load adm0 shapefile and rasterize
all_sf <- sf::st_read(get_admin_shapefile(version=modeling_shapefile_version))
all_sf <- all_sf[not(all_sf$ADM0_CODE %in% drop_codes),]
country_mask <- fasterize::fasterize(sf=all_sf, field='ADM0_CODE', raster=mask)
country_mask[ !is.na(country_mask) ] <- 1

# load shapefile
message('loading shapefile for plotting')
sad0 <- shapefile('<<<< FILEPATH REDACTED >>>>')


# loop over age groups
for(group in c('under5','infant','neonatal')){
  
  message(group)
  
  # set directory
  dirpath <- '<<<< FILEPATH REDACTED >>>>'
  
  # loop over raked and unraked
  for(rake in c('raked','unraked')){
    
    message(sprintf(' ... %s', rake))
    
    # loop over summstats
    message(' ... ... rasters ')
    for(stat in c('mean','lower','upper')){
      
      if(rake == 'raked' | (rake == 'unraked' & stat == 'mean')){
      
        message(sprintf(' ... ... ... %s', stat))
        
        # load all raster files for this rake / sumstat combination
        rastl <- list()
        for(f in list.files(path = dirpath, pattern = sprintf('%s_%s_%s_raster.tif',stat,rake,group)))
          rastl[[f]] <- brick(sprintf('%s/%s', dirpath ,f))
        
        # merge the rasters and save them
        if(length(rastl)>1){
          r <- do.call(raster::merge, unname(rastl))
        } else {
          r <- rastl[[1]]
        }
        # Mask to country mask
        country_mask_temp <- crop(country_mask, r)
        r <- mask(r, country_mask_temp)

        writeRaster(r, format='GTiff', overwrite = TRUE,
          file = sprintf('%s/died_%s_%s_%s_%i_%i.tif', dirpath, group, stat, rake, min(year_list), max(year_list)))
          
        ### do a nice lucas plot
        if(doplotting){
          pdf(sprintf('%s/died_%s_%s_%s_%i.pdf', dirpath, group, stat, rake, max(year_list)), height=8, width=12)
            rr <- r[[(length(year_list)-1)]]*mask*1000
            rr[rr>200] <- 200
            # define color breaks
            col.f1 <- colorRampPalette(c("#e58bba", "#f2e8b5"))
            col.f2 <- colorRampPalette(c("#f2e8b5", "#ed152e"))
            if(TRUE==FALSE){ #group == 'under5'){
              breaks <- c(0, 25,
                          26:50,
                          51:200,
                          max(c(200,max(as.vector(rr)))))
              col   <- c("#74039E",
                         col.f1(25),
                         col.f2(150),
                         "#ED152E")
              arg <- list(at=c(25,50,200), labels=c('<=25','50','200+'))
            } else {
              breaks <- c(0, 12,
                          13:25,
                          26:200,
                          max(c(200,max(as.vector(rr)))))
              col   <- c("#74039E",
                         col.f1(13),
                         col.f2(175),
                         "#ED152E")
              arg <- list(at=c(12,25,200), labels=c('<=12','25','200+'))
            }
            raster::plot(rr, 
                         axes      = FALSE,
                         breaks    = breaks,
                         col       = col,
                         maxpixels = length(rr),
                         legend    = FALSE,
                         main      = max(year_list))
            box(col='white')
            lines(sad0, lwd = 0.8)
            plot(rr, legend.only=TRUE, col=col, breaks= breaks, legend.width=1, legend.shrink=0.75,
                 smallplot = c(0.03,0.07, 0.3, 0.5),  axis.args = arg)
          dev.off()
          pdf(sprintf('%s/died_allyrs_%s_%s_%s_%i_%i.pdf', dirpath, group, stat, rake, min(year_list), max(year_list)), height=8, width=12)
            rr <- r*1000
            names(rr) <- year_list
            mn <- min(as.vector(rr),na.rm=TRUE)
            mx <- max(as.vector(rr),na.rm=TRUE)
            plot(rr, 
                 axes      = FALSE,
                 zlim      = c(mn,mx),
                 maxpixels = 1e7,
                 legend    = TRUE)
          dev.off()
        } # CLOSE do plotting
      } # CLOSE conidtional on rake
      
      # Do death rasters as well
      if(rake == 'raked'){
        
        message(' ... ... ... ... death counts')
        
        # load all raster files for this rake / sumstat combination
        rastl <- list()
        for(f in list.files(path = dirpath, pattern = sprintf('fr_%s_%s_deaths_raster',group,stat)))
          rastl[[f]] <- brick(sprintf('%s/%s', dirpath ,f))
        # merge the rasters and save them
        if(length(rastl)>1){
          r <- do.call(raster::merge, unname(rastl))
        } else {
          r <- rastl[[1]]
        }
        # Mask to country mask
        country_mask_temp <- crop(country_mask, r)
        r <- mask(r, country_mask_temp)

        writeRaster(r, format='GTiff', overwrite = TRUE,
                    file = sprintf('%s/died_%s_deathcounts_%s_%s_%i_%i.tif', dirpath, group, stat, rake, min(year_list), max(year_list)))
        
      } # CLOSE raked condtional for deaths rasters
      
      
    } # CLOSE stat loop
    
    # aggregate admin-level aggregate estimates
    for(ad in 0:2){
        
        message(sprintf(' ... ... admin %i ', ad))
        

        # load all admin draws files
        adml <- list()
        for(f in list.files(path = dirpath, pattern = sprintf('draws_%s_ad%i.csv',rake, ad))){
          adml[[f]] <- fread(sprintf('%s/%s', dirpath ,f))
          adml[[f]][, V1 := NULL] # bad first column from a csv, remove it
          drawcols <- colnames(adml[[f]])[grep('V',colnames(adml[[f]]))] # get column names
          newdrawcols <- paste0('V',1:length(drawcols)) # ensure draws are correctly named
          setnames(adml[[f]], drawcols, newdrawcols) # ensure draws are correctly named
          adml[[f]][, names(adml[[f]]) := lapply(.SD, as.numeric)]
          setcolorder(adml[[f]], c('name','year',newdrawcols))
        } 
        a <- rbindlist(adml, fill=TRUE, use.names=TRUE)
        
        # summarize
        drawcols <- colnames(a)[grep('V',colnames(a))] # get column names
        a[, names(a) := lapply(.SD, as.numeric)]
        aci <- t(a[,apply(.SD,1,quantile,probs=c(0.025,0.975),na.rm=TRUE),.SDcols=drawcols]) # Note this NA.RM=TRUE wont break the code but can lead to NA admins. 
        am  <- a[,apply(.SD,1,mean,na.rm=TRUE),.SDcols=drawcols]
        
        # clean up summary
        af <- setDT(cbind(a[,c('name','year'),with=FALSE],am,aci))
        names(af) <- c(sprintf('ADM%i_CODE',ad),'year','mean','lower','upper')
        # Drop particular countries
        ad_table <- foreign::read.dbf(get_admin_shapefile(
          version=modeling_shapefile_version,admin_level=ad,suffix='.dbf'
        ))
        drop_admins <- ad_table[ ad_table$ADM0_CODE %in% drop_codes, sprintf('ADM%i_CODE',ad)]
        keep_rows <- not(af[[sprintf('ADM%i_CODE',ad)]] %in% drop_admins)
        af <- af[keep_rows==TRUE,]
        
        # save a version for lucas 
        for(stat in c('mean','lower','upper')){
          tmp        <- af[,c(sprintf('ADM%i_CODE',ad),'year',stat),with=FALSE]
          names(tmp) <- c(sprintf('ADM%i_CODE',ad),'year','value')
          write.csv(tmp, file = sprintf('%s/died_%s_%s_%s_ad%i.csv', dirpath, group, stat, rake, ad))
        } # CLOSE stat loop
        
        # save a convenience version for me, with mean and ci
        # TODO, for this one load in admin names
        write.csv(af, file = sprintf('%s/died_%s_%s_ad%i_fullsummary.csv', dirpath, group, rake, ad))
        
    #  } # CLOSE conditional  !(ad==0 & rake=='unraked')
      
      # do deaths as well
      message(' ... ... ... DEATH COUNTS ')
      
      # load all admin draws files
      if(rake == 'raked'){
        adml <- list()
        for(f in list.files(path = dirpath, pattern = sprintf('draws_ADM%i', ad))){
          adml[[f]] <- readRDS(sprintf('%s/%s', dirpath ,f))
          drawcols <- colnames(adml[[f]])[grep('V',colnames(adml[[f]]))] # get column names
          adml[[f]][, names(adml[[f]]) := lapply(.SD, as.numeric)]
          keepcols <-c(sprintf('ADM%i_CODE',ad),'year',drawcols)
          adml[[f]] <- adml[[f]][, keepcols, with=FALSE]
          setcolorder(adml[[f]], keepcols)
        } 
        a <- rbindlist(adml, fill=TRUE, use.names=TRUE)
        
        # summarize
        drawcols <- colnames(a)[grep('V',colnames(a))] # get column names
        a[, names(a) := lapply(.SD, as.numeric)]
        aci <- t(a[,apply(.SD,1,quantile,probs=c(0.025,0.975),na.rm=TRUE),.SDcols=drawcols]) # Note this NA.RM=TRUE wont break the code but can lead to NA admins. 
        am  <- a[,apply(.SD,1,mean,na.rm=TRUE),.SDcols=drawcols]
        
        # clean up summary
        af <- setDT(cbind(a[,c(sprintf('ADM%i_CODE',ad),'year'),with=FALSE],am,aci))
        names(af) <- c(sprintf('ADM%i_CODE',ad),'year','mean','lower','upper')
        
        ad_table <- foreign::read.dbf(get_admin_shapefile(
          version=modeling_shapefile_version,admin_level=ad,suffix='.dbf'
        ))
        drop_admins <- ad_table[ ad_table$ADM0_CODE %in% drop_codes, sprintf('ADM%i_CODE',ad)]
        keep_rows <- not(af[[sprintf('ADM%i_CODE',ad)]] %in% drop_admins)
        af <- af[keep_rows==TRUE,]

        # save a version for lucas 
        for(stat in c('mean','lower','upper')){
          tmp        <- af[,c(sprintf('ADM%i_CODE',ad),'year',stat),with=FALSE]
          names(tmp) <- c(sprintf('ADM%i_CODE',ad),'year','value')
          write.csv(tmp, file = sprintf('%s/died_%s_deathcounts_%s_%s_ad%i.csv', dirpath, group, stat, rake, ad))
        } # CLOSE stat loop
        
        # save a convenience version for me, with mean and ci
        # TODO, for this one load in admin names
        write.csv(af, file = sprintf('%s/died_%s_deathcounts_%s_ad%i_fullsummary.csv', dirpath, group, rake, ad))
      } # CLOSE raked only for deaths (so we do it only)
        
      
    } # CLOSE ad loop
  } # CLOSE rake loop
} # CLOSE group loop


# NOTE potential reasons for larger ci than last time: nugget, data in annual vs 5yr, tmb not inla?, no stackers

