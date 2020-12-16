# HEADER ------------------------------------------------------------------
# Project: MBG
# Purpose: various R functions for use in vaccine data processing
# source("FILEPATH")

# FUNCTIONS ---------------------------------------------------------------

# Timer Functions ---------------------------------------------------------
#
## Functions to time within the child stacking regions
##
## General usage:
##    require(tictoc)
## 
##    tic("Step 1")
##    **your code here**
##    toc(log = T)
##
##    ticlog <- tic.log(format = F)
##    generate_time_log(ticlog)   
##
##  Returns: data table with two columns
##     "step": names of events (e.g. "Step 1")
##     "time": time elapsed (as text: Xh Xm Xs)
##
##  Note: can nest tic/toc pairs

generate_time_log <- function(ticlog) {
  
  # Set up packages
  require(magrittr)
  require(data.table)
  
  # Functions in functions
  strip_time <- function(x) {
    sec <- as.numeric(x$toc - x$tic)
    time <- format_time(sec)
    name <- x$msg
    
    df <- c(name, time) %>%
            t %>%
            as.data.table
    
    names(df) <- c("step", "time")
    
    return(df)
  }
  
  format_time <- function(run_time) {
    run_time <- round(as.numeric(run_time),0)
    
    hours <- run_time %/% 3600
    remainder <- run_time %% 3600
    minutes <- remainder %/% 60 
    seconds <- remainder %% 60
    
    run_time <- paste0(hours, "h ", minutes, "m ", seconds, "s")
    return(run_time)
  }
  
  df_out <- lapply(ticlog, strip_time) %>% rbindlist
  
  return(df_out)

}


# gaul_convert() ----------------------------------------------------------
# Now in prep_functions.R 
  
# combine_country_image_history -------------------------------------------

## Combine input data and covariates layers for models run by country
## Necessary for saving to Shiny tool
## Save "covs", "tv_*", and "df" to new combined snapshot in model_image_history

combine_country_image_history <- function(indicator, indicator_group, run_date, fixed_effects, countries) {

  # Combine non-varying covs
  load(paste0("FILEPATH"))
  nt_covs <- names(covs)
  combine_nt_cov <- function(nt_cov) {
    pull_nt_covs <- function(country) {
      load(paste0("FILEPATH"))
      cov_layer <- covs[[nt_cov]]
      return(cov_layer)
    }
    country_layers <- lapply(countries, pull_nt_covs)
    if (length(country_layers) > 1) {
        combined_layers <- do.call(raster::merge, country_layers)
      } else if (length(country_layers) == 1) {
        combined_layers <- country_layers
      }
      names(combined_layers) <- nt_cov
      return(combined_layers)
  }

  covs <- lapply(nt_covs, combine_nt_cov)
  covs <- do.call(raster::brick, covs)

  # Combine varying covs
  selected_covs <- strsplit(fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]
  for(c in selected_covs) {
    if(paste0('tv_',c) %in% grep('tv_*', ls(), value = TRUE)) {
      pull_tv_covs <- function(country) {
        load(paste0("FILEPATH"))
        tv_cov <- get(paste0('tv_',c))
        return(tv_cov)
      }
      country_layers <- lapply(countries, pull_tv_covs)
      
      if (length(country_layers) > 1) {
            combined_layers <- do.call(raster::merge, country_layers)
         } else if (length(country_layers) == 1) {
            combined_layers <- country_layers
         }

      names(combined_layers) <- gsub("layer", c, names(combined_layers))
      assign(paste0('tv_', c), combined_layers)
    }
  }

  # Combine input data
  pull_df <- function(region) {
     load(paste0("FILEPATH"))
    return(df)
  }
  df <- lapply(countries, pull_df)
  df <- do.call(rbind.fill, df)

  save(list = c('df','covs',grep('^tv_*', ls(), value = TRUE)), file = paste0("FILEPATH"))

}

# pushover_notify() -------------------------------------------------------

# Function to notify via pushover when job done on cluster

pushover_notify <- function(message, title = NULL, image = NULL) {
  
  ## Function to enable custom push notifications from within R on cluster
  ## Uses command line curl call 
  ## (the httr package would be a good alternative, but breaks on our cluster)
  
  ## Inputs: message (text you want to send) 
  #          title (optional)
  #          image (path to image to include)
  
  ## Outputs: command line curl call
  
  ## Setup: ---------------------------------------------------------------
  
  ## 1) Download the pushover app on your phone and sign up for an account
  ##    https://pushover.net
  
  ## 2) Find your user token and enter in the "keyfile" 
  ##    in your home directory (h:/lib/pushover_key.txt) 
  ##    (plain text, token only, no spaces)
  
  ## 3) Register an application (https://pushover.net/apps/build) with pushover
  ##    and find your app token, then enter that app token in the 'appfile'
  ##    in your home directory (h:/lib/pushover_app.txt)
  ##    (again, plain text, no spaces)
  
  ## 4) Enjoy custom push notifications from your R code with the syntax
  ##    `pushover_notify("my_message", "my_title")`.  You can even send images
  ##    if you point towards a saved .jpg or .png on the cluster!
  
  # -----------------------------------------------------------------------
  
  if (Sys.info()["sysname"] == "Linux") {  
    h_root <- paste0("FILEPATH")
  } else { 
    h_root <- "FILEPATH"
  }
  
  url <- "https://api.pushover.net/1/messages.json"

  keyfile <- paste0("FILEPATH")
  
  # Check if file exists to continue
  if (!file.exists(keyfile)) {
    return("Pushover not enabled")
  }
  
  key <- readChar(keyfile, file.info(keyfile)$size)
  
  appfile <- paste0("FILEPATH")
  app <- readChar(appfile, file.info(appfile)$size)
  
  if (is.null(image)) {
    img_addin <- "" 
  } else if (!is.null(image)) {
    img_addin <- paste0("-F \"attachment=@", image, "\" \\\n")
  }

  curl_cmd <- paste0("curl -s ",
                      "--form-string \"token=", app, "\" \\\n", 
                      "--form-string \"user=", key, "\" \\\n",
                      "--form-string \"message=", message, " \" \\\n",
                     "--form-string \"title=", title, " \" \\\n",
                     img_addin,
                      "https://api.pushover.net/1/messages.json")
  
  system(curl_cmd, ignore.stdout = T)

}


## Fork of combine_region_image_history
## Modified to work with INLA stacking code
## Combine input data and covariates layers for models run by region
## Necessary for saving to Shiny tool
## Save "covs", "tv_*", and "df" to new combined snapshot in model_image_history
combine_region_image_history_inlastack <- function(indicator, indicator_group, run_date, fixed_effects) {

  # Combine varying covs
  selected_covs <- strsplit(all_fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]
  for(c in selected_covs) {
    if(c %in% names(cov_list)) {
      pull_tv_covs <- function(region) {
        load(paste0("FILEPATH"))
        tv_cov <- cov_list[[c]]
        return(tv_cov)
      }
      region_layers <- lapply(Regions, pull_tv_covs)
      combined_layers <- do.call(raster::merge, region_layers)
      names(combined_layers) <- gsub("layer", c, names(combined_layers))
      assign(paste0('tv_', c), combined_layers)
    }
  }

  # Combine input data
  pull_df <- function(region) {
    load(paste0("FILEPATH"))
    return(df)
  }
  df <- lapply(Regions, pull_df)
  df <- do.call(rbind.fill, df)

  covs <- NULL  # No ntv covs

  save(list = c('df','covs',grep('^tv_*', ls(), value = TRUE)), file = paste0("FILEPATH"))

}


shiny_cov_layers_inla <- function(fixed_effects, gaul_list, run_date, indicator, indicator_group) {
  
  # Plot covariate layers
  library(ggplot2, lib.loc = package_lib)
  
  color_list <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
  
  # Make sure covs are loaded
  load(paste0("FILEPATH"))
  
  # Load actual data (df already in memory)
  if(time_stamp==TRUE) output_dir <- paste0("FILEPATH")
  if(time_stamp==FALSE) output_dir <- paste0("FILEPATH")
  plot_dir <- paste0(output_dir,"/plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Get templates
  if(exists("subset_shape")==FALSE) {
    message("Opening master shapefile because not found in global env...")
    master_shape <- shapefile(paste0("FILEPATH"))
    subset_shape <- master_shape[master_shape@data$GAUL_CODE %in% gaul_list, ]   
  }
  admin0.dt <- data.table(fortify(subset_shape)) 
  
  # Plot varying covariates
  gLegend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  selected_covs <- strsplit(fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]
  library(grid)
  tv_count <- 0
  for(c in selected_covs) {
    if(paste0('tv_',c) %in% grep('tv_*', ls(), value = TRUE)) {
      
      tv_count <- tv_count + 1
      tv_cov <- get(paste0('tv_',c))
      
      # Convert raster to SpatialPointsDataFrame
      preds.sp <- rasterToPoints(tv_cov, spatial=TRUE)
      projection <- proj4string(preds.sp)
      
      # reproject sp object  
      preds.sp <- spTransform(preds.sp, CRS(projection)) 
      preds.sp@data <- data.frame(preds.sp@data, long=coordinates(preds.sp)[,1],lat=coordinates(preds.sp)[,2]) 
      preds.dt <- data.table(preds.sp@data)
      
      ## Plot preds of proportion with 0 years of education
      names(preds.dt)[names(preds.dt) == "lat"] = "latitude"
      names(preds.dt)[names(preds.dt) == "long"] = "longitude" 
      
      # Plot predictions for all periods
      plot.preds <- function(x) {
        period <- gsub(paste0(c,'.'), "", x)
        loop.preds.gg <- ggplot(preds.dt,aes(longitude,latitude)) +
          geom_raster(aes(fill=get(x))) +
          coord_fixed() + 
          theme_minimal() +
          geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='white', lwd=.1) +
          scale_fill_gradientn(colours=rev(color_list), limits=c(min(minValue(tv_cov)), max(maxValue(tv_cov))), na.value = "grey") + 
          guides(fill=guide_colorbar(title=c, label=TRUE, ticks=FALSE)) +
          scale_x_continuous("", breaks=NULL) +
          scale_y_continuous("", breaks=NULL) +
          theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines")) + 
          ggtitle(paste0("Period ",period))
        return(loop.preds.gg)
      }
      for(i.period in 1:4) {
        assign(paste0("tv_cov_",tv_count,".gg.", i.period),plot.preds(paste0(c,'.',i.period)))
      }
    }    
  }
  
  # Plot all varying covariates
  png(paste0("FILEPATH"),width=1600,height=1600)
  total_rows <- tv_count*5
  # Initialize plot with master title
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(total_rows, 22, heights=c(.25,.25,.25,.25), widths=c(.25,.25,.25,.25))))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  # Plot all data coverage maps
  start_row <- 1
  for(i in 1:tv_count) {
    end_row <- start_row + 4
    print(get(paste0('tv_cov_',i,'.gg.1')) + theme(legend.position="none"), vp = vplayout(start_row:end_row, 1:5))
    print(get(paste0('tv_cov_',i,'.gg.2')) + theme(legend.position="none"), vp = vplayout(start_row:end_row, 6:10))
    print(get(paste0('tv_cov_',i,'.gg.3')) + theme(legend.position="none"), vp = vplayout(start_row:end_row, 11:15))
    print(get(paste0('tv_cov_',i,'.gg.4')) + theme(legend.position="none"), vp = vplayout(start_row:end_row, 16:20))
    p.legend <- gLegend(get(paste0('tv_cov_',i,'.gg.1')))
    p.legend$vp <- viewport(layout.pos.row = start_row:end_row, layout.pos.col = 21:22)
    grid.draw(p.legend)
    start_row <- start_row + 5
  }
  dev.off()
  
}

# Wrapper for easy plotting of run
plot_run_raster <- function(raster, 
                           vax_title,
                           summary_stat,
                           stat_title,
                           indicator = indicator, 
                           indicator_group = indicator_group, 
                           run_date = run_date,
                           layer_names = c("2000", "2005", "2010", "2015")) {

  out_dir <- paste0("FILEPATH")
  dir.create(out_dir)

  plot_raster(ras = raster,
              indicator = paste0(indicator, "_", summary_stat),
              region = 'africa', 
              return_map = F,
              out_dir = out_dir,
              highisbad = F,
              min_value = 0,
              mid_value = 0.5,
              max_value = 1,
              legend_title = paste0(vax_title, "\n", stat_title),
              plot_title = paste0(vax_title, " ", stat_title),
              layer_names = as.character(layer_names),
              cores = 16,
              individual_layers = T)

}

# Create raster of GBD estimates
make_gbd_rasterbrick <- function(gbd_data, 
                                 indicator, 
                                 gaul_list,
                                 years = c(2000, 2005, 2010, 2015)) {

  message(paste0("Creating GBD RasterBrick for ", indicator))

  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                        buffer = 0.4, shapefile_version=modeling_shapefile_version)

  subset_shape   <- simple_polygon_list[[1]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]

  # create template raster
  ext <- extent(simple_raster)
  ncol <- ncol(simple_raster)
  nrow <- nrow(simple_raster)
  r <- raster(ext, ncol, nrow)

  gbd_brick_list <- lapply(years, function(yr) {

    my_spdf <- merge(subset_shape, gbd_data[year == yr], by.x = "GAUL_CODE", by.y = "name", all.x = T, all.y = F)
    ras <- rasterize(my_spdf, r, "mean")

  })

  gbd_brick <- brick(gbd_brick_list)

  return(gbd_brick) 
  
}

make_rf_plot <- function(indicator_group, indicator, run_date, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking) {

  if (Sys.info()["sysname"] == "Linux") {
    j_root <- "FILEPATH"
  } else {
    j_root <- "FILEPATH"
  }

  require(data.table)
  require(ggplot2)
  require(magrittr)

  in_dir <- paste0("FILEPATH")
  out_dir <- paste0("FILEPATH")
  rf_file <- paste0("FILEPATH")

  table_file <- paste0("FILEPATH")
  gaul_table <- read.csv(table_file) %>% data.table
  gaul_table <- subset(gaul_table, select = c("short_name", "gaul"))

  # fix Sudan
  gaul_table[short_name == "Sudan", gaul := 6 ]

  rf <- read.csv(rf_file, stringsAsFactors = F) %>% as.data.table
  rf[, X := NULL]
  setnames(rf, "name", "gaul")

  rf <- merge(rf, gaul_table, all.x = T)

  regions <- c("cssa", "essa", "sssa", "wssa", "name")

  for (reg in regions) {
    
  rf_reg <- rf[gaul %in% get_adm0_codes(reg, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)]

  gg <- ggplot(rf_reg, aes(x = rake_to_mean, y = geo_mean)) +
    geom_point(aes(color = year)) +
    xlim(c(0,1)) + ylim(c(0,1)) +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face="bold", size=10)) +
    labs(x = "GBD predicted", y = "MBG predicted", color = "Year") +
    facet_wrap(~short_name)

  png(filename = paste0(out_dir, "rf_scatter_", reg, ".png"),
      width = 10, 
      height = 6,
      units = "in", 
      res = 300,
      type = "cairo")
    
  print(gg)
    
  dev.off()

  }

}


transform_spatial_effects <- function(params.mat,age){

    nr=nrow(params.mat)

    tau_mean =exp(params.mat[(nr-2),1])
    kappa_mean=exp(params.mat[(nr-1),1])

    tau_lb =exp(params.mat[(nr-2),3])
    kappa_lb=exp(params.mat[(nr-1),3])

    tau_ub =exp(params.mat[(nr-2),5])
    kappa_ub=exp(params.mat[(nr-1),5])

    res=data.table('parameter'=c("GPRandom Range in Decimal Degrees","GPRandom Nominal Variance"))
        res[,paste0('a',age,'_mean'):=c(sqrt(8)/kappa_mean,1/(4*pi*kappa_mean^2*tau_mean^2))]
        res[,paste0('a',age,'_lb')  :=c(sqrt(8)/kappa_ub,1/(4*pi*kappa_ub^2*tau_ub^2))]
        res[,paste0('a',age,'_ub')  :=c(sqrt(8)/kappa_lb,1/(4*pi*kappa_lb^2*tau_lb^2))]

    return(res)

  }


clean_model_results <- function(rd   = run_date,
                                regs = Regions,
                                ages = 0,
                                nm   = '',
                                indicator = indicator, 
                                indicator_group = indicator_group){

  require(magrittr)

  sharedir <- paste0("FILEPATH")

  # make loopvars
  lv <- expand.grid(regs,ages)

  # grab model fit objects
  mods <- model_fit_table(lv=lv,rd=rd,nullmodel=nm, indicator = indicator, indicator_group = indicator_group)

  # one table for each region
  tables=list()
  for(r in regs){

    # get parameter names
    tempbig <- lapply(ages, function(age) row.names(mods[[paste0(r, '_', age)]])) %>%
                  unlist %>%
                  unique %>%
                  data.table(parameter = .)

     tempbig = tempbig[!parameter%in%c('Theta1 for space','Theta2 for space','GroupRho for space'),]
     assign(r,tempbig)
      for(a in ages){
          tmp=data.table(cbind(parameter=row.names(mods[[paste0(r,'_',a)]]),
                                         mods[[paste0(r,'_',a)]]))[,c(1,2,4,6)]
          setnames(tmp, c("parameter",paste0("a",a,"_mean"),paste0("a",a,"_lb"),paste0("a",a,"_ub") ))

          # tranform spatial covariates to be more readbale
          tmp$parameter=as.character(tmp$parameter)
          tmp<-tmp[!parameter%in%c('Theta1 for space','Theta2 for space'),]
          tmp$parameter[tmp$parameter=="GroupRho for space"]="GPRandom Rho for time"
          tmp <- rbind(tmp,transform_spatial_effects(params.mat=mods[[paste0(r,'_',a)]],age=a))

          # round
          cols <- names(tmp)[2:4]
          tmp[,(cols) := round(.SD,4), .SDcols=cols]

          assign(r,merge(get(r),tmp,by='parameter',all=TRUE)) #,all.x=TRUE))

      }
      tables[[r]]<-get(r)

  }

  for(r in regs){
    write.csv(tables[[r]],sprintf('%s/model_results_table_%s.csv',sharedir,r))
  }
  return(tables)

}


# Load alternative raking targets
load_newest_gbd_vax  <- function(vaccine,
                                 gaul_list = NULL,
                                 return_cis = F,
                                 return_mode = "national", # national, subnational, or both
                                 gbd_year = 2020,
                                 years = c(2000:2019),
                                 gbd_date = "best", 
                                 return_field = "location_id") {

  source("FILEPATH")

  str_match <- stringr::str_match

  gaul_to_loc_id <- get_location_metadata(location_set_id = 22, gbd_round_id=6)
  names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"

  gbd_file <- paste0("FILEPATH")

  gbd <- readRDS(gbd_file)
  gbd <- merge(gaul_to_loc_id, gbd, by="location_id",  all.y = T)  

  if (return_mode == "national") {
    gbd <- subset(gbd, !grepl('_', ihme_loc_id))
  }

  if (return_field == "GAUL_CODE") {
    loc_map <- get_location_code_mapping(shapefile_version=modeling_shapefile_version)
    gbd <- merge(gbd, subset(loc_map, select = c("loc_id", "GAUL_CODE")), by.x = "location_id", by.y = "loc_id", all.x = T)
  }

  # subset to year & gaul lists
  gbd <- subset(gbd, year_id %in% years)
  if (!is.null(gaul_list)) {
    gbd <- subset(gbd, get(return_field) %in% gaul_list)
  }
  
  setnames(gbd, 
           c(return_field, "year_id", "gpr_mean", "gpr_lower", "gpr_upper"),
           c("name", "year", "value", "lower", "upper"))
  if (return_cis == F)  return_vars <- c("name", "year", "value")
  if (return_cis == T)  return_vars <- c("name", "year", "value", "lower", "upper")

  if (return_mode != "national") {
    gbd[, parent_ihme_loc_id := str_match(ihme_loc_id,"(.*)_")[,2]]

    if (return_mode == "subnational") {

      # Drop national estimates if we have subnational estimates
      ihme_loc_ids_to_drop <- c( NA , "CHN", "GBR", "ZAF", "USA", "MEX", "IDN", "BRA", "IND", "NOR", "SWE", "JPN",
                  "KEN" ,"NZL", "ETH", "IRN", "RUS")
     

      ihme_loc_ids_to_drop <- ihme_loc_ids_to_drop[!is.na(ihme_loc_ids_to_drop)]
      gbd <- subset(gbd, !(ihme_loc_id %in% ihme_loc_ids_to_drop))  
      gbd <- subset(gbd, !grepl('NGA_', ihme_loc_id))
      gbd <- subset(gbd, !grepl('PAK_', ihme_loc_id))
      gbd <- subset(gbd, !grepl('PHL_', ihme_loc_id))
    }
  }

  gbd <- subset(gbd, select = return_vars)

  return(gbd)
}

load_dropout_for_raking <- function(vaccine,
                                    gaul_list,
                                    years = c(2000:2016),
                                    abs_rel,
                                    gbd_date = "best") {

  vax1_gbd <- load_newest_gbd_vax(vaccine = paste0(vaccine, "1"),
                                  gaul_list = gaul_list,
                                  years = years,
                                  gbd_date = gbd_date)

  vax3_gbd <- load_newest_gbd_vax(vaccine = paste0(vaccine, "3"),
                                  gaul_list = gaul_list,
                                  years = years,
                                  gbd_date = gbd_date)

  setnames(vax1_gbd, "mean", "vax1")
  setnames(vax3_gbd, "mean", "vax3")
  vax1_3_dropout_gbd <- merge(vax1_gbd, vax3_gbd)
  if (abs_rel == "absolute") vax1_3_dropout_gbd[, mean := (vax1 - vax3)]
  if (abs_rel == "relative") vax1_3_dropout_gbd[, mean := ((vax1 - vax3)/vax1)]

  vax1_3_dropout_gbd <- subset(vax1_3_dropout_gbd, select = c("name", "year", "mean"))

   if (nrow(vax1_3_dropout_gbd[mean < 0,]) > 0) {
    n_fix <- nrow(vax1_3_dropout_gbd[mean < 0,])
    warning(paste0(n_fix, " country-year pairs have vax3 > vax1. Setting to 1e-4..."))
    message("Affected country-years: ")
    print(vax1_3_dropout_gbd[mean < 0,])
    vax1_3_dropout_gbd[mean < 0, mean := 1e-4]

  }

  return(vax1_3_dropout_gbd)

}


make_country_cell_pred <- function(indicator,
                                   indicator_group,
                                   run_date,
                                   reg,
                                   cell_pred,
                                   shapefile_version=modeling_shapefile_version,
                                   subnational_raking=subnational_raking) {

   # Getting the simple polygon and simple raster objects for this region alone
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  message("Getting the spatial objects associated with this region.")
    simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version=modeling_shapefile_version), buffer = 0.4, subset_only = TRUE, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)
    subset_shape   <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
  message("Building simple raster from subset_shape")
    raster_list    <- build_simple_raster_pop(subset_shape)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]
    rm(raster_list,simple_polygon_list,pop_raster);gc()
  message("All done loading spatial template files (subset_shape,simple_polygon,simple_raster,pop_raster, etc)")

  # Determining a list of the valid pixel indices based on the simple raster template
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  pixel_spatial<-data.table(pixel_id=pixel_id)
  message("Pixel ID that links cell_pred observations to raster locations created from simple raster, checking for sameness in dimensions compared to cell_pred.")
  if(!(length(pixel_id)==nrow(cell_pred)/16)){
    stop("Excuse me, but the number of valid pixels in your simple raster is not the same number of valid pixels in your cell_pred object. Look into this!")
  }else{
    message(paste0("Check passed: There are ",length(pixel_id),
                   " pixels in your simple_raster object, the same number of pixels in the cell_pred object for each year."))
  }

  region_adm0_list<-get_adm0_codes(reg, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking) 

  # Defining a function that will get the raster versions of each Admin level:
  GetAdmin0<-function(simple_raster, region_adm0_list){
    message("Loading admin level 0")
      admin_shp<-readOGR(dsn=paste0("FILEPATH"),
                         layer=paste0("g2015_2014_", admin_level))
             
    message("Rasterizing...")
      admin_rast<-rasterize(admin_shp,simple_raster,"ADM0_CODE")
    message("Converted to raster based on simple_raster template. Cropping and masking:")
      admin_rast  <- crop(admin_rast, extent(simple_raster))
      admin_rast  <- setExtent(admin_rast, simple_raster)
      admin_rast  <- mask(admin_rast, simple_raster)
    message("Subsetting polygon and point objects to only contain the relevant ADM0 codes; calculating centroids.")
      admin_shp<-admin_shp[admin_shp@data$ADM0_CODE %in% region_adm0_list,]
      admin_centroids<-SpatialPointsDataFrame(gCentroid(admin_shp, byid=TRUE), 
                                              admin_shp@data, match.ID=FALSE)
      message("Loaded in.")
    message("Compiling and returning results.")
      admin<-list()
      admin[["spdf"]]<-admin_shp
      admin[["centroids"]]<-admin_centroids
      admin[["rast"]]<-admin_rast
      admin[["attributes"]]<-copy(data.table(admin_shp@data))
    return(admin)
  }
  
  admin_levels <- list()
  message("Rasterizing shapefiles; this may take a while.")
  fieldname<-paste0("ADM0_CODE")
  admin_info<-GetAdmin0(simple_raster,region_adm0_list)
  pixel_spatial[[fieldname]]<-extract(admin_info[["rast"]],pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
  admin_levels[["0"]]<-admin_info # Add the admin info to the list
  if(sum(is.na(pixel_spatial[[fieldname]]))>0){ # Check to see if any of the pixels don't have a location assigned
    message(paste0("   Whoah, there are some pixels that are NA, and have not been assigned a location for level 0"))
    }

 # Merging together cell_pred and spatial information, population information.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cell_pred<-data.table(cell_pred) # Converting cell_pred to a data.table
    draw_colnames<-names(cell_pred)
     
  # Loading and Assigning Pixels to Admin Units (that are missing)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # We now need to figure out what admin units don't show up when we pull the locations of individual raster pixels of our results.
  # To do this, we need to know what Admin 1 and 2 units exist in the shapefiles we want to use, but don't show up in the results.
  
  message("Generating a raster of the pixel_id values")
    pixel_id_raster<-copy(simple_raster)
    for (pix in pixel_id){ # For each of the valid pixels, replace the value with the pixel_id, not whatever country ID is used in the simple_polygon/raster
      pixel_id_raster@data@values[pix]<-pix
    } 
      
  missing_admins<-list() # This will contain the cell_preds for otherwise missing admin units, by level.

  fieldname<-paste0("ADM0_CODE")
  # First, we discover what GAUL codes are missing from the pixel_spatial table compared to the shapefile:
  shpfile_attributes<-admin_levels[["0"]][["attributes"]] # Get the attribute table from the admin level's shapefile
  attribute_codes<-shpfile_attributes[[paste0("ADM0_CODE")]] # Get list of codes for that level from the admin level's attribute table
  admin_cols<-names(shpfile_attributes)[names(shpfile_attributes) == "ADM0_CODE"] # Any columns in ADM0 should be included in merge
  pixel_codes<-pixel_spatial[[fieldname]] # Get list of codes based on what's in the pixel_spatial object
  missing_codes<-attribute_codes[!(attribute_codes %in% pixel_spatial[[fieldname]])] # Get list of missing codes
  missing_table<-shpfile_attributes[shpfile_attributes[[fieldname]]%in%missing_codes,admin_cols,with=F]
  
  if(length(missing_codes)==0){
    message(paste0("No missing codes at level 0"))
  }else{
    message(paste0("Missing codes found for level 0: "))
    # Strategy 1: Assign a pixel location based on centroid location
    # Develop a raster of the pixel-IDs: 
    message("  Discovering centroid location")
      points<-admin_levels[["0"]][["centroids"]] # SpatialPointsDataFrame of that level
      missing_points<-points[(points@data[[fieldname]] %in% missing_codes),] # getting only the missing points
      missing_centroid_locs<-data.table(raster::extract(pixel_id_raster,missing_points)) # Extracting the missing location values as centroid points... 
      names(missing_centroid_locs)<-"point_method"
      missing_admins_centroids<-data.table(missing_codes,missing_centroid_locs)
      
    # Strategy 2: Assign a pixel location based on polygon extraction
    # Develop a raster of the pixel-IDs: 
    message("  Discovering first raster pixel touching polygon")
      polys<-admin_levels[["0"]][["spdf"]] # SpatialPointsDataFrame of that level
      missing_polys<-polys[(polys@data[[fieldname]] %in% missing_codes),] # getting only the missing points
      missing_poly_locs<-data.table(raster::extract(x=pixel_id_raster,y=missing_polys,small=T,fun=function(x,...)first(x))) # Extracting the missing location values as polygons, pulling the first raster cell that it touches... 
      names(missing_poly_locs)<-"poly_method"
      missing_admins_polys<-data.table(missing_codes,missing_poly_locs) # Extracting the 
    
    # Merging strategies together: centroids and polygons; adding to list.
      missing_locs<-merge(missing_admins_polys,missing_admins_centroids,by="missing_codes")
      setnames(missing_locs,"missing_codes",fieldname)
      missing_locs<-merge(missing_locs,missing_table,by=fieldname) # Add in admin 0, 1 levels if relevant
      missing_locs[,pixel_id:=point_method] # If centroid produced something, go with centroid
      missing_locs[is.na(point_method),pixel_id:=poly_method] # Otherwise, go with the polygon method.
    
    # Check to see if NAs still exist in missing locations, make a warning about missing areas.
      if(NA %in% missing_locs$pixel_id){
        message( "The following admin units appear to be getting lost forever:")
        print(missing_locs[is.na(pixel_id),c("pixel_id",admin_cols),with=F])
      }
    
    # Merging on locations with pixel IDs
      missing_locs<-missing_locs[!is.na(pixel_id),c("pixel_id",admin_cols),with=F]
      missing_locs<-merge(missing_locs,cell_pred[,c("pixel_id","pop","year",draw_colnames),with=F],by="pixel_id",allow.cartesian=T)
      missing_admins[["0"]]<-missing_locs
    }
  
}

# make_qsub, adapted to launch vaccine scripts
make_qsub_vaccine <- function(user = 'USERNAME',
                              cores = slots,
                              memory = 100,
                              proj  = 'PROJECT',
                              indic = indicator,  
                              vax = vaccine,                          
                              rd = run_date,
                              log_location = 'sgeoutput',
                              code,
                              shell = "r_shell.sh",
                              geo_nodes = FALSE){

   sharedir <- sprintf("FILEPATH")
   dir.create(sharedir, recursive = T, showWarnings = F)

   if(log_location=='sgeoutput')
     logloc = sprintf("FILEPATH")
   if(log_location=='sharedir'){
     logloc = sprintf('%s/output/%s',sharedir,rd)
     dir.create(sprintf('%s/errors',logloc), recursive = T, showWarnings = F)
     dir.create(sprintf('%s/output',logloc), recursive = T, showWarnings = F)
   }

  ## do we want to submit to geo nodes? if so, there are a few tweaks
  if(geo_nodes == TRUE){
    shell     <- "FIELPATH"   ## for the correct path to R
    proj      <- "PROJECT"    ## correct proj for geos nodes
    node.flag <- "NODES" ## send the job to geo nodes
  }else{
    node.flag <- "" ## don't do anything special
  }

  return(paste0("qsub -e ", logloc, "/errors -o ", logloc, "/output -cwd -l mem_free=", memory, "G -pe multi_slot ", cores, " -P ", proj, " ", node.flag, " -N launch_", indic, " ", shell, " ", code, ".R ", indic, " ", vax, " ", rd))

}

make_qsub_combine_vaccine <- function(user = 'USERNAME',
                                      cores = slots,
                                      memory = 100,
                                      proj  = 'PROJECT',
                                      ig = indicator_group,
                                      vax = vaccine,
                                      doses,
                                      region,                          
                                      rd = run_date,
                                      holdout = 0,
                                      log_location = 'sgeoutput',
                                      code,
                                      shell = "FILEPATH",
                                      geo_nodes = FALSE){

   indic <- paste0(vax, doses, "_cov")
   sharedir <- sprintf("FILEPATH")
   dir.create(sharedir, recursive = T, showWarnings = F)

   if(log_location=='sgeoutput')
     logloc = sprintf("FILEPATH")
   if(log_location=='sharedir'){
     logloc = sprintf('%s/output/%s',sharedir,rd)
     dir.create(sprintf('%s/errors',logloc), recursive = T, showWarnings = F)
     dir.create(sprintf('%s/output',logloc), recursive = T, showWarnings = F)
   }

  ## do we want to submit to geo nodes? if so, there are a few tweaks
  if(geo_nodes == TRUE){
    shell     <- "FILEPATH"  ## for the correct path to R
    proj      <- "PROJECT"    ## correct proj for geos nodes
    node.flag <- "NODES" ## send the job to geo nodes
  }else{
    node.flag <- "" ## don't do anything special
  }

  return(paste0("qsub -e ", logloc, "/errors -o ", logloc, "/output -cwd -l mem_free=", memory, "G -pe multi_slot ", cores, 
                " -P ", proj, " ", node.flag, " -N combine_", vax, "_",region, " ", ig, "/", shell, " ",ig, "/", code, ".R ", 
                vax, " ", doses, " ", rd, " ", region, " ", holdout))

}

make_qsub_mean_raster <- function(user = 'USERNAME',
                                  cores = slots,
                                  memory = 100,
                                  proj  = 'PROJECT',
                                  ig = indicator_group,
                                  vax = vaccine,
                                  doses,
                                  indic = indicator,                          
                                  rd = run_date,
                                  log_location = 'sgeoutput',
                                  code,
                                  shell = "r_shell.sh",
                                  geo_nodes = FALSE){

   sharedir <- sprintf("FILEPATH")
   dir.create(sharedir, recursive = T, showWarnings = F)

   if(log_location=='sgeoutput')
     logloc = sprintf("FILEPATH")
   if(log_location=='sharedir'){
     logloc = sprintf('%s/output/%s',sharedir,rd)
     dir.create(sprintf('%s/errors',logloc), recursive = T, showWarnings = F)
     dir.create(sprintf('%s/output',logloc), recursive = T, showWarnings = F)
   }

  ## do we want to submit to geo nodes? if so, there are a few tweaks
  if(geo_nodes == TRUE){
    shell     <- "FILEPATH"   ## for the correct path to R
    proj      <- "PROJECT"    ## correct proj for geos nodes
    node.flag <- "NODES" ## send the job to geo nodes
  }else{
    node.flag <- "" ## don't do anything special
  }

  return(paste0("qsub -e ", logloc, "/errors -o ", logloc, "/output -cwd -l mem_free=", memory, "G -pe multi_slot ", cores, 
                " -P ", proj, " ", node.flag, " -N mean_ras_", indic, " ", ig, "/", shell, " ",ig, "/", code, ".R ", 
                vax, " ", doses, " ", rd, " ", indic))

}


# make_qsub, adapted to launch vaccine scripts
make_qsub_postest_each_vax <- function(user = 'USERNAME',
                                       cores = slots,
                                       memory = 100,
                                       proj  = 'PROJECT',
                                       ig = indicator_group,
                                       indic = indicator,  
                                       vax = vaccine,                          
                                       rd = run_date,
                                       log_location = 'sgeoutput',
                                       code,
                                       shell = "r_shell.sh",
                                       geo_nodes = FALSE){

   sharedir <- sprintf("FILEPATH")
   dir.create(sharedir, recursive = T, showWarnings = F)

   if(log_location=='sgeoutput')
     logloc = sprintf("FILEPATH")
   if(log_location=='sharedir'){
     logloc = sprintf('%s/output/%s',sharedir,rd)
     dir.create(sprintf('%s/errors',logloc), recursive = T, showWarnings = F)
     dir.create(sprintf('%s/output',logloc), recursive = T, showWarnings = F)
   }

  ## do we want to submit to geo nodes? if so, there are a few tweaks
  if(geo_nodes == TRUE){
    shell     <- "FILEPATH"   ## for the correct path to R
    proj      <- "PROJECT"    ## correct proj for geos nodes
    node.flag <- "NODES" ## send the job to geo nodes
  }else{
    node.flag <- "" ## don't do anything special
  }

  return(paste0("qsub -e ", logloc, "/errors -o ", logloc, "/output -cwd -l mem_free=", memory, "G -pe multi_slot ", cores, 
                " -P ", proj, " ", node.flag, " -N pe_master_", indic, " ", ig, "/", shell, " ",ig, "/", code, ".R ", 
                indic, " ", vax, " ", rd))

}

waitforvaccinestofinish <- function(sleeptime=100,
                                    indicators,
                                    rd   = run_date,
                                    showfiles = TRUE){

      files =  paste0("FILEPATH")
      n_finished <- sum(file.exists(files)==T)

      lv <- data.table(indicator = indicators,
                       file = files)

      while(n_finished != nrow(lv)){

            n_finished <- sum(file.exists(lv$file) == T)

            message('\n====================================================================================')
            message(sprintf('=====================      Run Date: %s      ======================',rd))
            message(paste0('\nAt ',Sys.time(),',  ',n_finished,' indicators have written output.'))
            if(showfiles){
              message('\nCurrently missing indicators:')
              for(i in 1:nrow(lv))
                if(file.exists(lv[i, file]) == F)
                  message(paste('Indicator =',lv[i,1]))
            }
            message('\n====================================================================================')
            message('====================================================================================')
            message("\n")
            Sys.sleep(sleeptime)
      }

      unlink(lv$file) # Clean up by deleting extra files once done with the loop
}


waitforcombinetofinish <- function(sleeptime=100,
                                   combine_indicators,
                                   rd   = run_date,
                                   regions = Regions,
                                   showfiles = TRUE){

      lv <- data.table(expand.grid(combine_indicators, regions, stringsAsFactors = F))
      names(lv) <- c("combine_indicator", "region")
      setorder(lv, combine_indicator)
      lv[, file := paste0("FILEPATH")]
      n_finished <- sum(file.exists(lv$file)==T)

      while(n_finished != nrow(lv)){

            n_finished <- sum(file.exists(lv$file) == T)

            message('\n====================================================================================')
            message(sprintf('=====================      Run Date: %s      ======================',rd))
            message(paste0('\nAt ',Sys.time(),',  ',n_finished,' indicators have written output.'))
            if(showfiles){
              message('\nCurrently missing indicators:')
              for(i in 1:nrow(lv))
                if(file.exists(lv[i, file]) == F)
                  message(paste0('Indicator = ',lv[i,1]), " | Region = ", lv[i,2])
            }
            message('\n====================================================================================')
            message('====================================================================================')
            message("\n")
            Sys.sleep(sleeptime)
      }

      unlink(lv$file) # Clean up by deleting extra files once done with the loop
}

## check_input_data ################################################

#' Checks to see if input data has been saved into sharedir
#' If not, saves an input data object for the given region
#' This is useful if post-estimating indicators that are calculated
#' from other modeled indicators (e.g. if using continuation-ratio
#" approach for ordinal indicators)
#'
#' @param indicator Indicator
#' @param indicator_group Indicator group
#' @param age Age group
#' @param reg Region
#' @param run_date Run date
#' @param use_share Should the file look in share for master input_data csv?
#'
#' @return Does not return anything; saves input_data files by region
#'         into standard directory for this run
#'
#' @examples
#' 

check_input_data <- function(indicator,
                             indicator_group,
                             age,
                             reg,
                             holdout = 0,
                             run_date,
                             use_share = F,
                             ylist = year_list,
                             shapefile_version=modeling_shapefile_version,
                             sunmational_raking=subnational_raking) {

  message(paste0("Checking input data for indicator: ", indicator,
                 " | region: ", reg, 
                 " | holdout: ", holdout))
  # make a pathaddin
  pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
  input_file <- paste0("FILEPATH")
  
  if (!file.exists(input_file)) {
    ## Load simple polygon template to model over
    gaul_list           <- get_adm0_codes(reg, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T, shapefile_version=modeling_shapefile_version)
    subset_shape        <- simple_polygon_list[[1]]
    simple_polygon      <- simple_polygon_list[[2]]

    message('Loading in full dataset using load_input_data() and saving...')

    df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                          simple      = simple_polygon,
                          agebin      = age,
                          removeyemen = TRUE,
                          pathaddin   = pathaddin,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = ylist)
  }
}

## distribute_config ################################################

#' Distribute config file to directories for a set of indicators with the same run date
#'
#' After running this, use load_config() with post_est_only == T to pull this config
#' Allows the config to be loaded once only from the master script, e.g. in case
#' you want to change it and launch another model shortly afterwards
#'
#' @param cfg config object (from `load_config()`)
#' @param indicators vector of indicator names to distribute config to
#' @param indicator_group indicator group
#' @param run_date shared run date for all indicators
#'
#' @return does not return any objects; writes CSV with config to each of the 
#'         indicator directories on "FILEPATH"
#' @examples
#' 
# Distribute config file to directories for a set of indicators with the same run date.
# After running this, use `load_config()` with `post_est_only == T` to pull this config.
# Allows the config to be loaded once only from the master script, e.g. in case
#   you want to change it and launch another model shortly afterwards

distribute_config <- function(cfg = config,
                              indicators,
                              ig = indicator_group,
                              rd = run_date) {
  for (ind in indicators) {
    ind_dir <- paste0("FILEPATH")
    dir.create(ind_dir, recursive = T, showWarnings = F)
    write.csv(cfg, paste0(ind_dir,'/config.csv'), row.names = FALSE)
  }
}


load_gbd_covariates_newvax = function(vaccines, year_ids, gaul_list, template, use_subnationals = F){
  
  #load the analysis shapefile
  world_shape = readRDS("FILEPATH")
  
  #if we are not using subnationals, keep only national units; otherwise remove the parent units
  if(!use_subnationals){
    world_shape = world_shape[world_shape$level==3,]
  }else{
    world_shape = world_shape[!world_shape$loc_id %in% unique(world_shape$parent_id),]
  }
  
  world_shape = crop(world_shape, template)
  
  afras = rasterize(world_shape, template, 'GAUL_CODE', fun = 'last')
  
  #check to make sure all of the relevant gauls are in the file
  if(!all(world_shape$GAUL_CODE %in% unique(afras))){
    afras = raster::rasterize(world_shape[!world_shape$GAUL_CODE %in% unique(afras),], afras, 'GAUL_CODE', fun = 'first', update = T)
  }
    
  #loop over requested covariates
  fetch_gbd_cov = function(vax, afras){

    gbd <- load_newest_gbd_vax(vaccine = vax, gaul_list = unique(afras), years = year_list)
    
    #for each gaul code and year, update the values
    blank = brick(lapply(year_ids, function(x) afras * NA))
    
    for(yyy in 1:length(year_ids)){
      for(ggg in unique(afras)){
        blank[[yyy]][which(raster::getValues(afras)==ggg)] = gbd[name == ggg & year == year_list[yyy], mean]
      }
    }
    
    names(blank) = rep(vax, times = dim(blank)[3])
    
    return(blank)
  }
  
  all_gbd_layers = lapply(vaccines, function(x) fetch_gbd_cov(x, afras))
  names(all_gbd_layers) = vaccines
  
  return(all_gbd_layers)
  
}

## summarize_admins ################################################

#' Function to summarize admin_pred objects
#'
#' This is a wrapper for `make_admin_pred_summary()`
#'
#' @param ind indicator
#' @param ig indicator_group
#' @param summstats Summary statistics (functions) to compute.
#'                  Order will be the order they are written in csv
#'                  This is passed to `make_admin_pred_summary()`
#' @param raked Raked (T), unraked (F), or both (`c(T,F)`)?
#' @param ad_levels Admin levels to summarize (0, 1, and 2 by default)
#' @return Writes csv files to 'FILEPATH'
#' @examples
#' summarize_admins(summstats = c("mean", "lower", "upper", "cirange"), 
#'                  ad_levels = c(0,1,2), 
#'                  raked     = c(T,F))

summarize_admins_p_below <- function(ind = indicator,
                                     ig = indicator_group,
                                     raked = c(T,F),
                                     ad_levels = c(0,1,2),
                                     rd = run_date,
                                     target = 0.8) {

  sharedir       <- sprintf("FILEPATH")
  input_dir <- paste0("FILEPATH")
  output_dir <- paste0("FILEPATH")
  dir.create(output_dir, recursive = T, showWarnings = F)

  # Convert raked to character
  rr <- character()
  if (T %in% raked) rr <- c(rr, "raked")
  if (F %in% raked) rr <- c(rr, "unraked")

  # Summarize and save admin preds
  for (rake in rr) {
    load(paste0(input_dir, ind, "_", rake, "_admin_draws_eb_bin0_0.RData"))
    for (ad in ad_levels) {
      message(paste0("Summarizing ", ind, ": admin ", ad, " (", rake, ")"))
      ad_summary_table <- make_admin_pred_summary(admin_pred = get(paste0("admin_", ad)),
                                                  sp_hierarchy_list,
                                                  summary_stats = "p_below",
                                                  value = target,
                                                  equal_to = F)
      fwrite(ad_summary_table, 
             file = paste0(output_dir, ind, "_admin_", ad, "_", rake, "_p_below_", target, ".csv"))
    }
  }
}

## save_stratum_ho() ################################################

#' Save a stratum_ho object to a standardized location
#'
#' @param indic indicator for which holdouts have been generated
#' @param ig indicator group for `indic`
#' @param rd run_date for `indic`
#' @param stratum_obj the `stratum_ho` object to save
#' @return Saves an .rds file to the indicator's directory
#' @examples
#' 

save_stratum_ho <- function(indic = indicator,
                            ig = indicator_group,
                            rd = run_date,
                            stratum_obj = stratum_ho) {

  out_dir <- paste0("FILEPATH")
  saveRDS(stratum_obj, file = paste0(out_dir, "stratum.rds"))
}

## recreate_holdouts() ################################################
#' For indicator Y, recreate holdouts generated for indicator X using a 
#' row ID column. This is useful, for instance, in ordinal regression
#' settings if you want to generate holdouts once for the entire data
#' set, and apply those holdouts to each derivative data object (eg
#' the conditional data objects).  
#' 
#' In order for this to work, the holdouts must be generated on the 
#' most comprehensive data set (the one with all of the row_ids).  Those
#' holdouts can then be applied to other data sets with either the same
#' row_ids or a subset of the row_ids from the original data set
#'
#' @param data new data that you'd like to apply the holdouts to
#' @param row_id_col column in both `data` and the original data set
#'                   used to generate the master set of holdouts that
#'                   links the two data sets. Must be unique within
#'                   each data set
#' @param load_from_indic which indicator (within the same run_date and
#'                        indicator group) would you like to load the master
#'                        set of holdouts from?
#' @param rd run_date
#' @param ig indicator group
#' @return a `stratum_ho` style object containing the data from `data` but 
#'         divided up by the holdouts from `load_from_indic`
#' 

recreate_holdouts <- function(data = df,
                              row_id_col = "row_id",
                              load_from_indic,
                              rd = run_date,
                              ig = indicator_group) {
  data <- copy(data)
  n_before <- nrow(data)
  s_ho <- readRDS(paste0("FILEPATH"))

  stratum_df <- rbindlist(s_ho)

  # Peel off just the columns of interest
  stratum_df <- subset(stratum_df, select = c(row_id_col, "t_fold", "fold", "region"))

  # Replace any region columns
  if ("region" %in% names(data)) data[, region := NULL]

  # Merge with df; drop if no holdout assigned (buffer points)
  data <- merge(data, stratum_df, all.x = F, all.y = F, by = row_id_col)
  n_after <- nrow(data)
  
  # get regions
  regions <- unique(data$region)

  out_list <- lapply(regions, function(r) {
    return(data[region == r,])
  })

  names(out_list) <- paste0("region__", regions)

  return(out_list)

}


## check_for_tifs ################################################
#' Check for .tif files for a given set of indicators, summary stats, and (optional) regions
#' within a particular run date.
#'
#' @param indicators vector of indicators
#' @param rd run date
#' @param ig indicator group
#' @param regs vector of regions. If NULL will search for the already-combined files
#' @param rake either T/F (for all to take that value) or vector of T/F indicating whether indicator was raked or not
#' 
#' @return either gives a success message and returns NULL or gives a data table of missing tifs
#' @examples

check_for_tifs <- function(indicators,
                           rd   = run_date,
                           ig   = indicator_group,
                           regs = NULL,
                           rake = F,
                           summstats) {

  library(data.table)
  str_match <- stringr::str_match

  # Set up raked params
  if (length(rake) == 1) {
    rake <- rep(rake, length(indicators)) 
  } else if ((length(rake) > 1) & (length(rake) != length(indicators))) {
    stop("length(rake) must equal length(indicators)")
  }

  pasted_indicators <- paste(indicators, rake, sep = "_-_")

  if (!is.null(regs)) {
    tif_table <- data.table(expand.grid(regs, summstats, pasted_indicators), stringsAsFactors = F)
    names(tif_table) <- c("region", "summstat", "pasted_indicator")
  } else if (is.null(regs)) {
    tif_table <- data.table(expand.grid(summstats, pasted_indicators), stringsAsFactors = F)
    names(tif_table) <- c("summstat", "pasted_indicator")
    
    # Fix cirange for naming convention for merged rasters
    tif_table[summstat == "cirange", summstat := "range"]
  }

  tif_table[, indicator := str_match(pasted_indicator, "(.*)_-_(.*)")[,2]]
  tif_table[, raked := str_match(pasted_indicator, "(.*)_-_(.*)")[,3]]
  tif_table[, pasted_indicator := NULL]
  tif_table[raked == T, raked_addin := "_raked"]
  tif_table[raked == F, raked_addin := ""]

  # Make filenames
  if (!is.null(regs)) {
    tif_table[, filename := paste0("FILEPATH")]
  } else if (is.null(regs)) {
    tif_table[, filename := paste0("FILEPATH")]
  }

  tif_table[,raked_addin := NULL]

  # Check if files exist
  tif_table[, exists := file.exists(filename)]
  missing_files <- subset(tif_table, exists == F)

  # Return
  if (nrow(missing_files) == 0) {
    message("All files are present.")
    return(NULL)
  } else if (nrow(missing_files) > 0) {
    message(paste0(nrow(missing_files), " out of ", nrow(tif_table), " files are missing."))
    message("Returning data table of missing files.")
    return(missing_files)
  }
  
}

# Save TIFs to the standard directory / format for mapping
copy_tif_to_map_input_dir <- function(ind,
                                      ig = indicator_group, 
                                      measure, 
                                      rd = run_date, 
                                      raked, 
                                      yl = year_list) {

  # Fix measures
  if (measure == "cirange") measure <- "range"

  # Fix raked if boolean
  if (raked == T) raked <- "raked"
  if (raked == F) raked <- "unraked"

  # Define input dirs / files
  in_dir <- paste0("FILEPATH")
  in_file <- paste0("FILEPATH")
  
  # Define output dirs / files
  out_dir <- paste0("FILEPATH")
  dir.create(out_dir, recursive = T, showWarnings =F)

  out_filename <- paste0(ind, "_", measure, "_", raked, "_", min(yl), "_", max(yl), ".tif")
  out_file <- paste0(out_dir, out_filename)

  success <- file.copy(in_file, out_file, overwrite = T)

  if (success) message(paste0(out_filename, " copied successfully!"))
  if (!success) warning(paste0("Failed to copy ", out_filename, "!")) 

  return(invisible(success))

}

# Save CSVs to the standard directory / format for mapping
copy_admins_to_map_input_dir <- function(ind,
                                         ig = indicator_group, 
                                         measure, 
                                         rd = run_date, 
                                         raked, 
                                         yl = year_list, 
                                         ad_level) {

  str_match <- stringr::str_match

  # Fix raked if boolean
  if (raked == T) raked <- "raked"
  if (raked == F) raked <- "unraked"

  # Define basic dirs / files
  in_dir <- paste0("FILEPATH")
  out_dir <- paste0("FILEPATH")
  dir.create(out_dir, recursive = T, showWarnings =F)

  if (measure %in% c("mean", "lower", "upper", "cirange", "cfb")) {
    in_dir <- paste0("FILEPATH")
    in_file <- paste0(in_dir, ind, "_admin_", ad_level, "_", raked, "_summary.csv")  
    in_df <- as.data.table(read.csv(in_file, stringsAsFactors = F))
    setnames(in_df, measure, "value")

    if (measure == "cirange") measure <- "range" #naming convention
    out_df <- subset(in_df, select = c(paste0("ADM", ad_level, "_CODE"), "year", "value")) 
    out_filename <- paste0(ind, "_", measure, "_", raked, "_ad", ad_level, ".csv")
    message(paste0("Writing ", out_filename, "..."))
    write.csv(out_df, file = paste0(out_dir, out_filename))
  
  } else if  (grepl("psup", measure)) {

    pct <- str_match(measure, "psup(.*)")[,2]
    pct_div100 <- as.numeric(pct)/100

    in_dir <- paste0("FILEPATH")
    in_file <- paste0(in_dir, ind, "_admin_", ad_level, "_", raked, "_p_", pct_div100, "_or_better_summary.csv")  
    in_df <- as.data.table(read.csv(in_file, stringsAsFactors = F))
    setnames(in_df, "p_above", "value")

    out_df <- subset(in_df, select = c(paste0("ADM", ad_level, "_CODE"), "year", "value"))
    out_filename <- paste0(ind, "_psup", pct, "_", raked, "_ad", ad_level, ".csv")

    message(paste0("Writing ", out_filename, "..."))
    write.csv(out_df, file = paste0(out_dir, out_filename))

  } else if (grepl("diff", measure)) {
  
    in_dir <- paste0("FILEPATH")
    in_file <- paste0(in_dir, ind, "_admin_", ad_level, "_", raked, "_", measure, ".csv")

    in_df <- as.data.table(read.csv(in_file, stringsAsFactors = F))
    setnames(in_df, "mean", "value")

    out_df <- subset(in_df, select = c(paste0("ADM", ad_level, "_CODE"), "year", "value"))
    out_filename <- paste0(ind, "_diff_", min(year_list), "-", max(year_list), "_", raked, "_ad", ad_level, ".csv")

    message(paste0("Writing ", out_filename, "..."))
    write.csv(out_df, file = paste0(out_dir, out_filename))

  } else {

    warning(paste0("Measure ", measure, " not found for indicator ", ind, " | ", raked, " | admin level ", ad_level))

  }
}


# Copy tifs to geoviz directory
copy_tif_to_geoviz_dir <- function(ind,
                                   ig = indicator_group,
                                   measure,
                                   rd = run_date,
                                   raked,
                                   yl = year_list) {

  message(paste0(ind, " | ", measure))

  # Fix measures
  if (measure == "cirange") measure <- "range"

  # Fix raked if boolean
  if (raked == T) raked <- "raked"
  if (raked == F) raked <- "unraked"

  # Define input dirs / files
  in_dir <- paste0("FILEPATH")
  in_file <- paste0("FILEPATH")

  # Define output dirs / files
  out_dir <- paste0("FILEPATH")
  dir.create(out_dir, recursive = T, showWarnings =F)

  # Load raster and save each layer individually
  ras <- brick(in_file)

  # Adjust so that values are scaled to 0-100 for percents
  ras <- ras*100

  for (i in 1:length(year_list)) {
    yy <- year_list[i]
    ind_name <- ind
    out_filename <- paste0(ind_name, "_", measure, "_", yy)
    out_file <- paste0(out_dir, out_filename)

    writeRaster(ras[[i]],
                file = out_file,
                overwrite = TRUE,
                format='GTiff')
  }
}

# Copy admin csvs to geoviz dir
copy_admins_to_geoviz_dir <- function(ind,
                                      ig = indicator_group,
                                      measure,
                                      rd = run_date,
                                      raked,
                                      yl = year_list) {

  # Fix raked if boolean
  if (raked == T) raked <- "raked"
  if (raked == F) raked <- "unraked"

  # Define basic dirs / files
  in_dir <- paste0("FILEPATH")
  out_dir <- paste0("FILEPATH")
  dir.create(out_dir, recursive = T, showWarnings =F)

  all_admins <- lapply(c(0,1,2), function(ad_level) {
                  in_dir <- paste0("FILEPATH")
                  in_file <- paste0(in_dir, ind, "_admin_", ad_level, "_", raked, "_summary.csv")
                  in_df <- fread(in_file)
                  setnames(in_df, measure, "value")
                  setnames(in_df, paste0("ADM", ad_level, "_CODE"), "gaul_code")
                  return(subset(in_df, select = c("gaul_code", "value", "year")))
                })
  all_admins <- rbindlist(all_admins)
  all_admins <- all_admins[!is.na(value),]

  for (yy in yl) {
    year_admins <- subset(all_admins, year == yy)
    year_admins[, year := NULL]
    year_admins[, value := value * 100] # Adjust so that values scaled 0-100 for percents
    ind_name <- ifelse(ind=="dpt1_3_abs_dropout", "dropout", ind)
    fwrite(year_admins, file = paste0(out_dir, ind_name, "_", measure, "_", yy, ".csv"))
  }

}

  
#' Plot MBG vs GBD estimates along with data
#' @param ind indicator
#' @param ig indicator_group
#' @param rd run_date
#' @param gbd_ind gbd indicator name (without 'vacc_' prefix; i.e. "dpt3")
#' @param new_vax use new vaccines from load_newest_gbd_vax? 
#'          (only T for now; need to build in basic GBD defaults)   
#' @param yl year_list
#' @param vax_title title of vaccine for display
#' @param raked use raked estimates? (T/F)
#' @param include_biascorr include bias-corrected admin data 
#'                         (from MBG standard covariate folder)? 
#' @param regions NULL by default (use all regions); pass vector of regions
#'                to limit analysis to those regions only
#' @param pathaddin addin for filename to distinguish between subsets if needed
#' 
#' @return creates png images in a /comparisons/ folder in output directory
#' 
plot_mbg_vs_gbd <- function(ind = indicator,
                            ig = indicator_group,
                            rd = run_date,
                            gbd_ind,
                            new_vax = T,
                            yl = year_list,
                            vax_title,
                            raked = F,
                            include_biascorr = T,
                            gbd_date = "best",
                            regions = NULL,
                            pathaddin = "",
                            shapefile_version=modeling_shapefile_version,
                            subnational_raking=subnational_raking) {
  # Compare admin0s

  # Set up files and directories
  rd_dir <-paste0("FILEPATH")
  adm_dir <- paste0("FILEPATH")
  adm0_file <- paste0(adm_dir, ind, "_admin_0_", ifelse(raked, "raked", "unraked"), "_summary.csv")
  biascorr_file <- paste0("FILEPATH")

  # Logical check for bias correction
  if (include_biascorr == T) {
    if (!file.exists(biascorr_file)) {
      warning(paste0("No bias correction file found at \n", 
                     "file path: ", biascorr_file, "\n",
                     "Setting include_biascorr to FALSE and continuing..."))
    }
  }

  out_dir <- paste0(rd_dir, "comparisons/")
  dir.create(out_dir, showWarnings = F)

  adm0_df <- fread(adm0_file)

  gbd <- load_newest_gbd_vax(vaccine = gbd_ind,
                             gaul_list = NULL,
                             years = yl,
                             return_cis = T,
                             gbd_date = gbd_date,
                             return_mode = "national",
                             return_field = "location_id")

  loc_maps <- get_location_code_mapping(shapefile_version)
  loc_maps <- subset(loc_maps, select = c("loc_id", "ADM_CODE"))
  setnames(loc_maps, c("loc_id", "ADM_CODE"), c("name", "ADM0_CODE"))
  gbd <- merge(gbd, loc_maps)

  setnames(gbd, c("name", "value", "lower", "upper"), c("loc_id", "gbd", "gbd_lower", "gbd_upper"))
  adm0_df <- merge(adm0_df, gbd, by = c("ADM0_CODE", "year"), all.x = T)

  missing_countries <- unique(adm0_df[is.na(gbd), ADM0_NAME])
  if (length(missing_countries) > 0) {
    warning(paste0("The following countries have no GBD value and will be dropped: \n",
            paste(missing_countries, collapse = ", ")))
  }

  adm0_df <- subset(adm0_df, !(ADM0_NAME %in% missing_countries))

  # Truncate long country names
  adm0_df[ADM0_NAME == "Democratic Republic of the Congo", ADM0_NAME := "DRC"]
  adm0_df[ADM0_NAME == "Sao Tome and Principe", ADM0_NAME := "STP"]
  adm0_df[ADM0_NAME == "United Republic of Tanzania", ADM0_NAME := "Tanzania"]
  adm0_df[ADM0_NAME == "Central African Republic", ADM0_NAME := "CAR"]
  adm0_df[ADM0_NAME == "Equatorial Guinea", ADM0_NAME := "Eq. Guinea"]

  # Now pull in the input data
  output_draws <- fread(paste0(rd_dir, "output_draws_data.csv"))

  # Drop draw columns
  output_draws <- subset(output_draws, select = names(output_draws)[!grepl("draw[0-9]+", names(output_draws))])
  output_draws[, V1 := NULL]
  
  output_draws[, outcome := get(ind) / N]
  setnames(output_draws, "ad0", "ADM0_CODE")
  
  group_by <- c("svy_id", "source", "country", "ADM0_CODE", "year")
  data_est <- output_draws[, .(mean = weighted.mean(outcome, w = weighted_n),
                               data_tot_n = sum(weighted_n)), 
                 by = group_by]
  data_est[, estimate := "Data"]

  # Load in bias corrected administrative data if needed
  if (include_biascorr) {
    df_biascorr <- readRDS(biascorr_file)
    loc_mapping <- get_location_code_mapping(shapefile_version)
    df_biascorr <- merge(subset(df_biascorr, select = c("loc_id", "data", "year_id", "interpolated")),
                         subset(loc_mapping, select = c("loc_id", "ihme_lc_id", "GAUL_CODE")), 
                         by = "loc_id", all.x = T, all.y = F)
    setnames(df_biascorr, c("GAUL_CODE", "data", "ihme_lc_id", "year_id"), c("ADM0_CODE", "mean", "country", "year"))
    df_biascorr <- subset(df_biascorr, year %in% yl)
  }

  # Subset data if needed
  if (!is.null(regions)) {
    gaul_list <- get_adm0_codes(regions, shapefile_version = shapefile_version, subnational_raking=subnational_raking)
    data_est <- subset(data_est, ADM0_CODE %in% gaul_list)
    adm0_df <- subset(adm0_df, ADM0_CODE %in% gaul_list)
    if (include_biascorr) df_biascorr <- subset(adm0_df, ADM0_CODE %in% gaul_list)
  }

  # Make plots
  # Max = 16 per page
  n_pages <- ceiling(length(unique(adm0_df$ADM0_NAME))/16)
  adm0_df <- adm0_df[order(ADM0_NAME, year)]
  assign_table <- data.table(ADM0_NAME = unique(adm0_df$ADM0_NAME),
                             page = rep(1:n_pages, each=16, length.out = length(unique(adm0_df$ADM0_NAME))))
  adm0_df <- merge(adm0_df, assign_table, by = "ADM0_NAME")                           

  for (pg in 1:n_pages) {

    # Prepare time series data
    plot_df <- adm0_df[page == pg]
    plot_df_gbd <- subset(plot_df, select = !(names(plot_df) %in% c("mean", "lower", "upper")))
    plot_df_gbd[, estimate := "GBD Estimate"]
    setnames(plot_df_gbd, c("gbd", "gbd_lower", "gbd_upper"), c("mean", "lower", "upper"))
    plot_df_mbg <- subset(plot_df, select = !(names(plot_df) %in% c("gbd", "gbd_lower", "gbd_upper")))
    plot_df_mbg[, estimate := paste0(ifelse(raked, "Raked", "Unraked"), " MBG Estimate")]
    plot_df <- rbind(plot_df_mbg, plot_df_gbd)
    
    # Add data estimates
    data_df <- subset(data_est, ADM0_CODE %in% unique(plot_df$ADM0_CODE))
    setnames(data_df, "mean", "data_mean")
    data_df <- merge(data_df, unique(subset(plot_df, select = c("ADM0_CODE", "ADM0_NAME"))))

    all_df <- rbind(data_df, plot_df, fill=T, use.names=T)

    if (include_biascorr) {
      df_biascorr_subset <- subset(df_biascorr, ADM0_CODE %in% unique(plot_df$ADM0_CODE))
      setnames(df_biascorr_subset, "mean", "biascorr_mean")
      df_biascorr_subset[, estimate := "Bias-corr admin data"]
      df_biascorr_subset <- merge(df_biascorr_subset, 
                                  unique(subset(plot_df, select = c("ADM0_CODE", "ADM0_NAME"))),
                                  by = "ADM0_CODE")
      all_df <- rbind(all_df, df_biascorr_subset, fill = T, use.names = T)
    }

    out_file = paste0(out_dir, ind, "_compare_", ifelse(raked, "raked", "unraked"),
                      "_mbg_to_gbd_ad0_pg", pg, pathaddin, ".png")

    png(filename=out_file, 
        type="cairo",
        units="in", 
        width=14, 
        height=8, 
        pointsize=12, 
        res=300)

    gg <- ggplot(data = all_df, 
                 aes(x = year, color = estimate, group = estimate)) +
          geom_point(aes(y = data_mean, size = data_tot_n), alpha = 0.5, shape = 1)

    if (include_biascorr) {
      gg <- gg + 
            geom_point(aes(y = biascorr_mean), alpha = 0.5, shape = 2)
    }

    gg <- gg + 
          geom_errorbar(aes(ymin = lower, ymax = upper),
                        width = 0, alpha = 0.2) +
          geom_line(aes(y = mean), alpha = 0.8) +
          facet_wrap(~ADM0_NAME) +
          theme_bw() +
          labs(x = "Year", 
               y = paste0(vax_title, " Coverage"),
               color = "Estimate",
               shape = "Estimate",
               size = "Data: Weighted N") +
          scale_color_manual(values = c("#4daf4a", "#984ea3", "#e41a1c", "#377eb8")) +
          ylim(c(0,1))

    plot(gg)

    dev.off()
  }
}

#' Plot summary metrics by year and by region
#' Function to plot summary metrics by year and by region, depending on
#' what options are selected.  Produces time trends.  Useful as diagnostic
#' plots, for instance, or for publications
#'
#' @param indic indicator
#' @param ig indicator group
#' @param rd run date
#' @param use_oos use OOS predictions? T, F, or c(T,F)
#' @param by_regions plot by regions? Requires that summary
#'                   metric tables have been generated by region
#' @param ad_levels country, ad1, or ad2 - or a vector of these
#' @return writes a series of png files to [run date dir on share]/summary_metrics/
#' 
plot_summary_metrics <- function(indic,
                                 ig,
                                 rd,
                                 use_oos = c(T,F),
                                 by_regions = c(T,F),
                                 ad_levels = c("country", "ad1", "ad2"),
                                 facet_wrap = F) {
  
  
  # Define a color scheme for the regions
  get_color_scheme <- function(theme){
    
    # Set up categorical colors
    carto_discrete <- c("#7F3C8D","#11A579","#F2B701","#E73F74",
                        "#3969AC","#80BA5A","#E68310","#008695",
                        "#CF1C90","#f97b72","#4b4b8f","#A5AA99",
                        "#66C5CC","#F6CF71","#F89C74","#DCB0F2",
                        "#87C55F","#9EB9F3","#FE88B1","#C9DB74",
                        "#8BE0A4","#B497E7","#D3B484","#B3B3B3")
    return(get(theme))
  }
  
  in_dir <- paste0("FILEPATH")
  
  for (level in ad_levels) {
    
    message(paste0("Working on the ", level, " level..."))
    # Read in csvs
    df_reg <- fread(paste0(in_dir, level, "_metrics_by_year_region.csv"))
    # df_reg$Region <- toupper(df_reg$Region)
    df_all <- fread(paste0(in_dir, level, "_metrics_by_year.csv"))
    
    # Combine data tables
    df_all[, Region := "All"]
    df_all <- rbind(df_reg, df_all, use.names = T)
    rename_table <- data.table(rbind(c("Year", "year", "Year"),
                                     c("OOS", "oos", "OOS"),
                                     c("Mean Err.", "me", "Mean Error"),
                                     c("Mean Abs. Err.", "mae", "Mean Absolute Error"),
                                     c("Mean Pred.", "mean_pred", "Mean Prediction"),
                                     c("Mean Obs.", "mean_obs", "Mean Observed value"),
                                     c("RMSE", "rmse", "RMSE"),
                                     c("Median SS", "median_ss", "Median Sample Size"),
                                     c("Corr.", "corr", "Correlation"), 
                                     c("95% Cov.", "cov_95", "95% Coverage"),
                                     c("Region", "region", "Region")))
    
    names(rename_table) <- c("old", "new", "title")
    setnames(df_all, rename_table$old, rename_table$new)
    
    message(paste0("Writing images to ", in_dir))
    # Loop over variables & make plots
    for (is_oos in use_oos) {
      for (by_region in by_regions) {
        for (plot_obj in c("me", "mae", "rmse", "corr", "cov_95")) {
          
          if (level == "ad2") ad_caption = "admin 2"
          if (level == "ad1") ad_caption = "admin 1"
          if (level == "country") ad_caption = "country"
          
          plot_title <- rename_table[new == plot_obj, title]
          
          gg <- ggplot()
          
          if (by_region == T) {
            gg <- gg + geom_line(data = df_all[region != "All" & oos == is_oos,],
                                 aes(x = year, y = get(plot_obj), group = region, color = region), alpha = 0.35) +
              scale_color_manual(values = c("#000000", get_color_scheme("carto_discrete"))) +
              labs(color = "Region") +
              geom_line(data = df_all[region == "All" & oos == is_oos],
                        aes(x = year, y = get(plot_obj), group = region, color = region))
          } else {
            gg <- gg + geom_line(data = df_all[region == "All" & oos == is_oos],
                                 aes(x = year, y = get(plot_obj)), color = "black")
          }
          
          gg <- gg +
            labs(x = "Year", y = plot_title, title = plot_title, 
                 subtitle = paste0(ifelse(is_oos, "Out-of-sample", "In-sample"), " | Aggregated to the ", ad_caption, " level")) +
            theme_bw()
          
          # Get range
          range_obs <- df_all[oos == is_oos, get(plot_obj)]
          
          if (plot_obj == "me") gg <- gg + ylim(-max(abs(range_obs)), max(abs(range_obs))) + geom_abline(slope = 0, intercept = 0, color = "red", linetype = 2)
          if (plot_obj == "mae" | plot_obj == "rmse") gg <- gg + ylim(0, max(abs(range_obs))) + geom_abline(slope = 0, intercept = 0, color = "red", linetype = 2)
          if (plot_obj == "cov_95") gg <- gg + ylim(0,1) + geom_abline(slope = 0, intercept = 0.95, color = "red", linetype = 2)
          if (plot_obj == "corr") gg <- gg + ylim(0,1) + geom_abline(slope = 0, intercept = 1, color = "red", linetype = 2)
          
          message(paste0(" Writing file for ", plot_title, " | by_region = ", by_region, " | oos = ", is_oos))
          png_filename <- paste0(indic, "_summary_metric_plot_", level, "_", 
                                 ifelse(is_oos, "OOS" ,"IS"), "_", ifelse(by_region, "by_region", "all"), 
                                 "_", plot_obj, ".png")
          
          png(file = paste0(in_dir, png_filename), 
              width = 8, 
              height = 4, 
              units = "in",
              res = 300)
          plot(gg)
          dev.off()
          
          message(paste0("  File written to ", png_filename))
        }
      }
    }
  }
}


# Fix raster naming conventions -- rarely needed function
fix_raster_naming_convention <- function(run_date = run_date,
                                         indicator = indicator,
                                         indicator_group = indicator_group) {
  str_replace <- stringr::str_replace
  main_dir <- paste0("FILEPATH")
  files_to_fix <- list.files(path = main_dir,
                             pattern="_unraked_raster.tif", 
                             full.names = T)
  new_file_names <- str_replace(files_to_fix, "_unraked_", "_")
  file.rename(from = files_to_fix, to = new_file_names)

}


get_individual_countries <- function(gaul_list) {
  gaul_ref <- list()

  gaul_ref[["mwi"]] <- 152
  gaul_ref[["nga"]] <- 182
  gaul_ref[["egy"]] <- 40765
  gaul_ref[["gha"]] <- 94
  gaul_ref[["zaf"]] <- 227
  gaul_ref[["cod"]] <- 68
  gaul_ref[["ago"]] <- 8
  gaul_ref[["caf"]] <- 49
  gaul_ref[["cog"]] <- 59
  gaul_ref[["cod"]] <- 68
  gaul_ref[["gnq"]] <- 76
  gaul_ref[["gab"]] <- 89
  gaul_ref[["bwa"]] <- 35
  gaul_ref[["lso"]] <- 142
  gaul_ref[["nam"]] <- 172
  gaul_ref[["zaf"]] <- 227
  gaul_ref[["swz"]] <- 235
  gaul_ref[["zwe"]] <- 271
  gaul_ref[["ben"]] <- 29
  gaul_ref[["bfa"]] <- 42
  gaul_ref[["cmr"]] <- 45
  gaul_ref[["cpv"]] <- 47
  gaul_ref[["tcd"]] <- 50
  gaul_ref[["civ"]] <- 66
  gaul_ref[["gmb"]] <- 90
  gaul_ref[["gha"]] <- 94
  gaul_ref[["gin"]] <- 106
  gaul_ref[["gnb"]] <- 105
  gaul_ref[["lbr"]] <- 144
  gaul_ref[["mli"]] <- 155
  gaul_ref[["mrt"]] <- 159
  gaul_ref[["ner"]] <- 181
  gaul_ref[["nga"]] <- 182
  gaul_ref[["stp"]] <- 214
  gaul_ref[["sen"]] <- 217
  gaul_ref[["sle"]] <- 221
  gaul_ref[["tgo"]] <- 243
  gaul_ref[["bdi"]] <- 43
  gaul_ref[["com"]] <- 58
  gaul_ref[["dji"]] <- 70
  gaul_ref[["eri"]] <- 77
  gaul_ref[["eth"]] <- 79
  gaul_ref[["ken"]] <- 133
  gaul_ref[["mdg"]] <- 150
  gaul_ref[["mwi"]] <- 152
  gaul_ref[["moz"]] <- 170
  gaul_ref[["rwa"]] <- 205
  gaul_ref[["som"]] <- 226
  gaul_ref[["ssd"]] <- 74
  gaul_ref[["tza"]] <- 257
  gaul_ref[["uga"]] <- 253
  gaul_ref[["zmb"]] <- 270
  gaul_ref[["dza"]] <- 4
  gaul_ref[["egy"]] <- 40765
  gaul_ref[["lby"]] <- 145
  gaul_ref[["mar"]] <- 169
  gaul_ref[["sdn"]] <- 6
  gaul_ref[["tun"]] <- 248

  t <- c(unlist(gaul_ref)) %>% as.numeric
  t <- cbind(t, names(gaul_ref)) %>% as.data.table
  names(t) <- c("gaul", "iso3")

  return(t[gaul %in% gaul_list, iso3])
}



# Function to set up indicators for multiple-dose runs
set_up_indicators <- function(stem,
                              doses,
                              single_dose = T,
                              ir = indic_repo,
                              save_doses = F,
                              save_2_cov = F) {
  # Define objects
  vaccine            <- stem

  # Indicators to model in INLA
  model_indicators <- paste0(stem, doses, "_cov")
  if (single_dose == F & vaccine != "dpt") model_indicators <- c(paste0(stem, 1:(doses-1), "_cond"), model_indicators)
  if (single_dose == F & vaccine == "dpt") model_indicators <- c(paste0(stem, "12_cond"), model_indicators)

  # Dose specific indicators (i.e. p(exactly that # of doses))
  dose_indicators <- paste0(stem, doses, "_cov")
  if (single_dose == F & save_doses == T) dose_indicators <- c(paste0(stem, 0:(doses-1), "_dose"), dose_indicators)                          

  # Indicators to rake
  # This will need to be changed once we have full dose-specific information from GBD
  rake_indicators    <- c(unique(paste0(stem, doses, "_cov")))
  if (single_dose == F & vaccine == "dpt") rake_indicators <- c(rake_indicators, 
                                             paste0(stem, "12_cond") ) #,
                                             #paste0(stem, "1_", doses, "_", c("abs", "rel"), "_dropout"))
 

  postest_indicators <- rake_indicators  
  
  all_indicators     <- unique(c(model_indicators, rake_indicators, dose_indicators, postest_indicators))

  # Message outputs
  message("\nGenerating Indicator Lists:")
  message(paste0(" model_indicators:   ", paste(sort(model_indicators), collapse = ", ")))
  message(paste0(" dose_indicators:    ", paste(sort(dose_indicators), collapse = ", ")))
  message(paste0(" rake_indicators:    ", paste(sort(rake_indicators), collapse = ", ")))
  message(paste0(" postest_indicators: ", paste(sort(postest_indicators), collapse = ", ")))
  message(paste0(" all_indicators:     ", paste(sort(all_indicators), collapse = ", ")))

  # Check to ensure all files present
  message("\nChecking to ensure input data files present")

  config <- fread(paste0(ir, "config_", vaccine, ".csv"), header = F)
 
  use_share <- FALSE

  if(use_share==FALSE) load_dir <- paste0("FILEPATH")
  if(use_share==TRUE) load_dir  <- "FILEPATH"

  # For padding in output
  char_space <- max(nchar(all_indicators)) + 2

  # Loop over indicators for check
  for (ind in sort(unique(c(model_indicators)))) {
    check_exists <- file.exists(paste0(load_dir, ind, ".csv"))
    message(paste0("  ", ind, ":", 
                   paste0(rep(" ", char_space - nchar(ind)), collapse = ""), 
                   ifelse(check_exists, "SUCCESS", "FAILURE")))
    if (!check_exists) stop(paste0("CSV file for ", ind, " does not exist!"))
  }

  # Assign everything to global environment
  message("\nAssigning to global environment...")
  for (param in c("vaccine", "doses", "save_2_cov", "save_doses", 
                  "model_indicators", "dose_indicators", 
                  "rake_indicators", "postest_indicators", 
                  "all_indicators")) {
     assign(param, get(param), envir=globalenv())
  }

  message("  Completed successfully!")
}

# Function to load profiling information for a given vaccine, region, and draw combo

load_profiling <- function(vax,
                           draws,
                           reg,
                           ir = indic_repo,
                           default_mem = 200,
                           default_slots = 6) {

  profile_csv <- paste0(indic_repo, "profiling_", vax, ".csv")
  if (file.exists(profile_csv)) {
    profile <- fread(profile_csv)
    mem <- profile[region == reg & samples == draws, mem]
    slots <- profile[region == reg & samples == draws, slots]
  } else {
    mem <- NULL
    slots <- NULL
  }

  if (length(mem) == 0) {
    mem <- default_mem
    message(paste0("No profiling csv entry found for memory -- using default value of ", default_mem))
  }

  if (length(slots) == 0) {
    slots <- default_slots
    message(paste0("No profiling csv entry found for slots -- using default value of ", default_slots))
  }

  return(c(mem = mem, slots = slots))

}