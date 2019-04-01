###############################################################################
###############################################################################
## Covariate range plots
##
## Purpose: Produce plots of covariate distributions for locations with and
##          without data observations
###############################################################################
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<<<< FILEPATH REDACTED >>>>'
indic_repo         <- '<<<< FILEPATH REDACTED >>>>'

## sort some directory stuff
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific code begins here ##########################################

# Set up run options #########################################################

rd <- NULL # Change this to reflect the model's run date
ind <- "dpt3_cov"
ig <- "vaccine"
stacker_names <- c("gam", "gbm", "lasso")
year_list <- c(2000:2016)

out_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", rd, 
                  "/number_plugging/cov_range_plots/")
dir.create(out_dir, showWarnings = F)

# Loop over regions ##########################################################

for (reg in c("cssa", "name", "sssa", "essa", "wssa")) {

  message(paste0("Working on region ", reg, "..."))
  
  # Load stackers
  load(paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/model_image_history/",
              rd, "_bin0_", reg, "_0.RData"))  
  
  stacker_list <- cov_list[which(names(cov_list) %in% stacker_names)]
  
  # Load data
  df <- fread(paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", rd, 
                     "/input_data_bin0_", reg, "_0.csv"))
  
  extract_values <- function(stacker, yr, stack_list = stacker_list, yl = year_list, the_df = df) {
    df_year <- subset(the_df, year == yr)
    if (nrow(df_year) == 0) return(NULL)
    ras_year <- stack_list[[stacker]][[which(year_list == yr)]]
      
    # Determine locations where data is present
    ras_idx <- matrix(1:(nrow(ras_year)*ncol(ras_year)), nrow=nrow(ras_year), ncol=ncol(ras_year))
    ras_idx <- raster(x = ras_idx)
    crs(ras_idx) <- crs(ras_year)
    extent(ras_idx) <- extent(ras_year)  
    
    vals_idx <- na.omit(raster::extract(ras_idx, cbind(df_year$longitude, df_year$latitude)))
    cells_with_data <- as.vector(as.matrix(ras_year))[vals_idx]
    if (length(na.omit(cells_with_data)) == 0) return(NULL)
    cells_without_data <- as.vector(as.matrix(ras_year))[-na.omit(vals_idx)]
    
    range_with_data <- range(cells_with_data, na.rm = T)
    
    df_with_data <- data.table(val = na.omit(cells_with_data), 
                               data = "with_data",
                               stacker = stacker,
                               year = yr)
    
    df_without_data <- data.table(val = na.omit(cells_without_data),
                                  data = "without_data",
                                  stacker = stacker,
                                  year = yr)
    
    df_return <- rbind(df_with_data, df_without_data, use.names=T)
   
  }
  
  out_list <- list()
  
  for (ss in stacker_names) {
    for (yy in year_list) {
      message(ss, "_", yy)
      out_list[[paste0(ss, "_", yy)]] <- extract_values(stacker = ss, yr = yy)
    }
  }
  
  out_df <- rbindlist(out_list)
  out_df[stacker == "gam", stacker := "GAM"]
  out_df[stacker == "gbm", stacker := "GBM"]
  out_df[stacker == "lasso", stacker := "Lasso"]
  out_df[data == "with_data", data := "Observed coverage data"]
  out_df[data == "without_data", data := "No observed coverage data"]
  
  gg_all <- ggplot(data = out_df, aes(x = val, color = data)) +
    geom_density() +
    theme_bw() +
    facet_wrap(~stacker) +
    labs(x = "Value", y = "Density", color = "Pixel type") +
    theme(legend.position = "bottom")
  
  png(file = paste0(out_dir, "cov_range_plot_all_", reg, ".png"),
      height = 3, width = 8,
      units = "in", 
      res = 300)
  
  print(gg_all)
  dev.off()

}