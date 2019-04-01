###############################################################################
###############################################################################
## Create plots of mean error over space and time
##
## Purpose: Generate maps of residual mean error in space and time
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

indicator_group <- "vaccine"
run_date <- NULL # Change to model run date
regions <- c("wssa", "cssa", "essa", "sssa", "name")
ig <- indicator_group

# Load background shapes for map
background_shape <- readRDS("<<<< FILEPATH REDACTED >>>>/background_map_africa.rds")

background_shape@data$id <- rownames(background_shape@data)
background_shape_points <- fortify(background_shape, region = "id")
background_shape_df <- join(background_shape_points, background_shape@data, by = "id")
background_shape_df <- as.data.table(background_shape_df)

theme_empty <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_blank())

color_scheme <- c("#990021", "#CC51AB", "#CE89E5", "#6151CC", "#004999")

for (indicator in c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout")) {

  if (indicator == "dpt3_cov") indicator_title <- "DPT3 coverage"
  if (indicator == "dpt1_cov") indicator_title <- "DPT1 coverage"
  if (indicator == "dpt1_3_abs_dropout") indicator_title <- "DPT 1-3 absolute dropout"
  
  out_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", indicator, "/output/", run_date, 
                  "/number_plugging/residual_error_maps/")
  dir.create(out_dir, showWarnings = F)

  message(paste0("Working on ", indicator, "..."))
  
  sharedir <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", indicator, "/output/", run_date, "/")
  
  # Load draws
  message("   loading data...")
  input_df <- fread(paste0(sharedir, "output_draws_data.csv"))
  
  draw_cols <- names(input_df)[grepl("draw", names(input_df))]
  
  input_df[, pred_mean := rowMeans(.SD, na.rm = T), .SDcols = draw_cols]
  input_df[, me := pred_mean - (get(indicator) / N)]
  
  # Drop draws for memory
  input_df <- subset(input_df, select = !(names(input_df) %in% draw_cols))
  
  # Set years
  input_df[year > 1999 & year <= 2003, plot_year := "2000-2003"]
  input_df[year > 2003 & year <= 2007, plot_year := "2004-2007"]
  input_df[year > 2007 & year <= 2011, plot_year := "2008-2011"]
  input_df[year > 2011 & year <= 2016, plot_year := "2012-2016"]
  
  # Basic plot
  gg <- ggplot() +
    geom_polygon(data = background_shape_df,
                 aes(x = long, y = lat, group = group),
                 fill="white") +
    geom_path(data = background_shape_df,
              aes(x = long, y = lat, group = group),
              color="black", size = 0.2) +
    geom_point(data = input_df[weight != 1,], 
               aes(x = longitude, y = latitude, size = N, alpha = weight, color = me),
               pch = 16) +
    geom_point(data = input_df[weight == 1,], 
               aes(x = longitude, y = latitude, size = N, alpha = weight, color = me),
               pch = 16) +
    facet_wrap(~plot_year, ncol = 2) +
    scale_size_area(max_size = 1) +
    scale_color_gradientn(colors = color_scheme, limits = c(-1,1)) +
    coord_equal() +
    labs(title = paste0("Residual error: ", indicator_title), 
         color = "Residual error", alpha = "Weight", size = "N") +
    theme_empty
  
  message("    plotting...")
  png(file = paste0(out_dir, "error_plot_", indicator, ".png"),
      height = 7,
      width = 8,
      units = "in",
      res = 300)
  print(gg)
  dev.off()
  
}