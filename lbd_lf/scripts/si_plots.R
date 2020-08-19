### Making some plots for LF prevalence paper SI
######################################################################################

## Get current user name
user <- Sys.info()[["user"]]

## Set repo location and indicator group
core_repo <- paste0(<<<< FILEPATH REDACTED >>>>)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load libraries, packages, and miscellaneous MBG project functions.

commondir <- sprintf(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))
package_list <- c(package_list, "sf")

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

# These may not actually be necessary for this particular script.
# if these packages are not already installed in user's home directory. Consider trying to run the script first, before installing.
path <- paste0(<<<< FILEPATH REDACTED >>>>)
library(fasterize, lib.loc = path)
library(sf, lib.loc = path)
library(matrixStats, lib.loc = path)

if (TRUE) { ## Pull in custom LF scripts, some of which mask existing central functions.
  setwd(indic_repo)
  for (funk in list.files(indic_repo, recursive = TRUE, pattern = "functions")) {
    if (length(grep(<<<< FILEPATH REDACTED >>>>)) != 0) {
      message(funk)
      source(funk)
    }
  }
  source(<<<< FILEPATH REDACTED >>>>)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

### Total number of data points ############################################################################
if (FALSE) {
  all_data <- fread(<<<< FILEPATH REDACTED >>>>)
  all_data <- all_data[!(country %in% c("CHN"))]
  mbg_data <- all_data[!(country %in% c("BRA", "COK", "DOM", "FJI", "FSM", "HTI", "KIR", "MDV", "MHL", "NCL", "NIU", "NRU", "PLW", "PYF", "TON", "TUV", "VUT", "WLF", "WSM"))]
  nonmbg_data <- all_data[(country %in% c("BRA", "COK", "DOM", "FJI", "FSM", "HTI", "KIR", "MDV", "MHL", "NCL", "NIU", "NRU", "PLW", "PYF", "TON", "TUV", "VUT", "WLF", "WSM"))]
  
  mbg_data %>% nrow
  mbg_data$point %>% table
  mbg_data$diagnostic %>% table
}

### Polygon resampling (Simulated v. actual). Adjumani, Uganda 2015 TAS ########################################

### 1. Load relevant polygon. Copied from polygon resampling script/functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pull_poly_method <- "fast"

admin2_shapefile <- rgdal::readOGR(get_admin_shapefile(admin_level = 2, version = "current"))

shapes <- data.table(shapefile = "lf_g2015_2014_2", location_code = "28381")
shapes <- shapes[!is.na(shapes$shapefile), ]
n_shapes <- nrow(shapes)

# sort by shapefile name (faster if they're all clumped together)
o <- order(shapes$shapefile)
shapes <- shapes[o, ]

# empty, named list to store polygons
polys <- list()
polys[[n_shapes]] <- NA
names(polys) <- paste(shapes$shapefile, shapes$location_code, sep = "__")

# null first shapefile
shape <- ""

# report to the user
message(sprintf("extracting %i polygons", n_shapes))

# vector of shapefiles
shapefiles <- unique(shapes$shapefile)

if (pull_poly_method == "foreach") {
  polys <- pull_polys_foreach(cores, shapefiles, shapes, shp_path)
} else if (pull_poly_method == "mclapply") {
  polys <- pull_polys_mclapply(cores, shapefiles, shapes, shp_path)
} else if (pull_poly_method == "fast") {
  polys <- pull_polys_fast(shapefiles, shapes, shp_path)
} else {
  stop("pull_poly_method must be one of \"foreach\", \"mclapply\", or \"fast\"")
}

problem_shapes <- c()
for (i in 1:length(polys)) {
  problem_shapes <- c(problem_shapes, polys[[i]][["problem_shapes"]])
}
problem_shapes <- problem_shapes[is.na(problem_shapes) == F]

# find ones that didn't work
if (length(problem_shapes) != 0) {
  warning(sprintf("%i polygons could not be found:", length(problem_shapes), ". Dropping these by default. Please check and fix"))
  library(tidyr)
  problem_shapes <- as.data.table(problem_shapes)
  drop_shapes <- separate(problem_shapes, col = "problem_shapes", sep = "__", into = c("shapefile", "location_code")) %>% as.data.table()
  drop_shapes$location_code <- as.integer(drop_shapes$location_code)
  print(drop_shapes)
  if (!ignore_warnings) {
    stop("Since ignore_warnings==FALSE, the function is now breaking. Please fix your bad records or set ignore_warnings==TRUE to drop them.
")
  }
  warning("Since you have set ignore_warnings=TRUE, the function will drop all bad records and continue. Note: this is NOT recommended.
")
  data <- data[!(paste(data$shapefile, data$location_code, sep = "__") %in%
    problem_shapes$problem_shapes), ]
}

polys <- unlist(polys) # get to single-level list
polys[grep("problem_shapes", names(polys))] <- NULL

# Remake IDs to be unique across all of poly_shape_list (needed for rbind below)
makeID_shiny <- function(shape, shape_list) {
  shapefile <- shape_list[[shape]]
  shapefile <- spChFIDs(shapefile, paste0(shape, "_", row.names(shapefile@data)))
  return(shapefile)
}

working_shape_list <- lapply(names(polys), makeID_shiny, shape_list = polys) %>% unlist()
names(working_shape_list) <- names(polys)

## check to make sure that all the projections are the same
proj.strings <- unlist(lapply(working_shape_list, proj4string))
if (length(unique(proj.strings)) > 1) {
  ## take the most common one and assign it to the others
  tab.strings <- table(proj.strings)
  common.proj <- names(tab.strings[which(tab.strings == max(tab.strings))])
  for (ss in c(1:length(proj.strings))) {
    suppressWarnings(proj4string(working_shape_list[[ss]]) <- common.proj)
  }
}

## Combine all into one big SPDF
for (i in 1:nrow(shapes)) {
  working_shape_list[[i]]$shapefile <- shapes[i, shapefile]
  working_shape_list[[i]]$location_code <- shapes[i, location_code]

  working_shape_list[[i]] <- working_shape_list[[i]][which(names(working_shape_list[[i]]) %in% c("shapefile", "location_code"))]
}

poly_shapes_all <- do.call(rbind, working_shape_list)
rm(working_shape_list)

### 2. load and create pseudo-polygon data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_poly <- fread(<<<< FILEPATH REDACTED >>>>)
data_poly <- data_poly[country == "UGA" & point == 0 & location_code == "28381", ]
gaul_list <- get_adm0_codes("UGA", shapefile_version = "current")
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, subset_only = T, buffer = 1, tolerance = 0.4, use_premade = F, shapefile_version = "current")
subset_shape <- simple_polygon_list[[1]]
sub_shape_cropped <- crop(subset_shape, extent(poly_shapes_all) + 1)

## Load list of raster inputs (pop and simple)
raster_list <- build_simple_raster_pop(subset_shape)
simple_raster <- raster_list[["simple_raster"]]
pop_raster <- raster_list[["pop_raster"]]
pop_raster_cropped <- crop(pop_raster[[4]], extent(poly_shapes_all))
pop_raster_cropped <- as.data.frame(pop_raster_cropped, xy = TRUE)
admin2_shapefile <- crop(admin2_shapefile, sub_shape_cropped)

# for bug with NAs & polygon detection when read as character
data_poly$latitude <- as.numeric(data_poly$latitude)
data_poly$longitude <- as.numeric(data_poly$longitude)

### 3. Generate 3 permutations of polygon resampling when treating cluster surveys as a single polygon record ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_poly_data_pop_var1 <- resample_polygons(
  data = data_poly,
  cores = 1,
  indic = "had_lf_w_resamp",
  density = .01, # this is where we can change things up to determine the number of points that get sprinkled
  use_1k_popraster = TRUE,
  gaul_list = gaul_list
)

df_poly_data_pop_var2 <- resample_polygons(
  data = data_poly,
  cores = 1,
  indic = "had_lf_w_resamp",
  density = .01, # this is where we can change things up to determine the number of points that get sprinkled
  use_1k_popraster = TRUE,
  gaul_list = gaul_list
)

df_poly_data_pop_var3 <- resample_polygons(
  data = data_poly,
  cores = 1,
  indic = "had_lf_w_resamp",
  density = .01, # this is where we can change things up to determine the number of points that get sprinkled
  use_1k_popraster = TRUE,
  gaul_list = gaul_list
)

### 4. get corresponding espen data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
espen <- fread(<<<< FILEPATH REDACTED >>>>)
espen_adjumani <- espen[Country == "Uganda" & ADMIN2_NAME == "Adjumani", ]
setnames(espen_adjumani, c("Latitude", "Longitude", "Examined"), c("latitude", "longitude", "N"))
espen_adjumani$latitude <- as.numeric(espen_adjumani$latitude)
espen_adjumani$longitude <- as.numeric(espen_adjumani$longitude)
espen_adjumani <- espen_adjumani[latitude != min(latitude), ] # one misassigned point outside of admin shape. Dropping.

### 5. Make plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A blank theme
theme_empty <- theme_classic() +
  theme(
    axis.line = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# adding admin boundaries to plots
adm2_outline <- suppressMessages(
  geom_path(
    data = crop(admin2_shapefile, sub_shape_cropped),
    aes(
      x = long,
      y = lat,
      group = group
    ),
    color = "white",
    linetype = "dotted"
  )
)

eu_outline <- suppressMessages(
  geom_path(
    data = poly_shapes_all,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    color = "black",
    linetype = "dotted"
  )
)

background_map <- suppressMessages(
  geom_polygon(
    data = sub_shape_cropped,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "light gray"
  )
)

adm0_outline <- suppressMessages(
  geom_path(
    data = sub_shape_cropped,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    color = "black",
    linetype = "solid"
  )
)

diff_polygon <- suppressMessages(
  geom_polygon(
    data = gDifference(admin2_shapefile, poly_shapes_all),
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "light gray"
  )
)

diff_outline <- suppressMessages(
  geom_path(
    data = gDifference(admin2_shapefile, poly_shapes_all),
    aes(
      x = long,
      y = lat,
      group = group
    ),
    color = "black",
    linetype = "solid"
  )
)

map_actualdata <- ggplot() + theme_empty + eu_outline +
  geom_polygon(
    data = poly_shapes_all,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white"
  ) +
  geom_point(
    data = espen_adjumani,
    aes(
      x = as.numeric(longitude),
      y = as.numeric(latitude), color = N, size = N
    ), alpha = .75
  ) +
  scale_color_viridis_c(option = "magma", limits = c(0, max(espen_adjumani$N))) +
  scale_size_continuous(limits = c(0, 220)) +
  geom_path(data = poly_shapes_all, aes(x = long, y = lat, group = group), color = "black", linetype = "solid") +
  ggtitle("Actual survey sample sizes & locations") + coord_equal()

psuedodata <- list(df_poly_data_pop_var1, df_poly_data_pop_var2, df_poly_data_pop_var3)

map_psuedodata <- lapply(psuedodata, FUN = function(dat) {
  map <- ggplot() + theme_empty + eu_outline +
    geom_polygon(
      data = poly_shapes_all,
      aes(
        x = long,
        y = lat,
        group = group
      ),
      fill = "white"
    ) +
    geom_point(
      data = dat,
      aes(
        x = as.numeric(longitude),
        y = as.numeric(latitude), color = N * weight, size = N * weight
      ), alpha = .75
    ) +
    scale_color_viridis_c(option = "magma", limits = c(0, max(espen_adjumani$N))) +
    scale_size_continuous(limits = c(0, 220)) +
    geom_path(data = poly_shapes_all, aes(x = long, y = lat, group = group), color = "black", linetype = "solid") + guides(color = "none", size = "none") + coord_equal()

  return(map)
})

for (i in 1:length(map_psuedodata)) {
  map_psuedodata[[i]] <- map_psuedodata[[i]] + ggtitle(paste0("Resample version ", i))
}

map_pop <- ggplot() + theme_empty + background_map +
  geom_raster(
    data = pop_raster_cropped,
    aes(
      x = x,
      y = y, fill = WorldPop_total_global_stack.4
    )
  ) +
  geom_path(data = poly_shapes_all, aes(x = long, y = lat, group = group), color = "black", linetype = 2) + scale_fill_continuous(type = "viridis") + diff_polygon + diff_outline +
  adm0_outline + ggtitle("Underlying population distribution (WorldPop)") + guides(fill = guide_legend("Population (total)", reverse = T)) + coord_equal()

g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(map_actualdata)

map_actualdata <- ggplot() + theme_empty + eu_outline +
  geom_polygon(
    data = poly_shapes_all,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white"
  ) +
  geom_point(
    data = espen_adjumani,
    aes(
      x = as.numeric(longitude),
      y = as.numeric(latitude), color = N, size = N
    ), alpha = .75
  ) +
  scale_color_viridis_c(option = "magma", limits = c(0, max(espen_adjumani$N))) +
  scale_size_continuous(limits = c(0, 220)) +
  geom_path(data = poly_shapes_all, aes(x = long, y = lat, group = group), color = "black", linetype = "solid") +
  ggtitle("Actual survey sample sizes & locations") + coord_equal() + guides(color = "none", size = "none")

# arrange and plot
grid.arrange(
  grobs = list(map_pop, map_actualdata, mylegend, map_psuedodata[[1]], map_psuedodata[[2]], map_psuedodata[[3]]), widths = c(1, 1, 1, 1, 1, 1),
  layout_matrix = rbind(c(1, 1, 1, NA, 2, 3), c(1, 1, 1, 4, 5, 6)), top = textGrob("Example: Polygon resampling of TAS polygon in Adjumani, Uganda",
    gp = gpar(fontsize = 18)
  )
)
