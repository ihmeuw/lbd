#####################################################################################################################################
### Script objective: Polygon resampling using population rasters
#
#####################################################################################################################################

### Setup #####################################################################################################################

source(<<<< FILEPATH REDACTED >>>>)
load_from_parallelize()
message(paste0("Using ", core_repo))

user <- Sys.info()[["user"]]

## Set repo locations
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- paste0(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

config <- set_up_config_focal_3(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator, 
                                config_file = paste0(<<<< FILEPATH REDACTED >>>>), run_tests = FALSE)

### Some set up
source(paste0(<<<< FILEPATH REDACTED >>>>))

# -------------------------------------------- getting polygon dataset --------------

# Load simple polygon template to model over
gaul_list <- get_adm0_codes(region_list, shapefile_version = modeling_shapefile_version)
message(print(gaul_list))

simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 2, tolerance = 0.4, use_premade = use_premade, shapefile_version = modeling_shapefile_version)
subset_shape <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]

raster_list <- build_simple_raster_pop(subset_shape, link_table = modeling_shapefile_version)

simple_raster <- raster_list[["simple_raster"]]
pop_raster <- raster_list[["pop_raster"]]

year_list <- eval(parse(text = year_list))

## Load input data and make sure it is cropped to modeling area for points
df_pts <- load_input_data(
  indicator = indicator,
  simple = simple_polygon,
  removeyemen = FALSE,
  years = "annual",
  update_run_date = FALSE,
  withtag = as.logical(withtag),
  datatag = datatag,
  use_share = as.logical(use_share)
)

df_pts <- df_pts[, -("keep")]

adm0 <- load_adm0_lookup_table()[gadm_geoid %in% gaul_list]$iso3

# Get polygon data
df_poly <- load_polygon_data(
  indicator = indicator,
  adm0 = adm0,
  withtag = as.logical(withtag),
  datatag = datatag,
  use_share = use_share
)

df_poly[is.na(sampling), sampling := 1]

## locate the path to the pop raster
root <- ifelse(Sys.info()[1] == "Windows", <<<< FILEPATH REDACTED >>>>, <<<< FILEPATH REDACTED >>>>)
if (pop_ras == "default") {
  path.to.pop.rast <- paste0(<<<< FILEPATH REDACTED >>>>)
  p_ras <- NULL
} else {
  pop_ras <- as.character(pop_ras)
  path.to.pop.rast <- pop_ras
  p_ras <- brick(path.to.pop.rast)
}
message(path.to.pop.rast)

# -------------------------------------------- Defining density, specified in config file ----------------
dens <- as.numeric(poly_resamp_density)
message(paste0("Resampling density set to: ", dens))

if (exists("minimum_clusters_per_polygon")) {
  min_clusters <- as.numeric(minimum_clusters_per_polygon)
  message(paste0("You're specified a minimum of ", min_clusters, " clusters for each resampled polygons."))
}

# -------------------------------------------- Resample here ------------------------

df_poly[!is.na(df_poly$latitude) & df_poly$latitude == "", "latitude"] <- NA
df_poly[!is.na(df_poly$longitude) & df_poly$longitude == "", "longitude"] <- NA

df_poly_data_pop <- resample_polygons(
  data = df_poly,
  cores = 1,
  indic = indicator,
  density = dens, # this is where we can change things up to determine the number of points that get sprinkled
  use_1k_popraster = FALSE,
  gaul_list = gaul_list,
  pull_poly_method = pull_poly_method,
  shapefile_version = modeling_shapefile_version,
  ignore_warnings = FALSE
)

if (exists("minimum_clusters_per_polygon")) {
  if (!is.null(minimum_clusters_per_polygon)) {
    # check for number of clusters per polygon
    check_clusters <- table(df_poly_data_pop$Master_UID) %>% as.data.table()
    check_clusters <- check_clusters[order(check_clusters$N, decreasing = T), ]
    setnames(check_clusters, "V1", "Master_UID")
    df_poly_goodenough <- df_poly_data_pop[!(Master_UID %in% check_clusters[N < min_clusters, Master_UID]), ]
    c <- 0
    while (nrow(check_clusters[N < min_clusters, ]) > 0) {
      # subset to any polygons that do not exceed min_clusters
      # rerun polygon resampling iteratively with change in parameters until min_clusters is reached
      
      c <- c + 1
      message(paste0("Iteration ", c, ": changing dens parameter to make sure all polygons are resampled to a minimum of ", min_clusters, " clusters."))
      dens <- dens * 10
      
      message(paste0("Density: ", dens))
      df_poly_data_pop_try <- resample_polygons(
        data = df_poly[!(Master_UID %in% df_poly_goodenough$Master_UID), ],
        cores = 1,
        indic = indicator,
        density = dens, # this is where we can change things up to determine the number of points that get sprinkled
        use_1k_popraster = TRUE,
        gaul_list = gaul_list,
        perpixel = F,
        shapefile_version = modeling_shapefile_version
      )
      check_clusters <- table(df_poly_data_pop_try$Master_UID) %>% as.data.table()
      check_clusters <- check_clusters[order(check_clusters$N, decreasing = T), ]
      setnames(check_clusters, "V1", "Master_UID")
      
      if (c == 5) { # Set maximum iterations at which point, break out of while loop and move on
        df_poly_goodenough <- rbind(df_poly_goodenough, df_poly_data_pop_try)
        
        message(paste0("Maximum iterations have been reached. ", nrow(check_clusters[N < min_clusters, ]), " polygons below min_clusters (printed below). Moving forward."))
        print(unique(df_poly[Master_UID %in% check_clusters[N < min_clusters, Master_UID], .(shapefile, location_code)]))
        break
      }
      
      df_poly_goodenough <- rbind(df_poly_goodenough, df_poly_data_pop_try[!(Master_UID %in% check_clusters[N < min_clusters, Master_UID]), ])
      
      message(paste0(nrow(check_clusters[N < min_clusters, ]), " polygons still below min_clusters threshold."))
    }
  }
}

df_poly_data <- df_poly_goodenough

run_date_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
dir.create(run_date_dir, showWarnings = FALSE)

write.csv(df_poly_data, file = paste0(<<<< FILEPATH REDACTED >>>>))

# --------------------------------------- Join re-sampled polygons to point data ----------------------
# get point data
df_pts <- df_pts[, pseudocluster := FALSE]

# bind together
pts_polys_comb <- rbind.fill(df_pts, df_poly_data)

# subset to only necessary mbg fields
mbg_ready <- copy(pts_polys_comb)
new_indic <- indicator_pts_poly
setnames(mbg_ready, indicator, new_indic)
mbg_ready <- as.data.table(mbg_ready)
if (is.null(mbg_ready$cohort_id)) {
  mbg_ready$cohort_id <- NA
}

if (exists("other_weight")) {
  if (other_weight != "") {
    mbg_ready <- mbg_ready[, c("nid", "Master_UID", "data_collect_method", "country", "source", "cluster_id", "point", "latitude", "longitude", "N", new_indic, "weight", "year", "diagnostic", "shapefile", "location_code", "age_start", "age_end", "cohort_id", other_weight), with = FALSE]
  }
} else {
  mbg_ready <- mbg_ready[, c("nid", "Master_UID", "data_collect_method", "country", "source", "cluster_id", "point", "latitude", "longitude", "N", new_indic, "weight", "year", "diagnostic", "shapefile", "location_code", "age_start", "age_end", "cohort_id"), with = FALSE]
}

# define archive folder path
arch_path <- as.character(archive_data_path)

# define indicator to save as
if (withtag == T) {
  id_save <- paste0(new_indic, datatag)
} else {
  id_save <- new_indic
}

dir.create(arch_path, showWarnings = FALSE)

# archive previous datasets and save
if (use_share == FALSE) {
  if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
    file.copy(paste0(<<<< FILEPATH REDACTED >>>>))
    write.csv(mbg_ready, file = paste0(<<<< FILEPATH REDACTED >>>>))
  } else {
    write.csv(mbg_ready, file = paste0(<<<< FILEPATH REDACTED >>>>))
  }
} else {
  if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
    write.csv(mbg_ready, file = paste0(<<<< FILEPATH REDACTED >>>>))
  } else {
    write.csv(mbg_ready, file = paste0(<<<< FILEPATH REDACTED >>>>))
  }
}

check_resampling(check_resamp_rd = run_date, indicator_group = indicator_group, indicator_poly = indicator)
