###############################################################
# Purpose: Generate Quick and Dirty Maps using MapSuite
# based on RasterStack objects
# Note: Must be run where MapSuite library exists
# for maps to be generated.
##############################################################

#####################################################################
# Prepare Workspace
#####################################################################

# Setting the root, loading in libraries, and functions:
setwd("<<<< FILEPATH REDACTED >>>>")
for(function_script in list.files(getwd(),pattern="*_functions.R")){message(function_script);source(function_script)};message("Central Functions Loaded.")

load_libs(c('data.table',"MapSuite","stringr","raster"))

# Getting parameters from arguments passed in the qsub call
indicator_group<-commandArgs()[4]; message(indicator_group)
indicator<-commandArgs()[5]; message(indicator)
run_date<-commandArgs()[6]; message(run_date)
rake<-commandArgs()[7]; message(rake)
results_timestamp<-commandArgs()[8]; message(results_timestamp)

results<-fread("<<<< FILEPATH REDACTED >>>>")

indicator_name<-results[indicator_name==indicator &
                        raked==rake]$indicator_longname[1]
#####################################################################
# Load in raster data (results raster, mask)
#####################################################################

# Pulling results .tif into R
results_tif<-raster::stack("<<<< FILEPATH REDACTED >>>>")

# Load in raster mask
mask_layer<-raster::raster("<<<< FILEPATH REDACTED >>>>")

# Mask results raster based on Geospatial's mask layer
results_tif<-raster::mask(results_tif,mask_layer,maskvalue=1)

# Converting raster to map and data tables objects:
results_tables<-MapSuite::PixelsToTables(results_tif,dimensions = seq(2000,2015),dim_name="year",var_name="mean")
map<-results_tables$map
data<-results_tables$data

################################################################
# Load in Polygon boundaries
################################################################

# Prepare the shapefile you want to use as an outline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #               <---------------------- Change these parameters to alter the Admin unit boundaries
subset_shp<-"gauls" #' This can be either "gauls" (which will take a
#' list of guals, (in this case, Africa), and
#' subset based on the ADM0 gauls in that list,
#' or this can be set to "crop", which will
#' include all boundaries within the raster's extent,
#'  although they may not be included in the analysis

admin_level<-0      #' You can set this to whatever level you want (0,1,2...)-
#' MapSuite will only allow 1 shapefile as an outline as part
#' of the RasterMap function (there are fancy ways to have more than
#' 1 outline, or other shapes, but these go beyond the simple
#' RasterMap function-- see MapSuite documentation)

# Load in an Admin shapefile (level of your choice)
admin_shp<-rgdal::readOGR(dsn="<<<< FILEPATH REDACTED >>>>",
                          layer=paste0("g2015_2014_", admin_level))

if(subset_shp=="gauls"){
  # List of Africa GAULS
  afr_adm0_gaul<-c(4,6,8,29,35,42,43,45,47,49,50,58,59,66,68,70,
                   74,76,77,79,89,90,94,95,105,106,142,144,145,150,
                   152,155,159,169,170,172,181,182,205,214,217,221,
                   226,235,243,248,253,268,270,271,40762,40765,
                   227,257,133) # taken from load_simple_polygon mbg central code; no Yemen.
  admin_shp<-admin_shp[admin_shp@data$ADM0_CODE %in% afr_adm0_gaul,]
}else if(subset_shp=="crop"){
  admin_shp<-raster::crop(admin_shp,results_tif)
}

################################################################
# Map the Data
################################################################

# Prepare the data to be mapped
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Manipulating numbers to be percentages
data[,mean:=mean*100]

# You can make a  PDF of your maps as well, with a histogram at the bottom to show trends over time:
pdf_path<-"<<<< FILEPATH REDACTED >>>>"
pdf(file=pdf_path,width=6,height=10) # Starting a PDF of the results

message("Mapping time series to PDF.")

MapSuite::RasterMap(coords=map,data=data,
          id="id", xcol="x", ycol="y",
          outline = admin_shp,
          series_dimension = "year",
          series_sequence = c(2000,2005,2010,2015),
          histogram=T,
          outline_color = "grey",
          variable="mean",
          map_title=indicator_name,
          map_colors = wpal("sky"),
          legend_title="Mean Estimate")

dev.off()

message("Maps should be done; go look for the PDF!")
