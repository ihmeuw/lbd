library(rgl)
library(INLA)
library(rgdal)

lonlat3D <- function(lon, lat){
  cbind(cos((lon/180)*pi)*cos((lat/180)*pi),
        sin((lon/180)*pi)*cos((lat/180)*pi),
        sin((lat/180)*pi))
}

## draw.mesh  ---------------------------------------------------
# this function allows 3d drawing and imaging of a mesh.inla
# object. it can plot either a mesh constructed on the R2 or S2
# manifold. NOTE that it requires the 'rgl' R package and it spawns
# interactive graphics windows. it is untested on the cluster and is
# meant for use on local machines
#
#' @description Plots, in 3d, a s2-manifold or r2-manifold mesh 
#'
#' @param mesh an inla.mesh object 
#'
#' @param draw.edges Logical. Draw the edges between the vertices?
#'
#' @param draw.segments Logical. Draw the segments that bound the mesh
#'   object?
#'
#' @param draw.plane Logical. Draw a planar shape to aid in displaying
#'   curvature of mesh?
#'
#' @param node.cols Numeric vector with length equal to the number of
#'   mesh vertices (mesh$n). e.g. pass in the posterior mean of the
#'   spatial random effects heights to visualize the fitted GP.
#'
#' @param col.type String taking value of either 'bw' or 'col' and
#'   determining whether the color of the surface should be drawn in
#'   black and white or in color. Only used if a node.cols vector is
#'   passed in.
#'
#' @param window.dims 2d numeric vector describing the width and
#'   height of the plotting in pixels
#'
#' @param plot.mirror Logical. Should a mirror image of the mesh be
#'   added to the plot (could help visualize mesh if printing to a
#'   static image)
#'
#' @param returns nothing but spawns an interactive plot window
#' 
#' @examples
#'
#' ## plot a mesh. add color to the background that results from
#' ## the linear interpolation of randomly generated nodal (basis
#' ## height) values
#' draw.s2.mesh(mesh_s, draw.edges = T, draw.segments = T, col.type = 'col',
#'              node.cols = rnorm(n = mesh_s$n), draw.plane = F)
#'
#' ## take a snapshot to save to file
#' fig.path <- <<<< FILEPATH REDACTED >>>>
#' rgl.snapshot(file.path(fig.path, "mesh.png"), top=TRUE)
#' 
#' ## shut down the graphics window
#' rgl.close()

draw.mesh <- function(mesh, draw.edges=TRUE, draw.segments=TRUE,
                      draw.plane = F, node.cols = NULL, 
                      col.type = 'bw', window.dims = c(840, 840),
                      plot.mirror = FALSE, shape=NULL){
  
  require('rgl') ## this is an interactive R graphics. won't work on the cluster
  
  window.dims = c(50, 50, 50 + window.dims[1], 50 + window.dims[2])
  
  if(is.null(node.cols)){
    node.cols <- rep(1, mesh$n)
  }
  
  if(col.type == 'bw')  cp <- function (n, ...) { return (grey.colors(n, start=0.05, end=0.95, alpha=0.05, ...))}
  if(col.type == 'col') cp <- colorRampPalette(c("darkblue", "blue", "cyan",
                                                 "yellow", "red", "darkred"))
  
  mesh0 = inla.mesh.create(loc=cbind(0,0), extend=list(offset=1.1,n=4))
  
  mesh01 = mesh0
  mesh02 = mesh0
  mesh1 = mesh
  mesh2 = mesh
  mesh02$loc[,1] = mesh02$loc[,1]*(-1)
  mesh02$loc[,3] = mesh02$loc[,3]*(-1)
  mesh2$loc[,1] = mesh2$loc[,1]*(-1)
  mesh2$loc[,3] = mesh2$loc[,3]*(-1)
  
  mesh01$loc[,1] = mesh01$loc[,1]-1.1
  mesh02$loc[,1] = mesh02$loc[,1]+1.1
  mesh1$loc[,1] = mesh1$loc[,1]-1.1
  mesh2$loc[,1] = mesh2$loc[,1]+1.1
  
  open3d(windowRect=window.dims)
  if(draw.plane){
    plot(mesh01, rgl=TRUE, col="white", color.palette=cp,
         draw.vertices=FALSE, draw.edges=FALSE, add=TRUE)
    if(plot.mirror){
      plot(mesh02, rgl=TRUE, col="white", color.palette=cp,
           draw.vertices=FALSE, draw.edges=FALSE, add=TRUE)
    }
  }
  
  plot(mesh1, rgl=TRUE, col=node.cols, color.palette=cp,
       draw.vertices=FALSE, draw.edges=draw.edges, add=TRUE,
       draw.segments=draw.segments, edge.color="#999999")
  if(plot.mirror){
    plot(mesh2, rgl=TRUE, col=node.cols, color.palette=cp,
         draw.vertices=FALSE, draw.edges=draw.edges, add=TRUE,
         draw.segments=draw.segments)
  }
  
  for (a in 1:length(shape@polygons)) {
    for (b in 1:length(shape@polygons[[a]]@Polygons)) {
      polygon3d(x=shape@polygons[[a]]@Polygons[[b]]@coords[, 1] - 1.1, y=shape@polygons[[a]]@Polygons[[b]]@coords[, 2], z=shape@polygons[[a]]@Polygons[[b]]@coords[, 3], fill=FALSE, random=FALSE, col="black", lwd=2.5)
    }
  }
  
  view3d(0,0,fov=0,zoom=0.8)
  rgl.bringtotop()
}

#### Load admin shapefile with national borders
adm0 <- readOGR(<<<< FILEPATH REDACTED >>>>) # Replace this path as needed

#### Load mesh object from INLA run
load(<<<< FILEPATH REDACTED >>>>) # Replace this path as needed

#### Subset admin0 shapefile to desired countries
## Replace list as appropriate for the model region
shape <- adm0[adm0@data$ADM0_CODE %in% c(3, 18, 20, 22, 39, 45, 46, 47, 48, 72, 79, 83, 85, 87, 88, 89, 118, 127, 146, 152, 157, 162, 164, 191, 193, 194, 200, 203, 206, 207, 217, 218, 230, 231, 253),] # oncho_endem_afr

#### Project shapefile coordinates onto a sphere (this can be very slow)
for (a in 1:length(shape@polygons)) {
  for (b in 1:length(shape@polygons[[a]]@Polygons)) {
    print(paste0(a, " ", b))
    shape@polygons[[a]]@Polygons[[b]]@coords <- cbind(shape@polygons[[a]]@Polygons[[b]]@coords, NA)
    for (c in 1:nrow(shape@polygons[[a]]@Polygons[[b]]@coords)) {
      shape@polygons[[a]]@Polygons[[b]]@coords[c,] <- lonlat3D(shape@polygons[[a]]@Polygons[[b]]@coords[c, 1], shape@polygons[[a]]@Polygons[[b]]@coords[c, 2])
    }
  }
}

draw.mesh(mesh=mesh_s, col.type="bw", shape=shape)

rgl.postscript("mesh", fmt = "pdf", drawText = TRUE)

