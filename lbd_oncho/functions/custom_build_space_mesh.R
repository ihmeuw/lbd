## build_space_mesh() ------------------------------------------------------------
#'
#' @title Buils spatial FEM mesh
#'
#' @description Build a finite-elements triangulation mesh to use in
#'   SPDE GP approximation
#'
#' @param d data.frame or data.table of observational data with
#'   'longitude' and 'latitude' columns
#'
#' @param simple simple_polygon defining modeling domain
#'
#' @param max_edge string of 2 element numeric vector in R
#'   notation. e.g. "c(0.25, 5)". First entry decribes the maximum
#'   allowed triangle edge length within the simple_polygon modeling
#'   boundary. Second entry describes the maximum allowed triaqngle
#'   edge length in the extended buffer region. Units are in
#'   lat-long distances. Used if s2mesh=FALSE.
#'
#' @param mesh_offset string of 2 element numeric vector in R
#'   notation. e.g. "c(1, 5)". Describes the automatic extension
#'   distance from the 'simple' boundary. The entries are the inner
#'   and outer extension distances, respectively. Units are in
#'   lat-long distance. Used if s2mesh=FALSE.
#'
#' @param plot_mesh Logical. Should a plot of the mesh be generated?
#'   Not currently implemented for s2mesh
#'
#' @param s2mesh Logical. Should the mesh be created on the surface of
#'   a sphere? If TRUE, s2params is used to specify mesh parameters
#'   instead of max_edge and mesh_offset
#'
#' @param s2params string of 3 element numeric vector in R
#'   notation. e.g. "c(25, 500, 1000)". The entries describe the
#'   minimum triangle edge length allowed, hos far to extend the mesh
#'   beyond the 'simple' boundary, and the maximum allowed triangle
#'   edge length, respectively. Units are in kilometers. Used only if
#'   s2mesh=TRUE.
#'
#' @return an 'inla.mesh' object
#'

build_space_mesh <- function(d, simple, max_edge, mesh_offset,
                             plot_mesh = F, s2mesh = FALSE,
                             s2params = NULL, cutoff){

  if(!as.logical(s2mesh)){ ## build mesh on R2 plane
    message(paste0('Creating spatial mesh, max edge parameter: ', max_edge))
    max.edge <- eval(parse(text=max_edge))
    mesh_offset <- eval(parse(text=mesh_offset))
    mesh_s <- inla.mesh.2d(
      boundary = inla.sp2segment(simple),
      loc = cbind(d$longitude,d$latitude),
      max.edge = max.edge,
      offset = mesh_offset,
      cutoff = cutoff
    )

    if (plot_mesh) {
      plot(mesh_s, asp=1)
      points(d$longitude, d$latitude, col=d$year)
    }

  } else{ ## build mesh on sphere surface

    if(is.null(s2params)){
      stop("You've chosen to build an s2 mesh but haven't provided parameters to do so (i.e. s2params=NULL)!")
    }

    message(paste0('Creating spatial SPHERICAL mesh, min edge, max edge, extension kms are: ', s2params))
    s2params <- eval(parse(text = s2params))

    ## convert data locs to 3d coords on the unit-sphere
    true.radius.of.earth = 6371
    s3 <- lonlat3D(d$longitude, d$latitude)

    ## convert boundary of simple_polygon to 3d coords on unit-sphere
    boundary <- inla.sp2segment(simple)
    boundary.loc <- lonlat3D(boundary$loc[, 1], boundary$loc[, 2])

    ## make a s2 domain mesh using data locs and boundary locs
    all.loc <- rbind(s3, boundary.loc)
    mesh_s <- inla.mesh.create(loc=all.loc,
                               cutoff = s2params[1] / true.radius.of.earth, ## minimum triangle edge allowed
                               extend = list(offset = s2params[2] / true.radius.of.earth), ## how far to extend mesh
                               refine=list(max.edge = s2params[3] / true.radius.of.earth)) ## max triangle edge allowed

    ## TODO add 3d mesh plotting
    ## an example of how to do this can be found in the code examples here:
    ## http://www.r-inla.org/examples/case-studies/simpson2011
    if(plot_mesh){
      ## plot 3d mesh
    }
  }

  return(mesh_s)

}
