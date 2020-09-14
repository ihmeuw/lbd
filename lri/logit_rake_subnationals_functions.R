################################################
#functions for logit raking
################################################

#install faraway if on cluster
personal_lib <- ifelse(grepl('health_fin', Sys.getenv("SINGULARITY_NAME")) ,
                       '<<<< FILEPATH REDACTED >>>>',
                       '<<<< FILEPATH REDACTED >>>>')
if(!dir.exists(personal_lib)) dir.create(personal_lib, recursive = TRUE)
Sys.setenv(R_LIBS_USER = personal_lib)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()) )
install.packages('faraway')

require(faraway) #for ilogit

G <- function(x){
	val <- logit(x)
	return(val)
}

G_Inv <- function(x){
	val <- ilogit(x)
	return(val)
}

Mult <- FALSE

Agg <- function(p_i, N_i){
	return(sum(p_i * N_i)/sum(N_i))
}


EvalK <- function(K, p_i, p_N, N_i){
	if (Mult == TRUE){
		p_tilde <- G_Inv(G(p_i) * K)
	} else {
		p_tilde <- G_Inv(G(p_i) + K)
	}
	LHS <- Agg(p_tilde, N_i)
	RHS <- p_N
	SE <- (LHS-RHS)^2
	return(SE)
}

FindK <- function(p_i, p_N, N_i, a) {
	if (Mult == TRUE){
		Limits <- c(0,a)
	} else {
		Limits <- c(-a,a)
	}
	iter <- 1
	Boundary = TRUE
	while (Boundary & iter < 10){
		Limits <- Limits * iter
		val <- optimize(EvalK, Limits, p_i = p_i, p_N = p_N, N_i = N_i, tol = 1e-20)$min
		Boundary <- (round(abs(val)) == Limits[2])
		iter <- iter + 1
	}
	return(val)
}

## Load list of raster inputs (pop and simple), with all years of population
build_simple_raster_pop_all_yrs <- function(subset_shape, u5m=FALSE, field=NULL, raking=F, link_table=modeling_shapefile_version) {
  
  if (is.null(field)) {
    if ('GAUL_CODE' %in% names(subset_shape@data)) field <- 'GAUL_CODE'
    if ('ADM0_CODE' %in% names(subset_shape@data)) field <- 'ADM0_CODE'
  }
  
  if(raking){
    field <- 'loc_id'
    # no 'loc_id' field in the link table, so we can't use it
    link_table <- NULL
  }
  
  if(u5m==FALSE){
    master_pop <- brick('<<<< FILEPATH REDACTED >>>>')
  } else {
    master_pop <- brick('<<<< FILEPATH REDACTED >>>>')
  }
  
  cropped_pop <- crop(master_pop, extent(subset_shape), snap="out")
  ## Fix rasterize
  initial_raster <- rasterize_check_coverage(subset_shape, cropped_pop, field = field, link_table = link_table)
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster))))!=0) {
    rasterized_shape <- 
      raster::merge(
        rasterize_check_coverage(subset(subset_shape, !(get(field) %in% unique(initial_raster))),
                                 cropped_pop,
                                 field = field,
                                 link_table = link_table),
        initial_raster)
  }
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster))))==0) {
    rasterized_shape <- initial_raster
  }
  masked_pop <- raster::mask(x=cropped_pop, mask=rasterized_shape)
  
  raster_list <- list()
  raster_list[['simple_raster']] <- rasterized_shape
  raster_list[['pop_raster']] <- masked_pop
  
  return(raster_list)
  
}