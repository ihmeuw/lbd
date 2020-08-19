build_spde <- function(mesh_s, sig0=0.5, rho0=0.3) {
  
    # construct an SPDE model with a Matern kernel
    message('Building SPDE...')

    spde <- inla.spde2.pcmatern(mesh=mesh_s, prior.range=c(rho0, 0.5), prior.sigma=c(sig0, 0.5))
    
    return(spde)
}
