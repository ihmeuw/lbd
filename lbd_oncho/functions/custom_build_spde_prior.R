#' @title Build the spde and spde prior
#' 
#' @description construct the spde and relevant priors for Matern GP so that it 
#'   can be used with INLA or TMB
#'   
#' @param spde_prior String containing a list. Specifies the type of prior used
#'   for the Matern model parameters. If type = "nonpc" then use non pc 
#'   prior. User can specify nominal prior means: spde_prior$prior$range.nominal 
#'   and spde_prior$prior$variance.nominal. Defaults are INLA defaults, that is 
#'   20\% of the range of the mesh and 1, respectively. If type = "pc" then use 
#'   pc prior.  User can specify spde_prior$prior$range = (range0, Prange) i.e. 
#'   P(range < range0) = Prange and spde_prior$prior$sigma = (sigma0, Psigma) 
#'   i.e. P(sigma > sigma0) = Psigma. Defaults for Prange and Psigma are 0.05. 
#'   Default for range0 is 5\% of the range of the mesh and sigma0=3.
#'   
#' @param mesh_s An 'inla.mesh' object used to fit the spatial SPDE GP
#'   approximation
#'   
#' @param st_gp_int_zero Logical. Should the space-time GP be forced to interate to 0?
#' 
#' @return List containing the spde and the prior filled in with defaults if no
#'   other values are specified
build_spde_prior <- function(spde_prior, mesh_s, st_gp_int_zero, spde_alpha = 2) {
  if(spde_prior$type=="pc"){ # PC prior
    if(is.null(spde_prior$prior$sigma)) {
      spde_prior$prior$sigma <- c(3, 0.05) # P(sigma > 3) = 0.05
    }
    if(is.null(spde_prior$prior$range)) {
      mesh.range <- max(c(diff(range(mesh_s$loc[, 1])), 
                          diff(range(mesh_s$loc[, 2])), 
                          diff(range(mesh_s$loc[, 3]))))
      spde_prior$prior$range <- c(mesh.range*0.05, 0.05) # P(range < 5% max extent of mesh) = 0.05
    }
    message(paste("Building spde with pc prior,",
                  spde_prior$prior$range[2]*100,
                  "% probability that the range is lower than",
                  spde_prior$prior$range[1], 
                  "and a",
                  spde_prior$prior$sigma[2]*100,
                  "% probability that sigma is greater than",
                  spde_prior$prior$sigma[1]))
    spde <- inla.spde2.pcmatern(mesh = mesh_s, 
                                alpha = spde_alpha,
                                prior.range = spde_prior$prior$range,
                                prior.sigma = spde_prior$prior$sigma,
                                constr = st_gp_int_zero)
  } else { # Non PC prior
    if(is.null(spde_prior$prior$variance.nominal)) {
      spde_prior$prior$variance.nominal <- 1
    }
    spde <- inla.spde2.matern(mesh = mesh_s,  alpha = spde_alpha, constr = st_gp_int_zero,
                              prior.range.nominal = spde_prior$prior$range.nominal,
                              prior.variance.nominal = spde_prior$prior$variance.nominal)
  }
  return(list(spde=spde, spde_prior=spde_prior))
}