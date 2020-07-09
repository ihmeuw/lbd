require(INLA)
require(data.table)
require(ggplot2)

plot_hyperparameters <- function(indicator, indicator_group, run_date, age, holdout, save_file = NULL, regs) {
  
  # get regions
  outputdir <- '<<<< FILEPATH REDACTED >>>>'
  regions <- regs
  
  # extract prior and posterior distributions from INLA model objects
  message('Load models & extract priors and posteriors for hyper-parameters')
  dist <- rbindlist(lapply(regions, function(r) {
    
    # load model
    message(paste0('...', r))
    result = tryCatch({
      load(paste0(outputdir, indicator, '_model_eb_bin', age, '_', r, '_', holdout, '.RData'))
    }, error = function(e) {
      return()
    })
    
    # extract hyper-priors from INLA (based on plot.inla() code)
    all.hyper <- INLA:::inla.all.hyper.postprocess(res_fit$all.hyper)
    hyper <- res_fit$marginals.hyperpar
    id <- strsplit(sapply(hyper, attr, 'hyperid'), split = '\\|')
    prior <- rbindlist(lapply(names(id), function(x) {
      print(x)
      if (grepl('Theta. for', x)) range <- c(-5, 5)
      if (grepl('GroupRho for', x)) range <- c(-0.999, 0.999)
      if (grepl('Group PACF. for', x)) range <- c(-0.999, 0.999)
      if (grepl('Precision for', x)) range <- c(1, 1000)
      p <- INLA:::inla.get.prior.xy(section = tolower(id[[x]][2]), hyperid = id[[x]][1], all.hyper = all.hyper, range = range, intern = F)
      if (grepl('Precision for', x)) {
        p <- inla.tmarginal(function(x) sqrt(1/x), p, method = 'linear')
        x <- gsub('Precision for', 'SD for', x)
      }
      data.table(region = r, type = 'prior', name = x, x = p$x, y = p$y)
    }))
    
    # extract corresponding posteriors from INLA
    post <- rbindlist(lapply(names(hyper), function(x) {
      p <- hyper[[x]]
      if (grepl('Precision for', x)) {
        p <- inla.tmarginal(function(x) sqrt(1/x), p, method = 'linear')
        x <- gsub('Precision for', 'SD for', x)
      }
      data.table(region = r, type = 'posterior', name = x, x = p[, 'x'], y = p[, 'y'])
    }))
    
    # combine
    all <- rbind(prior, post)
    all[, name := factor(name, unique(name))]
    return(all)
  }))
  
  # make plots
  message('Plotting hyper-parameters')
  
  if (is.null(save_file)) save_file <- paste0(outputdir, '/diagnostic_plots/inla_hyperparameters.pdf')
  pdf(save_file, width = 14, height = 8)
  gg <- ggplot(dist[y > 1e-8,], aes(x = x, y = y, color = region, linetype = type)) +
    facet_wrap(~ name, scales = 'free') +
    geom_line() +
    labs(x = '', y = '', title = paste0('Hyper-parameter prior and posterior distributions for ', indicator)) +
    theme_bw()
  print(gg)
  dev.off()
  
  return(dist)
}