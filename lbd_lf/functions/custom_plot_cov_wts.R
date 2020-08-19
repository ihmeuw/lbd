###########################################################################################
# Getting covariate weights for model runs (w/ custom bar graph output)
###########################################################################################

# mask over central plot code with mine (for bar plot)
plot.cov.wts <- function(rd, ## run_date
                         ind, ## indicator
                         ind_gp, ## indicator_group
                         reg,
                         plot.inla.col = TRUE,
                         age = 0,
                         holdout = 0,
                         vallevel = "", ## validation level if used in qsub
                         type = "bar") {

  ## load a cov.wts object
  load(paste0(<<<< FILEPATH REDACTED >>>>))
  cov.wts <- imp.mat
  waic <- waic

  ## remove the final column which has inla weight
  if (!plot.inla.col) {
    cov.wts <- cov.wts[, -ncol(cov.wts)]
  } else {
    ## rescale to be -1:1
    cov.wts[, ncol(cov.wts)] <- cov.wts[, ncol(cov.wts)] / sum(cov.wts[, ncol(cov.wts)], na.rm = TRUE)
  }

  ## melt
  cw.m <- as.data.table(reshape2::melt(as.matrix(cov.wts)))
  colnames(cw.m)[1:3] <- c("Model", "Covar", "Imp")

  ## reorder factors
  cw.m$Model <- factor(cw.m$Model, levels = c("INLA COMBINED", sort(setdiff(unique(cw.m$Model), "INLA COMBINED"))))

  ## setup the plot
  if (type == "bar") {
    require(gridExtra)
    cw.m$Imp <- cw.m$Imp * 100
    cw.m[Model == "INLA COMBINED", Model := "INLA"]
    cw.m$Imp <- round(cw.m$Imp, 1)

    p <- ggplot(data = cw.m[Covar != "SCALED.INLA.COEFS"], aes(x = factor(Covar, levels = unique(sort(Covar))))) +
      geom_histogram(aes(y = Imp, fill = Imp, color = "black"), stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 10)) +
      facet_grid(Model ~ ., scales = "free_y") +
      scale_fill_gradient(high = "red", low = "white") +
      ggtitle(sprintf("%s: covariate importance for %s. WAIC = %s", ind, reg, waic)) +
      xlab("Covariates") +
      guides(color = "none") +
      ylab("Importance")

    p_weights <- ggplot(cw.m[Covar == "SCALED.INLA.COEFS" & Model != "INLA"], aes(x = factor(Model, levels = sort(Model)), y = Imp)) +
      geom_histogram(aes(y = Imp, fill = Imp, color = "black"), stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 10)) + ylab("Importance") +
      guides(color = "none") + xlab("Models") +
      ggtitle(sprintf("%s: model contributions for %s", ind, reg)) +
      scale_fill_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0)

    ## save
    ggsave(p, filename = paste0(<<<< FILEPATH REDACTED >>>>)), width = 10, height = 10, units = "in")
    ggsave(p_weights, filename = paste0(<<<< FILEPATH REDACTED >>>>)), width = 5, height = 5, units = "in")
  } else {
    base_size <- 9
    p <- ggplot(cw.m, aes(Covar, Model)) +
      geom_tile(aes(fill = Imp), colour = "white") +
      scale_fill_gradient2(low = "red", mid = "white", high = "steelblue", midpoint = 0) +
      theme_grey(base_size = base_size) +
      labs(x = "Covariate", y = "Model") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme(axis.text.x = element_text(
        size = base_size * 0.8,
        angle = 270,
        hjust = 0,
        colour = "grey50"
      ))
  }
  return(p)
}
