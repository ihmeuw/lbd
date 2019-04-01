###############################################################################
###############################################################################
## Create plots comparing model specifications (master script)
##
## Purpose: Create plots comparing OOS validation metrics for various 
##          combinations of stacking vs raw covariates with and without 
##          Gaussian process. Calls `plot_of_plots_parallel.R`
###############################################################################
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<<<< FILEPATH REDACTED >>>>'
indic_repo         <- '<<<< FILEPATH REDACTED >>>>'

## sort some directory stuff
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific code begins here ##########################################

library(ggrepel)

indicator_group <- 'vaccine'
indicators <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout")
years_to_plot <- c(2000, 2005, 2010, 2016)
drop_models <- NULL
samples <- 1000

## Set up run dates (one for each model specification)

raw <- NULL         # Raw covariates only
stack <- NULL       # Stacked covariates only
gp <- NULL          # Gaussian process only
raw_gp <- NULL      # Raw covariates with GP
stack_gp <- NULL    # Stacked covariates with GP

run_date <- stack_gp  # A run date is needed for parallelize()

# All the relevant run_dates
all.rd <- c(raw, stack, gp, raw_gp, stack_gp)

pp_lv <- list(indicator = indicators,
              run_date = all.rd, 
              samples = samples)

pp_qsub_output <- parallelize(script = "plot_of_plots_parallel",
                              script_dir = paste0(indic_repo, "share_scripts"),
                              log_location = "sgeoutput",
                              expand_vars = pp_lv,
                              save_objs = c("indicator_group"),
                              prefix = "pp",
                              slots = 4,
                              memory = 20,
                              geo_nodes = TRUE, 
                              singularity = 'default')

monitor_jobs(pp_qsub_output)

## ##################################################################################################
## ##################################################################################################

for(ind in c('dpt3_cov', 'dpt1_cov', 'dpt1_3_abs_dropout')){

  indicator <- ind

  message(paste0("Working on indicator ", ind))

  ## get ready to pull stats
  pull_OOS_stats <- function(i, rd, model_type, oos = "TRUE", agg = 'ad2', ind_group = indicator_group) {
    stats <- fread(paste0('<<<< FILEPATH REDACTED >>>>/', ind_group, '/', i, '/output/', rd, '/si_figs/', agg, '_metrics.csv'))
    stats <- stats[OOS == oos, ]
    stats <- stats[, Model := model_type]
    names(stats)[names(stats)=='95% Cov.'] <- 'Coverage'
    names(stats)[names(stats)=='Mean Err.'] <- 'Bias'
    
    return(stats)
  }

##### In-sample (IS) statistics
  raw_stats <- pull_OOS_stats(ind, raw, 'Raw covs', "FALSE", 'ad2')
  stack_stats <- pull_OOS_stats(ind, stack, 'Stacked covs', "FALSE", 'ad2')
  gp_stats <- pull_OOS_stats(ind, gp, 'GP', "FALSE", 'ad2')
  raw_gp_stats <- pull_OOS_stats(ind, raw_gp, 'Raw + GP', "FALSE", 'ad2')
  stack_gp_stats <- pull_OOS_stats(ind, stack_gp, 'Stacked + GP', "FALSE", 'ad2')

  IS_stats <- rbind(raw_stats, 
                    stack_stats, 
                    gp_stats, 
                    raw_gp_stats, 
                    stack_gp_stats,
                    use.names = T)

  IS_stats <- subset(IS_stats, !(Model %in% drop_models))

  ## Plot
  IS_stats$Model <- factor(IS_stats$Model, levels = c('Raw covs', 'Stacked covs', 'GP', 'Raw + GP', 'Stacked + GP'))

  annotations <- data.frame(
    xpos = c(-Inf , -Inf, 0   , Inf),
    ypos =  c(0.95, 0   , 0.95, Inf),
    annotateText = rep("Best\nPerformance",4), 
    hjustvar = c(0.0, 0.0, 0.0, 1) ,
    vjustvar = c(1.0, 1.0, 1.0, 1))
  
  x_lim <- c(min(IS_stats$RMSE), max(IS_stats$RMSE))
  y_lim <- c(min(IS_stats$Coverage), max(IS_stats$Coverage))

  gg.is.1 <- ggplot(data = IS_stats,
                aes(x = RMSE,
                    y = Coverage)) +
    geom_point(size = 5, shape = 21, aes(fill = Model)) +
    theme_light() + 
    scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')) + 
    ggtitle('In-sample statistics aggregated to Admin 2') +
    geom_text_repel(aes(RMSE, Coverage, label = Model)) +
    geom_text(data = annotations[1,], aes(x=xpos,y=ypos,hjust=hjustvar,
                                          vjust=vjustvar,label=annotateText, col = 'red'),
               show.legend = FALSE) +
    geom_hline(yintercept = 0.95, color = 'red') +
    ylim(c(min(y_lim[1], 0.9), max(y_lim[2], 1)))


  x_lim <- c(min(IS_stats$RMSE), max(IS_stats$RMSE))
  y_lim <- c(min(IS_stats$Bias), max(IS_stats$Bias))
    
  gg.is.2 <- ggplot(data = IS_stats,
                aes(x = RMSE,
                    y = Bias)) +
    geom_point(size = 5, shape = 21, aes(fill = Model)) +
    theme_light() + 
    scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')) + 
    ggtitle('In-sample statistics aggregated to Admin 2') +
    geom_text_repel(aes(RMSE, Bias, label = Model)) +
    geom_text(data = annotations[2, ], aes(x=xpos,y=ypos,hjust=hjustvar,
                                           vjust=vjustvar,label=annotateText, col = 'red'),
               show.legend = FALSE) +
    geom_hline(yintercept = 0.00, color = 'red') + 
    ylim(c(min(y_lim[1], -0.05), max(y_lim[2], 0.05)))

  x_lim <- c(min(IS_stats$Bias), max(IS_stats$Bias))
  y_lim <- c(min(IS_stats$Coverage), max(IS_stats$Coverage))

  gg.is.3 <- ggplot(data = IS_stats,
                aes(x = Bias,
                    y = Coverage)) +
    geom_point(size = 5, shape = 21, aes(fill = Model)) +
    theme_light() + 
    scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')) + 
    ggtitle('In-sample statistics aggregated to Admin 2') +
    geom_text_repel(aes(Bias, Coverage, label = Model)) +
    geom_text(data = annotations[3, ], aes(x=xpos,y=ypos,hjust=hjustvar,
                                           vjust=vjustvar,label=annotateText, col = 'red'),
               show.legend = FALSE) +
    geom_vline(xintercept = 0.00, color = 'red') + geom_hline(yintercept = 0.95, color = 'red') +
    xlim(c(min(x_lim[1], -0.05), max(x_lim[2], 0.05))) +
    ylim(c(min(y_lim[1], 0.9), max(y_lim[2], 1)))
  

##### Out-of-sample (OOS) statistics
  raw_stats <- pull_OOS_stats(ind, raw, 'Raw covs', "TRUE", "ad2")
  stack_stats <- pull_OOS_stats(ind, stack, 'Stacked covs', "TRUE", "ad2")
  gp_stats <- pull_OOS_stats(ind, gp, 'GP', "TRUE", "ad2")
  raw_gp_stats <- pull_OOS_stats(ind, raw_gp, 'Raw + GP', "TRUE", "ad2")
  stack_gp_stats <- pull_OOS_stats(ind, stack_gp, 'Stacked + GP', "TRUE", "ad2")

  OOS_stats <- rbind(raw_stats, 
                     stack_stats, 
                     gp_stats, 
                     raw_gp_stats, 
                     stack_gp_stats)

  OOS_stats <- subset(OOS_stats, !(Model %in% drop_models))

  ## Plot
  OOS_stats$Model <- factor(OOS_stats$Model, levels = c('Raw covs', 'Stacked covs', 'GP', 'Raw + GP', 'Stacked + GP'))
  annotations <- data.frame(
    xpos = c(-Inf , -Inf, 0   , Inf),
    ypos =  c(0.95, 0   , 0.95, Inf),
    annotateText = rep("Best\nPerformance",4), 
    hjustvar = c(0.0, 0.0, 0.0, 1) ,
    vjustvar = c(1.0, 1.0, 1.0, 1))
  
  x_lim <- c(min(OOS_stats$RMSE), max(OOS_stats$RMSE))
  y_lim <- c(min(OOS_stats$Coverage), max(OOS_stats$Coverage))

  gg.oos.1 <- ggplot(data = OOS_stats,
                aes(x = RMSE,
                    y = Coverage)) +
    geom_point(size = 5, shape = 21, aes(fill = Model)) +
    theme_light() + 
    scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')) + 
    ggtitle('Out-of-sample statistics aggregated to Admin 2') +
    geom_text_repel(aes(RMSE, Coverage, label = Model)) +
    geom_text(data = annotations[1,], aes(x=xpos,y=ypos,hjust=hjustvar,
                                          vjust=vjustvar,label=annotateText, col = 'red'),
              show.legend = FALSE) +
    geom_hline(yintercept = 0.95, color = 'red') +
    ylim(c(min(y_lim[1], 0.9), max(y_lim[2], 1)))

  x_lim <- c(min(OOS_stats$RMSE), max(OOS_stats$RMSE))
  y_lim <- c(min(OOS_stats$Bias), max(OOS_stats$Bias))

  gg.oos.2 <- ggplot(data = OOS_stats,
                     aes(x = RMSE,
                         y = Bias)) +
    geom_point(size = 5, shape = 21, aes(fill = Model)) +
    theme_light() + 
    scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')) + 
    ggtitle('Out-of-sample statistics aggregated to Admin 2') +
    geom_text_repel(aes(RMSE, Bias, label = Model)) +
    geom_text(data = annotations[2, ], aes(x=xpos,y=ypos,hjust=hjustvar,
                                           vjust=vjustvar,label=annotateText, col = 'red'),
               show.legend = FALSE) +
    geom_hline(yintercept = 0.00, color = 'red') +
    ylim(c(min(y_lim[1], -0.05), max(y_lim[2], 0.05)))

  x_lim <- c(min(OOS_stats$Bias), max(OOS_stats$Bias))
  y_lim <- c(min(OOS_stats$Coverage), max(OOS_stats$Coverage))

  gg.oos.3 <- ggplot(data = OOS_stats,
                     aes(x = Bias,
                         y = Coverage)) +
    geom_point(size = 5, shape = 21, aes(fill = Model)) +
    theme_light() + 
    scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')) + 
    ggtitle('Out-of-sample statistics aggregated to Admin 2') +
    geom_text_repel(aes(Bias, Coverage, label = Model)) +
    geom_text(data = annotations[3, ], aes(x=xpos,y=ypos,hjust=hjustvar,
                                           vjust=vjustvar,label=annotateText, col = 'red'),
               show.legend = FALSE) +
    geom_vline(xintercept = 0.00, color = 'red') + geom_hline(yintercept = 0.95, color = 'red') +
    xlim(c(min(x_lim[1], -0.05), max(x_lim[2], 0.05))) +
    ylim(c(min(y_lim[1], 0.9), max(y_lim[2], 1)))
  

  pdf(paste0('<<<< FILEPATH REDACTED >>>', user, '/plot_of_plots_cov_', ind, '.pdf'), width = 30, height = 10)
  grid.arrange(gg.is.1, gg.is.2, gg.is.3, ncol = 3)
  grid.arrange(gg.oos.1, gg.oos.2, gg.oos.3, ncol = 3)
  dev.off()

  png(paste0('<<<< FILEPATH REDACTED >>>', user, '/plot_of_plots_cov_', ind, '_is.png'), width = 3000, height = 1000, res = 120)
    grid.arrange(gg.is.1, gg.is.2, gg.is.3, ncol = 3)
  dev.off()

  png(paste0('<<<< FILEPATH REDACTED >>>', user, '/plot_of_plots_cov_', ind, '_oos.png'), width = 3000, height = 1000, res = 120)
    grid.arrange(gg.oos.1, gg.oos.2, gg.oos.3, ncol = 3)
  dev.off()

}