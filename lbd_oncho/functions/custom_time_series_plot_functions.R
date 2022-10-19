## time_series plots  #############################

# Functions to create time series plots at the specified admin level, faceting by the number of admin levels in that geography
#
# This function is called in the subnational_ts_function and is a way to generalize making time
# series for each admin unit. These functions are all very similar but come in four flavors:
#
#
# 1) time_series:  This plots the model predictions along with the uncertainty, faceted by the admin unit
# 2) time_series_data: This plots the model predictions and uncertainty but adds the data aggregated to that admin level
# (using input_aggregated_admin function), and it colors the points by the data source, as well as differentiating between
# survey size and point and polygon with size and shape aesthetics respectively
# 3) time_series_multiple: This plots multiple model predictions (either different model runs or multiple indicators, specified by model_runs = T in subnational_ts_plots)
# on each time trend at the specified admin unit (0/1/2)
# 4) time_series_multiple_data: This plots multiple model predictions and also includes the aggregated input data to that admin unit.
# I only recommend using this if comparing different model runs of the same indicator, as plotting the data and multiple indicators makes the plots too busy
#
# @author Michael Cork, \email{mcork23@uw.edu}
#
# @param model_data admin level data table from the standard aggregation code (typically ad0_df/ad1_df/ad2_df specified to certain geographic area)
# @param input_data input data collapsed to admin level (typically ad0_data/ad1_data/ad2_data specified to certain geographic area)
# @param admin specified admin level (0/1/2)
# @param val_range Range of values for y axis
# @param title_plot_size Size of plot title
#
# @return
# This function returns a ggplot object
#
# @examples
#
# ad0_df_reg is the ad0_df specified to specific region, say for example 'cssa' (Central Sub-Saharan Africa)
# if (plot_data == F & multiple_runs == F) ad0_time_series <- time_series(ad0_df_reg, admin = 0, val_range = c(0,1))

time_series <- function(model_data,
                        admin = 0,
                        val_range = c(0, 1),
                        title_plot_size = 30,
                        ind_title = "",
                        mda = NULL,
                        plot_stkrs = plot_stacker_trends,
                        stkrs = stackers,
                        s_fixed_effects = stacked_fixed_effects) {

  # Color palette
  carto_discrete <- c(
    "#7F3C8D", "#11A579", "#F2B701", "#E73F74",
    "#3969AC", "#80BA5A", "#E68310", "#008695",
    "#CF1C90", "#f97b72", "#4b4b8f", "#A5AA99"
  )

  if (is.null(mda) == T) {
    gg_admin <-
      ggplot(data = model_data, aes(x = year)) +
      geom_ribbon(alpha = 0.3, aes(ymin = lower, ymax = upper)) +
      geom_line(aes(y = mean)) +
      theme_bw(base_size = 16) +
      coord_cartesian(ylim = val_range) +
      facet_wrap(~plot_name) +
      theme(
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0.5),
        plot.title = element_text(size = title_plot_size, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        y = ind_title,
        x = "Year"
      ) + {
        if (admin == 0) labs(title = paste0(ind_title, " by Country"))
      } + {
        if (admin == 1) {
          labs(
            title = paste0(ind_title, " by First-level Administrative Unit"),
            caption = paste0(
              "Time series depict first-level administrative units except for NATIONAL,\n",
              "which shows the time series for the entire country"
            )
          )
        }
      } + {
        if (admin == 2) {
          labs(
            title = paste0(ind_title, " by Second-level Administrative Unit"),
            caption = paste0(
              "Vertical dashed line denotes inital year of \n",
              "Mass Drug Administration (MDA)"
            )
          )
        }
      }
  } else {
    gg_admin <-
      ggplot(data = model_data, aes(x = year, y = mean, ymin = lower, ymax = upper)) +
      geom_ribbon(alpha = 0.3) +
      geom_line() +
      theme_bw(base_size = 16) +
      coord_cartesian(ylim = val_range) +
      facet_wrap(~plot_name) +
      geom_vline(data = mda, aes(xintercept = mda_start), linetype = "dotted") +
      theme(
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0.5),
        plot.title = element_text(size = title_plot_size, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        y = ind_title,
        x = "Year"
      ) + {
        if (admin == 0) labs(title = paste0(ind_title, " by Country"))
      } + {
        if (admin == 1) {
          labs(
            title = paste0(ind_title, " by First-level Administrative Unit"),
            caption = paste0(
              "Time series depict first-level administrative units except for NATIONAL,\n",
              "which shows the time series for the entire country"
            )
          )
        }
      } + {
        if (admin == 2) {
          labs(
            title = paste0(ind_title, " by Second-level Administrative Unit"),
            caption = paste0(
              "Vertical dashed line denotes inital year of \n",
              "Mass Drug Administration (MDA)"
            )
          )
        }
      }
  }

  if (plot_stkrs) {
    stacker_colors <- data.table("stacker" = c("lasso", "gbm", "gam"), "color" = c("orange", "green", "blue"))
    model_data_df <- as.data.frame(model_data)
    for (s in stkrs) {
      gg_admin <- gg_admin + geom_line(aes_string(y = model_data_df[, s], color = paste0("\"", stacker_colors[stacker == s, "color"], "\"")))
    }
    gg_admin <- gg_admin + scale_color_discrete(name = "Stackers", labels = stkrs)
  }
  return(gg_admin)
}

# Plot time series trends with aggregated input_data
# Maintain consistency between plots in shapes and color graph attributes
time_series_data <- function(model_data,
                             input_data,
                             admin = 0,
                             val_range = c(0, 1),
                             title_plot_size = 30,
                             ind_title = "",
                             mda = NULL,
                             diag_vals = diagnostic_vals,
                             collect_method_vals = data_collect_method_vals,
                             plot_stkrs = plot_stacker_trends,
                             stkrs = stackers,
                             s_fixed_effects = stacked_fixed_effects) {

  # Color palette
  carto_discrete <- c(
    "#7F3C8D", "#11A579", "#F2B701", "#E73F74",
    "#3969AC", "#80BA5A", "#E68310", "#008695",
    "#CF1C90", "#f97b72", "#4b4b8f", "#A5AA99"
  )

  carto_discrete <- carto_discrete[1:length(collect_method_vals)]
  names(carto_discrete) <- collect_method_vals

  if (length(diag_vals) == 2) {
    shapes <- c(22, 21)
  }
  if (length(diag_vals == 3)) {
    shapes <- c(22, 21, 25)
  }
  if (length(diag_vals == 4)) {
    shapes <- c(22, 21, 25, 23)
  }

  names(shapes) <- diag_vals

  if (is.null(mda) == T) {
    gg_admin <-
      ggplot() +
      geom_ribbon(data = model_data, aes(x = year, ymin = lower, ymax = upper), alpha = 0.3) +
      geom_line(data = model_data, aes(x = year, y = mean)) +
      theme_bw(base_size = 16) +
      geom_point(
        data = input_data,
        aes(x = year, y = outcome, size = N, shape = diagnostic, fill = data_collect_method), alpha = 0.7
      ) +
      scale_fill_manual("Stage", values = carto_discrete, drop = F) +
      scale_shape_manual("Diagnostic", values = shapes) +
      facet_wrap(~plot_name) +
      coord_cartesian(ylim = val_range) +
      scale_size(range = c(1, 7)) +
      theme(
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0.5),
        plot.title = element_text(size = title_plot_size, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.justification = "top"
      ) +
      guides(
        fill = guide_legend(order = 1, override.aes = list(size = 5, shape = 21)),
        shape = guide_legend(order = 2, override.aes = list(size = 5))
      ) +
      labs(
        y = ind_title,
        x = "Year"
      ) + {
        if (admin == 0) labs(title = paste0(ind_title, " by Country"))
      } + {
        if (admin == 1) {
          labs(
            title = paste0(ind_title, " by First-level Administrative Unit"),
            caption = paste0(
              "Time series depict first-level administrative units except for NATIONAL,\n",
              "which shows the time series for the entire country"
            )
          )
        }
      } + {
        if (admin == 2) {
          labs(
            title = paste0(ind_title, " by Second-level Administrative Unit"),
            caption = paste0(
              "Vertical dashed line denotes inital year of \n",
              "Mass Drug Administration (MDA)"
            )
          )
        }
      }
  } else {
    gg_admin <-
      ggplot() +
      geom_ribbon(data = model_data, aes(x = year, ymin = lower, ymax = upper), alpha = 0.3) +
      geom_line(data = model_data, aes(x = year, y = mean)) +
      theme_bw(base_size = 16) +
      geom_point(
        data = input_data,
        aes(x = year, y = outcome, size = N, shape = diagnostic, fill = data_collect_method), alpha = 0.7
      ) +
      scale_fill_manual("Stage", values = carto_discrete, drop = F) +
      scale_shape_manual("Diagnostic", values = shapes) +
      facet_wrap(~plot_name) +
      geom_vline(data = mda, aes(xintercept = mda_start), linetype = "dotted") +
      coord_cartesian(ylim = val_range) +
      scale_size(range = c(1, 7)) +
      theme(
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0.5),
        plot.title = element_text(size = title_plot_size, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.justification = "top"
      ) +
      guides(
        fill = guide_legend(order = 1, override.aes = list(size = 5, shape = 21)),
        shape = guide_legend(order = 2, override.aes = list(size = 5))
      ) +
      labs(
        y = ind_title,
        x = "Year"
      ) + {
        if (admin == 0) labs(title = paste0(ind_title, " by Country"))
      } + {
        if (admin == 1) {
          labs(
            title = paste0(ind_title, " by First-level Administrative Unit"),
            caption = paste0(
              "Time series depict first-level administrative units except for NATIONAL,\n",
              "which shows the time series for the entire country"
            )
          )
        }
      } + {
        if (admin == 2) {
          labs(
            title = paste0(ind_title, " by Second-level Administrative Unit"),
            caption = paste0(
              "Vertical dashed line denotes inital year of \n",
              "Mass Drug Administration (MDA)"
            )
          )
        }
      }
  }

  if (plot_stkrs) {
    stacker_colors <- data.table("stacker" = c("lasso", "gbm", "gam"), "color" = c("orange", "green", "blue"))
    model_data_df <- as.data.frame(model_data)
    for (s in stkrs) {
      gg_admin <- gg_admin + geom_line(data = model_data_df, aes_string(y = model_data_df[, s], x = model_data_df$year, color = paste0("\"", stacker_colors[stacker == s, "color"], "\"")))
    }
    gg_admin <- gg_admin + scale_color_discrete(name = "Stackers", labels = stkrs)
  }

  return(gg_admin)
}

# Plot time series plot with multiple indicators or multiple model runs, with run column specified
time_series_multiple <- function(model_data,
                                 admin = 0,
                                 val_range = c(0, 1),
                                 title_plot_size = 30,
                                 ind_title = "",
                                 mda = NULL) {

  # Color palette
  carto_discrete <- c(
    "#7F3C8D", "#11A579", "#F2B701", "#E73F74",
    "#3969AC", "#80BA5A", "#E68310", "#008695",
    "#CF1C90", "#f97b72", "#4b4b8f", "#A5AA99"
  )
  if (is.null(mda) == T) {
    gg_admin <-
      ggplot(data = model_data, aes(x = year, y = mean, ymin = lower, ymax = upper)) +
      geom_ribbon(aes(fill = run), alpha = 0.1) +
      geom_line(aes(color = run, linetype = run), alpha = 0.9) +
      theme_bw(base_size = 16) +
      scale_color_manual("", values = carto_discrete) +
      scale_fill_manual("", values = carto_discrete) +
      scale_linetype_discrete("") +
      coord_cartesian(ylim = val_range) +
      facet_wrap(~plot_name) +
      theme(
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0.5),
        plot.title = element_text(size = title_plot_size, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.justification = "top"
      ) +
      labs(
        y = ind_title,
        x = "Year"
      ) +
      guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 1),
        linetype = guide_legend(order = 1)
      ) + {
        if (admin == 0) labs(title = paste0(ind_title, " by Country"))
      } + {
        if (admin == 1) {
          labs(
            title = paste0(ind_title, " by First-level Administrative Unit"),
            caption = paste0(
              "Time series depict first-level administrative units except for NATIONAL,\n",
              "which shows the time series for the entire country"
            )
          )
        }
      } + {
        if (admin == 2) {
          labs(
            title = paste0(ind_title, " by Second-level Administrative Unit"),
            caption = paste0(
              "Vertical dashed line denotes inital year of \n",
              "Mass Drug Administration (MDA)"
            )
          )
        }
      }
  } else {
    gg_admin <-
      ggplot(data = model_data, aes(x = year, y = mean, ymin = lower, ymax = upper)) +
      geom_ribbon(aes(fill = run), alpha = 0.1) +
      geom_line(aes(color = run, linetype = run), alpha = 0.9) +
      theme_bw(base_size = 16) +
      scale_color_manual("", values = carto_discrete) +
      scale_fill_manual("", values = carto_discrete) +
      scale_linetype_discrete("") +
      coord_cartesian(ylim = val_range) +
      facet_wrap(~plot_name) +
      geom_vline(data = mda, aes(xintercept = mda_start), linetype = "dotted") +
      theme(
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0.5),
        plot.title = element_text(size = title_plot_size, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.justification = "top"
      ) +
      labs(
        y = ind_title,
        x = "Year"
      ) +
      guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 1),
        linetype = guide_legend(order = 1)
      ) + {
        if (admin == 0) labs(title = paste0(ind_title, " by Country"))
      } + {
        if (admin == 1) {
          labs(
            title = paste0(ind_title, " by First-level Administrative Unit"),
            caption = paste0(
              "Time series depict first-level administrative units except for NATIONAL,\n",
              "which shows the time series for the entire country"
            )
          )
        }
      } + {
        if (admin == 2) {
          labs(
            title = paste0(ind_title, " by Second-level Administrative Unit"),
            caption = paste0(
              "Vertical dashed line denotes inital year of \n",
              "Mass Drug Administration (MDA)"
            )
          )
        }
      }
  }

  return(gg_admin)
}

# Plot time series with multiple model runs, including data
# Maintain consistency between plots in shapes and color graph attributes
time_series_multiple_data <- function(model_data,
                                      input_data,
                                      admin = 0,
                                      val_range = c(0, 1),
                                      title_plot_size = 30,
                                      ind_title,
                                      mda = NULL) {

  # Color palette
  carto_discrete <- c(
    "#7F3C8D", "#11A579", "#F2B701", "#E73F74",
    "#3969AC", "#80BA5A", "#E68310", "#008695",
    "#CF1C90", "#f97b72", "#4b4b8f", "#A5AA99"
  )

  carto_discrete <- carto_discrete[1:length(collect_method_vals)]
  names(carto_discrete) <- collect_method_vals

  if (length(diag_vals) == 2) {
    shapes <- c(22, 21)
  } else if (length(diag_vals == 3)) {
    shapes <- c(22, 21, 25)
  }
  names(shapes) <- diag_vals

  if (is.null(mda) == T) {
    gg_admin <-
      ggplot() + {
        if (length(unique(model_data$run)) <= 3) {
          geom_ribbon(
            data = model_data,
            aes(x = year, ymin = lower, ymax = upper, fill = run), alpha = 0.1
          )
        }
      } + {
        if (length(unique(model_data$run)) > 3) {
          geom_ribbon(
            data = model_data,
            aes(x = year, ymin = lower, ymax = upper, fill = run), alpha = 0.03
          )
        }
      } +
      geom_line(data = model_data, aes(x = year, y = mean, color = run, linetype = run), alpha = 0.8) +
      geom_point(
        data = input_data,
        aes(x = year, y = outcome, size = N, shape = diagnostic), fill = "grey", alpha = 0.7
      ) +
      theme_bw(base_size = 16) +
      scale_color_manual("", values = carto_discrete) +
      scale_fill_manual("", values = carto_discrete) +
      scale_shape_manual("Diagnostic", values = shapes) +
      scale_linetype_manual("", values = rep(c("solid", "dotdash", "dashed"), 5)) +
      coord_cartesian(ylim = val_range) +
      facet_wrap(~plot_name) +
      theme(
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0.5),
        plot.title = element_text(size = title_plot_size, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.justification = "top"
      ) +
      labs(
        y = ind_title,
        x = "Year"
      ) +
      guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 1),
        linetype = guide_legend(order = 1),
        shape = guide_legend(order = 2, override.aes = list(size = 5, fill = "grey"))
      ) + {
        if (admin == 0) labs(title = paste0(ind_title, " by Country"))
      } + {
        if (admin == 1) {
          labs(
            title = paste0(ind_title, " by First-level Administrative Unit"),
            caption = paste0(
              "Time series depict first-level administrative units except for NATIONAL,\n",
              "which shows the time series for the entire country"
            )
          )
        }
      } + {
        if (admin == 2) {
          labs(
            title = paste0(ind_title, " by Second-level Administrative Unit"),
            caption = paste0(
              "Vertical dashed line denotes inital year of \n",
              "Mass Drug Administration (MDA)"
            )
          )
        }
      }
  } else {
    gg_admin <-
      ggplot() + {
        if (length(unique(model_data$run)) <= 3) {
          geom_ribbon(
            data = model_data,
            aes(x = year, ymin = lower, ymax = upper, fill = run), alpha = 0.1
          )
        }
      } + {
        if (length(unique(model_data$run)) > 3) {
          geom_ribbon(
            data = model_data,
            aes(x = year, ymin = lower, ymax = upper, fill = run), alpha = 0.03
          )
        }
      } +
      geom_line(data = model_data, aes(x = year, y = mean, color = run, linetype = run), alpha = 0.8) +
      geom_point(
        data = input_data,
        aes(x = year, y = outcome, size = N, shape = diagnostic), fill = "grey", alpha = 0.7
      ) +
      theme_bw(base_size = 16) +
      scale_color_manual("", values = carto_discrete) +
      scale_fill_manual("", values = carto_discrete) +
      scale_shape_manual("Diagnostic", values = shapes) +
      scale_linetype_manual("", values = rep(c("solid", "dotdash", "dashed"), 5)) +
      coord_cartesian(ylim = val_range) +
      facet_wrap(~plot_name) +
      geom_vline(data = mda, aes(xintercept = mda_start), linetype = "dotted") +
      theme(
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0.5),
        plot.title = element_text(size = title_plot_size, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.justification = "top"
      ) +
      labs(
        y = ind_title,
        x = "Year"
      ) +
      guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 1),
        linetype = guide_legend(order = 1),
        shape = guide_legend(order = 2, override.aes = list(size = 5, fill = "grey"))
      ) + {
        if (admin == 0) labs(title = paste0(ind_title, " by Country"))
      } + {
        if (admin == 1) {
          labs(
            title = paste0(ind_title, " by First-level Administrative Unit"),
            caption = paste0(
              "Time series depict first-level administrative units except for NATIONAL,\n",
              "which shows the time series for the entire country"
            )
          )
        }
      } + {
        if (admin == 2) {
          labs(
            title = paste0(ind_title, " by Second-level Administrative Unit"),
            caption = paste0(
              "Vertical dashed line denotes inital year of \n",
              "Mass Drug Administration (MDA)"
            )
          )
        }
      }
  }
  return(gg_admin)
}
