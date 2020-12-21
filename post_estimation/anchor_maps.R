# ------------------------------------------------------------------------------------------------------------
# make_anchor_maps()
#
# Function that calculate changes in treatment at the admin 0, 1, and 2 draw level and save outputs
#
# Inputs:
# run_date - run date for current model
# year_start - start date for analysis (start of modeling period, pre/post intervention, etc.)
# year_end - end date for analysis (end of modeling period, pre/post intervention, etc.)
# holdout - which holdout we're running this for (0 is the full model)
# country - iso3 for single country to run mapping for
# country_name - full names of the country to run mapping for
# shp - shapefile to use (doesn't need subsetted)
# change_limits - upper and lower limits for the treatment change per year plots
# holdout - which holdout you're running this for
#
# Outputs (saved in /results_maps/ of no_ort output folder):
# - PDFs containing the treatment change and efficacy maps, saved separately by country
# ------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Start function
make_anchor_maps <- function(run_date,
                             year_start,
                             year_end,
                             country,
                             country_name,
                             shp,
                             change_limits = c(-0.025, 0.025),
                             holdout = 0) {
  # -----------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------------
  # Set-up
  
  # set indicator group and indicators
  indicator_group <- 'ort'
  indicators <- c('no_ort', 'rhf_only', 'any_ors')
  
  # define share directories
  share_dirs <- paste0('<<<< FILEPATH REDACTED >>>>')
  names(share_dirs) <- indicators
  maindir <- paste0('<<<< FILEPATH REDACTED >>>>')
  outdir <- paste0('<<<< FILEPATH REDACTED >>>>')
  # -----------------------------------------------------------------------------------
  
  
  # ----------------------------------------------------------------------------------------------
  # Plotting functions
  
  # plotting function 1
  plot_map1 <- function(indat, measure, measure_name, invert = F, axis_limits = c(-0.02, 0.02)) {
    inv <- ifelse(invert, -1, 1)
    ggplot(indat) +
      geom_sf(aes(fill = get(measure)*inv)) +
      scale_fill_scico(palette = 'roma',
                       limits = axis_limits,
                       oob = scales::squish) +
      labs(fill = measure_name) +
      theme_classic() +
      coord_sf(datum = NA)
  }
  
  # plotting function 2
  plot_map2 <- function(indat, measure, measure_name, labs, custom_pal, legend_title = NULL) {
    ggplot(indat) +
      geom_sf(aes(fill = get(measure))) +
      scale_fill_manual(values = custom_pal, 
                        labels = labs,
                        name = legend_title) +
      labs(fill = measure_name) +
      theme_classic() +
      coord_sf(datum = NA)
  }
  # ----------------------------------------------------------------------------------------------
  
  
  # -------------------------------------------------------------------------------------------------------
  # Load and clean data
  
  # load and clean treatment change data
  filepaths <- list()
  for (i in indicators) filepaths[[i]] <- list.files(share_dirs[[i]], pattern = paste0(year_start, '_', year_end, '_', country, '_change_summary_ad2'))
  dt1 <- lapply(paste0(share_dirs, filepaths), fread)
  names(dt1) <- indicators
  
  # clean treatment change data
  for (i in indicators) dt1[[i]][, indicator := i]
  dt1 <- rbindlist(dt1)
  dt1 <- dcast(dt1, ADM0_NAME + ADM2_CODE + ADM2_NAME ~ indicator, value.var = c('mean', 'upper', 'lower'))
  dt1[, mean_rhf_not_replaced := mean_no_ort > 0]
  dt1[, mean_ors_declined := mean_any_ors < 0]
  
  # load and clean ORS efficacy data
  dt2 <- fread(paste0(maindir, 'efficacy_ratio_needed_', year_start, '_', year_end, '_', country, '_ad2.csv'))[, -1]
  dt2 <- merge(dt2, dt1[, c('ADM2_CODE', 'mean_rhf_not_replaced', 'mean_ors_declined')], by = 'ADM2_CODE')
  
  # load and merge draw summary table
  summary <- fread(paste0(maindir, 'draw_summary_', year_start, '_', year_end, '_',  country, '.csv'))
  dt2 <- merge(dt2, summary[agg_level == 'ADM2', c('code', 'ors_declined_prop', 'rhf_replaced_prop', 'efficacy_100_prop')], 
               by.x = 'ADM2_CODE', by.y = 'code')
  
  # change to numeric
  dt1[, mean_rhf_not_replaced := as.numeric(mean_rhf_not_replaced)]
  if (length(unique(dt1$mean_rhf_not_replaced)) == 2) {
    scenario <- 2
  } else if (unique(dt1$mean_rhf_not_replaced) == 1) {
    scenario <- 1
  } else {
    scenario <- 3
  }

  # create an RHF replaced column with uncertainty incorporated
  dt2[mean_rhf_not_replaced == TRUE, repl_uncert := rhf_replaced_prop -1]
  dt2[mean_rhf_not_replaced == FALSE, repl_uncert := rhf_replaced_prop]
  
  # create a mean efficacy column to use for mapping
  dt2[, efficacy_map := mean]
  dt2[efficacy_map >= 10 & efficacy_map < 100, efficacy_map := 10]
  dt2[mean_rhf_not_replaced == FALSE, efficacy_map := 0]
  dt2[mean_ors_declined == TRUE, efficacy_map := 999]
  dt2[, efficacy_map := as.factor(floor(efficacy_map))]
  # -------------------------------------------------------------------------------------------------------
  
  
  # -------------------------------------------------------------------------------------------------------
  # Make the maps
  
  # purple green diverging palette
  palg <- rev(brewer.pal(n = 9, name = 'Greens'))
  palp <- brewer.pal(n = 9, name = 'Purples')
  get_palm <- colorRampPalette(brewer.pal(n = 9, name = 'RdPu')[2:9])
  palm <- get_palm(11)
  
  # pdf to write to
  pdf(paste0(outdir, country_name, '_ors_rhf_efficacy_maps_', year_start, '_', year_end, '.pdf'),
      width = 6, height = 4)
  
  # load and subset shapefile
  shp_cty <- shp[which(shp@data$ADM0_NAME == country_name), ]
  shp_cty@data$ADM2_CODE <- as.numeric(shp_cty@data$ADM2_CODE)
  shp_cty <- st_as_sf(shp_cty)
  
  # map differences in ORS, RHF, and ORT
  map <- left_join(shp_cty, dt1, by = 'ADM2_CODE')
  print(plot_map1(map, 'mean_any_ors', paste0('Rate of change \nin ORS per year\n', year_start, ' to ', year_end),
                  axis_limits = change_limits))
  print(plot_map1(map, 'mean_rhf_only', paste0('Rate of change \nin RHF per year\n', year_start, ' to ', year_end),
                  axis_limits = change_limits))
  print(plot_map1(map, 'mean_no_ort', paste0('Rate of change \nin ORT per year\n', year_start, ' to ', year_end),
                  axis_limits = change_limits, invert = T))
  
  # map places where fewer kids are getting treated in some places
  map$mean_rhf_not_replaced <- as.factor(map$mean_rhf_not_replaced)
  if (scenario == 1) {
    print(plot_map2(map, 'mean_rhf_not_replaced', '', 
                    labs = c('RHF not replaced'),
                    custom_pal = c(palp[6])))
  }
  if (scenario == 2) {
    print(plot_map2(map, 'mean_rhf_not_replaced', '', 
                    labs = c('RHF replaced', 'RHF not replaced'),
                    custom_pal = c(palg[4], palp[6])))
  }
  if (scenario == 3) {
    print(plot_map2(map, 'mean_rhf_not_replaced', '', 
                    labs = c('RHF replaced'),
                    custom_pal = c(palg[4], palp[6])))
  }
  
  # prep to map uncertainty
  get_category <- function(x) {
    if (x == 0) {
      warning('There is a 0 value in your RHF replacement percent of draws analysis. This is suspicious.')
    }
    pos <- ifelse(x > 0, T, F)
    x <- abs(x)
    if (x <= 0.5) c <- 1
    if (x > 0.5 & x <= 0.65) c <- 2
    if (x > 0.65 & x <= 0.75) c <- 3
    if (x > 0.75 & x <= 0.85) c <- 4
    if (x > 0.85 & x <= 0.95) c <- 5
    if (x > 0.95) c <- 6
    if (!pos) c <- c*-1
    return(c)
  }
  dt2[, category := get_category(repl_uncert), 
       by = c(names(dt2)[1:9])]
  dt2$category <- as.factor(dt2$category)
  map_thr <- left_join(shp_cty, dt2, by = 'ADM2_CODE')
  
  # create labels
  get_labels <- function(x) {
    label_list <- list(
      '1' = '< 50',
      '2' = '51-65',
      '3' = '66-75',
      '4' = '76-85',
      '5' = '86-95',
      '6' = '> 95',
      '-1' = '< 50',
      '-2' = '51-65',
      '-3' = '66-75',
      '-4' = '76-85',
      '-5' = '86-95',
      '-6' = '> 95'
    )
    return(label_list[[x]])
  }
  nums <- as.numeric(levels(dt2$category))
  uncert_pal <- c(palp[c(3:7,9)][nums[which(nums < 0)]*-1],
                  rev(palg)[c(3:7,9)][nums[which(nums > 0)]],
                  if(999 %in% nums) 'darkgrey')
  uncert_labels <- sapply(as.character(nums), get_labels)
  
  # make the same map, but have the colors scale from white with uncertainty
  print(plot_map2(map_thr, 'category', '',
                  labs = uncert_labels,
                  custom_pal = uncert_pal,
                  legend_title = 'Percent of draws'))
  
  # create labels and palette for mapping efficacy
  nums <-  as.numeric(as.character(sort(unique(map_thr$efficacy_map))))
  efficacy_labels <- paste0(nums, '-', nums+1)
  efficacy_labels <- efficacy_labels[which(nums < 10)]
  efficacy_pal <- palm[nums[which(nums < 10)]]
  if (0 %in% nums) {
    efficacy_labels <- c('RHF replaced', efficacy_labels[-1])
    efficacy_pal <- c(palg[4], efficacy_pal)
  }
  if (10 %in% nums) {
    efficacy_labels <- c(efficacy_labels, '10-99')
    efficacy_pal <- c(efficacy_pal, palm[10])
  }
  if (100 %in% nums) {
    efficacy_labels <- c(efficacy_labels, '>99')
    efficacy_pal <- c(efficacy_pal, palm[11])
  }
  if (999 %in% nums) {
    efficacy_labels <- c(efficacy_labels, 'ORS declined')
    efficacy_pal <- c(efficacy_pal, 'darkgrey')
  }
  
  # map efficacy ratios
  print(plot_map2(map_thr, 'efficacy_map', '', 
                  labs = efficacy_labels,
                  custom_pal = efficacy_pal,
                  legend_title = 'Relative ORS:RHF \nefficacy needed \nto have saved lives'))
  
  # end plotting
  dev.off()
  # -------------------------------------------------------------------------------------------------------
  
  
  # ------------------------------------------------------
  # End function
  return(paste0('Treatment change and efficacy maps made for ', country, '. Plots saved.'))
  
}
# --------------------------------------------------------
