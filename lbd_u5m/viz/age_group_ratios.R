rm(list=ls())

# Set J directory path
if(Sys.info()[1]=='Windows') j_head <- '<<<< FILEPATH REDACTED >>>>' else j_head <- '<<<< FILEPATH REDACTED >>>>'
# Imports
library(data.table)
library(foreign)
library(ggplot2)
library(ggrepel)
library(ggtern, lib.loc='<<<< FILEPATH REDACTED >>>>')
library(manipulate, lib.loc='<<<< FILEPATH REDACTED >>>>')
library(RColorBrewer)
# Source LBD core functions
source(paste0(j_head,"<<<< FILEPATH REDACTED >>>>"))


## DEFINE INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_date <- "<<<< REDACTED >>>>"
adm_level <- 1
out_dir <- paste0('<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()) )
dir.create(out_dir, showWarnings = FALSE)


## FUNCTION DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function to prep administrative data for ternary plotting ~~~~~~~~~~~~~~~~~~~
prep_ternary_data <- function(run_date, adm_level=2){
  # List to store individual data types
  data_types <- list(
    'neonatal'           = NA,
    'infant'             = NA,
    'under5'             = NA,
    'under5_deathcounts' = NA
  )
  for (ag in names(data_types)){
    if(ag=='under5_deathcounts'){
      # Handle death counts a bit differently
      ag_lab <- 'under5'
      dc <- '_deathcounts'
    } else {
      ag_lab <- ag
      dc <- ''
    }
    # Define input filepath
    fp <- paste0('<<<< FILEPATH REDACTED >>>>',ag_lab,'<<<< FILEPATH REDACTED >>>>',run_date,
                 '<<<< FILEPATH REDACTED >>>>',ag_lab,dc,'<<<< FILEPATH REDACTED >>>>',adm_level,'<<<< FILEPATH REDACTED >>>>')
    # Read file
    dt <- fread(fp)
    setnames(dt, paste0("ADM",adm_level,"_CODE"), "adm_code")
    dt <- dt[,.(adm_code,year,mean)]
    setnames(dt,'mean','q')
    dt[, ag := ag]
    data_types[[ag]] <- dt
  }
  # Combine into a single data.table for all data types
  data_full <- rbindlist(data_types)
  # Reshape wide by age group
  data_wide <- dcast(data_full, adm_code + year ~ ag,
                     value.var = 'q', fun.aggregate = sum)
  # Rename death counts
  setnames(data_wide, 'under5_deathcounts', 'd')
  # Get RELATIVE probaility of dying in a given age group
  data_wide[, prob_nn  := neonatal / under5 ]
  # Keep only needed columns
  prepped <- data_wide[, .(adm_code, year, prob_nn, under5, infant, neonatal, d)]

  # Fix for now because not much 2017 data is coming in
  in_data <- prepped[year < 2017,]

  # Merge on admin metadata
  hierarchy <- read.dbf(get_admin_shapefile(adm_level, suffix=".dbf"))
  hierarchy <- as.data.table(hierarchy)
  setnames(hierarchy, "ADM0_NAME", "country")
  setnames(hierarchy, paste0("ADM",adm_level,"_CODE"), "adm_code")
  setnames(hierarchy, paste0("ADM",adm_level,"_NAME"), "adm_name")
  
  hierarchy <- hierarchy[, .(adm_code, adm_name, country)]
  
  with_country <- merge(
    x = prepped,
    y = hierarchy,
    by = c('adm_code')
  )
  with_country[, level := adm_level]
  # Return data with country metadata
  return(with_country)
}


## Function to find particular exemplars for plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~
find_exemplars <- function(in_data){
  ## Create a test dataset that contains the first and last years of data
  ##  as well as ARCs to better identify exemplars
  measure_cols <- c('prob_nn','under5')
  start_year <- min(in_data[,year])
  end_year   <- max(in_data[,year])
  in_data_sub <- copy(in_data)[year %in% c(start_year, end_year),
                               c('adm_code', 'year', measure_cols), with=F]
  ex_dt <- dcast(
    in_data_sub,
    formula   = adm_code ~ year,
    value.var = measure_cols)
  # Calculate ARC for each measure variable
  for(measure_col in measure_cols){
    arc_col   <- paste0(measure_col,'_arc')
    start_col <- paste0(measure_col,'_',start_year)
    end_col   <- paste0(measure_col,'_',end_year)
    ex_dt[, (arc_col) := log(get(end_col)/get(start_col)) / (end_year-start_year)]
  }
  
  # Helper function to concatenate text together for the following functions
  safe_concat <- function(text1, text2){
    if(is.na(text1)){
      return(text2)
    } else {
      return(paste(text1, text2, sep=';\n'))
    } 
  }

  ## SEEK EXAMPLARS:
  ex_dt[, exemplar := character(0)]
  # - One with the biggest swing in nn/5q0 ratio over the time period
  ex_dt[, abs_prob_nn_arc := abs(prob_nn_arc)]
  ex_dt[ abs_prob_nn_arc == max(abs_prob_nn_arc, na.rm=TRUE),
        exemplar := safe_concat(exemplar, '1. Biggest swing')]
  # - one with the smallest swing in nn/5q0 ratio over the time period
  ex_dt[ abs_prob_nn_arc == min(abs_prob_nn_arc, na.rm=TRUE),
        exemplar := safe_concat(exemplar, '2. Smallest swing')]  
  # - One with a 'backwards' swing. if this happened anywhere
  ex_dt[ prob_nn_arc == min(prob_nn_arc, na.rm=TRUE),
        exemplar := safe_concat(exemplar, '3. Backwards swing')]
  # - One with the most extreme (high) ratio
  end_ratio <- paste0('prob_nn_',end_year)
  ex_dt[ get(end_ratio) == max(get(end_ratio), na.rm=TRUE),
        exemplar := safe_concat(exemplar, paste0('4. Highest NN/U5\nratio in ',end_year))]
  # - One with the most extreme (low) ratio
  ex_dt[ get(end_ratio) == min(get(end_ratio), na.rm=TRUE),
         exemplar := safe_concat(exemplar, paste0('5. Lowest NN/U5\nratio in ',end_year))]
  # - One with high 5q0 but also high NN proportion
  end_5q0 <- paste0('under5_',  end_year)
  end_nn  <- paste0('prob_nn_', end_year)
  ex_6_nn_prop <- ex_dt[ get(end_5q0) >= quantile(get(end_5q0), .9, na.rm=TRUE),
                         max(get(end_nn), na.rm=TRUE)]
  ex_dt[ (get(end_5q0) >= quantile(get(end_5q0), .9, na.rm=TRUE)) & 
         (get(end_nn)  >= ex_6_nn_prop), 
        exemplar := safe_concat(exemplar, '6. High 5q0 and\nhigh NN proportion')]
  # - One with low 5q0 and with low NN proportion
  ex_7_nn_prop <- ex_dt[ get(end_5q0)<=quantile(get(end_5q0), .1, na.rm=TRUE),
                         min(get(end_nn), na.rm=TRUE)]
  ex_dt[ (get(end_5q0) <= quantile(get(end_5q0), .1, na.rm=TRUE)) & 
         (get(end_nn)  <= ex_7_nn_prop), 
         exemplar := safe_concat(exemplar, '7. High 5q0 and\nhigh NN proportion')]
  
  ## Merge back onto the main dataset, keeping only the exemplar locations
  ex_dt_sub <- ex_dt[ !is.na(exemplar), .(adm_code, exemplar)]
  exemplars_only <- merge(
    x  = in_data,
    y  = ex_dt_sub,
    by = c('adm_code')
  )
  # Return the dataset of exemplars only
  return(exemplars_only)
}


## Function to create ternary line plots for a given country ~~~~~~~~~~~~~~~~~~~
plot_ternary <- function(data, country_name){
  # Restrict to the given country
  data_sub <- data[ country == country_name, ]
  # PLOT BY U5M
  fig <- ggtern(data=data_sub,
                aes(x=prob_nn, y=prob_inf, z=prob_u5, group=ADM2_CODE,
                    color=under5)
                ) +
    geom_line(alpha=.4) +
    scale_color_gradientn(
      colors=rev(brewer.pal(name="RdYlBu",9)),
      guide = guide_colorbar(
        title.position='top',
        barwidth=unit(18,'line')
      )
    ) +
    labs(
      title = 'Evolution of mortality distribution over time by administrative unit',
      subtitle = paste0('Administrative units: Admin',2,' in ',country_name),
      x = 'Prob. of dying\n(Neonatal)',
      y = 'Prob. of dying\n(1-12 months)',
      z = 'Prob. of dying\n(1-4 years)',
      color = 'Probability of dying under age 5 in location-year'
    ) +
    tern_limit(.8, .8, .8) +
    theme_bw() +
    theme(legend.position='bottom')
  # PLOT BY YEAR
  fig2 <- ggtern(data=data_sub,
                aes(x=prob_nn, y=prob_inf, z=prob_u5, group=ADM2_CODE,
                    color=year)
    ) +
    geom_line(alpha=.2) +
    scale_color_gradientn(
      colors=brewer.pal(name="Spectral",11),
      guide = guide_colorbar(
        title.position='top',
        barwidth=unit(18,'line')
      )
    ) +
    labs(
      title = 'Evolution of mortality distribution over time by administrative unit',
      subtitle = paste0('Administrative units: Admin',2,' in ',country_name),
      x = 'Prob. of dying\n(Neonatal)',
      y = 'Prob. of dying\n(1-12 months)',
      z = 'Prob. of dying\n(1-4 years)',
      color = 'Year'
    ) +
    tern_limit(.8, .8, .8) +
    theme_bw() +
    theme(legend.position='bottom')

  return(list(fig1=fig, fig2=fig2))
}


################################################################################
## PROGRAM EXECUTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################################

## Prep input data
in_data <- prep_ternary_data(
  run_date  = run_date,
  adm_level = adm_level
)
# Get probability of death from 1 year to 4 years
in_data[, q_child := 1 - (1-under5)/(1-infant)]

in_data <- in_data[ year < 2017,  ]
in_data <- in_data[order(country, adm_name, year)]


## Prepare a plot showing change in overall 5q0 versus NN/5q0 probability
in_data_final <- in_data[year == max(year), ]


pdf(
  paste0(out_dir,'<<<< FILEPATH REDACTED >>>>'),
  width=10,
  height=10
)

for(loc in unique(in_data[,country])){
  fig_1 <- ggplot() +
    geom_line(data=in_data[!(country%in%loc)], 
              aes(x=under5, y=prob_nn, group=adm_name),lwd=.1, alpha=.2, color='#888888') +
    geom_point(data=in_data_final[!(country%in%loc)],
               aes(x=under5, y=prob_nn, size=d), lwd=.1, alpha=.2, color='#888888') +
    geom_line(data=in_data[(country%in%loc)], 
              aes(x=under5, y=prob_nn, group=adm_name, color=adm_name)) +
    geom_point(data=in_data_final[(country%in%loc)],
               aes(x=under5, y=prob_nn, size=d, color=adm_name)) +
    geom_text_repel(
      data=in_data_final[(country%in%loc)],
      aes(x=under5, y=prob_nn, label=adm_name, color=adm_name),
          box.padding = unit(.3, 'lines')
    ) +
    labs(
      title = 'Changes in NN mortality probability versus 5q0 over time',
      subtitle = paste0("Country selected:", loc),
      x = "Probability of death before age 5 (5q0)",
      y = "Neonatal mortality probability / Under-5 mortality probability",
      size = "Number\nof deaths\n(2016)"
    ) +
    xlim(0,.3) + ylim(.15,.80) +
    theme_bw() + theme(legend.position='none')
  print(fig_1)
}

dev.off()



## Prepare a plot that shows the relationship between 4q1 and NNq0 over time

fig_2 <- ggplot(data=in_data, aes(x=q_child, y=neonatal, group=adm_name, color=country)) +
  geom_line() + theme(legend.position='none')
fig_2
  

## PREPARE A PLOT SHOWING NN versus 5q0 in 2010 versus in 2017

for(out_year in c(2010, 2016)){
  in_data_thisyear <- in_data[year == out_year, ]
  pdf(
    paste0(out_dir,'<<<< FILEPATH REDACTED >>>>',adm_level,'_',out_year,'.pdf'),
    width=10,
    height=10
  )
  fig_1 <- ggplot() +
    geom_point(data=in_data_thisyear,
               aes(x=under5, y=prob_nn, size=d, color=adm_name)) +
    labs(
      title = paste0('NN mortality probability versus 5q0 in ',out_year),
      x = "Probability of death before age 5 (5q0)",
      y = "Neonatal mortality probability / Under-5 mortality probability",
      size = paste0("Number\nof deaths\n(",out_year,")")
    ) +
    xlim(0,.3) + ylim(.15,.80) +
    theme_bw() + theme(legend.position='none')
  print(fig_1)

  dev.off()
}