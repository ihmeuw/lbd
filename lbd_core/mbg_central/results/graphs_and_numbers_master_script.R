##########
# Master script to run markdown, and any other plots and numbers you want regularly pulled from your model runs.
##########

# Setting the root, loading in libraries, and functions:
setwd("<<<< FILEPATH REDACTED >>>>")
for(function_script in list.files(getwd(),pattern="*_functions.R")){message(function_script);source(function_script)};message("Central Functions Loaded.")
load_libs(c('data.table','rmarkdown'))

indicator_group<-commandArgs()[4]; message(indicator_group)
results_pull_time<-commandArgs()[5]; message(results_pull_time)

results<-fread("<<<< FILEPATH REDACTED >>>>")


# Running the Markdown
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rmarkdown::render("<<<< FILEPATH REDACTED >>>>",
                  params= list(indicator_group= indicator_group,
                               results_pull_time=results_pull_time),
                  envir = new.env(),
                  output_file="<<<< FILEPATH REDACTED >>>>")


# Discover the correlation coefficient between the raw and estimated data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Getting the correlation between the raw data and the estimates (unraked):")

corrs<-list()
for (i in seq(1,nrow(results[raked==F]))){
  message(i)
  corrs[[results[raked==F]$indicator_longname[i]]] <- get_weighted_correlation_raw_estimated(indicator_group=results[raked==F]$indicator_group[i],
                                                                                             indicator=results[raked==F]$indicator_name[i],
                                                                                             rundate=results[raked==F]$date_string[i])
}
corrs<-as.data.table(corrs)
write.csv(corrs, file="<<<< FILEPATH REDACTED >>>>")

