#####################################################################
# Summarize gam models and output model summaries
#####################################################################

# (1) Setup -------------------------------------------------------------
rm(list = ls())

run_date <- '2019_09_12_17_11_54'

library(rmarkdown)

# (2) Render HTML in diagnostics plots folder ----------------------------

code_dir <- '<<<< FILEPATH REDACTED >>>>'
report_filename <- "vet_gam.Rmd"
report_filename <- file.path(code_dir, report_filename)

output_dir <- '<<<< FILEPATH REDACTED >>>>'
print(output_dir)
if (!dir.exists(output_dir)) dir.create(output_dir)
output <- file.path("..",output_dir)

render(report_filename, output_dir = output_dir, params = list(output_dir = output))
