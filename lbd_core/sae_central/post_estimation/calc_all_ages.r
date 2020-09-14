####################################################################################################
## Description: Define function for adding crude and age-standardized rates to age-specific rates
##
## Inputs:      draws - a data.table of age-specific rate draws plus the corresponding
##                population counts to be used for generating crude rates.
##              std_wt - a data.table of the weights for generating age-standardized rates.
##
## Outputs:     data.table of age-specific PLUS crude and age-standardized rates.
####################################################################################################

calc_all_ages <- function(draws, std_wt) {

  if (!setequal(std_wt$age, unique(draws$age))) stop("age mismatch")

  # calculate weights for generating crude rates based on observed populations
  draws[pop < 1e-5, pop := 0] # change effective 0s to 0s to avoid some annoying edge cases where one age gets all the weight even if it also has effectively no population
  draws[, crude_wt := pop/sum(pop), by='level,area,year,sex,sim']

  # for areas with no population, calculate weights based on total populations across all areas
  draws[, total := sum(pop), by='level,year,sex,age,sim']
  draws[is.na(crude_wt), crude_wt := total/sum(total), by='level,area,year,sex,sim']

  # collapse to get crude draws
  crude <- draws[, list(value = sum(value * crude_wt)), by='level,area,year,sex,sim']
  crude[, age := 98L]
  draws[, c("pop", "total", "crude_wt") := NULL]

  # merge on standard weights
  draws <- merge(draws, std_wt, by='age')

  # collapse to get age-standardized rates
  std <- draws[, list(value = sum(value*wt)), by='level,area,year,sex,sim']
  std[, age := 99L]
  draws[, wt := NULL]

  # combine age-specific rates with crude and age-standardized rates
  draws <- rbind(draws, crude, std, use.names=T)

  # resort, reorder, and return
  setkeyv(draws, c("level", "area", "year", "sex", "age", "sim"))
  setcolorder(draws, c("level", "area", "year", "sex", "age", "sim", "value"))
  return(draws)
}
