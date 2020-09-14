###############################################################
# Make mapping inputs for deaths averted figures
# LRI Stg 2
###############################################################
prep_deaths_averted_figures <- function(map_date,
                                        run_date,
                                        end_year){
  #(1) Setup #############################################################
  map_outdir <- '<<<< FILEPATH REDACTED >>>>'
  
  #(2) Data prep ##############################################################
  #read in deaths averted data and aggregated pop values
  results <- fread('<<<< FILEPATH REDACTED >>>>')
  load('<<<< FILEPATH REDACTED >>>>')
  pop <- admin_2 %>%
    select(year, ADM2_CODE, pop) %>%
    filter(year == end_year)
  rm(admin_0, admin_1, admin_2)
  
  #add on a column for total_paf_2017, total_paf_2000, total cgf paf 2000, total cgf paf 2017
  results <- mutate(results, 
                    total_paf_2000 = 1 - (1 - stunting_paf00) * (1 - wasting_paf00) * (1 - underweight_paf00) * (1 - tap_paf00),
                    total_paf_2017 = 1 - (1 - stunting_paf17) * (1 - wasting_paf17) * (1 - underweight_paf17) * (1 - tap_paf17),
                    total_cgf_paf_2000 = 1 - (1 - stunting_paf00) * (1 - wasting_paf00) * (1 - underweight_paf00),
                    total_cgf_paf_2017 = 1 - (1 - stunting_paf17) * (1 - wasting_paf17) * (1 - underweight_paf17))
  
  #calculate total deaths averted & merge on pop column
  results <- mutate(results,
                    total_deaths_averted = lri_deaths_17 *((1-total_paf_2017)/(1-total_paf_2000)*total_paf_2000-total_paf_2017),
                    cgf_deaths_averted = lri_deaths_17 *((1-total_cgf_paf_2017)/(1-total_cgf_paf_2000)*total_cgf_paf_2000-total_cgf_paf_2017)) %>%
    merge(pop, by = 'ADM2_CODE')
  
  #add rate columns for deaths averted due to cgf and ors
  results <- mutate(results,
                    cgf_deaths_averted_rate = cgf_deaths_averted / pop,
                    tap_deaths_averted_rate = tap_deaths_averted / pop,
                    total_deaths_averted_rate = total_deaths_averted / pop)
  
  #make column for dominant cause (rate and count)
  # dominant rate total
  results <- as.data.table(results)
  results[, dominant_rate_total := ifelse(cgf_deaths_averted_rate / total_deaths_averted_rate > 0.5 & tap_deaths_averted_rate / total_deaths_averted_rate < 0.5, 'cgf', 'none')]
  results[, dominant_rate_total :=  ifelse(cgf_deaths_averted_rate / total_deaths_averted_rate < 0.5 & tap_deaths_averted_rate / total_deaths_averted_rate >  0.5, 'tap', dominant_rate_total)]
  
  # dominant count total
  results[, dominant_count_total := ifelse(cgf_deaths_averted / total_deaths_averted > 0.5 & tap_deaths_averted / total_deaths_averted < 0.5, 'cgf', 'none')]
  results[, dominant_count_total :=  ifelse(cgf_deaths_averted / total_deaths_averted < 0.5 & tap_deaths_averted / total_deaths_averted >  0.5, 'tap', dominant_count_total)]
  
  #(3) Prep csvs ##############################################################
  #need year, adm2 code, value
  
  #Fig 6: wash + cgf deaths averted
  #deaths averted total (cgf + tap) rate
  da_total_rate <- select(results, total_deaths_averted_rate, ADM2_CODE, year) %>%
    rename(value = total_deaths_averted_rate)
  
  #deaths averted total (cgf + tap) count
  da_total_count <- select(results, total_deaths_averted, ADM2_CODE, year) %>%
    rename(value = total_deaths_averted)
  
  #deaths averted total where cgf dominant/wash dominant/none dominant rate
  da_total_rate_cgf <- select(results, total_deaths_averted_rate, ADM2_CODE, dominant_rate_total, year) %>%
    filter(dominant_rate_total == 'cgf') %>%
    rename(value = total_deaths_averted_rate)
  
  da_total_rate_tap <- select(results, total_deaths_averted_rate, ADM2_CODE, dominant_rate_total, year) %>%
    filter(dominant_rate_total == 'tap') %>%
    rename(value = total_deaths_averted_rate)
  
  da_total_rate_none <- select(results, total_deaths_averted_rate, ADM2_CODE, dominant_rate_total, year) %>%
    filter(dominant_rate_total == 'none') %>%
    rename(value = total_deaths_averted_rate)
  
  #deaths averted total where cgf dominant/wash dominant/none dominant count
  da_total_count_cgf <- select(results, total_deaths_averted, ADM2_CODE, dominant_count_total, year) %>%
    filter(dominant_count_total == 'cgf') %>%
    rename(value = total_deaths_averted)
  
  da_total_count_tap <- select(results, total_deaths_averted, ADM2_CODE, dominant_count_total, year) %>%
    filter(dominant_count_total == 'tap') %>%
    rename(value = total_deaths_averted)
  
  da_total_count_none <- select(results, total_deaths_averted, ADM2_CODE, dominant_count_total, year) %>%
    filter(dominant_count_total == 'none') %>%
    rename(value = total_deaths_averted)
  
  #(4) Write out results ##############################################################
  write.csv(da_total_rate, '<<<< FILEPATH REDACTED >>>>')
  write.csv(da_total_count, '<<<< FILEPATH REDACTED >>>>')
  
  write.csv(da_total_rate_cgf, '<<<< FILEPATH REDACTED >>>>')
  write.csv(da_total_rate_tap, '<<<< FILEPATH REDACTED >>>>')
  write.csv(da_total_rate_none, '<<<< FILEPATH REDACTED >>>>')
  
  write.csv(da_total_count_cgf, '<<<< FILEPATH REDACTED >>>>')
  write.csv(da_total_count_tap, '<<<< FILEPATH REDACTED >>>>')
  write.csv(da_total_count_none, '<<<< FILEPATH REDACTED >>>>')
}
