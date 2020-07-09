# Oral rehydration therapies - low and middle-income countries - 2020

This repo contains cleaned code used to generate estimates of oral rehydration solution (ORS) and recommended home fluids (RHF) coverage in children under the age of five with diarrhoea from 'Mapping geographic inequalities in oral rehydration therapy coverage in low-income and middle-income countries, 2000–17' published in _Lancet Global Health_ on July 22, 2020.

The repo contains six directories:
1) covariate_selection - containing code to perform covariate selection by variance inflation factors
2) custom_functions - containing data cleaning and analysis functions
3) gbm_optim - containing code to perform gradient boosted machine (GBM) optimization
4) mbg_central - containing model-based geostatistical code similar to code used across multiple projects at IHME
5) ors - containing code for data scoping, extraction, and cleaning, as well as code to launch model-based geostatistical models for ORS, RHF, and ORT (ORS or RHF)
6) post_estimation - containing code to clean, vet, and analyze modeled estimates

# Other LBD work

To see what work our team has published, please visit the [main IHME LBD team page](http://www.healthdata.org/lbd).

To see other published LBD code repositories, please visit our [LBD code homepage](https://github.com/ihmeuw/lbd).