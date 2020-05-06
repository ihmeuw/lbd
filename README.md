# Diarrhea morbidity and mortality - low and middle-income countries - 2020

This repo contains cleaned code used to generate estimates of diarrhea prevalence, incidence, and mortality in children under the age of five from 'Mapping geographic inequalities in childhood diarrhoeal morbidity and mortality in low-income and middle-income countries, 2000–2017: analysis for the Global Burden of Disease Study 2017'  published in _Lancet_ on May 6, 2020.

The repo contains seven directories:
1) covariate_selection - containing code to perform covariate selection by variance inflation factors
2) custom_functions - containing data cleaning and analysis functions
3) data - containing code for data extraction, cleaning, and vetting
4) gbm_optim - containing code to perform gradient boosted machine (GBM) optimization
5) had_diarrhea - containing code to launch model-based geostatistical models for diarrhea
6) mbg_central - containing model-based geostatistical code similar to code used across multiple projects at IHME
7) post_estimation - containing code to clean, vet, and analyze modeled estimates

# Other LBD work

To see what work our team has published, please visit the [main IHME LBD team page](http://www.healthdata.org/lbd).

To see other published LBD code repositories, please visit our [LBD code homepage](https://github.com/ihmeuw/lbd).