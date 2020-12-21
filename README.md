# Oral rehydration therapies - Senegal, Mali, and Sierra Leone - 2020

This repo contains cleaned code used to generate coverage estimates of any oral rehydration solution (any ORS), only recommended home fluids (only RHF), and no oral rehydration therapy (no ORT) for children under the age of five with diarrhoea from the manuscript 'Oral rehydration therapies in Senegal, Mali, and Sierra Leone: a spatial analysis of changes over time and implications for policy'.

The repo contains six directories:
1) covariate_selection - containing code to perform covariate selection by variance inflation factors
2) custom_functions - containing data cleaning and analysis functions
3) gbm_optim - containing code to perform gradient boosted machine (GBM) optimization
4) mbg_central - containing model-based geostatistical code similar to code used across multiple projects at IHME
5) ors - containing code for data extraction and cleaning, as well as code to launch model-based geostatistical models for any ORS (ORS alone or in combination with RHF), only RHF (RHF alone), and no ORT (no ORS or RHF).
6) post_estimation - containing code to clean, vet, analyze, and plot modeled estimates, including analyses of changes over time before and after events and policy changes that may have impacted diarrhea treatment

# Other LBD work

To see what work our team has published, please visit the [main IHME LBD team page](http://www.healthdata.org/lbd).

To see other published LBD code repositories, please visit our [LBD code homepage](https://github.com/ihmeuw/lbd).