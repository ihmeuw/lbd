
#  Documentation for estimating HIV mortality using vital registration and small area estimation in Latin American countries

## Data
This section describes the input data required for modeling and the existing code to prepare and format this data correctly

#### Data sources
* Population data: Pulled from Worldpop raster using VR functions `vr_get_population`
* Mortality data: Pull from VR database using `vr_pull_cod`
* Covariate data: Pulled from MBG standard covariates. Uses central function `frag_agg_covs`
* Local shapefile: Pull from VR database using `vr_pull_shp`

### Data preparation
#### Code
Below is the location and a description of the R scripts required to prepare VR data for small area estimation, the process describes how the code is run and what is accomplished in each script. The data preparation consists of running the run_all prep functions, which also generates priors on VR completeness by launching the completeness code prep functions. These prep functions are seperated by VR completeness prep for Colombia, Costa Rica and Guatemala, and prep for Mexico and Brazil.

Prep data code: 

* `lbd_hiv/sae/vr/1_data_processing/run_all_data_prep.r`

Prep completeness code:

* `lbd_hiv/sae/vr/1_data_processing/vr_completeness_prep.r`
* `lbd_hiv/sae/vr/functions/vr_mex_bra_completeness.r`
* `lbd_hiv/sae/functions/vr_completeness_functions.r`

Extra functions: 

* `lbd_hiv/sae/functions/vr_gbd_comparison_functions.r` -- these functions are used for comparing unraked country estimates to GBD and loading GBD estimates
* `lbd_hiv/sae/functions/vr_model_prep_functions.r` -- these funtions are used for prepping HIV and all cause mortality from VR data 
* `lbd_hiv/sae/functions/vr_model_viz_functions.r` -- these functions are used in post estimation to create visualizations of model fit


#### Process
##### Prep all input data (`run_all_data_prep.r`)

* Run line by line through `run_all_data_prep.r` code, and add correct country (by iso3, bra = Brazil, col = Colombia, cri = Costa Rica, ecu = Ecuador, gtm = Guatemala, mex = Mexico). This script accomplishes the following steps:
* Define global parameters and pull in input data from VR database
	* Assign country iso3 and admin level for VR data
	* Determine rake level for GBD (National for most countries, Admin-1 for Mexico/Brazil and other countries that GBD estimates at admin1 level; as determined by `get_gbd_locs function`)
	* Pull the years of data available from formatted VR data
	* Load shapefile and shapefile path
	* Mark date of data prep and create a directory to save inputs
* Prepare population estimates from Worldpop for every modeled area-age-year
	* Process uses function `vr_get_population` to fractionally aggregated population to each admin unit in 5 year age bins. This function uses Worldpop data for the appropriate age group that is then raked to GBD population for the country or first administrative unit
	* Function returns the raking factors as well as population counts, which are saved in the output directory 
	* Additional checks are performed to ensure no modeled units are missing population data or have a population of zero.  A plot of raking factors by year and population changes by area is saved to the output directory. 
* Create a geographic aggregation file for post-estimation aggregation to administrative units from modeled areas 
	* For national geographic aggregation file use population output and add the country location id
	* For admin-1 aggregation file the code pulls the location metadata and creates data table assigning each modeled area to an admin-1
*  Prepares covariates: MBG standard covariates can also be fractionally aggregated, and these can be used as inputs to the model. 
	* Currently the code is set up to use the MBG standard covariates: access, urbanicity, nightlights, fractionally aggregated and population-weighted by stable admin unit
	* Creates a plot of covariate values for each modeled area over the time period, saved to the input data folder
* Prepares VR data: VR mortality data is specified using function `vr_prep_hiv_mortality`
* Saves shapefile, adjacency matrices and age weights from GBD 2017.
* Launch completeness code for generating priors for all countries except Costa Rica
* After final checks, dated input data folder is created with SAE settings at `<<<< FILEPATH REDACTED >>>>`


## SAE modeling
This section provides an overview of how to launch small area estimation model after data is formatted correctly using the processes outlined above. A more detailed explanation of the code and methods used for using a small-area model on administrative count data can be found at `/lbd_core/sae_central/readme.md`. This describes the required settings for modeling, as well as the input files and code overview. Code described below complements the code found in the `sae_central` folder. Currently, completeness priors can only be used for model 4 and 6 (though this could be modified to include more models). 

#### Code
Generic settings file (for defining model run parameters and location of input data)

*  `lbd_hiv/sae/vr/2_modeling/settings.csv`

List of covariates to use in modeling framework

* `lbd_hiv/sae/vr/2_modeling/cov_list.csv`

Launch script (for setting models and estimation code running)

* `lbd_hiv/sae/vr/2_modeling/launch.r`

Run all scripts (for running the launch script for multiple sets of models and multiple countries)

* `lbd_hiv/sae/vr/2_modeling/run_all_vr.r`

Post-estimation code (primarily for visualizing results)

* `lbd_hiv/sae/vr/3_post_estimation/launch_post_estimation.r`
* `lbd_hiv/sae/functions/*` -- includes functions for plotting results, completeness prep, and comparison to GBD

#### Process
##### Launch models (`run_all_vr.r` & `launch.r`)
* Currently this is set up to launch models by country and model type, where model type refers to the suite of small area models defined in sae_central folder. To launch a model, pass a data table of country name and model number and the `launch.r` script will run with selected model and country. Information on process will be saved to `<<<< FILEPATH REDACTED >>>>` for posteriority
* Run all script calls launch script, which creates an output directory for models located at `<<<< FILEPATH REDACTED >>>>`

##### Post estimation 
* Currently code is set up to make several different diagnostic plots in post estimation. After the model is completed, you can run `lauch_post_estimation.r` which will load in the table saved to `<<<< FILEPATH REDACTED >>>>`
* Diagnostics include comparisons to GBD model results at the national level, and time-trend plots at the first administrative level, as well as additional model parameter plots. Code for creating these diagnostics is found at `lbd_hiv/sae/functions/*`

### Model output

Model runs are stored here `<<<< FILEPATH REDACTED >>>>`. In this folder you can find:
* `/temp_dir/` folder which contains and errors and outputs folder for each job, and draws and estimates for the SAE model at each geographic resolution
* `/time_trend/` folder that has a time-trend plot of age-standardized HIV mortality for all admin-1 that fall within a country, as well as plots of age-specific mortality rate over time by each admin-1 unit
* `/rake_plots/` folder that has a comparison between the unraked model fit and GBD national mortality rate, faceted by year and by age group
* `est_all.rdata` which contains the model estimates (unraked)
* `est_all_raked.rdata` which contains the raked model estimates
* `model_fit.pdf` plots the hyperparameters of the SAE model
* `model_results.pdf` makes several plots of the unraked model results
* `model_results_raked.pdf` makes several plots of the raked model results 
