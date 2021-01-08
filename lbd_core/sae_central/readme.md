This folder contains code for running small area models on administrative count data (e.g., deaths registration data, notification data) to estimate rates by area, year, sex, and age. 

## Settings

The code in this folder is intended to be generic; different datasets and model options are provided via a settings CSV file called 'settings.csv'. This settings file should contain the following arguments: 

| Argument | Description | Required/Optional |
| :------ | :--------------------------------------------------------------| :-------------------- | 
| model | The SAE model to run. This matches the code in the 'models' folder. Current options are 1, 1b, 2, 2b, 3, 4, 5, 6, and 7. All models are described in 'models/models.tex' | Required |
| area_var | The name of the variable in all input folders denoting the areas to be modeled | Required |
| years | The range of years to model | Required |
| ages | The age groups to model | Required |
| sexes | The sexes to model | Required |
| covars | Area- and year-level covariates to include in the model | Optional |
| covars_as | Area-, year-, age-, and sex-level covariates to include in the model | Optional | 
| covars_trans | Transformations to apply to the covariates in covars and/or covars_as | Optional |
| n.sims | The number of posterior draws to generate and use | Required |
| events_file | File path to the event counts data | Required |
| pop_file | File path to the population data | Required |
| covar_file | File path to the area- and year-level covariates data | Required only if 'covars' is specified |
| covar_as_file | File path to the area-, year-, age-, and sex-level covariates data | Required only if 'covars_as' is specified |
| geoagg_files | File path(s) to population files relating the modeled areas to geographic aggregates | Optional |
| raked | If raking, named character vector that contains the country iso3 code, cause_id, gbd_round_id and rake_level. The rake level must be one of the geographic aggregates specified in geoagg_files. | Optional |
| adjmat_file | File path to the adjacency matrix | Required |
| shape_file | File path to the area-level shape file | Required |
| age_std_file | File path to the weights used for generating age-standardized estimates | Required |
| temp_dir | File path to the directory where intermediate and draws files are stored | Required |

An example settings file is located in this directory. 

## Input files 

All input files must be formatted properly in terms of the object type, object name, variables, and dimensions. All data should be stored as .rdata files except the age standard weights (age_std_file), which should be a .csv. Examples are available in '<<<< FILEPATH REDACTED >>>>/inputs/' corresponding to the settings.csv file in this directory. In summary: 

1. Events data (events_file) 
    - A data.table named 'events' 
    - One row per area, year, sex, and age (no duplicates allowed, excluded rows are assumed to correspond to 0 counts) 
    - Columns for area, year, sex, age, and events 

2. Population data (pop_file)
    - A data.table named 'pop'
    - One row per area, year, sex, and age (no duplicates allowed, must be square)
    - Columns for area, year, sex, age, and pop

3. Covariates data by area and year (covar_file)
    - A data.table named 'covar'
    - One row per area and year (no duplicates allowed, must be square)
    - Columns for area, year, and each covariate specified in 'covars'

4. Covariates data by area, year, sex, and age (covar_as_file)
    - A data.table named 'covar'
    - One row per area, year, sex, and age (no duplicates allowed, must be square)
    - Columns for area, year, sex, age, and each covariate specified in 'covars_as'

5. Geographic aggregation files (geoagg_files)
    - A data.table named 'weights'
    - One row per aggregated area, corresponding base area, year, sex, and age (no duplicates allowed) 
    - Columns for aggregated area (level), base area (area_var), year, sex, age, and pop

6. Adjacency matrix 
    - A dgTMatrix named 'adjmat'
    - One row and column for each area, with 0's denoting non-neighbors and 1's denoting neighbors. 

7. Shape file
    - A SpatialPolygonsDataFrame (object name doesn't matter)
    - One polygon and one row in the attribute table per area

8. Age standard weights (age_std_file)
    - One row per age group (no duplicates allowed)
    - Columns for age and wt

Some notes about the variables in these various files: 

- The area variable can be named anything and is specified via the 'area_var' argument in the settings file. 
  The area variable must be encoded such that areas are numbered sequentially from 0 and this must be consistent across all of the input files. 
  In the adjacency matrix, the rows and columns must be ordered in this same way. 
  The shape file must also be sorted according to this ordering, and the 'ID' slots for each polygon must use this encoding for the labels. 

- Variables for year, sex, and age must be named 'year', 'sex', and 'age'. 
  Year should be coded as an integer calendar year. 
  Sex should always be coded '1' for males and '2' for females (estimates for both sexes combined are added in the code and will be labeled '3'). 
  Age should be coded using the beginning of each age group. 
  Year, sex, and age coding must be consistent across all input files and all input files must contain all years, sexes, and ages as needed. 
  Additional years, sexes, and ages in the input files will generally be ignored. 

- Variables for events (in events_file) and population (in pop_file and geoagg_files) must be called 'events' and 'pop', respectively. 

- Variables for the covariates in covars_file or covars_as_file must match the covars or covars_as arguments, respectively. 

Basic formatting and consistency checks are carried out by the `check_settings()` function in `settings.r`. 
This is run by default when new models are launched, but can also be run manually to check a new settings file. 

## Code overview

### Main code 
- `submit.r` -- submit all jobs required to fit models and generate predictions based on the settings provided in settings.csv. 
- `models/prep_inputs.r` -- load, merge, and format the inputs files in preparation for modeling 
- `models/fit_mod_#.r` (and `models/mod_#.CPP`) -- fit the SAE model. The .r files are run directly and call on the corresponding .CPP files. 
- `models/plot_mod.r` -- plot the fitted SAE model parameters 
- `models/pred.r` -- generate posterior draws and summary measures of the area-year-sex-age-level rates based on the fitted SAE model 
- `post_estimation/agg_geos.r` -- generate posterior draws and summary measures of the area-year-sex-age-level rates for aggregated geographies based on population-weighted averages. 
- `post_estimation/agg_sex.r` -- generate posterior draws and summary measures of the area-year-age-level rates for both sexes combined based on population-weighted averages. 
- `post_estimation/compile_estimates.r` -- compile all summary estimates into a single output file. 

### Functions
- `qsub.r` -- a function for submitting jobs, used in `submit.r` 
- `settings.r` -- functions for loading settings, used in all code, and checking the settings and input files, used in `submit.r`
- `models/lcar_strmat.hpp` -- function for generating a LCAR structure matrix, used in `models/mod_#.CPP` 
- `models/get_params.r` -- function for extracting estimated model parameters from a saved model object, used in `models/plot_mod.r` 
- `post_estimation/calc_all_ages.r` -- function for calculating crude and standardized estimates for all ages combined, used in `pred.r`, `agg_geos.r`, and `agg_sex.r` 
- `post_estimationg/collapse.r` -- function for calculating summary measures from draws, used in `pred.r`, `agg_geos.r`, and `agg_sex.r `

## Running a model 

To run a model: 

1. Create (or update) a personal clone of the mbg repo containing the sae_central code on the cluster. Generally this will go in the '<<<< FILEPATH REDACTED >>>>' directory in a folder labeled with your netid, e.g.: '<<<< FILEPATH REDACTED >>>>'

2. On the cluster, navigate to the sae_central folder in this directory, e.g.: 

    `cd '<<<< FILEPATH REDACTED >>>>'`

3. qsub `submit.r`, providing four arguments: (1) the main directory, where 'settings.csv' is located and where final output will be saved; (2) the type of run, i.e. 'models' (only the models and initial predictions), 'post_estimation' (everything that comes after models and initial prediction), or 'all' (all steps); (3) a logical indicating whether or not this is a resubmission (i.e., if only jobs where output is missing should be submitted); and (4) the name of the cluster project these jobs should run under. For example: 

    `qsub -N sub -cwd r_shell.sh submit.r <<<< FILEPATH REDACTED >>>> all F proj_geospatial` 

Assuming everything runs as intended, intermediate output (the combined, prepped data; draws files; non-compiled summary estimates) will be saved in the temp directory indicated by the 'temp_dir' argument in the settings file and final output (the plots of the model parameters and the compiled summary estimates) will be saved in the main directory specified in the qsub (i.e., where your settings.csv file also is). Additionally, a record of the settings and jobs submitted will be saved in the temp directory. Usually the main directory will be on the J drive and the temp directory will be on the scratch space, though this is not required. Dating both folders is often useful to keep track of model runs, but this is also not required. 
