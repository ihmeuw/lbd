# HIV covariates for Sub-Saharan Africa

## Data

### Input datasets
Input datasets are often grouped by risky sexual practices (rsp) displayed below:

| Indicator | Types | Indicator Definition| Subset | Age Range | Archived |
| ---------------- | ----------------------------- | ---- |---------- |  ---------- | ---- |
| in_union | BOTH | True if marital status is currently married or living with partner; False for all other marital status values. |none| 15-49 | |
| partner_away | WN | partner_away | Keep only currently married or living with partner | 15-49 |  |
| condom_last_time | BOTH | condom_last_time | Keep only had sex in last 12 months | 15-49 | MN, WN |
| multiple_partners_year | MN, WN | True if num_partners_year is 0 or 1 | Drop all surveys that required constructing num_partners_year using spousal adjust. | 15-49 | BOTH |
| had_intercourse | WN | had_intercourse | none | 15-24 | BOTH, MN |
male circumcision (defined in men from 15-49 years old), and presence of STI symptoms in men and women (15-49). Therefore in many parts of the code there are three different groups for covariates: rsp, male circumcision, and sti, which we refer to as *covariate categories*. When refering to the exaustive list of covariates (in_union, partner_away, condom_last_time, multiple_partners_year, had_intercourse, male_circumcision and sti_symptoms) we refer to as *covariate indicator*. 
* Covariate survey microdata (to be collapsed): 
    * Extraction code: "lbd_hiv/data/covariate_surveys/microdata/[*Covariate Category*]/"
    * Data: "<<<< FILEPATH REDACTED >>>>".
* Male circumcision survey reports (appended onto microdata): 
    * Data: "lbd_hiv/data/covariate_surveys/reports/circumcision_extraction.csv"

### Data collapse and compilation

#### Code
* Collapse code: "lbd_hiv/mbg/covariates/1_collapse/1_[*covariate indicator*]_collapse.r"
* Collapse checks: "lbd_hiv/mbg/covariates/1_collapse/2_[*covariate indicator*]_collapse_checks.r"

#### Process

Each indicator has at least one collapse code and often as many as 3 (one for both men and women combined (BOTH), one for just men (MN), and one for just women (WN) - this is referred to as the collapse type which is only used in file names etc. if an indicator has multiple types). These collapse codes can be found in /lbd_hiv/2_collapse/rsp/1_\<indicator name>\_collapse_\<type>.R. The general steps in each collapse code are described below. 

1. Read in the most recent extracted rsp data from <<<< FILEPATH REDACTED >>>> and drop unnecessary variables.
2. Drop all surveys pre-1998.
3. All observations with missing pweights or missing indicator are dropped.
4. If type is MN or WN: keep only male or female respondents respectively. If type is BOTH: only keep surveys that asked both men and women the indicator question.
5. Drop observations from respondents outside the age range.
6. Drop surveys that don't include the full age range (determined by whether a survey includes at least one age_min-year-old respondent AND at least one age_max-year-old respondent).
7. Drop point data with missing latitude/longitude.
8. Drop data from countries outside of Africa.
9. Save the fully subset but pre-collapsed data to <<<< FILEPATH REDACTED >>>>.
10. Microdata is collapsed by survey and location. The median interview year is used as year. 
11. Create data coverage maps. Output: <<<< FILEPATH REDACTED >>>>.
12. Convert to counts if prevalence indicator. Resample polygons.
13. Compare the list of surveys and their sample sizes at this point to that saved before the data were collapsed (i.e., saved after step 6).
14. Compare the data to that produced by the last run of the collapse code and highlight surveys that have been added, dropped, or where data changed
15. Add columns to final data indicating the date of the geomatched input data and the current date for collapse. Write the geomatched data date and the collapse date to the date log at <<<< FILEPATH REDACTED >>>>.
16. Save final data to <<<< FILEPATH REDACTED >>>> and a version to the archive at <<<< FILEPATH REDACTED >>>>

Note: We decided to drop condom_every_time_3_partner and condom_every_time_2_partner before even writing collapse codes due to lack of data. After writing collapse codes for client_sex_worker_3_partner and client_sex_worker_2_partner and looking at the coverage plots, we decided to also drop these two indicators due to lack of data and such low prevalence levels. These collapse codes are now in /lbd_hiv/2_collapse/rsp/archive/.

Subsetting, indicator creation, etc. specific to each indicator are summarized in the following table. As we decide which types to model for each indicator, the rejected types are moved to the 'Archived' column and the corresponding collapse codes are moved to /lbd_hiv/2_collapse/rsp/archive/.


## MBG Modeling

### Code

Config files (for defining model run parameters)

 - Each modeled indicator shares a similar config file, found at 'lbd_hiv/mbg/covariates/2_modeling/config_*.csv
	 -  config_default_pc.csv -- standard penalized complexity prior (used for final results)
	 - config_default_pc_strict_rho.csv -- more strict parameterization of time correlation
	 - config_loose_pc.csv -- less restrictive PC prior specification
	 - config_loose_pc_strict_rho.csv -- less restrictive PC prior specification with strict parameterization of time correlation
	 - config_inla_default.csv -- INLA default inla.spde2.matern() parameterization
	 - config_inla_default_strict_rho -- INLA default inla.spde2.matern() parameterization with strict time correlation

- Covariates list (though ignored): 'lbd_hiv/mbg/covariates/2_modeling/cov_list.csv'

Launch script (for setting models and estimation code running)

-   "lbd_hiv/mbg/covariates/2_modeling/launch.r"

Run all scripts (for running the launch script for multiple sets of models, ie, using multiple config files)

-   "lbd_hiv/mbg/covariates/2_modeling/run_all.r" (runs in-sample models to compare different hyper-prior specifications)

Post-estimation code (primarily for visualizing results; in most cases these are functions that are called automatically from the launch script for each new model run)

-   Shared code (across indicators): 'lbd_hiv/mbg/functions/*'
-   Covariate code (typically from determining which sex to include in covariates): 'lbd_hiv/mbg/covariates/3_post_estimation.r






