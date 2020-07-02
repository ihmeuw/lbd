

# Male Circumcision documentation

## Data identification strategy 
- STATCompiler to identify DHS, AIS with relevant variables
- Review surveys extracted for other HIV topics
  - Found via ubCov codebooks (hiv_anc, hiv_biomarker, hiv)
- Review surveys identified for rsp by [Cooper](http://internal-ghdx.healthdata.org/hivaids-geospatial-sexual-activity)
- Cooper project - [HIV/AIDS Geospatial - Male Circumcision](http://internal-ghdx.healthdata.org/hivaids-geospatial-male-circumcision)
  - Keywords: Circumcision

Initial phase just extracted microdata, although kept an ongoing list of sources including male circumcision for which only reports were available that were extracted later (described in Report Extraction).

### Inclusion criteria
- INCLUDE any source asking male respondents if they personally are circumcised
- EXCLUDE sources only asking about circumcision for persons other than the respondent (e.g., do you plan to circumcise your son)


## Microdata extraction - UbCov
### Indicators
The sole indicator for this topic (described below) is extracted from ubCov via the [male_circumcision codebook](https://docs.google.com/spreadsheets/d/1T9j7mV9zYNcz05UZunIk6dqiHhAranyFwpAduh0eaT0/edit?usp=drive_web&ouid=108770465896075084680). There is no custom code for male circumcision - the indicator is generated using ubCov's standard indicator generation code. The standard set of indicators (e.g., sex, age, interview dates, geographic information, etc.) is also extracted. Because this set is common to all extractions, only topic specific indicators are described in the following table.  
 
| Indicator | Definition | Coding | Source |
| ---------------- | ----------------------------- | ---------- | ---------- |
| male_circumcision | Has the respondent been circumcised (medical or traditional) | 0 False 1 True | male_circumcision cb via standard indicator generation |

Some surveys include both men and women, with the circumcision variable representing either male or female circumcision depending on the sex of the respondent. All responses are pulled via ubCov and the subsetting to just males is done in postprocessing (described below). 

### Extraction code
The extraction code for this topic is found at lbd_hiv/data/covariate_surveys/microdata/male_circumcision/1_male_circumcision_extraction.do and saves the extracted datafiles to <<<< FILEPATH REDACTED >>>>. The extraction code automatically calculates the current date, but the folder in ubCov_extractions must be manually created first. 

### Extracted data
For each module extracted via ubCov, a corresponding standardized .dta file is extracted to  <<<< FILEPATH REDACTED >>>> is the date the extraction was run. Each dataset includes the following variables in addition to the circumcision indicator: nid, nid_n, survey_name, ihme_loc_id, year_start, year_n, yeard_end, end_year_n, survey_module, file_path, strata, psu, pweight, geospatial_id, hh_id, line_id, sex_id, age_year, age_month, int_month, int_year, urban, admin_1, admin_2. 

### Special cases
Ethiopia surveys require year and month adjustments due to Ethiopia's use of a non-Gregorian calendar. The year adjustment is done in ubCov (8 years are added to dates via basic) but the true difference is 7 years and 8 months so an additional 4 months must be subtracted from dates. This month adjustment is done in postprocessing (described below). NOTE: This may at some point be moved centrally to ubCov, at which point the manual adjustment will need to be removed from postprocessing.  


## Microdata extraction - Post-processing
Postprocessing code is found at lbd_hiv/lbd_hiv/data/covariate_surveys/microdata/male_circumcision/2_male_circumcision_postprocessing.R.

Each dataset from the most recent extraction is read into R and bound together into a single dataset. That dataset is then put through a series of formatting and subsetting steps as follows:

1. Specific survey issues are addressed 
  - ZMB Sexual Behavior Survey 2005 (nid 27987) has a series of additional clusters appended onto the data that don't correspond to actual interviews (clusters greater than 105). These additional clusters are dropped. (In other topics, these are dropped in custom code, but since male_circumcision doesn't have any custom code, this must be done in postprocessing.)
  - Interview dates are adjusted for Ethiopia surveys: an additional 4 months must be subtracted from interview dates to compensate for the fact that only whole year adjustments can be made in ubCov (see note in above).

2. The data is subset
  - Keep only male responses. Because some surveys combine male and female responses about circumcision into one variable, we must restrict the data to only male responses (sex_id == 1).

Once the data has been thus transformed, it is geomatched. The full set of geography codebooks are read in and bound onto the data. Any observations missing geographical information are written to a 'geography_matching_to_do.csv' file that should be manually reviewed and the necessary changes should be made to the geography codebooks to geomatch as many of these observations as possible. 


## Microdata extraction - checks
### Checks code
#### Checks report
Checks code is found at lbd_hiv/data/covariate_surveys/microdata/male_circumcision/3_male_circumcision_checks_report.Rmd. 

Knitting the checks_report file will produce an html report that summarizes the results of a series of checks on the most recent extracted postprocessed dataset. It takes four parameters:

1. \<comments_path\>: A filepath to a csv where comments about each survey can be recorded by nid. These comments will be read-in and displayed for the corresponding surveys. 
Default: "<<<< FILEPATH REDACTED >>>>"

2. \<vetted_path\>: A filepath to a csv where the nids of previously vetted surveys can be recorded. Once a survey has been reviewed and all issues addressed or explained, its nid should be recorded in this csv. This list of vetted surveys will be read-in and the status for each survey will be displayed.
Default: "<<<< FILEPATH REDACTED >>>>"

3. \<input_version\>: The date of the most recent postprocessed dataset. Format: "mm-dd-yyyy". 

4. \<fast_checks\>: A boolean indicating whether or not the report should be created using 'fast-checks'. If TRUE, diagnostics will be printed only for unvetted surveys. If FALSE, diagnostics will be printed for all surveys. 
Default: TRUE

The checks included in the report are summarized below, along with the corresponding sections of the report: 

1. Check geography information:
    - For all point data, check that each point is within the borders of the country. Print a table of all points outside the correct country borders, with points more than 5 km distant highlighted in orange and points more than 10 km distant in red. 
    - For polygon data, check whether any observations are missing shapefiles or gaul codes. Print tables of all observations that are. Check whether any observations use invalid shapefiles (i.e., shapefiles that are not stored in the shapefile directory). Print a table of all observations that are, with the invalid shapefile names highlighted in red. Check that the gaul code each observation references can be found in the corresponding referenced shapefile. Print a table of all surveys using mismatched gaul codes.

2. Check data for each survey:

    (If fast checks is turned on, these diagnostics will only be printed for unvetted surveys)
    - A warning is raised for:
        - any missing geography information
        - any missing individual weights
        - any missing interview years
        - any interview years outside of the survey range
        - any missing sex ids
        - any missing age years
        - an incomplete age range (i.e. not the full 15-49)
    - If the survey will fail any of the conditions at the beginning of collapse code and be dropped, that fact is highlighted at the beginning of the survey diagnostics
    - A plot of the male circumcision data for that survey is printed to be verified visually
    - A table of the sample size and missingness for male circumcision is included, with any missingness over 5% highlighted
    - Any comments for the survey are printed (with the idea being that as each survey is vetted, any remaining warnings, oddities in the plot, or missingness in the indicator will be explained in comments to allow for quick future review). 


## Report comparison
Comparison code can be found at /lbd_hiv/data/covariate_surveys/microdata/male_circumcision/4_male_circumcision_report_comparison.R.

Report data was pulled from statcompiler for DHS/AIS surveys and manually pulled for other surveys where available. It is stored in the same repo folder. A description of the various data files follows:

- 4_statcompiler_mc_db.xlsx: download of sample sizes for male circumcision indicator from statcompiler
- 4_statcompiler_mc_indicators.xlsx: download of male circumcision prevalences from statcompiler
-4_statcompiler_mc_prev_n.xlsx: download of prevalence and sample size from statcompiler for Non-SSA and SSA countries.
- 4_non_statcompiler_report_data_mc.xlsx: sample size and indicator prevalence where available for other surveys


This data is read into the code and compared to numbers calculated from the extracted, postprocessed data for these indicators. Sample size is also compared. Plots of the results are output for manual review to <<<< FILEPATH REDACTED >>>>.

##Before Collapse Checks Comparison and Time Series Plots
Comparison code can be found at /lbd_hiv/data/covariate_surveys/microdata/circumcision_before_collapse_check.R.
Input:
The geomatched dataset is the one created in post-processing.  
The pre-collapsed dataset is the cleaned geo-matched dataset before the resampling polygons function in the collapse code. 

-Enter the  old and new input dates for pre-collapsed and geomatched datasets. 
-Data is loaded from previous and latest pre-collapsed dataset.
	- Checks which NIDS were added, dropped, and changed in pre-collapsed dataset.
- Data is loaded from previous and latest geomatched dataset.
	- Checks geomatched dataset for NID changes, to validate and narrow down pre-collapsed NID changes(if there are any).
- Create time-series-plot



## Report extraction 
### Data identification strategy 
- Any sources via the microdata identification that were marked as report only
- Review all reports extracted for HIV test
- Review all sources identified in the [HIV Blood Test Cooper project](http://internal-ghdx.healthdata.org/hivaids-geospatial-blood-tests-or-medical-tests)

### Inclusion criteria
- INCLUDE any source reporting the prevalence of male circumcision from population surveys or other sources intended to be representative of the general adult population at a geographic level below admin 0 for any age grouping (ideally 15-49 but all age ranges were extracted)


## Extracted data
Extracted report data for male circumcision is located in /lbd_hiv/data/covariate_surveys/reports/circumcision_extraction.csv. Included in the same folder is a readme detailing the description of each variable included in the extraction csv.


##Before Collapse Checks Comparison and Time Series Plots
Comparison code can be found at '/lbd_hiv/data/covariate_surveys/microdata/male_circumcision/5_pre_collapse_checks_time_series.r'
Input:
The folder-dates for new geomatched dataset and old geomatched dataset. 

1. Enter the old and new input dates for geomatched datasets. 
2. Data is loaded from previous and latest pre-collapsed dataset.
3. Checks which NIDS were added, dropped, and changed in pre-collapsed dataset.
4. Checks for difference in prevalence and sample size and diagnostic information.
5. Creates time-series plots for new surveys, and all surveys. 


## Data collapse and compilation 
Collapse code for circumcision is found at /lbd_hiv/mbg/covariates/1_collapse/male_circumcision/1_circumcision_collapse.R. The code does the following:

1. Read in the most recent extracted circumcision data from <<<< FILEPATH REDACTED >>>> and drop any unnecessary variables.
2. Drop observations from respondents outside the 15-49 age range.
3. Drop surveys that don't include the full age range 15-49 (determined by whether a survey includes at least one 15-year-old respondent AND at least one 49-year-old respondent)
4. Drop any observations that are missing pweight or male_circumcision
5. Drop point data with missing latitude/longitude
6. Drop data from countries outside of Africa
7. Collapse microdata by survey and location. The median interview year is used as year. 
8. Add in the survey report data
  - Keep only data with 15-49 age group. 
  - Drop national estimates
  - Drop any nids for which we also have microdata
  - Approximate the effective sample size to ensure report is not unfairly masking microdata
9. Create data coverage maps. Output: <<<< FILEPATH REDACTED >>>>
10. Convert to counts and resample
11. Compare the list of surveys and their sample sizes at this point to that saved before the data were collapsed (i.e., saved after step 6)
12. Compare the data to that produced by the last run of the collapse code and highlight surveys that have been added, dropped, or where data changed

#### Checks
1. Some checks are built into  the collapse code itself as described in steps 11 and 12 above. The output should be carefully reviewed to verify any survey additions, drops, or changes were intentional.
2. Collapse checks code (/lbd_hiv/mbg/covariates/1_collapse/male_circumcision/2_circumcision_collapse_checks.R) produces histograms of sample size and prevalence, scatter plots of sample size and prevalence, and maps of cluster locations for review. 

#### Exclusions
Dropping Nigeria National HIV/AIDS and Reproductive Health Survey 2007 (nid 325046) and 2012 (nid 324443) and Nigeria Living Standards Survey (nid 151719) because the data are too different from DHS.
