CSV files in this folder contain data extracted from survey reports and published literature. 
These are generally data sources where microdata are not available, but the estimates reported are similar in nature to what we would calculate from the microdata, if we had access. 
We extract these data in this form so that they can be combined with collapsed microdata from other sources. 

# Male Circumcision (circumcision\_extraction.csv)

### Data: 
Self reported male circumcision prevalence from population surveys or other sources intended to be representative of the general adult population. 
Ideally, we are looking for information on the age group 15-49. 

### Variables: 

| Variable | Description | Required? | 
| ---------------    | ------------------------------------ | ------- | 
| nid                | Source NID.| Required. |
| file_path          | J drive location of the file from which data were extracted. | Required. | 
| survey_name        | Name of survey. | Required. | 
| survey_series      | Name of series. For major series, this follows normal conventions (e.g., MACRO_DHS), but otherwise should usually be COUNTRY_SPECIFIC. | Required. | 
| country            | Country iso3 code. | Required. | 
| start_year         | Year data collection began. | Required. | 
| end_year           | Year data collection ended. | Required. | 
| int_year           | Reference year for data collection (usually this will be the midpoint of a reported range). | Required. | 
| int_month          | Reference month for data collection (usually this will be the midpoint of a reported range). | Required. | 
| start_age          | Beginning of the age range data refer to, in years. | Required. | 
| end_age            | End of the age range data refer to, in years. | Required. | 
| sex_id             | Sex ID for the data, coded 1 = males, 2 = females, 3 = both combined. | Required. | 
| location           | Location information, as presented in the data source (preferably cut and paste, though you may possibly need to remove special characters). | Required. | 
| location_type      | Descriptor of the location type, as presented in the data source (e.g., village, city, district, region, etc.). | Required. | 
| point              | Indicator for point vs polygon data, coded 0 = polygon and 1 = point. | Required. |
| shapefile          | Shape file containing the polygon this data refers to.  | Required if point = 0. | 
| location_code      | Location code in the specified shape file of the polygon this data refers to. | Required if point = 0. |
| admin_level        | Administrative level the data refer to. | Required if point = 0 and the polygons correspond to actual administrative boundaries. | 
| latitude           | Latitude of the specified location. | Required if point = 1. | 
| longitude          | Longitude of the specified location. | Required if point = 1. | 
| male_circumcision  | Reported HIV prevalence (%). | Required. | 
| mc_unknown         | *This column is to be filled out for PHIA reports for transparency in prevalence total.
| N                  | Reported sample size. | Required. | 
| date_source        | Source of information for start_year, end_year, int_month, and int_year. Usually this will be a quote (with corresponding page number\*) and any information about how this was transformed into the values entered. | Required. | 
| hiv_test_source    | Source of information for male_circumcision. Usually this will be a table or figure number (with corresponding page number\*) and any information about how this was transformed into the values entered. | Required. | 
| N_source           | Source of information for N. Usually this will be a table or figure number (with corresponding page number\*) and any information about how this was transformed into the values entered. | Required. | 
| data_notes         | Any notes about the underlying data or extraction process. | Optional. | 
| extractor          | netid for the person who extracted this data. This could be more than one person if data are updated or added to a particular row. | Required. | 

\**Note: page numbers should be actual page numbers, as displayed in the document, and not PDF page numbers.* 
