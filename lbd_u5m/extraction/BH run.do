/***********************************************************************************************************
 Project: ubCov
 Purpose: Run Script
																					
***********************************************************************************************************/


//////////////////////////////////
// Setup
//////////////////////////////////

clear all
set more off
set obs 1


// Settings
// Set locals
local central_root "<<<< FILEPATH REDACTED >>>>"								// This is set to the folder which contains the ubcov code - shouldnt need changing
global outpath "<<<< FILEPATH REDACTED >>>>"							// Set this as the folder you wish to save your extracts in.
local topics birthhistories																	// List of this topics you want to extract data for - shouldnt need changing
local array <<<< REDACTED >>>>
						 													// List the UbCov ID's you wish to extract here

// Set working directory
cd `central_root'

// Load functions
do "<<<< FILEPATH REDACTED >>>>"

// Make sure you're in central
cd `central_root'

// Initialize the system
/* 
	Brings in the databases, after which you can run
	extraction or sourcing functions like: new_topic_rows

	You can view each of the loaded databases by running: get, *db* (eg. get, codebook)
*/

ubcov_path
init, topics(`topics')

// list of ubcov ids to process


foreach number in `array'{
    local i `number'
    run_extract `i', bypass   
	cap gen pweight = .
	cap gen admin_2 = .
	
	// renaming of vars to comply with current model reqs
	rename geospatial_id cluster_number
	rename hh_id household_number
	rename year_end year
	rename ihme_loc_id country
	rename survey_name source
	rename pweight weight
	rename survey_module survey
	gen caseid = mother_id
	rename mother_id mid
	rename child_no childs_line

	// drop module ('survey'), and rename 'source' as 'survey'
	drop survey
	gen survey = source
	
	if regexm(source, "/") {  
	replace source = "COUNTRY_SPECIFIC"
	}
	
	// generate a unique survey name for output
	tostring nid, gen(temp_nid)
	tostring year, gen(temp_year)
	gen temp_id = country + " " + temp_year + " " + source + " " + temp_nid
	
	// create required variables and ensure existing ones are strings
	foreach var in admin_1 admin_2 strata mothers_line mid{
	capture confirm new variable `var'
	if _rc == 0 {
		gen `var' = .
		}
	tostring `var', replace	
	}
	
	capture confirm variable ced_male
	if _rc {
		gen ced_male = .
	}
	capture confirm variable ced_female
	if _rc {
		gen ced_female = .
	}
	capture confirm variable age_of_death_units
	if _rc {
		gen age_of_death_units = .
	}
	capture confirm variable age_of_death_number
	if _rc {
		gen age_of_death_number = .
	}
	capture confirm variable children_alive_male
	if _rc {
		gen children_alive_male = .
	}
	capture confirm variable children_alive_female
	if _rc {
		gen children_alive_female = .
	}

	keep nid survey country year cluster_number household_number mothers_line childs_line mothers_age mothers_age_group child_alive birthtointerview_cmc child_age_at_death_months ceb ced weight mid admin_1 admin_2 interview_date_cmc child_dob_cmc temp_id age_of_death_units age_of_death_number strata temp_id ced_male ced_female children_alive_male children_alive_female

	// reorder variables to match previous formatting
	cap order nid survey country year cluster_number household_number mothers_line childs_line mothers_age mothers_age_group child_alive birthtointerview_cmc child_age_at_death_months ceb ced weight mid admin_1 admin_2 interview_date_cmc child_dob_cmc temp_id age_of_death_units age_of_death_number strata ced_male ced_female children_alive_male children_alive_female
	
	local filename = "$outpath"+temp_id+".DTA"
	
	drop temp_id
	
	save "`filename'", replace emptyok
	
}
STOP






