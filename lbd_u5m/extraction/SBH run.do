/***********************************************************************************************************
 Project: ubCov
 Purpose: CBH Run Script																					
***********************************************************************************************************/


//////////////////////////////////
// Setup
//////////////////////////////////

clear all
set more off
set obs 1

// Set locals
local central_root "<<<< FILEPATH REDACTED >>>>"														// This is set to the folder which contains the ubcov code - shouldnt need changing
global outpath "<<<< FILEPATH REDACTED >>>>"										// Set this as the folder you wish to save your extracts in. In this folder you require one folder called 'WN' and one called 'Aggregated'
local geographies "<<<< FILEPATH REDACTED >>>>"								//Set this filepath so be that where you have a combined geography codebook dta file saved
local topics sbh																									// List of this topics you want to extract data for - shouldnt need changing
local array <<<< REDACTED >>>>

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

////////////////////
// Run extraction //
////////////////////


foreach number in `array'{
    local i `number'
    run_extract `i', bypass bypass_map 
	rename geospatial_id cluster_number
	rename year_end year
	gen country = substr(ihme_loc_id,1,3)
	rename ihme_loc_id iso3
	rename survey_name source
	cap gen pweight = .
	rename pweight weight
	rename survey_module survey
	tostring year, gen(temp_year)
	if regexm(source, "/") {  
	replace source = "COUNTRY_SPECIFIC"
	}
	tostring nid, gen (temp_nid)
	gen temp_id = iso3 + " " + temp_year + " " + survey + " " + source + " " + temp_nid
		
	//format variabels to %20.0f
	tostring cluster_number, replace format("%20.0f")
	cap tostring line_id, replace format("%20.0f")
	cap tostring hh_id, replace format("%20.0f")
	
	global filename = "$outpath"+temp_id+".DTA"
	
	
	// Drop unwanted variables in datathe variables you want to keep in the dataset
	foreach var in psu hhweight iso3 int_year admin_1_mapped admin_1_id admin_2_mapped admin_2_id smaller_site_unit urban int_month file_path temp_nid temp_year sex_id year_start _merge admin_1 admin_2 admin_3 admin_4 buffer uncertain_point iso3 temp_id admin_1_mapped admin_5 age_day geospatial_id int_year age_month{
		cap drop `var'
	}
	
	// if age_year isnt in the dataset (i.e only categories available) gen a missing variable
	cap gen age_year = .
	
	save "$filename", replace emptyok
	
}
STOP

