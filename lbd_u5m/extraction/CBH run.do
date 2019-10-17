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
local central_root "<<<< FILEPATH REDACTED >>>>"								// This is set to the folder which contains the ubcov code - shouldnt need changing
global outpath "<<<< FILEPATH REDACTED >>>>"							// Set this as the folder you wish to save your extracts in.
local geographies "<<<< FILEPATH REDACTED >>>>"		//Set this filepath so be that where you have a combined geography codebook dta file saved
local topics cbh 																			// List of this topics you want to extract data for - shouldnt need changing
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
	gen country = substr(ihme_loc_id,1,3)
	rename year_end year
	rename ihme_loc_id iso3
	rename survey_name source
	cap gen pweight = .
	rename pweight weight
	rename survey_module survey
	tostring year, gen(temp_year)
	if regexm(source, "/") {  
	replace source = "COUNTRY_SPECIFIC"
	}
	tostring nid, gen(temp_nid)
	gen temp_id = iso3 + "_" + temp_year + "_" + source + "_" + survey +"_"+ temp_nid
	
	cap tostring geospatial_id, replace
	rename geospatial_id cluster_number
	cap tostring hh_id, replace
	cap tostring mother_id, replace
	cap tostring child_id, replace
	
	// Drop unwanted variables in datathe variables you want to keep in the dataset
	foreach var in age_month psu urban admin_1_urban_id admin_1_urban_mapped admin_2_mapped admin_2_id file_path hhweight smaller_site_unit year_start temp_year temp_nid admin_1_mapped admin_1_id iso3 urban age_day file_path year_start temp_year temp_nid _merge temp_id child_age_at_death_raw admin_1 admin_2 admin_3 admin_4 buffer uncertain_point iso3 latitude longitude age_of_death_units age_of_death_number hh_id, child_id, mother_id ever_married_weight int_month int_year geospatial_id mother_age_years mother_dob_cmc hh_resident{
		cap drop `var'
	}
	
	save "`filename'", replace emptyok
	
}

STOP

///////////
/// END ///
///////////
