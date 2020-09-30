/***********************************************************************************************************
// Redacted                                                              
 Date: 02/10/2020
 Project: ubCov
 Purpose: Run Script
***********************************************************************************************************/


//////////////////////////////////
// Setup
//////////////////////////////////

// Redacted

clear all
set more off
set obs 1

///////////////////////////////////////////////////////////////////////
/*  Arguments:
		- topic: your research topic
        - array: Ubcov_id. The id of the codebook row
// Redacted
    Optional:
        - keep: Keeps both raw data and extracted data, allows manual extraction check before final output.
        - bypass: Skips the extraction check, output the data 
        - run_all: Loops through all ubcov_ids in the codebook.
*/
////////////////////////////////////////////////////////////////////////

local topics hap
local array 1075
// Redacted
local options bypass //leave it blank if you don't want to keep, bypass, or run_all. If you are running on cluster, you have to use bypass as an option otherwise it will error.



///////////////////////////////////////////////////////////////////////
//Extraction
//////////////////////////////////////////////////////////////////////


// Load functions
cd "`central_root'"
// Redacted

// Load the base code for ubCov
// Redacted
foreach path in `paths' {
    local files : dir "`path'" files "*.do"
    foreach file in `files' {
        if "`file'" != "run.do" do "`path'/`file'"
    }
}

// Make sure you're in central
cd `central_root'

// Initialize the system
ubcov_path
init, topics(`topics')

// Launches extract
foreach number in `array'{
    local uid `number'
    run_extract `uid', `options'
    tostring year_start, gen(year)
    tostring year_end, gen(end_year)
    tostring nid, gen(nid_n)
    local filename = ihme_loc_id + "_" + survey_name + "_" + year + "_" + end_year + "_" + nid_n
    local filename = subinstr("`filename'", "/", "_",.)
    drop year end_year nid_n
	cd  "`outpath'"
    export delimited using "`filename'.csv", replace
}
