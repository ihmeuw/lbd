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
local central_root "<<<< FILEPATH REDACTED >>>>"
local topics ort //ENTER YOUR TOPIC HERE

// Load functions
cd "`central_root'"
do "`central_root'<<<< FILEPATH REDACTED >>>>load.do"

// Load the base code for ubCov
local paths  `central_root'<<<< FILEPATH REDACTED >>>> `central_root'<<<< FILEPATH REDACTED >>>>
foreach path in `paths' {
    local files : dir "`path'" files "*.do"
    foreach file in `files' {
        if "`file'" != "run.do" do "`path'/`file'"
    }
}

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
// Run extraction
/* Launches extract
    Arguments:
        - ubcov_id: The id of the codebook row
    Optional:
        - keep: Keeps 
        - bypass: Skips the extraction check
        - run_all: Loops through all ubcov_ids in the codebook.
*/
//Enter path where you'd like your extracts saved below between the quotes
local outpath = "<<<< FILEPATH REDACTED >>>>"

///////////////////////////////////////////////////////////////////////////////////////////////

//Enter the ubCov ID of your survey after "array"
local array 5404 5553 5944 6014 6099 6183 6245 6329 // 10470 10469// Insert ubCov ID here
foreach number in `array'{
    run_extract `number', bypass //keep  
    tostring year_start, gen(year)
    tostring year_end, gen(end_year)
    tostring nid, gen(nid_n)
    local filename = survey_name + "_" + nid_n + "_" + survey_module + "_" + ihme_loc_id + "_" + year + "_" + end_year
    local filename = subinstr("`filename'", "/", "_",.)
    drop year end_year nid_n
    cd  `outpath'
    save "`filename'", replace //if you upgrade to Stata 14 will need to change this to saveold
}
