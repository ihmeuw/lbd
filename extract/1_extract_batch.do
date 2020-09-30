/***********************************************************************************************************
// Redacted                                                              
 Date: 7/13/2015
 Project: ubCov
 Purpose: Run Script
***********************************************************************************************************/
// Redacted

//////////////////////////////////
// Setup
//////////////////////////////////

// Redacted

clear all
set more off
set obs 1

// Settings
// Redacted
local topics hap //ENTER YOUR TOPIC HERE
// Redacted

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
// Redacted
local outpath = // Redacted
local thisvar = "cooking_fuel_mapped"

   get, vars
    levelsof var_name if topic_name == "`topics'", l(vars) clean
    local n : list sizeof vars
    get, codebook
    egen keep = rowmiss(`vars')
    keep if keep < `n'

     batch_extract, topics(`topics') ubcov_ids(`ubcov_ids') /// 
                 central_root(`central_root') ///
                 cluster_project(proj_geospatial) ///
                 output_path("`outpath'") ///
// Redacted
