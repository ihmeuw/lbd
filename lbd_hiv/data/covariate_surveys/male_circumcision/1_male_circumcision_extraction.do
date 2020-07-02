/***********************************************************************************************************
 Author: Patrick Liu (pyliu@uw.edu) edited by Manny Garcia (gmanny@uw.edu) and Ruby Liu (hcrliu@uw.edu)                                                               
 Date: 02/10/2020
 Project: ubCov
 Purpose: Run Script
***********************************************************************************************************/


//////////////////////////////////
// Setup
//////////////////////////////////

if c(os) == "Unix" {
    local j "<<<< FILEPATH REDACTED >>>>"
	local h "<<<< FILEPATH REDACTED >>>>"
	local l "<<<< FILEPATH REDACTED >>>>"
    set odbcmgr unixodbc
	local central_root "<<<< FILEPATH REDACTED >>>>"
}
else if c(os) == "Windows" {
    local j "<<<< FILEPATH REDACTED >>>>"
	local h "<<<< FILEPATH REDACTED >>>>"
	local l "<<<< FILEPATH REDACTED >>>>"
	local central_root "\\dengue-fs.ihme.washington.edu\ubcov\ubcov_central"
}

clear all
set more off
set obs 1

///////////////////////////////////////////////////////////////////////
/*  Arguments:
		- topic: your research topic
        - ubcov_id: The id of the codebook row
		- outpath_L: output file path for limited drive files
			(GBD: <<<< FILEPATH REDACTED >>>>)
			(LBD: <<<< FILEPATH REDACTED >>>>)
		- outpath_J: output file path for general files
    Optional:
        - keep: Keeps both raw data and extracted data, allows manual extraction check before final output.
        - bypass: Skips the extraction check, output the data 
        - run_all: Loops through all ubcov_ids in the codebook.
*/
////////////////////////////////////////////////////////////////////////

local topics male circumcision
local array 16358
local outpath_L "`l'<<<< FILEPATH REDACTED >>>>"
//local outpath_J "`j'<<<< FILEPATH REDACTED >>>>"
local options bypass //leave it blank if you don't want to keep, bypass, or run_all. If you are running on cluster, you have to use bypass as an option otherwise it will error.



///////////////////////////////////////////////////////////////////////
//Extraction
//////////////////////////////////////////////////////////////////////


// Load functions
cd "`central_root'"
do "`central_root'/modules/extract/core/load.do"

// Load the base code for ubCov
local paths  `central_root'/modules/extract/core/ `central_root'/modules/extract/core/addons/
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
    tostring year_start, gen(start_year)
    tostring year_end, gen(end_year)
    tostring nid, gen(nid_n)
    local filename = survey_name + "_" + nid_n + "_" + survey_module + "_" + ihme_loc_id + "_" + start_year + "_" + end_year
    local filename = subinstr("`filename'", "/", "_",.)
    drop year_start year_end nid_n
	if (strpos("$file_path", "LIMITED_USE")|strpos("$file_path", "IDENT")){
		local outpath = "`outpath_L'"
	}
	else{
		local outpath = "`outpath_J'"
	}
	cd  `outpath'
    save "`filename'", replace
}
