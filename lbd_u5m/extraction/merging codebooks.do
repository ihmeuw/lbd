/////////////////////////////////////////////////////////////
/// Combined all geography codebooks into one DTA file //////
/// Uses an R script - merging codebooks.R first then 
/// finishes the cleaning in stata

// Set working directory to where you want to save the combined codebook
cd "<<<< FILEPATH REDACTED >>>>"

use "combined_codebook.dta"
//convert labels vars to string
foreach var in iso3 cluster_number location_name location_code shapefile admin_level uncertain_point buffer{
	decode `var', gen(temp)
	drop `var'
	rename temp `var'
	}

order nid iso3 cluster_number cluster_number point latnum longnum uncertain_point buffer location_name location_code admin_level shapefile
duplicates drop nid iso3 cluster_number, force
save "combined_codebook.dta", replace
