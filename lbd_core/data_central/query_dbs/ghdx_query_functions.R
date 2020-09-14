################################################################################
##
## TOOLS FOR QUERYING THE GHDX
##
## Created: Sept 11, 2018
## Purpose: Functions for querying the GHDx to get information about source
##   availability metadata.
##
################################################################################



## ghdx_strip_formatting ---------------------------------------------------->
#'
#' @title Strip HTML formatting from text
#' @description Strip formatting from queried text from the GHDx
#'
#' @details This function takes a character vector as input and removes any
#'   HTML-style formatting tags from each character in the vector. The output
#'   of this function should be plain text.
#'
#' @param formatted_text Character vector of text that includes HTML formatting
#'
#' @return Character vector stripped of formatting tags
#'
ghdx_strip_formatting <- function(formatted_text){
  # The following matches all HTML text tags from <b> to </strong> and
  #  strips them from the text
  library(XML)
  stripped <- sapply(
    formatted_text,
    function(x){
      doc <- htmlParse(x, asText=TRUE)
      plain_text <- xpathSApply(
        doc,
        paste0("//text()[not(ancestor::script)][not(ancestor::style)]",
               "[not(ancestor::noscript)][not(ancestor::form)]"),
        xmlValue
      )
      return(paste(plain_text, collapse=""))
    }
  )
  return(unname(stripped))
}



## ghdx_execute_query --------------------------------------------------------->
#'
#' @title Execute a GHDx query
#' @description Query the GHDx about a set of NIDs given a query template
#'
#' @details Query the GHDx about a set of NIDs given a query template.
#'   \strong{This function requires that the user has set up a valid ODBC file
#'   in their H:/ drive.}
#'
#' @param nids Vector of all NIDs to query
#' @param query_string_template Character vector containing the full SQL query
#'   containing the section "WHERE <nid field> IN (%s);". This template will
#'   be replaced with the actual NIDs and used to run the query
#'
#' @return a data.table of results from the GHDx query
#'
ghdx_execute_query <- function(nids, query_string_template){
  ## Chunk NIDs vector into a list, where each sub-vector contains no more than
  ##  100 NIDs
  nid_chunks <- split(nids, ceiling(seq_along(nids)/100) )
  ## For each vector of 100 NIDs, execute the query and return as a data.table
  result_chunks <- lapply(
    nid_chunks,
    function(nids_sub){
      # Define the full query string
      query_string <- sprintf(
        query_string_template,
        paste(nids_sub, collapse=',')
      )
      # Execute query
      query_results_sub <- execute_sql(
        statement=query_string,
        conn_def='ghdx',
        odbc_filepath='~/.odbc.ini',
        return_output=TRUE
      )
      return(query_results_sub)
    }
  )
  ## Combine the results into a single data.table and return
  full_results <- rbindlist(result_chunks)
  return(full_results)
}


## ghdx_get_all_nid_citations ------------------------------------------------->
#'
#' @title Get NID citations
#' @description Pull citation metadata from the GHDx, INCLUDING CONFIDENTIAL NIDs
#'
#' @details Given a vector of NIDs (the unique "Node IDentifiers" representing
#'   a source in the GHDx), pull citation and confidentiality metadata about
#'   each identified source. NIDs not contained within the GHDx will have
#'   metadata assigned as "missing".
#'
#'
#' @param nids A vector of GHDx NIDs
#' @param core_repo [default="<<<< FILEPATH REDACTED >>>>"] Path to
#'   the core repository, used to load database reading packages.
#' @param strip_html_tags Should HTML tags such as '<i>' and '<p>' be removed?
#'
#'
#' @return A data.table of formatted NIDs
#'
#' @seealso \code{\link{ghdx_execute_query}}
#'
ghdx_get_all_nid_citations <- function(
  nids,
  core_repo="<<<< FILEPATH REDACTED >>>>",
  strip_html_tags=TRUE
){
  # Import needed functions
  supp_source <- function(fp) suppressMessages(source(fp))
  supp_source(paste0(core_repo,'/mbg_central/setup.R'))
  suppressMessages(load_R_packages(c('data.table','XML')))
  supp_source(paste0(core_repo,'/data_central/query_dbs/db_base_functions.R'))
  # Define SQL query template
  query_string_template <- (
    "
    SELECT
    node.nid as nid,
    fdfc.field_citation_value as citation,
    fdfc.field_citation_format as citation_format,
    ttd.name as confidentiality,
    ttd.description as confidentiality_details
    FROM
    ghdx.node node
    LEFT JOIN
    ghdx.field_data_field_private_data fdfpd
    ON node.nid = fdfpd.entity_id
    LEFT JOIN
    ghdx.taxonomy_term_data ttd
    ON fdfpd.field_private_data_tid = ttd.tid
    LEFT JOIN
    ghdx.field_data_field_citation fdfc
    ON fdfpd.entity_id = fdfc.entity_id
    WHERE fdfpd.entity_id IN (%s);
    "
  )
  # Get results from the query
  citations_all <- ghdx_execute_query(
    nids = unique(nids),
    query_string_template = query_string_template
  )
  if(nrow(citations_all)==0){
    # Add columns if none exist
    citations_all <- data.table(
      nid=integer(0),
      citation = character(0),
      citation_format = character(0),
      confidentiality = character(0),
      confidentiality_details = character(0)
    )
  }
  citations_all[is.na(citation), citation:="**CITATION MISSING**" ]
  # If specified, strip HTML tags
  if (strip_html_tags){
    citations_all$citation <- ghdx_strip_formatting(citations_all$citation)
  }
  # Return filled results
  return(citations_all)
}


## get_public_nid_citations --------------------------------------------------->
#'
#' @title Get Public NID citations
#' @description Pull public-facing citation metadata from the GHDx
#'
#' @details Not all NIDs in the GHDx should be cited publicly; some NIDs should
#'   not be cited at all, while other NIDs should have information from their
#'   parent NID cited if it exists. This function wraps the function
#'   \code{get_all_nid_citations} to iteratively pull citation metadata for a
#'   group of NIDs, then handle that data appropriately by either keeping it,
#'   dropping it (for internal-only NIDs), or citing parent NID data if it
#'   exists.
#'
#'
#' @param nids A vector of GHDx NIDs
#' @param core_repo [default="<<<< FILEPATH REDACTED >>>>"] Path to
#'   the core repository, used to load database reading packages.
#'
#' @return A data.table of NIDs with public-facing citation information
#'
#' @seealso \code{\link{ghdx_get_all_nid_citations}}
#'
get_public_nid_citations <- function(
  nids,
  core_repo="<<<< FILEPATH REDACTED >>>>"
  ){
  # Import needed functions
  supp_source <- function(fp) suppressMessages(source(fp))
  supp_source(paste0(core_repo,'/mbg_central/setup.R'))
  suppressMessages(load_R_packages(c('data.table')))
  supp_source(paste0(core_repo,'/data_central/query_dbs/db_base_functions.R'))

  # Define vector determining which sources should not be cited
  drop_types <- c('Research','Composite source, do not cite','Private')

  # Sub-function for pulling NIDs and dropping missing/private sources
  pull_sub_citations <- function(nids_sub){
    # Pull a single set of citations
    citations <- ghdx_get_all_nid_citations(nids=nids_sub, core_repo=core_repo)
    # Drop any private sources
    citations <- citations[!(confidentiality %in% drop_types),]
    # Return sub-dataset
    return(citations)
  }

  ## Pull citations the first time
  working_citations <- pull_sub_citations(nids)
  # Create a field, 'no_parents', which indicates that a citation definitively
  #  does not have a parent ID
  working_citations[, no_parents := 0]

  ## This function runs iteratively in a "while" loop - while the set of NIDs
  ##  still needs any parent IDs checked, keep running
  parent_iterations <- 1
  while( any(working_citations[,no_parents]==0) & parent_iterations <= 5 ){
    # Check for parent NIDs
    parent_query_template <- (
      "
      SELECT
        entity_id as nid,
        field_parent_projects_target_id as parent_id
      FROM
        ghdx.field_data_field_parent_projects fdfp
      WHERE
        fdfp.entity_type = 'node'
        AND
        fdfp.entity_id IN (%s)
      "
    )
    parent_ids_dt <- ghdx_execute_query(
      nids = working_citations[ no_parents==0, nid ],
      query_string_template = parent_query_template
    )
    # Keep only data for NIDs without parents
    nids_with_parents <- unique(parent_ids_dt[,nid])
    without_parents <- working_citations[!(nid %in% nids_with_parents),]
    without_parents[, no_parents := 1]
    # Query citation data for the parent NIDs and add it to the dataset
    parent_citations <- pull_sub_citations( unique(parent_ids_dt[,parent_id]) )
    parent_citations[, no_parents := 0]
    working_citations <- rbindlist( list(without_parents, parent_citations) )
    # Iterate the number of times parents have been searched
    parent_iterations <- parent_iterations + 1
  }
  if(parent_iterations > 5){
    warning(paste0("There is an issue with the parent IDs of the NIDs you are ",
                   "pulling; maximum parent iterations exceeded."))
  }
  # All NIDs should now have citation information
  # All parents have also been checked
  working_citations[, c('no_parents','confidentiality',
                        'confidentiality_details') := NULL ]
  return(working_citations)
}


## ghdx_resolve_urls ---------------------------------------------------------->
#'
#' @title Resolve GHDx URL paths
#' @description Find public GHDx URLs for a set of NIDs
#'
#' @details Users can access GHDx records using the URL shortcut
#'   `https://ghdx.healthdata.org/node/<NID>`. However, this format is not
#'   appropriate for usage in publication-ready citations. This function
#'   programmatically follows links to GHDx websites using the `curl` function,
#'   returning the resolved links matched to each NID. If an NID does not
#'   match to a resolved link, it will be matched to `NA`.
#'
#' @param nids A vector of GHDx NIDs
#'
#' @return A data.frame of NIDs with matching GHDx URLs
#'
ghdx_resolve_urls <- function(nids){
  # Attempt to resolve the URL for each NID
  url_vector <- sapply(nids, function(n) {
    curl_cmd <- paste0("curl -Ls -o /dev/null -w %{url_effective} ",
                       "http://ghdx.healthdata.org/node/", n, "/")
    curl_output <- system(curl_cmd, intern = T)
    if (grepl("/node", curl_output)) curl_output <- NA
    return(curl_output)
  })
  # Convert to a data frame
  url_dt <- data.frame(
    nid      = nids,
    ghdx_url = url_vector
  )
  return(url_dt)
}


## ghdx_locations_to_iso ------------------------------------------------------>
#'
#' @title Correct GHDx locations
#' @description Correct GHDx locations
#'
#' @details The GHDx does not always return a standard set of locations. This
#'   function takes a character vector containing location names and matches
#'   them to standard ISO3 codes by country
#'
#' @param locations_dt Data.table containing the field 'location_name'
#'
#' @return Standardized character vector of ISO3 codes
#'
ghdx_locations_to_iso <- function(locations_dt){
  ## Location name cleaning
  ## Function to strip special characters from a vector of location names
  strip_special_chars <- function(vec){
    vec <- iconv(vec, to="ASCII//TRANSLIT")
    vec <- gsub("\\'", '', vec)
    return(vec)
  }

  ## This function fixes common location-name to ISO code errors
  apply_iso_fixes <- function(locs_df){
    locs_df[ is.na(location_name), iso:='Unknown' ]
    locs_df[ iso=="NA", iso := NA ]
    # Location-specific fixes
    start_fix_list = list(
      'Sudan'='SDN',
      '710'  ='TWN'
    )
    end_fix_list = list(
      'Ivoire'    ='CIV',
      'Kosovo'    ='KSV',
      'Zaire'     ='COD',
      'Equatoria' ='SUD'

    )
    contains_fix_list = list(
      'Yugoslavia'              ='SRB',
      'Serbia'                  ='SRB',
      'San Marino'              ='SMR',
      'Palau'                   ='PLW',
      "South Sudan"             ='SSD',
      "Serbia"                  ='SRB',
      "United States"           ='USA',
      "Niue"                    ='NIU',
      "Nauru"                   ='NRU',
      "Saint Kitts"             ='KNA',
      'Greenland'               ='GRL',
      'Bermuda'                 ='BMU',
      'Puerto Rico'             ='PRI',
      'Virgin Islands, U.S.'    ='VIR',
      'American Samoa'          ='ASM',
      'Guam'                    ='GUM',
      'Northern Mariana Islands'='MNP',
      'Al Gharbiyah'            ='EGY',
      'Al Minya'                ='EGY',
      'Al Qalyubiyah'           ='EGY',
      'Asyut'                   ='EGY',
      'Qina'                    ='EGY',
      'Suhaj'                   ='EGY',
      'Accra'                   ='GHA',
      'Bungoma'                 ='KEN',
      'Garissa'                 ='KEN',
      'Kakamega'                ='KEN',
      'Mandera'                 ='KEN',
      'Nyanza'                  ='KEN',
      'Turkana'                 ='KEN',
      'Wajir'                   ='KEN',
      'Fianarantsoa'            ='MDG',
      'Toliara'                 ='MDG',
      'Bari'                    ='SOM',
      'Mudug'                   ='SOM',
      'Nugaal'                  ='SOM',
      'Somaliland'              ='SOM'
    )

    ## Make the replacements
    for(k in names(start_fix_list)){
      locs_df[ startsWith(location_name,k), iso:=start_fix_list[[k]] ]
    }
    for(k in names(end_fix_list)){
      locs_df[ endsWith(location_name,k), iso:=end_fix_list[[k]] ]
    }
    for(k in names(contains_fix_list)){
      locs_df[ grepl(k, location_name), iso:=contains_fix_list[[k]]]
    }
    return(locs_df)
  }

  ## This function pulls all vector-ISO code combinations as a data.table from
  ##  the GHDx
  get_iso_map <- function(){
    query_str = "SELECT location_name, map_id as iso FROM shared.location;"
    loc_table <- execute_sql(statement = query_str, conn_def  = 'ghdx')
    loc_table <- loc_table[ !is.na(iso) & !is.na(location_name), ]
    ## CLEANING
    # - Strip special character from location names
    # - Keep only the first three characters from each ISO code
    # - Apply fixes for ISO codes
    # - Drop USA locations (to avoid Georgia mixup)
    loc_table[, location_name := strip_special_chars(location_name) ]
    loc_table[, iso := substr(iso, 1, 3) ]
    loc_table <- apply_iso_fixes(loc_table)
    loc_table <- unique(loc_table[ !(iso %in% c('USA','Non')), ])
    # Return cleaned data
    return(loc_table)
  }

  ## This function converts ISO codes back to countries
  get_iso_to_country <- function(){
    query_str <- "SELECT location_name as country, map_id as iso FROM shared.location;"
    loc_table <- execute_sql(statement = query_str, conn_def  = 'ghdx')
    loc_table <- loc_table[ !is.na(iso) & !is.na(country), ]
    # Keep only Iso codes with three digits
    loc_table <- loc_table[ nchar(iso)==3, ]
    # Basic corrections for some countries
    additions <- data.table(
      country = c("Taiwan"),
      iso     = c("TWN")
    )
    loc_table <- rbindlist(list(loc_table, additions))
    return(loc_table)
  }

  ## ------------------------- MAIN FUNCTION EXECUTION -------------------------

  # Create a data.table of the locations
  locations_dt[, location_name := strip_special_chars(location_name)]
  # Drop the ISO field if it already exists
  if('iso' %in% names(locations_dt)){
    locations_dt[,iso:=NULL]
  }
  # Map to ISO codes
  iso_map <- get_iso_map()
  locs_joined <- merge(
    x     = locations_dt,
    y     = iso_map,
    by    = c('location_name'),
    all.x = T
  )
  if(nrow(locs_joined) != nrow(locations_dt)){
    stop("There is an issue of non-unique location name <-> ISO translation.")
  }
  # Fix some ISO codes by location names
  locs_joined <- apply_iso_fixes(locs_joined)
  # Add on country names
  locs_joined <- merge(
    x     = locs_joined,
    y     = get_iso_to_country(),
    by    = c('iso'),
    all.X = TRUE
  )
  # Return the vector of ISO codes
  return(locs_joined)
}


## ghdx_construct_pub_table --------------------------------------------------->
#'
#' @title Construct source publication table
#' @description Pull metadata of GHDx source for a LBD publication appendix
#'
#' @details Pull relevant GHDx metadata that could be included in the source
#'   citations for a Local Burden of Disease publication.
#'
#'
#' @param nids A vector of GHDx NIDs
#' @param core_repo [default="<<<< FILEPATH REDACTED >>>>"] Path to
#'   the core repository, used to load database reading packages.
#' @param resolve_urls [default=TRUE] Should GHDx links to each source be pulled?
#'   (Specifying `TRUE` in this argument will slow down the function a bit)
#'
#' @return A data.table of NIDs with public-facing citation information
#'
#' @seealso \code{\link{get_all_nid_citations}}, \code{\link{ghdx_resolve_urls}}
#'
ghdx_construct_pub_table <- function(
  nids,
  core_repo="<<<< FILEPATH REDACTED >>>>",
  resolve_urls=TRUE
  ){
  ## Import needed functions
  supp_source <- function(fp) suppressMessages(source(fp))
  supp_source(paste0(core_repo,'/mbg_central/setup.R'))
  suppressMessages(load_R_packages(c('data.table')))
  supp_source(paste0(core_repo,'/data_central/query_dbs/db_base_functions.R'))

  ## Pull metadata
  # Pull all public citations for the given set of NIDs
  base_citations <- get_public_nid_citations(nids=nids, core_repo=core_repo)
  # Define query to pull additional metadata:
  pub_meta_template <- (
    "
    SELECT
      node.nid   as nid,
      node.title as title,
      ttd1.name  as location_name,
      ttd2.name  as geoprecision,
      SUBSTRING(fdft.field_time_value, 1, 4)  as start_year,
      SUBSTRING(fdft.field_time_value2, 1, 4) as end_year
    FROM
    	 ghdx.node node
    LEFT JOIN
    	ghdx.field_data_field_geography fdfg
      ON node.nid = fdfg.entity_id
    LEFT JOIN
    	ghdx.taxonomy_term_data ttd1
      ON fdfg.field_geography_tid = ttd1.tid
    LEFT JOIN
    	ghdx.field_data_field_time fdft
      ON node.nid = fdft.entity_id
    LEFT JOIN
    	ghdx.field_data_field_geoprecision fdfgp
      ON node.nid = fdfgp.entity_id
    LEFT JOIN
    	ghdx.taxonomy_term_data ttd2
    	ON fdfgp.field_geoprecision_target_id = ttd2.tid
    WHERE (fdfg.delta=0) AND (node.nid IN (%s));
    "
  )
  # Pull additional metadata
  pub_metadata <- ghdx_execute_query(
    nids = unique(base_citations[,nid]),
    query_string_template = pub_meta_template
  )
  # If no rows were pulled, create a data.table with the appropriate columns
  if(nrow(pub_metadata)==0){
    pub_metadata <- data.table(
      nid           = integer(0),
      title         = character(0),
      location_name = character(0),
      geoprecision  = character(0),
      start_year    = integer(0),
      end_year      = integer(0)
    )
  }
  # Add a field for ISO code; keep only unique survey metadata by ISO
  with_isos <- ghdx_locations_to_iso(copy(pub_metadata))
  with_isos[, location_name := NULL]
  with_isos <- unique(with_isos)

  # Function to keep only the most granular spatial data
  get_most_granular_gp <- function(gp_group){
    gp_group <- tolower(gp_group)
    gran_ranking <- c('Latitude/longitude','Precise place names','Admin4',
                      'Admin3','Admin2','Admin1','Country')
    for(gr in gran_ranking){
      if(tolower(gr) %in% gp_group){
        return(gr)
      }
    }
    # There was no match
    return("Unknown")
  }
  deduped <- with_isos[, .(geoprecision=get_most_granular_gp(geoprecision)),
                        by = .(nid, title, iso, country, start_year, end_year)]
  ## Merge citations, GHDx metadata, and URL
  all_merged <- merge(
    x     = base_citations,
    y     = deduped,
    by    = c('nid'),
    all.x = TRUE
  )

  # If specified, pull link information for each NID
  if(resolve_urls==TRUE){
    link_metadata <- ghdx_resolve_urls( nids=unique(base_citations[,nid]) )
    all_merged <- merge(
      x     = all_merged,
      y     = link_metadata,
      by    = c('nid'),
      all.x = TRUE
    )
  }

  ## Clean up merged dataset, then return
  all_merged[, years := ifelse(start_year==end_year,
                               as.character(start_year),
                               paste0(start_year,'-',end_year))]
  all_merged[, c('start_year','end_year','citation_format') := NULL]
  # Remove improper geoprecision levels
  improper_gp <- c('Country','None','Unknown','Not applicable')
  all_merged[ (geoprecision %in% improper_gp), geoprecision:=NA ]
  all_merged <- all_merged[order(iso, years)]
  return(all_merged)
}
