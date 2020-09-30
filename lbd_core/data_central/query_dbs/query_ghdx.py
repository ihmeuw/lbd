'''
FUNCTIONS FOR QUERYING THE GHDx

Created: August 27, 2018
Purpose: Functions for querying the GHDx for metadata about source NIDs.
  Examples of useful data stored in the GHDx include citations, start and end
  years, countries, and data types related to each NID.

NOTE: These functions must be run in a Python environment (preferably Anaconda)
  with the GBD Shared Functions installed.

'''

# Imports (general)
import numpy as np
import pandas as pd
import os
import unidecode
from datetime import datetime

# Imports (GBD)
from db_queries import get_location_metadata, get_covariate_estimates
from db_tools.ezfuncs import query


################################################################################
## FUNCTIONS: DIRECT QUERIES OF THE GHDX
################################################################################

def get_nids_natl_representative():
    '''
    Get a dataframe of NIDs for all nationally-representative surveys in the GHDx
    '''
    # Define query
    q = '''
        SELECT
            representative.nid
        FROM
            (SELECT
                DISTINCT(entity_id) as nid
            FROM
                ghdx.field_data_field_secondary
            WHERE field_secondary_tid IN (
                SELECT
                    tid
                FROM
                    ghdx.ihme_terms_synonyms
                WHERE
                    name='Nationally representative')
                ) representative
        INNER JOIN
            (SELECT
                DISTINCT(entity_id) as nid
            FROM
                ghdx.field_data_field_type fdft
            INNER JOIN (
                SELECT
                    tid, name
                FROM
                    ghdx.ihme_terms_synonyms
                WHERE
                    name LIKE "%%survey%%"
                ) its
                ON fdft.field_type_tid = its.tid) surveys
            ON surveys.nid = representative.nid;
    '''
    # Query GHDx database
    natl_rep = query(q, conn_def='ghdx')
    # Return all nationally-representative surveys by NID
    return natl_rep


def get_chunks(l, n):
    """
    Helper function to chunk data.
    Return a list of successive n-sized chunks from l.
    """
    out_list = list()
    for i in range(0, len(l), n):
        out_list.append(l[i:i+n])
    return out_list


def get_survey_start_end_year(nids_list):
    '''
    Given a list of NIDs, get the start and end year for each.
    '''
    # Split list of NIDs into chunks of 100
    chunks = get_chunks(nids_list, 100)
    # For each chunk, get associated time info
    chunked_results = list()
    for chunk in chunks:
        # Define the search query
        q = '''
        SELECT
            entity_id as nid,
            field_time_value as start_date,
            field_time_value2 as end_date
        FROM
            ghdx.field_data_field_time
        WHERE
            entity_id IN ({nids});
        '''.format( nids=','.join([str(i) for i in chunk]) )
        # Query the DB and append to list
        chunked_result = query(q, conn_def='ghdx')
        chunked_results.append(chunked_result)
    # Concatenate all results into a single dataframe
    all_results = pd.concat(chunked_results)
    # FORMATTING
    all_results['start_year'] = (
        all_results['start_date'].apply(lambda x: str(x)[0:4]).astype(np.int32)
    )
    all_results['end_year'] = (
        all_results['end_date'].apply(lambda x: str(x)[0:4]).astype(np.int32)
    )
    all_results = all_results.drop(labels=['start_date','end_date'], axis=1)
    return all_results


def get_location_name(nids_list):
    '''
    Given a list of NIDs, get the location name associated with each. In cases where
    more than one location name is associated with a single NID, keep only the first
    location name.
    '''
    # Split list into chunks of 100
    chunks = get_chunks(nids_list, 100)
    # For each chunk, get associated location info
    chunked_results = list()
    for chunk in chunks:
        # Define search query
        q = '''
        SELECT
            field_data_field_geography.entity_id as nid,
            ihme_terms_synonyms.name as location_name
        FROM
            ghdx.field_data_field_geography
        LEFT JOIN
            ghdx.ihme_terms_synonyms
            ON field_data_field_geography.field_geography_tid = ihme_terms_synonyms.tid
        WHERE
            field_data_field_geography.entity_id IN ({nids});
        '''.format( nids=','.join([str(i) for i in chunk]) )
        # Query the DB and append to list
        chunked_result = query(q, conn_def='ghdx')
        chunked_results.append(chunked_result)
    # Concatenate all results into a single dataframe
    all_results = pd.concat(chunked_results)
    # In cases where multiple locations are associated with one NID, keep only
    #  the first
    all_results = all_results.drop_duplicates(subset=['nid'], keep='first')
    # Ensure the proper column formatting for NIDs
    all_results['nid'] = all_results['nid'].astype(np.int64)
    return all_results



################################################################################
## FUNCTIONS FOR CLEANING GHDx LOCATION NAMES AND MATCHING TO ISO3 CODES
################################################################################

def clean_loc_names(loc_series):
    '''
    Given a series of location names, convert all special characters to ASCII
    '''
    # Remove special characters
    loc_series = loc_series.apply(
        lambda x: unidecode.unidecode(str(x))
    )
    return loc_series


def apply_iso_fixes(isos_df):
    '''
    Takes in a dataFrame with location names and ISO3 codes
    Fixes some known issues with missing ISO codes based on location name
    '''
    # Clean location names to avoid missingness
    isos_df.loc[isos_df['location_name'].isnull(), 'location_name'] = 'Unknown'
    # Make sure that no NaNs are coded as strings
    isos_df.loc[isos_df['iso']=="NaN",'iso'] = np.nan
    # Location-specific fixes
    start_fix_dict = {'Sudan': 'SDN',
                      '710'  : 'TWN'
                      }
    end_fix_dict   = {'Ivoire': 'CIV',
                      'Kosovo': 'KSV',
                      'Zaire' : 'COD'
                      }
    contains_fix_dict = {'Yugoslavia'              : 'SRB',
                         'Serbia'                  : 'SRB',
                         'San Marino'              : 'SMR',
                         'Palau'                   : 'PLW',
                         "South Sudan"             : 'SSD',
                         "Serbia"                  : 'SRB',
                         "United States"           : 'USA',
                         "Niue"                    : 'NIU',
                         "Nauru"                   : 'NRU',
                         "Saint Kitts"             : 'KNA',
                         'Greenland'               : 'GRL',
                         'Bermuda'                 : 'BMU',
                         'Puerto Rico'             : 'PRI',
                         'Virgin Islands, U.S.'    : 'VIR',
                         'American Samoa'          : 'ASM',
                         'Guam'                    : 'GUM',
                         'Northern Mariana Islands': 'MNP'
                         }
    # Make the replacements
    for k, v in start_fix_dict.items():
        isos_df.loc[isos_df['location_name'].str.startswith(k),'iso'] = v
    for k, v in end_fix_dict.items():
        isos_df.loc[isos_df['location_name'].str.endswith(k),'iso'] = v
    for k, v in contains_fix_dict.items():
        isos_df.loc[isos_df['location_name'].str.contains(k),'iso'] = v
    # Fill all missing codes with '999' for now
    isos_df.loc[isos_df['iso'].isnull(),'iso'] = '999'
    # Taiwan ISO code fix
    isos_df.loc[isos_df['iso'].str.startswith('710'),'iso']='TWN'
    # Revert to NaNs
    isos_df.loc[isos_df['iso']=='999','iso'] = np.nan
    return isos_df


def get_iso_map():
    '''
    Create a dataframe that maps all location names in the GBD database to ISO3
    codes.
    '''
    # Query the location table
    q = '''
    SELECT location_name, map_id
    FROM shared.location;
    '''
    locs = query(q, conn_def='ghdx')
    # Define ISO field
    locs = locs.loc[ locs['map_id'].notnull(), : ]
    locs['iso'] = locs['map_id'].apply(lambda x: str(x)[:3])
    # Clean location names to remove special characters
    locs['location_name'] = clean_loc_names(locs['location_name'])
    # Drop USA locations (to avoid Georgia mixup)
    # Get unique combinations of location names to ISOs
    locs = apply_iso_fixes(locs)
    locs = (locs.loc[(locs['iso'] != 'USA') & (locs['iso'] != 'Non'),
                     ['location_name','iso']]
                .drop_duplicates())
    return(locs)


def full_nids_to_isos(nids_list):
    '''
    Combination of the functions above: given a list of NIDs, return a cleaned
    DF containing each NID and the corresponding ISO code.
    '''
    # Get location names and clean
    loc_names = get_location_name(nids_list)
    loc_names['location_name'] = clean_loc_names(loc_names['location_name'])
    # Match to ISO codes
    iso_map = get_iso_map()
    isos = pd.merge(
        left  = loc_names,
        right = iso_map,
        on    = ['location_name'],
        how   = 'left'
    )
    # Fix some ISO codes by location names
    isos = apply_iso_fixes(isos)
    return isos.loc[:,['nid','iso']]
