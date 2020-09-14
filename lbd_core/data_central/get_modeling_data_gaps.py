'''
GET DATA GAPS FROM MODELED DATA

Created: August 27, 2018
Purpose: Show the number of surveys available to a project by year, country, and
  types of data available.

NOTE: This code must be run in a Python environment (preferably Anaconda) with
  the GBD shared functions installed. For more information, see the forthcoming
  Hub documentation.

'''

# Imports (general)
import numpy as np
import os
import pandas as pd
import re
import sys
# Imports (shared functions)
from db_queries import get_population
# Imports (project-specific)
from query_dbs import query_ghdx as ghdx


################################################################################
# FUNCTIONS TO GENERATE DATA GAPS INFORMATION FROM MODEL INPUT DATA
################################################################################

def iso_to_region(iso_df, merge_col='mbg_reg'):
    '''
    Given a df that has an 'iso' field, add a 'region' field containing standard
    MBG region.
    '''
    # Read in the ISO lookup table
    country_table_file = '/home/j/WORK/11_geospatial/10_mbg/stage_master_list.csv'
    country_table = pd.read_csv(country_table_file)
    # Subset to columns to be added
    country_table_sub = (
        country_table.loc[:,['iso3','loc_id',merge_col]]
                     .rename(columns={
                        'iso3'   : 'iso',
                        'loc_id' : 'location_id'
                        })
    )
    # Fill missing regions with "other"
    country_table_sub[merge_col] = country_table_sub[merge_col].fillna('other')
    # Merge the lookup table onto the main dataframe
    merged = pd.merge(
        left  = iso_df,
        right = country_table_sub,
        on    = ['iso'],
        how   = 'left'
    )
    # Get TOTAL population of country from GBD database
    unique_locations = merged.loc[merged['location_id'].notnull(),
                                  'location_id'].unique().tolist()
    pops = get_population(
        gbd_round_id = 5, # GBD 2017
        location_id  = unique_locations, # Locations in DF
        year_id      = 2017,
        sex_id       = 3, # Both sexes combined
        age_group_id = 22 # All ages combined: 22, under-5s only: 1
    ).loc[:,['location_id','population']]
    merged = pd.merge(
        left  = merged,
        right = pops,
        on    = ['location_id'],
        how   = 'left'
    )
    return merged


def process_mbg_in_data(in_df):
    '''
    This function takes a dataframe, reduces to one row per NID, and then
    merges on ISO and region
    '''
    # If 'nid' is not in the columns, rename it
    if 'nid' not in in_df.columns:
        in_df = in_df.rename(columns={'svy_id':'nid'})
    # Keep only necessary fields
    in_df = in_df.loc[:,['nid','year','point']]
    ## STANDARDIZE YEAR BINS
    # Drop anything that wouldn't be binned into 2000, 2005, 2010, or 2015
    in_df = in_df.loc[(in_df['year'] > 1997) & (in_df['year'] < 2018),:]
    #  Round to the nearest five years
    in_df['year'] = np.round(in_df['year']/5,0) * 5
    # If any 'point' field is 1 by NID, then it's point data
    in_df = in_df.groupby(by=['nid','year']).max().reset_index()
    # Keep only the first row from each NID
    by_nid = in_df.drop_duplicates(subset=['nid'], keep='first')
    # Drop rows where NID is missing
    by_nid = by_nid.loc[by_nid['nid'].notnull(),:]
    # Get a list of ISO codes for all unique NIDs in the data
    # Merge back onto the main data
    iso_merge_table = ghdx.full_nids_to_isos( in_df['nid'].unique().tolist() )
    nids_with_isos = pd.merge(
        left  = by_nid,
        right = iso_merge_table,
        on    = ['nid'],
        how   = 'left'
    )
    # Check that there are no missing ISO codes
    missing_isos = (
        nids_with_isos.loc[nids_with_isos['iso'].isnull(),'nid'].unique().tolist()
    )
    if len(missing_isos) > 0:
        print("The following NIDs have missing ISO codes: {}".format(
            ','.join(missing_isos)
        ))
    # Add MBG regions to the standard ISO table
    with_regions = iso_to_region(nids_with_isos)
    # Return the formatted data
    return with_regions


def mbg_modeling_data_gaps(indicator, out_file):
    '''
    Given an MBG indicator, read in the modeling input data from a standard folder,
    process for data availability plots, and then write to file.
    '''
    # Read file from standard location
    save_loc = "<<<< FILEPATH REDACTED >>>>"
    in_data = pd.read_csv( save_loc )
    # Process the file for data availability
    processed = process_mbg_in_data(in_data)
    # Save to output file
    processed.to_csv(out_file, index=False, encoding='utf8')
    return None
