"""Save population xarrays into CSV and RDS

This script will take in the specified population run version from the FHS outputs
in `"<<<< FILEPATH REDACTED >>>>"`, create means from the draws and
save out a CSV and RDS file in `"<<<< FILEPATH REDACTED >>>>"`

The only argument is `pop_version` corresponding to the folder names in the FHS
population folder mentioned above.

Necessary Python packages (with versions used when developing this script):
`xarray`: 0.12.1
`pandas`: 0.24.2
`rpy2`: 2.9.5*

NOTE*: It is very important to use this version of rpy2 and nothing above; rpy2's
API changed a lot from version 3 and will NOT work with this script.
Use `pip install rpy2==2.9.5` to ensure the correct version of rpy2 is installed
in your (conda) environment.
"""
import xarray as xr
import pandas as pd
import argparse
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()


def get_xarray_mean(pop_version):
    """ Loads the population xarray from "<<<< FILEPATH REDACTED >>>>"
    particularly the population_agg file, and creates the mean value across draws

    Parameters
    ----------
    pop_version : str
        The version of population from FHS outputs

    Returns
    -------
    pd.DataFrame
        A pandas dataframe with mean population
    """
    filepath = "<<<< FILEPATH REDACTED >>>>"
    print(f"Getting population version: {pop_version}")
    print(f"Filepath: {filepath}")
    xar_open = xr.open_dataset(f"{filepath}")
    print(xar_open)
    print("Get mean of draws")
    xar_open_mean = xar_open.mean('draw')
    print("Converting to pd dataframe and prepping for save")
    xar_open_mean = xar_open_mean.to_dataframe().reset_index()
    xar_open_mean.rename(columns={"value": "population"}, inplace=True)
    return xar_open_mean


def save_out(pop_version, mean_file):
    """Saves out a CSV and RDS file with the mean population from `get_xarray_mean()`

    Parameters
    ----------
    pop_version : str
        The version of population from FHS outputs
    mean_file : pd.DataFrame
        The pandas dataframe with population (output of `get_xarray_mean()`)

    Returns
    -------
    int
       0 if everything is good
    """

    print("Saving csv...")
    mean_file.to_csv("<<<< FILEPATH REDACTED >>>>")
    print("Move the dataframe into R and save RDS...")
    r_dataframe = pandas2ri.py2ri(mean_file)
    ro.globalenv['population'] = r_dataframe
    ro.r('population <- data.table::data.table(population)')
    ro.r(f"saveRDS(population, file = "
         f"\"<<<< FILEPATH REDACTED >>>>"")")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--pop_version',
        type=str,
        help="The population version by FHS"
    )
    args = parser.parse_args()
    mean_pop = get_xarray_mean(pop_version=args.pop_version)
    save_out(pop_version=args.pop_version, mean_file=mean_pop)


if __name__ == '__main__':
    main()
