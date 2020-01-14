# Link Table Generation in Parallel

## Overview

This script generates the link table in parallel. Most of the parameters are set in `launch.sh`.

The following pattern is followed:

1. Partition the shapefile by ADM0 geometry area.
2. Build a link table for each partition.
3. Combine the partitioned link tables together.
4. Cleanup/ensure script ran successfully.

## Use

    ./launch.sh <directory>

Where `directory` is a child directory within the base admin shapefile directory.

Example:

    ./launch.sh 2019_08_01

will generate a link table for `/home/j/WORK/11_geospatial/admin_shapefiles/2019_08_01` and place the new link table and id raster there.
