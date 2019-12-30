#!/bin/bash

# Author: Tim Whitson (twhitson@uw.edu)
# Date: July 2019
# Details:
#   This script utilizes the following model for building the link table in parallel:
#       1. Partition the shapefile
#       2. Build the link table for each partition
#       3. Combine the link tables
#       4. Cleanup process (check if all jobs completed successfully, remove temp files)
#
# Arguments:
#   $1 - shapefile directory (directory to pull shapefile from/save link table to)
#
# Notes:
#   You must QLOGIN before running this script.

shape_dir=$1
script_dir="$(cd "$( dirname "${BASH_SOURCE[0]}" )" || exit ; pwd -P)"
core_dir="$(cd "$script_dir/../../.." || exit ; pwd -P)"

export core_dir # export for use in R scripts

# r shell script
r_shell="$script_dir"/../shell_sing.sh

# get base shapefile directory from mbg scripts
# then parse the output from singR, which will print numerous lines of
# code followed by, for example:
# [1] "<<< FILEPATH REDACTED >>>
# so we just need to find [1] and trim the quotes
echo "Getting base shapefile directory ..."
base_shapefile_dir=$("<<< FILEPATH REDACTED >>>")

# ensure shapefile exists in given directory
if [[ ! -f "$base_shapefile_dir/$shape_dir/lbd_standard_admin_2.shp" ]]; then
    echo -e "$base_shapefile_dir/$shape_dir is not a valid shapefile directory.\\nExiting ..."
    exit 1
fi

# ensure user is not providing "current" directory
if [ "$shape_dir" == "current" ]; then
    echo -e "'current' should not be modified. Please provide a different directory.\\nExiting ..."
    exit 1
fi

# Partition the shapefile
n=50
temp_dir=$('<<< FILEPATH REDACTED >>>>')
job_name="partition_lt"
q="geospatial.q"
p="proj_geo_nodes"
mem="10G"
threads=4
runtime="4:00:00"
sing_image="default"
r_script="$script_dir/run"
output="<<< FILEPATH REDACTED >>>"
errors="<<< FILEPATH REDACTED >>>"

echo "Logging output to $output."
echo "Logging errors to $errors."

pjid=$(qsub -terse \
		-N $job_name \
		-q $q \
		-P $p \
		-l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
		-v sing_image=$sing_image,core_dir="$core_dir"\
		-o "$output" \
        -e "$errors" \
		"$r_shell" \
		"$r_script" --args -c partition -n $n -d "$temp_dir" -s "$shape_dir")

echo "Partition job $pjid queued."

# Build link tables
job_name="build_lt"
mem="30G"

bjid=$(qsub -terse \
	-N $job_name \
	-t 1:$n \
	-q $q \
	-P $p \
	-l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
	-hold_jid "$pjid" \
	-v sing_image=$sing_image,core_dir="$core_dir" \
	-o "$output" \
    -e "$errors" \
	"$r_shell" \
	"$r_script" --args -c build -d "$temp_dir" -s "$shape_dir" \
	| cut -d "." -f 1) # cut to remove decimal from job id

echo "Build job $bjid queued."

# Combine link tables
job_name="combine_lt"
mem="20G"

cjid=$(qsub -terse \
    -N $job_name \
    -q $q \
    -P $p \
    -l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
    -hold_jid "$bjid" \
    -v sing_image=$sing_image,core_dir="$core_dir" \
    -o "$output" \
    -e "$errors" \
    "$r_shell" \
    "$r_script" --args -c combine -d "$temp_dir" -s "$shape_dir" -n $n)

echo "Combine job $cjid queued."

# Check if job was successful and cleanup
job_name="cleanup_lt"
mem="1G"
threads=2
runtime="1:00:00"
shell_script="$script_dir/cleanup.sh"

cljid=$(qsub -terse \
    -N $job_name \
    -q $q \
    -P $p \
    -l m_mem_free=$mem,fthread=$threads,h_rt=$runtime \
    -hold_jid "$cjid" \
    -v sing_image=$sing_image,core_dir="$core_dir" \
    -o "$output" \
    -e "$errors" \
    "$shell_script" "$script_dir" $n "$temp_dir" "$shape_dir")

echo "Cleanup job $cljid queued."
