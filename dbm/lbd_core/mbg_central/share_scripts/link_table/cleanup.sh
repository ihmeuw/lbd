#!/bin/bash

# This script runs a final cleanup (R) script and notifies the user via email
#
# Arguments:
#
# $1 run_dir - directory the parent script is running in (for referencing)
# $2 n - number of jobs run
# $3 tmp_dir - temporary directory to check for link table files
# $4 shape_dir - shapefile directory

script_dir=$1
n=$2
tmp_dir=$3
shape_dir=$4

# image for sing_shell.sh
export sing_image="default"

# run cleanup R script
msg=$("$script_dir/../shell_sing.sh" "$script_dir/run" --args -c cleanup -n "$n" -d "$tmp_dir" -s "$shape_dir")

# email user
echo "$msg" | grep -o '\[1].*' | cut -d '"' -f2 | xargs -r0 printf | mail -s "Link Table Status" "$USER"@uw.edu

# remove temp directory
rm -r "$tmp_dir"
