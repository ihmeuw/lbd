#!/bin/bash
#
# Summary     :
# This script will determine the appropriate number of threads to use for R
# which has been compiled against Intel's Math Kernel Library (MKL) based on
# the number of slots requested for the job and how many total cores are
# available on the host machine the job lands on through UGE via qsub. After
# this determination, a container is launched from a specified Singularity
# image containing this version of R and the specified R script is executed.
#
# It is assumed that this script is launched by an R script (such as
# 'parallelize' or 'make_qsub_share' found in
# 'lbd_core/mbg_central/misc_functions.R') as a qsub job.
#
# This following information is used in this script:
#   * host            : Geos node or prod node (cluster2, intel, or AMD) as
#                       determined from $host value.
#   * NSLOTS          : number of slots requested with the job. UGE sets this
#                       environmental variable with slots requested at the time
#                       of the qsub.
#   * sing_image      : The Singularity image to launch a container from
#                       containing the R MKL version. Default is hardcode in
#                       this script to:
#                       <<<< FILEPATH REDACTED >>>>/r_pkgs3.4.3gcc7mkl.simg
#   * SET_OMP_THREADS : Number of threads allowed for packages with OpenMP
#                       parallel regions (user defined [default=1])
#   * SET_MKL_THREADS : Number of threads allowed for packages linked against
#                       Intel's MKL (user defined
#                       [default=${max_threads} calculated by script])
#   * JOB_NAME        : Job name passed to UGE with qsub
#   * JOB_ID          : Job ID passed to UGE with qsub
#
# Details     :
# This script first determines what node it is has landed on. It has the number
# of cores/slot hardcoded, so it can calculate how many threads are available
# ${max_threads} based on cores/slot and $NSLOTS passed in by UGE without
# oversubscribing the machine. If $NSLOTS is unset, it will exit with an error.
# Since ${max_threads} may not be an integer, this script takes some care not
# to aggressively round up the number of cores based on the number slots
# requested and potentially oversubscribe the machine. If this script does not
# recognize the machine, and therefore the number of cores/slot, it will print
# a warning and set ${max_threads} to 1/2 $NSLOTS.
#
# After ${max_threads} is determined, the user settings for SET_OMP_THREADS and
# SET_MKL_THREADS (if any) are considered to deterimine how to set
# OMP_NUM_THREADS and MKL_NUM_THREADS following suggestions from Intel here:
# https://software.intel.com/en-us/articles/recommended-settings-for-calling-intel-mkl-routines-from-multi-threaded-applications
# R with the MKL works best when both OMP_NUM_THREADS as well as
# MKL_NUM_THREADS are set in such as way as to not oversubscribe the machine.
# In the event that both OMP_NUM_THREADS and MKL_NUM_THREADS are both unset,
# the MKL will grab all available cores on the machine. This is probably not
# desired on a shared cluster node. Since most R functions are not SMP parallel
# the default set up below is to set OMP_NUM_THREADS=1 and
# MKL_NUM_THREADS=${max_threads}. There are some packages (like TMB) which
# itself has OpenMP parallel regions and each parallel region may also make use
# of the MKL. In such a case, it may make sense for the user to define how many
# threads they wish to set for TMB (with OMP_NUM_THREADS, such as the number of
# physical slots on the machine for example), and a separate number of MKL
# threads (with MKL_NUM_THREADS) which each OpenMP parallel region will launch.
# Below shows the different configurations that may be set by the user:
#
# Scenario  | User Setting           | Result Set by This Script
# ----------|-----------------------------------------------------------
# A         | SET_OMP_THREADS unset  | OMP_NUM_THREADS=1
# [default] | SET_MKL_THREADS unset  | MKL_NUM_THREADS={max_threads}
# ----------------------------------------------------------------------
# B         | SET_OMP_THREADS=N      | OMP_NUM_THREADS=N
#           | SET_MKL_THREADS unset  | MKL_NUM_THREADS={max_threads}/N
# ----------------------------------------------------------------------
# C         | SET_OMP_THREADS unset  | OMP_NUM_THREADS=1
#           | SET_MKL_THREADS=M      | MKL_NUM_THREADS=M
# ----------------------------------------------------------------------
# D         | SET_OMP_THREADS=N      | OMP_NUM_THREADS=N
#           | SET_MKL_THREADS=M      | MKL_NUM_THREADS=M
#
# Note that in scenarios A and B, the script goes to some care to check the
# machine is not oversubscribed, but in scenario's C & D the assumption is that
# your settings are for good reason, and the user still may oversubcribe the
# node.
#
# After a determination is made for SET_OMP_THREADS and SET_MKL_THREADS, these
# values are used to set OMP_NUM_THREADS and MKL_NUM_THREADS.
#
# Finally, the script will launch a Singularity container from a user provided
# image location or, if the 'default' keyword is given, from a default image.
# The default image location is maintained in this script so that the default
# can be updated in this script without changing any upstream R code.
#
# The current default image is:
# <<<< FILEPATH REDACTED >>>>/r_pkgs3.4.3gcc7mkl.simg
#
# This script is dependent on being launched through qsub where UGE sets the
# $JOB_NAME, $JOB_ID, and $NSLOTS variables. These are used here, as well as
# specific machine name prefixes on the prod and geos clusters which correspond
# to a certain hardcoded cores/slot value.
#
# This script is also dependent on the upstream R script that launches it to
# pass $sing_image variable (qsub -v sing_image=...) with either the full path
# to the Singularity image to launch a container from or the 'default' keyword
# to use the current default image. Optionally, the SET_OMP_THREADS and
# SET_MKL_THREADS values can be passed with the '-v' option in the qsub as
# well.
#
# Updated after "geos" nodes were renamed to "lbd-cluster-p*"
#
# -----------------------------------------------------------------------------
#
# <-------------------------------------------------------------------->
# <------------------- Define Some Functions We Need ------------------>
# <-------------------------------------------------------------------->
#
# Use 'tr' to make a string all lower case to make our test easier
function make_lower_case() {
  echo $(echo $1 | tr '[:upper:]' '[:lower:]')
}

# Pull the vendor name out of '<<<< FILEPATH REDACTED >>>>'
# Should be either "GenuineIntel" or "AuthenticAMD"
function get_cpu_vendor() {
  local vendor=$(cat <<<< FILEPATH REDACTED >>>> | grep "vendor_id" | uniq | cut -d ':' -f 2)
  echo $(make_lower_case ${vendor})
}

# Unnecessary function to do float math with 'awk' because 'bc' is not
# installed on the geos nodes.
# Nice solution for rounding with 'awk' used below found here:
# http://bits.mdminhazulhaque.io/linux/round-number-in-bash-script.html
# Returns a rounded integer value of threads is careful not to round up
# to aggressively.  Will only round up if decimal is > 0.8 of a thread
function calc_threads() {
  # calculate the raw number of cores based on slots and cores
  local raw_threads=$(echo $1 $2 | awk '{print $1 * $2}')
  raw_threads=$(echo ${raw_threads} | awk '{print ($0-int($0)<0.799)?int($0):int($0)+1}')
  if [[ "${raw_threads}" -le 0 ]]; then
    raw_threads=1
  fi
  echo "${raw_threads}"
}

# Based on the machine name, we know the number of cores/slot and UGE
# has set the 'NSLOTS' environmental variable when this script was
# qsub'ed. This function calculates the maximum number of threads for
# each machine, limiting it to the number of hyperthreaded cores
# (physical cores x 2)

function get_threads() {
  # set some variables
  local host=$(hostname)
  local threads=0
  # geos nodes
  if [[ $host == *"lbd-cluster"* ]]; then
    threads=$(calc_threads 1.04 ${NSLOTS})
    if [[ "${threads}" -gt 56 ]]; then
      threads=56
    fi
  # c2-6f ("cluster 2") prod nodes
  elif [[ $host == *"c2-6f"* ]]; then
    threads=$(calc_threads 1.02 ${NSLOTS})
    if [[ "${threads}" -gt 56 ]]; then
      threads=56
    fi
  # Older prod nodes have different cores/slot if they are Intel or AMD
  elif [[ $host == *"cn"* ]] && [[ "$(get_cpu_vendor)" == *"intel"* ]]; then
    threads=$(calc_threads 0.83 ${NSLOTS})
    if [[ "${threads}" -gt 40 ]]; then
      threads=40
    fi
  elif [[ $host == *"cn"* ]] && [[ "$(get_cpu_vendor)" == *"amd"* ]]; then
    threads=$(calc_threads 0.64 ${NSLOTS})
    if [[ "${threads}" -gt 64 ]]; then
      threads=64
    fi
  else
    threads="unknown"
  fi
  echo "${threads}"
}

# <-------------------------------------------------------------------->
# <---------------------- Calculate Threads to Use -------------------->
# <-------------------------------------------------------------------->

# If this script was launched through `qsub` (as it should have been),
# $NSLOTS should have been set by UGE. If this is not set, exit, since
# it is sort of the whole point
if [[ -z "${NSLOTS}" ]]; then
  echo "'NSLOTS' environmental not set.  Did you qsub this script?"
  echo "Exiting..."
  exit 1
fi

# First, get the max threads based on $NSLOTS requested and machine:
max_threads=$(get_threads)
# Warn the user if the node is unknown and set the threads to 1/2 slots
if [[ "${max_threads}" == "unknown" ]]; then
  echo "----> WARNING: Machine $(hostname) Unknown"
  echo "----> max_threads set to 1/2 requested slots"
  max_threads=$(calc_threads 0.5 ${NSLOTS})
fi

# Set OMP_NUM_THREADS and MKL_NUM_THREADS based on user settings and
# ${max_threads} calculated from slots and type of machine
#
# If both SET_OMP_THREADS and SET_MKL_THREADS are unset (default)
# This is scenario A as described above
if [[ -z "${SET_OMP_THREADS}" ]] && [[ -z "${SET_MKL_THREADS}" ]]; then
  SET_OMP_THREADS=1
  SET_MKL_THREADS=${max_threads}
# If the user has supplied SET_OMP_THREADS but SET_MKL_THREADS is unset
# This is scenario B as described above
elif [[ -n "${SET_OMP_THREADS}" ]] && [[ -z "${SET_MKL_THREADS}" ]]; then
  if [[ "${SET_OMP_THREADS}" -ge "${max_threads}" ]]; then
    echo "----> WARNING: SET_OMP_THREADS requested > available cores. Setting as follows:"
    echo "----> OMP_NUM_THREADS=${max_threads}"
    echo "----> MKL_NUM_THREADS=1"
    SET_OMP_THREADS=${max_threads}
    SET_MKL_THREADS=1
  elif [[ "${SET_OMP_THREADS}" -le 0 ]]; then
    echo "----> WARNING: SET_OMP_THREADS requested <= 0. Setting as follows:"
    echo "----> OMP_NUM_THREADS=1"
    echo "----> MKL_NUM_THREADS=${max_threads}"
    SET_OMP_THREADS=1
    SET_MKL_THREADS=${max_threads}
  else
    SET_MKL_THREADS=$((${max_threads} / ${SET_OMP_THREADS}))
  fi
# If the user has supplied SET_MKL_THREADS but SET_OMP_THREADS is unset
# This is scenario C as described above
elif [[ -z "${SET_OMP_THREADS}" ]] && [[ -n "${SET_MKL_THREADS}" ]]; then
  if [[ "${SET_MKL_THREADS}" -ge "${max_threads}" ]]; then
    echo "----> WARNING: SET_MKL_THREADS requested > available cores."
    echo "----> You may be oversubscribing the node!"
    SET_OMP_THREADS=1
  elif [[ "${SET_MKL_THREADS}" -le 0 ]]; then
    echo "----> WARNING: SET_MKL_THREADS requested <= 0. Setting as follows:"
    echo "----> OMP_NUM_THREADS=1"
    echo "----> MKL_NUM_THREADS=${max_threads}"
    SET_OMP_THREADS=1
    SET_MKL_THREADS=${max_threads}
  else
    SET_OMP_THREADS=1
  fi
# If the user has supplied both SET_MKL_THREADS and SET_OMP_THREADS
# This is scenario D as described above
else
  if [[ "$((${SET_OMP_THREADS} * ${SET_MKL_THREADS}))" -ge "${max_threads}" ]]; then
    echo "----> WARNING: SET_OMP_THREADS * SET_MKL_THREADS > available cores."
    echo "----> [$((${SET_OMP_THREADS} * ${SET_MKL_THREADS})) > ${max_threads}]"
    echo "----> You may be oversubscribing the node!"
  fi
fi

# Let's tell the user what is going on:
echo "--------------------------------------------------------------------------------"
echo "Job Name           : ${JOB_NAME}"
echo "Job Id             : ${JOB_ID}"
echo "Execution Date     : $(date '+%c')"
echo "Node               : $(hostname)"
echo "Slots              : ${NSLOTS}"
echo "OMP_NUM_THREADS    : ${SET_OMP_THREADS}"
echo "MKL_NUM_THREADS    : ${SET_MKL_THREADS}"
# Make sure image is defined or the default is specified
if [[ -z "${sing_image}" ]]; then
  echo "No Singularity image supplied"
  echo "Exiting..."
  exit 1
elif [[ "${sing_image}" == "default" ]]; then
  echo "Launching Singularity container from default image:"
  echo "<<<< FILEPATH REDACTED >>>>/r_pkgs3.4.3gcc7mkl.simg"
  sing_image="<<<< FILEPATH REDACTED >>>>/r_pkgs3.4.3gcc7mkl.simg"
else
  echo "Launching Singularity container from specified image:"
  echo "${sing_image}"
fi
echo "--------------------------------------------------------------------------------"

# Launch the container
SINGULARITYENV_OMP_NUM_THREADS=${SET_OMP_THREADS}                 \
SINGULARITYENV_MKL_NUM_THREADS=${SET_MKL_THREADS}                 \
SINGULARITYENV_OMP_NESTED=true SINGULARITYENV_MKL_DYNAMIC=false   \
<<<< FILEPATH REDACTED >>>>/singularity exec --cleanenv    \
${sing_image} <<<< FILEPATH REDACTED >>>>/bin/R <$1 --no-save $@
