#!/bin/bash
#
# -----------------------------------------------------------------------------
# Author      : Ian M. Davis
# Group       : LBD
# Date        : 25 Oct 2018
# Summary     :
# This script will determine how many cores are available for R which has been
# compiled against Intel's Math Kernel Library (MKL) in a qsub job on IHME's
# compute cluster. This script will honor resource requests as either slots
# ("NSLOTS") or threads ("fthread") and attempts to set the appropriate values
# for environmental variables controlling both OpenMP (OMP_NUM_THREADS) and the
# MKL (MKL_NUM_THREADS) so as not to oversubscribe the node or use more cores
# than was granted for the job.
# 
# After the determination of the number of total cores to use, a container is
# launched from the latest Singularity image found in the "image_dir"
# (unless a specific image is given), the version of R in the image is launched,
# and the specified R script is executed.
#
# It is assumed that this script is launched by an R script (such as
# 'parallelize' or 'make_qsub_share' found in
# 'lbd_core/mbg_central/misc_functions.R') as a qsub job.
#
# This following information is used in this script:
#   * host               : lbd (lbd-cluster-p*) node or prod node (cluster2,
#                          intel, or AMD) as determined from $host value.
#   * NSLOTS             : number of slots requested with the job set by UGE.
#   * SGE_HGR_fthread    : The specific number of threads requested for the job
#                          set by UGE.
#   * SGE_HGR_m_mem_free : The max amount of memory requested for the job set by
#                          UGE. We don't do anything with this except echo it
#                          for the user if it is set.
#   * PE                 : Parallel environment environmental set by UGE. If
#                        : is set to "multi_slot" the "NSLOTS" is used to
#                          determine maximum cores available for the job,
#                          otherwise "SGE_HGR_fthread" is used. 
#   * sing_image         : The Singularity image to launch a container from
#                          containing the R MKL version. Default is to located
#                          the latest image found in "image_dir" with the
#                          "get_latest_image()" function
#   * SET_OMP_THREADS    : Number of threads allowed for packages with OpenMP
#                          parallel regions (user defined [default=1])
#   * SET_MKL_THREADS    : Number of threads allowed for packages linked against
#                          Intel's MKL (user defined
#                          [default=${max_threads} calculated by script])
#   * JOB_NAME           : Job name passed to UGE with qsub
#   * JOB_ID             : Job ID passed to UGE with qsub
#   * /proc/cpuinfo      : Used to find the total number of processors on the
#                          host machine in order to calculate cores/slot
#   * qstat              : Used to grep out total number of slots the machine
#                          has been partitioned into in order to calculate
#                          cores/slot
#   * singularity        : This command must be present in order to spin up a
#                          container
#
# Details:
# This script first determines what node it is has landed on (host machine)
# after a `qsub` to execute the R script. It then makes a determination if the
# number of cores allocated to the job is given by "SGE_HGR_fthread" or by 
# "NSLOTS". In the past, all jobs were given a number of slots. As of October
# 29th, 2018, a portion of the cluster began to schedule jobs enforcing that
# both "fthread" (setting the "SGE_HGR_fthread" environmental variable) and
# "m_mem_free" (setting the "SGE_HGR_m_mem_free" environmental variable) be
# specified instead. Since a portion of the cluster still will use "NSLOTS" and
# another users "fthread"/"m_mem_free" over a period of several months, this 
# script supports both types of resource requests in the following way:
#
# fthread/m_mem_free:
# If the script finds that the "SGE_HGR_fthread" and "SGE_HGR_m_mem_free"
# environmental variables are not empty, the "NSLOTS" environmental variable is
# 1, and the "PE" environmental variable does not equal "multi_slot", then the
# value of "SGE_HGR_fthread" is used as the maximum cores (${max_threads})
# allowed for the job. The value of maximum memory from "SGE_HGR_m_mem_free" is
# echoed for the user.
#
# Slots:
# If the script finds that the "PE" environmental variable is "multi_slot" and
# the "NSLOTS" environmental variable is not empty, it will use slots to
# ultimately determine the maximum number of cores allowed for the job. It first
# determines the number of total cores on the machine and total slots allocated
# to the machine so it can calculate the cores/slot ratio. It then determines
# how many slots have been allocated to the job ($NSLOTS) and, using the
# cores/slot ratio, calculates how many cores it actually has to use
# in the calculation ${max_threads}. Most machines have a cores/slot ratio = 1,
# but there are some machines which do not, making the above necessary.
# Further some types of machines (the AMD's) can have a variable number of slots
# they have been partitioned into (51, 64, or 100), giving a variable cores/slot
# ratio. Since ${max_threads} may not be an integer, this script takes some care
# not to aggressively round up the number of cores based on the number slots
# requested and potentially oversubscribe the machine.
#
# After ${max_threads} is determined by either "NSLOTS" or "SGE_HGR_fthread",
# the user settings for SET_OMP_THREADS and SET_MKL_THREADS (if any) are
# considered to determine how to set OMP_NUM_THREADS and MKL_NUM_THREADS
# following suggestions from Intel here:
# https://software.intel.com/en-us/articles/recommended-settings-for-calling-
# intel-mkl-routines-from-multi-threaded-applications
# R with the MKL works best when both OMP_NUM_THREADS as well as MKL_NUM_THREADS
# are set in such as way as to not oversubscribe the machine. In the event that
# both OMP_NUM_THREADS and MKL_NUM_THREADS are both unset, the MKL will grab all
# available cores on the machine. This is probably undesired behavior on a
# shared cluster node. Since most R functions are not SMP parallel the default
# for this script is to set OMP_NUM_THREADS=1 and
# MKL_NUM_THREADS=${max_threads}. There are some packages (like TMB and INLA)
# which itself have OpenMP parallel regions and each parallel region may also
# make use of the MKL. In such a case, it may make sense for the user to define
# how many threads they wish to set for TMB (with OMP_NUM_THREADS, such as the
# number of physical slots on the machine for example), and a separate number of
# MKL threads (with MKL_NUM_THREADS) which each OpenMP parallel region will
# launch.
#
# Below shows the different configurations that may be set by the user:
#
# Scenario  | User Setting           | Result Set by This Script
# ----------|-----------------------------------------------------------
# A         | SET_OMP_THREADS unset  | OMP_NUM_THREADS=1
# [default] | SET_MKL_THREADS unset  | MKL_NUM_THREADS={max_threads}
# ----------------------------------------------------------------------
# B         | SET_OMP_THREADS=N      | OMP_NUM_THREADS=N (N <= {max_threads})
#           | SET_MKL_THREADS unset  | MKL_NUM_THREADS={max_threads}/N
# ----------------------------------------------------------------------
# C         | SET_OMP_THREADS unset  | OMP_NUM_THREADS={max_threads}/M
#           | SET_MKL_THREADS=M      | MKL_NUM_THREADS=M (M <= {max_threads})
# ----------------------------------------------------------------------
# D         | SET_OMP_THREADS=N      | OMP_NUM_THREADS=N
#           | SET_MKL_THREADS=M      | MKL_NUM_THREADS=M
#
# Note that in scenarios A, B, and C, the script goes to the trouble to check
# the machine is not using more cores than have been allocated for the job. In
# all three of these scenarios, the user may not oversubscribe the machine, i.e.
# use more than ${max_threads}. In scenario D, the assumption is that since the
# user is explicitly setting both the SET_OMP_THREADS and SET_MKL_THREADS this
# is for good reason, and the user may oversubscribe the node with a warning.
#
# After a determination is made for SET_OMP_THREADS and SET_MKL_THREADS, these
# values are used to set OMP_NUM_THREADS and MKL_NUM_THREADS.
#
# The script then gives the user some details about the host node, the job, and
# the Singularity image it is making use of.
#
# Finally, the script will launch a Singularity container from a user provided
# image location or, if the 'default' keyword is given, from the latest image
# found by the "get_latest_image()" function in the "image_dir" directory.
#
# This script is dependent on being launched through qsub where UGE sets the
# $JOB_NAME, $JOB_ID, and $NSLOTS variables, and optionally the $PE,
# $SGE_HGR_fthread, $SGE_HGR_m_mem_free variables.
#
# This script is also dependent on the upstream R script that launches it to
# pass the ${sing_image} variable (qsub -v sing_image=...) with either the full
# path to the Singularity image to launch a container from or the 'default'
# keyword indicating to use the most current image. Optionally, the
# SET_OMP_THREADS and SET_MKL_THREADS values may be passed with the '-v' option
# in the qsub as well.
#
# -----------------------------------------------------------------------------
#
# Default directory containing LBD Singularity images
image_dir="/share/singularity-images/lbd"
shell_sing_version="2.0.0"
#
# Set up some empty strings to catch warning messages
slot_warnings=""
cores_warnings=""
#
# <---------------------------------------------------------------------------->
# <---------------------------- Print Script Header --------------------------->
# <---------------------------------------------------------------------------->
echo "--------------------------------------------------------------------------------"
# Let's first print out the head info for the script so if there are errors, the
# user will know what they came from.
echo "R Job Submission Script Info"
echo "----------------------------------------"
echo "Script Name        : 'shell_sing'"
echo -e "Version            : ${shell_sing_version}\n"
#
# <---------------------------------------------------------------------------->
# <----------------------- Define Some Functions We Need ---------------------->
# <---------------------------------------------------------------------------->
#
# Simply function to echo an error message then exit the script
function throw_error() {
 echo -e "ERROR: ${1}"
 echo "Exiting..."
 exit 1
}

# Use 'tr' to make a string all lower case to make our test easier
function make_lower_case() {
  echo $(echo $1 | tr '[:upper:]' '[:lower:]')
}

# Pull the vendor name out of '/proc/cpuinfo'
# Should be either "GenuineIntel" or "AuthenticAMD"
function get_cpu_vendor() {
  local vendor=$(cat /proc/cpuinfo | grep "vendor_id" | uniq | cut -d ':' -f 2)
  echo $(make_lower_case ${vendor})
}

# Unnecessary function to do float math with 'awk' because 'bc' was not
# originally installed on the LBD nodes when they were Ubuntu. It is now under
# CentOS, but no need to change it up now.
#
# Nice solution for rounding with 'awk' used below found here:
# http://bits.mdminhazulhaque.io/linux/round-number-in-bash-script.html
# Returns a rounded integer value of threads is careful not to round up to
# aggressively.  Will only round up if decimal is > 0.8 of a thread. If for some
# reason, the calculated number of threads exceed the number of cores available
# on the machine, we set the number of threads to the total number of cores on
# the machine.
function calc_threads() {
  local cores_per_slot=$1 # cores_per_slot
  local job_slots=$2      # slots obtained for the job
  local max_cores=$3      # total cores on the machine
  # calculate the raw number of cores based on slots and cores
  local raw_threads=$(echo ${cores_per_slot} ${job_slots} | awk '{print $1 * $2}')
  raw_threads=$(echo ${raw_threads} | awk '{print ($0-int($0)<0.799)?int($0):int($0)+1}')
  if [[ ${raw_threads} -le 0 ]]; then
    raw_threads=1
  elif [[ ${raw_threads} -gt ${max_cores} ]]; then
    raw_threads=${max_cores}
  fi
  echo "${raw_threads}"
}

# Some functions that I copied this from "sing_functs.sh" used in launching a
# Singularity container (or R within a Singularity container) for interactive
# use. I wanted to have a single, independent script in the lbd_core repo. These
# functions are simple enough to have repeated here.

# Function to determine if the requested thread value makes sense. Will throw
# errors for non-integer, negative, or zero values.
function check_threads() {
  # use regex to determine if the arg is an integer and exit if it is not
  if [[ ! $1 =~ ^[0-9]+$ ]] || [[ $1 -eq 0 ]]; then
    throw_error "'$1' not an integer number of threads > 0 for $2."
  fi
}

# Get the latest Singularity image file in the default LBD Singularity images
# directory.
function get_latest_image() {
  # Throw an error if we can't access the image directory
  if [[ ! -d "${image_dir}" ]]; then
    throw_error "Could access image directory: '${image_dir}'"
  fi

  # Use `ls` to find the latest Singularity image (with *.simg extension)
  # Will give the complete path to the file
  local latest_image=$(ls -t "${image_dir}"/*.simg 2> /dev/null | head -1)
  if [[ -z "${latest_image}" ]]; then
    throw_error "No *.simg images found in '${image_dir}'"
  fi

  # Return the image name and complete path to file
  echo "${latest_image}"
}

# Echo's out some information on the Singularity image file that is nice for the
# user to know
function echo_sing_info() {
  # Expecting on a single argument
  if [[ $# -ne 1 ]]; then
    throw_error "function 'echo_sing_info' expects the Singularity image name as the only argument"
  fi
  if [[ ! -f "${1}" ]]; then
    throw_error "Can't access Singularity image: '${1}'"
  fi

  # Get the name of the Singularity image file only as well as the last
  # modification date.
  local latest_image_name=$(basename ${1})
  local latest_image_dir=$(dirname ${1})
  local image_creation_date=$(stat -c %y ${1} | cut -f1 -d ".")

  # Echo some details about the container we are attempting to launch and then
  # launch the container with Singularity run
  echo "Image name         : '${latest_image_name}'"
  echo "Image directory    : ${latest_image_dir}"
  echo "Creation date      : ${image_creation_date}"
  echo "Singularity version: $(singularity --version)"
}

# <---------------------------------------------------------------------------->
# <------------------- Determine threads for MKL and OpenMP ------------------->
# <---------------------------------------------------------------------------->
# Grab the hostname as a global since we use it all over, and error if it comes
# up empty since we need it to determine the number of slots
host=$(hostname -s)
if [[ -z "${host}" ]]; then
  throw_error "Could not obtain host name"
fi

# Throw an error if we can't find the "singularity" executable
if [[ -z $(command -v singularity) ]]; then
  throw_error "Could not find 'singularity' executable."
fi

# We always want to know the number of cores on the machine so we can inform the
# user, but we need it in the case where NSLOTS are specified for the machine
# resource.
#
# First Make sure that /proc/cpuinfo exists, because we use it to determine how
# many cores are on the machine (faster than grepping this out of the `lscpu`
# output) and we also use it to determine the CPU vendor. If it exists, use it
# to find the number of cores on the host machine.
if [[ ! -f "/proc/cpuinfo" ]]; then
  throw_error "no /proc/cpuinfo exists to find number of cores"
else
  machine_cores=$(cat /proc/cpuinfo | grep "processor" | wc -l)
fi

# Make sure we got something for machine_cores and make sure that the value
# actually makes sense
if [[ -z $machine_cores ]]; then
  throw_error "Was not able to determine number of cores on machine"
elif [[ ! $machine_cores =~ ^[0-9]+$ ]] || [[ $machine_cores -eq 0 ]]; then
  throw_error "Invalid value found for machine cores: '${machine_cores}'"
fi

# Now let's determine if this is a "slots" based request or an "fthread" based
# request. If this is an "fthread" based request, we don't need to go to all
# the trouble to determine how many cores are available to us on the host
# machine based on slots.
#
# Jobs on the current scheduler will be specified with both "fthread" and
# "m_mem_free", and NSLOTS will be specifically set to 1. PE should not be 
# set at all, but definitely not to "multi_slot" if "SGE_HGR_fthread" is set.
# In this case, our "max_threads" is just equal the "SGE_HGR_fthread".
if [[ -n "${SGE_HGR_fthread}" ]] && [[ -n "${SGE_HGR_m_mem_free}" ]] &&
   [[ "${PE}" != "multi_slot" ]] && [[ "${NSLOTS}" == 1 ]]; then
  request="fthread"
  max_threads="${SGE_HGR_fthread}"
  check_threads "${max_threads}" "SGE_HGR_fthread"
# Jobs on the old scheduler will have the environmental variable "PE" set to
# "multi_slot", if not it will be empty and this test will fail. Also "NSLOTS"
# should be set to something.
elif [[ "${PE}" == "multi_slot" ]] && [[ -n "${NSLOTS}" ]]; then
  request="slots"
else
  throw_error "Could not determine if resources requested with 'NSLOTS' or 'fthread'"
fi

# If the resources are requested with slots, then we need to do some more work
# to calculate max_threads from the number slots allocated to the machine
# ("machine_slots") and total cores on the machine "machine_cores". 
if [[ "${request}" == "slots" ]]; then
  # Determine the total slots the machine has been partitioned to in order to
  # calculate cores/slot
  machine_slots=$(qstat -f -l hostname=${host} | sed -n 3p | awk -F '( *|/)' '{print $5}')
  # Use a default value if `qstat` gave us trouble
  if [[ -z $machine_slots ]] || [[ ! $machine_slots =~ ^[0-9]+$ ]] || [[ $machine_slots -eq 0 ]]; then
    slot_warnings+="----> WARNING: Could not determine machine slots from 'qstat'.\n"
    # AMD nodes are either 100 or 51 slots, so being conservative if we can't find
    # the true total slot value.
    if [[ "${host}" == *"cn"* ]] && [[ "$(get_cpu_vendor)" == *"amd"* ]]; then
      machine_slots=100
      slot_warnings+="---->          Using default of: ${machine_slots}"
    else
      machine_slots=${machine_cores}
      slot_warnings+="---->          Using default of machine cores: ${machine_slots}"
    fi
  fi

  # Now that we know the total slots and the number of cores on the machine, we
  # can calculate cores/slot:
  cores_per_slot="$(echo ${machine_cores} ${machine_slots} | awk '{print $1 / $2}')"

  # Get the max threads based on $NSLOTS requested and machine:
  max_threads=$(calc_threads ${cores_per_slot} ${NSLOTS} ${machine_cores})

  # Make sure max_threads is an integer > 0:
  check_threads "${max_threads}" "maximum threads calculated from 'NSLOTS'"
fi

# Now let's give a warning if the maximum threads is greater than the total
# machine cores
if [[ "${max_threads}" -gt "${machine_cores}" ]]; then
  cores_warnings+="\n----> WARNING: Maximum threads > Total machine cores\n"
  cores_warnings+="---->          Max threads   = ${max_threads}\n"
  cores_warnings+="---->          Machine cores = ${machine_cores}\n"
fi

# Set OMP_NUM_THREADS and MKL_NUM_THREADS based on user settings and number of
# ${max_threads} determined above
#
# ------------------------------------------------------------------------------
# | Scenario A |
# --------------
# Both SET_OMP_THREADS and SET_MKL_THREADS are unset (default)
# ------------------------------------------------------------------------------
if [[ -z "${SET_OMP_THREADS}" ]] && [[ -z "${SET_MKL_THREADS}" ]]; then
  SET_OMP_THREADS=1
  SET_MKL_THREADS=${max_threads}

# ------------------------------------------------------------------------------
# | Scenario B |
# --------------
# SET_OMP_THREADS set but SET_MKL_THREADS is unset
# ------------------------------------------------------------------------------
elif [[ -n "${SET_OMP_THREADS}" ]] && [[ -z "${SET_MKL_THREADS}" ]]; then
  # Check to make sure that the requested thread value make sense.
  check_threads "${SET_OMP_THREADS}" "SET_OMP_THREADS"
  # Check for oversubscription
  if [[ "${SET_OMP_THREADS}" -gt "${max_threads}" ]]; then
    cores_warnings+="\n----> WARNING: SET_OMP_THREADS requested > available cores. Setting as follows:\n"
    cores_warnings+="---->          OMP_NUM_THREADS=${max_threads}\n"
    cores_warnings+="---->          MKL_NUM_THREADS=1\n"
    SET_OMP_THREADS=${max_threads}
    SET_MKL_THREADS=1
  else
    SET_MKL_THREADS=$((${max_threads} / ${SET_OMP_THREADS}))
  fi

# ------------------------------------------------------------------------------
# | Scenario C |
# --------------
# SET_OMP_THREADS unset but SET_MKL_THREADS is set
# ------------------------------------------------------------------------------
elif [[ -z "${SET_OMP_THREADS}" ]] && [[ -n "${SET_MKL_THREADS}" ]]; then
  # Check to make sure that the requested thread value make sense.
  check_threads "${SET_MKL_THREADS}" "SET_MKL_THREADS"
  # Check for oversubscription
  if [[ "${SET_MKL_THREADS}" -gt "${max_threads}" ]]; then
    cores_warnings+="\n----> WARNING: SET_MKL_THREADS requested > available cores. Setting as follows:\n"
    cores_warnings+="---->          OMP_NUM_THREADS=1\n"
    cores_warnings+="---->          MKL_NUM_THREADS=${max_threads}\n"
    SET_OMP_THREADS=1
    SET_MKL_THREADS=${max_threads}
  else
    SET_OMP_THREADS=$((${max_threads} / ${SET_MKL_THREADS}))
  fi

# ------------------------------------------------------------------------------
# | Scenario D |
# --------------
# BOTH SET_OMP_THREADS and SET_MKL_THREADS set
# ------------------------------------------------------------------------------
else
  # Check to make sure that the requested thread values make sense.
  check_threads "${SET_MKL_THREADS}" "SET_MKL_THREADS"
  check_threads "${SET_OMP_THREADS}" "SET_OMP_THREADS"
  # Check for oversubscription
  if [[ "$((${SET_OMP_THREADS} * ${SET_MKL_THREADS}))" -gt "${max_threads}" ]]; then
    cores_warnings+="\n----> WARNING: SET_OMP_THREADS * SET_MKL_THREADS > available cores.\n"
    cores_warnings+="---->          [$((${SET_OMP_THREADS} * ${SET_MKL_THREADS})) > ${max_threads}]\n"
    cores_warnings+="---->          You may be oversubscribing the node!\n"
  fi
fi

# At this point, if either SET_MKL_THREADS or SET_OMP_THREADS is empty, there's
# a problem
if [[ -z "${SET_OMP_THREADS}" ]] || [[ -z "${SET_MKL_THREADS}" ]]; then
  error_message="Either 'SET_OMP_THREADS' or 'SET_MKL_THREADS' empty:\n"
  error_message+="       SET_OMP_THREADS=${SET_OMP_THREADS}\n"
  error_message+="       SET_MKL_THREADS=${SET_MKL_THREADS}"
  throw_error "${error_message}"
else
  # If both are assigned, let's double check that the values make sense, since
  # calcs may have been done to determine them
  check_threads "${SET_MKL_THREADS}" "SET_MKL_THREADS"
  check_threads "${SET_OMP_THREADS}" "SET_OMP_THREADS"
fi

# <---------------------------------------------------------------------------->
# <------------------------ Print Remaining User Info ------------------------->
# <---------------------------------------------------------------------------->
echo "Machine Info"
echo "----------------------------------------"
echo "Node               : ${host}"
echo "Total Cores        : ${machine_cores}"
# Throw any warnings related to slots that have accumulated
if [[ -n "${slot_warnings}" ]]; then
  echo -e "${slot_warnings}"
fi
# Echo slot related info if resources request was based on "NSLOTS"
if [[ "${request}" == "slots" ]]; then
  echo "Total Slots        : ${machine_slots}"
  echo -e "Cores/Slot         : ${cores_per_slot}\n"
else
  echo ""
fi
echo "Job Info"
echo "----------------------------------------"
echo "Job Name           : ${JOB_NAME}"
echo "Job Id             : ${JOB_ID}"
echo "Execution Date     : $(date '+%c')"
if [[ "${request}" == "slots" ]]; then
  echo "Slots              : ${NSLOTS}"
fi
# Echo cores/memory related info if resources request was based on "fthread"
if [[ "${request}" == "fthread" ]]; then
  echo "Max Threads        : ${SGE_HGR_fthread}"
  echo "Max Memory         : ${SGE_HGR_m_mem_free}"
fi
# Throw any warnings related to cores that have accumulated
if [[ -n "${cores_warnings}" ]]; then
  echo -e "${cores_warnings}"
fi
echo "OMP_NUM_THREADS    : ${SET_OMP_THREADS}"
echo -e "MKL_NUM_THREADS    : ${SET_MKL_THREADS}\n"
echo "Singularity Image Info"
echo "----------------------------------------"

# Make sure image is defined or the default is specified
if [[ -z "${sing_image}" ]]; then
  throw_error "No Singularity image supplied"
elif [[ "${sing_image}" == "default" ]]; then
  sing_image=$(get_latest_image)
  echo "Launching Singularity container from default image:"
else
  echo "Launching Singularity container from user specified image:"
fi
echo_sing_info "${sing_image}"
echo "--------------------------------------------------------------------------------"

# <---------------------------------------------------------------------------->
# <----------------------------- Launch Container ----------------------------->
# <---------------------------------------------------------------------------->


SINGULARITYENV_R_EXEC="/usr/local/bin/R < $1 --no-save $@"        \
SINGULARITYENV_UMASK_ORIG=$(umask)                                \
SINGULARITYENV_OMP_NUM_THREADS=${SET_OMP_THREADS}                 \
SINGULARITYENV_MKL_NUM_THREADS=${SET_MKL_THREADS}                 \
SINGULARITYENV_OMP_NESTED=true                                    \
SINGULARITYENV_MKL_DYNAMIC=false                                  \
singularity exec --bind /tmp:/tmp "${sing_image}" /bin/sh -c      \
    'umask "${UMASK_ORIG}" && eval ${R_EXEC}'
