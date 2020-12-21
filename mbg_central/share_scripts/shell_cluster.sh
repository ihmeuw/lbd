#!/bin/bash

MEMINFO=<<<< FILEPATH REDACTED >>>>
MEMINFO_ROOT=<<<< FILEPATH REDACTED >>>>
if [[ -z ${JOB_ID+x} ]]; then
    # JOB_ID unset. This is not a job, so record by PID as a fallback.
    F="$MEMINFO_ROOT/pids/$$"
else # JOB_ID set
    F="$MEMINFO_ROOT/jobs/$JOB_ID"
    # also handle array jobs
    if [ -z $SGE_TASK_ID ]; then
        F="$f.$SGE_TASK_ID"
    fi
fi

# Determine R path
if [[ $(hostname) == *geos* ]]; then
    RPATH=<<<< FILEPATH REDACTED >>>>
else
    RPATH=<<<< FILEPATH REDACTED >>>>
fi

# Run meminfo and report usage every second
$MEMINFO --loop-interval=1 > $F &
MEMINFO_PID=$!

$RPATH <$1 --no-save $@

RETCODE=$?

# meminfo does not die normally. we must kill it
kill $MEMINFO_PID

exit $RETCODE