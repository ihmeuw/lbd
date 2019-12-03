#!/bin/bash

MEMINFO="<<< FILEPATH REDACTED >>>"
MEMINFO_ROOT="<<< FILEPATH REDACTED >>>"
if [[ -z ${JOB_ID+x} ]]; then
    # JOB_ID unset. This is not a job, so record by PID as a fallback.
    F="<<< FILEPATH REDACTED >>>"
else # JOB_ID set
    F="<<< FILEPATH REDACTED >>>"
    # also handle array jobs
    if [ -z $SGE_TASK_ID ]; then
        F="$f.$SGE_TASK_ID"
    fi
fi

# Run meminfo and report usage every second
$MEMINFO --loop-interval=1 > $F &
MEMINFO_PID=$!

/usr/bin/R <$1 --no-save $@

RETCODE=$?

# meminfo does not die normally. we must kill it
kill $MEMINFO_PID

exit $RETCODE
