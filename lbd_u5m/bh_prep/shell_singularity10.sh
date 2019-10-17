#!/bin/sh
SINGULARITYENV_OMP_NUM_THREADS=10 <<<< FILEPATH REDACTED >>>> exec --cleanenv <<<< FILEPATH REDACTED >>>> <<<< FILEPATH REDACTED >>>> <$1 --no-save $@
