#!/bin/sh
SINGULARITYENV_OMP_NUM_THREADS=10 <<<< FILEPATH REDACTED >>>>/singularity exec --cleanenv <<<< FILEPATH REDACTED >>>>/r_pkgs3.4.3gcc7mkl.simg <<<< FILEPATH REDACTED >>>>/bin/R <$1 --no-save $@
