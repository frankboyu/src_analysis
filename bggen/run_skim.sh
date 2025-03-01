#!/bin/bash

REACTION=$1

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml

python launch.py configs/jobs_analysis_bggen_${REACTION}.cfg 90034 90034
swif2 run -workflow src_analysis_bggen -maxconcurrent 1000