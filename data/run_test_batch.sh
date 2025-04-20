#!/bin/bash

REACTION=phi_d_2H
RUNMIN=90212
RUNMAX=90215

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

python launch.py configs/jobs_analysis_${REACTION}.cfg ${RUNMIN} ${RUNMAX}
swif2 run -workflow src_analysis_data