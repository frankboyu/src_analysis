#!/bin/bash

REACTION=phi_d_2H
WRAPPER_VERSION=gluex_MCwrapper_srcct-v2.11.0
RUN=90213
EVENTS=10000

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml
export MCWRAPPER_CENTRAL=/work/halld2/home/boyu/src_software_builds/${WRAPPER_VERSION}
export PATH=/work/halld2/home/boyu/src_software_builds/${WRAPPER_VERSION}:${PATH}

gluex_MC.py configs/wrapper_${REACTION}.cfg ${RUN} ${EVENTS} batch=2 per_file=1000000