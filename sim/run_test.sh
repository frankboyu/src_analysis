#!/bin/bash

REACTION=phi_d_2H
RUN=90213
EVENTS=10000
TAG=1_initial

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-dev.xml

mkdir -p output/test/${REACTION}_${TAG}
cd output/test/${REACTION}_${TAG}
gluex_MC.py ../../../configs/wrapper_${REACTION}.cfg ${RUN} ${EVENTS} per_file=1000000 batch=0