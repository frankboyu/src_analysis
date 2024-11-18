#!/bin/bash

REACTION=piminus_p_2H
XML_VERSION=ver01_1.1
RUN=90213
EVENTS=10000

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-${XML_VERSION}.xml

cd output

gluex_MC.py ../../configs/wrapper_${REACTION}.cfg ${RUN} ${EVENTS} batch=0 per_file=1000000

rm -r 90*