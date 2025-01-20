#!/bin/bash

REACTION=phi_d_2H
XML_VERSION=dev
RUN=90213
EVENTS=1000
OPTIONS='cleangenerate=0 cleangeant=0 cleanmcsmear=0'

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-${XML_VERSION}.xml

cd output

gluex_MC.py ../configs/wrapper_${REACTION}.cfg ${RUN} ${EVENTS} ${OPTIONS} batch=0 per_file=1000000

rm -r 90*