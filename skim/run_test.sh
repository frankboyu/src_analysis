#!/bin/bash

INPUTFILE="/cache/halld/RunPeriod-2021-11/recon/ver01/REST/090213/dana_rest_090213_000.hddm"
REACTION="phi_c_2H"
TAG="2_version511"

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.11.0.xml

mkdir output/${REACTION}_${TAG}
cd output/${REACTION}_${TAG}
hd_root --config=../../configs/jana_analysis_${REACTION}.cfg ${INPUTFILE}