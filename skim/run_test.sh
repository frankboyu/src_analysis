#!/bin/bash

INPUTFILE="/cache/halld/RunPeriod-2021-11/recon/ver01/REST/090213/dana_rest_090213_000.hddm"
JANA_CONFIG_FILE="phi_c_2H"
TAG="1_initial"

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.19.0.xml

mkdir output/${JANA_CONFIG_FILE}_${TAG}
cd output/${JANA_CONFIG_FILE}_${TAG}
hd_root --config=../../configs/jana_analysis_${JANA_CONFIG_FILE}.cfg ${INPUTFILE}