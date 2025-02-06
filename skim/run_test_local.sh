#!/bin/bash

INPUTFILE=/cache/halld/RunPeriod-2021-11/recon/ver01/REST/090213/dana_rest_090213_000.hddm
REACTION=phi_d_2H

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/analysis-2021_11-ver10.xml

mkdir -p output/${REACTION}_test
cd output/${REACTION}_test
hd_root --config=../../configs/jana_analysis_${REACTION}.cfg ${INPUTFILE}