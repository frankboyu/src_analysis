#!/bin/bash

INPUTFILE=/cache/halld/RunPeriod-2021-11/recon/ver04/REST/090213/dana_rest_090213_000.hddm
REACTION=phi_d_2H

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /group/halld/www/halldweb/html/halld_versions/version_6.3.0.xml

mkdir -p output/${REACTION}_test
cd output/${REACTION}_test
hd_root --loadconfigs /work/halld2/home/boyu/src_analysis/data/configs/jana_analysis_${REACTION}.cfg ${INPUTFILE}