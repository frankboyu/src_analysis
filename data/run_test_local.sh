#!/bin/bash

INPUTFILE=/cache/halld/RunPeriod-2021-11/rawdata/Run090213/hd_rawdata_090213_000.evio
REACTION=recon_2021-11

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /group/halld/www/halldweb/html/halld_versions/version_6.1.3.xml

mkdir -p output/${REACTION}_test
cd output/${REACTION}_test
hd_root --loadconfigs /work/halld2/home/boyu/src_analysis/data/configs/jana_${REACTION}.cfg ${INPUTFILE}