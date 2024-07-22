#!/bin/bash

REACTION=$1
TARGET=$2
EVENTS=$3

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv builds/halld_versions_srcct/recon_srcct.xml

if   [ ${TARGET} == "2H" ]
then
    gluex_MC.py configs/wrapper/wrapper_test.cfg 90213 ${EVENTS} per_file=1000000 batch=0
elif [ ${TARGET} == "4He" ]
then
    gluex_MC.py configs/wrapper/wrapper_test.cfg 90061 ${EVENTS} per_file=1000000 batch=0
elif [ ${TARGET} == "12C" ]
then
    gluex_MC.py configs/wrapper/wrapper_test.cfg 90290 ${EVENTS} per_file=1000000 batch=0
fi
