#!/bin/bash

REACTION=$1
VERSION=$2
EVENTS=$3

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv builds/halld_versions_srcct/recon_srcct.xml

if  [ ${REACTION} == "test" ]
then
    rm -r /volatile/halld/home/boyu/src_analysis/sim/test/
    mkdir /volatile/halld/home/boyu/src_analysis/sim/test/

    if   [ ${VERSION} == "2H" ]
    then
        gluex_MC.py configs/wrapper/wrapper_test.cfg 90213 ${EVENTS} per_file=1000000 batch=2 logdir=/farm_out/boyu/src_analysis/sim
    elif [ ${VERSION} == "4He" ]
    then
        gluex_MC.py configs/wrapper/wrapper_test.cfg 90061 ${EVENTS} per_file=1000000 batch=2 logdir=/farm_out/boyu/src_analysis/sim
    elif [ ${VERSION} == "12C" ]
    then
        gluex_MC.py configs/wrapper/wrapper_test.cfg 90290 ${EVENTS} per_file=1000000 batch=2 logdir=/farm_out/boyu/src_analysis/sim
    fi
else
    mkdir /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/
    mkdir /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/

    python get_configs.py ${REACTION} ${VERSION}  # Generate the config files for each individual run.
    python get_number.py ${REACTION} ${EVENTS}    # Generate a text file of the number of events for each individual run

    while read -r run number;
    do
        echo ${run} ${number}
        gluex_MC.py configs/wrapper/wrapper_${REACTION}_ver${VERSION}_${run}.cfg ${run} ${number} per_file=1000000 batch=1  logdir=/farm_out/boyu/src_analysis/sim
    done < list.txt

    swif2 run -workflow src_analysis_sim

    rm list.txt
    rm configs/wrapper/wrapper_${REACTION}_ver${VERSION}_*.cfg
fi