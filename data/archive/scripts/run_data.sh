#!/bin/bash

REACTION=$1
VERSION=$2

source env.sh

rm -r /volatile/halld/home/boyu/src_analysis/data/single
rm -r /volatile/halld/home/boyu/src_analysis/data/merged
mkdir /volatile/halld/home/boyu/src_analysis/data/single
mkdir /volatile/halld/home/boyu/src_analysis/data/merged

if  [ ${REACTION} == "test" ]
then
    if   [ ${VERSION} == "2H" ]
    then
        python launch.py configs/jobs_analysis_test.cfg 90213 90213
    elif [ ${VERSION} == "4He" ]
    then
        python launch.py configs/jobs_analysis_test.cfg 90061 90061
    elif [ ${VERSION} == "12C" ]
    then
        python launch.py configs/jobs_analysis_test.cfg 90290 90290
    fi
else
    if   [ ${REACTION:(-2)} == "2H" ]
    then
        python launch.py configs/jobs_analysis_${REACTION}_ver${VERSION}.cfg 90213 90213
    elif [ ${REACTION:(-3)} == "4He" ]
    then
        python launch.py configs/jobs_analysis_${REACTION}_ver${VERSION}.cfg 90061 90061
    elif [ ${REACTION:(-3)} == "12C" ]
    then
        python launch.py configs/jobs_analysis_${REACTION}_ver${VERSION}.cfg 90290 90290
    fi
fi

swif2 run -workflow src_analysis_data