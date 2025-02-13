#!/bin/bash

start=`date +%s`

CHANNEL=$1
REACTION=$2
INPUTMODE=$3
OUTPUTMODE=$4

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

root -b -q -l "/work/halld2/home/boyu/src_analysis/filter/configs/filter_$CHANNEL.C(\"$REACTION\", \"$INPUTMODE\", \"$OUTPUTMODE\")"

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"