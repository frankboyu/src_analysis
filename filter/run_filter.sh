#!/bin/bash

start=`date +%s`

CHANNEL=$1
REACTION=$2
OUTPUTMODE=$3

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

root -b -q -l "/work/halld2/home/boyu/src_analysis/filter/configs/filter_$CHANNEL.C(\"$REACTION\", \"$OUTPUTMODE\")"

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"