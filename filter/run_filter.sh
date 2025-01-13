#!/bin/bash
echo "Start time"
date

CHANNEL=$1
REACTION=$2
INPUTMODE=$3
OUTPUTMODE=$4

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

root -b -q -l "/work/halld2/home/boyu/src_analysis/filter/filters/filter_$CHANNEL.C(\"$REACTION\", \"$INPUTMODE\", \"$OUTPUTMODE\")"

echo "End time"
date