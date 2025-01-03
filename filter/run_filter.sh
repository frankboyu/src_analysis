#!/bin/bash
echo "Start time"
date

CHANNEL=$1
TAG=$2
INPUTMODE=$3
OUTPUTMODE=$4

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

root -b -q -l "filters/filter_$CHANNEL.C(\"$TAG\", \"$INPUTMODE\", \"$OUTPUTMODE\")"

echo "End time"
date