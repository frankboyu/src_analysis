#!/bin/bash

CHANNEL=$1
REACTION=$2
OBSERVABLE=$3
TAG=$4

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

root -b -q -l "get_yield.C(\"$CHANNEL\", \"$REACTION\", \"$OBSERVABLE\", \"$TAG\")"