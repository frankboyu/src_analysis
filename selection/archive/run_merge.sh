#!/bin/bash

REACTION=$1
TAG=$2

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

hadd output/flattree_${REACTION}_${TAG}.root output/${REACTION}/**/*.root
rm -r output/${REACTION}