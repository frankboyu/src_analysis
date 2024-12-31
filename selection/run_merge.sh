#!/bin/bash

REACTION=$1

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

hadd output/selectedhist_${REACTION}.root output/selectedhist_${REACTION}/**/*.root
hadd output/selectedtree_${REACTION}.root output/selectedtree_${REACTION}/**/*.root

rm -r output/selectedhist_${REACTION}/
rm -r output/selectedtree_${REACTION}/