#!/bin/bash

REACTION=$1

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml

if [[ $REACTION == *"data"* || $REACTION == *"bggen"*]]; then
    hadd -f output/selectedhist_${REACTION}.root output/selectedhist_${REACTION}/*.root
    hadd -f output/selectedtree_${REACTION}.root output/selectedtree_${REACTION}/*.root
elif [[ $REACTION == *"sim"* ]]; then
    hadd -f output/selectedhist_${REACTION}.root output/selectedhist_${REACTION}/**/*.root
    hadd -f output/selectedtree_${REACTION}.root output/selectedtree_${REACTION}/**/*.root
fi

rm -r output/selectedhist_${REACTION}/
rm -r output/selectedtree_${REACTION}/