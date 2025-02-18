#!/bin/bash

REACTION=$1

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml

if [[ $REACTION == *"data"* ]]; then
    hadd -f output/selectedhist_${REACTION}.root output/selectedhist_${REACTION}/*.root
    hadd -f output/selectedtree_${REACTION}_all.root output/selectedtree_${REACTION}/*.root
    if [[ $REACTION == *"2H"* ]]; then
        hadd -f output/selectedtree_${REACTION}_one.root output/selectedtree_${REACTION}/*090213.root
    elif [[ $REACTION == *"4He"* ]]; then
        hadd -f output/selectedtree_${REACTION}_one.root output/selectedtree_${REACTION}/*090061.root
    elif [[ $REACTION == *"12C"* ]]; then
        hadd -f output/selectedtree_${REACTION}_one.root output/selectedtree_${REACTION}/*090291.root
    fi
elif [[ $REACTION == *"sim"* ]]; then
    hadd -f output/selectedhist_${REACTION}.root output/selectedhist_${REACTION}/**/*.root
    hadd -f output/selectedtree_${REACTION}_all.root output/selectedtree_${REACTION}/**/*.root
    if [[ $REACTION == *"2H"* ]]; then
        hadd -f output/selectedtree_${REACTION}_all.root output/selectedtree_${REACTION}/090213/*.root
    elif [[ $REACTION == *"4He"* ]]; then
        hadd -f output/selectedtree_${REACTION}_all.root output/selectedtree_${REACTION}/090061/*.root
    elif [[ $REACTION == *"12C"* ]]; then
        hadd -f output/selectedtree_${REACTION}_all.root output/selectedtree_${REACTION}/090291/*.root
    fi
fi

rm -r output/selectedhist_${REACTION}/
rm -r output/selectedtree_${REACTION}/