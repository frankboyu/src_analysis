#!/bin/bash

REACTION=$1
TAG=$2
SUBTAG=$3

if [[ $REACTION == *"recon"* ]]; then
    mv output/selectedhist_${REACTION}_${TAG}.root output/selectedhist_${REACTION}_${TAG}_${SUBTAG}.root
    mv output/selectedtree_${REACTION}_${TAG}.root output/selectedtree_${REACTION}_${TAG}_${SUBTAG}.root
elif [[ $REACTION == *"thrown"* ]]; then
    mv output/selectedtree_${REACTION}.root output/selectedtree_${REACTION}_${SUBTAG}.root
else
    echo "Invalid REACTION type. Please use 'recon' or 'thrown'."
    exit 1
fi