#!/bin/bash

REACTION=$1
TAG=$2

source env.sh

hadd output/flattree_${REACTION}_${TAG}.root output/flattree_${REACTION}/*.root

rm -r output/flattree_${REACTION}/
rm -r output/log/