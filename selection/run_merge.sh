#!/bin/bash

REACTION=$1
TAG=$2

source env.sh

hadd output/${REACTION}/flattree_${REACTION}_${TAG}.root output/${REACTION}/flattree/*.root

rm -r output/${REACTION}/flattree/
rm -r output/${REACTION}/log/