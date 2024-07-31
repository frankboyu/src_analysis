#!/bin/bash

REACTION=$1
TAG=$2

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

cd output/
mkdir ${REACTION}_${TAG}
tree_list=$(find . -type d -name "tree_${REACTION}_*")
for tree in ${tree_list};
do
    hadd ${REACTION}_${TAG}/${tree}.root ${tree}/**/*.root
    rm -r ${tree}
done