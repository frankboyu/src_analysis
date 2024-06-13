#!/bin/bash

REACTION=$1
VERSION=$2

source env.sh

if   [ ${VERSION} == "2H" ] || [ ${REACTION:(-2)} == "2H" ]
then
    run_number=90213
elif [ ${VERSION} == "4He" ] || [ ${REACTION:(-3)} == "4He" ]
then
    run_number=90061
elif [ ${VERSION} == "12C" ] || [ ${REACTION:(-3)} == "12C" ]
then
    run_number=90290
fi

for tree_dir in /volatile/halld/home/boyu/src_analysis/data/single/tree_*;
do
    tree_name=$(basename "$tree_dir")
    echo ${tree_dir}
    hadd /volatile/halld/home/boyu/src_analysis/data/merged/${tree_name}_0${run_number}.root ${tree_dir}/0${run_number}/*.root
    rm -r ${tree_dir}
done