#!/bin/bash

echo "Initializing..."

# INPUT ARGUMENTS
REACTION=$1
VERSION=$2

# SET UP ENVIRONMENT
source env.sh

# GET THE TARGET NAME
if  [ ${REACTION:(-2)} == "2H" ]
then
    TARGET="deuterium"
elif [ ${REACTION:(-3)} == "4He" ]
then
    TARGET="helium"
elif [ ${REACTION:(-3)} == "12C" ]
then
    TARGET="carbon"
fi

# GET THE LIST OF TREES TO MERGE
tree_list=()
tree_list+=("genOut_gen_MF")       # output trees from generator
tree_list+=("tree_thrown_gen_MF")  # output trees from mcthrown_tree plugin

for file in /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/root/trees/*;  # output trees from reaction filter plugin
do
    tree_name=$(basename "$file")
    tree_name=${tree_name:0:${#tree_name}-16}  # remove the part of run number and file number, "_09XXXX_XXX.root"

    if [ ${tree_name} != ${tree_list[-1]} ];
    then
        tree_list+=(${tree_name})
    fi
done

# MAKE NEW DIRECTORIES
echo "Making new directories..."
mkdir /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/REST
for tree_name in ${tree_list[@]};
do
    mkdir /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/${tree_name}
    mkdir /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/${tree_name}/merged
done

# MOVE THE HDDM FILES
while read -r run flux;
do
    echo "Moving HDDM files for run ${run}..."
    mkdir /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/REST/0${run}
    cp /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/hddm/dana_rest_gen_MF_0${run}_*.hddm /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/REST/0${run}
    jcache put /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/REST/0${run}/*.hddm
done < /work/halld2/home/boyu/src_analysis/flux/output/${TARGET}/flux_total_${TARGET}.txt

# MERGE THE ROOT TREES
for tree_name in ${tree_list[@]};
do
    echo "Merging ${tree_name}.root..."
    if   [ ${tree_name} == "genOut_gen_MF" ]
    then
        hadd /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/${tree_name}/merged/${tree_name}_090000.root /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/root/generator/${tree_name}_*.root
        jcache put /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/${tree_name}/merged/${tree_name}_090000.root
    elif [ ${tree_name} == "tree_thrown_gen_MF" ]
    then
        hadd /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/${tree_name}/merged/${tree_name}_090000.root /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/root/thrown/${tree_name}_*.root
        jcache put /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/${tree_name}/merged/${tree_name}_090000.root
    else
        hadd /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/${tree_name}/merged/${tree_name}_090000.root /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/root/trees/${tree_name}_*.root
        jcache put /cache/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/${tree_name}/merged/${tree_name}_090000.root
    fi
done

# REMOVE THE INDIVIDUAL FILES
echo "Removing individual files..."
rm /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/hddm/*
rm /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/root/generator/*
rm /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/root/thrown/*
rm /volatile/halld/home/boyu/src_analysis/sim/${REACTION}/ver${VERSION}/root/trees/*
