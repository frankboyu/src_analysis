#!/bin/bash

REACTION=$1

tree_list=()

for file in output/${REACTION}/root/trees/*;
do
    tree_name=$(basename "$file")
    tree_name=${tree_name:0:${#tree_name}-16}  # remove the part of run number and file number, "_09XXXX_XXX.root"

    if [ ${#tree_list[@]} -eq 0 ] || [ ${tree_name} != ${tree_list[-1]} ];
    then
        tree_list+=(${tree_name})
    fi
done

for tree_name in ${tree_list[@]};
do
    echo "Moving files for ${REACTION} ${tree_name}"
    mkdir -p output/${REACTION}/root/${tree_name}
    mv output/${REACTION}/root/trees/${tree_name}_*.root output/${REACTION}/root/${tree_name}
done

rm -r output/${REACTION}/root/trees/