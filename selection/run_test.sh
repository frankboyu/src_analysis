#!/bin/bash

INPUTFILE=/cache/halld/RunPeriod-2021-11/analysis/ver05/tree_ghe_kpkmmisshe__B4_F4_T1_S4/merged/tree_ghe_kpkmmisshe__B4_F4_T1_S4_090061.root
TREE_NAME=ghe_kpkmmisshe__B4_F4_T1_S4_Tree
SELECTOR_FILE=phi_c_4He_data
TAG=1_initial

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.17.0.xml

root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C process_tree.C'("'$INPUTFILE'", "'$TREE_NAME'", "'selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root output/flattree_${SELECTOR_FILE}_${TAG}.root