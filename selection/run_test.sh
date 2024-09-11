#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/skim/output/phi_c_4He_ver01/tree_ghe_kpkmmisshe__B4_F4_T1_S4/090656.root
TREE_NAME=ghe_kpkmmisshe__B4_F4_T1_S4_Tree
SELECTOR_FILE=phi_c_4He_data
TAG=090656

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.19.0.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREE_NAME'", "'../../selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root flattree_${SELECTOR_FILE}_${TAG}.root