#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/skim/output/phi_c_2H_ver01/tree_gd_kpkmd__B4_F4_T0_S4/090213/*.root
TREE_NAME=gd_kpkmd__B4_F4_T0_S4_Tree
SELECTOR_FILE=phi_d_2H_data
TAG=1_initial

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.19.0.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREE_NAME'", "'../../selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root flattree_${SELECTOR_FILE}_${TAG}.root