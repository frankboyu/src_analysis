#!/bin/bash

INPUTFILE=/cache/halld/RunPeriod-2021-11/analysis/ver09/tree_gd_kpkmd__B4_F4_T0_S4/merged/*.root
TREE_NAME=gd_kpkmd__B4_F4_T0_S4_Tree
SELECTOR_FILE=phi_c_2H_kpkmd
TAG=2_all

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.17.0.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREE_NAME'", "'../../selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root flattree_${SELECTOR_FILE}_${TAG}.root