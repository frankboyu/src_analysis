#!/bin/bash

INPUTFILE=/cache/halld/RunPeriod-2021-11/analysis/ver06/tree_gd_kpkmprotmissn__B4_F4_T2_S5/merged/tree_gd_kpkmprotmissn__B4_F4_T2_S5_090207.root
TREE_NAME=gd_kpkmprotmissn__B4_F4_T2_S5_Tree
SELECTOR_FILE=phi_p_2H_data
TAG=3_test

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.17.0.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREE_NAME'", "'../../selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root flattree_${SELECTOR_FILE}_${TAG}.root