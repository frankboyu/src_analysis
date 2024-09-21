#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/phi_p_2H_ver01/root/trees/tree_gd_kpkmprotmissn__B4_F4_T2_S5_gen_MF_090208_000.root
TREE_NAME=gd_kpkmprotmissn__B4_F4_T2_S5_Tree
SELECTOR_FILE=phi_p_2H_data
TAG=4_simtest

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.17.0.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREE_NAME'", "'../../selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root flattree_${SELECTOR_FILE}_${TAG}.root