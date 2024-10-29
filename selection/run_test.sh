#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/test/phi_d_2H_3_newcut/root/trees/tree_gd_kpkmmissd__B4_F4_T1_S4_gen_MF_090213_000.root
TREE_NAME=gd_kpkmmissd__B4_F4_T1_S4_Tree
SELECTOR_FILE=phi_c_2H_kpkmmissd
TAG=sim_test

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.17.0.xml

cd output/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../process_chain.C'("'$INPUTFILE'", "'$TREE_NAME'", "'../selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root flattree_${SELECTOR_FILE}_${TAG}.root