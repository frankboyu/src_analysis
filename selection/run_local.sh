#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver01/root/trees/*090262_000.root
TREE_NAME=gc12_pimprotinc__B4_F4_T2_S5_Tree
SELECTOR_FILE=piminus_p_recon
TAG=test

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

cd output/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../process_chain.C'("'$INPUTFILE'", "'$TREE_NAME'", "'../selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root flattree_${SELECTOR_FILE}_${TAG}.root