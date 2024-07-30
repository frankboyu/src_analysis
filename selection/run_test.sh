#!/bin/bash

INPUTFILE="/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver01/root/thrown/tree_thrown_gen_MF_090208_000.root"
# INPUTFILE="/cache/halld/RunPeriod-2021-11/analysis/ver06/tree_gd_kpkmprotmissn__B4_F4_T2_S5/merged/tree_gd_kpkmprotmissn__B4_F4_T2_S5_090213.root"
TREE_NAME="Thrown_Tree"
# TREE_NAME="gd_kpkmprotmissn__B4_F4_T2_S5_Tree"
SELECTOR_FILE="piminus_p_2H_thrown"
# NUM_THREADS=1
TAG="test"

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.17.0.xml

root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C process_tree.C'("'$INPUTFILE'", "'$TREE_NAME'", "'selectors/DSelector_${SELECTOR_FILE}.C+'", '1')'

mv flattree_${SELECTOR_FILE}.root output/flattree_${SELECTOR_FILE}_${TAG}.root