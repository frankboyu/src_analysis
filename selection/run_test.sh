#!/bin/bash

# INPUTFILE="/cache/halld/RunPeriod-2021-11/analysis/ver05/tree_gd_kpkmmissd__B4_F4_T1_S4/merged/tree_gd_kpkmmissd__B4_F4_T1_S4_090213.root"
INPUTFILE="/cache/halld/RunPeriod-2021-11/analysis/ver06/tree_gd_kpkmprotmissn__B4_F4_T2_S5/merged/tree_gd_kpkmprotmissn__B4_F4_T2_S5_090213.root"
# TREE_NAME="gd_kpkmmissd__B4_F4_T1_S4_Tree"
TREE_NAME="gd_kpkmprotmissn__B4_F4_T2_S5_Tree"
SELECTOR_FILE="phi_c_2H_data"
NUM_THREADS=1
TAG="4_measured_incoherent"

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.17.0.xml

root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C process_tree.C'("'$INPUTFILE'", "'$TREE_NAME'", "'selectors/DSelector_${SELECTOR_FILE}.C+'", '${NUM_THREADS}')'

mv flattree_${SELECTOR_FILE}.root output/flattree_${SELECTOR_FILE}_${TAG}.root