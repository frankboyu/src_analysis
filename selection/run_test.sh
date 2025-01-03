#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_test/root/trees/tree_gc12_pimprotinc*.root
TREENAME=gc12_pimprotinc__B4_F4_T1_S4_Tree
SELECTOR=piminus_p_recon

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C+'", '1')'