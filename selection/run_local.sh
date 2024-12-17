#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver01/root/trees/tree_gd_pimprotinc*.root
TREENAME=gd_pimprotinc__B4_F4_T1_S4_Tree
SELECTOR=piminus_p_recon
TAG=sim_2H_inc

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

rm -r output/flattree_${SELECTOR}_${TAG}_backup.root
mv output/flattree_${SELECTOR}_${TAG}.root output/flattree_${SELECTOR}_${TAG}_backup.root

cd output/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../configs/DSelector_${SELECTOR}.C+'", '1')'