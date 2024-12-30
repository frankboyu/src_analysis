#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver01/root/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF/*90208_000.root
TREENAME=gd_pimprotmissprot__B4_F4_T1_S4_Tree
SELECTOR=piminus_p_recon
TAG=sim_2H_missprot

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

rm -r output/selectedtree_${SELECTOR}_${TAG}_backup.root
mv output/selectedtree_${SELECTOR}_${TAG}.root output/selectedtree_${SELECTOR}_${TAG}_backup.root

cd output/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../configs/DSelector_${SELECTOR}.C+'", '1')'