#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver01/root/trees/*090262_000.root
TREENAME=gc12_pimprotinc__B4_F4_T2_S5_Tree
SELECTOR=piminus_p_recon
TAG=sim_2H_missprot

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

rm -r output/flattree_${SELECTOR}_${TAG}_backup.root
mv output/flattree_${SELECTOR}_${TAG}.root output/flattree_${SELECTOR}_${TAG}_backup.root

cd output/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../selectors/DSelector_${SELECTOR}.C+'", '1')'