#!/bin/bash

echo "Start time"
date

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF/*.root
TREENAME=gd_pimprotmissprot__B4_F4_T1_S4_Tree
SELECTOR=piminus_p_recon

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C+'", '1')'

echo "Start time"
date