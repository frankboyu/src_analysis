#!/bin/bash

INPUTFILE=/work/halld2/home/boyu/src_analysis/skim/output/piminus_p_2H_test/tree_gd_pimprotmissprot__B4_F4_T1_S4/090213/*000.root
TREENAME=gd_pimprotmissprot__B4_F4_T1_S4_Tree
SELECTOR=piplus_n_recon

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C+'", '1')'