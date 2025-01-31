#!/bin/bash

INPUTFILE=/cache/halld/RunPeriod-2021-11/analysis/ver10/tree_ghe_kpkmdinc__B4_F4/merged/tree_ghe_kpkmdinc__B4_F4_090660.root
TREENAME=ghe_kpkmdinc__B4_F4_Tree
SELECTOR=phi_d_recon

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C+'", '1')'