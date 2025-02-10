#!/bin/bash

echo "Start time"
date

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/thrown/*.root
TREENAME=Thrown_Tree
SELECTOR=piminus_p_thrown

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C+'", '1')'

echo "End time"
date