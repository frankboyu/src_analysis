#!/bin/bash

start=`date +%s`

INPUTFILE=/work/halld2/home/boyu/src_analysis/bggen/output/4He_p_3736/tree_ghe_pimprotinc__B4_F4_T1_S4/*.root
TREENAME=ghe_pimprotinc__B4_F4_T1_S4_Tree
SELECTOR=piminus_p_recon

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C+'", '1')'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"