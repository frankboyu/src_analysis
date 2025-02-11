#!/bin/bash

start=`date +%s`

INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_test/root/trees/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF_090213_000.root
TREENAME=gd_pimprotmissprot__B4_F4_T1_S4_Tree
SELECTOR=piminus_p_recon

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C+'", '1')'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"