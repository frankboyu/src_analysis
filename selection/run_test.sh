#!/bin/bash

start=`date +%s`

INPUTFILE=/work/halld2/home/oscarlin/analysis/rho_deuteron/simulation/output/model_piecewise/root/trees/tree_gd_pi0pimprotmissprot__B4_F4_T2_S2_gen_gcf_rhoMinus_model_noBKG_090221_000.root
# INPUTFILE=/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/tree_gd_kpkmd__B4_F4_gen_coherent/*90213_000.root
TREENAME=gd_pi0pimprotmissprot__B4_F4_T2_S2_Tree
SELECTOR=rhoMinus_p

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml
# gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C++'", '1')'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"