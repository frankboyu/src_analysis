#!/bin/bash

start=`date +%s`

INPUTFILE=/volatile/halld/home/psharp/sim_D_1p_0T_inc_AV18/root/trees/inc/tree_gd_pippimprotinc__B4_F4_T2_S4_gen_gcf_090561_000.root
TREENAME=gd_pippimprotinc__B4_F4_T2_S4_Tree
SELECTOR=1p_0T_inc

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml
# gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DS_${SELECTOR}.C++'", '1')'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"