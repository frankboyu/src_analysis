#!/bin/bash

start=`date +%s`

INPUTFILE=/cache/halld/RunPeriod-2021-11/analysis/ver06/tree_ghe_pippimprotinc__B4_F4_T2_S5/merged/tree_ghe_pippimprotinc__B4_F4_T2_S5_090661.root
TREENAME=ghe_pippimprotinc__B4_F4_T2_S5_Tree
SELECTOR=1p_0T_inc

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml
echo gxenv $HALLD_VERSIONS/version_5.21.0.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DS_${SELECTOR}.C++'", '1')'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"