#!/bin/bash

start=`date +%s`

INPUTFILE=/work/halld2/home/boyu/src_analysis/data/output/phi_d_2H_ver10/*90213.root
TREENAME=gd_kpkmd__B4_F4_Tree
SELECTOR=phi_d_recon

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml
# gxenv $HALLD_VERSIONS/version.xml

cd output/test/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../../configs/DSelector_${SELECTOR}.C++'", '1')'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"