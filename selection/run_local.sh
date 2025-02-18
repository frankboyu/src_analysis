#!/bin/bash

INPUTFILE=$1
TREENAME=$2
SELECTOR=$3

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version_5.21.0.xml

cd output/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../process_chain.C'("'$INPUTFILE'", "'$TREENAME'", "'../configs/DSelector_${SELECTOR}.C+'", '1')'