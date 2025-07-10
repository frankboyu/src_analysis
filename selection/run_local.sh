#!/bin/bash

start=`date +%s`

REACTION=$1
FILENAME=$2
TREENAME=$3

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
if [[ $FILENAME == *"data"* ]]; then
    if [[ $FILENAME == *"ver10"* ]]; then
        echo gxenv $HALLD_VERSIONS/version_5.21.0.xml
        gxenv $HALLD_VERSIONS/version_5.21.0.xml
    elif [[ $FILENAME == *"ver11"* ]]; then
        echo gxenv $HALLD_VERSIONS/version_6.3.0.xml
        gxenv $HALLD_VERSIONS/version_6.3.0.xml
    fi
elif [[ $FILENAME == *"sim"* ]]; then
    if [[ $FILENAME == *"ver01"* || $FILENAME == *"ver02"* ]]; then
        echo gxenv $HALLD_VERSIONS/version_5.21.0.xml
        gxenv $HALLD_VERSIONS/version_5.21.0.xml
    elif [[ $FILENAME == *"ver03"* ]]; then
        echo gxenv $HALLD_VERSIONS/version_6.3.0.xml
        gxenv $HALLD_VERSIONS/version_6.3.0.xml
    fi
fi

cd output/
root -b -q $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C ../process_chain.C'("'$FILENAME'", "'$TREENAME'", "'../configs/DSelector_${REACTION}.C++'", '1')'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"