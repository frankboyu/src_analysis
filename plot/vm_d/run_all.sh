#!/bin/bash

start=`date +%s`

REACTION_LIST=()
REACTION_LIST+=("recon_data_2H" "recon_sim_2H" "thrown_tagged_2H")

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml
for REACTION in "${REACTION_LIST[@]}"
do
    root -b -q -l "get_yield.C(\"$REACTION\")"
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"