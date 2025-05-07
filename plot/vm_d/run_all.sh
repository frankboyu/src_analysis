#!/bin/bash

start=`date +%s`

CHANNEL_LIST=()
CHANNEL_LIST+=("phi_d")

REACTION_LIST=()
REACTION_LIST+=("recon_data_2H" "recon_sim_2H" "thrown_tagged_2H")

OBSERVABLE_LIST=()
OBSERVABLE_LIST+=("ds_dt" "W_costheta" "W_phi")

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml
for OBSERVABLE in "${OBSERVABLE_LIST[@]}"
do
    for CHANNEL in "${CHANNEL_LIST[@]}"
    do
        for REACTION in "${REACTION_LIST[@]}"
        do
            root -b -q -l "get_yield.C(\"$CHANNEL\", \"$REACTION\", \"$OBSERVABLE\")"
        done
    done
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"