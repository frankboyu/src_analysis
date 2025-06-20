#!/bin/bash

start=`date +%s`

CHANNEL_LIST=()
# CHANNEL_LIST+=("phi_d")
CHANNEL_LIST+=("rho_d")

REACTION_LIST=()
REACTION_LIST+=("recon_data_2H_exc" "recon_sim_2H_exc" "thrown_tagged_2H")

OBSERVABLE_LIST=()
# OBSERVABLE_LIST+=("dsdt" "Wcostheta" "Wphi" "WPhi" "Wpsi")
OBSERVABLE_LIST+=("Wcostheta")

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

python get_plot.py

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"