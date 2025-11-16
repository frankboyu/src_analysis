#!/bin/bash

start=`date +%s`

CHANNEL_LIST=()
CHANNEL_LIST+=("phi_d")
# CHANNEL_LIST+=("rho_d")

REACTION_LIST=()
# REACTION_LIST+=("recon_exc_data_2H" "recon_exc_sim_2H" "thrown_exc_tagged_2H")
REACTION_LIST+=("recon_exc_data_2H")

OBSERVABLE_LIST=()
OBSERVABLE_LIST+=("dsdt" "Wcostheta" "Wphi" "WPhi" "Wpsi")
# OBSERVABLE_LIST+=("dsdt")

TAG_LIST=()
# TAG_LIST+=("nominal")
# TAG_LIST+=("chisquared_3" "chisquared_4" "chisquared_6" "chisquared_7")
# TAG_LIST+=("momentum_0.35" "momentum_0.45" "theta_1.0" "theta_3.0")
TAG_LIST+=("vertexZ_13" "vertexZ_15" "vertexR_0.5" "vertexR_1.5")
# TAG_LIST+=("fitfunc_quadratic" "fitfunc_none")

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml
for OBSERVABLE in "${OBSERVABLE_LIST[@]}"
do
    for CHANNEL in "${CHANNEL_LIST[@]}"
    do
        for REACTION in "${REACTION_LIST[@]}"
        do
            for TAG in "${TAG_LIST[@]}"
            do
                if [[ "$REACTION" == *"thrown"* && "$TAG" != "nominal" ]]; then
                    continue
                fi
                if [[ "$REACTION" == *"sim"* && "$TAG" == *"fitfunc"* ]]; then
                    continue
                fi
                echo "Processing CHANNEL: $CHANNEL, REACTION: $REACTION, OBSERVABLE: $OBSERVABLE, TAG: $TAG"
                root -b -q -l "get_yield.C(\"$CHANNEL\", \"$REACTION\", \"$OBSERVABLE\", \"$TAG\")"
            done
        done
    done
done

# python get_plot.py

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"