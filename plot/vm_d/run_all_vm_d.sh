#!/bin/bash

start=`date +%s`

CHANNEL_LIST=()
CHANNEL_LIST+=("phi_d")
# CHANNEL_LIST+=("rho_d")

REACTION_LIST=()
REACTION_LIST+=("recon_exc_data_2H" "recon_exc_sim_2H" "thrown_exc_tagged_2H")
# REACTION_LIST+=("recon_exc_data_2H")

OBSERVABLE_LIST=()
# OBSERVABLE_LIST+=("dsdt" "Wcostheta" "Wphi" "WPhi" "Wpsi")
OBSERVABLE_LIST+=("dsdt")

TAG_LIST=()
# TAG_LIST+=("nominal")
# TAG_LIST+=("dEdx_1.0" "dEdx_1.5" "dEdx_2.5" "dEdx_3.0")
# TAG_LIST+=("misspminus_0.010_chisquared_3.5" "misspminus_0.015_chisquared_3.5" "misspminus_0.020_chisquared_3.5" "misspminus_0.025_chisquared_3.5" "misspminus_0.030_chisquared_3.5")
# TAG_LIST+=("misspminus_0.010_chisquared_4.0" "misspminus_0.015_chisquared_4.0" "misspminus_0.020_chisquared_4.0" "misspminus_0.025_chisquared_4.0" "misspminus_0.030_chisquared_4.0")
# TAG_LIST+=("misspminus_0.010_chisquared_5.0" "misspminus_0.015_chisquared_5.0" "misspminus_0.020_chisquared_5.0" "misspminus_0.025_chisquared_5.0" "misspminus_0.030_chisquared_5.0")
# TAG_LIST+=("misspminus_0.010_chisquared_6.0" "misspminus_0.015_chisquared_6.0" "misspminus_0.020_chisquared_6.0" "misspminus_0.025_chisquared_6.0" "misspminus_0.030_chisquared_6.0")
# TAG_LIST+=("misspminus_0.010_chisquared_7.0" "misspminus_0.015_chisquared_7.0" "misspminus_0.020_chisquared_7.0" "misspminus_0.025_chisquared_7.0" "misspminus_0.030_chisquared_7.0")
# TAG_LIST+=("momentum_0.400" "momentum_0.425" "momentum_0.475" "momentum_0.500")
# TAG_LIST+=("theta_1.0" "theta_1.5" "theta_2.5" "theta_3.0")
# TAG_LIST+=("vertexZ_13.0" "vertexZ_13.5" "vertexZ_14.5" "vertexZ_15.0")
# TAG_LIST+=("vertexR_0.50" "vertexR_0.75" "vertexR_1.25" "vertexR_1.50")
# TAG_LIST+=("fitfunc_quadratic" "fitfunc_phenomenological" "fitfunc_fulllinear" "fitfunc_fullquadratic")
# TAG_LIST+=("beamaccid_3" "beamaccid_5" "beamaccid_4out")
# TAG_LIST+=("comboaccid_all" "comboaccid_none")
# TAG_LIST+=("fitmax_1.06" "fitmax_1.07" "fitmax_1.09" "fitmax_1.10")
# TAG_LIST+=("fitwidth_0.0040" "fitwidth_0.0048" "fitwidth_0.0060" "fitwidth_0.0075")
# TAG_LIST+=("fitbkg_fulllinear" "fitbkg_quadratic" "fitbkg_fullquadratic" "fitbkg_phenomenological")
# TAG_LIST+=("fitsig_noBL" "fitsig_nonrel" "fitsig_relBWsim")
# TAG_LIST+=("simweight_pass0" "simweight_pass1" "simweight_pass2" "simweight_pass3")
TAG_LIST+=("simweight_pass0")


source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

# root -b -q -l get_num_combo.C

for CHANNEL in "${CHANNEL_LIST[@]}"
do
    for REACTION in "${REACTION_LIST[@]}"
    do
        for OBSERVABLE in "${OBSERVABLE_LIST[@]}"
        do
            for TAG in "${TAG_LIST[@]}"
            do
                if [[ "$REACTION" == *"thrown"* && "$TAG" != "nominal" && "$TAG" != *"simweight"* ]]; then
                    continue
                fi
                if [[ "$REACTION" == *"sim"* && "$TAG" == *"fit"* ]]; then
                    continue
                fi
                if [[ "$REACTION" == *"data"* && "$TAG" == *"simweight"* ]]; then
                    continue
                fi
                # echo "Processing: CHANNEL=$CHANNEL, REACTION=$REACTION, OBSERVABLE=$OBSERVABLE, TAG=$TAG"
                root -b -q -l "get_yield.C(\"$CHANNEL\", \"$REACTION\", \"$OBSERVABLE\", \"$TAG\")"
            done
        done
    done
done

# python get_plot.py

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"