#!/bin/bash

start=`date +%s`

RUN_MODE=$1

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

CHANNEL_LIST=()
CHANNEL_LIST+=("phi_d")

REACTION_LIST=()
# REACTION_LIST+=("exc_recon_data_ver12")
REACTION_LIST+=("exc_recon_sim_ver12_07")
REACTION_LIST+=("exc_thrown_tagged_ver12_07")

OBSERVABLE_LIST=()
# OBSERVABLE_LIST+=("dsdt")
OBSERVABLE_LIST+=("dsdt" "Wcostheta" "Wdecayphi" "Wpolphi" "Wsumpsi" "Wdiffpsi")

TAG_LIST=()
# TAG_LIST+=("nominal")
# TAG_LIST+=("dEdx_1.50"              "dEdx_1.75"                 "dEdx_2.50"             "dEdx_3.00")
# TAG_LIST+=("misspminus_0.0150"      "misspminus_0.0175"         "misspminus_0.0250"     "misspminus_0.0300")
# TAG_LIST+=("chisquared_4.50"        "chisquared_4.75"           "chisquared_5.50"       "chisquared_6.00")
# TAG_LIST+=("momentum_0.350"         "momentum_0.375"            "momentum_0.425"        "momentum_0.450")
# TAG_LIST+=("theta_1.90"             "theta_1.95"                "theta_2.05"            "theta_2.10")
# TAG_LIST+=("vertexZ_13.50"          "vertexZ_13.75"             "vertexZ_14.25"         "vertexZ_14.50")
# TAG_LIST+=("vertexR_0.50"           "vertexR_0.75"              "vertexR_1.25"          "vertexR_1.50")
# TAG_LIST+=("fitfunc_quadratic"      "fitfunc_phenomenological"  "fitfunc_fulllinear"    "fitfunc_fullquadratic")
# TAG_LIST+=("beamaccid_3"            "beamaccid_5"               "beamaccid_4out")
# TAG_LIST+=("comboaccid_all"         "comboaccid_none")
# TAG_LIST+=("fitmax_1.06"            "fitmax_1.07"               "fitmax_1.09"           "fitmax_1.10")
# TAG_LIST+=("fitwidth_0.0040"        "fitwidth_0.0048"           "fitwidth_0.0060"       "fitwidth_0.0075")
# TAG_LIST+=("fitbkg_fulllinear"      "fitbkg_quadratic"          "fitbkg_fullquadratic" "fitbkg_phenomenological")
# TAG_LIST+=("fitsig_noBL"            "fitsig_nonrel"             "fitsig_relBWsim")
# TAG_LIST+=("simweight_syst_a1_-1.0" "simweight_syst_a1_-0.5" "simweight_syst_a1_0.5" "simweight_syst_a1_1.0")
# TAG_LIST+=("simweight_syst_b1_-1.0" "simweight_syst_b1_-0.5" "simweight_syst_b1_0.5" "simweight_syst_b1_1.0")
# TAG_LIST+=("simweight_syst_a2_-1.0" "simweight_syst_a2_-0.5" "simweight_syst_a2_0.5" "simweight_syst_a2_1.0")
# TAG_LIST+=("simweight_syst_b2_-1.0" "simweight_syst_b2_-0.5" "simweight_syst_b2_0.5" "simweight_syst_b2_1.0")
TAG_LIST+=("rungroup_0_90" "rungroup_45_135" "rungroup_amo")
# TAG_LIST+=("sideband")

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
                if [[ "$REACTION" == *"sim"* && "$TAG" == *"fit"* && "$TAG" != "fitsig_relBWsim" ]]; then
                    continue
                fi
                if [[ "$REACTION" == *"data"* && "$TAG" == *"simweight"* ]]; then
                    continue
                fi
                if [[ "$REACTION" != *"data"* && "$TAG" == "sideband" ]]; then
                    continue
                fi
                if [[ "$OBSERVABLE" == "dsdt" && "$TAG" == "sideband" ]]; then
                    continue
                fi
                if [[ "$REACTION" == *"model"* && "$TAG" != "nominal" ]]; then
                    continue
                fi
                if [[ "$REACTION" == *"data"* && "$TAG" == "fitsig_relBWsim" ]]; then
                    continue
                fi

                if [[ "$RUN_MODE" == "echo" ]]; then
                    echo "Dry run: CHANNEL=$CHANNEL, REACTION=$REACTION, OBSERVABLE=$OBSERVABLE, TAG=$TAG"
                elif [[ "$RUN_MODE" == "local" ]]; then
                    root -b -q -l "get_yield.C(\"$CHANNEL\", \"$REACTION\", \"$OBSERVABLE\", \"$TAG\")"
                elif [[ "$RUN_MODE" == "batch" ]]; then
                    JOB_WORKFLOW="-workflow src_analysis_plot"
                    JOB_NAME="-name yield_${CHANNEL}_${REACTION}_${OBSERVABLE}_${TAG}_$(date '+%Y-%m-%d-%H-%M')"
                    JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -ram 1GB -disk 4GB -time 24hrs"
                    JOB_OUT="-stdout /farm_out/boyu/src_analysis/plot/yield_${CHANNEL}_${REACTION}_${OBSERVABLE}_${TAG}_$(date '+%Y-%m-%d').out"
                    JOB_ERR="-stderr /farm_out/boyu/src_analysis/plot/yield_${CHANNEL}_${REACTION}_${OBSERVABLE}_${TAG}_$(date '+%Y-%m-%d').err"
                    JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/plot/vm_d/run_yield_single.sh \"$CHANNEL\" \"$REACTION\" \"$OBSERVABLE\" \"$TAG\""
                    swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND
                else
                    echo "Error: Unknown RUN_MODE '$RUN_MODE'. Please set RUN_MODE to 'echo', 'local' or 'batch'."
                fi
            done
        done
    done
done

if [[ "$RUN_MODE" == "batch" ]]; then
    swif2 run src_analysis_plot
fi

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"