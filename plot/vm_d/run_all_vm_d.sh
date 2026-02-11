#!/bin/bash

start=`date +%s`

# RUN_MODE="batch"
RUN_MODE="local"
# RUN_MODE="echo"

CHANNEL_LIST=()
CHANNEL_LIST+=("phi_d")
# CHANNEL_LIST+=("rho_d")

REACTION_LIST=()
REACTION_LIST+=("recon_exc_data_2H_ver12")
REACTION_LIST+=("recon_exc_sim_2H_ver12_flat" "thrown_exc_tagged_2H_ver12_flat")
REACTION_LIST+=("recon_exc_sim_2H_ver12_model" "thrown_exc_tagged_2H_ver12_model")

OBSERVABLE_LIST=()
# OBSERVABLE_LIST+=("Wcostheta" "dsdt")
OBSERVABLE_LIST+=("dsdt" "Wcostheta" "Wdecayphi" "Wpolphi" "Wpsi")

TAG_LIST=()
TAG_LIST+=("nominal")
TAG_LIST+=("dEdx_1.50" "dEdx_1.75" "dEdx_2.50" "dEdx_3.00")
TAG_LIST+=("misspminus_0.0150" "misspminus_0.0175" "misspminus_0.0250" "misspminus_0.0300")
TAG_LIST+=("chisquared_4.50" "chisquared_4.75" "chisquared_5.50" "chisquared_6.00")
TAG_LIST+=("momentum_0.350" "momentum_0.375" "momentum_0.425" "momentum_0.450")
TAG_LIST+=("theta_1.90" "theta_1.95" "theta_2.05" "theta_2.10")
TAG_LIST+=("vertexZ_13.50" "vertexZ_13.75" "vertexZ_14.25" "vertexZ_14.50")
TAG_LIST+=("vertexR_0.50" "vertexR_0.75" "vertexR_1.25" "vertexR_1.50")
TAG_LIST+=("fitfunc_quadratic" "fitfunc_phenomenological" "fitfunc_fulllinear" "fitfunc_fullquadratic")
TAG_LIST+=("beamaccid_3" "beamaccid_5" "beamaccid_4out")
TAG_LIST+=("comboaccid_all" "comboaccid_none")
TAG_LIST+=("fitmax_1.06" "fitmax_1.07" "fitmax_1.09" "fitmax_1.10")
TAG_LIST+=("fitwidth_0.0040" "fitwidth_0.0048" "fitwidth_0.0060" "fitwidth_0.0075")
TAG_LIST+=("fitbkg_fulllinear" "fitbkg_quadratic" "fitbkg_fullquadratic" "fitbkg_phenomenological")
TAG_LIST+=("fitsig_noBL" "fitsig_nonrel" "fitsig_relBWsim")
# TAG_LIST+=("simweight_iter0" "simweight_iter1" "simweight_iter2" "simweight_iter3" "simweight_iter4" "simweight_iter5" "simweight_iter6")
# for i in {-1..1}
# do
#     for j in {-1..1}
#     do
#         for k in {-1..1}
#         do
#             for l in {-1..1}
#             do
#                 TAG_LIST+=("simweight_syst_a1_${i}_b1_${j}_a2_${k}_b2_${l}")
#             done
#         done
#     done
# done
# TAG_LIST+=("sideband")

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

for CHANNEL in "${CHANNEL_LIST[@]}"
do
    for REACTION in "${REACTION_LIST[@]}"
    do
        # root -b -q -l "get_num_combo.C(\"$CHANNEL\", \"$REACTION\")"

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
                if [[ "$REACTION" == *"data"* && "$TAG" == "sideband" ]]; then
                    continue
                fi
                if [[ "$OBSERVABLE" == "dsdt" && "$TAG" == "sideband" ]]; then
                    continue
                fi

                if [[ "$RUN_MODE" == "echo" ]]; then
                    echo "Dry run: CHANNEL=$CHANNEL, REACTION=$REACTION, OBSERVABLE=$OBSERVABLE, TAG=$TAG"
                elif [[ "$RUN_MODE" == "local" ]]; then
                    root -b -q -l "get_yield.C(\"$CHANNEL\", \"$REACTION\", \"$OBSERVABLE\", \"$TAG\")"
                elif [[ "$RUN_MODE" == "batch" ]]; then
                    JOB_WORKFLOW="-workflow src_analysis_plot"
                    JOB_NAME="-name yield_phi_d_${REACTION}_$(date '+%Y-%m-%d-%H-%M')"
                    JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -ram 1GB -disk 4GB -time 24hrs"
                    JOB_OUT="-stdout /farm_out/boyu/src_analysis/plot/yield_phi_d_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d').out"
                    JOB_ERR="-stderr /farm_out/boyu/src_analysis/plot/yield_phi_d_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d').err"
                    JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/plot/vm_d/run_yield_extraction.sh \"$CHANNEL\" \"$REACTION\" \"$OBSERVABLE\" \"$TAG\""
                    swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND
                else
                    echo "Error: Unknown RUN_MODE '$RUN_MODE'. Please set RUN_MODE to 'echo', 'local' or 'batch'."
                fi
            done
        done
    done
done

swif2 run src_analysis_plot
# python get_plot.py

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"