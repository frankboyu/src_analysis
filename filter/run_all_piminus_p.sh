#!/bin/bash
RUNMODE=$1
INPUTMODE=$2
OUTPUTMODE=$3

start=`date +%s`

REACTION_LIST=()
# REACTION_LIST+=("data_2H_inc"      "data_2H_missprot"       "data_4He_inc"      "data_4He_misshe3"      "data_12C_inc"      "data_12C_missb11")
# REACTION_LIST+=("sim_2H_inc_flat"  "sim_2H_missprot_flat"   "sim_4He_inc_flat"  "sim_4He_misshe3_flat"  "sim_12C_inc_flat"  "sim_12C_missb11_flat")
# REACTION_LIST+=("sim_2H_inc_model" "sim_2H_missprot_model"  "sim_4He_inc_model" "sim_4He_misshe3_model" "sim_12C_inc_model" "sim_12C_missb11_model")
REACTION_LIST+=("bggen_4He_n_inc"   "bggen_4He_n_misshe3"   "bggen_4He_p_inc"   "bggen_4He_p_misshe3")
REACTION_LIST+=("bggen_12C_n_inc"   "bggen_12C_n_missb11"   "bggen_12C_p_inc"   "bggen_12C_p_missb11")

if [ "$RUNMODE" == "local" ]; then
    for REACTION in "${REACTION_LIST[@]}"
    do
        sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh piminus_p_recon $REACTION $INPUTMODE $OUTPUTMODE
    done
elif [ "$RUNMODE" == "batch" ]; then
    for REACTION in "${REACTION_LIST[@]}"
    do
        JOB_WORKFLOW="-workflow src_analysis_filter"
        JOB_NAME="-name piminus_p_recon_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d-%H-%M')"
        JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -ram 4GB -disk 8GB -time 24hrs"
        JOB_OUT="-stdout /farm_out/boyu/src_analysis/filter/piminus_p_recon_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d').out"
        JOB_ERR="-stderr /farm_out/boyu/src_analysis/filter/piminus_p_recon_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d').err"
        JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh piminus_p_recon $REACTION $INPUTMODE $OUTPUTMODE"
        swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND
    done
    swif2 run -workflow src_analysis_filter
fi

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"