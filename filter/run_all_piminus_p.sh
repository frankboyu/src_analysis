#!/bin/bash
RUNMODE=$1
INPUTMODE=$2
OUTPUTMODE=$3

start=`date +%s`

RECON_LIST=()
# RECON_LIST+=("data_2H_inc"      "data_2H_missprot"      "data_4He_inc"      "data_4He_misshe3"      "data_12C_inc"      "data_12C_missb11")
# RECON_LIST+=("sim_2H_inc_flat"  "sim_2H_missprot_flat"  "sim_4He_inc_flat"  "sim_4He_misshe3_flat"  "sim_12C_inc_flat"  "sim_12C_missb11_flat")
# RECON_LIST+=("sim_2H_inc_model" "sim_2H_missprot_model" "sim_4He_inc_model" "sim_4He_misshe3_model" "sim_12C_inc_model" "sim_12C_missb11_model")
RECON_LIST+=("bggen_4He_inc"    "bggen_4He_misshe3"     "bggen_12C_inc"     "bggen_12C_missb11")

THROWN_LIST=()
# THROWN_LIST+=("gen_2H_model"    "gen_4He_model"     "gen_12C_model")
# THROWN_LIST+=("gen_2H_flat"     "gen_4He_flat"      "gen_12C_flat")
# THROWN_LIST+=("tagged_2H_model" "tagged_4He_model"  "tagged_12C_model")
# THROWN_LIST+=("tagged_2H_flat"  "tagged_4He_flat"   "tagged_12C_flat")

if [ "$RUNMODE" == "local" ]; then
    for REACTION in "${RECON_LIST[@]}"
    do
        sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh piminus_p_recon $REACTION $INPUTMODE $OUTPUTMODE
    done

    for REACTION in "${THROWN_LIST[@]}"
    do
        sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh piminus_p_thrown $REACTION $INPUTMODE $OUTPUTMODE
    done
elif [ "$RUNMODE" == "batch" ]; then
    for REACTION in "${RECON_LIST[@]}"
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

    for REACTION in "${THROWN_LIST[@]}"
    do
        JOB_WORKFLOW="-workflow src_analysis_filter"
        JOB_NAME="-name piminus_p_thrown_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d-%H-%M')"
        JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -ram 4GB -disk 8GB -time 24hrs"
        JOB_OUT="-stdout /farm_out/boyu/src_analysis/filter/piminus_p_thrown_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d').out"
        JOB_ERR="-stderr /farm_out/boyu/src_analysis/filter/piminus_p_thrown_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d').err"
        JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh piminus_p_thrown $REACTION $INPUTMODE $OUTPUTMODE"
        swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND
    done
    swif2 run -workflow src_analysis_filter
fi

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"