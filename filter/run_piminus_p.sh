#!/bin/bash
RUNMODE=$1
INPUTMODE=$2
OUTPUTMODE=$3

echo "Start time"
date

REACTION_LIST=("data_2H_missprot")
# REACTION_LIST+=("data_2H_missprot" "data_2H_inc" "data_4He_inc" "data_12C_inc")
# REACTION_LIST+=("sim_2H_missprot_flat" "sim_2H_inc_flat" "sim_4He_inc_flat" "sim_12C_inc_flat")
# REACTION_LIST+=("sim_2H_missprot_model" "sim_2H_inc_model" "sim_4He_inc_model" "sim_12C_inc_model")

if [ "$RUNMODE" == "local" ]; then
    source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
    gxenv $HALLD_VERSIONS/version.xml
    for REACTION in "${REACTION_LIST[@]}"
    do
        root -b -q -l "filters/filter_piminus_p_recon.C(\"$REACTION\", \"$INPUTMODE\", \"$OUTPUTMODE\")"
    done
elif [ "$RUNMODE" == "batch" ]; then
    for REACTION in "${REACTION_LIST[@]}"
    do
        JOB_WORKFLOW="-workflow src_analysis_filter"
        JOB_NAME="-name piminus_p_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d-%H-%M')"
        JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -disk 8GB -ram 4GB -time 24hrs"
        JOB_OUT="-stdout /farm_out/boyu/src_analysis/filter/piminus_p_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d').out"
        JOB_ERR="-stderr /farm_out/boyu/src_analysis/filter/piminus_p_${REACTION}_${INPUTMODE}_${OUTPUTMODE}_$(date '+%Y-%m-%d').err"
        JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh piminus_p_recon $REACTION $INPUTMODE $OUTPUTMODE"
        swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND
    done
    swif2 run -workflow src_analysis_filter
fi

echo "End time"
date
