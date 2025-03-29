#!/bin/bash
RUNMODE=$1
OUTPUTMODE=$2

start=`date +%s`

RECON_LIST=()
RECON_LIST+=("data_2H" "data_4He" "data_12C")
RECON_LIST+=("sim_2H" "sim_4He" "sim_12C")
THROWN_LIST=()
THROWN_LIST+=("tagged_2H" "tagged_4He" "tagged_12C")
THROWN_LIST+=("gen_2H" "gen_4He" "gen_12C")

if [ "$RUNMODE" == "local" ]; then
    for REACTION in "${RECON_LIST[@]}"
    do
        sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh phi_d_recon $REACTION $OUTPUTMODE
    done
elif [ "$RUNMODE" == "batch" ]; then
    for REACTION in "${RECON_LIST[@]}"
    do
        JOB_WORKFLOW="-workflow src_analysis_filter"
        JOB_NAME="-name phi_d_recon_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d-%H-%M')"
        JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -ram 4GB -disk 8GB -time 24hrs"
        JOB_OUT="-stdout /farm_out/boyu/src_analysis/filter/phi_d_recon_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d').out"
        JOB_ERR="-stderr /farm_out/boyu/src_analysis/filter/phi_d_recon_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d').err"
        JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh phi_d_recon $REACTION $OUTPUTMODE"
        swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND
    done
    swif2 run -workflow src_analysis_filter
fi

if [ "$RUNMODE" == "local" ]; then
    for REACTION in "${THROWN_LIST[@]}"
    do
        sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh phi_d_thrown $REACTION $OUTPUTMODE
    done
elif [ "$RUNMODE" == "batch" ]; then
    for REACTION in "${THROWN_LIST[@]}"
    do
        JOB_WORKFLOW="-workflow src_analysis_filter"
        JOB_NAME="-name phi_d_thrown_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d-%H-%M')"
        JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -ram 4GB -disk 8GB -time 24hrs"
        JOB_OUT="-stdout /farm_out/boyu/src_analysis/filter/phi_d_thrown_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d').out"
        JOB_ERR="-stderr /farm_out/boyu/src_analysis/filter/phi_d_thrown_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d').err"
        JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh phi_d_thrown $REACTION $OUTPUTMODE"
        swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND
    done
    swif2 run -workflow src_analysis_filter
fi

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"