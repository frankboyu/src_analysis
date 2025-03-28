#!/bin/bash
RUNMODE=$1
OUTPUTMODE=$2

start=`date +%s`

REACTION_LIST=()
REACTION_LIST+=("data_2H" "data_4He" "data_12C")

if [ "$RUNMODE" == "local" ]; then
    for REACTION in "${REACTION_LIST[@]}"
    do
        sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh phi_d_recon $REACTION $OUTPUTMODE
    done
elif [ "$RUNMODE" == "batch" ]; then
    for REACTION in "${REACTION_LIST[@]}"
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

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"