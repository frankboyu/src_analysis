#!/bin/bash

start=`date +%s`

RUN_MODE=$1

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

TAG_LIST=()
# TAG_LIST+=("two_para_nominal")
# TAG_LIST+=("two_para_sgamman_12.0_bgamman_5.0")
TAG_LIST+=("two_para_sgamman_12.0_bgamman_5.5")

sh run_theory_compile.sh

for TAG in "${TAG_LIST[@]}"
do
    if [[ "$RUN_MODE" == "echo" ]]; then
        echo "Dry run: TAG=$TAG"
    elif [[ "$RUN_MODE" == "local" ]]; then
        root -b -q -l "get_theory_run.C(\"$TAG\")"
    elif [[ "$RUN_MODE" == "batch" ]]; then
        JOB_WORKFLOW="-workflow src_analysis_plot"
        JOB_NAME="-name theory_${TAG}_$(date '+%Y-%m-%d-%H-%M')"
        JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -ram 1GB -disk 4GB -time 24hrs"
        JOB_OUT="-stdout /farm_out/boyu/src_analysis/plot/theory_${TAG}_$(date '+%Y-%m-%d').out"
        JOB_ERR="-stderr /farm_out/boyu/src_analysis/plot/theory_${TAG}_$(date '+%Y-%m-%d').err"
        JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/plot/vm_d/run_theory_single.sh \"$TAG\""
        swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND
    else
        echo "Error: Unknown RUN_MODE '$RUN_MODE'. Please set RUN_MODE to 'echo', 'local' or 'batch'."
    fi
done

if [[ "$RUN_MODE" == "batch" ]]; then
    swif2 run src_analysis_plot
fi

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"