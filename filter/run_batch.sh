#!/bin/bash

CHANNEL=$1
REACTION=$2
OUTPUTMODE=$3

JOB_WORKFLOW="-workflow src_analysis_filter"
JOB_NAME="-name ${CHANNEL}_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d-%H-%M')"
JOB_RESOURCES="-account halld -partition production -os el9 -cores 1 -ram 4GB -disk 8GB -time 24hrs"
JOB_OUT="-stdout /farm_out/boyu/src_analysis/filter/${CHANNEL}_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d').out"
JOB_ERR="-stderr /farm_out/boyu/src_analysis/filter/${CHANNEL}_${REACTION}_${OUTPUTMODE}_$(date '+%Y-%m-%d').err"
JOB_COMMAND="sh /work/halld2/home/boyu/src_analysis/filter/run_filter.sh ${CHANNEL} ${REACTION} ${OUTPUTMODE}"

swif2 add-job $JOB_WORKFLOW $JOB_NAME $JOB_RESOURCES $JOB_OUT $JOB_ERR $JOB_COMMAND