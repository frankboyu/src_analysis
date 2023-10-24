#!/bin/bash

REACTION=$1
RUNMIN=$2
RUNMAX=$3

source env.sh
python launch.py jobs_data_${REACTION}.cfg ${RUNMIN} ${RUNMAX}
swif2 run -workflow src_analysis_data -maxconcurrent 1000