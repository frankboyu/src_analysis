#!/bin/bash

REACTION=$1

source env.sh
python launch.py configs/jobs_root_analysis_${REACTION}.cfg 90000 90662
swif2 run -workflow src_analysis_selection -maxconcurrent 1000