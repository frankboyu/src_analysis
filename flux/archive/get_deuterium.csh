#!/bin/bash

python select_runs.py "@is_src_production and @status_approved and target_type=='FULL & Ready Deuterium'"

python2.7 calc_flux.py -f runs.txt -e 9.4276605900

rm runs.txt




