#!/bin/bash

TAG=$1

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

root -b -q -l "/work/halld2/home/boyu/src_analysis/plot/vm_d/get_edved_wkng_pol_run.C(\"$TAG\")"