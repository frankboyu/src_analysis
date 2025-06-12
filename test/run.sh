#!/bin/bash
set -e

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/psharp/env/version_2022_09_28.xml

RUN=$1

INPUT=/volatile/halld/home/psharp/src_rho0_C12_1p_0T_inc/flattree_1p_0T/flattree_1p_0T_${RUN}.root

TREE=Tree_1p
OUTPUT=/volatile/halld/home/psharp/deltapp_Gary/filter1_C12_1p_0T_inc/trees/filter1_rho0_C12_1p_0T_inc${RUN}.root
OUTPUTHIST=/volatile/halld/home/psharp/deltapp_Gary/filter1_C12_1p_0T_inc/hists/filter1_rho0_C12_1p_0T_inc_hist_${RUN}.root

root -b -q '/work/halld2/home/psharp/rho-dyssey/scripts/filter1/Gary/filter_rho0_1_Gary.C+("'${INPUT}'","'${TREE}'","'${OUTPUT}'","'${OUTPUTHIST}'")'