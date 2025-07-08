#!/bin/bash

REACTION=$1
PINDAYS=$2

DIR_LIST=()

if   [ "${REACTION}" == "phi_d_2H_ver10" ]
then
    DIR_LIST+=("/mss/halld/RunPeriod-2021-11/analysis/ver10/tree_gd_kpkmd__B4_F4/merged")
elif [ "${REACTION}" == "phi_d_2H_ver11" ]
then
    DIR_LIST+=("/mss/halld/RunPeriod-2021-11/analysis/ver11/tree_gd_kpkmd__B4_F4/merged")
    DIR_LIST+=("/mss/halld/RunPeriod-2021-11/analysis/ver11/tree_gd_kpkmdinc__B4_F4_T1_S2/merged")
elif [ "${REACTION}" == "rho_d_2H_ver10" ]
then
    DIR_LIST+=("/mss/halld/RunPeriod-2021-11/analysis/ver10/tree_gd_pippimd__B4_F4/merged")
elif [ "${REACTION}" == "rho_d_2H_ver11" ]
then
    DIR_LIST+=("/mss/halld/RunPeriod-2021-11/analysis/ver11/tree_gd_pippimd__B4_F4/merged")
    DIR_LIST+=("/mss/halld/RunPeriod-2021-11/analysis/ver11/tree_gd_pippimdinc__B4_F4_T1_S2/merged")
elif [ "${REACTION}" == "omega_d_2H_ver11" ]
then
    DIR_LIST+=("/mss/halld/RunPeriod-2021-11/analysis/ver11/tree_gd_pi0pippimd__B4_F4/merged")
    DIR_LIST+=("/mss/halld/RunPeriod-2021-11/analysis/ver11/tree_gd_pi0pippimdinc__B4_F4_T1_S2/merged")
fi

for DIR in "${DIR_LIST[@]}"
do
    jcache get ${DIR}/*.root -D ${PINDAYS}
done