#!/bin/bash

CHANNEL=$1
VERSION=$2

source env.sh

echo "Making directory"
mkdir output_volatile/${CHANNEL}/${TAG}/

echo "Merging ROOT files"
if [ ${CHANNEL} == "piminus_p_2H_MF" ]
then
    hadd output_volatile/${CHANNEL}/${TAG}/genOut_gen_MF.root                               output_volatile/${CHANNEL}/root/generator/genOut_gen_MF_*.root
    hadd output_volatile/${CHANNEL}/${TAG}/tree_thrown_gen_MF.root                          output_volatile/${CHANNEL}/root/thrown/tree_thrown_gen_MF_*.root
    hadd output_volatile/${CHANNEL}/${TAG}/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF.root      output_volatile/${CHANNEL}/root/trees/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF_*.root
    hadd output_volatile/${CHANNEL}/${TAG}/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF.root output_volatile/${CHANNEL}/root/trees/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF_*.root  
elif [ ${CHANNEL} == "phi_p_2H_MF" ]
then
    hadd output_volatile/${CHANNEL}/${TAG}/genOut_gen_MF.root                               output_volatile/${CHANNEL}/root/generator/genOut_gen_MF_*.root
    hadd output_volatile/${CHANNEL}/${TAG}/tree_thrown_gen_MF.root                          output_volatile/${CHANNEL}/root/thrown/tree_thrown_gen_MF_*.root
    hadd output_volatile/${CHANNEL}/${TAG}/tree_gd_kpkmprotinc__B4_F4_T2_S5_gen_MF.root     output_volatile/${CHANNEL}/root/trees/tree_gd_kpkmprotinc__B4_F4_T2_S5_gen_MF_*.root
    hadd output_volatile/${CHANNEL}/${TAG}/tree_gd_kpkmprotmissn__B4_F4_T2_S5_gen_MF.root   output_volatile/${CHANNEL}/root/trees/tree_gd_kpkmprotmissn__B4_F4_T2_S5_gen_MF_*.root
fi

echo "Removing individual files"
rm -r output_volatile/${CHANNEL}/configurations
rm -r output_volatile/${CHANNEL}/hddm
rm -r output_volatile/${CHANNEL}/root