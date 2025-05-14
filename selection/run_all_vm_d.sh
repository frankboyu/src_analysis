#!/bin/bash

MODE=$1

start=`date +%s`

if [[ $MODE == "run" ]]; then
    # data
    sh run_batch.sh 'phi_d_recon_data_2H_exc'
    # sh run_batch.sh 'phi_d_recon_data_2H_inc'
    # sh run_batch.sh 'phi_d_recon_data_4He_IA'
    # sh run_batch.sh 'phi_d_recon_data_4He_inc'
    # sh run_batch.sh 'phi_d_recon_data_12C_IA'
    # sh run_batch.sh 'phi_d_recon_data_12C_inc'
    sh run_batch.sh 'rho_d_recon_data_2H_exc'
    # sh run_batch.sh 'rho_d_recon_data_2H_inc'
    # sh run_batch.sh 'rho_d_recon_data_4He_IA'
    # sh run_batch.sh 'rho_d_recon_data_4He_inc'
    # sh run_batch.sh 'rho_d_recon_data_12C_IA'
    # sh run_batch.sh 'rho_d_recon_data_12C_inc'

    # sim
    sh run_local.sh     'phi_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/tree_gd_kpkmd__B4_F4_gen_coherent/*.root'             'gd_kpkmd__B4_F4_Tree'
    # sh run_local.sh     'phi_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/tree_gd_kpkmdinc__B4_F4_gen_coherent/*.root'          'gd_kpkmdinc__B4_F4_Tree'
    # sh run_local.sh     'phi_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_4He_ver01/root/tree_ghe_kpkmdmissd__B4_F4_gen_coherent/*.root'      'ghe_kpkmdmissd__B4_F4_Tree'
    # sh run_local.sh     'phi_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_4He_ver01/root/tree_ghe_kpkmdinc__B4_F4_gen_coherent/*.root'        'ghe_kpkmdinc__B4_F4_Tree'
    # sh run_local.sh     'phi_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_12C_ver01/root/tree_gc12_kpkmdmissb10__B4_F4_gen_coherent/*.root'   'gc12_kpkmdmissb10__B4_F4_Tree'
    # sh run_local.sh     'phi_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_12C_ver01/root/tree_gc12_kpkmdinc__B4_F4_gen_coherent/*.root'       'gc12_kpkmdinc__B4_F4_Tree'
    sh run_local.sh     'rho_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/tree_gd_pippimd__B4_F4_gen_coherent/*.root'           'gd_pippimd__B4_F4_Tree'
    # sh run_local.sh     'rho_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/tree_gd_pippimdinc__B4_F4_gen_coherent/*.root'        'gd_pippimdinc__B4_F4_Tree'
    # sh run_local.sh     'rho_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_4He_ver01/root/tree_ghe_pippimdmissd__B4_F4_gen_coherent/*.root'    'ghe_pippimdmissd__B4_F4_Tree'
    # sh run_local.sh     'rho_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_4He_ver01/root/tree_ghe_pippimdinc__B4_F4_gen_coherent/*.root'      'ghe_pippimdinc__B4_F4_Tree'
    # sh run_local.sh     'rho_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_12C_ver01/root/tree_gc12_pippimdmissb10__B4_F4_gen_coherent/*.root' 'gc12_pippimdmissb10__B4_F4_Tree'
    # sh run_local.sh     'rho_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_12C_ver01/root/tree_gc12_pippimdinc__B4_F4_gen_coherent/*.root'     'gc12_pippimdinc__B4_F4_Tree'

    # tagged
    sh run_local.sh     'phi_d_thrown' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/thrown/*.root'   'Thrown_Tree'
    sh run_rename.sh    'phi_d_thrown' ''  'tagged_2H'
    # sh run_local.sh     'phi_d_thrown' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_4He_ver01/root/thrown/*.root'  'Thrown_Tree'
    # sh run_rename.sh    'phi_d_thrown' ''  'tagged_4He'
    # sh run_local.sh     'phi_d_thrown' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_12C_ver01/root/thrown/*.root'  'Thrown_Tree'
    # sh run_rename.sh    'phi_d_thrown' ''  'tagged_12C'
    sh run_local.sh     'rho_d_thrown' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/thrown/*.root'   'Thrown_Tree'
    sh run_rename.sh    'rho_d_thrown' ''  'tagged_2H'
    # sh run_local.sh     'rho_d_thrown' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_4He_ver01/root/thrown/*.root'  'Thrown_Tree'
    # sh run_rename.sh    'rho_d_thrown' ''  'tagged_4He'
    # sh run_local.sh     'rho_d_thrown' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_12C_ver01/root/thrown/*.root'  'Thrown_Tree'
    # sh run_rename.sh    'rho_d_thrown' ''  'tagged_12C'

    # gen
    hadd output/selectedtree_phi_d_thrown_gen_2H.root   /work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/generator/*.root
    # hadd output/selectedtree_phi_d_thrown_gen_4He.root  /work/halld2/home/boyu/src_analysis/sim/output/phi_d_4He_ver01/root/generator/*.root
    # hadd output/selectedtree_phi_d_thrown_gen_12C.root  /work/halld2/home/boyu/src_analysis/sim/output/phi_d_12C_ver01/root/generator/*.root
    hadd output/selectedtree_rho_d_thrown_gen_2H.root   /work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/generator/*.root
    # hadd output/selectedtree_rho_d_thrown_gen_4He.root  /work/halld2/home/boyu/src_analysis/sim/output/rho_d_4He_ver01/root/generator/*.root
    # hadd output/selectedtree_rho_d_thrown_gen_12C.root  /work/halld2/home/boyu/src_analysis/sim/output/rho_d_12C_ver01/root/generator/*.root
elif [[ $MODE == "merge" ]]; then
    # data
    sh run_merge.sh 'phi_d_recon_data_2H_exc'
    # sh run_merge.sh 'phi_d_recon_data_2H_inc'
    # sh run_merge.sh 'phi_d_recon_data_4He_IA'
    # sh run_merge.sh 'phi_d_recon_data_4He_inc'
    # sh run_merge.sh 'phi_d_recon_data_12C_IA'
    # sh run_merge.sh 'phi_d_recon_data_12C_inc'
    sh run_merge.sh 'rho_d_recon_data_2H_exc'
    # sh run_merge.sh 'rho_d_recon_data_2H_inc'
    # sh run_merge.sh 'rho_d_recon_data_4He_IA'
    # sh run_merge.sh 'rho_d_recon_data_4He_inc'
    # sh run_merge.sh 'rho_d_recon_data_12C_IA'
    # sh run_merge.sh 'rho_d_recon_data_12C_inc'
else
    echo "Invalid mode. Please use 'run' or 'merge'."
    exit 1
fi

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"