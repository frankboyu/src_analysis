#!/bin/bash

start=`date +%s`

# data
sh run_local.sh     'phi_d_recon'   '/work/halld2/home/boyu/src_analysis/data/output/phi_d_2H_ver10/tree_gd_kpkmd__B4_F4/*.root'        'gd_kpkmd__B4_F4_Tree'
sh run_rename.sh    'phi_d_recon'   'data_2H_exc'  'jana1'
sh run_local.sh     'phi_d_recon'   '/work/halld2/home/boyu/src_analysis/data/output/phi_d_2H_ver11/tree_gd_kpkmd__B4_F4/*.root'        'gd_kpkmd__B4_F4_Tree'
sh run_rename.sh    'phi_d_recon'   'data_2H_exc'  'jana2'
sh run_local.sh     'rho_d_recon'   '/work/halld2/home/boyu/src_analysis/data/output/rho_d_2H_ver10/tree_gd_pippimd__B4_F4/*.root'      'gd_pippimd__B4_F4_Tree'
sh run_rename.sh    'rho_d_recon'   'data_2H_exc'  'jana1'
sh run_local.sh     'rho_d_recon'   '/work/halld2/home/boyu/src_analysis/data/output/rho_d_2H_ver11/tree_gd_pippimd__B4_F4/*.root'      'gd_pippimd__B4_F4_Tree'
sh run_rename.sh    'rho_d_recon'   'data_2H_exc'  'jana2'
sh run_local.sh     'omega_d_recon' '/work/halld2/home/boyu/src_analysis/data/output/omega_d_2H_ver11/tree_gd_pi0pippimd__B4_F4/*.root' 'gd_pi0pippimd__B4_F4_Tree'
sh run_rename.sh    'omega_d_recon' 'data_2H_exc'  'jana2'

# sim
sh run_local.sh     'phi_d_recon'   '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/tree_gd_kpkmd__B4_F4_gen_coherent/*.root'     'gd_kpkmd__B4_F4_Tree'
sh run_rename.sh    'phi_d_recon'   'sim_2H_exc'  'jana1_uniform'
sh run_local.sh     'phi_d_recon'   '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver02/root/tree_gd_kpkmd__B4_F4_gen_coherent/*.root'     'gd_kpkmd__B4_F4_Tree'
sh run_rename.sh    'phi_d_recon'   'sim_2H_exc'  'jana1_helicity'
sh run_local.sh     'rho_d_recon'   '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/tree_gd_pippimd__B4_F4_gen_coherent/*.root'   'gd_pippimd__B4_F4_Tree'
sh run_rename.sh    'rho_d_recon'   'sim_2H_exc'  'jana1_uniform'
sh run_local.sh     'rho_d_recon'   '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver02/root/tree_gd_pippimd__B4_F4_gen_coherent/*.root'   'gd_pippimd__B4_F4_Tree'
sh run_rename.sh    'rho_d_recon'   'sim_2H_exc'  'jana1_helicity'

# tagged
sh run_local.sh     'phi_d_thrown'  '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/thrown/*.root'  'Thrown_Tree'
sh run_rename.sh    'phi_d_thrown'  '' 'tagged_2H_jana1_uniform'
sh run_local.sh     'phi_d_thrown'  '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver02/root/thrown/*.root'  'Thrown_Tree'
sh run_rename.sh    'phi_d_thrown'  '' 'tagged_2H_jana1_helicity'
sh run_local.sh     'rho_d_thrown'  '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/thrown/*.root'  'Thrown_Tree'
sh run_rename.sh    'rho_d_thrown'  '' 'tagged_2H_jana1_uniform'
sh run_local.sh     'rho_d_thrown'  '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver02/root/thrown/*.root'  'Thrown_Tree'
sh run_rename.sh    'rho_d_thrown'  '' 'tagged_2H_jana1_helicity'

# gen
hadd output/selectedtree_phi_d_thrown.root  /work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/generator/*.root
sh run_rename.sh    'phi_d_thrown'  '' 'gen_2H_jana1_uniform'
hadd output/selectedtree_phi_d_thrown.root  /work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver02/root/generator/*.root
sh run_rename.sh    'phi_d_thrown'  '' 'gen_2H_jana1_helicity'
hadd output/selectedtree_rho_d_thrown.root  /work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/generator/*.root
sh run_rename.sh    'rho_d_thrown'  '' 'gen_2H_jana1_uniform'
hadd output/selectedtree_rho_d_thrown.root  /work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver02/root/generator/*.root
sh run_rename.sh    'rho_d_thrown'  '' 'gen_2H_jana1_helicity'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"