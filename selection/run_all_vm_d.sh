#!/bin/bash

start=`date +%s`

# data
sh run_local.sh     'phi_d_recon' '/work/halld2/home/boyu/src_analysis/data/output/phi_d_2H_ver10/*.root'   'gd_kpkmd__B4_F4_Tree'
sh run_local.sh     'rho_d_recon' '/work/halld2/home/boyu/src_analysis/data/output/rho_d_2H_ver10/*.root'   'gd_pippimd__B4_F4_Tree'

# sim
sh run_local.sh     'phi_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/tree_gd_kpkmd__B4_F4_gen_coherent/*.root'     'gd_kpkmd__B4_F4_Tree'
sh run_local.sh     'rho_d_recon' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/tree_gd_pippimd__B4_F4_gen_coherent/*.root'   'gd_pippimd__B4_F4_Tree'

# tagged
sh run_local.sh     'phi_d_thrown' '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/thrown/*.root'   'Thrown_Tree'
sh run_rename.sh    'phi_d_thrown' ''  'tagged_2H'
sh run_local.sh     'rho_d_thrown' '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/thrown/*.root'   'Thrown_Tree'
sh run_rename.sh    'rho_d_thrown' ''  'tagged_2H'

# gen
hadd output/selectedtree_phi_d_thrown_gen_2H.root   /work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/generator/*.root
hadd output/selectedtree_rho_d_thrown_gen_2H.root   /work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/generator/*.root

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"