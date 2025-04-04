#!/bin/bash

start=`date +%s`

# # data
# sh run_batch.sh rho_d_recon_data_2H
# sh run_batch.sh rho_d_recon_data_4He
# sh run_batch.sh rho_d_recon_data_12C

# sim
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/tree_gd_pippimd__B4_F4_gen_coherent/*.root' 'gd_pippimd__B4_F4_Tree' 'rho_d_recon'
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_4He_ver01/root/tree_ghe_pippimdinc__B4_F4_gen_coherent/*.root' 'ghe_pippimdinc__B4_F4_Tree' 'rho_d_recon'
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_12C_ver01/root/tree_gc12_pippimdinc__B4_F4_gen_coherent/*.root' 'gc12_pippimdinc__B4_F4_Tree' 'rho_d_recon'

# thrown
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/thrown/*.root' 'Thrown_Tree' 'rho_d_thrown'
mv output/selectedtree_rho_d_thrown.root output/selectedtree_rho_d_thrown_tagged_2H.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_4He_ver01/root/thrown/*.root' 'Thrown_Tree' 'rho_d_thrown'
# mv output/selectedtree_rho_d_thrown.root output/selectedtree_rho_d_thrown_tagged_4He.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/rho_d_12C_ver01/root/thrown/*.root' 'Thrown_Tree' 'rho_d_thrown'
# mv output/selectedtree_rho_d_thrown.root output/selectedtree_rho_d_thrown_tagged_12C.root

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"