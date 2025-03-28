#!/bin/bash

start=`date +%s`

# # data
# sh run_batch.sh phi_d_recon_data_2H
# sh run_batch.sh phi_d_recon_data_4He
# sh run_batch.sh phi_d_recon_data_12C

# sim
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/tree_gd_kpkmd__B4_F4_gen_coherent/*.root' 'gd_kpkmd__B4_F4_Tree' 'phi_d_recon'
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_4He_ver01/root/tree_ghe_kpkmdinc__B4_F4_gen_coherent/*.root' 'ghe_kpkmdinc__B4_F4_Tree' 'phi_d_recon'
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_12C_ver01/root/tree_gc12_kpkmdinc__B4_F4_gen_coherent/*.root' 'gc12_kpkmdinc__B4_F4_Tree' 'phi_d_recon'

# thrown
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver01/root/thrown/*.root' 'Thrown_Tree' 'phi_d_thrown'
# mv output/selectedtree_phi_d_thrown.root output/selectedtree_phi_d_thrown_tagged_2H.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_4He_ver01/root/thrown/*.root' 'Thrown_Tree' 'phi_d_thrown'
mv output/selectedtree_phi_d_thrown.root output/selectedtree_phi_d_thrown_tagged_4He.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_12C_ver01/root/thrown/*.root' 'Thrown_Tree' 'phi_d_thrown'
mv output/selectedtree_phi_d_thrown.root output/selectedtree_phi_d_thrown_tagged_12C.root

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"