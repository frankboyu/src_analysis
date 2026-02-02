#!/bin/bash

start=`date +%s`

# data
sh run_local.sh     'phi_d_recon'   '/work/halld2/home/boyu/src_analysis/data/output/phi_d_2H_ver12/tree_gd_kpkmd__B5_F4/*.root'        'gd_kpkmd__B5_F4_Tree'
sh run_rename.sh    'phi_d_recon'   'exc_data_2H' 'ver12'

# sim
sh run_local.sh     'phi_d_recon'   '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_01/root/tree_gd_kpkmd__B5_F4_gen_coherent/*.root'     'gd_kpkmd__B5_F4_Tree'
sh run_rename.sh    'phi_d_recon'   'exc_sim_2H' 'ver12_flat'
sh run_local.sh     'phi_d_recon'   '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_02/root/tree_gd_kpkmd__B5_F4_gen_coherent/*.root'     'gd_kpkmd__B5_F4_Tree'
sh run_rename.sh    'phi_d_recon'   'exc_sim_2H' 'ver12_model'

# gen
hadd output/selectedtree_phi_d_thrown.root  /work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_01/root/generator/*.root
sh run_rename.sh    'phi_d_thrown'  '' 'gen_2H_ver12_flat'
hadd output/selectedtree_phi_d_thrown.root  /work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_02/root/generator/*.root
sh run_rename.sh    'phi_d_thrown'  '' 'gen_2H_ver12_model'

# tagged
sh run_local.sh     'phi_d_thrown'  '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_01/root/thrown/*.root'  'Thrown_Tree'
sh run_rename.sh    'phi_d_thrown'  '' 'tagged_2H_ver12_flat'
sh run_local.sh     'phi_d_thrown'  '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_02/root/thrown/*.root'  'Thrown_Tree'
sh run_rename.sh    'phi_d_thrown'  '' 'tagged_2H_ver12_model'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"