#!/bin/bash

start=`date +%s`

# data
sh run_local.sh     'phi_d_exc_recon'   '/work/halld2/home/boyu/src_analysis/data/output/phi_d_2H_ver12/tree_gd_kpkmd__B5_F4/*.root'        'gd_kpkmd__B5_F4_Tree'
sh run_rename.sh    'phi_d_exc_recon'   'data' 'ver12'

# sim
sh run_local.sh     'phi_d_exc_recon'   '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_05/root/tree_gd_kpkmd__B5_F4_gen_coherent/*.root'     'gd_kpkmd__B5_F4_Tree'
sh run_rename.sh    'phi_d_exc_recon'   'sim' 'ver12'

# gen
hadd output/selectedtree_phi_d_exc_thrown.root  /work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_05/root/generator/*.root
sh run_rename.sh    'phi_d_exc_thrown'  '' 'gen_ver12'

# tagged
sh run_local.sh     'phi_d_exc_thrown'  '/work/halld2/home/boyu/src_analysis/sim/output/phi_d_2H_ver12_05/root/thrown/*.root'  'Thrown_Tree'
sh run_rename.sh    'phi_d_exc_thrown'  '' 'tagged_ver12'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"