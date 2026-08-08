#!/bin/bash

start=`date +%s`

# data
sh run_local.sh     'jpsi_d_exc_recon'   '/work/halld2/home/boyu/src_analysis/data/output/jpsi_d_2H_ver13/tree_gd_epemd__B5_F4/*.root'        'gd_epemd__B5_F4_Tree'
sh run_rename.sh    'jpsi_d_exc_recon'   'data' 'ver13'

# sim
sh run_local.sh     'jpsi_d_exc_recon'   '/work/halld2/home/boyu/src_analysis/sim/output/jpsi_d_2H_ver13_01/root/tree_gd_epemd__B5_F4_gen_coherent/*.root'     'gd_epemd__B5_F4_Tree'
sh run_rename.sh    'jpsi_d_exc_recon'   'sim' 'ver13_01'

# gen
hadd output/selectedtree_jpsi_d_exc_thrown.root  /work/halld2/home/boyu/src_analysis/sim/output/jpsi_d_2H_ver13_01/root/generator/*.root
sh run_rename.sh    'jpsi_d_exc_thrown'  '' 'gen_ver13_01'

# tagged
sh run_local.sh     'jpsi_d_exc_thrown'  '/work/halld2/home/boyu/src_analysis/sim/output/jpsi_d_2H_ver13_01/root/thrown/*.root'  'Thrown_Tree'
sh run_rename.sh    'jpsi_d_exc_thrown'  '' 'tagged_ver13_01'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"