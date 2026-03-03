#!/bin/bash

start=`date +%s`

# data
sh run_local.sh 'phi_d_recon_exc' 'data_2H_ver12'           'hist'

# sim
sh run_local.sh 'phi_d_recon_exc' 'sim_2H_ver12_flat'       'hist'
sh run_local.sh 'phi_d_recon_exc' 'sim_2H_ver12_model'      'hist'

# tagged
sh run_local.sh 'phi_d_thrown_exc' 'tagged_2H_ver12_flat'   'hist'
sh run_local.sh 'phi_d_thrown_exc' 'tagged_2H_ver12_model'  'hist'

# gen
sh run_local.sh 'phi_d_thrown_exc' 'gen_2H_ver12_flat'      'hist'
sh run_local.sh 'phi_d_thrown_exc' 'gen_2H_ver12_model'     'hist'

# hists
# python configs/hist_vm_d.py

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"