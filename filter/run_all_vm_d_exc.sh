#!/bin/bash

start=`date +%s`

# data
sh run_local.sh 'phi_d_exc_recon'   'data_2H_ver12'     'both'

# sim
sh run_local.sh 'phi_d_exc_recon'   'sim_2H_ver12'      'both'

# tagged
sh run_local.sh 'phi_d_exc_thrown'  'tagged_2H_ver12'   'both'

# gen
sh run_local.sh 'phi_d_exc_thrown'  'gen_2H_ver12'      'both'

# hists
# python configs/hist_vm_d.py

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"