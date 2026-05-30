#!/bin/bash

start=`date +%s`

# data
sh run_local.sh 'phi_d_exc_recon'   'data_ver12'     'both'

# sim
sh run_local.sh 'phi_d_exc_recon'   'sim_ver12'      'both'

# tagged
sh run_local.sh 'phi_d_exc_thrown'  'tagged_ver12'   'both'

# gen
sh run_local.sh 'phi_d_exc_thrown'  'gen_ver12'      'both'

# hists
# python configs/hist_phi_d_exc.py

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"