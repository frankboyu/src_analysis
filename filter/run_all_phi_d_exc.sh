#!/bin/bash

start=`date +%s`

# data
sh run_local.sh 'phi_d_exc_recon'   'data_ver12'     'both'

# sim
sh run_local.sh 'phi_d_exc_recon'   'sim_ver12_07'      'both'

# tagged
sh run_local.sh 'phi_d_exc_thrown'  'tagged_ver12_07'   'both'

# gen
sh run_local.sh 'phi_d_exc_thrown'  'gen_ver12_07'      'both'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"