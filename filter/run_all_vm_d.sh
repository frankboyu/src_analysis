#!/bin/bash

start=`date +%s`

# data
sh run_local.sh 'phi_d_recon_exc' 'data_2H'     'both'
sh run_local.sh 'rho_d_recon_exc' 'data_2H'     'both'

# sim
sh run_local.sh 'phi_d_recon_exc' 'sim_2H'      'both'
sh run_local.sh 'rho_d_recon_exc' 'sim_2H'      'both'

# tagged
sh run_local.sh 'phi_d_thrown' 'tagged_2H'  'both'
sh run_local.sh 'rho_d_thrown' 'tagged_2H'  'both'

# gen
sh run_local.sh 'phi_d_thrown' 'gen_2H'     'both'
sh run_local.sh 'rho_d_thrown' 'gen_2H'     'both'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"