#!/bin/bash

start=`date +%s`

# data
sh run_local.sh 'phi_d_recon' 'data_2H_exc'     'both'
# sh run_local.sh 'phi_d_recon' 'data_2H_inc'     'both'
# sh run_local.sh 'phi_d_recon' 'data_4He_IA'     'both'
# sh run_local.sh 'phi_d_recon' 'data_4He_inc'    'both'
# sh run_local.sh 'phi_d_recon' 'data_12C_IA'     'both'
# sh run_local.sh 'phi_d_recon' 'data_12C_inc'    'both'
sh run_local.sh 'rho_d_recon' 'data_2H_exc'     'both'
# sh run_local.sh 'rho_d_recon' 'data_2H_inc'     'both'
# sh run_local.sh 'rho_d_recon' 'data_4He_IA'     'both'
# sh run_local.sh 'rho_d_recon' 'data_4He_inc'    'both'
# sh run_local.sh 'rho_d_recon' 'data_12C_IA'     'both'
# sh run_local.sh 'rho_d_recon' 'data_12C_inc'    'both'

# sim
sh run_local.sh 'phi_d_recon' 'sim_2H_exc'      'both'
# sh run_local.sh 'phi_d_recon' 'sim_2H_inc'      'both'
# sh run_local.sh 'phi_d_recon' 'sim_4He_IA'      'both'
# sh run_local.sh 'phi_d_recon' 'sim_4He_inc'     'both'
# sh run_local.sh 'phi_d_recon' 'sim_12C_IA'      'both'
# sh run_local.sh 'phi_d_recon' 'sim_12C_inc'     'both'
sh run_local.sh 'rho_d_recon' 'sim_2H_exc'      'both'
# sh run_local.sh 'rho_d_recon' 'sim_2H_inc'      'both'
# sh run_local.sh 'rho_d_recon' 'sim_4He_IA'      'both'
# sh run_local.sh 'rho_d_recon' 'sim_4He_inc'     'both'
# sh run_local.sh 'rho_d_recon' 'sim_12C_IA'      'both'
# sh run_local.sh 'rho_d_recon' 'sim_12C_inc'     'both'

# tagged
sh run_local.sh 'phi_d_thrown' 'tagged_2H'  'both'
# sh run_local.sh 'phi_d_thrown' 'tagged_4He' 'both'
# sh run_local.sh 'phi_d_thrown' 'tagged_12C' 'both'
sh run_local.sh 'rho_d_thrown' 'tagged_2H'  'both'
# sh run_local.sh 'rho_d_thrown' 'tagged_4He' 'both'
# sh run_local.sh 'rho_d_thrown' 'tagged_12C' 'both'

# gen
sh run_local.sh 'phi_d_thrown' 'gen_2H'     'both'
# sh run_local.sh 'phi_d_thrown' 'gen_4He'    'both'
# sh run_local.sh 'phi_d_thrown' 'gen_12C'    'both'
sh run_local.sh 'rho_d_thrown' 'gen_2H'     'both'
# sh run_local.sh 'rho_d_thrown' 'gen_4He'    'both'
# sh run_local.sh 'rho_d_thrown' 'gen_12C'    'both'

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"