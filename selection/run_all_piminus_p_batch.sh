#!/bin/bash

start=`date +%s`

# data
sh run_batch.sh piminus_p_recon_data_2H_inc
sh run_batch.sh piminus_p_recon_data_2H_missprot
sh run_batch.sh piminus_p_recon_data_4He_inc
sh run_batch.sh piminus_p_recon_data_4He_misshe3
sh run_batch.sh piminus_p_recon_data_12C_inc
sh run_batch.sh piminus_p_recon_data_12C_missb11

# sim
# sh run_batch.sh piminus_p_recon_sim_2H_inc_flat
# sh run_batch.sh piminus_p_recon_sim_2H_inc_model
# sh run_batch.sh piminus_p_recon_sim_2H_missprot_flat
# sh run_batch.sh piminus_p_recon_sim_2H_missprot_model
# sh run_batch.sh piminus_p_recon_sim_4He_inc_flat
# sh run_batch.sh piminus_p_recon_sim_4He_inc_model
# sh run_batch.sh piminus_p_recon_sim_4He_misshe3_flat
# sh run_batch.sh piminus_p_recon_sim_4He_misshe3_model
# sh run_batch.sh piminus_p_recon_sim_12C_inc_flat
# sh run_batch.sh piminus_p_recon_sim_12C_inc_model
# sh run_batch.sh piminus_p_recon_sim_12C_missb11_flat
# sh run_batch.sh piminus_p_recon_sim_12C_missb11_model

# thrown
# sh run_batch.sh piminus_p_thrown_tagged_2H_flat
# sh run_batch.sh piminus_p_thrown_tagged_2H_model
# sh run_batch.sh piminus_p_thrown_tagged_4He_flat
# sh run_batch.sh piminus_p_thrown_tagged_4He_model
# sh run_batch.sh piminus_p_thrown_tagged_12C_flat
# sh run_batch.sh piminus_p_thrown_tagged_12C_model

# bggen
# sh run_batch.sh piminus_p_recon_bggen_4He_n_inc
# sh run_batch.sh piminus_p_recon_bggen_4He_n_misshe3
# sh run_batch.sh piminus_p_recon_bggen_4He_p_inc
# sh run_batch.sh piminus_p_recon_bggen_4He_p_misshe3
# sh run_batch.sh piminus_p_recon_bggen_12C_n_inc
# sh run_batch.sh piminus_p_recon_bggen_12C_n_missb11
# sh run_batch.sh piminus_p_recon_bggen_12C_p_inc
# sh run_batch.sh piminus_p_recon_bggen_12C_p_missb11

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"