#!/bin/bash

start=`date +%s`

# data
sh run_batch.sh phi_p_recon_data_2H_inc
sh run_batch.sh phi_p_recon_data_2H_missn
sh run_batch.sh phi_p_recon_data_4He_inc
sh run_batch.sh phi_p_recon_data_4He_misstri
sh run_batch.sh phi_p_recon_data_12C_inc
sh run_batch.sh phi_p_recon_data_12C_missb11

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"