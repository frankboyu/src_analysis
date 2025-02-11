#!/bin/bash

start=`date +%s`

# data
sh run_batch.sh phi_d_recon_data_2H
sh run_batch.sh phi_d_recon_data_4He
sh run_batch.sh phi_d_recon_data_12C

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"