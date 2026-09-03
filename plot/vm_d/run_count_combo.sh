#!/bin/bash

start=`date +%s`

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

root -b -q -l "get_num_combo.C(\"phi_d\", \"exc_recon_data_ver12\")"
root -b -q -l "get_num_combo.C(\"phi_d\", \"exc_recon_sim_ver12_07\")"

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"