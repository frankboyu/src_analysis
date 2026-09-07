#!/bin/bash

# ITERATION=$1

start=`date +%s`

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

for ITERATION in {0..10}
do
    echo "Iteration: ${ITERATION}"
    if [ ${ITERATION} == 0 ]; then
        root -b -q -l "get_yield.C(\"phi_d\", \"exc_recon_data_ver12\",       \"dsdt\", \"nominal\")"
    fi
    root -b -q -l     "get_yield.C(\"phi_d\", \"exc_recon_sim_ver12_07\",     \"dsdt\", \"simweight_iter${ITERATION}\")"
    root -b -q -l     "get_yield.C(\"phi_d\", \"exc_thrown_tagged_ver12_07\", \"dsdt\", \"simweight_iter${ITERATION}\")"
    output=$(python get_simweight.py ${ITERATION})
    echo "$output"
    if [ "$output" == "Converged!" ]; then
        break
    fi
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"