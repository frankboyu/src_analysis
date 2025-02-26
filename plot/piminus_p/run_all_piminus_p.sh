#!/bin/bash

start=`date +%s`

REACTION_LIST=()
# REACTION_LIST+=("recon_data_2H_missprot" "recon_data_4He_misshe3" "recon_data_12C_missb11")
REACTION_LIST+=("recon_sim_2H_missprot_model" "recon_sim_4He_misshe3_model" "recon_sim_12C_missb11_model")
# REACTION_LIST+=("thrown_gen_2H_model" "thrown_gen_4He_model" "thrown_gen_12C_model")
# REACTION_LIST+=("thrown_tagged_2H_model" "thrown_tagged_4He_model" "thrown_tagged_12C_model")

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml
for REACTION in "${REACTION_LIST[@]}"
do
    root -b -q -l "get_yield.C(\"$REACTION\")"
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"