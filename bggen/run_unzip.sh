#!/bin/bash

# for f in /cache/halld/gluex_simulations/REQUESTED_MC/bggen_upd-2021-11-nucleus-Helium-target-proton_3736/dana*
# do
#     tar -xvf $f
# done
# mv work/test-xrootd/gluex/mcwrap/REQUESTEDMC_OUTPUT/bggen_upd-2021-11-nucleus-Helium-target-proton_3736/hddm/ output/bggen_4He_p_3736/hddm/

# for f in /cache/halld/gluex_simulations/REQUESTED_MC/bggen_upd-2021-11-nucleus-Helium-target-neutron_3739/dana*
# do
#     tar -xvf $f
# done
# mv work/test-xrootd/gluex/mcwrap/REQUESTEDMC_OUTPUT/bggen_upd-2021-11-nucleus-Helium-target-neutron_3739/hddm/ output/bggen_4He_n_3739/hddm/

# for f in /cache/halld/gluex_simulations/REQUESTED_MC/bggen_upd-2021-11-nucleus-Carbon-target-proton_3738/dana*
# do
#     tar -xvf $f
# done
# mv work/test-xrootd/gluex/mcwrap/REQUESTEDMC_OUTPUT/bggen_upd-2021-11-nucleus-Carbon-target-proton_3738/hddm/ output/bggen_12C_p_3738/hddm/

for f in /cache/halld/gluex_simulations/REQUESTED_MC/bggen_upd-2021-11-nucleus-Carbon-target-neutron_3742/dana*
do
    tar -xvf $f
done
mv work/test-xrootd/gluex/mcwrap/REQUESTEDMC_OUTPUT/bggen_upd-2021-11-nucleus-Carbon-target-neutron_3742/hddm/ output/bggen_12C_n_3742/hddm/

rm -r work/