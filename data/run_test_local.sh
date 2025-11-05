#!/bin/bash

INPUTFILE=/cache/halld/home/gxproj2/simulations/rp-2021-11-h2-27082025-rdm-bkg-flux-based-test-bggen-v2/bggen_upd-to-all_mesons/deuteron/neutron/rnb-90213/hddm/bggen_upd_090213_000_geant4_smeared.hddm
REACTION=phi_d_2H

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /group/halld/www/halldweb/html/halld_versions/version_6.5.0.xml

mkdir -p output/${REACTION}_test
cd output/${REACTION}_test
hd_root --loadconfigs /work/halld2/home/boyu/src_analysis/data/configs/jana_analysis_${REACTION}.cfg ${INPUTFILE}