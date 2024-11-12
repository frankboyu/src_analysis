#!/bin/bash

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-ver01_4.xml

cd /work/halld2/home/boyu/src_software_builds/halld_sim_srcct-4.52.0^rec211111/src
scons install -j32