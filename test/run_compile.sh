#!/bin/bash

XML_VERSION=recon_srcct-2021_11-dev
SIM_VERSION=halld_sim_srcct

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv analysis-2021_11-dev.xml

cd /work/halld2/home/boyu/src_software_builds/halld_recon_srcct/src
scons install -j32