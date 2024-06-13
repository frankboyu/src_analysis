#!/bin/bash

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv halld_versions_srcct/recon_srcct-2021_11-ver01_2.xml

cd halld_sim_srcct/halld_sim_srcct-4.49.0^hdr438/src
scons install -j32

