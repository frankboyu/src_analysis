#!/bin/bash

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-ver01_1.xml

cd /work/halld2/home/boyu/src_software_builds/halld_sim_srcct-4.48.0^hdr438/src
scons install -j32

