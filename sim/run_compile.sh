#!/bin/bash

XML_VERSION=ver01_4.1
SIM_VERSION=halld_sim_srcct-4.52.0.1^rec211111

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-${XML_VERSION}.xml

cd /work/halld2/home/boyu/src_software_builds/${SIM_VERSION}/src
scons install -j32