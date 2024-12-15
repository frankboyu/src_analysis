#!/bin/bash

XML_VERSION=dev
SIM_VERSION=halld_sim_srcct

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/recon_srcct-2021_11-${XML_VERSION}.xml

cd /work/halld2/home/boyu/src_software_builds/${SIM_VERSION}/src
scons install -j32