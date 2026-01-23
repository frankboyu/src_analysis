#!/bin/bash

XML_VERSION=recon_srcct-2021_11-ver04_0.2
SIM_VERSION=halld_sim_srcct-5.1.0.2^hdr512

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/${XML_VERSION}.xml

cd /work/halld2/home/boyu/src_software_builds/${SIM_VERSION}/src
scons install -j32