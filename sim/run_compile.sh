#!/bin/bash

XML_VERSION=analysis-2021_11-ver10
SIM_VERSION=halld_recon_srcct-4.51.0

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/${XML_VERSION}.xml

cd /work/halld2/home/boyu/src_software_builds/${SIM_VERSION}/src
scons install -j32