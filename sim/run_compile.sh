#!/bin/bash

XML_VERSION=recon_srcct-2021_11-ver04_0.5
SOFTWARE_VERSION=halld_sim_srcct-5.1.0.4^hdr512

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv /work/halld2/home/boyu/src_software_builds/halld_versions_srcct/${XML_VERSION}.xml

if [[ $SOFTWARE_VERSION == *"halld_sim"* ]]; then
    cd /work/halld2/home/boyu/src_software_builds/${SOFTWARE_VERSION}/src
    scons install -j32
elif [[ $SOFTWARE_VERSION == *"HDGeant4"* ]]; then
    cd /work/halld2/home/boyu/src_software_builds/${SOFTWARE_VERSION}/build
    cmake -DCMAKE_INSTALL_PREFIX=/work/halld2/home/boyu/src_software_builds/${SOFTWARE_VERSION}/Linux_Alma9-x86_64-gcc11.5.0 -DGeant4_DIR=$G4ROOT/lib64/Geant4-10.7.4 -DDiracxx_DIR=$DIRACXX_DIR/lib/cmake/Diracxx -DJANA_DIR=$JANA_HOME/lib/JANA/cmake ..
    cmake --build .
    cmake --install .
fi