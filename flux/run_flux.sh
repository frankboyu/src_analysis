#!/bin/bash

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

python get_flux.py deuterium
python get_flux.py helium
python get_flux.py carbon
python get_flux.py empty