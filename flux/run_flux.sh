#!/bin/bash

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

python3.6 get_flux.py 2H
python3.6 get_flux.py 4He
python3.6 get_flux.py 12C
python3.6 get_flux.py empty