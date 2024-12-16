#!/bin/bash

REACTION=piminus_p_2H
RUN=90213

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

python launch.py configs/jobs_analysis_${REACTION}.cfg ${RUN} ${RUN}
swif2 run -workflow src_analysis_skim