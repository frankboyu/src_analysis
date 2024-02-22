#!/usr/bin/env bash

VERSION=$1

source env.sh

cd /work/halld2/home/boyu/src_analysis/sim/builds/halld_sim_srcct/halld_sim_srcct-${VERSION}/src/programs/Simulation/gen_MF
scons -u install