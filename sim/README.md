## src_analysis_sim

Simulation suite for the SRC/CT analysis.

JLab: /work/halld2/home/boyu/src_analysis/sim

GitHub: https://github.com/frankboyu/src_analysis/tree/master/sim

### Usage

1.  To test the simulation or run small sample size of simulation interactively on the command line

    `sh run_test.sh`

    If any of the software is not installed on the current system, run the singularity container first by

    `sh run_singularity.sh`

2.  Run simulaiton of large sample size on the batch system

    `sh run_sim.sh REACTION EVENTS`
