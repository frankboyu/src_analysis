## src_analysis_sim

Simulation suite for the SRC/CT analysis.

### Usage

1.  To test the simulation or run small sample size of simulation interactively on the command line

    `sh run_test_local.sh`

    If any of the software is not installed on the current system, run the singularity container first by

    `sh run_singularity.sh`

    If the simulation softwares need to be compiled, run

    `sh run_compile.sh`

2.  To test the simulation or run small sample size of simulation on the batch system

    `sh run_test_batch.sh`

3.  Run simulaiton of large sample size on the batch system

    `sh run_sim.sh REACTION EVENTS`

4.  To organize the output trees for further processing

    `sh run_postprocess.sh`
