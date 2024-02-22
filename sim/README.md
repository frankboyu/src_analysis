## src_analysis_sim

Simulation suite for the SRC/CT analysis.

JLab: /work/halld2/home/boyu/src_analysis/sim

GitHub: https://github.com/frankboyu/src_analysis/tree/master/sim

### Usage

1.  Install the softwares in the builds directory, including halld_sim_srcct and gluex_MCwrapper_srcct.

2.  Recompile halld_sim_srcct if needed.

    `sh run_compile.sh VERSION`

    VERSION is the version number of the halld_sim_srcct software. Only the gen_MF generator is re-compiled for efficiency.

3.  Edit the config files (configs/gen_REACTION.cfg, jana_REACTION.cfg, wrapper_REACTION.cfg).

4.  Run the simulaiton

    `sh run_sim.sh REACTION VERSION EVENTS`

    REACTION can be one of the physics channel names or test.

    VERSION is the version number for the physics channel or the target name for the test files.

    EVENTS is the number of events to be run.

6.  Merge the output root files and back up all the files on tape.

    `sh run_merge.sh REACTION VERSION`

    Root files are merged by run number, while REST files are untouched. Test files are not applicable for this step.