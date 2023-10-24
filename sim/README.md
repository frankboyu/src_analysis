## src_analysis_sim

Simulation suite for the SRC/CT analysis.

JLab: /work/halld2/home/boyu/src_analysis/sim

GitHub:

### Usage

1.  Install the softwares in the build directory.

    `cd build`

    `git clone https"//github.com/JacksonPybus/halld_sim_srcct.git` Install the halld_sim_srcct package. (derivative from halld_sim with custom generator, private repository)

    `git clone https://github.com/JeffersonLab/gluex_MCwrapper.git` Install the gluex_MCwrapper package.

    `cd gluex_MCwrapper`

    `git checkout -b src-ct` Switch to the src-ct branch.

2.  Recompile halld_sim_srcct after installation, checkout, or personal edits.

    `sh run_compile.sh`

3.  Edit the config files (configs/CHANNEL/gen_CHANNEL.cfg, jana_CHANNEL.cfg, wrapper_CHANNEL.cfg) and the environment files (env.csh, env.sh).

4.  Test the simulation

     `sh run_test.sh CHANNEL EVENTS`
    
     It will run the simulation on the ifarm using wrapper_CHANNEL.cfg with specified number of events.

5.  Run the simulaiton

    `sh run_sim.sh CHANNEL EVENTS`

    which contains 3 steps:

    First, get_configs.py will generate the config files in the configs/ directory for each individual run.

    Then, get_number.py will generate a text files containing the number of events to run for each individual run based on the total luminosity.

    Finally, gluex_MC.py will submit the simulation jobs to the batch system.

6.  Merge the output root files and back up all the files on tape.

    `sh run_merge.sh CHANNEL TAG`

    where TAG is used to specify the configurations of the simulation.

    The output files of the simulation will be relocated to output_cache/CHANNEL/TAG. Root files are merged by run number, while REST files are untouched.