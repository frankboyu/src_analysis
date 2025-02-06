## src_analysis_selection

Select the events with tighter cuts and flatten the trees, using DSelector scheme from gluex_root_analysis software

### Usage

1.  To run the DSelector over simulation or a small sample of data interactively on the command line

    `sh run_test.sh`

2.  To run the DSelector over the entire data set on the batch system

    `sh run_selection.sh REACTION`

    REACTION is one of the suffix in the configs folder

3.  To merge the output files

    `sh merge.sh REACTION`