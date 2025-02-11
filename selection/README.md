## src_analysis_selection

Select the events with tighter cuts and flatten the trees, using DSelector scheme from gluex_root_analysis software

### Usage

1.  To test the DSelector interactively on the command line

    `sh run_test.sh`

2.  To run the DSelector locally

    `sh run_local.sh INPUTFILE TREENAME SELECTOR`

3.  To run the DSelector on the batch system

    `sh run_selection.sh REACTION`

4.  To merge the output files from step 3

    `sh run_merge.sh REACTION`