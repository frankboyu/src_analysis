## src_analysis_selection

Select the events with tighter cuts, using DSelector scheme from gluex_root_analysis software

JLab: /work/halld2/home/boyu/src_analysis/selection

GitHub: https://github.com/frankboyu/src_analysis/tree/master/selection

### Usage

1.  Edit the configuration and DSelector files

2.  Run the selection

    `sh run_selection.sh REACTION`

    REACTION is one of the suffix in the configs folder

3.  Merge the output flat trees

    `sh run_merge.sh REACTION TAG`

    TAG is the suffix added to the end of the output flat tree.