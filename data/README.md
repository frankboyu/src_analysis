## src_analysis_data

Use reaction filter plugin of hd_root to skim the REST files.

JLab: /work/halld2/home/boyu/src_analysis/data

GitHub: https://github.com/frankboyu/src_analysis/tree/master/data

### Usage

1. Edit the config files jana_data_REACTION.cfg and jobs_data_REACTION.cfg with the reactions to run and job resources

2. To run the skim:

`sh run_skim.sh REACTION RUNMIN RUNMAX`

### Notes

1. This is more intended for testing. To skim the entire dataset, submit the reactions to the offcial analysis launch.