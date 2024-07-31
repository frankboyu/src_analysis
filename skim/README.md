## src_analysis_skim

Run ReactionFilter plugin of hd_root to skim the reconstructed data

JLab: /work/halld2/home/boyu/src_analysis/skim

GitHub: https://github.com/frankboyu/src_analysis/tree/master/skim

### Usage

1.  To run a test skim or run the skim on a small number of files interactively on the command line

    ```sh run_test.sh```

    If you want to run with an older version of software that is not compiled in the current environment, try to run inside the container

    ```sh run_sigularity.sh```

2.  To run the skim on one run

    ```sh run_skim.sh <REACTION> <RUN>```

### Note
1.  This part of code should be used for testing the reaction configs on a small number of files. To run the skim on the entire dataset, submit the reaction to the analysis launch.