## src_analysis_skim

Run ReactionFilter plugin of hd_root to skim the reconstructed data

### Usage

1.  To run a test skim or run the skim on a small number of files interactively on the command line

    ```sh run_test_local.sh```

    If you want to run with an older version of software that is not compiled in the current environment, try to run inside the container first

    ```sh run_sigularity.sh```

2.  To run a test skim on one run

    ```sh run_test_batch.sh```

3.  To run through the entire dataset (not recommended, submit to the analysis launch instead)

    ```sh run_skim.sh <REACTION>```