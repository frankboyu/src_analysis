#!/bin/bash

python select_runs.py "@is_src_production and @status_approved"

while read run;
do
  echo $run
  python2.7 plot_flux.py -r $run -e 6.0 -o output/run_$run.root
done < runs.txt

rm runs.txt




