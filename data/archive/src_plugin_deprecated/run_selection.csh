#!/bin/bash

while read run;
do
  echo $run
  ls /volatile/halld/home/sns/1p2pi_hddmver03/0${run}/ | wc -l
done < runs.dat
