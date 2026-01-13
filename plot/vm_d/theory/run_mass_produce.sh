#!/bin/bash

start=`date +%s`

LABEL=case3_4.0_alpha0
mkdir output/${LABEL}/

gfortran edved_wkng_pol.f

for EPHIN in 2.1
do
    for SPHIN in $(seq 0 1 120)
    do
        for BPHIN in $(seq 0 1 30)
        do
            echo "Running with E=$EPHIN GeV, sigma=$SPHIN mb, b=$BPHIN GeV^-2"
            echo "$EPHIN" > input.txt
            echo "$SPHIN" >> input.txt
            echo "$BPHIN" >> input.txt
            ./a.out > output/${LABEL}/E_${EPHIN}_sigma_${SPHIN}_b_${BPHIN}.txt
        done
    done
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"