#!/bin/bash

start=`date +%s`

LABEL=case3_6.0
mkdir output/${LABEL}/

gfortran edved_wkng_pol.f

for EPHIN in 2.1
do
    for SPHIN in $(seq 0 0.5 80)
    do
        for BPHIN in $(seq 0 20 0.5)
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