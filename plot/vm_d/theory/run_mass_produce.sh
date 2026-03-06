#!/bin/bash

start=`date +%s`

LABEL=nominal
mkdir output/${LABEL}/

gfortran edved_wkng_pol.f

for EPHIN in 8.3 6.9 9.7 2.1 3.1
do
    for SPHIN in $(seq 30 1 40)
    do
        for BPHIN in $(seq 10 0.25 15)
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