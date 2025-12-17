#!/bin/bash

start=`date +%s`

LABEL=case5_4.5_5.0
mkdir output/${LABEL}/

gfortran edved_wkng_pol.f

for EPHIN in 8.5
do
    for SPHIN in {20..40..1}
    do
        for BPHIN in {5..25..1}
        do
            echo "Running with E=$EPHIN, sigma=$SPHIN mb, b=$BPHIN"
            echo "$EPHIN" > input.txt
            echo "$SPHIN" >> input.txt
            echo "$BPHIN" >> input.txt
            ./a.out > output/${LABEL}/E_${EPHIN}_sigma_${SPHIN}_b_${BPHIN}.txt
        done
    done
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"