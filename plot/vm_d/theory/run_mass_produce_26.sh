#!/bin/bash

start=`date +%s`

LABEL=four_para
mkdir output/${LABEL}/

gfortran -ffixed-line-length-none -o run_26.out edved_wkng_pol.f

for BEAMENERGY in 3.1
do
    mkdir output/${LABEL}/E_${BEAMENERGY}/
    for SGAMMAN in $(seq 10 0.2 12)
    do
        for BGAMMAN in $(seq 2 0.2 6)
        do
            mkdir output/${LABEL}/E_${BEAMENERGY}/E_${BEAMENERGY}_gamma_s_${SGAMMAN}_b_${BGAMMAN}/
            for SPHIN in $(seq 20 2.0 40)
            do
                for BPHIN in $(seq 5 1.0 15)
                do
                    echo "Running with E=$BEAMENERGY GeV, sigma_gn=$SGAMMAN mb, b_gn=$BGAMMAN GeV^-2, sigma_vn=$SPHIN mb, b_vn=$BPHIN GeV^-2"
                    echo "$BEAMENERGY" > input_paras_26.txt
                    echo "$SGAMMAN" >> input_paras_26.txt
                    echo "$BGAMMAN" >> input_paras_26.txt
                    echo "$SPHIN" >> input_paras_26.txt
                    echo "$BPHIN" >> input_paras_26.txt
                    ./run_26.out > output/${LABEL}/E_${BEAMENERGY}/E_${BEAMENERGY}_gamma_s_${SGAMMAN}_b_${BGAMMAN}/E_${BEAMENERGY}_gamma_s_${SGAMMAN}_b_${BGAMMAN}_phi_s_${SPHIN}_b_${BPHIN}.txt
                done
            done
        done
    done
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"