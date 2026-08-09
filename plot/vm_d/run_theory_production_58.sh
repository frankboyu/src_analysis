#!/bin/bash

start=`date +%s`

LABEL=four_para
mkdir output/${LABEL}/

gfortran -ffixed-line-length-none -o run_58.out edved_wkng_pol_58.f

for BEAMENERGY in 6.9
do
    mkdir output/${LABEL}/E_${BEAMENERGY}/
    for SGAMMAN in $(seq 13.8 0.2 14.2)
    do
        for BGAMMAN in $(seq 4 0.2 6)
        do
            mkdir output/${LABEL}/E_${BEAMENERGY}/E_${BEAMENERGY}_gamma_s_${SGAMMAN}_b_${BGAMMAN}/
            for SPHIN in $(seq 20 2.0 40)
            do
                for BPHIN in $(seq 5 1.0 15)
                do
                    echo "Running with E=$BEAMENERGY GeV, sigma_gn=$SGAMMAN mb, b_gn=$BGAMMAN GeV^-2, sigma_vn=$SPHIN mb, b_vn=$BPHIN GeV^-2"
                    echo "$BEAMENERGY" > input_paras_58.txt
                    echo "$SGAMMAN" >> input_paras_58.txt
                    echo "$BGAMMAN" >> input_paras_58.txt
                    echo "$SPHIN" >> input_paras_58.txt
                    echo "$BPHIN" >> input_paras_58.txt
                    ./run_58.out > output/${LABEL}/E_${BEAMENERGY}/E_${BEAMENERGY}_gamma_s_${SGAMMAN}_b_${BGAMMAN}/E_${BEAMENERGY}_gamma_s_${SGAMMAN}_b_${BGAMMAN}_phi_s_${SPHIN}_b_${BPHIN}.txt
                done
            done
        done
    done
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"