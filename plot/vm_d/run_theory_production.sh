#!/bin/bash

start=`date +%s`

LABEL=theory_two_para_nominal
mkdir output/${LABEL}/

gfortran -ffixed-line-length-none -o exe_theory.out get_edved_wkng_pol.f

for BEAMENERGY in 6.9 8.3 9.7
do
    for SGAMMAN in $(seq 11 0.1 11)
    do
        for BGAMMAN in $(seq 4.0 0.2 4.0)
        do
            for SPHIN in $(seq 20.0 0.2 40.0)
            do
                for BPHIN in $(seq 8.0 0.1 14.0)
                do
                    echo "Running with E=$BEAMENERGY GeV, sigma_gn=$SGAMMAN mb, b_gn=$BGAMMAN GeV^-2, sigma_vn=$SPHIN mb, b_vn=$BPHIN GeV^-2"
                    echo "$BEAMENERGY" > input/theory_paras.txt
                    echo "$SGAMMAN" >> input/theory_paras.txt
                    echo "$BGAMMAN" >> input/theory_paras.txt
                    echo "$SPHIN" >> input/theory_paras.txt
                    echo "$BPHIN" >> input/theory_paras.txt
                    ./exe_theory.out > output/${LABEL}/E_${BEAMENERGY}_gamma_s_${SGAMMAN}_b_${BGAMMAN}_phi_s_${SPHIN}_b_${BPHIN}.txt
                done
            done
        done
    done
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"