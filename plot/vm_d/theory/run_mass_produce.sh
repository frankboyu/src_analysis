#!/bin/bash

start=`date +%s`

for EPHIN in 7.0 8.5 10.0
do
    for SPHIN in {20..40..1}
    do
        for BPHIN in {5..25..1}
        do
            echo "Running with E=$EPHIN, sigma=$SPHIN mb, b=$BPHIN"
            echo "$EPHIN" > input.txt
            echo "$SPHIN" >> input.txt
            echo "$BPHIN" >> input.txt
            # ./a.out > output/theory_${EPHIN}_${SPHIN}_${BPHIN}.txt
        done
    done
done

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"