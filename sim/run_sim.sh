#!/bin/bash

CHANNEL=$1
EVENTS=$2

source env.sh





python get_configs.py ${CHANNEL}
python get_number.py ${CHANNEL} ${EVENTS}

while read -r run number;
do
    echo ${run} ${number}
    gluex_MC.py configs/wrapper_${CHANNEL}_${run}.cfg ${run} ${number} per_file=1000000 batch=1
done < list.txt

swif2 run -workflow src_analysis_sim

rm list.txt
rm configs/${CHANNEL}/wrapper_${CHANNEL}_*.cfg








if   [ ${CHANNEL:(-11):2} == "2H" ]
then
    gluex_MC.py configs/wrapper_${CHANNEL}.cfg 90213 ${EVENTS} per_file=1000000 batch=0
elif [ ${CHANNEL:(-12):3} == "4He" ]
then
    gluex_MC.py configs/wrapper_${CHANNEL}.cfg 90061 ${EVENTS} per_file=1000000 batch=0
elif [ ${CHANNEL:(-12):3} == "12C" ]
then
    gluex_MC.py configs/wrapper_${CHANNEL}.cfg 90290 ${EVENTS} per_file=1000000 batch=0
fi

rm -r 90*