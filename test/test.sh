#!/bin/bash
#set -e

#i= 090206
# source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
# gxenv /work/halld2/home/psharp/env/version_2022_09_28.xml

#for i in $(cat /work/halld2/home/psharp/boat_your_rho/2_preselection/c12runnumbs.txt);

#cat /work/halld2/home/psharp/boat_your_rho/2_preselection/c12runnumbs.txt | while read i;
#for i in 090252

while read i
do
    echo $i
	# swif2 add-job filter1_C12_1p_0T_inc_Gary -account 'halld' -partition production -disk 1TB -ram 1GB -time 4hrs -cores 1 -shell /bin/sh -stdout /farm_out/psharp/filter1_C12_1p_0T_inc_Gary/${i}_out.txt -stderr /farm_out/psharp/filter1_C12_1p_0T_inc_Gary/${i}_err.txt sh /work/halld2/home/psharp/rho-dyssey/scripts/filter1/Gary/swif_scripts/run_filter1_C12_1p_0T_Gary.sh $i
done < /work/halld2/home/psharp/howto_ifarm/c12runnumbs.txt

echo 'done'