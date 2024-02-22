#!/bin/bash

source env.sh

echo $CCDB_USER
echo $CCDB_CONNECTION
echo $JANA_CALIB_URL

while read RUN;
do
    ccdb dump /PHOTON_BEAM/endpoint_energy:${RUN}:default -> files/endpoint_${RUN}.txt
    
    ccdb add /PHOTON_BEAM/endpoint_energy -v mc -r ${RUN}-${RUN} files/endpoint_${RUN}.txt
done < list.txt