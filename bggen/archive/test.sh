#!/bin/bash

for run in output/4He_n_3739/hists/* ; do
    mv output/4He_n_3739/hists/${run: -6}/hd_root_${run: -6}_-1.root output/4He_n_3739/hists/hd_root_${run: -6}.root
    rm -r output/4He_n_3739/hists/${run: -6}
done