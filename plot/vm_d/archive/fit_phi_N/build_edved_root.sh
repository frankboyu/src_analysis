#!/bin/sh
set -eu

# Build only the Fortran backend. ROOT macros load this shared library at runtime.
gfortran -O2 -fPIC -ffixed-line-length-none -shared \
    get_edved_wkng_pol_callable.f \
    -o libedved.so

echo "Built $(pwd)/libedved.so"
