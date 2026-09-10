#!/bin/sh
set -eu

# Build only the Fortran backend. ROOT macros load this shared library at runtime.
gfortran -O2 -fPIC -ffixed-line-length-none -shared \
    get_edved_wkng_pol_callable.f \
    -o exe_theory_callable.out

echo "Built $(pwd)/exe_theory_callable.out"
