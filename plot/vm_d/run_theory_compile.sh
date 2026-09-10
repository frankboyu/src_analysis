#!/bin/sh
set -eu

# Build only the Fortran backend. ROOT macros load this shared library at runtime.
gfortran -O2 -fPIC -ffixed-line-length-none -shared \
    get_theory_callable.f \
    -o lib_edved.so

echo "Built $(pwd)/lib_edved.so"
