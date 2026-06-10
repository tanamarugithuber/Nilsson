#!/bin/bash

rm -f *.mod
rm -f *.exe
rm -f *.o

export OMP_NUM_THREADS=4

gfortran -fopenmp -c nilsson_basis_mod.f90

gfortran -fopenmp -c test.f90

gfortran -fopenmp \
    nilsson_basis_mod.o \
    test.o \
    -o test.exe