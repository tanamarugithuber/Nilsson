#!/bin/bash

rm -f *.mod
rm -f *.exe
rm -f *.o

export OMP_NUM_THREADS=4

gfortran -c nilsson_basis_mod.f90
gfortran -c nilsson_matrix_mod.f90
gfortran -c test.f90

gfortran \
    nilsson_basis_mod.o \
    nilsson_matrix_mod.o \
    test.o \
    -llapack \
    -lblas \
    -o test.exe