# C++ directives
CC=g++
CFLAGS=-std=c++11

# Fortran directives
FC=gfortran
FFLAGS=-std=f2003

all: clean
    $(FC) $(FFLAGS) -c -fcheck=all ConvertUnitsLib.f03
    $(CC) $(CFLAGS) -c CFort.cpp
    $(FC) $(FFLAGS) ConvertUnitsLib.o CFort.o -o convert

clean:
    rm -f *.o
    rm -f *.mod
