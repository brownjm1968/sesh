
gfortran -Wall -c mc-routines.f90

gfortran -Wall sesh.f90 mc-routines.o -o sesh