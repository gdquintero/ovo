rm -f main

gfortran -w -fcheck=all -g main.f90 sort.o subset.o -o main

./main
