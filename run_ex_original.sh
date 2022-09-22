export ALGENCAN=/home/gustavo/algencan-3.1.1

rm -f ex_original

gfortran -O3 -w -fcheck=all -g ex_original.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o ex_original

./ex_original