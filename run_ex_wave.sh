export ALGENCAN=/home/gustavo/algencan-3.1.1

rm -f ex_zika

gfortran -O3 -w -fcheck=all -g main.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o zika

./zika
