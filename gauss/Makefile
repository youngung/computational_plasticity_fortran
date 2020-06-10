#FC=gfortran
FC=gfortran
FFLAGS= -Werror -O3

## objects
guass.o: gauss.f Makefile
main.o: main.f Makefile
main_mixed.o: main_mixed.f Makefile

## executables
main: gauss.o main.o Makefile
	$(FC) -o main gauss.o main.o
main_mixed: gauss.o main_mixed.o Makefile
	$(FC) -o main_mixed gauss.o main_mixed.o
