FFLAGS=-Ofast -march=native -mtune=native -ffixed-line-length-none
#FFLAGS= -O3 -m 3 
precon=-I/usr/local/include -L/usr/local/lib -cpp
#FFLAGS=-O0 -g -Wall -ffixed-line-length-none -finit-real=nan -finit-integer=666 -fbounds-check -fbacktrace

F90=ftn

default: all

OBJS = system.o global.o params.o grid.o writeout.o allocation.o operate.o solver.o stochastic.o timeav.o


%.o: src/%.F90
	$(F90) $(precon) $(FFLAGS) -c $< -o src/$@

all: clean main


clean:
	rm -f *.mod $(addsuffix .mod, $(basename $(OBJS)))
	rm -f src/*.o $(OBJS) main.o
	rm -f sw.out

main: $(OBJS) main.o
	$(F90) $(FFLAGS) src/*.o -o sw.out

params.o: system.o
global.o: params.o

