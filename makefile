#FFLAGS=-O2 -march=native -mtune=native -ffixed-line-length-none
#FFLAGS= -O3 -m 3 
precon=-I/usr/local/include -L/usr/local/lib -cpp
FFLAGS=-O0 -g -Wall -ffixed-line-length-none -finit-real=nan -finit-integer=666 -fbounds-check

F90=mpif90

default: all

OBJS = system.o global.o params.o grid.o operate.o solver.o writeout.o


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

