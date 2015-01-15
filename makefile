#FFLAGS=-I/usr/local/include -L/usr/local/lib -cpp -O2 -march=native -mtune=native -ffixed-line-length-none
#FFLAGS= -O3 -m 3 
FFLAGS=-I/usr/local/include -L/usr/local/lib -cpp -O0 -g -Wall -ffixed-line-length-none -finit-real=nan -finit-integer=666 -fbounds-check

F90=mpif90

default: all

OBJS = system.o global.o params.o grid.o 


%.o: src/%.F90
	$(F90) $(FFLAGS) -c $< -o src/$@

all: clean main


clean:
	rm -f *.mod $(addsuffix .mod, $(basename $(OBJS)))
	rm -f src/*.o $(OBJS) main.o
	rm -f sw.out

main: $(OBJS) main.o
	$(F90) $(FFLAGS) src/*.o -o sw.out

params.o: system.o
global.o: params.o

