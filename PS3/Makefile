ODIR=./src/obj
IDIR=./include
SDIR=./src
CFLAGS+= -std=c99 -g -O3 -I$(IDIR)
CC=gcc

_DEPS = bitmap.h CA.h RPS.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_UTIL_OBJ = CA.o bitmap.o
_PTHREAD_OBJ = RPS_pthread.o
_OMP_OBJ = RPS_omp.o
UTIL_OBJ = $(patsubst %,$(ODIR)/%,$(_UTIL_OBJ))
PTHREAD_OBJ = $(patsubst %,$(ODIR)/%,$(_PTHREAD_OBJ))
OMP_OBJ = $(patsubst %,$(ODIR)/%,$(_OMP_OBJ))

_SRC = bitmap.c RPS_pthread.c RPS_omp.c CA.c
SRC = $(patsubst %,$(SDIR)/%,$(_SRC))

$(ODIR)/%.o: ./src/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

#SHELL := /bin/bash

omp : $(UTIL_OBJ) $(OMP_OBJ)
	$(CC) $(CFLAGS) $^ -fopenmp -o ompRPS.out

pth : $(UTIL_OBJ) $(PTHREAD_OBJ)
	$(CC) $(CFLAGS) $^ -pthread -o pthreadRPS.out


remake : clean all

# We add .PHONY when a target doesn't actually create any output. In this case we just run a shell
# command, removing all object files, i.e files ending on .o
# the * syntax means [anything].o
.PHONY : clean
clean :
	rm -f src/obj/*.o && rm -f *.gch && rm -f *.bmp && rm -f *.mp4 && rm -f *.out

