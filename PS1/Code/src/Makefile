IDIR =../include
CC=mpicc

ifeq ($(s),)
	s=1
endif

CFLAGS=-I$(IDIR) -std=c99 -D S=$(s)

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = bitmap.h julia_mpi.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = bitmap.o julia_mpi.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

julia_mpi: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
