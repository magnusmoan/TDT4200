IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -std=c99

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = bitmap.h julia.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = bitmap.o julia.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

julia: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
