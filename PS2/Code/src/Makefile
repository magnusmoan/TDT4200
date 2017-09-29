CC_PAR=mpicc
CC=gcc
ODIR=obj
IDIR=../include
CFLAGS=-I$(IDIR) -std=c99

_DEPS = bitmap.h CA.h RPS.h RPS_MPI.h
DEPS = $(patsubst %,$(IDIR)/%,$(_UTIL_DEPS))

_SER_OBJ = RPS.o
_PAR_OBJ = RPS_MPI.o 
_UTIL_OBJ = bitmap.o CA.o
SER_OBJ = $(patsubst %,$(ODIR)/%,$(_SER_OBJ))
PAR_OBJ = $(patsubst %,$(ODIR)/%,$(_PAR_OBJ))
UTIL_OBJ = $(patsubst %,$(ODIR)/%,$(_UTIL_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC_PAR) -c -o $@ $< $(CFLAGS)

RPSparallel: $(PAR_OBJ) $(UTIL_OBJ)
	$(CC_PAR) -o $@ $^ $(CFLAGS) 

serial: $(SER_OBJ) $(UTIL_OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: all remake clean video
all: parallel serial

remake: clean all

clean:
	rm -f $(ODIR)/*.o && rm -f *.gch && rm -f *.bmp && rm -f *.mp4

video: clean all
	./serial && ffmpeg -framerate 60 -i data/CA-%000d.bmp -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
