CC=gcc
ODIR=src/obj
IDIR=include
SRCDIR=src
CFLAGS=-I$(IDIR) -std=c99 -Wall

_DEPS = bitmap.h CA.h RPS.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = bitmap.o CA.o RPS.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_SRC = bitmap.c CA.c RPS.c
SRC= $(patsubst %,$(SRCDIR)/%,$(_SRC))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


serial: $(OBJ)
	$(CC) -std=c99 -g -O3 $^ $(CFLAGS) -v -o $@

omp: $(OBJ)
	$(CC) -fopenmp -o ompRPS.out $^ $(CFLAGS)

pth: $(OBJ)
	$(CC) -pthread -o pthreadRPS.out $^ $(CFLAGS)

all: serial omp pth

remake: clean all

# We add .PHONY when a target doesn't actually create any output. In this case we just run a shell
# command, removing all object files, i.e files ending on .o
# the * syntax means [anything].o
.PHONY : clean
clean :
	rm -f src/obj/*.o && rm -f src/*.gch && rm -f src/data/*.bmp && rm -f src/*.mp4

# Finally, the test target. Builds the 'all' target, then runs the test script on the output
.PHONY : video
video : clean all
	./myProgram && ffmpeg -framerate 60 -i data/CA-%000d.bmp -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

.PHONY : run
run : clean omp
	./meme.out
