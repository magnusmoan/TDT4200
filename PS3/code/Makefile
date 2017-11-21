CFLAGS+= -std=c99 -g -O3

# by uncommenting this line the preprocessor will see #ifdef DEBUG as true
# CFLAGS+= -DDEBUG
util_objects = CA.o bitmap.o 
serial_objects = RPS_serial.o
parallel_pth = RPS_pthread.o
omp_objects = RPS_omp.o

SHELL := /bin/bash

serial : $(util_objects) $(serial_objects)
	gcc -std=c99 -g -O3 $(util_objects) $(serial_objects) -o serialRPS.out

omp : $(util_objects) $(omp_objects)
	gcc -std=c99 -g  -fopenmp $(util_objects) $(omp_objects) -o ompRPS.out

pth : $(util_objects) $(parallel_pth)
	gcc -std=c99 -g -O3 -pthread $(util_objects) $(parallel_pth) -o pthreadRPS.out


remake : clean all

# We add .PHONY when a target doesn't actually create any output. In this case we just run a shell
# command, removing all object files, i.e files ending on .o
# the * syntax means [anything].o
.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f data/*.bmp && rm -f *.mp4 && rm -f *.out

# Finally, the test target. Builds the 'all' target, then runs the test script on the output
.PHONY : video
video : clean all
	./myProgram && ffmpeg -framerate 60 -i data/CA-%000d.bmp -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

.PHONY : run
run : clean omp
	./meme.out
