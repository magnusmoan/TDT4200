all: nbody nbody_cuda

nbody: nbody.c
	gcc -std=c99 -O0 -g -o nbody nbody.c -lm

nbody_cuda: nbody_cuda.cu
	nvcc -O0 -arch=sm_61 -o nbody_cuda nbody_cuda.cu -lm -lcuda -lcudart
	
clean:
	-rm -f nbody nbody_cuda
