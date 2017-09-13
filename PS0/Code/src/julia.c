#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <julia.h>

double x_start=-2.01;
double x_end=1;
double yupper;
double ylower;

double ycenter=1e-6;
double step;

int pixel[XSIZE*YSIZE];


// I suggest you implement these, however you can do fine without them if you'd rather operate
// on your complex number directly.
complex_t square_complex(complex_t c){
	complex_t squared;
	squared.real = c.real*c.real - c.imag*c.imag;
	squared.imag = 2*c.real*c.imag;

	return squared;
}

complex_t add_complex(complex_t a, complex_t b){
	complex_t sum;
	sum.real = a.real + b.real;
	sum.imag = a.imag + b.imag;

	return sum;
}

complex_t add_real(complex_t a, int b){
	complex_t sum;
	sum.real = a.real + ((float) b);
	sum.imag = a.imag;

	return sum;
}

complex_t add_imag(complex_t a, int b) {
	complex_t sum;
	sum.real = a.real;
	sum.imag = a.imag + ((float) b);

	return sum;
}



// add julia_c input arg here?
void calculate(complex_t julia_C) {
	for(int i=0;i<XSIZE;i++) {
		for(int j=0;j<YSIZE;j++) {
			/* Calculate the number of iterations until divergence for each pixel.
			   If divergence never happens, return MAXITER */
			complex_t c;
      			complex_t z;
      			complex_t temp;
			int iter=0;

      			// find our starting complex number c
			c.real = (x_start + step*i);
			c.imag = (ylower + step*j);

      			// our starting z is c
			z = c;

      			// iterate until we escape
			while(z.real*z.real + z.imag*z.imag < 4) {
				// Each pixel in a julia set is calculated using z_n = (z_n-1)² + C
				// C is provided as user input, so we need to square z and add C until we
				// escape, or until we've reached MAXITER

				z = add_complex(square_complex(z), julia_C);
				if(++iter==MAXITER) break;
			}
			pixel[PIXEL(i,j)]=iter;
		}
	}
}


int main(int argc,char **argv) {
	if(argc==1) {
		puts("Usage: JULIA\n");
		puts("Input real and imaginary part. ex: ./julia 0.0 -0.8");
		return 0;
	}

	/* Calculate the range in the y-axis such that we preserve the
	   aspect ratio */
	step=(x_end-x_start)/XSIZE;
	yupper=ycenter+(step*YSIZE)/2;
	ylower=ycenter-(step*YSIZE)/2;

	// Unlike the mandelbrot set where C is the coordinate being iterated, the
	// julia C is the same for all points and can be chosed arbitrarily
	complex_t julia_C;

	// Get the command line args
	julia_C.real = strtod(argv[1], NULL);
	julia_C.imag = strtod(argv[2], NULL);

	calculate(julia_C);

	/* create nice image from iteration counts. take care to create it upside
	down (bmp format) */
	unsigned char *buffer=calloc(XSIZE*YSIZE*3,1);
	for(int i=0;i<XSIZE;i++) {
		for(int j=0;j<YSIZE;j++) {
			int p=((YSIZE-j-1)*XSIZE+i)*3;
			fancycolour(buffer+p,pixel[PIXEL(i,j)]);
		}
	}

	/* write image to disk */
    char file_name[35];
    sprintf(file_name, "../bmp/julia_%.3f_%.3fi.bmp", julia_C.real, julia_C.imag);
	savebmp(file_name,buffer,XSIZE,YSIZE);

    printf("Suksess. Se filen %s for resultat\n", file_name);
	return 0;
}
