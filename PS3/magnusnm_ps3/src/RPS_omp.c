#include "RPS.h"
#include <time.h>
#include <omp.h>

void swap_petris();
void iterate();

cell* petri_A;
cell* petri_B;

int rows_per_thread;

int main(int argc, char** argv){

    if(argc==1) {
        puts("Expected one command line argument, got 0");
        return 0;
    }

    int thread_count = strtol(argv[1], NULL, 10);

    printf("running %d iterations\n",ITERATIONS);

    srand(time(NULL));
    petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
    petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));

    rows_per_thread = (IMG_X-2) / thread_count;

    int seed = rand();

    // Seed some CAs
    for(int ii = 0; ii < 100; ii++){
        int rx = rand() % (IMG_X - 1);
        int ry = rand() % (IMG_Y - 1);
        int rt = rand() % 4;

        petri_A[TRANS(rx,ry)].color = rt;
        petri_A[TRANS(rx,ry)].strength = 1;
    }

    #pragma omp parallel num_threads(thread_count)
    { 
        iterate(); 
    }

    char filename[50] = "RPS_omp.bmp";
    make_bmp(petri_A, 0, filename);
}


void swap_petris(){
    cell* tmp1 = petri_A;
    petri_A = petri_B;
    petri_B = tmp1;
}

void iterate(void) {
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    int start = my_rank * rows_per_thread + 1;
    int end = start + rows_per_thread; 
    
    for(int ii = 0; ii < ITERATIONS; ii++) {
        iterate_image(petri_A, petri_B, start, end);
        # pragma omp barrier
        # pragma omp master
        { swap_petris(); }
        # pragma omp barrier
    }
}
