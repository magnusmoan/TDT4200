#include "RPS.h"
#include <time.h>
#include <pthread.h>
#include <semaphore.h>

void swap_petris();
void* iterate(void* rank);

cell* petri_A;
cell* petri_B;

int rows_per_thread;

int thread_count;
int counter;
sem_t count_sem;
sem_t barrier_sem;

int main(int argc, char** argv){

    if(argc==1) {
        puts("Expected one command line argument, got 0");
        return 0;
    }

    printf("running %d iterations\n",ITERATIONS);
    
    long thread;
    pthread_t* thread_handles;
    thread_count = strtol(argv[1], NULL, 10);
    
    sem_init(&count_sem, 0, 1);
    sem_init(&barrier_sem, 0, 0);
    counter = 0;

    rows_per_thread = (IMG_X-2) / thread_count;

    thread_handles = malloc (thread_count*sizeof(pthread_t));

    srand(time(NULL));
    petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
    petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));

    int seed = rand();

    // Seed some CAs
    for(int ii = 0; ii < 100; ii++){
        int rx = rand() % (IMG_X - 1);
        int ry = rand() % (IMG_Y - 1);
        int rt = rand() % 4;

        petri_A[TRANS(rx,ry)].color = rt;
        petri_A[TRANS(rx,ry)].strength = 1;
    }

    for(thread = 0; thread < thread_count; thread++) {
        pthread_create(&thread_handles[thread], NULL, iterate, (void*) thread);
    }

    for(thread = 0; thread < thread_count; thread++) {
        pthread_join(thread_handles[thread], NULL);
    }
    free(thread_handles);

    sem_destroy(&count_sem);
    sem_destroy(&barrier_sem);
    char filename[50] = "RPS_pthread.bmp";
    make_bmp(petri_A, 0, filename);

    free(petri_A);
    free(petri_B);
    return 0;

}


void swap_petris(){
    cell* tmp1 = petri_A;
    petri_A = petri_B;
    petri_B = tmp1;
}

void* iterate(void* rank) {
    long my_rank = (long) rank;

    int start = my_rank * rows_per_thread + 1;
    int end = start + rows_per_thread; 
    
    for(int ii = 0; ii < 3000; ii++) {
        iterate_image(petri_A, petri_B, start, end);
        sem_wait(&count_sem);
        if(counter == thread_count-1) {
            counter = 0;
            swap_petris();
            sem_post(&count_sem);
            for(int jj = 0; jj < thread_count-1; jj++) {
                sem_post(&barrier_sem);
            }
        } else {
            counter++;
            sem_post(&count_sem);
            sem_wait(&barrier_sem);
        }
    }
   return NULL; 
}
