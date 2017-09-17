#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>


int X = 2560;
int Y = 2048;

void calc(int rank, int comm_sz, int pixel[]) {
	for(int i=rank; i < X; i+=comm_sz) {
		for(int j=0; j < Y; j++) {
			/*if(rank == 1) {
				printf("%d %d %d\n", i, j, (i/comm_sz)*Y+j);
			}*/
			pixel[(i/comm_sz)*Y + j] = rank;
		}
	}
}

int main(){
	int comm_sz;
	int rank;
	int size = X*Y;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	int proc_rows = ((X - rank + comm_sz - 1) / comm_sz)*Y;
	int pixel[proc_rows];
	int result[size];

	calc(rank, comm_sz, pixel);
	if(1) {
		printf("FOO");
		MPI_Finalize();
		return 0;
	}

	if(rank == 0) {
		for(int i = 0; i < X; i+= comm_sz) {
			for(int j = 0; j < Y; j++) {
				result[i+(j*X)] = pixel[(i/comm_sz)*X+j];
			}
		}

		for(int r = 1; r < comm_sz; r++) {
			int r_size = ((X - r + comm_sz - 1) / comm_sz)*Y;
			int recieved[r_size];
			MPI_Recv(recieved, r_size, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Received from: %d\n", r);
			for(int i = 0; i < r_size; i++) {
				printf("%d ", recieved[i]);
			}
			printf("\n");
			for(int i = r; i < X; i+= comm_sz) {
				for(int j = 0; j < Y; j++) {
					printf("Result index: %d, Recieved index: %d, Value: %d\n", i+(j*X), (i/comm_sz)*Y+j, recieved[(i/comm_sz)*Y+j]);
					result[i+(j*X)] = recieved[(i/comm_sz)*Y+j];
				}
			}
			printf("\n");
		}
		for(int i = 0; i < size; i++) {
			printf("%d", result[i]);
		}

	} else {
		MPI_Send(pixel, proc_rows, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	
	return 0;
}
