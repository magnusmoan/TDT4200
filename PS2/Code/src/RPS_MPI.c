#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

#include "RPS_MPI.h"

void initialize();
void exchange_borders(cell* petri);
void iterate_CA(cell* old, cell* new);
void gather_petri(cell* local_petri);
void create_types();

// 2D to 1D Helpers
int get_1d_index(int x, int y);
int get_north_row_index(int offset);
int get_south_row_index(int offset);
int get_west_column_index(int offset);
int get_east_column_index(int offset);

// CA methods adjusted for a 1D world
bool mpi_beats(cell me, cell other);
void flat_iterate_image(cell* old_petri, cell* new_petri);
cell flat_pick_neighbor(int x, int y, cell* petri);
cell flat_next_cell(int x, int y, cell* petri);

int rank;
int size;

// I denote mpi process specific values with hungarian notation, adding a p

// The dimensions of the processor grid. Same for every process
int p_x_dims;
int p_y_dims;

// The location of a process in the process grid. Unique for every process
int p_my_x_dim;
int p_my_y_dim;

int p_north, p_south, p_east, p_west;

// The dimensions for the process local petri
int p_local_petri_x_dim;
int p_local_petri_y_dim;

// All local petris have an extra border around it representing the borders to its neighbors
int p_local_petri_x_dim_with_border;
int p_local_petri_y_dim_with_border;

MPI_Comm cart_comm;

// some datatypes, useful for sending data with somewhat less primitive semantics
MPI_Datatype border_row_t;  // TODO: Implement this
MPI_Datatype border_col_t;  // TODO: Implement this
MPI_Datatype local_petri_t; // Already implemented
MPI_Datatype mpi_cell_t;    // Already implemented

// Each process is responsible for one part of the petri dish.
// Since we can't update the petri-dish in place each process actually
// gets two petri-dishes which they update in a lockstep fashion.
// dish A is updated by writing to dish B, then next step dish B updates dish A.
// (or you can just swap them inbetween iterations)
cell* local_petri_A;
cell* local_petri_B;


// Storing the entire petri
cell** petri;


// For sending and receiving in the border exchange
cell* send;
cell* recv;

int main(int argc, char** argv){

	// Ask MPI what size (number of processors) and rank (which process we are)
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Timing the program
	double local_start, local_finish, local_elapsed, elapsed;
	MPI_Barrier(MPI_COMM_WORLD);
	local_start = MPI_Wtime();

	// Making sure each processor have a different seed so we dont work with identical local petris
	srand(rank * 1000);

	////////////////////////////////
	// Create cartesian communicator
	int dims[2];
	dims[0] = p_x_dims;
	dims[1] = p_y_dims;

	int periods[2]; // we set these to 0 because we are not interested in wrap-around
	periods[0] = 0;
	periods[1] = 0;

	int coords[2];
	coords[0] = p_my_x_dim;
	coords[1] = p_my_y_dim;

	MPI_Dims_create(size, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
	MPI_Cart_coords(cart_comm, rank, 2, coords);

	MPI_Cart_shift(cart_comm, 0, 1, &p_north, &p_south);
	MPI_Cart_shift(cart_comm, 1, 1, &p_west, &p_east);

	p_x_dims = dims[0];
	p_y_dims = dims[1];

	p_my_x_dim = coords[0];
	p_my_y_dim = coords[1];
	////////////////////////////////
	////////////////////////////////


	initialize();
	create_types();

	int step;
	int it = ITERATIONS;
	for(step = 0; step < it; step++) {
		if(step % 2 == 0) {
			iterate_CA(local_petri_A, local_petri_B);;
		  	exchange_borders(local_petri_A);
	  	} else {
		  	iterate_CA(local_petri_B, local_petri_A);
		  	exchange_borders(local_petri_B);
	  	}
	}

	gather_petri(local_petri_A);
	if(rank == 0) {
	  	make_bmp(petri, step);
	}

	local_finish = MPI_Wtime();
	local_elapsed = local_finish - local_start;

	MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank == 0) {
		printf("Elapsed time: %lf\n seconds", elapsed);
	}

	MPI_Finalize();


	// free_stuff()
	if(rank==0){
	  	for(int ii = 0; ii < IMG_X; ii++){
		  	free(petri[ii]);
	  	}
	  	free(petri);
	}

	free(local_petri_A);
	free(local_petri_B);
	free(send);
	free(recv);

  exit(0);
}


void create_types(){

	// cell type
	const int    nitems=2;
	int          blocklengths[2] = {1,1};
	MPI_Datatype types[2] = {MPI_INT, MPI_INT};
	MPI_Aint     offsets[2];

	offsets[0] = offsetof(cell, color);
	offsets[1] = offsetof(cell, strength);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
	MPI_Type_commit(&mpi_cell_t);

	// A message for a local petri-dish
	MPI_Type_contiguous(p_local_petri_x_dim_with_border * p_local_petri_y_dim_with_border,
		      mpi_cell_t,
		      &local_petri_t);
	MPI_Type_commit(&local_petri_t);

	// A border row
	MPI_Type_contiguous(p_local_petri_x_dim,
		      mpi_cell_t,
		      &border_row_t);
	MPI_Type_commit(&border_row_t);

	// A border column
	MPI_Type_contiguous(p_local_petri_y_dim,
		      mpi_cell_t,
		      &border_col_t);
	MPI_Type_commit(&border_col_t);

}


void initialize(){
	p_local_petri_x_dim = IMG_X/p_x_dims;
	p_local_petri_y_dim = IMG_Y/p_y_dims;
	p_local_petri_x_dim_with_border = p_local_petri_x_dim + 2;
	p_local_petri_y_dim_with_border = p_local_petri_y_dim + 2;

	// Allocating space for the petri dishes. Adding extra rows and columns on the edges to store the borders
	local_petri_A = (cell*)malloc(p_local_petri_x_dim_with_border*p_local_petri_y_dim_with_border*sizeof(mpi_cell_t));
	local_petri_B = (cell*)malloc(p_local_petri_x_dim_with_border*p_local_petri_y_dim_with_border*sizeof(mpi_cell_t));

	// Allocating space for send and recv used in the border exchange
	int cells_to_send = 2*(p_local_petri_y_dim) + 2*(p_local_petri_x_dim);	
	send = malloc(cells_to_send*sizeof(cell));
	recv = malloc(cells_to_send*sizeof(cell));

	// Setting all cells in the A buffer to be white
	for(int i = 0; i < p_local_petri_x_dim_with_border*p_local_petri_y_dim_with_border; i++) {
		local_petri_A[i].color = 0;
		local_petri_A[i].strength = 0;
	}

	// Assigning random values to 100 random cells
	for(int ii = 0; ii < 100; ii++) {
	  int rx = rand() % p_local_petri_x_dim;
	  int ry = rand() % p_local_petri_y_dim;
	  int rt = rand() % 4;

	  local_petri_A[get_1d_index(rx,ry)].color = rt;
	  local_petri_A[get_1d_index(rx,ry)].strength = 1;
	}

}


void exchange_borders(cell* local_petri){
	int sizes[4] = {1, 1, 1, 1};
	long send_dis[4] = {0, (p_local_petri_x_dim)*sizeof(cell), 2*(p_local_petri_x_dim)*sizeof(cell), 
		(2*(p_local_petri_x_dim) + p_local_petri_y_dim)*sizeof(cell)};

	MPI_Datatype send_types[4] = {border_row_t, border_row_t, border_col_t, border_col_t};
		
	long recv_dis[4] = {0, (p_local_petri_y_dim)*sizeof(cell), 2*(p_local_petri_y_dim)*sizeof(cell), 
		(2*(p_local_petri_y_dim) + p_local_petri_x_dim)*sizeof(cell)};

	MPI_Datatype recv_types[4] = {border_col_t, border_col_t, border_row_t, border_row_t};
	

	// Adding the borders and rows that are going to be sent to the send pointer
	int send_index = 0;
	for(int ii = 1; ii < p_local_petri_x_dim + 1; ii++) {
		send[send_index] = local_petri[get_north_row_index(ii)];
		send[send_index + p_local_petri_x_dim] = local_petri[get_south_row_index(ii)];
		send_index++;
	}
	send_index += p_local_petri_x_dim;

	for(int ii = 1; ii < p_local_petri_y_dim + 1; ii++) {
		send[send_index] = local_petri[get_west_column_index(ii)];
		send[send_index + p_local_petri_y_dim] = local_petri[get_east_column_index(ii)];
		send_index++;
	}
	////////////////////////////////
	////////////////////////////////
	
	// Sending and receiving
	MPI_Neighbor_alltoallw(send, sizes, send_dis, send_types, recv, sizes, recv_dis, recv_types, cart_comm);

	// Handling the received data
	if(p_north >= 0) {
		for(int ii = 0; ii < p_local_petri_y_dim; ii++) {
			local_petri[ii+1] = recv[ii];
		}
	}
	if(p_south >= 0) {
		int bottom_row_index = p_local_petri_x_dim_with_border*(p_local_petri_y_dim_with_border-1)+1;
		for(int ii = 0; ii < p_local_petri_x_dim; ii++) {
			local_petri[bottom_row_index + ii] = recv[(p_local_petri_y_dim) + ii];
		}
	}
	if(p_west >= 0) {
		for(int ii = 0; ii < p_local_petri_y_dim; ii++) {
			local_petri[(ii+1)*p_local_petri_x_dim_with_border] = recv[2*p_local_petri_y_dim + ii];
		}
	}
	if(p_east >= 0) {
		int bottom_row_index = p_local_petri_x_dim_with_border*(p_local_petri_y_dim_with_border-1)+1;
		for(int ii = 0; ii < p_local_petri_y_dim; ii++) {
			local_petri[(ii+2)*p_local_petri_x_dim_with_border - 1] = 
				recv[2*p_local_petri_y_dim + p_local_petri_x_dim + ii];
		}
	}
}

void iterate_CA(cell* old_petri, cell* new_petri){
	for(int ii = 0; ii < p_local_petri_y_dim; ii++){
		for(int jj = 0; jj < p_local_petri_x_dim; jj++){
			new_petri[get_1d_index(jj,ii)]  = flat_next_cell(jj, ii, old_petri);
    		}
 	}
}

void gather_petri(cell* local_petri){
	if(rank == 0) {
		// Initializing the full petri/img
		petri = (cell**)malloc(IMG_X*sizeof(cell*));

		for(int ii = 0; ii < IMG_X; ii++) {
			cell* col = (cell*)malloc(IMG_Y*sizeof(cell));
			petri[ii] = col;
		}

		// Inserting the local petri of processor with rank 0 into the full petri
		for(int ii = 0; ii < p_local_petri_y_dim; ii++) {
			for(int jj = 0; jj < p_local_petri_x_dim; jj++) {
				petri[jj][ii] = local_petri[get_1d_index(ii,jj)];
			}
		}

		// Inserting petri from all processor except for rank 0 into the full petri
		for(int r = 0; r < size-1; r++) {
			cell* petri_part = (cell*)malloc(p_local_petri_x_dim_with_border*
					p_local_petri_y_dim_with_border*sizeof(mpi_cell_t));
			MPI_Status status;
			MPI_Recv(petri_part, 1, local_petri_t, MPI_ANY_SOURCE, 0, cart_comm, &status);

			int curr_coords[2];
			MPI_Cart_coords(cart_comm, status.MPI_SOURCE, 2, curr_coords);
			int curr_x = curr_coords[0]*p_local_petri_x_dim;
			int curr_y = curr_coords[1]*p_local_petri_y_dim;
			for(int ii = 0; ii < p_local_petri_y_dim; ii++) {
				for(int jj = 0; jj < p_local_petri_x_dim; jj++) {
					petri[curr_x + jj][curr_y + ii] = petri_part[get_1d_index(ii,jj)];
				}
			}
			free(petri_part);
		}

	} else {
		MPI_Send(local_petri, 1, local_petri_t, 0, 0, cart_comm);
	}
}


// Helpers for calculating indexes in a 1d world
int get_1d_index(int x, int y) {
	return p_local_petri_x_dim_with_border * (y + 1) +  (x + 1);
}

int get_north_row_index(int offset) {
	return p_local_petri_x_dim_with_border + offset;
}

int get_south_row_index(int offset) {
	return p_local_petri_x_dim_with_border*p_local_petri_y_dim + offset;
}

int get_west_column_index(int offset) {
	return p_local_petri_x_dim_with_border*offset + 1;
}

int get_east_column_index(int offset) {
	return p_local_petri_x_dim_with_border*(offset+1) - 2;
}


// Methods used in the CA. Similar to the methods found in CA.c, but adjusted to work on 1 dimension instead of 2.
cell flat_pick_neighbor(int x, int y, cell* petri) {
	int chosen = rand() % 8;

	if(chosen >= 4) { chosen++; }
	int c_x = chosen % 3;
	int c_y = chosen / 3;
	int index  = get_1d_index(x,y) + (c_x-1) + (c_y-1)*p_local_petri_x_dim_with_border;

	return petri[index];
}

cell flat_next_cell(int ii, int jj, cell* petri) {
	cell neighbor_cell = flat_pick_neighbor(ii, jj, petri);
	cell my_cell = petri[get_1d_index(ii,jj)];
	if(neighbor_cell.color == WHITE){
		return my_cell;
	}
	cell next_cell = my_cell;

	if(my_cell.color == WHITE){
		next_cell.strength = 1;
		next_cell.color = neighbor_cell.color;
		return next_cell;
	} else {
		if(mpi_beats(my_cell, neighbor_cell)){
			next_cell.strength++;
		} else {
			next_cell.strength--;
		}
	}

	if(next_cell.strength == 0){
		next_cell.color = neighbor_cell.color;
		next_cell.strength = 1;
	}

	if(next_cell.strength > 4){
		next_cell.strength = 4;
	}

	return next_cell;
}

bool mpi_beats(cell me, cell other){
  return
    (((me.color == SCISSOR) && (other.color == PAPER)) ||
     ((me.color == PAPER) && (other.color == ROCK))    ||
     ((me.color == ROCK) && (other.color == SCISSOR))  ||
     (me.color == other.color));
}
