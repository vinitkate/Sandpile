#include <stdio.h>
#include <stdlib.h>
// #include "sandpile.h"
#include "sandpile_context.h"

#define PRINT_FREQ 5000
#define ROOT_RANK 0

//Helper function to update the old_height matrix
static int update_heights(int **oldh, int **newh, int x_start, int x_end, int y_start, int y_end){
	int i, j, h;
	int change = 0;
	for(i = x_start; i <= x_end; i++){
		for(j = y_start; j <= y_end; j++){
			h = oldh[i][j];
			if(h > MAXH) h -= 4;
			if(oldh[i-1][j] > MAXH) h++;
			if(oldh[i+1][j] > MAXH) h++;
			if(oldh[i][j-1] > MAXH) h++;
			if(oldh[i][j+1] > MAXH) h++;
			if(h != oldh[i][j]) change++;
			newh[i][j] = h;
		}
	}
	return change;
}

/*
exchange_halos(SandpileContext *, MPI_Requests[]): exchaning halos

---------------------------------------------------------------------------------------
run_simulation_step(SandpileContext *): Public API function: runs a single step of simulation
---------------------------------------------------------------------------------------
*/

//Helper function to perform non-blocking halo exchange
static void exchange_halos(SandpileContext *sand_ctx, MPI_Request requests[]){
    int i;
    // ### 1. Halo Exchange (Communication) ---

	// --- X-Direction Exchange (Rows, since X is the first index in array) ---
	// Top row to neighbor[0] (x-1)
	MPI_Isend(&sand_ctx->old_height[1][1], 1, sand_ctx->ROW_DATA_TYPE, sand_ctx->neighbours_ranks[0], 0, sand_ctx->comm_2d, &requests[0]);
	// Receive halo from neighbor[0] into old_height[0][1] (the top halo row)
	MPI_Irecv(&sand_ctx->old_height[0][1], 1, sand_ctx->ROW_DATA_TYPE, sand_ctx->neighbours_ranks[0], 1, sand_ctx->comm_2d, &requests[1]);
		
	// Bottom row for neighbor[1] (x+1)
	MPI_Isend(&sand_ctx->old_height[sand_ctx->LX][1], 1, sand_ctx->ROW_DATA_TYPE, sand_ctx->neighbours_ranks[1], 1, sand_ctx->comm_2d, &requests[2]);
	// Receive halo from neighbor[1] into old_height[LX+1][1] (the bottom halo row)
	MPI_Irecv(&sand_ctx->old_height[sand_ctx->LX+1][1], 1, sand_ctx->ROW_DATA_TYPE, sand_ctx->neighbours_ranks[1], 0, sand_ctx->comm_2d, &requests[3]);
		
    //Point to point comms to exchange columns 
	// --- Y-Direction Exchange :(columns, periodic)
    // Pack the internal columns into send buffers 
	for(i = 1; i <= sand_ctx->LX; i++){
		sand_ctx->left_col_send[i-1] = sand_ctx->old_height[i][1];
		sand_ctx->right_col_send[i-1] = sand_ctx->old_height[i][sand_ctx->LY];
	}
		
	//Send Left column to neighbor[2] (y-1)
	MPI_Isend(sand_ctx->left_col_send, 1, sand_ctx->CONT_DATA_FOR_COLUMN, sand_ctx->neighbours_ranks[2], 2, sand_ctx->comm_2d, &requests[4]);
	// Receive halo from neighbor[2] into left_col_recv buffer
	MPI_Irecv(sand_ctx->left_col_recv, 1, sand_ctx->CONT_DATA_FOR_COLUMN, sand_ctx->neighbours_ranks[2], 3, sand_ctx->comm_2d, &requests[5]);
		
	// send Right column to neighbor[3] (y+1)
	MPI_Isend(sand_ctx->right_col_send, 1, sand_ctx->CONT_DATA_FOR_COLUMN, sand_ctx->neighbours_ranks[3], 3, sand_ctx->comm_2d, &requests[6]);
	// Receive halo from neighbor[3] into right_col_recv buffer
	MPI_Irecv(sand_ctx->right_col_recv, 1, sand_ctx->CONT_DATA_FOR_COLUMN, sand_ctx->neighbours_ranks[3], 2, sand_ctx->comm_2d, &requests[7]);
		
}

//Function to run a simulation step---------------------------------------------------
int run_simulation_step(SandpileContext *sand_ctx){
    int i, j;
    int local_change = 0;
    int total_change = 0;

    MPI_Request requests[8];
    double tloops;

    MPI_Barrier(sand_ctx->comm_2d);
	tloops = MPI_Wtime(); // Start time for the current iteration loop

    // --- 1. Halo Exchange Start ---
    exchange_halos(sand_ctx, requests);

    // --- 2. Computation/Communication Overlap ---
    // Compute heights for the core area (2 to LX-1, 2 to LY-1) 
    local_change += update_heights(sand_ctx->old_height, sand_ctx->new_height, 2, sand_ctx->LX - 1, 2, sand_ctx->LY - 1);

    // --- 3. Synchronize and update Halos ---
    MPI_Waitall(8, requests, MPI_STATUSES_IGNORE);

    // Update received column buffers into the actual halo cells (column 0 and LY+1)
	for(i = 1; i <= sand_ctx->LX; i++){
		sand_ctx->old_height[i][0] = sand_ctx->left_col_recv[i-1];
		sand_ctx->old_height[i][sand_ctx->LY+1] = sand_ctx->right_col_recv[i-1];
	}

    // ### 4. Finish Remaining Boundary Computation ---
    // Compute height changes for the boundary regions that rely on the received halo data.

	// Boundary columns (Col 1 and Col LY)
	local_change += update_heights(sand_ctx->old_height, sand_ctx->new_height, 1, sand_ctx->LX, 1, 1); // Leftmost inner column (col 1)
	local_change += update_heights(sand_ctx->old_height, sand_ctx->new_height, 1, sand_ctx->LX, sand_ctx->LY, sand_ctx->LY); // Rightmost inner column (col LY)

	// Boundary rows (Row 1 and Row LX, excluding corners already handled by above calls)
	local_change += update_heights(sand_ctx->old_height, sand_ctx->new_height, 1, 1, 2, sand_ctx->LY-1); // Topmost inner row (row 1, inner cells)
	local_change += update_heights(sand_ctx->old_height, sand_ctx->new_height, sand_ctx->LX, sand_ctx->LX, 2, sand_ctx->LY-1); // Bottommost inner row (row LX, inner cells)

    // ### 5. Global Termination Check ---
    // Sum the local changes across all processors using MPI_Allreduce
	MPI_Allreduce(&local_change, &total_change, 1, MPI_INT, MPI_SUM, sand_ctx->comm_2d);

    // ### 6. Update State ---
    // Copy new_height back to old_height for the next iteration (only inner cells)
	for(i = 1; i <= sand_ctx->LX; i++){
		for(j = 1; j <= sand_ctx->LY; j++){
			sand_ctx->old_height[i][j] = sand_ctx->new_height[i][j];
		}
	}

    sand_ctx->num_loops++;

    MPI_Barrier(sand_ctx->comm_2d);
	sand_ctx->tloop += MPI_Wtime() - tloops; // Accumulate time spent in the a simulation loop
    // Output and Termination Condition ---
	if(sand_ctx->rank == 0 ){
		if(sand_ctx->num_loops % PRINT_FREQ == 0){
    		printf("Total change %d -> step %d\n", total_change, sand_ctx->num_loops);
		}
			
		if(total_change == 0){
			printf("Total change %d -> step %d\n", total_change, sand_ctx->num_loops);
			printf("Done !!!\n");
		}
	}

	return total_change;
}
