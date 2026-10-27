#ifndef SANDPILE_CONTEXT_H
#define SANDPILE_CONTEXT_H

#include <mpi.h>
#include <stdbool.h>

//global variables to hold color code and height constraint
#define MAXH 3
#define NRGB 3

//---Data structure to hold all local variables

typedef struct {
    //Global and comm details
    int L;               // Global side length of the grid
    int size;            // Total number of processes (MPI_COMM_WORLD)
    int rank;            // Rank in MPI_COMM_WORLD
    int new_rank;        // Rank in comm_2d (Cartesian)
    MPI_Comm comm;      // Original communicator
    MPI_Comm comm_2d;    // 2D cartesian communicator

    //Domain decomposition
    int PX, PY;          // Processor counts in X and Y dimensions
    int LX, LY;          // Local domain size (excluding halos, includes remainder)
    int extra_x, extra_y; // Remainders for non-uniform decomposition
    int coords[2];       // Processor coordinates (x, y)
    int neighbours_ranks[4]; // [Top, Bottom, Left, Right]

    // Simulation Arrays (LX+2 x LY+2 including halos)
    int **old_height;
    int **new_height;

    // Halo Exchange Buffers
    int *left_col_send, *left_col_recv;
    int *right_col_send, *right_col_recv;
    
    // Custom MPI Datatypes
    MPI_Datatype ROW_DATA_TYPE;         // For X-direction halos (rows)
    MPI_Datatype CONT_DATA_FOR_COLUMN;  // For Y-direction halos (columns)
    
    // Initial State & Loop Control
    int offset_x, offset_y;
    int start_rank;
    int max_iterations;
    int num_loops;
    int start_height;
    
    // Timing variables
    double tstart, tloop; 
} SandpileContext;


//functions to write output
void sandwritedynamic(char *cellfile, int **cell, int l);

//Initializes MPI, sets up the 2D Cartesian topology, and calculates domain sizes.
void sandpile_initialization(SandpileContext *sand_ctx, int argc, char *argv[]);

/*
Executes one simulation step halo exchange, computation, and reduction.
returns : total change
*/
int run_simulation_step(SandpileContext *sand_ctx);

//Gather result at rank = 0 and write output and timing data
void finalize_and_cleanup(SandpileContext *sand_ctx, double tend);

#endif