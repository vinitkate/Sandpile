#include <stdio.h>
#include <stdlib.h>

#include "sandpile_context.h"
#include "arraymalloc.h"  

#define MAX_LOOPS 500000
/*
    Functions to parse the input, set up cartesian topology,
    set up the data types for contiguous ROW and COLUMN
    parse_bcast_input(SandpileContext *, int, char **): checks input and 
                    broadcasts parameters to all processors
    setup_topology_and_domain(SandpileContext *): 2D cartesian topology,
                    initialize the correct LX and LY
    setup_datatypes(SandpileContext *): datatypes for exchanging halos

    ------------------------------------------------------------------------------
    sandpile_initialization(SandpileContext *, int, char **): API call function 
    ------------------------------------------------------------------------------
*/


// Helper function to handle input parsing and broadcasting=========================================================
static void parse_bcast_input(SandpileContext *sand_ctx, int argc, char *argv[]){
    int L = 512, offset_x = 256, offset_y = 256;

    if(sand_ctx->rank == 0){
        if(argc != 2 && argc != 4){
            fprintf(stderr, "1. Usage: <sandpile> <offset>\n");
            fprintf(stderr, "2. Usage: <sandpile> <sandpile side_len> <offset_x> <offset_y>\n");
            // Standard error handling for incorrect usage
            MPI_Abort(sand_ctx->comm, EXIT_FAILURE);
        }
        if(argc == 2){
            offset_x = atoi(argv[1]);
            offset_y = offset_x;
        }
        else if(argc == 4){
            L = atoi(argv[1]);
            offset_x = atoi(argv[2]);
            offset_y = atoi(argv[3]);
        }

        if(offset_x < 0 || offset_x >= L || offset_y < 0 || offset_y >= L){
            fprintf(stderr, "0 <= <offset_x> < L and 0 <= <offset_y> < L\n");
            MPI_Abort(sand_ctx->comm, EXIT_FAILURE);
        }
    }

    MPI_Comm_size(sand_ctx->comm, &sand_ctx->size);

    // --- Broadcast Configuration ---
	// Broadcast the grid side length and the offset coordinates to all ranks.
	int temp[3] = {L, offset_x, offset_y};
	MPI_Bcast(temp, 3, MPI_INT, 0, sand_ctx->comm);

	sand_ctx->L = temp[0];
	sand_ctx->offset_x = temp[1];
	sand_ctx->offset_y = temp[2];
}

// Helper function to set up the 2D Cartesian topology and domain sizes
static void setup_topology_and_domain(SandpileContext *sand_ctx){
    // --- 2D Cartesian Topology Creation ---
    // MPI determines the optimal 2D dimensions (P = PX * PY)
	int dims[2] = {0, 0};
	MPI_Dims_create(sand_ctx->size, 2, dims);

    // Non-periodic in X-dimension assuming -> (0), Periodic in the Y-dimension assuming -> (1)
	int periods[2] = {false, true};
	int reorder = true;
	
	//Creating new comm using the dims periodic in second dim
	MPI_Cart_create(sand_ctx->comm, 2, dims, periods, reorder, &sand_ctx->comm_2d);

    MPI_Comm_rank(sand_ctx->comm_2d, &sand_ctx->new_rank);
	
	sand_ctx->PX = dims[0]; // Number of processors in X-dimension (Rows)
	sand_ctx->PY = dims[1]; // Number of processors in Y-dimension (Columns)

    // Initial calculation for local sizes (LX, LY) based on even division
	//Update later based on LX, PX
	sand_ctx->LX = sand_ctx->L / sand_ctx->PX;
	sand_ctx->LY = sand_ctx->L / sand_ctx->PY;
    
    // Calculate remainders for non-uniform decomposition (for LX+1 or LY+1 domains)
    sand_ctx->extra_x = sand_ctx->L % sand_ctx->PX;
    sand_ctx->extra_y = sand_ctx->L % sand_ctx->PY;

    MPI_Cart_coords(sand_ctx->comm_2d, sand_ctx->new_rank, 2, sand_ctx->coords);

    //// Update the local size (LX, LY) for the current processor based on remainder
	if(sand_ctx->coords[0] < sand_ctx->extra_x)
		(sand_ctx->LX)++;
	if(sand_ctx->coords[1] < sand_ctx->extra_y)
		(sand_ctx->LY)++;

    // ### Neighbor Determination 
	// X-dir (index 0) shift: neighbours_ranks[0] = Top (x-1), [1] = Bottom (x+1)
	MPI_Cart_shift(sand_ctx->comm_2d, 0, 1, &sand_ctx->neighbours_ranks[0], &sand_ctx->neighbours_ranks[1]);
	// Y-dir (index 1) shift: neighbours_ranks[2] = Left (y-1), [3] = Right (y+1) (Periodic)
	MPI_Cart_shift(sand_ctx->comm_2d, 1, 1, &sand_ctx->neighbours_ranks[2], &sand_ctx->neighbours_ranks[3]);
}

//Helper function for setting up datatypes =========================================================================
static void setup_datatypes(SandpileContext *sand_ctx){
    // CONT_DATA_FOR_COLUMN: Contiguous block of LX integers (used for Y-halo (column) exchange buffers)
    // The buffer size is LX (height of the column). 
	MPI_Type_contiguous(sand_ctx->LX, MPI_INT, &sand_ctx->CONT_DATA_FOR_COLUMN);
	MPI_Type_commit(&sand_ctx->CONT_DATA_FOR_COLUMN);

	// ROW_DATA_TYPE: Contiguous block of LY integers (used for Y-halo exchange)
    // The buffer size is LY (width of the row).
    MPI_Type_contiguous(sand_ctx->LY, MPI_INT, &sand_ctx->ROW_DATA_TYPE);
    MPI_Type_commit(&sand_ctx->ROW_DATA_TYPE);

}

//Helper function for array allocation and initialization ============================================================
static void setup_arrays_and_sandpile(SandpileContext *sand_ctx){
    int i, j;
    int start_x, start_y;
    int lx, ly;

    lx = sand_ctx->L / sand_ctx->PX;
    ly = sand_ctx->L / sand_ctx->PY;

    //Initializing arrays to store local heights and updates
	sand_ctx->new_height = (int **)arraymalloc2d(sand_ctx->LX+2, sand_ctx->LY+2, sizeof(int));
	sand_ctx->old_height = (int **)arraymalloc2d(sand_ctx->LX+2, sand_ctx->LY+2, sizeof(int));

    for(i = 0; i < sand_ctx->LX+2; i++){
		for(j = 0; j < sand_ctx->LY+2; j++){
			sand_ctx->new_height[i][j] = 0;
			sand_ctx->old_height[i][j] = 0;
		}
	}

    // Allocation for the temporary column buffers (used to pack/unpack X-direction halos)
	sand_ctx->left_col_send = (int *)malloc(sizeof(int) * sand_ctx->LX);
	sand_ctx->left_col_recv = (int *)calloc(sand_ctx->LX, sizeof(int));
	sand_ctx->right_col_send = (int *)malloc(sizeof(int) * sand_ctx->LX);
	sand_ctx->right_col_recv = (int *)calloc(sand_ctx->LX, sizeof(int));

    // --- Determine Local Domain Size and Initial Grain Processor 
    // Calculate the coordinate (start_x, start_y) of the processor containing the initial offset.
    // This logic handles the non-uniform distribution to find the correct processor index.
	if(sand_ctx->offset_x < sand_ctx->extra_x*(lx+1)){
		start_x = sand_ctx->offset_x/(lx+1);
	}else{
		start_x = sand_ctx->extra_x + (sand_ctx->offset_x-(lx + 1) * sand_ctx->extra_x)/lx;
	}

	if(sand_ctx->offset_y < sand_ctx->extra_y*(ly+1)){
		start_y = sand_ctx->offset_y/(ly+1);
	}else{
		start_y = sand_ctx->extra_y+(sand_ctx->offset_y-(ly + 1)* sand_ctx->extra_y)/(ly);
	}

    // ### Determine the global rank (start_rank) that receives the initial grain 
	int start_coords[2] = {start_x, start_y};
    MPI_Cart_rank(sand_ctx->comm_2d, start_coords, &sand_ctx->start_rank);

    //Initial configuration: height and max iterations
	sand_ctx->start_height = 2*sand_ctx->L*sand_ctx->L;
	sand_ctx->max_iterations = 5*sand_ctx->start_height/4;
	sand_ctx->max_iterations = (sand_ctx->max_iterations > MAX_LOOPS)?MAX_LOOPS:sand_ctx->max_iterations;
    //Initialize starting height
	if(sand_ctx->new_rank == sand_ctx->start_rank){
        int local_offset_x = sand_ctx->offset_x;
        int local_offset_y = sand_ctx->offset_y;

		// Calculate the local offset relative to this processor's starting point
        // recalculate the global offset based on the uneven boundaries
		if(sand_ctx->extra_x != 0 && sand_ctx->coords[0] >= sand_ctx->extra_x){
			local_offset_x -= (sand_ctx->LX + 1)*sand_ctx->extra_x;
		}
		if(sand_ctx->extra_y != 0 && sand_ctx->coords[1] >= sand_ctx->extra_y){
			local_offset_y -= (sand_ctx->LY+1)*sand_ctx->extra_y;
		}

		 // Setting the starting height in the local grid +1 for halos
		sand_ctx->old_height[local_offset_x % sand_ctx->LX + 1][local_offset_y % sand_ctx->LY + 1] = sand_ctx->start_height;
	}
}


// Setup initialization function: ====================================================================================
void sandpile_initialization(SandpileContext *sand_ctx, int argc, char *argv[]){
    sand_ctx->comm = MPI_COMM_WORLD;
    MPI_Comm_rank(sand_ctx->comm, &sand_ctx->rank);
    MPI_Comm_size(sand_ctx->comm, &sand_ctx->size);

    // Input Parsing and Bcast
    parse_bcast_input(sand_ctx, argc, argv);

    //Settingup topology
    setup_topology_and_domain(sand_ctx);

    //datatypes
    setup_datatypes(sand_ctx);

    setup_arrays_and_sandpile(sand_ctx);

    sand_ctx->num_loops = 0;
    sand_ctx->tloop = 0.0;
}















