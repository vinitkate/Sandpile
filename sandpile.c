#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<stdbool.h>
#include<unistd.h>

#include "sandpile_context.h"

#define ROOT_RANK 0

/*
main sandpile function:
Calls subroutines 
- initialize the sandpile
- Run the simulation steps to update the sandpile
- Writing the output 
*/

int main(int argc, char *argv[]){
    SandpileContext sand_ctx;

    double tstart, tend;
    int loop_index;

    // MPI startup
    MPI_Init(&argc, &argv);

    //Start timing
    MPI_Barrier(MPI_COMM_WORLD);
    tstart = MPI_Wtime();

    // 1.--- initialization: sandpile_init.c
    sandpile_initialization(&sand_ctx, argc, argv);

    // 2.--- Main simulation loop: sandpile_simulation.c
    for(loop_index = 0; loop_index < sand_ctx.max_iterations; loop_index++){
        int total_change;

        total_change = run_simulation_step(&sand_ctx);
		
	if(total_change == 0){
            break;
        }
    }

    //Timing final
    MPI_Barrier(MPI_COMM_WORLD);
    tend = MPI_Wtime()-tstart;

    // 3.---Write the output to file and done: sandpile_io.c
    finalize_and_cleanup(&sand_ctx, tend);

    MPI_Finalize();

    return EXIT_SUCCESS;

}
