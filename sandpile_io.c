#include <stdio.h>
#include <stdlib.h>
#include "sandpile_context.h"
#include "arraymalloc.h" 

#define HLINE "-----------------------------------------------------------"
#define ROOT_RANK 0

/*
I/O functions
- samdwritedynamic(char *, int **, int): writes output to file .ppm
- gather_result(SandpileContext *): Collects result to rank 0 from all processors
- write_time_and_cleanup(SandpileContext *, double): Write timing to a file time.txt
													clean up the dynamically allocated arrays
---------------------------------------------------------------------------------------------------------
- finalize_and_cleanup(SandpileContext *, double): Public API - function for I/O  
--------------------------------------------------------------------------------------------------------
*/

//Function to write the output to file
void sandwritedynamic(char *sandfile, int **height, int l)
{
  FILE *fp;

  int i, j, npix, c;
  static int pixperline = 12;

  // Encoding is black, green, magenta, yellow

  int col[MAXH+1][NRGB] = { {0,0,0}, {0,1,0}, {1,0,1}, {1,1,0} };

  /*
   *  Write the file
   */

  printf("sandwritedynamic: opening file <%s>\n", sandfile);

  fp = fopen(sandfile, "w");

  printf("sandwritedynamic: writing data ...\n");

  /*
   *  Start with the PPM header
   */

  fprintf(fp, "P3\n");
  fprintf(fp, "# Written by sandwritedynamic\n");
  fprintf(fp, "%d %d\n", l, l) ;
  fprintf(fp, "%d\n", 1) ;

  /*
   *  Now write the heights to file so that height[0][0] is in the
   *  bottom-left-hand corner and height[l-1][l-1] is in the
   *  top-right-hand corner
   */

  npix = 0;

  for (j=l-1; j >= 0; j--)
    {
      for (i=0; i < l; i++)
	{
	  npix++;

          c = height[i][j];

          if (c > MAXH) c = MAXH;

	  // Make sure lines wrap after "pixperline" pixels

	  if (npix == 1)
	    {
              fprintf(fp, "%1d %1d %1d", col[c][0], col[c][1], col[c][2]);
	    }
	  else if (npix < pixperline)
	    {
              fprintf(fp, " %1d %1d %1d", col[c][0], col[c][1], col[c][2]);
	    }
	  else
	    {
              fprintf(fp, " %1d %1d %1d\n", col[c][0], col[c][1], col[c][2]);
	      npix = 0;
	    }
	}
    }

  if (npix != 0) fprintf(fp, "\n");

  printf("sandwritedynamic: ... done\n");

  fclose(fp);
  printf("sandwritedynamic: file closed\n");
}

static void write_matrix(int **matrix, int rows, int cols, const char *fileName){
	FILE *fp = fopen(fileName, "w");
	int i, j;
	if(fp == NULL){
		fprintf(stderr, "failed to write matrix to file\n");
		return;
	}

	fprintf(fp, "%d %d\n", rows, cols);
	for(i = 0; i < rows; i++){
		for(j = 0; j < cols; j++){
			fprintf(fp, "%d", matrix[i][j]);
			if(j < cols-1)fprintf(fp, " ");
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	return;
}

//Helper function to gather results
static void gather_results(SandpileContext *sand_ctx){
    if(sand_ctx->new_rank == ROOT_RANK){
        int i, j;
        int stx, sty, lx, ly;
        int **height;
        int coords[2];

        // Datatype for receiving rows in the gathering phase (LY-1 width, for smaller size due non uniform distribution)
		MPI_Datatype ROW_DATA_TYPE_2;
		MPI_Type_contiguous(sand_ctx->LY-1, MPI_INT, &ROW_DATA_TYPE_2);
    	MPI_Type_commit(&ROW_DATA_TYPE_2);

        // Allocate array of requests for all incoming row receives (L global rows x PY processors)
		MPI_Request **reqs = (MPI_Request**)arraymalloc2d(sand_ctx->L, sand_ctx->PY, sizeof(MPI_Request));
        // Allocate global result array
		height = (int **)arraymalloc2d(sand_ctx->L, sand_ctx->L, sizeof(int));

        // Copy rank 0's local old_height to the global array's (height) starting position
		for(i = 0; i < sand_ctx->LX; i++){
			for(j = 0; j < sand_ctx->LY; j++){
				height[i][j] = sand_ctx->old_height[i+1][j+1];
			}
		}

        for(i = 1; i < sand_ctx->size; i++){
			MPI_Cart_coords(sand_ctx->comm_2d, i, 2, coords); // Get coordinates of source rank i

			// Calculate starting position (stx, sty) and local size (lx, ly) of the source domain
			lx = sand_ctx->LX;
			if(sand_ctx->extra_x == 0 || coords[0] < sand_ctx->extra_x){
				stx = (sand_ctx->LX) * coords[0];
			}else{
				stx =(sand_ctx->LX)*(sand_ctx->extra_x) + (coords[0]- sand_ctx->extra_x)*(sand_ctx->LX-1);
				lx--; // Reduce local X size i.e. LX for these ranks
			}

			ly = sand_ctx->LY;
			if(sand_ctx->extra_y == 0 || coords[1] < sand_ctx->extra_y){
				sty = sand_ctx->LY * coords[1];
			}else{
				sty = (sand_ctx->LY)*(sand_ctx->extra_y) + (coords[1]-sand_ctx->extra_y)*(sand_ctx->LY-1);
				ly--; // Reduce local Y size i.e. LY for these ranks
			}

			//Label for receiving message
			int tag = 1;

			// Receive rows into the global height array
			for(j = stx; j < stx+lx; j++){
				// Check which datatype to use based on the source rank's domain size
				if(sand_ctx->extra_y == 0 || coords[1] < sand_ctx->extra_y){
					// Use the type for LY elements (for larger domains)
					MPI_Irecv(&height[j][sty], 1, sand_ctx->ROW_DATA_TYPE, i, tag+i, sand_ctx->comm_2d, &reqs[j][coords[1]]);
				}else{
					// Use the type for LY-1 elements (for smaller domains)
					MPI_Irecv(&height[j][sty], 1, ROW_DATA_TYPE_2, i, tag+i, sand_ctx->comm_2d, &reqs[j][coords[1]]);
				}
				tag++;
			}
		}
        // Wait for all receives to complete. This logic ensures all row segments are received.
		//For first LX rows processors 0 is not included
		for(i = 0; i < sand_ctx->L; i++){
			if(i < sand_ctx->LX){
				MPI_Waitall(sand_ctx->PY - 1, &reqs[i][1], MPI_STATUSES_IGNORE);
			}
			else{
				MPI_Waitall(sand_ctx->PY, &reqs[i][0], MPI_STATUSES_IGNORE);
			}
		}
	char fileName[200];
	snprintf(fileName, 200, "sand_%d_%d.mat", sand_ctx->size, sand_ctx->L);
    	write_matrix(height, sand_ctx->L, sand_ctx->L, fileName);
	
	snprintf(fileName, 200, "sand_%d_%d.ppm", sand_ctx->size, sand_ctx->L);
        sandwritedynamic(fileName, height, sand_ctx->L);

	// Clean up root-specific allocations
		free(height);
		free(reqs);
        MPI_Type_free(&ROW_DATA_TYPE_2);
    }else{
        // ### Sender Ranks (new_rank != 0) ---
		// Allocate requests for sending local data back to rank 0 of size num rows
		int i;
		MPI_Request *req = (MPI_Request *)malloc(sand_ctx->LX*sizeof(MPI_Request)); 

		// Send local inner rows (1 to LX) to rank 0
		for(i = 1; i <= sand_ctx->LX; i++){
			// Using ROW_DATA_TYPE for the send 
			//(this is based on the assumption that rank 0 handles 
			// the type matching based on receiving processor's dimension)
			MPI_Isend(&sand_ctx->old_height[i][1], 1, sand_ctx->ROW_DATA_TYPE, 0, sand_ctx->new_rank + i, sand_ctx->comm_2d, &req[i-1]);
		}
		MPI_Waitall(sand_ctx->LX, req, MPI_STATUSES_IGNORE);

		free(req);
    }
}

//Helper function to write final timings and free dynamically allocated space
static void write_time_and_cleanup(SandpileContext *sand_ctx, double tend){
    if(sand_ctx->new_rank == ROOT_RANK){
        FILE *fp;
        if((fp = fopen("time.txt", "a+")) != NULL){
			double time_per_loop = sand_ctx->tloop / (1.0*sand_ctx->num_loops);
			fprintf(fp, "Processors %d, Problem size: %d\n", sand_ctx->size, sand_ctx->L);
			fprintf(fp, "Total time: %.20f \n", tend);
			fprintf(fp, "Total loops time: %.20f \n", sand_ctx->tloop);
			fprintf(fp, "Total loops: %d \n", sand_ctx->num_loops);
			fprintf(fp, "one loop time %.20f \n", time_per_loop);
			fprintf(fp, HLINE);
			putc('\n', fp);
			fclose(fp);
		}
    }


    // --- Cleanup Shared Allocations ---
	free(sand_ctx->left_col_recv);
	free(sand_ctx->left_col_send);
	free(sand_ctx->right_col_recv);
	free(sand_ctx->right_col_send);

    free(sand_ctx->old_height);
	free(sand_ctx->new_height);

	MPI_Type_free(&sand_ctx->ROW_DATA_TYPE);
    MPI_Type_free(&sand_ctx->CONT_DATA_FOR_COLUMN);
}

//Function finalize_and_cleanup
void finalize_and_cleanup(SandpileContext *sand_ctx, double tend){
    gather_results(sand_ctx);

    write_time_and_cleanup(sand_ctx, tend);
}
