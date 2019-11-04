/* Mitchell Cassaday & William (Luke) Benton
 * ECE4730
 * Fall 2019
 * 4 NOV 2019
 * Goal: Apply floyd's algorithm using mpi
 */

#include "graph.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>


void distribute (
   MPI_Datatype dtype,   /* IN - Element type */
   int n,               /* IN - Array cols */
   int ***local_array,  /* OUT - 2D Array */
   int **full_matrix,   /* IN - Full 2D Array from File */
   MPI_Comm grid_comm
   );

void compute(   MPI_Comm grid,      // Grid comm
                MPI_Comm row,       // Row comm
                MPI_Comm cols,      // Cols comm
                int rank,           // Process rank
                int *coords, 
                int *dim, 
                int n, 
                int size, 
                int **loc_matrix
                );

int main(int argc, char* argv[]) {
    int i,j,k,n;
    int dims[2] = {0,0};          // [# of rows, # of cols]
    
    int periodic[2];
    periodic[0] = 0;
    periodic[1] = 0;
    int coords[2] = {0,0};
    int grid_coord[2] = {0,0};  /* Process coords */
    int grid_period[2] = {0,0}; /* Wraparound */
    int grid_size[2] = {0,0};   /* Dimensions of grid */
    int **A;
    int **loc_matrix;
    char *file_in = calloc(100,sizeof(char));
    char *file_out = calloc(100,sizeof(char));
    
    // Initialize the MPI environment
	MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
	// Find out rank, size
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    

    /* Create the grid of processes to communicate over */
    MPI_Comm comm_grid;
    MPI_Dims_create(world_size,2,dims);
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periodic,0,&comm_grid);
    MPI_Cart_coords(comm_grid,world_rank,2,coords);

    /* Divide the grid into rows and cols */
    MPI_Comm row_comm,col_comm;
    MPI_Comm_split(comm_grid,coords[0],coords[1],&row_comm);
    MPI_Comm_split(comm_grid,coords[1],coords[0],&col_comm);

    if (world_rank == 0) {
        // Parse out arguemnets
        if (argc != 3) {
            printf("Usage:  print-graph {file_name}\n");
            exit(1);
        }
        file_in = argv[1];
        file_out = argv[2];

        // Read in data
        read_graph(file_in,&n,&A);
        MPI_Bcast (&n, 1, MPI_INT, 0, comm_grid);
    }
    distribute(MPI_INT,n,&loc_matrix,A,comm_grid);


    /* Begin Parallel Portion of the program */
    compute(comm_grid,row_comm,col_comm,world_rank,coords,dims,n,world_size,loc_matrix);

    // 
    MPI_Cart_get(comm_grid, 2, grid_size, grid_period,grid_coord);
    MPI_Datatype rowType,arrType;
    MPI_Type_contiguous(n, MPI_INT, &rowType);
    MPI_Type_commit(&rowType);
    MPI_Type_vector(n*n,1,0, rowType, &arrType);
    MPI_Type_commit(&arrType);
    
    if(world_rank == 0) {

        MPI_Gather(&loc_matrix, 1, arrType, A, 1, A, 0,comm_grid);
    } else {
        MPI_Gather(&loc_matrix, 1, arrType, NULL, 1, A, 0,comm_grid);
    }
    print_graph(n,A);
    //write_graph(file_out,n,A);
    return 0;
}

void gather(    MPI_Comm grid,      // Grid comm
                MPI_Comm row,       // Row comm
                MPI_Comm cols,      // Cols comm
                int rank,           // Process rank
                int *coords, 
                int *dim, 
                int n, 
                int size, 
                int **loc_matrix
                ) {

}

void compute(   MPI_Comm grid,      // Grid comm
                MPI_Comm row,       // Row comm
                MPI_Comm cols,      // Cols comm
                int rank,           // Process rank
                int *coords, 
                int *dim, 
                int n, 
                int size, 
                int **loc_matrix
                ) {
    int k,i,j,l;      // Loop values
    int sender;     // Rank of Broadcaster
    int grid_coord[2] = {0,0};  /* Process coords */
    int grid_period[2] = {0,0}; /* Wraparound */
    int grid_size[2] = {0,0};   /* Dimensions of grid */
    int *k_row;
    int k_col;
    int rel_r, rel_c;

    MPI_Cart_get(grid, 2, grid_size, grid_period,grid_coord);
    
    k_row = (int *) calloc(BLOCK_SIZE(grid_coord[0],grid_size[0],n),sizeof(int));
    
    for(k=0;k<n;k++) {
        
        sender = k/(n/size)*grid_size[0]+coords[1];
        rel_r = k % (grid_size[0]*BLOCK_SIZE(grid_coord[0],grid_size[0],n));
        if (sender == rank) {
            k_row = loc_matrix[rel_r];
        }
        
        MPI_Bcast((void *) k_row, n, MPI_INT, sender, cols);
        for (i=0;i<BLOCK_SIZE(grid_coord[0],grid_size[0],n);i++) {
            //printf("%d\t %d %d %d %d %d\n",rank,k,n,size,grid_size[0],coords[1]);
            sender = k/(n/size)*grid_size[0]+coords[1];
            rel_c = k % ((coords[1]+1)*BLOCK_SIZE(grid_coord[0],grid_size[0],n));
            if (sender == rank) {
                k_col = loc_matrix[i][rel_c];
            }
            //printf("%d %d %d\t %d %d %d\n",rel_r,rel_c,k_col,k,i,j);
            MPI_Bcast(&k_col, n, MPI_INT, sender, row);
            //printf("\t%d\n",k_col);
            for (j=0;j<BLOCK_SIZE(grid_coord[1],grid_size[1],n);j++) {
                if (k_col > 0 && k_row[j] > 0) {
                    if ((k_col + k_row[j]) < loc_matrix[i][j] || loc_matrix[i][j] == -1) {
                        loc_matrix[i][j] = k_col + k_row[j]; 
                    }
                }
            }
        } 
    }


}

void distribute (
   MPI_Datatype dtype,   /* IN - Element type */
   int n,               /* IN - Array cols */
   int ***local_array,  /* OUT - 2D Array */
   int **full_matrix,   /* IN - Full 2D Array from File */
   MPI_Comm grid_comm)   /* IN - Communicator */
{
    void      *buffer;         /* File buffer */
    int        coords[2];      /* Coords of proc receiving
                                    next row of matrix */
    int        datum_size;     /* Bytes per elements */
    int        dest_id;        /* Rank of receiving proc */
    int        grid_coord[2];  /* Process coords */
    int        grid_id;        /* Process rank */
    int        grid_period[2]; /* Wraparound */
    int        grid_size[2];   /* Dimensions of grid */
    int        i, j, k;
    void      *laddr;          /* Used when proc 0 gets row */
    int        local_cols;     /* Matrix cols on this proc */
    int        local_rows;     /* Matrix rows on this proc */
    void     **lptr;           /* Pointer into 'subs' */
    int        p;              /* Number of processes */
    void      *raddr;          /* Address of first element
                                    to send */
    void      *rptr;           /* Pointer into 'storage' */
    MPI_Status status;         /* Results of read */

    MPI_Comm_rank (grid_comm, &grid_id);
    MPI_Comm_size (grid_comm, &p);
    datum_size = sizeof(int);


    /* Each process determines the size of the submatrix
        it is responsible for. */

    MPI_Cart_get (grid_comm, 2, grid_size, grid_period,
        grid_coord);
    local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],n);
    local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],n);

    /* Dynamically allocate two-dimensional matrix 'subs' */

    local_array[0] = (int **) calloc(local_rows,sizeof(int *));
    for (i=0;i<local_rows;i++) {
        local_array[0][i] = (int *) calloc(local_cols,sizeof(int));
    }


   /* Grid process 0 reads in the matrix one row at a time
      and distributes each row among the MPI processes. */

   if (grid_id == 0)
      buffer = (int *) calloc(n,sizeof(int));

   /* For each row of processes in the process grid... */
   for (i = 0; i < grid_size[0]; i++) {
      coords[0] = i;

      /* For each matrix row controlled by this proc row...*/
      for (j = 0; j < BLOCK_SIZE(i,grid_size[0],n); j++) {

         /* Read in a row of the matrix */

         if (grid_id == 0) {
             buffer = full_matrix[i*grid_size[0]+j];
            //fread (buffer, datum_size, *n, infileptr);
         }

         /* Distribute it among process in the grid row */

         for (k = 0; k < grid_size[1]; k++) {
            coords[1] = k;

            /* Find address of first element to send */
            raddr = buffer +
               BLOCK_LOW(k,grid_size[1],n) * datum_size;

            /* Determine the grid ID of the process getting
               the subrow */
            MPI_Cart_rank (grid_comm, coords, &dest_id);

            /* Process 0 is responsible for sending...*/
            if (grid_id == 0) {

               /* It is sending (copying) to itself */
               if (dest_id == 0) {
                    memcpy (local_array[0][j], raddr,
                        local_cols * datum_size);

               /* It is sending to another process */
               } else {
                  MPI_Send (raddr,
                     BLOCK_SIZE(k,grid_size[1],n), dtype,
                  dest_id, 0, grid_comm);
               }

            /* Process 'dest_id' is responsible for
               receiving... */
            } else if (grid_id == dest_id) {
               MPI_Recv (local_array[0][j], local_cols, dtype, 0,
                  0, grid_comm,&status);
            }
         }
      }
   }
   if (grid_id == 0) free (buffer);
}