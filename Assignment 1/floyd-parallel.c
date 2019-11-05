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


int findSource(int *grid_size,int *grid_coords,int k,int n);
int hasRow(int k,int *coords,int n,int *grid_size,int *grid_coord);
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
    MPI_Status status;
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
        //print_graph(n,A);
        
    }
    MPI_Bcast (&n, 1, MPI_INT, 0, comm_grid);
    //printf("%d %d\n",world_rank,n);
    distribute(MPI_INT,n,&loc_matrix,A,comm_grid);
    

    /* Begin Parallel Portion of the program */
    compute(comm_grid,row_comm,col_comm,world_rank,coords,dims,n,world_size,loc_matrix);

    // 
    MPI_Cart_get(comm_grid, 2, grid_size, grid_period,grid_coord);
    int local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],n);
    int local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],n);
    int send_coords[2] = {0,0};
    int dest_id;
    int **temp_arr;
    int *red_arr;
    int isColZ;
    int t;
    /* Loop through every row */
    for (k=0;k<n;k++) {
        /* If process contains row, then gatherv to the root of the row */
        isColZ = hasRow(k,coords,n,grid_size,grid_coord);
        if (isColZ) {
            /* If task is root of row, gather */
            if (0 == coords[1]) {
                //temp_arr = calloc(n,sizeof(int));
                red_arr = (int *) calloc(n,sizeof(int));
                MPI_Gather(&loc_matrix[i],local_cols,MPI_INT,temp_arr,1,MPI_INT,world_rank,row_comm);
                t = 0;
                for(i=0;i<grid_size[1];i++) {
                    t += BLOCK_SIZE(i,grid_size[1],n);
                    for (j=0;j<BLOCK_SIZE(i,grid_size[1],n);j++) {
                        red_arr[t+j] = temp_arr[i][j];
                    }
                }
            } else {
                /* Else send to gatherv */
                MPI_Gather(loc_matrix[i],local_cols,MPI_INT,temp_arr,1,MPI_INT,(world_rank - world_rank % grid_size[1]),row_comm);
            }
        }
        /* If task is world 0 RECV with tag k and store in A[k] */
        if (world_rank == 0) {
            if (isColZ) {
                memcpy(A[k],temp_arr,sizeof(int)*n);
            } else {
                MPI_Recv(A[k],n,MPI_INT,findSource(grid_size,grid_coord,k,n),k,comm_grid,&status);
            }
        /* If task has the gathered row, SEND with tag k to root */
        } else if (coords[1] == 0 && isColZ) {
            MPI_Send(temp_arr,n,MPI_INT,0,k,comm_grid);
        }
    }


    if (world_rank == 0) {
        print_graph(n,A);
    }
    int error = 1; 
    MPI_Bcast(&error, 1, MPI_INT, 0, comm_grid);
    if (error != 0) {
        if( world_rank == 0 ) {
            
            
            MPI_Finalize();
            //print_graph(n,loc_matrix);
            //write_graph(file_out,n,A);
        }
    }
    return 0;
}

int findSource(int *grid_size,int *grid_coords,int k,int n) {
    int i=0,ret;
    while(ret < k) {
        ret += BLOCK_SIZE(i++,grid_size[0],n);
    }
    return i*grid_size[1];
}
int hasRow(int k,int *coords,int n,int *grid_size,int *grid_coord) {
    int rows = 0;
    int i,j;
    int *loc_coord;
    memcpy(loc_coord,coords,2*sizeof(int));
    for(i=coords[0-1];i>=0;i--) {
        rows += BLOCK_SIZE(i,grid_size[0],n);
    }
    if (k >= rows && k < (rows+BLOCK_SIZE(grid_coord[0],grid_size[0],n)))
        return k-rows;
    return -1;
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
    int send_coords[2] = {0,0};
    int grid_period[2] = {0,0}; /* Wraparound */
    int grid_size[2] = {0,0};   /* Dimensions of grid */
    int *temp_k;
    int *k_row;
    int k_col = 0;
    int rel_r, rel_c;

    MPI_Cart_get(grid, 2, grid_size, grid_period,grid_coord);
    
    
    int local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],n);
    int local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],n);
    k_row = (int *) calloc(local_cols,sizeof(int));
    temp_k = (int *) calloc(local_cols,sizeof(int));

    printf("Rank %d %d %d\n",rank,local_rows,local_cols);
    printf("\n");
    for(k=0;k<n;k++) {
        send_coords[0] = k/(n/grid_size[0]);
        send_coords[1] = coords[1];
        
        MPI_Cart_rank (grid, send_coords, &sender);
        
        rel_r = k % (local_rows);
        
        if (sender == rank) {
            //k_row = loc_matrix[rel_r];
            memcpy(temp_k,loc_matrix[rel_r],local_cols*sizeof(int));//k_row[l] = 0;
            printf("SENDER %d %d\n",sender,rel_r);
            //printf("Sender: %d %d\n",rank,rel_r);
            MPI_Allreduce((void *) temp_k,(void *) k_row,local_cols,MPI_INT,MPI_SUM,cols);
        } else {
            for (l=0;l<local_cols;l++) {
                temp_k[l] = 0;
            }
            MPI_Allreduce((void *) temp_k,(void *) k_row,local_cols,MPI_INT,MPI_SUM,cols);
        }

        if (sender != rank) {
            /*printf("K_row (rank %d): ",rank);
            for (l=0;l<local_cols;l++) {
                printf("%d ", k_row[l]);
            }
            printf("\n");*/
        }
        //MPI_Bcast((void *) k_row, n, MPI_INT, sender, cols);
        
        for (i=0;i<BLOCK_SIZE(grid_coord[0],grid_size[0],n);i++) {
            //printf("%d\t %d %d %d %d %d\n",rank,k,n,size,grid_size[0],coords[1]);
            //sender = k/(n/size)*grid_size[0]+coords[1];
            send_coords[0] = coords[0];
            send_coords[1] = k/(n/grid_size[1]);
            MPI_Cart_rank(grid,send_coords,&sender);

            rel_c = k % (local_cols);
            if (sender == rank) {
                memcpy(&l,&(loc_matrix[i][rel_c]),sizeof(int));
                //printf("Rank: %d %d %d %d\t %d %d %d\n",rank,rel_r,rel_c,l,k,i,j);
            } else {
                l = 0;
            }
            
            //MPI_Bcast((void *) &k_col, 1, MPI_INT, sender, row);
            MPI_Allreduce(&l,&k_col,1,MPI_INT,MPI_SUM,row);

            //if (sender != rank) {
                //printf("k_col: \t%d ",k_col);
            //}
            for (j=0;j<BLOCK_SIZE(grid_coord[1],grid_size[1],n);j++) {
                if (k_col > 0 && k_row[j] > 0) {
                    if ((k_col + k_row[j]) < loc_matrix[i][j] || loc_matrix[i][j] == -1) {
                        loc_matrix[i][j] = k_col + k_row[j]; 
                    }
                }
            }
        } printf("\n");
    }


}

void distribute (
   MPI_Datatype dtype,   /* IN - Element type */
   int n,               /* IN - Array cols */
   int ***local_array,  /* OUT - 2D Array */
   int **full_matrix,   /* IN - Full 2D Array from File */
   MPI_Comm grid_comm) {
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