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

int main(int argc, char* argv[]) {
    int i,j,k,n;
    int *dims;      // [# of rows, # of cols]
    int periodic[2];
    periodic[0] = 0;
    periodic[1] = 0;
    int *coords;
    int **A;
    int **floyd;
    char *file_in = calloc(100,sizeof(char));
    char *file_out = calloc(100,sizeof(char));
    
    // Initialize the MPI environment
	MPI_Init(NULL, NULL);
	// Find out rank, size
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    /* Create the grid of processes to communicate over */
    MPI_Comm comm_grid;
    MPI_Comm row_comm,col_comm;
    MPI_Dims_create(world_size,2,dims);
    MPI_Cart_create(MPI_WORLD_COMM,2,dims,periodic,0,&comm_grid);
    MPI_Cart_coords(comm_grid,world_rank,2,coords);
    MPI_Comm_split(comm_grid,coords[0],coords[1],&row_comm);
    

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
        distribute(world_size,A);
    }
    //print_graph(n,A);


    /* Begin Parallel Portion of the program */
    
    



    for(k=0;k<n;k++) {
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                if (A[i][k] > 0 && A[k][j] > 0) {
                    if ((A[i][k]+A[k][j]) < A[i][j] || A[i][j] == -1) {
                        A[i][j] = A[i][k]+A[k][j]; 
                    }
                }
            }
        }
    }
    printf("\n\n\n");
    print_graph(n,A);
    write_graph(file_out,n,A);
    return 0;
}

/* 
    Create the Cartesian space for 1, 2, 4, 8, 16, 32, and 64 processes
*/
void distribute(int size, int **matrix) {
    int ndim;
    int *dim = (int *) calloc(2,sizeof(int));
    if (size < 4) {
        ndim = 1;
        dim = (int *) calloc(1,sizeof(int));
        dim[0] = size;
        dim[1] = 0;
    } else {
        ndim = 2; 
        switch (size) {

            case 4:
                dim[0] = 2;
                dim[1] = 2;
                break;

            case 8:
                dim[0] = 4;
                dim[1] = 2;
                break;

            case 16:
                dim[0] = 4;
                dim[1] = 4;
                break;
            
            case 32:
                dim[0] = 8;
                dim[1] = 4;
                break;

            case 64:
                dim[0] = 8;
                dim[1] = 8;
                break;

            default:
                printf("Please use 1,2,4,8,16,32, or 64 processors\n");
                exit(1);
                break;
        }
    }
    
}

void compute(MPI_Comm comm, int rank, int dim[], int n, int size, int locMatrix[]) {
    int k,i,j;      // Loop values
    int sender;     // Rank of Broadcaster
    int *k_row = (int *) calloc(n*n/p,sizeof(int));
    int k_col;

    
    for(k=0;k<n;k++) {
        sender = k/n/p;
    }


}

