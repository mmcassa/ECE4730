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
    int *dim;
    if (size < 4) {
        ndim = 1;
        dim = (int *) calloc(1,sizeof(int));
    } else {
        ndim = 2;
        dim = (int *) calloc(2,sizeof(int));
    }
}

void compute() {

}

