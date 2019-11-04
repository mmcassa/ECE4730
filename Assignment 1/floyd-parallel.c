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


    // Parse out arguemnets
    if (argc != 3) {
        printf("Usage:  print-graph {file_name}\n");
        exit(1);
    }

    file_in = argv[1];
    file_out = argv[2];

    // Read in data
    read_graph(file_in,&n,&A);
    print_graph(n,A);

    /*floyd = (int **) calloc(n,sizeof(int *));
    for (i=0;i<n;i++) {
        floyd[i] = (int *) calloc(n,sizeof(int));
    }*/

    /* Begin Parallel Portion of the program */
    // https://www.palmetto.clemson.edu/jupyterhub/hub/spawn
    // Initialize the MPI environment
	MPI_Init(NULL, NULL);
	// Find out rank, size
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);



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

